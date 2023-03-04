#include "Includes.h"

#ifdef ENABLE_CUSTOM_FORMULA
#define _USE_MATH_DEFINES
#include <math.h>
#include <ranges>
#include <stack>

#define and &&
#define or ||
#define not !
#include <symengine/basic.h>
#include <symengine/visitor.h>
#include <symengine/solve.h>
#include <symengine/eval_mpc.h>
#undef and
#undef or
#undef not

#include "Encoder/Encoder.h"
#include "Encoder/XMMComplex.h"

class FormulaCompileError : public std::exception {
public:
	explicit FormulaCompileError(const std::string &Message) : std::exception(Message.c_str()) {}
	explicit FormulaCompileError(const char *Message) : std::exception(Message) {}
};

struct Operation {
	enum OpEnum : uint16_t {
		Nop,

		LoadConst,
		LoadZ,
		LoadC,

		Add,
		Sub,
		Mul,
		Div,
		
		Exp,
		PowUint,
		Log,

		Pow,

		Sin,
		Cos,
		Tan,
	};

	OpEnum operation;
	bool reference = true;
	uint32_t parameter = 0;

	Operation(OpEnum operation) : operation(operation) {}
	Operation(OpEnum operation, bool reference) : operation(operation), reference(reference) {}
	Operation(OpEnum operation, uint32_t parameter) : operation(operation), parameter(parameter) {}
	Operation(OpEnum operation, uint32_t parameter, bool reference) : operation(operation), reference(reference), parameter(parameter) {}
};

struct Constant {
	SymEngine::RCP<const SymEngine::Basic> Value;

	Constant(const SymEngine::RCP<const SymEngine::Number> &x) : Value(x) {}
	Constant(const SymEngine::RCP<const SymEngine::Constant> &x) : Value(x) {}
	Constant(const SymEngine::Number &x) : Value(x.rcp_from_this()) {}
	Constant(const SymEngine::Constant &x) : Value(x.rcp_from_this()) {}
};

struct CompiledFormula {
	std::vector<Operation> Operations;
	std::vector<Constant> Constants;

	size_t AddOperation(Operation operation) {
		Operations.push_back(operation);
		return Operations.size() - 1;
	}
	size_t AddOperation(Operation::OpEnum operation) {
		Operations.emplace_back(operation);
		return Operations.size() - 1;
	}
	size_t AddOperation(Operation::OpEnum operation, uint32_t parameter) {
		Operations.emplace_back(operation, parameter);
		return Operations.size() - 1;
	}

	size_t AddConstant(Constant constant) {
		Constants.push_back(constant);
		return Constants.size() - 1;
	}
};

class EncodeVisitor : public SymEngine::BaseVisitor<EncodeVisitor> {
protected:
	XMMComplexEncoder &encoder;
	uint8_t &stack;

public:
	RegisterXMM Top() {
		return RegisterXMM(stack - 1);
	}

	RegisterXMM Pop() {
		stack--;
		return RegisterXMM(stack);
	}

	void EmitPow(uint32_t exp) {
		if (exp == 1) {
		} else if (exp == 2) {
			encoder.Csqr(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
		} else if (exp == 3) {
			encoder.MovPD(RegisterXMM(stack), RegisterXMM(stack - 1));
			encoder.Csqr(RegisterXMM(stack), RegisterXMM(stack + 1), RegisterXMM(stack + 2));
			encoder.Cmul(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
		} else if (exp == 4) {
			encoder.Csqr(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
			encoder.Csqr(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
		} else throw;
	}

	void EmitMul() {
		stack--;
		encoder.Cmul(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
	}

	void EmitDiv() {
		stack--;
		encoder.Cdiv(RegisterXMM(stack - 1), RegisterXMM(stack), RegisterXMM(stack + 1));
	}

	void EmitCall1(void *func) {
		using enum Register64;
		encoder.Sub(rsp, 32 + 2 * sizeof(__m128));
		encoder.Lea(rcx, rsp + 32);
		encoder.Lea(rdx, rsp + (32 + sizeof(__m128)));
		encoder.MovPD(*rdx, Top());

		encoder.Call(func);

		encoder.MovPD(Top(), *rax);
		encoder.Add(rsp, 32 + 2 * sizeof(__m128));
	}

	void EmitCall2(void *func) {
		using enum Register64;
		encoder.Sub(rsp, 32 + 3 * sizeof(__m128));
		encoder.Lea(rcx, rsp + 32);
		encoder.Lea(rdx, rsp + (32 + sizeof(__m128)));
		encoder.Lea(r8, rsp + (32 + sizeof(__m128) * 2));
		encoder.MovPD(*r8, Pop());
		encoder.MovPD(*rdx, Top());

		encoder.Call(func);

		encoder.MovPD(Top(), *rax);
		encoder.Add(rsp, 32 + 3 * sizeof(__m128));
	}

	EncodeVisitor(XMMComplexEncoder &encoder, uint8_t &stack) : encoder(encoder), stack(stack) {}

	void bvisit(const SymEngine::Basic &x) {
		throw;
	}

	void bvisit(const SymEngine::Number &x) {
		const auto value = SymEngine::eval_complex_double(x);
		if (stack >= 14) throw;
		encoder.MovPD(RegisterXMM(stack), *encoder.Const(_mm_set_pd(value.imag(), value.real())));
		stack++;
	}
	void bvisit(const SymEngine::Constant &x) {
		const auto value = SymEngine::eval_complex_double(x);
		if (stack >= 14) throw;
		encoder.MovPD(RegisterXMM(stack), *encoder.Const(_mm_set_pd(value.imag(), value.real())));
		stack++;
	}
	void bvisit(const SymEngine::Symbol &x) {
		const auto name = x.get_name();
		if (name == "Z" || name == "z") {
			if (stack >= 14) throw;
			encoder.MovPD(RegisterXMM(stack), *Register64::rdi);
			stack++;
		} else if (name == "C" || name == "c") {
			if (stack >= 14) throw;
			encoder.MovPD(RegisterXMM(stack), *Register64::rsi);
			stack++;
		} else throw;
	}
	void bvisit(const SymEngine::Add &x) {
		const auto args = x.get_args();

		apply(args.front());

		for (const auto &p : args | std::views::drop(1)) {
			apply(p);
			stack--;
			encoder.AddPD(RegisterXMM(stack - 1), RegisterXMM(stack));
		}
	}
	void bvisit(const SymEngine::Mul &x) {
		apply(x.get_coef());

		for (const auto &p : x.get_dict()) {
			apply(p.first);

			if (neq(*p.second, *SymEngine::one)) {
				if (!SymEngine::is_a<SymEngine::Integer>(*p.second)) throw;

				auto exp = down_cast<const SymEngine::Integer &>(*p.second).as_int();
				if (exp < 0) {
					EmitPow(-exp);
					EmitDiv();
				} else {
					EmitPow(exp);
					EmitMul();
				}
			} else {
				EmitMul();
			}
		}
	}
	void bvisit(const SymEngine::Pow &x) {
		if (SymEngine::is_a<SymEngine::Integer>(*x.get_exp())) {
			apply(x.get_base());
			auto exp = down_cast<const SymEngine::Integer &>(*x.get_exp()).as_uint();
			EmitPow(exp);
		} else if (eq(*(x.get_base()), *SymEngine::E)) {
			apply(x.get_exp());
			EmitCall1((void *)(std::complex<double>(*)(const std::complex<double>&))std::exp<double>);
		} else {
			apply(x.get_base());
			apply(x.get_exp());
			EmitCall2((void *)(std::complex<double>(*)(const std::complex<double>&, const std::complex<double>&))std::pow<double>);
		}
	}
	void bvisit(const SymEngine::Log &x) {
		apply(x.get_arg());
		EmitCall1((void *)(std::complex<double>(*)(const std::complex<double>&))std::log<double>);
	}
	void bvisit(const SymEngine::Sin &x) {
		apply(x.get_arg());
		EmitCall1((void *)(std::complex<double>(*)(const std::complex<double>&))std::sin<double>);
	}
	void bvisit(const SymEngine::Cos &x) {
		apply(x.get_arg());
		EmitCall1((void *)(std::complex<double>(*)(const std::complex<double>&))std::cos<double>);
	}
	void bvisit(const SymEngine::Tan &x) {
		apply(x.get_arg());
		EmitCall1((void *)(std::complex<double>(*)(const std::complex<double>&))std::tan<double>);
	}

	void apply(const SymEngine::Basic &b) {
		apply(b.rcp_from_this());
	}
	void apply(const SymEngine::RCP<const SymEngine::Basic> &b) {
		b->accept(*this);
	}
};

class CompileVisitor : public SymEngine::BaseVisitor<CompileVisitor> {
public:
	CompiledFormula Compiled;

	void bvisit(const SymEngine::Basic &x) {
		throw FormulaCompileError(std::string("Unsupported operation: ") + SymEngine::type_code_name(x.get_type_code()) + "\nIn subexpression: " + SymEngine::str(x));
	}

	void bvisit(const SymEngine::FunctionSymbol &x) {
		throw FormulaCompileError(std::string("Unsupported function: ") + x.get_name());
	}

	void bvisit(const SymEngine::Number &x) {
		Compiled.AddOperation(Operation::LoadConst, Compiled.AddConstant(x));
	}
	void bvisit(const SymEngine::Constant &x) {
		Compiled.AddOperation(Operation::LoadConst, Compiled.AddConstant(x));
	}
	void bvisit(const SymEngine::Symbol &x) {
		const auto name = x.get_name();
		if (name == "Z" || name == "z") {
			Compiled.AddOperation(Operation::LoadZ);
		} else if (name == "C" || name == "c") {
			Compiled.AddOperation(Operation::LoadC);
		} else throw FormulaCompileError("Invalid symbol");
	}
	void bvisit(const SymEngine::Add &x) {
		const auto args = x.get_args();

		apply(args.front());

		for (const auto &p : args | std::views::drop(1)) {
			apply(p);
			Compiled.AddOperation(Operation::Add);
		}
	}
	void bvisit(const SymEngine::Mul &x) {
		apply(x.get_coef());

		for (const auto &p : x.get_dict()) {
			apply(p.first);

			if (neq(*p.second, *SymEngine::one)) {
				if (!SymEngine::is_a<SymEngine::Integer>(*p.second)) throw;

				auto exp = down_cast<const SymEngine::Integer &>(*p.second).as_int();
				if (exp < 0) {
					Compiled.AddOperation(Operation::PowUint, -exp);
					Compiled.AddOperation(Operation::Div);
				} else {
					Compiled.AddOperation(Operation::PowUint, exp);
					Compiled.AddOperation(Operation::Mul);
				}
			} else {
				Compiled.AddOperation(Operation::Mul);
			}
		}
	}
	void bvisit(const SymEngine::Pow &x) {
		if (SymEngine::is_a<SymEngine::Integer>(*x.get_exp())) {
			apply(x.get_base());
			auto exp = down_cast<const SymEngine::Integer &>(*x.get_exp()).as_uint();
			Compiled.AddOperation(Operation::PowUint, exp);
		} else if (eq(*(x.get_base()), *SymEngine::E)) {
			apply(x.get_exp());
			Compiled.AddOperation(Operation::Exp);
		} else {
			apply(x.get_base());
			apply(x.get_exp());
			Compiled.AddOperation(Operation::Pow);
		}
	}
	void bvisit(const SymEngine::Log &x) {
		apply(x.get_arg());
		Compiled.AddOperation(Operation::Log);
	}
	void bvisit(const SymEngine::Sin &x) {
		apply(x.get_arg());
		Compiled.AddOperation(Operation::Sin);
	}
	void bvisit(const SymEngine::Cos &x) {
		apply(x.get_arg());
		Compiled.AddOperation(Operation::Cos);
	}
	void bvisit(const SymEngine::Tan &x) {
		apply(x.get_arg());
		Compiled.AddOperation(Operation::Tan);
	}

	void apply(const SymEngine::Basic &b) {
		apply(b.rcp_from_this());
	}
	void apply(const SymEngine::RCP<const SymEngine::Basic> &b) {
		b->accept(*this);
	}
};

void call1(Encoder &encoder, void *func, RegisterXMM dstsrc) {
	using enum Register64;
	encoder.Sub(rsp, 32 + 2 * sizeof(__m128));
	encoder.Lea(rcx, rsp + 32);
	encoder.Lea(rdx, rsp + (32 + sizeof(__m128)));
	encoder.MovPD(*rdx, dstsrc);

	encoder.Call(func);

	encoder.MovPD(dstsrc, *rax);
	encoder.Add(rsp, 32 + 2 * sizeof(__m128));
}

void call2(Encoder &encoder, void *func, RegisterXMM dstsrc1, RegisterXMM src2) {
	using enum Register64;
	encoder.Sub(rsp, 32 + 3 * sizeof(__m128));
	encoder.Lea(rcx, rsp + 32);
	encoder.Lea(rdx, rsp + (32 + sizeof(__m128)));
	encoder.Lea(r8, rsp + (32 + sizeof(__m128) * 2));
	encoder.MovPD(*rdx, dstsrc1);
	encoder.MovPD(*r8, src2);

	encoder.Call(func);

	encoder.MovPD(dstsrc1, *rax);
	encoder.Add(rsp, 32 + 3 * sizeof(__m128));
}

void EncodeDouble(XMMComplexEncoder &encoder, const CompiledFormula &Formula) {
	using enum RegisterXMM;
	std::vector<MemoryAddress<>> ConstantAddresses;
	ConstantAddresses.reserve(Formula.Constants.size());
	for (auto &constant : Formula.Constants) {
		const auto value = SymEngine::eval_complex_double(*constant.Value);
		ConstantAddresses.push_back(encoder.Const(_mm_set_pd(value.imag(), value.real())));
	}
	uint8_t stack = 6; // start at xmm6, the first non volatile xmm register
	auto Top = [&]() { return RegisterXMM(stack - 1); };
	auto Push = [&]() { return RegisterXMM(stack++); };

	for (auto &operation : Formula.Operations) {
		switch (operation.operation) {
			case Operation::Nop: break;
			case Operation::LoadConst: {
				encoder.MovPD(Push(), *ConstantAddresses[operation.parameter]);
				break;
			}
			case Operation::LoadZ: {
				encoder.MovPD(Push(), *Register64::rdi);
				break;
			}
			case Operation::LoadC: {
				encoder.MovPD(Push(), *Register64::rsi);
				break;
			}
			case Operation::Add: {
				stack--;
				encoder.AddPD(RegisterXMM(stack - 1), RegisterXMM(stack));
				break;
			}
			case Operation::Sub: {
				stack--;
				encoder.SubPD(RegisterXMM(stack - 1), RegisterXMM(stack));
				break;
			}
			case Operation::Mul: {
				stack--;
				encoder.Cmul(RegisterXMM(stack - 1), RegisterXMM(stack), xmm0);
				break;
			}
			case Operation::Div: {
				stack--;
				encoder.Cdiv(RegisterXMM(stack - 1), RegisterXMM(stack), xmm0);
				break;
			}
			case Operation::Exp: {
				std::complex<double>(*func)(const std::complex<double>&) = std::exp<double>;
				call1(encoder, (void *)func, Top());
				break;
			}
			case Operation::PowUint: {
				encoder.Cpowui(operation.parameter, Top(), xmm0, xmm1, xmm2);
				break;
			}
			case Operation::Log: {
				std::complex<double>(*func)(const std::complex<double>&) = std::log<double>;
				call1(encoder, (void *)func, Top());
				break;
			}
			case Operation::Pow: {
				std::complex<double>(*func)(const std::complex<double>&, const std::complex<double>&) = std::pow<double>;
				stack--;
				call2(encoder, (void *)func, RegisterXMM(stack - 1), RegisterXMM(stack));
				break;
			}
			case Operation::Sin: {
				std::complex<double>(*func)(const std::complex<double>&) = std::sin<double>;
				call1(encoder, (void *)func, Top());
				break;
			}
			case Operation::Cos: {
				std::complex<double>(*func)(const std::complex<double>&) = std::cos<double>;
				call1(encoder, (void *)func, Top());
				break;
			}
			case Operation::Tan: {
				std::complex<double>(*func)(const std::complex<double>&) = std::tan<double>;
				call1(encoder, (void *)func, Top());
				break;
			}
			default: throw;
		}
	}
	if (stack != 7) throw;
}

void EncodeMpc(Encoder &encoder, const CompiledFormula &Formula) {
	using enum Register32;
	using enum Register64;

	for (auto &operation : Formula.Operations) {
		switch (operation.operation) {
			case Operation::Nop: break;
			case Operation::LoadConst: {
				throw;
				break;
			}
			case Operation::LoadZ: {
				encoder.Mov(rcx, r13);
				encoder.Mov(rdx, rdi);
				encoder.Mov(r8d, MPC_RNDNN);
				encoder.Call(mpc_set);
				encoder.Add(r13, sizeof(mpc_t));
				break;
			}
			case Operation::LoadC: {
				encoder.Mov(rcx, r13);
				encoder.Mov(rdx, rsi);
				encoder.Mov(r8d, MPC_RNDNN);
				encoder.Call(mpc_set);
				encoder.Add(r13, sizeof(mpc_t));
				break;
			}
			case Operation::Add: {
				encoder.Sub(r13, sizeof(mpc_t));
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8, r13);
				encoder.Mov(r9d, MPC_RNDNN);
				
				encoder.Call(mpc_add);
				break;
			}
			case Operation::Sub: {
				encoder.Sub(r13, sizeof(mpc_t));
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8, r13);
				encoder.Mov(r9d, MPC_RNDNN);

				encoder.Call(mpc_sub);
				break;
			}
			case Operation::Mul: {
				encoder.Sub(r13, sizeof(mpc_t));
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8, r13);
				encoder.Mov(r9d, MPC_RNDNN);

				encoder.Call(mpc_mul);
				break;
			}
			case Operation::Div: {
				encoder.Sub(r13, sizeof(mpc_t));
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8, r13);
				encoder.Mov(r9d, MPC_RNDNN);

				encoder.Call(mpc_div);
				break;
			}
			case Operation::Exp: {
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8d, MPC_RNDNN);

				encoder.Call(mpc_exp);
				break;
			}
			case Operation::PowUint: {
				encoder.Lea(rcx, r13 - sizeof(mpc_t));
				encoder.Mov(rdx, rcx);
				encoder.Mov(r8d, operation.parameter);
				encoder.Mov(r9d, MPC_RNDNN);

				encoder.Call(mpc_pow_ui);
				break;
			}
			default: throw;
		}
	}
}

void EncodePerturbation(XMMComplexEncoder &encoder, Encoder &RefEncoder, const CompiledFormula &Formula, size_t &ReferenceElementSize, bool NonZeroCriticalPoint) {
	using enum Register32;
	using enum Register64;
	using enum RegisterXMM;
	//using namespace RegisterAsMemAddress;

	std::vector<MemoryAddress<>> ConstantAddresses;
	ConstantAddresses.reserve(Formula.Constants.size());
	for (auto &constant : Formula.Constants) {
		const auto value = SymEngine::eval_complex_double(*constant.Value);
		ConstantAddresses.push_back(encoder.Const(_mm_set_pd(value.imag(), value.real())));
	}

	uint8_t stack = 0;
	auto Top = [&]() { return RegisterXMM(6 + stack - 1); }; // start at xmm6, the first non volatile xmm register
	auto Pop = [&]() { return RegisterXMM(6 + --stack); };
	auto Push = [&]() { return RegisterXMM(6 + stack++); };

	std::stack<Int> RefValueIndexes;
	ReferenceElementSize = NonZeroCriticalPoint ? 2 : 1;

	const Register64 MpcStackPointer = r13, ReferenceElementPointer = r12;

	auto SaveRefValue = [&]() {
		RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t) + offsetof(__mpc_struct, re));
		RefEncoder.Mov(edx, MPC_RNDNN);
		RefEncoder.Call(mpfr_get_d);
		RefEncoder.MovSD(*(ReferenceElementPointer + ReferenceElementSize * 16), xmm0);

		RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t) + offsetof(__mpc_struct, im));
		RefEncoder.Mov(edx, MPC_RNDNN);
		RefEncoder.Call(mpfr_get_d);
		RefEncoder.MovSD(*(ReferenceElementPointer + (ReferenceElementSize * sizeof(__m128) + sizeof(double))), xmm0);

		RefValueIndexes.push(ReferenceElementSize);
		ReferenceElementSize++;
	};

	auto LoadReferenceValue = [&](RegisterXMM dst, Int index) {
		if (index >= 0) {
			encoder.MovPD(dst, *(ReferenceElementPointer + index * sizeof(__m128)));
		} else if (index == -1) {
			encoder.MovPD(dst, *rsi);
		} else {
			encoder.MovPD(dst, *ConstantAddresses[-2 - index]);
		}
	};

	auto PopReferenceValue = [&](RegisterXMM dst) {
		LoadReferenceValue(dst, RefValueIndexes.top());
		RefValueIndexes.pop();
	};

	for (auto &operation : Formula.Operations) {
		switch (operation.operation) {
			case Operation::Nop: break;
			case Operation::LoadConst: {
				const auto &Const = *Formula.Constants[operation.parameter].Value;
				if (SymEngine::is_a<SymEngine::Integer>(Const)) {
					auto ConstInt = down_cast<const SymEngine::Integer &>(Const).as_int();

					RefEncoder.Mov(rcx, MpcStackPointer);
					RefEncoder.Mov(edx, ConstInt);
					RefEncoder.Mov(r8d, MPC_RNDNN);
					RefEncoder.Call(mpc_set_si);
					RefEncoder.Add(MpcStackPointer, sizeof(mpc_t));
				} else if (SymEngine::is_a<SymEngine::Rational>(Const)) {
					auto &ConstRational = down_cast<const SymEngine::Rational &>(Const);

					RefEncoder.Mov(rcx, MpcStackPointer);
					RefEncoder.Mov(edx, ConstRational.get_num()->as_int());
					RefEncoder.Mov(r8d, MPC_RNDNN);
					RefEncoder.Call(mpc_set_si);

					RefEncoder.Mov(rcx, MpcStackPointer);
					RefEncoder.Mov(rdx, rcx);
					RefEncoder.Mov(r8d, ConstRational.get_den()->as_uint());
					RefEncoder.Mov(r9d, MPC_RNDNN);
					RefEncoder.Call(mpc_div_ui);

					RefEncoder.Add(MpcStackPointer, sizeof(mpc_t));
				} else {
					RefEncoder.Mov(rcx, MpcStackPointer);
					const auto value = SymEngine::eval_complex_double(Const);
					RefEncoder.MovSD(xmm1, *RefEncoder.Const(value.real()));
					RefEncoder.MovSD(xmm2, *RefEncoder.Const(value.imag()));
					RefEncoder.Mov(r9d, MPC_RNDNN);
					RefEncoder.Call(mpc_set_d_d);
					RefEncoder.Add(MpcStackPointer, sizeof(mpc_t));
					throw;
				}
				Push();
				encoder.XorPD(Top(), Top());

				RefValueIndexes.push(Int(-2) - operation.parameter);

				break;
			}
			case Operation::LoadZ: {
				RefEncoder.Mov(rcx, MpcStackPointer);
				RefEncoder.Mov(rdx, rdi);
				RefEncoder.Mov(r8d, MPC_RNDNN);
				RefEncoder.Call(mpc_set);
				RefEncoder.Add(MpcStackPointer, sizeof(mpc_t));

				RefValueIndexes.push(0);

				encoder.MovPD(Push(), *rdi);
				break;
			}
			case Operation::LoadC: {
				RefEncoder.Mov(rcx, MpcStackPointer);
				RefEncoder.Mov(rdx, rsi);
				RefEncoder.Mov(r8d, MPC_RNDNN);
				RefEncoder.Call(mpc_set);
				RefEncoder.Add(MpcStackPointer, sizeof(mpc_t));

				RefValueIndexes.push(-1);

				encoder.MovPD(Push(), *(rdi + sizeof(__m128d)));
				break;
			}
			case Operation::Add: {
				RefEncoder.Sub(MpcStackPointer, sizeof(mpc_t));
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8, MpcStackPointer);
				RefEncoder.Mov(r9d, MPC_RNDNN);
				
				RefEncoder.Call(mpc_add);

				stack--;
				RefValueIndexes.pop();
				RefValueIndexes.pop();
				encoder.AddPD(RegisterXMM(6 + stack - 1), RegisterXMM(6 + stack));

				SaveRefValue();

				break;
			}
			case Operation::Sub: {
				RefEncoder.Sub(MpcStackPointer, sizeof(mpc_t));
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8, MpcStackPointer);
				RefEncoder.Mov(r9d, MPC_RNDNN);

				RefEncoder.Call(mpc_sub);

				stack--;
				RefValueIndexes.pop();
				RefValueIndexes.pop();
				encoder.SubPD(RegisterXMM(6 + stack - 1), RegisterXMM(6 + stack));

				SaveRefValue();

				break;
			}
			case Operation::Mul: {
				RefEncoder.Sub(MpcStackPointer, sizeof(mpc_t));
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8, MpcStackPointer);
				RefEncoder.Mov(r9d, MPC_RNDNN);

				RefEncoder.Call(mpc_mul);

				const RegisterXMM A = xmm0, B = xmm1;
				const RegisterXMM a = Pop();
				const RegisterXMM b = Top();

				PopReferenceValue(A);
				PopReferenceValue(B);
				encoder.AddPD(A, a);
				encoder.Cmul(b, A, xmm2);
				encoder.Cmul(B, a, xmm2);
				encoder.AddPD(b, B);

				SaveRefValue();
				break;
			}
			case Operation::Div: {
				RefEncoder.Sub(MpcStackPointer, sizeof(mpc_t));
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8, MpcStackPointer);
				RefEncoder.Mov(r9d, MPC_RNDNN);

				RefEncoder.Call(mpc_div);

				const RegisterXMM A = xmm0, B = xmm1;
				const RegisterXMM b = Pop();
				const RegisterXMM a = Top();

				PopReferenceValue(B);
				PopReferenceValue(A);
				encoder.MovPD(xmm2, B);
				encoder.Cmul(a, xmm2, xmm3);

				encoder.MovPD(xmm2, B);
				encoder.AddPD(xmm2, b);
				encoder.Cmul(B, xmm2, xmm3);

				encoder.Cmul(A, b, xmm2);

				encoder.SubPD(a, A);
				encoder.Cdiv(a, B, xmm2);

				SaveRefValue();

				break;
			}
			case Operation::Exp: {
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8d, MPC_RNDNN);

				RefEncoder.Call(mpc_exp);
				throw;
				break;
			}
			case Operation::PowUint: {
				if (operation.parameter == 1) break;
				RefEncoder.Lea(rcx, MpcStackPointer - sizeof(mpc_t));
				RefEncoder.Mov(rdx, rcx);
				RefEncoder.Mov(r8d, operation.parameter);
				RefEncoder.Mov(r9d, MPC_RNDNN);

				RefEncoder.Call(mpc_pow_ui);

				Int RefValueIndex = RefValueIndexes.top();
				RefValueIndexes.pop();

				uint32_t Power = operation.parameter;
				if (Power == 2) {
					LoadReferenceValue(xmm0, RefValueIndex);
					encoder.AddPD(xmm0, xmm0);
					encoder.AddPD(xmm0, Top());
					encoder.Cmul(Top(), xmm0, xmm1);
				} else if (Power == 3) {
					encoder.MovPD(xmm0, Top());
					LoadReferenceValue(xmm1, RefValueIndex);
					encoder.MovPD(xmm2, xmm1);
					encoder.MulPD(xmm1, *encoder.Const(_mm_set1_pd(3.0)));
					encoder.AddPD(xmm0, xmm1);
					encoder.Cmul(xmm1, xmm2, xmm3);
					encoder.MovPD(xmm2, Top());
					encoder.Cmul(xmm0, xmm2, xmm3);
					encoder.AddPD(xmm0, xmm1);
					encoder.Cmul(Top(), xmm0, xmm1);
				} else if (Power == 4) {
					LoadReferenceValue(xmm0, RefValueIndex);
					encoder.MovPD(xmm1, xmm0);
					encoder.AddPD(xmm0, xmm0);
					encoder.MovPD(xmm2, xmm0);
					encoder.Cmul(xmm1, xmm2, xmm3);
					encoder.AddPD(xmm0, Top());
					encoder.Cmul(Top(), xmm0, xmm2);
					encoder.AddPD(xmm1, Top());
					encoder.Cmul(Top(), xmm1, xmm2);
				} else {
					LoadReferenceValue(xmm0, RefValueIndex);
					while ((Power & 1) == 0) {
						encoder.MovPD(xmm1, xmm0);
						encoder.AddPD(xmm1, xmm1);
						encoder.AddPD(xmm1, Top());
						encoder.Cmul(Top(), xmm1, xmm2);
						Power >>= 1;
						if (Power > 1) encoder.Csqr(xmm0, xmm1, xmm2);
					}
					if (Power == 1) {
						SaveRefValue();
						break;
					}
					const RegisterXMM A = xmm0, B = xmm1, a = xmm2, b = Top();

					encoder.MovPD(a, Top());
					encoder.MovPD(B, xmm0);

					encoder.MovPD(xmm3, A); // A = A^2
					encoder.AddPD(xmm3, xmm3);
					encoder.AddPD(xmm3, a);
					encoder.Cmul(a, xmm3, xmm4);

					encoder.Csqr(A, xmm3, xmm4);
					Power >>= 1;

					while (Power > 1) {
						if ((Power & 1) != 0) { // B *= A
							encoder.MovPD(xmm3, A);
							encoder.AddPD(xmm3, a);
							encoder.Cmul(b, xmm3, xmm4);
							encoder.MovPD(xmm3, B);
							encoder.MovPD(xmm4, a);
							encoder.Cmul(xmm3, xmm4, xmm5);
							encoder.AddPD(b, xmm3);

							encoder.MovPD(xmm3, A);
							encoder.Cmul(B, xmm3, xmm4);
						}
						encoder.MovPD(xmm3, A); // A = A^2
						encoder.AddPD(xmm3, xmm3);
						encoder.AddPD(xmm3, a);
						encoder.Cmul(a, xmm3, xmm4);

						encoder.Csqr(A, xmm3, xmm4);
						Power >>= 1;
					}

					encoder.AddPD(A, a);
					encoder.Cmul(b, A, xmm3);
					encoder.Cmul(B, a, xmm3);
					encoder.AddPD(b, B);
				}
				SaveRefValue();

				break;
			}
			default: throw;
		}
	}
	ReferenceElementSize--;
}

void CompileDoublePrecisionFunction(XMMComplexEncoder &encoder) {
	using enum Register64;
	using enum RegisterXMM;
	using namespace RegisterAsMemAddress;

	encoder.Push(rbx);
	encoder.Push(rsi);
	encoder.Push(rdi);

	encoder.Push(r12);

	encoder.Push(rbp);
	encoder.Mov(rbp, rsp);
	encoder.And(rsp, ~(16ull - 1));
	encoder.Sub(rsp, 10 * 16);

	encoder.MovPD(Rsp[0 * sizeof(__m128)], xmm6);
	encoder.MovPD(Rsp[1 * sizeof(__m128)], xmm7);
	encoder.MovPD(Rsp[2 * sizeof(__m128)], xmm8);
	encoder.MovPD(Rsp[3 * sizeof(__m128)], xmm9);
	encoder.MovPD(Rsp[4 * sizeof(__m128)], xmm10);
	encoder.MovPD(Rsp[5 * sizeof(__m128)], xmm11);
	encoder.MovPD(Rsp[6 * sizeof(__m128)], xmm12);
	encoder.MovPD(Rsp[7 * sizeof(__m128)], xmm13);
	encoder.MovPD(Rsp[8 * sizeof(__m128)], xmm14);
	encoder.MovPD(Rsp[9 * sizeof(__m128)], xmm15);

	encoder.Mov(rdi, rcx);
	encoder.Mov(rsi, rdx);
	encoder.Mov(r12, r8);
	encoder.Xor(rbx, rbx);

	bool ConvergentBailout = false, DivergentBailout = true;

	encoder.MovPD(xmm6, *rsi);
	encoder.MovPD(*rdi, xmm6);

	if (ConvergentBailout) {
		encoder.MovPD(xmm15, *rdi);
	}

	size_t LoopPtr = encoder.GetPointer();

	try {
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression(Global::CustomFormula));
		EncodeDouble(encoder, visitor.Compiled);
	} catch (const std::exception &exception) {
		ErrorMessage(exception.what());
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodeDouble(encoder, visitor.Compiled);
	} catch (...) {
		ErrorMessage("Unknown exception");
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodeDouble(encoder, visitor.Compiled);
	}

	encoder.MovPD(*rdi, xmm6);

	size_t BreakConvergentPtr;
	if (ConvergentBailout) {
		encoder.MovPD(xmm0, xmm6);
		encoder.SubPD(xmm0, xmm15);
		encoder.DPPD(xmm0, xmm0, 0xFF);
		encoder.ComISD(xmm0, *encoder.Const(0x1p-24));
		BreakConvergentPtr = encoder.GetPointer(); // TODO: Improve this
		encoder.Jb(BreakConvergentPtr); // Placeholder

		encoder.MovPD(xmm15, xmm6);
	}

	size_t BreakDivergentPtr;
	if (DivergentBailout) {
		encoder.MovPD(xmm0, xmm6);
		encoder.DPPD(xmm0, xmm0, 0xFF);
		encoder.ComISD(xmm0, *encoder.Const(0x1p256));
		BreakDivergentPtr = encoder.GetPointer(); // TODO: Improve this
		encoder.Ja(BreakDivergentPtr); // Placeholder
	}

	encoder.Inc(rbx);
	encoder.Cmp(rbx, r12);
	encoder.Jb(LoopPtr);
	size_t LoopEndPtr = encoder.GetPointer();

	if (ConvergentBailout) {
		encoder.SetPointer(BreakConvergentPtr);
		encoder.Jb(LoopEndPtr);
	}
	if (DivergentBailout) {
		encoder.SetPointer(BreakDivergentPtr);
		encoder.Ja(LoopEndPtr);
	}
	encoder.SetPointer(LoopEndPtr);

	encoder.Mov(rax, rbx);

	encoder.MovPD(xmm6, Rsp[0 * sizeof(__m128)]);
	encoder.MovPD(xmm7, Rsp[1 * sizeof(__m128)]);
	encoder.MovPD(xmm8, Rsp[2 * sizeof(__m128)]);
	encoder.MovPD(xmm9, Rsp[3 * sizeof(__m128)]);
	encoder.MovPD(xmm10, Rsp[4 * sizeof(__m128)]);
	encoder.MovPD(xmm11, Rsp[5 * sizeof(__m128)]);
	encoder.MovPD(xmm12, Rsp[6 * sizeof(__m128)]);
	encoder.MovPD(xmm13, Rsp[7 * sizeof(__m128)]);
	encoder.MovPD(xmm14, Rsp[8 * sizeof(__m128)]);
	encoder.MovPD(xmm15, Rsp[9 * sizeof(__m128)]);

	encoder.Mov(rsp, rbp);
	encoder.Pop(rbp);

	encoder.Pop(r12);

	encoder.Pop(rdi);
	encoder.Pop(rsi);
	encoder.Pop(rbx);

	encoder.Ret();
}

void CompileHighPrecisionFunction(Encoder &encoder) {
	using enum Register32;
	using enum Register64;

	encoder.Push(rsi);
	encoder.Push(rdi);
	//encoder.Push(r12);
	encoder.Push(r13);

	encoder.Sub(rsp, 32);

	encoder.Mov(rdi, rcx);
	encoder.Mov(rsi, rdx);
	encoder.Mov(r13, r8);

	try {
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression(Global::CustomFormula));
		EncodeMpc(encoder, visitor.Compiled);
	} catch (const std::exception &exception) {
		ErrorMessage(exception.what());
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodeMpc(encoder, visitor.Compiled);
	} catch (...) {
		ErrorMessage("Unknown exception");
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodeMpc(encoder, visitor.Compiled);
	}

	encoder.Mov(rcx, rdi);
	encoder.Lea(rdx, r13 - sizeof(mpc_t));
	encoder.Mov(r8d, MPC_RNDNN);
	encoder.Call(mpc_set);

	encoder.Add(rsp, 32);

	encoder.Pop(r13);
	encoder.Pop(rdi);
	encoder.Pop(rsi);

	encoder.Ret();
}

void CompilePerturbationFunction(const SymEngine::Expression &Formula, Encoder &RefEncoder, XMMComplexEncoder &encoder, size_t &ReferenceElementSize, bool NonZeroCriticalPoint) {
	using enum Register32;
	using enum Register64;
	using enum RegisterXMM;
	using namespace RegisterAsMemAddress;

	RefEncoder.Push(rsi);
	RefEncoder.Push(rdi);
	RefEncoder.Push(r12);
	RefEncoder.Push(r13);

	RefEncoder.Sub(rsp, 32 + 8);

	RefEncoder.Mov(rdi, rcx);
	RefEncoder.Mov(rsi, rdx);
	RefEncoder.Mov(r13, r8);
	RefEncoder.Mov(r12, r9);


	encoder.Push(rbx);
	encoder.Push(rsi);
	encoder.Push(rdi);

	encoder.Push(r12);
	encoder.Push(r13);

	encoder.Push(rbp);
	encoder.Mov(rbp, rsp);
	encoder.Sub(rsp, 10 * 16);

	encoder.MovPD(Rsp[0 * sizeof(__m128)], xmm6 );
	encoder.MovPD(Rsp[1 * sizeof(__m128)], xmm7 );
	encoder.MovPD(Rsp[2 * sizeof(__m128)], xmm8 );
	encoder.MovPD(Rsp[3 * sizeof(__m128)], xmm9 );
	encoder.MovPD(Rsp[4 * sizeof(__m128)], xmm10);
	encoder.MovPD(Rsp[5 * sizeof(__m128)], xmm11);
	encoder.MovPD(Rsp[6 * sizeof(__m128)], xmm12);
	encoder.MovPD(Rsp[7 * sizeof(__m128)], xmm13);
	encoder.MovPD(Rsp[8 * sizeof(__m128)], xmm14);
	encoder.MovPD(Rsp[9 * sizeof(__m128)], xmm15);

	encoder.Mov(rdi, rcx);
	encoder.Mov(rsi, rdx);
	encoder.Mov(r12, r8);

	try {
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)Formula);

		EncodePerturbation(encoder, RefEncoder, visitor.Compiled, ReferenceElementSize, NonZeroCriticalPoint);
	} catch (const std::exception &exception) {
		ErrorMessage(exception.what());
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodePerturbation(encoder, RefEncoder, visitor.Compiled, ReferenceElementSize, false);
	} catch (...) {
		ErrorMessage("Unknown exception");
		CompileVisitor visitor;
		visitor.apply((const SymEngine::Basic &)SymEngine::Expression("z^2+c"));
		EncodePerturbation(encoder, RefEncoder, visitor.Compiled, ReferenceElementSize, false);
	}

	RefEncoder.Mov(rcx, rdi);
	RefEncoder.Lea(rdx, r13 - sizeof(mpc_t));
	RefEncoder.Mov(r8d, MPC_RNDNN);
	RefEncoder.Call(mpc_set);

	RefEncoder.Add(rsp, 32 + 8);

	RefEncoder.Pop(r13);
	RefEncoder.Pop(r12);
	RefEncoder.Pop(rdi);
	RefEncoder.Pop(rsi);

	RefEncoder.Ret();


	encoder.MovPD(*rdi, xmm6);

	encoder.MovPD(xmm6,  Rsp[0 * sizeof(__m128)]);
	encoder.MovPD(xmm7,  Rsp[1 * sizeof(__m128)]);
	encoder.MovPD(xmm8,  Rsp[2 * sizeof(__m128)]);
	encoder.MovPD(xmm9,  Rsp[3 * sizeof(__m128)]);
	encoder.MovPD(xmm10, Rsp[4 * sizeof(__m128)]);
	encoder.MovPD(xmm11, Rsp[5 * sizeof(__m128)]);
	encoder.MovPD(xmm12, Rsp[6 * sizeof(__m128)]);
	encoder.MovPD(xmm13, Rsp[7 * sizeof(__m128)]);
	encoder.MovPD(xmm14, Rsp[8 * sizeof(__m128)]);
	encoder.MovPD(xmm15, Rsp[9 * sizeof(__m128)]);

	encoder.Mov(rsp, rbp);
	encoder.Pop(rbp);

	encoder.Pop(r13);
	encoder.Pop(r12);

	encoder.Pop(rdi);
	encoder.Pop(rsi);
	encoder.Pop(rbx);

	encoder.Ret();
}

JITEvaluator::JITReference::~JITReference() {
	if (Ref) delete[]Ref;
}

JITEvaluator::JITEvaluator() {
	using enum Register16;
	using enum Register32;
	using enum Register64;
	using enum RegisterXMM;
	using namespace RegisterAsMemAddress;

	Buffer = (uint8_t *)VirtualAlloc(nullptr, 16384, MEM_COMMIT | MEM_RESERVE, PAGE_EXECUTE_READWRITE);
	if (!Buffer) throw;

	Encoder RefEncoder(Buffer, 8192, 2048);
	XMMComplexEncoder encoder(Buffer + 8192, 8192, 2048);

	size_t Begin = encoder.GetPointer();

	size_t BeginRef = RefEncoder.GetPointer();
	size_t BeginPt = encoder.GetPointer();

	SymEngine::Expression Formula = SymEngine::Expression(Global::CustomFormula);
	SymEngine::Expression Derivative = Formula.diff(SymEngine::symbol("z"));
	auto CriticalPoints = SymEngine::solve(Derivative, SymEngine::symbol("z"));
	if (SymEngine::is_a<SymEngine::FiniteSet>(*CriticalPoints)) {
		CriticalPoint = CriticalPoints->get_args().front();
		NonZeroCriticalPoint = CriticalPoint != 0;
	} else {
		CriticalPoint = 0;
		NonZeroCriticalPoint = false;
	}
	CompilePerturbationFunction(Formula, RefEncoder, encoder, ReferenceElementSize, NonZeroCriticalPoint);

	funcRef = FuncPointerRef(Buffer + BeginRef);
	funcPt = FuncPointerPt(Buffer + 8192 + BeginPt);
}

JITEvaluator::~JITEvaluator() {
	VirtualFree(Buffer, 0, MEM_RELEASE);
}

Reference *JITEvaluator::GenerateReference(const EvaluationParameters &parameters) {
	IntIter MaxIt = Global::ItLim;
	JITReference &Ref = *new JITReference;
	Ref.Refcr = convert<SRReal>(parameters.CenterCoordinate.X);
	Ref.Refci = convert<SRReal>(parameters.CenterCoordinate.Y);
	Ref.CenterCoordinate = parameters.CenterCoordinate;
	Ref.Ref = new __m128d[MaxIt * ReferenceElementSize + 2];

	mpc_t Z, C, CriticalPointMPC, temp[10];
	mp_bitcnt_t Precision = parameters.CenterCoordinate.X.get_prec();
	mpc_init2(Z, Precision);
	mpc_init2(C, Precision);
	mpc_init2(CriticalPointMPC, Precision);
	SymEngine::eval_mpc(CriticalPointMPC, CriticalPoint.subs({{ SymEngine::symbol("c"), SymEngine::complex_mpc(SymEngine::mpc_class(C)) }}), MPFR_RNDN);
	for (size_t i = 0; i < 10; i++) {
		mpc_init2(temp[i], Precision);
	}

	mpc_set_f_f(C, parameters.CenterCoordinate.X.get_mpf_t(), parameters.CenterCoordinate.Y.get_mpf_t(), MPC_RNDNN);
	mpc_set(Z, CriticalPointMPC, MPC_RNDNN);

	Ref.Ref[0] = _mm_set_pd(mpfr_get_d(Z->im, MPFR_RNDN), mpfr_get_d(Z->re, MPFR_RNDN));//_mm_setzero_pd();
	if (NonZeroCriticalPoint) {
		Ref.Ref[1] = _mm_setzero_pd();
	}
	size_t i = 0;
	while (i < MaxIt) {
		funcRef(Z, C, temp, Ref.Ref + i * ReferenceElementSize);

		i++;

		if (NonZeroCriticalPoint) {
			mpc_sub(temp[0], Z, CriticalPointMPC, MPC_RNDNN);
			double re = mpfr_get_d(temp[0]->re, MPFR_RNDN);
			double im = mpfr_get_d(temp[0]->im, MPFR_RNDN);
			Ref.Ref[i * ReferenceElementSize + 1] = _mm_set_pd(im, re);
		}
		__m128d z = Ref.Ref[i * ReferenceElementSize];
		double Magnitude = _mm_cvtsd_f64(_mm_dp_pd(z, z, 0xFF));
		if (Magnitude > 0x1p256) break;
	}
	Ref.RefIteration = i;

	mpc_clear(Z);
	mpc_clear(C);
	mpc_clear(CriticalPointMPC);
	for (size_t i = 0; i < 10; i++) {
		mpc_clear(temp[i]);
	}

	return &Ref;
}

void JITEvaluator::Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *reference) {
	JITReference &Ref = *static_cast<JITReference *>(reference);
	IntIter MaxIt = Global::ItLim;
	Global::MaxIt = MaxIt;

	RasterizingInterface &RI = rasterizer.GetInterface();

	HRReal dcr, dci;
	while (RI.GetCoordinate(dcr, dci)) {
		size_t i = 0, j = 0;
		__m128d Variables[2] = { _mm_set_pd(0.0, 0.0), _mm_set_pd(double(dci), double(dcr)) }; // dz, dc
		__m128d Constants[1] = { _mm_set_pd(double(Ref.Refci), double(Ref.Refcr)) }; // C
		__m128d PrevZ = Variables[0], diff;

		while (i < MaxIt) {
			funcPt(Variables, Constants, Ref.Ref[j * ReferenceElementSize]);
			i++;
			j++;

			__m128d z = _mm_add_pd(Ref.Ref[j * ReferenceElementSize], Variables[0]);
			__m128d zmcp;

			if (NonZeroCriticalPoint) {
				zmcp = _mm_add_pd(Ref.Ref[j * ReferenceElementSize + 1], Variables[0]);
			} else {
				zmcp = z;
			}

			double Magnitude = _mm_cvtsd_f64(_mm_dp_pd(z, z, 0xFF));
			if (Magnitude > 0x1p16) {
				break;
			}
			diff = _mm_sub_pd(z, PrevZ);
			Magnitude = _mm_cvtsd_f64(_mm_dp_pd(diff, diff, 0xFF));
			if (Magnitude < 0x1p-24) {
				break;
			}
			double MagnitudeDz = _mm_cvtsd_f64(_mm_dp_pd(Variables[0], Variables[0], 0xFF));
			double MagnitudeZmcp = _mm_cvtsd_f64(_mm_dp_pd(zmcp, zmcp, 0xFF));
			if (j >= Ref.RefIteration || MagnitudeZmcp < MagnitudeDz) {
				Variables[0] = zmcp;
				j = 0;
			}
			PrevZ = z;
		}

		RI.WriteResults(i);
	}
	rasterizer.FreeInterface(RI);
}
#endif