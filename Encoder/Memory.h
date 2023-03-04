#pragma once
#include "Register.h"

inline bool IsPowerOfTwo(int64_t x) {
	return !(x & (x - 1));
}

uint8_t ScaleToS(uint8_t Scale);

template <Register32Or64 RegType>
struct MemoryOperand;

template <Register32Or64 RegType = Register64>
struct MemoryAddress {
	//using RegType = Register64;
	uint8_t S = 1;
	RegType I = RegType(), B = RegType();
	int32_t Disp = 0;
	bool HasI = false, HasB = false;

	constexpr MemoryAddress(int32_t Disp) : Disp(Disp) {}
	constexpr MemoryAddress(RegType reg) : B(reg), HasB(true) {}

	MemoryAddress &operator+=(const MemoryAddress &x);
	MemoryAddress &operator-=(int32_t x) {
		Disp -= x;
		return *this;
	}
	MemoryAddress &operator*=(int32_t x);
	inline MemoryOperand<RegType> operator[](const MemoryAddress &x) const;
	inline MemoryOperand<RegType> operator[](int32_t x) const;
};

namespace RegisterAsMemAddress {
	constexpr MemoryAddress Rax = Register64::rax;
	constexpr MemoryAddress Rcx = Register64::rcx;
	constexpr MemoryAddress Rdx = Register64::rdx;
	constexpr MemoryAddress Rbx = Register64::rbx;
	constexpr MemoryAddress Rsp = Register64::rsp;
	constexpr MemoryAddress Rbp = Register64::rbp;
	constexpr MemoryAddress Rsi = Register64::rsi;
	constexpr MemoryAddress Rdi = Register64::rdi;
	constexpr MemoryAddress R8  = Register64::r8;
	constexpr MemoryAddress R9  = Register64::r9;
	constexpr MemoryAddress R10 = Register64::r10;
	constexpr MemoryAddress R11 = Register64::r11;
	constexpr MemoryAddress R12 = Register64::r12;
	constexpr MemoryAddress R13 = Register64::r13;
	constexpr MemoryAddress R14 = Register64::r14;
	constexpr MemoryAddress R15 = Register64::r15;

	constexpr MemoryAddress Eax  = Register32::eax;
	constexpr MemoryAddress Ecx  = Register32::ecx;
	constexpr MemoryAddress Edx  = Register32::edx;
	constexpr MemoryAddress Ebx  = Register32::ebx;
	constexpr MemoryAddress Esp  = Register32::esp;
	constexpr MemoryAddress Ebp  = Register32::ebp;
	constexpr MemoryAddress Esi  = Register32::esi;
	constexpr MemoryAddress Edi  = Register32::edi;
	constexpr MemoryAddress R8d  = Register32::r8d;
	constexpr MemoryAddress R9d  = Register32::r9d;
	constexpr MemoryAddress R10d = Register32::r10d;
	constexpr MemoryAddress R11d = Register32::r11d;
	constexpr MemoryAddress R12d = Register32::r12d;
	constexpr MemoryAddress R13d = Register32::r13d;
	constexpr MemoryAddress R14d = Register32::r14d;
	constexpr MemoryAddress R15d = Register32::r15d;
}

template <Register32Or64 RegType> MemoryAddress<RegType> operator+(RegType x, const RegType &y) { return MemoryAddress(x) += y; }

template <Register32Or64 RegType> MemoryAddress<RegType> operator+(RegType x, int32_t y) { return MemoryAddress(x) += y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator-(RegType x, int32_t y) { return MemoryAddress(x) -= y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator*(RegType x, int32_t y) { return MemoryAddress(x) *= y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator+(int32_t x, RegType y) { return MemoryAddress(y) += x; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator*(int32_t x, RegType y) { return MemoryAddress(y) *= x; }

template <Register32Or64 RegType> MemoryAddress<RegType> operator+(MemoryAddress<RegType> x, const MemoryAddress<RegType> &y) { return x += y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator+(MemoryAddress<RegType> x, int32_t y) { return x += y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator-(MemoryAddress<RegType> x, int32_t y) { return x -= y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator*(MemoryAddress<RegType> x, int32_t y) { return x *= y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator+(int32_t x, MemoryAddress<RegType> y) { return y += x; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator*(int32_t x, MemoryAddress<RegType> y) { return y *= x; }

template <Register32Or64 RegType> MemoryAddress<RegType> operator+(MemoryAddress<RegType> x, RegType y) { return x += y; }
template <Register32Or64 RegType> MemoryAddress<RegType> operator+(RegType x, MemoryAddress<RegType> y) { return y += x; }

template <Register32Or64 RegType>
struct MemoryOperand {
	int32_t Disp = 0;
	uint8_t Mod, RM, S = 0b00, I = 0b100, B = 0b101;
	uint8_t DispSize = 0;
	bool HasSIB = false;

	void SetDisp(const MemoryAddress<RegType> &Address);
	void SetSIB(const MemoryAddress<RegType> &Address);
	void SetSI(const MemoryAddress<RegType> &Address);
	void SetB(const MemoryAddress<RegType> &Address);

	explicit MemoryOperand(const MemoryAddress<RegType> &Address);
};
