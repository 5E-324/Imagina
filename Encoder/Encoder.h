#pragma once
#include <stdint.h>
//#include <concepts>

#include "Register.h"
#include "Memory.h"

template <size_t> struct IntOfSizeStruct { using Type = void; };
template <> struct IntOfSizeStruct<8> { using Type = int8_t; };
template <> struct IntOfSizeStruct<16> { using Type = int16_t; };
template <> struct IntOfSizeStruct<32> { using Type = int32_t; };
template <> struct IntOfSizeStruct<64> { using Type = int64_t; };

template <size_t size> using IntOfSize = typename IntOfSizeStruct<size>::Type;

template <Register T> using ImmType = IntOfSize<RegisterSize<T>>;
template <Register T> using ImmMost32Type = IntOfSize<std::min(size_t(32), RegisterSize<T>)>;

template <Register32Or64 RegType> MemoryOperand<RegType> operator*(const MemoryAddress<RegType> &x) {
	return MemoryOperand<RegType>(x);
}

template <Register32Or64 RegType> MemoryOperand<RegType> operator*(const RegType &x) {
	return MemoryOperand<RegType>(x);
}

template <Register32Or64 RegType> MemoryOperand<RegType> MemoryAddress<RegType>::operator[](const MemoryAddress<RegType> &x) const {
	return *(*this + x);
}

template <Register32Or64 RegType> MemoryOperand<RegType> MemoryAddress<RegType>::operator[](int32_t x) const {
	return *(*this + x);
}

class Encoder {
	uint8_t *Buffer;
	size_t Size, Pointer = 0;
	size_t ConstSize = 0, ConstPointer = 0;

public:
	Encoder(uint8_t *Buffer, size_t Size) : Buffer(Buffer), Size(Size) {}
	Encoder(uint8_t *Buffer, size_t Size, size_t ConstSize) : Buffer(Buffer), Size(Size), Pointer(ConstSize), ConstSize(ConstSize) {}
	Encoder(const Encoder &) = delete;

	template<typename T>
	MemoryAddress<Register64> Const(const T &x) {
		ConstPointer = (((ConstPointer - 1) / alignof(T)) + 1) * alignof(T); // Round up to multiple of alignof(T)
		if (ConstPointer + sizeof(T) > ConstSize) throw;
		MemoryAddress Address = MemoryAddress(ConstPointer);
		(T &)Buffer[ConstPointer] = x;
		ConstPointer += sizeof(T);
		return Address;
	}

	size_t GetPointer() { return Pointer; }
	void SetPointer(size_t pointer) { Pointer = pointer; }
	MemoryAddress<Register64> GetAddress() { return MemoryAddress<Register64>(Pointer); }

	uint8_t HighBitIfIsReg(Register auto reg) { return uint8_t(reg) >> 3; }
	uint8_t HighBitIfIsReg(uint8_t x) { return x; }

	uint8_t HighBit(Register auto reg) { return uint8_t(reg) >> 3; }
	uint8_t Low3Bits(Register auto reg) { return uint8_t(reg) & 0b0111; }

	Encoder &Byte(uint8_t byte) {
		Buffer[Pointer++] = byte;
		return *this;
	}

	Encoder &Bytes(uint8_t byte) {
		Buffer[Pointer++] = byte;
		return *this;
	}
	template <typename... Args>
	Encoder &Bytes(uint8_t byte0, Args... bytes) {
		Buffer[Pointer++] = byte0;
		return Bytes(bytes...);
	}

	Encoder &Bytes(uint8_t byte0, uint8_t byte1, uint8_t byte2) {
		Buffer[Pointer++] = byte0;
		Buffer[Pointer++] = byte1;
		Buffer[Pointer++] = byte2;
		return *this;
	}

	Encoder &Imm8(int8_t  imm) { (int8_t &)Buffer[Pointer] = imm; Pointer += 1; return *this; }
	Encoder &Imm16(int16_t imm) { (int16_t &)Buffer[Pointer] = imm; Pointer += 2; return *this; }
	Encoder &Imm32(int32_t imm) { (int32_t &)Buffer[Pointer] = imm; Pointer += 4; return *this; }
	Encoder &Imm64(int64_t imm) { (int64_t &)Buffer[Pointer] = imm; Pointer += 8; return *this; }

	Encoder &Lock() { return Byte(0xF0); }
	Encoder &OperandSize() { return Byte(0x66); }
	Encoder &AddressSize() { return Byte(0x67); }
	Encoder &AddressSize(const MemoryOperand<Register32> &) { return Byte(0x67); }
	Encoder &AddressSize(const MemoryOperand<Register64> &) { return *this; }

	Encoder &Rex(uint8_t W, uint8_t R, uint8_t X, uint8_t B) {
		if (!W && !R && !X && !B) { // Not needed
			return *this;
		}
		return Byte(0b0'0100'0000 | (W << 3) | (R << 2) | (X << 1) | B);
	}

	Encoder &RexW(uint8_t	 R, uint8_t	   B) { return Rex(1, R, 0, B); }
	Encoder &RexW(Register64 R, uint8_t	   B) { return RexW(HighBit(R), B); }
	Encoder &RexW(uint8_t	 R, Register64 B) { return RexW(R, HighBit(B)); }
	Encoder &RexW(Register64 R, Register64 B) { return RexW(HighBit(R), HighBit(B)); }

	Encoder &Rex(RegisterOrInt8 auto R, RegisterOrInt8 auto B) { return Rex(0, HighBitIfIsReg(R), 0, HighBitIfIsReg(B)); }

	Encoder &Rex(RegisterOrInt8 auto R, RegisterOrInt8 auto X, RegisterOrInt8 auto B) {
		return Rex(0, HighBitIfIsReg(R), HighBitIfIsReg(X), HighBitIfIsReg(B));
	}

	Encoder &RexW(Register64OrInt8 auto R, Register64OrInt8 auto X, Register64OrInt8 auto B) {
		return Rex(1, HighBitIfIsReg(R), HighBitIfIsReg(X), HighBitIfIsReg(B));
	}

	template <Register32Or64 RegType>
	Encoder &Rex(RegisterOrInt8 auto R, const MemoryOperand<RegType> &Mem) {
		if (Mem.HasSIB) {
			return Rex(R, Mem.I >> 3, Mem.B >> 3);
		} else {
			return Rex(R, Mem.RM >> 3);
		}
	}

	template <Register32Or64 RegType>
	Encoder &RexW(Register64OrInt8 auto R, const MemoryOperand<RegType> &Mem) {
		if (Mem.HasSIB) {
			return RexW(R, Mem.I >> 3, Mem.B >> 3);
		} else {
			return RexW(R, Mem.RM >> 3);
		}
	}

	template <Register32Or64 RegType> Encoder &RexAndOpSize(Register16 Reg, const MemoryOperand<RegType> &Mem) { return OperandSize().Rex(Reg, Mem); }
	template <Register32Or64 RegType> Encoder &RexAndOpSize(Register32 Reg, const MemoryOperand<RegType> &Mem) { return Rex(Reg, Mem); }
	template <Register32Or64 RegType> Encoder &RexAndOpSize(Register64 Reg, const MemoryOperand<RegType> &Mem) { return RexW(Reg, Mem); }
	//template <Register32Or64 RegType> Encoder &RexAndOpSize(RegisterXMM Reg, const MemoryOperand<RegType> &Mem) { return Rex(Reg, Mem); }

	template <Register RegType, Register32Or64 AddrRegType> Encoder &RexAndAddrSize(RegType Reg, const MemoryOperand<AddrRegType> &Mem) {
		return AddressSize(Mem).Rex(Reg, Mem);
	}

	template <Register RegType, Register32Or64 AddrRegType> Encoder &RexAndSizes(RegType Reg, const MemoryOperand<AddrRegType> &Mem) {
		return AddressSize(Mem).RexAndOpSize(Reg, Mem);
	}

	Encoder &RexAndOpSize(uint8_t R, Register16 B) { return OperandSize().Rex(R, B); }
	Encoder &RexAndOpSize(uint8_t R, Register32 B) { return Rex(R, B); }
	Encoder &RexAndOpSize(uint8_t R, Register64 B) { return RexW(R, B); }

	Encoder &RexAndOpSize(Register16 R, Register16 B) { return OperandSize().Rex(R, B); }
	Encoder &RexAndOpSize(Register32 R, Register32 B) { return Rex(R, B); }
	Encoder &RexAndOpSize(Register64 R, Register64 B) { return RexW(R, B); }

	Encoder &ModRM(uint8_t Mod, uint8_t		  Reg, uint8_t		 RM) { return Byte((Mod << 6) | (Reg << 3) | RM); }
	Encoder &ModRM(uint8_t Mod, Register auto Reg, uint8_t		 RM) { return ModRM(Mod, Low3Bits(Reg), RM); }
	Encoder &ModRM(uint8_t Mod, uint8_t		  Reg, Register auto RM) { return ModRM(Mod, Reg, Low3Bits(RM)); }
	Encoder &ModRM(uint8_t Mod, Register auto Reg, Register auto RM) { return ModRM(Mod, Low3Bits(Reg), Low3Bits(RM)); }

	Encoder &SIB(uint8_t Scale, RegisterOrInt8 auto Index, RegisterOrInt8 auto Base) { return ModRM(Scale, Index, Base); } // SIB has the same layout as ModRM

	template <Register32Or64 RegType>
	Encoder &RegMem(RegisterOrInt8 auto Reg, const MemoryOperand<RegType> &Mem) {
		if (Mem.HasSIB) {
			ModRM(Mem.Mod, Reg, Mem.RM).SIB(Mem.S, Mem.I & 0b111, Mem.B & 0b111);
		} else {
			ModRM(Mem.Mod, Reg, Mem.RM & 0b111);
		}
		if (Mem.DispSize == 32) {
			if (Mem.Mod == 0b00 && Mem.RM == 0b101) {
				Imm32(Mem.Disp - (GetPointer() + 4));
			} else {
				Imm32(Mem.Disp);
			}
		} else if (Mem.DispSize == 8) {
			Imm8(uint8_t(Mem.Disp));
		}
		return *this;
	}

	template <Register T> void Imm(ImmType<T> imm) {
		(ImmType<T> &)Buffer[Pointer] = imm;
		Pointer += sizeof(ImmType<T>);
	}

	template <Register T> void ImmMost32(ImmMost32Type<T> imm) {
		(ImmMost32Type<T> &)Buffer[Pointer] = imm;
		Pointer += sizeof(ImmMost32Type<T>);
	}

	template <Register T> void Op(uint8_t Opcode, T dst, T src) {
		RexAndOpSize(dst, src).Byte(Opcode).ModRM(0b11, dst, src);
	}

	template <Register RegType, Register32Or64 AddrRegType>
	void Op(uint8_t Opcode, RegType dst, const MemoryOperand<AddrRegType> &src) {
		RexAndSizes(dst, src).Byte(Opcode).RegMem(dst, src);
	}

	template <Register RegType, Register32Or64 AddrRegType>
	void Op(uint8_t Opcode, const MemoryOperand<AddrRegType> &dst, RegType src) {
		RexAndSizes(src, dst).Byte(Opcode).RegMem(src, dst);
	}

	void OpImm8(uint8_t Opcode, uint8_t Ext, Register auto dst, int8_t imm) {
		RexAndOpSize(0, dst).Byte(Opcode).ModRM(0b11, Ext, dst).Imm8(imm);
	}

	template <Register T> void OpImm(uint8_t Opcode, uint8_t Ext, T dst, ImmMost32Type<T> imm) {
		RexAndOpSize(0, dst).Byte(Opcode).ModRM(0b11, Ext, dst).template ImmMost32<T>(imm);
	}

	template <Register T> void OpGroup1(uint8_t op, T dst, ImmMost32Type<T> imm) {
		if (imm == int8_t(imm)) {
			return OpImm8(0x83, op, dst, int8_t(imm));
		} else {
			return OpImm(0x81, op, dst, imm);
		}
	}

	template <Register32Or64 AddrRegType>
	void Call(MemoryOperand<AddrRegType> Mem) {
		RexW(0, 0).Byte(0xFF).RegMem(2, Mem);
	}

	void Call(void *func) {
		int64_t Disp = (uint8_t *)func - (Buffer + Pointer + 5);
		if (Disp == int32_t(Disp)) {
			Byte(0xE8).Imm32(Disp);
		} else {
			Call(*Const(func));
		}
	}

	void Push(Register64 reg) {
		Rex(0, reg).Byte(0x50 | Low3Bits(reg));
	}

	void Pop(Register64 reg) {
		Rex(0, reg).Byte(0x58 | Low3Bits(reg));
	}

	void Inc(Register auto Reg) {
		RexAndOpSize(0, Reg).Byte(0xFF).ModRM(0b11, 0, Reg);
	}

	void Dec(Register auto Reg) {
		RexAndOpSize(0, Reg).Byte(0xFF).ModRM(0b11, 1, Reg);
	}

	template <Register T> void Add(T dst, T src) { Op(0x03, dst, src); }
	template <Register T> void Or(T dst, T src) { Op(0x0B, dst, src); }
	template <Register T> void Adc(T dst, T src) { Op(0x13, dst, src); }
	template <Register T> void Sbb(T dst, T src) { Op(0x1B, dst, src); }
	template <Register T> void And(T dst, T src) { Op(0x23, dst, src); }
	template <Register T> void Sub(T dst, T src) { Op(0x2B, dst, src); }
	template <Register T> void Xor(T dst, T src) { Op(0x33, dst, src); }
	template <Register T> void Cmp(T dst, T src) { Op(0x3B, dst, src); }

	template <Register32Or64 AddrRegType> void Add(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x01, dst, src); }
	template <Register32Or64 AddrRegType> void Or(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x09, dst, src); }
	template <Register32Or64 AddrRegType> void Adc(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x11, dst, src); }
	template <Register32Or64 AddrRegType> void Sbb(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x19, dst, src); }
	template <Register32Or64 AddrRegType> void And(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x21, dst, src); }
	template <Register32Or64 AddrRegType> void Sub(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x29, dst, src); }
	template <Register32Or64 AddrRegType> void Xor(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x31, dst, src); }
	template <Register32Or64 AddrRegType> void Cmp(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x39, dst, src); }

	template <Register32Or64 AddrRegType> void Add(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x03, dst, src); }
	template <Register32Or64 AddrRegType> void Or(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x0B, dst, src); }
	template <Register32Or64 AddrRegType> void Adc(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x13, dst, src); }
	template <Register32Or64 AddrRegType> void Sbb(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x1B, dst, src); }
	template <Register32Or64 AddrRegType> void And(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x23, dst, src); }
	template <Register32Or64 AddrRegType> void Sub(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x2B, dst, src); }
	template <Register32Or64 AddrRegType> void Xor(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x33, dst, src); }
	template <Register32Or64 AddrRegType> void Cmp(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x3B, dst, src); }

	template <Register T> void Add(T dst, ImmMost32Type<T> imm) { OpGroup1(0x00, dst, imm); }
	template <Register T> void Or (T dst, ImmMost32Type<T> imm) { OpGroup1(0x01, dst, imm); }
	template <Register T> void Adc(T dst, ImmMost32Type<T> imm) { OpGroup1(0x02, dst, imm); }
	template <Register T> void Sbb(T dst, ImmMost32Type<T> imm) { OpGroup1(0x03, dst, imm); }
	template <Register T> void And(T dst, ImmMost32Type<T> imm) { OpGroup1(0x04, dst, imm); }
	template <Register T> void Sub(T dst, ImmMost32Type<T> imm) { OpGroup1(0x05, dst, imm); }
	template <Register T> void Xor(T dst, ImmMost32Type<T> imm) { OpGroup1(0x06, dst, imm); }
	template <Register T> void Cmp(T dst, ImmMost32Type<T> imm) { OpGroup1(0x07, dst, imm); }

	template <Register T> void Mov(T dst, T src) { Op(0x8B, dst, src); }

	template <Register T> void Mov(T dst, ImmType<T> imm) { RexAndOpSize(0, dst).Byte(0xB8 | uint8_t(dst)).template Imm<T>(imm); }

	template <> void Mov<Register64>(Register64 dst, int64_t imm) {
		if (uint64_t(imm) == uint32_t(imm)) {
			return Mov(Register32(dst), int32_t(imm));
		} else if (imm == int32_t(imm)) {
			return OpImm(0xC7, 0, dst, int32_t(imm));
		}
		RexAndOpSize(0, dst).Byte(0xB8 | uint8_t(dst)).Imm64(imm);
	}

	template <Register32Or64 AddrRegType> void Mov(const MemoryOperand<AddrRegType> &dst, Register auto src) { Op(0x89, dst, src); }
	template <Register32Or64 AddrRegType> void Mov(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x8B, dst, src); }

	template <Register32Or64 AddrRegType> void Lea(Register auto dst, const MemoryOperand<AddrRegType> &src) { Op(0x8D, dst, src); }
	template <Register32Or64 AddrRegType> void Lea(Register auto dst, const MemoryAddress<AddrRegType> &src) { Op(0x8D, dst, MemoryOperand(src)); }

	void MovSD(RegisterXMM dst, RegisterXMM src) { Byte(0xF2).Rex(dst, src).Bytes(0x0F, 0x10).ModRM(0b11, dst, src); }
	void AddSD(RegisterXMM dst, RegisterXMM src) { Byte(0xF2).Rex(dst, src).Bytes(0x0F, 0x58).ModRM(0b11, dst, src); }
	void SubSD(RegisterXMM dst, RegisterXMM src) { Byte(0xF2).Rex(dst, src).Bytes(0x0F, 0x5C).ModRM(0b11, dst, src); }
	void MulSD(RegisterXMM dst, RegisterXMM src) { Byte(0xF2).Rex(dst, src).Bytes(0x0F, 0x59).ModRM(0b11, dst, src); }
	void DivSD(RegisterXMM dst, RegisterXMM src) { Byte(0xF2).Rex(dst, src).Bytes(0x0F, 0x5E).ModRM(0b11, dst, src); }

	void ComISD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x2F).ModRM(0b11, dst, src); }

	void MovPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x10).ModRM(0b11, dst, src); }
	void AddPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x58).ModRM(0b11, dst, src); }
	void SubPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x5C).ModRM(0b11, dst, src); }
	void MulPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x59).ModRM(0b11, dst, src); }
	void DivPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x5E).ModRM(0b11, dst, src); }

	void AndPD (RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x54).ModRM(0b11, dst, src); }
	void AndNPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x55).ModRM(0b11, dst, src); }
	void OrPD  (RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x56).ModRM(0b11, dst, src); }
	void XorPD (RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x57).ModRM(0b11, dst, src); }

	void AddSubPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0xD0).ModRM(0b11, dst, src); }

	void UnpckLPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x14).ModRM(0b11, dst, src); }
	void UnpckHPD(RegisterXMM dst, RegisterXMM src) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x15).ModRM(0b11, dst, src); }

	void ShufPD(RegisterXMM dst, RegisterXMM src, uint8_t imm8) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0xC6).ModRM(0b11, dst, src).Imm8(imm8); }
	void DPPD(RegisterXMM dst, RegisterXMM src, uint8_t imm8) { Byte(0x66).Rex(dst, src).Bytes(0x0F, 0x3A, 0x41).ModRM(0b11, dst, src).Imm8(imm8); }

	void AddPS(RegisterXMM dst, RegisterXMM src) { Rex(dst, src).Bytes(0x0F, 0x58).ModRM(0b11, dst, src); }
	void AddSS(RegisterXMM dst, RegisterXMM src) { Byte(0xF3).Rex(dst, src).Bytes(0x0F, 0x58).ModRM(0b11, dst, src); }


	template <Register32Or64 AddrRegType> void MovSD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0xF2).RexAndAddrSize(dst, src).Bytes(0x0F, 0x10).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void AddSD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0xF2).RexAndAddrSize(dst, src).Bytes(0x0F, 0x58).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void SubSD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0xF2).RexAndAddrSize(dst, src).Bytes(0x0F, 0x5C).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void MulSD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0xF2).RexAndAddrSize(dst, src).Bytes(0x0F, 0x59).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void DivSD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0xF2).RexAndAddrSize(dst, src).Bytes(0x0F, 0x5E).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void ComISD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x2F).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void MovPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x10).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void AddPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x58).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void SubPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x5C).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void MulPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x59).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void DivPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x5E).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void MovDDup(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x12).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void AndPD (RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x54).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void AndNPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x55).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void OrPD  (RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x56).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void XorPD (RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x57).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void AddSubPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0xD0).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void UnpckLPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x14).RegMem(dst, src); }
	template <Register32Or64 AddrRegType> void UnpckHPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x15).RegMem(dst, src); }

	template <Register32Or64 AddrRegType> void ShufPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src, uint8_t imm8) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0xC6).RegMem(dst, src).Imm8(imm8); }
	template <Register32Or64 AddrRegType> void DPPD(RegisterXMM dst, const MemoryOperand<AddrRegType> &src, uint8_t imm8) { Byte(0x66).RexAndAddrSize(dst, src).Bytes(0x0F, 0x3A, 0x41).RegMem(dst, src).Imm8(imm8); }

	template <Register32Or64 AddrRegType> void MovSD(const MemoryOperand<AddrRegType> &dst, RegisterXMM src) { Byte(0xF2).RexAndAddrSize(src, dst).Bytes(0x0F, 0x11).RegMem(src, dst); }
	template <Register32Or64 AddrRegType> void MovPD(const MemoryOperand<AddrRegType> &dst, RegisterXMM src) { Byte(0x66).RexAndAddrSize(src, dst).Bytes(0x0F, 0x11).RegMem(src, dst); }

	void Jcc(uint8_t Opcode, size_t Dst) {
		int64_t diff = int64_t(Dst) - int64_t(Pointer);

		if (diff - 2 == int8_t(diff - 2)) {
			Byte(Opcode).Imm8(int8_t(diff - 2));
		} else if (diff - 6 == int32_t(diff - 6)) {
			Bytes(0x0F, Opcode + 0x10).Imm32(int32_t(diff - 6));
		} else throw;
	}

	void Ja(size_t Dst) { Jcc(0x77, Dst); }
	void Jb(size_t Dst) { Jcc(0x72, Dst); }
	void Jg(size_t Dst) { Jcc(0x7F, Dst); }
	void Jl(size_t Dst) { Jcc(0x7C, Dst); }

	void Jae(size_t Dst) { Jcc(0x73, Dst); }
	void Jbe(size_t Dst) { Jcc(0x76, Dst); }
	void Jge(size_t Dst) { Jcc(0x7D, Dst); }
	void Jle(size_t Dst) { Jcc(0x7E, Dst); }

	void Ret() {
		Byte(0xC3);
	}
};
