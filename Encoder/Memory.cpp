#include "Memory.h"

uint8_t ScaleToS(uint8_t Scale) {
	switch (Scale) {
		case 1: return 0;
		case 2: return 1;
		case 4: return 2;
		case 8: return 3;
		default: throw;
	}
}

template <Register32Or64 RegType>
MemoryAddress<RegType> &MemoryAddress<RegType>::operator+=(const MemoryAddress<RegType> &x) {
	if (HasI) {
		if (HasB) {
			if (x.HasI || x.HasB) throw;
		} else {
			if (x.HasI) {
				if (x.HasB) throw;
				if (S == 1) {
					B = I;
					I = x.I;
					S = x.S;
				} else if (x.S == 1) {
					B = x.I;
				} else {
					throw;
				}
				HasB = true;
			} else if (x.HasB) {
				B = x.B;
				HasB = true;
			}
		}
	} else if (HasB) {
		if (x.HasI) {
			if (x.HasB) throw;
			I = x.I;
			S = x.S;
			HasI = true;
		} else if (x.HasB) {
			I = x.B;
			S = 1;
			HasI = true;
		}
	}
	Disp += x.Disp;
	return *this;
}

template <Register32Or64 RegType>
MemoryAddress<RegType> &MemoryAddress<RegType>::operator*=(int32_t x) {
	if (x == 0) {
		HasI = false;
		HasB = false;
		Disp = 0;
		return *this;
	}
	if (HasB) {
		if (HasI) throw;
		if (x > 8 || !IsPowerOfTwo(x)) throw;
		I = B;
		S = x;
		HasI = true;
		HasB = false;
	} else if (HasI) {
		if (x > 8 || !IsPowerOfTwo(x)) throw;
		S *= x;
		if (S > 8) throw;
	}
	Disp *= x;
	return *this;
}

template<Register32Or64 RegType>
void MemoryOperand<RegType>::SetDisp(const MemoryAddress<RegType> &Address) {
	if (Address.Disp != 0 || Address.B == RegType(Register64::rbp)) {
		if (Address.Disp == int8_t(Address.Disp)) {
			Mod = 0b01;
			DispSize = 8;
		} else {
			Mod = 0b10;
			DispSize = 32;
		}
	} else {
		Mod = 0b00;
		DispSize = 0;
	}
}

template<Register32Or64 RegType>
void MemoryOperand<RegType>::SetSIB(const MemoryAddress<RegType> &Address) {
	if (Address.I == RegType(Register64::rsp)) throw;
	SetDisp(Address);

	RM = 0b100; // Has SIB
	HasSIB = true;
	S = ScaleToS(Address.S);
	I = uint8_t(Address.I);
	B = uint8_t(Address.B);
}

template<Register32Or64 RegType>
void MemoryOperand<RegType>::SetSI(const MemoryAddress<RegType> &Address) {
	if (Address.I == RegType(Register64::rsp)) throw;
	Mod = 0b00;
	RM = 0b100; // Has SIB
	HasSIB = true;
	S = ScaleToS(Address.S);
	I = uint8_t(Address.I);
	B = 0b101; // Base == 0b101 && Mod == 0b00 -> No base
	DispSize = 32;
}

template<Register32Or64 RegType>
inline void MemoryOperand<RegType>::SetB(const MemoryAddress<RegType> &Address) {
	SetDisp(Address);
	if ((uint8_t(Address.B) & 0b111) == 0b100) { //RegType(Register64::rsp)) {
		RM = 0b100;
		HasSIB = true;
		S = 0b00;
		I = 0b100;
		B = uint8_t(Address.B);
	} else {
		HasSIB = false;
		RM = uint8_t(Address.B);
	}
}

template<Register32Or64 RegType>
MemoryOperand<RegType>::MemoryOperand(const MemoryAddress<RegType> &Address) : Disp(Address.Disp) {
	if (Address.HasI) {
		if (Address.HasB) {
			SetSIB(Address);
		} else {
			SetSI(Address);
		}
	} else if (Address.HasB) {
		SetB(Address);
	} else { // RIP Based!
		Mod = 0b00;
		RM = 0b101;
		HasSIB = false;
		DispSize = 32;
	}
}

template MemoryAddress<Register64> &MemoryAddress<Register64>::operator+=(const MemoryAddress<Register64> &x);
template MemoryAddress<Register32> &MemoryAddress<Register32>::operator+=(const MemoryAddress<Register32> &x);
template MemoryAddress<Register64> &MemoryAddress<Register64>::operator*=(int32_t x);
template MemoryOperand<Register64>::MemoryOperand(const MemoryAddress<Register64> &Address);
template MemoryOperand<Register32>::MemoryOperand(const MemoryAddress<Register32> &Address);