#pragma once

#include "Encoder.h"

class XMMComplexEncoder : public Encoder {
public:
	XMMComplexEncoder(uint8_t *Buffer, size_t Size) : Encoder(Buffer, Size) {}
	XMMComplexEncoder(uint8_t *Buffer, size_t Size, size_t ConstSize) : Encoder(Buffer, Size, ConstSize) {}

	void Csqr(RegisterXMM dst, RegisterXMM temp1, RegisterXMM temp2);
	void Cmul(RegisterXMM dst, RegisterXMM src, RegisterXMM temp); // src is clobbered
	void Cdiv(RegisterXMM dst, RegisterXMM src, RegisterXMM temp); // src is clobbered
	void Cpowui(uint32_t exp, RegisterXMM dstsrc, RegisterXMM temp1, RegisterXMM temp2, RegisterXMM temp3); // src is clobbered
};