#include "XMMComplex.h"

void XMMComplexEncoder::Csqr(RegisterXMM dst, RegisterXMM temp1, RegisterXMM temp2) {
	MovPD(temp1, dst);
	MovPD(temp2, dst);

	UnpckHPD(temp1, temp1);
	UnpckLPD(temp2, temp2);
	MulPD(temp1, dst);
	MulPD(dst, temp2);
	ShufPD(temp1, temp1, 1);

	AddSubPD(dst, temp1);
}

void XMMComplexEncoder::Cmul(RegisterXMM dst, RegisterXMM src, RegisterXMM temp) {
	MovPD(temp, src);

	UnpckHPD(src, src);
	UnpckLPD(temp, temp);
	MulPD(src, dst);
	MulPD(dst, temp);
	ShufPD(src, src, 1);

	AddSubPD(dst, src);
}

// src is clobbered

void XMMComplexEncoder::Cdiv(RegisterXMM dst, RegisterXMM src, RegisterXMM temp) {
	MovPD(temp, dst);

	UnpckLPD(temp, temp);
	UnpckHPD(dst, dst);
	MulPD(dst, src);
	MulPD(temp, src);

	DPPD(src, src, 0xFF);
	ShufPD(temp, temp, 1);
	AddSubPD(dst, temp);

	ShufPD(dst, dst, 1);

	DivPD(dst, src);
}

void XMMComplexEncoder::Cpowui(uint32_t exp, RegisterXMM dstsrc, RegisterXMM temp1, RegisterXMM temp2, RegisterXMM temp3) {
	if (exp == 0) throw;

	while ((exp & 1) == 0) {
		Csqr(dstsrc, temp1, temp2);
		exp = exp >> 1;
	}
	if (exp == 1) {
		return;
	}
	MovPD(temp1, dstsrc);
	Csqr(temp1, temp2, temp3);
	exp >>= 1;
	while (exp > 1) {
		if ((exp & 1) != 0) {
			MovPD(temp2, temp1);
			Cmul(dstsrc, temp2, temp3);
		}
		Csqr(temp1, temp2, temp3);
		exp >>= 1;
	}
	Cmul(dstsrc, temp1, temp2);
}
