#ifndef _ARBITRARY_PRECISION_H_
#define _ARBITRARY_PRECISION_H_

#include <mpir.h>
#include <mpirxx.h>

#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
#include <stdint.h>
#include <algorithm>
#include "FloatExp.h"

template <typename T> T Power2(const int64_t &Power);

template <> inline double Power2<double>(const int64_t &Power) {
	uint64_t Exponent = (0x3FF + Power) << 52;
	return reinterpret_cast<const double &>(Exponent);
}
template <> inline long double Power2<long double>(const int64_t &Power) {
	struct { uint64_t Mantissa; uint16_t SignExponent; } x;
	x.SignExponent = uint16_t(0x3FFF + Power);
	x.Mantissa = 0x8000000000000000;
	return reinterpret_cast<const long double &>(x);
}

template <typename T>
struct OMITemp;

template <typename T>
struct OMIMTTemp;
constexpr int OMI_PROLOGUE = -1;
constexpr int OMI_EPILOGUE = -2;

inline void SetDefaultPrecision(size_t precision) {
	mpf_set_default_prec(precision);
}

template<>
struct OMITemp<mpf_class> {
	mpf_class temp[2];
};

inline void OMI(mpf_class &ZR, mpf_class &ZI, const mpf_class &CR, const mpf_class &CI, OMITemp<mpf_class> &T) {
	T.temp[0] = ZR + ZI;
	T.temp[1] = ZR - ZI;
	ZI *= ZR;
	ZI += ZI;
	ZI += CI;
	
	ZR = T.temp[0] * T.temp[1] + CR;
}

template<>
struct OMIMTTemp<mpf_class> {
	mpf_class temp[3];
};
template <int ThreadIndex>
inline void OMIMT(mpf_class &ZR, mpf_class &ZI, const mpf_class &CR, const mpf_class &CI, OMIMTTemp<mpf_class> &T) {
	if constexpr (ThreadIndex == OMI_PROLOGUE) {
		T.temp[0] = ZR + ZI;
		T.temp[1] = ZR - ZI;
		T.temp[2] = ZR;
	} else if constexpr (ThreadIndex == 0) {
		ZI = ZI * T.temp[2] * 2 + CI;
	} else if constexpr (ThreadIndex == 1) {
		ZR = T.temp[0] * T.temp[1] + CR;
	} else if constexpr (ThreadIndex == OMI_EPILOGUE) {
		// No epilogue
	}
	static_assert(ThreadIndex >= -2 && ThreadIndex < 2);
}

#endif