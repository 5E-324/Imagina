#pragma once
#include <stdint.h>
#include <intrin.h>
#include <type_traits>
#include <bit>
#include <mpirxx.h>
#include <complex>

struct FExpDouble;
struct FExpSingle;
struct FExpI64;

using FloatExp = FExpDouble;

struct FExpDouble {
	union {
		double Mantissa;
		uint64_t MantissaI64;
	};
	int64_t Exponent;

	static constexpr int64_t SignMask = 1ll << 63;
	static constexpr int64_t ExponentMask = 0x7FFll << 52;
	static constexpr int64_t ExponentMSBMask = 0x400ll << 52;
	static constexpr int64_t ExponentZeroOffset = 0x3FFll << 52;
	static constexpr int64_t MantissaMask = (1ll << 52) - 1;

	FExpDouble() : Mantissa(1.0), Exponent(-0x10000000000000) {}

	inline void constexpr Normalize() {
		uint64_t MantissaI64 = std::bit_cast<uint64_t>(Mantissa);
		Exponent = Exponent + int64_t(MantissaI64 << 1 >> 53) - 0x3FF;
		if (Mantissa == 0.0) Exponent = -0x10000000000000;
		MantissaI64 = (MantissaI64 & ~ExponentMask) | ExponentZeroOffset;
		Mantissa = std::bit_cast<double>(MantissaI64);
	}

	inline constexpr FExpDouble(const float			&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const double		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const long double	&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const int8_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const uint8_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const int16_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const uint16_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const int32_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const uint32_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const int64_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	inline constexpr FExpDouble(const uint64_t		&n)	: Mantissa(double(n)), Exponent(0) { Normalize(); }
	
	inline constexpr FExpDouble(const double &mantissa, const int64_t &exponent) : Mantissa(mantissa), Exponent(exponent) { Normalize(); }
	inline constexpr FExpDouble(const double &mantissa, const int64_t &exponent, int) : Mantissa(mantissa), Exponent(exponent) {}

	inline FExpDouble(const mpf_class& n) {
#ifdef mpf_get_2exp_d
		Exponent = mpf_get_2exp_d(&Mantissa, n.get_mpf_t());
#else
		Exponent = 0;
		Mantissa = mpf_get_d_2exp((signed long *)&Exponent, n.get_mpf_t());
#endif
		Normalize();
	}

	inline operator mpf_class() const {
		if (Exponent <= -0x100000000000) return 0.0;
		mpf_class result(Mantissa);
		if (Exponent >= 0) {
			mpf_mul_2exp(result.get_mpf_t(), result.get_mpf_t(), Exponent);
		} else {
			mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), -Exponent);
		}
		return result;
	}

	inline mpf_class to_mpf_class(mp_bitcnt_t prec) const {
		if (Exponent <= -0x100000000000) return 0.0;
		mpf_class result(Mantissa, prec);
		if (Exponent >= 0) {
			mpf_mul_2exp(result.get_mpf_t(), result.get_mpf_t(), Exponent);
		} else {
			mpf_div_2exp(result.get_mpf_t(), result.get_mpf_t(), -Exponent);
		}
		return result;
	}

	inline constexpr explicit operator double() const {
		if (Exponent >= 0x400) {
			int64_t temp = (MantissaI64 | ExponentMask) & ~MantissaMask;
			return reinterpret_cast<double &>(temp);
		} else if (Exponent <= -0x3FF) {
			return 0;
		} else {
			int64_t temp = MantissaI64 + (Exponent << 52);
			return reinterpret_cast<double &>(temp);
		}
	}

	inline constexpr explicit operator int() const {
		return double(*this);
	}

	inline constexpr FExpDouble operator+() const {
		return *this;
	}

	inline constexpr FExpDouble operator-() const {
		return FExpDouble(-Mantissa, Exponent, 0);
	}

	inline FExpDouble operator+(const FExpDouble &n) const {
		__m128i ExponentA = _mm_cvtsi64_si128(Exponent);
		__m128i ExponentB = _mm_cvtsi64_si128(n.Exponent);
		__m128d MantissaA = _mm_load_sd(&Mantissa);
		__m128d MantissaB = _mm_load_sd(&n.Mantissa);
		__m128i Mask = _mm_cmpgt_epi64(ExponentA, ExponentB);
		__m128d MaskD = _mm_castsi128_pd(Mask);
		__m128i ExponentX = _mm_blendv_epi8(ExponentB, ExponentA, Mask);
		__m128i ExponentY = _mm_blendv_epi8(ExponentA, ExponentB, Mask);
		__m128d MantissaX = _mm_blendv_pd(MantissaB, MantissaA, MaskD);
		__m128d MantissaY = _mm_blendv_pd(MantissaA, MantissaB, MaskD);

		__m128i ExponentDiff = _mm_sub_epi64(ExponentX, ExponentY);

		Mask = _mm_cmpgt_epi64(_mm_cvtsi64_si128(0x100), ExponentDiff);
		MantissaY = _mm_castsi128_pd(_mm_and_si128(_mm_sub_epi64(_mm_castpd_si128(MantissaY), _mm_slli_epi64(ExponentDiff, 52)), Mask));
		MantissaX = _mm_add_sd(MantissaX, MantissaY);

		__m128d Zero = _mm_cmpeq_sd(MantissaX, _mm_setzero_pd());
		ExponentX = _mm_sub_epi64(_mm_add_epi64(ExponentX, _mm_srli_epi64(_mm_slli_epi64(_mm_castpd_si128(MantissaX), 1), 53)), _mm_cvtsi64_si128(0x3FF));
		MantissaX = _mm_castsi128_pd(_mm_or_si128(_mm_and_si128(_mm_castpd_si128(MantissaX), _mm_cvtsi64_si128(~ExponentMask)), _mm_cvtsi64_si128(ExponentZeroOffset)));
		ExponentX = _mm_blendv_epi8(ExponentX, _mm_set1_epi64x(-0x10000000000000), _mm_castpd_si128(Zero));

		FExpDouble Result;
		_mm_store_sd(&Result.Mantissa, MantissaX);
		Result.Exponent = _mm_cvtsi128_si64(ExponentX);
		return Result;
	}

	inline FExpDouble operator+=(const FExpDouble &n) {
		return *this = *this + n;
	}

	inline FExpDouble operator-(const FExpDouble &n) const {
		return operator+(-n);
	}

	inline FExpDouble operator-=(const FExpDouble &n) {
		return *this = *this - n;
	}

	constexpr FExpDouble abs() const {
		return FExpDouble(std::abs(Mantissa), Exponent, 0);
	}

	inline constexpr FExpDouble operator*(const FExpDouble &n) const {
		return FExpDouble(Mantissa * n.Mantissa, Exponent + n.Exponent);
	}

	inline constexpr FExpDouble operator*=(const FExpDouble &n) {
		Mantissa *= n.Mantissa;
		Exponent += n.Exponent;
		Normalize();
		return *this;
	}

	inline constexpr FExpDouble operator/(const FExpDouble &n) const {
		return FExpDouble(Mantissa / n.Mantissa, Exponent - n.Exponent);
	}

	inline constexpr FExpDouble operator/=(const FExpDouble &n) {
		Mantissa /= n.Mantissa;
		Exponent -= n.Exponent;
		return *this;
	}

	inline constexpr bool operator==(const FExpDouble &n) const {
		return Exponent == n.Exponent && Mantissa == n.Mantissa;
	}

	inline constexpr bool operator!=(const FExpDouble &n) const {
		return Exponent != n.Exponent || Mantissa != n.Mantissa;
	}

	inline constexpr bool operator>(const FExpDouble &n) const {
		FExpDouble a = *this;
		FExpDouble b = n;
		a.Normalize();
		b.Normalize();
		bool Positive = a.Mantissa > 0;
		if (Positive == b.Mantissa > 0) {
			return Positive == (a.Exponent > b.Exponent || (a.Exponent == b.Exponent && a.Mantissa > b.Mantissa));
		} else {
			return Positive;
		}
	}

	inline constexpr bool operator>=(const FExpDouble &n) const {
		FExpDouble a = *this;
		FExpDouble b = n;
		a.Normalize();
		b.Normalize();
		bool Positive = a.Mantissa > 0;
		if (Positive == b.Mantissa > 0) {
			return Positive == (a.Exponent > b.Exponent || (a.Exponent == b.Exponent && a.Mantissa >= b.Mantissa));
		} else {
			return Positive;
		}
	}

	inline constexpr bool operator<(const FExpDouble &n) const {
		return n > *this;
	}

	inline constexpr bool operator<=(const FExpDouble &n) const {
		return n >= *this;
	}

	FExpDouble operator+(const double &n) const { return operator+(FExpDouble(n)); }
	FExpDouble operator-(const double &n) const { return operator-(FExpDouble(n)); }
	constexpr FExpDouble operator*(const double &n) const { return operator*(FExpDouble(n)); }
	constexpr FExpDouble operator/(const double &n) const { return operator/(FExpDouble(n)); }

	FExpDouble operator+(const int &n) const { return operator+(FExpDouble(n)); }
	FExpDouble operator-(const int &n) const { return operator-(FExpDouble(n)); }
	constexpr FExpDouble operator*(const int &n) const { return operator*(FExpDouble(n)); }
	constexpr FExpDouble operator/(const int &n) const { return operator/(FExpDouble(n)); }

	FExpDouble operator+(const size_t &n) const { return operator+(FExpDouble(n)); }
	FExpDouble operator-(const size_t &n) const { return operator-(FExpDouble(n)); }
	constexpr FExpDouble operator*(const size_t &n) const { return operator*(FExpDouble(n)); }
	constexpr FExpDouble operator/(const size_t &n) const { return operator/(FExpDouble(n)); }

	FExpDouble operator+=(const double &n) { return operator+=(FExpDouble(n)); }
	FExpDouble operator-=(const double &n) { return operator-=(FExpDouble(n)); }
	constexpr FExpDouble operator*=(const double &n) { return operator*=(FExpDouble(n)); }
	constexpr FExpDouble operator/=(const double &n) { return operator/=(FExpDouble(n)); }

	inline constexpr bool operator==(const double &n) const { return operator==(FExpDouble(n)); }
	inline constexpr bool operator!=(const double &n) const { return operator!=(FExpDouble(n)); }
	inline constexpr bool operator>(const double &n) const { return operator>(FExpDouble(n)); }
	inline constexpr bool operator>=(const double &n) const { return operator>=(FExpDouble(n)); }
	inline constexpr bool operator<(const double &n) const { return operator<(FExpDouble(n)); }
	inline constexpr bool operator<=(const double &n) const { return operator<=(FExpDouble(n)); }

	inline constexpr bool operator!=(const int &n) const { return operator!=(FExpDouble(n)); }
	inline constexpr bool operator==(const int &n) const { return operator==(FExpDouble(n)); }
	inline constexpr bool operator>(const int &n) const { return operator>(FExpDouble(n)); }
	inline constexpr bool operator<(const int &n) const { return operator<(FExpDouble(n)); }
	inline constexpr bool operator>=(const int &n) const { return operator>=(FExpDouble(n)); }
	inline constexpr bool operator<=(const int &n) const { return operator<=(FExpDouble(n)); }
};

inline FExpDouble operator+(const double &a, const FExpDouble &b) { return (FExpDouble(a)) + b; }
inline FExpDouble operator-(const double &a, const FExpDouble &b) { return (FExpDouble(a)) - b; }
inline constexpr FExpDouble operator*(const double &a, const FExpDouble &b) { return (FExpDouble(a)) * b; }
inline constexpr FExpDouble operator/(const double &a, const FExpDouble &b) { return (FExpDouble(a)) / b; }

inline FExpDouble operator+(const int &a, const FExpDouble &b) { return (FExpDouble(a)) + b; }
inline FExpDouble operator-(const int &a, const FExpDouble &b) { return (FExpDouble(a)) - b; }
inline constexpr FExpDouble operator*(const int &a, const FExpDouble &b) { return (FExpDouble(a)) * b; }
inline constexpr FExpDouble operator/(const int &a, const FExpDouble &b) { return (FExpDouble(a)) / b; }

inline mpf_class& operator+=(mpf_class& a, const FExpDouble& b) { return a += b.operator mpf_class(); }
inline mpf_class& operator-=(mpf_class& a, const FExpDouble& b) { return a -= b.operator mpf_class(); }
inline mpf_class& operator*=(mpf_class& a, const FExpDouble& b) { return a *= b.operator mpf_class(); }
inline mpf_class& operator/=(mpf_class& a, const FExpDouble& b) { return a /= b.operator mpf_class(); }

inline mpf_class operator+=(const mpf_class& a, const FExpDouble& b) { return a + b.operator mpf_class(); }
inline mpf_class operator-=(const mpf_class& a, const FExpDouble& b) { return a - b.operator mpf_class(); }
inline mpf_class operator*=(const mpf_class& a, const FExpDouble& b) { return a * b.operator mpf_class(); }
inline mpf_class operator/=(const mpf_class& a, const FExpDouble& b) { return a / b.operator mpf_class(); }

inline FExpDouble sqrt(const FExpDouble &n) {
	FExpDouble Result;
	Result.MantissaI64 = n.MantissaI64 + (uint64_t(n.Exponent << 63) >> 11);
	Result.Exponent = n.Exponent >> 1;
	Result.Mantissa = sqrt(Result.Mantissa);
	return Result;
}

inline double log10(const FExpDouble &n) {
	constexpr double Log10_2 = 0.30102999566398119521373889472449;
	return log10(n.Mantissa) + n.Exponent * Log10_2;
}

inline double log(const FExpDouble &n) {
	constexpr double Ln_2 = 0.693147180559945309417;
	return log(n.Mantissa) + n.Exponent * Ln_2;
}

inline std::complex<double> log(const std::complex<FExpDouble> &x) {
	constexpr double Ln_2 = 0.693147180559945309417;
	int64_t CommonExponent;
	double real, imag;
	real = x.real().Mantissa;
	imag = x.imag().Mantissa;

	if (x.real().Exponent >= x.imag().Exponent) {
		CommonExponent = x.real().Exponent;
		if (CommonExponent - x.imag().Exponent <= 0x100) {
			reinterpret_cast<int64_t &>(imag) += (x.imag().Exponent - CommonExponent) << 52;
		} else {
			imag = 0.0;
		}
	} else {
		CommonExponent = x.imag().Exponent;
		if (CommonExponent - x.real().Exponent <= 0x100) {
			reinterpret_cast<int64_t &>(real) += (x.real().Exponent - CommonExponent) << 52;
		} else {
			real = 0.0;
		}
	}

	std::complex<double> Result = { real, imag };
	Result = log(Result);

	Result.real(Result.real() + CommonExponent * Ln_2);

	return Result;
}

inline double log2(const FExpDouble &n) {
	return log2(n.Mantissa) + n.Exponent;
}

template <typename T> T HRExp(double);

template<>
inline FExpDouble HRExp<FExpDouble>(double x) {
	constexpr double Log2_e = 1.44269504088896340736;
	constexpr double Ln_2 = 0.693147180559945309417;
	FExpDouble Result;
	Result.Exponent = (int64_t)floor(x * Log2_e);
	Result.Mantissa = exp(x - Result.Exponent * Ln_2);
	return Result;
}

inline FExpDouble pow2(const double &n) {
	FExpDouble Result;
	Result.Exponent = (int64_t)floor(n);
	Result.Mantissa = pow(2.0, n - Result.Exponent);
	return Result;
}

inline FloatExp abs(const FloatExp &n) {
	return n.abs();
}

inline FloatExp fabs(const FloatExp &n) {
	return n.abs();
}

inline FloatExp copysign(const FloatExp &mag, const FloatExp &sgn) {
	FloatExp Result;
	Result.Exponent = mag.Exponent;
	Result.Mantissa = std::copysign(mag.Mantissa, sgn.Mantissa);
	return Result;
}

inline FloatExp scalbn(const FloatExp &x, int e) {
	FloatExp Result;
	Result.Exponent = x.Exponent + e; // FIXME handle overflow
	Result.Mantissa = x.Mantissa;
	return Result;
}

inline bool isfinite(const FloatExp &x) {
	return isfinite(x.Mantissa);
}

inline bool isnan(const FloatExp &x) {
	return isnan(x.Mantissa);
}

inline bool isinf(const FloatExp &x) {
	return isinf(x.Mantissa);
}

inline FloatExp fmax(const FloatExp &a, const FloatExp &b) {
	return a > b ? a : b;
}

inline FloatExp logb(const FloatExp &x) {
	return logb(x.Mantissa) + x.Exponent;
}
