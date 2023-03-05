#pragma once

#include <intrin.h>
#include <limits>
#include <stdint.h>
#include "FloatExp.h"

#include "Vector4.h"
struct DExpVec;

template <> struct Vector4Selector<FExpDouble> { using type = DExpVec4; };

struct DExpVec4 {
	union {
		struct {
			dvec4 Mantissa;
			i64vec4 Exponent;
		};
		struct {
			dvec4 Mantissa;
			i64vec4 Exponent;

		} AssumePositive;
		struct {
			dvec4 Mantissa;
			i64vec4 Exponent;

			__m256d operator==(const DExpVec4 &x) {
				return (Mantissa == x.Mantissa) && (Exponent == x.Exponent);
			}
		} AssumeNormalized;
		struct {
			dvec4 Mantissa;
			i64vec4 Exponent;

			mask64x4 operator> (const DExpVec4 &x) const { return (Exponent > x.Exponent) || ((Exponent == x.Exponent) && (Mantissa >  x.Mantissa)); }
			mask64x4 operator>=(const DExpVec4 &x) const { return (Exponent > x.Exponent) || ((Exponent == x.Exponent) && (Mantissa >= x.Mantissa)); }
			mask64x4 operator< (const DExpVec4 &x) const { return (Exponent < x.Exponent) || ((Exponent == x.Exponent) && (Mantissa <  x.Mantissa)); }
			mask64x4 operator<=(const DExpVec4 &x) const { return (Exponent < x.Exponent) || ((Exponent == x.Exponent) && (Mantissa <= x.Mantissa)); }
		} AssumeNormalizedPositive;
	};
	
	mask64x4 operator==(const DExpVec4 &x) const { return (Exponent == x.Exponent) && (Mantissa == x.Mantissa); }
	mask64x4 operator!=(const DExpVec4 &x) const { return (Exponent != x.Exponent) || (Mantissa != x.Mantissa); }
	mask64x4 operator> (const DExpVec4 &x) const { return (Exponent > x.Exponent) || ((Exponent == x.Exponent) && (Mantissa >  x.Mantissa)); }
	mask64x4 operator>=(const DExpVec4 &x) const { return (Exponent > x.Exponent) || ((Exponent == x.Exponent) && (Mantissa >= x.Mantissa)); }
	mask64x4 operator< (const DExpVec4 &x) const { return (Exponent < x.Exponent) || ((Exponent == x.Exponent) && (Mantissa <  x.Mantissa)); }
	mask64x4 operator<=(const DExpVec4 &x) const { return (Exponent < x.Exponent) || ((Exponent == x.Exponent) && (Mantissa <= x.Mantissa)); }

	static constexpr double SignMask = -0.0;										//0x8000000000000000
	static constexpr double ExponentMask = std::numeric_limits<double>::infinity();	//0x7FF0000000000000
	static constexpr double ExponentZeroOffset = 1.0;								//0x3FF0000000000000
	static constexpr int64_t ExponentZeroOffsetI64 = 0x3FF;

	inline void Normalize() noexcept {
		i64vec4 MantissaI64 = Mantissa.casti64();
		i64vec4 MantissaExponent = _mm256_srli_epi64(_mm256_slli_epi64(MantissaI64, 1), 53);

		Exponent = Exponent + MantissaExponent - i64vec4(ExponentZeroOffsetI64);
		Exponent = Select(Mantissa == dvec4(0.0), i64vec4(-0x100000000000), Exponent);
		Mantissa = _mm256_or_pd(_mm256_andnot_pd(_mm256_set1_pd(ExponentMask), Mantissa), _mm256_set1_pd(ExponentZeroOffset));
	}

	inline DExpVec4 Normalized() noexcept {
		DExpVec4 Result;
		i64vec4 MantissaI64 = Mantissa.casti64();
		i64vec4 MantissaExponent = _mm256_srli_epi64(_mm256_slli_epi64(MantissaI64, 1), 53);

		Result.Exponent = Exponent + MantissaExponent - i64vec4(ExponentZeroOffsetI64);
		Result.Exponent = Select(Mantissa == dvec4(0.0), i64vec4(-0x100000000000), Result.Exponent);
		Result.Mantissa = _mm256_or_pd(_mm256_andnot_pd(_mm256_set1_pd(ExponentMask), Mantissa), _mm256_set1_pd(ExponentZeroOffset));

		return Result;
	}

	inline DExpVec4 NormalizedSimplified() noexcept {
		DExpVec4 Result;
		__m256i MantissaI64 = _mm256_castpd_si256(Mantissa);
		__m256i MantissaExponent = _mm256_srli_epi64(_mm256_slli_epi64(MantissaI64, 1), 53);
		Result.Exponent = _mm256_sub_epi64(_mm256_add_epi64(Exponent, MantissaExponent), _mm256_set1_epi64x(ExponentZeroOffsetI64));
		Result.Mantissa = _mm256_or_pd(_mm256_andnot_pd(_mm256_set1_pd(ExponentMask), Mantissa), _mm256_set1_pd(ExponentZeroOffset));
		return Result;
	}

	DExpVec4() : Mantissa(_mm256_set1_pd(1.0)), Exponent(-0x100000000000) {};

	inline DExpVec4(double x) noexcept : Mantissa(_mm256_set1_pd(x)), Exponent(_mm256_setzero_si256()) { Normalize(); }
	inline DExpVec4(int x) noexcept : DExpVec4(double(x)) {}
	inline DExpVec4(double x, double y, double z, double w) noexcept
		: Mantissa(_mm256_set_pd(w, z, y, x)), Exponent(_mm256_setzero_si256()) {
		Normalize();
	}
	inline DExpVec4(FExpDouble x) noexcept : Mantissa(_mm256_set1_pd(x.Mantissa)), Exponent(_mm256_set1_epi64x(x.Exponent)) { Normalize(); }
	inline DExpVec4(FExpDouble x, FExpDouble y, FExpDouble z, FExpDouble w) noexcept :
		Mantissa(_mm256_set_pd(w.Mantissa, z.Mantissa, y.Mantissa, x.Mantissa)),
		Exponent(_mm256_set_epi64x(w.Exponent, z.Exponent, y.Exponent, x.Exponent)) {}
	inline DExpVec4(FExpDouble *x) noexcept :
		Mantissa(_mm256_set_pd(x[3].Mantissa, x[2].Mantissa, x[1].Mantissa, x[0].Mantissa)),
		Exponent(_mm256_set_epi64x(x[3].Exponent, x[2].Exponent, x[1].Exponent, x[0].Exponent)) {
	}

	inline DExpVec4(dvec4 x) noexcept : Mantissa(x), Exponent(0) { Normalize(); }

	inline explicit DExpVec4(dvec4 mantissa, i64vec4 exponent) noexcept : Mantissa(mantissa), Exponent(exponent) {}

	inline explicit operator dvec4() {
		Normalize();
		return ((Mantissa.casti64() + (Exponent << 52)) & (Exponent > i64vec4(-0x3FF))).castd();
	}

	inline FExpDouble operator[](const size_t &i) { return FExpDouble(Mantissa[i], Exponent[i]); }

	inline DExpVec4 operator+(const DExpVec4 &x) const noexcept {
		mask64x4 Mask = _mm256_cmpgt_epi64(Exponent, x.Exponent);
		i64vec4 Exponent1 = Select(Mask, Exponent, x.Exponent);
		i64vec4 Exponent2 = Select(Mask, x.Exponent, Exponent);
		dvec4 Mantissa1 = Select(Mask, Mantissa, x.Mantissa);
		dvec4 Mantissa2 = Select(Mask, x.Mantissa, Mantissa);

		i64vec4 ExponentDiff = Exponent1 - Exponent2;

		Mask = ExponentDiff < i64vec4(100);
		Mantissa2 = (Mantissa2.casti64() - (ExponentDiff << 52)).castd();
		Mantissa2 = _mm256_and_pd(Mantissa2, Mask);
		Mantissa1 = _mm256_add_pd(Mantissa1, Mantissa2);

		return DExpVec4(Mantissa1, Exponent1).Normalized();
	}

	inline DExpVec4 operator-(const DExpVec4 &x) const noexcept {
		return operator+(DExpVec4(_mm256_xor_pd(x.Mantissa, _mm256_set1_pd(SignMask)), x.Exponent));
	}

	inline DExpVec4 operator*(const DExpVec4 &x) const noexcept {
		return DExpVec4(Mantissa * x.Mantissa, Exponent + x.Exponent);
	}

	inline DExpVec4 operator/(const DExpVec4 &x) const noexcept {
		return DExpVec4(Mantissa / x.Mantissa, Exponent - x.Exponent);
	}

	inline DExpVec4 operator*(const FExpDouble &x) const noexcept {
		return DExpVec4(Mantissa * dvec4(x.Mantissa), Exponent + i64vec4(x.Exponent));
	}

	inline DExpVec4 operator/(const FExpDouble &x) const noexcept {
		return DExpVec4(Mantissa / dvec4(x.Mantissa), Exponent - i64vec4(x.Exponent));
	}

	inline DExpVec4 operator*=(const DExpVec4 &x) noexcept {
		Mantissa *= x.Mantissa;
		Exponent += x.Exponent;
		return *this;
	}

	inline DExpVec4 operator/=(const DExpVec4 &x) noexcept {
		Mantissa /= x.Mantissa;
		Exponent -= x.Exponent;
		return *this;
	}

	inline DExpVec4 operator*=(const FExpDouble &x) noexcept {
		Mantissa *= dvec4(x.Mantissa);
		Exponent += i64vec4(x.Exponent);
		return *this;
	}

	inline DExpVec4 operator/=(const FExpDouble &x) noexcept {
		Mantissa /= dvec4(x.Mantissa);
		Exponent -= i64vec4(x.Exponent);
		return *this;
	}

	inline DExpVec4 abs() const noexcept { return DExpVec4(_mm256_andnot_pd(_mm256_set1_pd(-0.0), Mantissa), Exponent); }
	inline DExpVec4 operator+() const { return *this; }
	inline DExpVec4 operator-() const { return DExpVec4(-Mantissa, Exponent); }

	static inline DExpVec4 absmax(const DExpVec4 &x, const DExpVec4 &y) noexcept {
		DExpVec4 a = x.abs();
		DExpVec4 b = y.abs();

		return Select(a.AssumeNormalizedPositive > b, a, b);
	}
};

inline DExpVec4 Select(mask64x4 sel, DExpVec4 t, DExpVec4 f) {
	return DExpVec4(Select(sel, t.Mantissa, f.Mantissa), Select(sel, t.Exponent, f.Exponent));
}

inline DExpVec4 sqrt(const DExpVec4& n) {
	DExpVec4 Result;
	Result.Mantissa = (n.Mantissa.casti64() + (n.Exponent << 63).srl(11)).castd();
	Result.Exponent = n.Exponent >> 1;
	Result.Mantissa = sqrt(Result.Mantissa);
	return Result;
}
#ifndef M_LN2
#define M_LN2      0.693147180559945309417  // ln(2)
#endif

inline dvec4 log(const DExpVec4& n) {
	return log(n.Mantissa) + n.Exponent.ToDouble() * dvec4(M_LN2);
}

inline dvec4 log2(const DExpVec4 &n) {
	return log2(n.Mantissa) + n.Exponent.ToDouble();
}

inline DExpVec4 min(const DExpVec4& a, const DExpVec4& b) {
	return Select(a.AssumeNormalizedPositive < b, a, b);
}

inline DExpVec4 max(const DExpVec4& a, const DExpVec4& b) {
	return Select(a.AssumeNormalizedPositive > b, a, b);
}

template <typename T> T HRExp(dvec4);

template<>
inline DExpVec4 HRExp<DExpVec4>(dvec4 x) {
	constexpr double Log2_e = 1.44269504088896340736;
	constexpr double Ln_2 = 0.693147180559945309417;
	DExpVec4 Result;
	Result.Exponent = floor(x * Log2_e).ToInt64();
	Result.Mantissa = exp(x - Result.Exponent.ToDouble() * Ln_2);
	return Result;
}

inline DExpVec4 ArrayToVec4(FExpDouble *x) {
	return DExpVec4(x);
}