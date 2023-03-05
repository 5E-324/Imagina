#pragma once
#include <intrin.h>
#include <stdint.h>
#include <complex>

struct mask64x4;
struct dvec4;
struct DExpVec4;

struct i64vec4;

template <typename T> struct Vector4Selector { using type = void; };
template <> struct Vector4Selector<double> { using type = dvec4; };
template <> struct Vector4Selector<int64_t> { using type = i64vec4; };
template <> struct Vector4Selector<bool> { using type = mask64x4; };

template <typename T>
using vec4 = typename Vector4Selector<T>::type;

template <typename T>
using complexvec4 = std::complex<typename Vector4Selector<T>::type>;

struct mask64x4 {
	union {
		__m256d ymm;
		uint64_t data[4];
	};

	inline uint64_t &operator[](const size_t &i) { return data[i]; }
	inline uint64_t operator[](const size_t &i) const { return data[i]; }

	inline mask64x4() = default;
	inline mask64x4(__m256d n) : ymm(n) {}
	inline mask64x4(__m256i n) : ymm(_mm256_castsi256_pd(n)) {}
	inline mask64x4(size_t mask) {
		__m256i ymmi = _mm256_set1_epi64x(mask);
		ymmi = _mm256_and_si256(ymmi, _mm256_set_epi64x(0x8, 0x4, 0x2, 0x1));
		ymmi = _mm256_cmpeq_epi64(ymmi, _mm256_set_epi64x(0x8, 0x4, 0x2, 0x1));
		ymm = _mm256_castsi256_pd(ymmi);
	}
	inline mask64x4(bool x) : ymm(_mm256_castsi256_pd(_mm256_set1_epi64x(x ? 0xFFFFFFFFFFFFFFFFull : 0x0ull))) {}
	inline mask64x4(bool x, bool y, bool z, bool w) : ymm(_mm256_castsi256_pd(_mm256_sub_epi64(_mm256_set1_epi64x(0), _mm256_set_epi64x(w, z, y, x)))) {}

	inline operator __m256d() const { return ymm; }
	inline operator __m256i() const { return _mm256_castpd_si256(ymm); }
	inline operator bool() const { return !_mm256_testz_pd(ymm, ymm); }

	inline mask64x4 operator||(const mask64x4& x) const { return mask64x4(_mm256_or_pd(ymm, x.ymm)); }
	inline mask64x4 operator&&(const mask64x4& x) const { return mask64x4(_mm256_and_pd(ymm, x.ymm)); }
};

inline mask64x4 AndNot(mask64x4 x, mask64x4 y) {
	return _mm256_andnot_pd(x, y);
}

inline std::ostream &operator<< (std::ostream &stream, const mask64x4 &mask) {
	stream << '[';
	for (auto x : mask.data) std::cout << (x ? '*' : ' ');
	stream << "]";
	return stream;
}

struct dvec4 {
	union {
		__m256d ymm;
		double data[4];
		struct {
			double x, y, z, w;
		};
	};
	struct mul_result;
	inline dvec4() = default;
	inline dvec4(__m256d n) : ymm(n) {}
	inline dvec4(double n) : ymm(_mm256_set1_pd(n)) {}
	inline dvec4(double x, double y, double z, double w) : ymm(_mm256_set_pd(w, z, y, x)) {}
	inline explicit dvec4(const double *x) : ymm(_mm256_set_pd(x[3], x[2], x[1], x[0])) {}

	inline operator __m256d() const { return ymm; }

	inline i64vec4 casti64() const;
	inline i64vec4 ToInt64() const;

	inline double &operator[](const size_t &i) { return data[i]; }

	inline dvec4 operator+(const dvec4 &n) const { return dvec4(_mm256_add_pd(ymm, n.ymm)); }
	inline dvec4 operator-(const dvec4 &n) const { return dvec4(_mm256_sub_pd(ymm, n.ymm)); }
	inline dvec4 operator*(const dvec4 &n) const;
	inline dvec4 operator*(const double &n) const;
	inline dvec4 operator/(const dvec4 &n) const { return dvec4(_mm256_div_pd(ymm, n.ymm)); }
	inline dvec4 operator/(const double &n) const { return dvec4(_mm256_div_pd(ymm, _mm256_set1_pd(n))); }

	inline dvec4 operator+=(const dvec4 &n) { ymm = _mm256_add_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator-=(const dvec4 &n) { ymm = _mm256_sub_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator*=(const dvec4 &n) { ymm = _mm256_mul_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator/=(const dvec4 &n) { ymm = _mm256_div_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator*=(const double &n) { ymm = _mm256_mul_pd(ymm, _mm256_set1_pd(n)); return *this; }
	inline dvec4 operator/=(const double &n) { ymm = _mm256_div_pd(ymm, _mm256_set1_pd(n)); return *this; }

	inline dvec4 operator+(const mul_result &n) const;
	inline dvec4 operator-(const mul_result &n) const;

	inline mask64x4 operator<(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_LT_OQ)); }
	inline mask64x4 operator>(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_GT_OQ)); }
	inline mask64x4 operator<=(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_LE_OQ)); }
	inline mask64x4 operator>=(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_GE_OQ)); }
	inline mask64x4 operator==(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_EQ_OQ)); }
	inline mask64x4 operator!=(const dvec4 &n) const { return mask64x4(_mm256_cmp_pd(ymm, n.ymm, _CMP_NEQ_OQ)); }

	inline dvec4 operator&(const mask64x4 &n) const { return dvec4(_mm256_and_pd(ymm, n.ymm)); }
	inline dvec4 operator|(const mask64x4 &n) const { return dvec4(_mm256_or_pd(ymm, n.ymm)); }
	inline dvec4 operator^(const mask64x4 &n) const { return dvec4(_mm256_xor_pd(ymm, n.ymm)); }

	inline dvec4 operator&(const dvec4 &n) const { return dvec4(_mm256_and_pd(ymm, n.ymm)); }
	inline dvec4 operator|(const dvec4 &n) const { return dvec4(_mm256_or_pd(ymm, n.ymm)); }
	inline dvec4 operator^(const dvec4 &n) const { return dvec4(_mm256_xor_pd(ymm, n.ymm)); }

	inline dvec4 operator&=(const dvec4 &n) { ymm = _mm256_and_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator|=(const dvec4 &n) { ymm = _mm256_or_pd(ymm, n.ymm); return *this; }
	inline dvec4 operator^=(const dvec4 &n) { ymm = _mm256_xor_pd(ymm, n.ymm); return *this; }

	inline dvec4 abs() const { return dvec4(_mm256_andnot_pd(_mm256_set1_pd(-0.0), ymm)); }
	inline dvec4 operator+() const { return *this; }
	inline dvec4 operator-() const { return dvec4(_mm256_xor_pd(ymm, _mm256_set1_pd(-0.0))); }
};

struct dvec4::mul_result {
	struct {
		__m256d Operand1, Operand2;
		inline operator __m256d() { return _mm256_mul_pd(Operand1, Operand2); }
	} ymm;
	inline mul_result(const __m256d &a, const __m256d &b) : ymm{ a, b } {}
	inline operator dvec4() { return dvec4(_mm256_mul_pd(ymm.Operand1, ymm.Operand2)); }

	inline dvec4 operator+(const dvec4 &n) const { return dvec4(_mm256_fmadd_pd(ymm.Operand1, ymm.Operand2, n.ymm)); }
	inline dvec4 operator-(const dvec4 &n) const { return dvec4(_mm256_fmsub_pd(ymm.Operand1, ymm.Operand2, n.ymm)); }
	inline dvec4 operator*(const dvec4 &n) { return dvec4(*this) * n; }
	inline dvec4 operator/(const dvec4 &n) { return dvec4(*this) / n; }
	inline dvec4 operator*(const double &n) { return dvec4(*this) * n; }
	inline dvec4 operator/(const double &n) { return dvec4(*this) / n; }
};

inline dvec4 dvec4::operator*(const dvec4 &n) const { return _mm256_mul_pd(ymm, n.ymm); }
inline dvec4 dvec4::operator*(const double &n) const { return _mm256_mul_pd(ymm, _mm256_set1_pd(n)); }

inline dvec4 dvec4::operator+(const dvec4::mul_result &n) const { return _mm256_fmadd_pd(n.ymm.Operand1, n.ymm.Operand2, ymm); }
inline dvec4 dvec4::operator-(const dvec4::mul_result &n) const { return _mm256_fnmadd_pd(n.ymm.Operand1, n.ymm.Operand2, ymm); }

inline dvec4 sqrt(const dvec4& n) { return _mm256_sqrt_pd(n); }
inline dvec4 log(const dvec4& n) {
#if _MSC_VER >= 1920
#ifndef __clang__
	return _mm256_log_pd(n);
#else
	register __m256d ymm0 asm("ymm0") = n.ymm;

	asm(
		"call __vdecl_log4 \t\n"
		: "=v" (ymm0)
		: "0" (ymm0)
		: "%ymm1", "%ymm2", "%ymm3", "%ymm4", "%ymm5", "%rax", "%rcx", "%rdx", "%r8", "%r9", "%r10", "%r11"
	);
	return ymm0;
#endif
#else
	return dvec4(log(n.x), log(n.y), log(n.z), log(n.w));
#endif
}
inline dvec4 log2(const dvec4 &n) {
#if _MSC_VER >= 1920
#ifndef __clang__
	return _mm256_log2_pd(n);
#else
	register __m256d ymm0 asm("ymm0") = n.ymm;

	asm(
		"call __vdecl_log24 \t\n"
		: "=v" (ymm0)
		: "0" (ymm0)
		: "%ymm1", "%ymm2", "%ymm3", "%ymm4", "%ymm5", "%rax", "%rcx", "%rdx", "%r8", "%r9", "%r10", "%r11"
	);
	return ymm0;
#endif
#else
	return dvec4(log2(n.x), log2(n.y), log2(n.z), log2(n.w));
#endif
}
inline dvec4 exp(const dvec4 &n) {
#if _MSC_VER >= 1920
#ifndef __clang__
	return _mm256_exp_pd(n);
#else
	register __m256d ymm0 asm("ymm0") = n.ymm;

	asm(
		"call __vdecl_exp4 \t\n"
		: "=v" (ymm0)
		: "0" (ymm0)
		: "%ymm1", "%ymm2", "%ymm3", "%ymm4", "%ymm5", "%rax", "%rcx", "%rdx", "%r8", "%r9", "%r10", "%r11"
	);
	return ymm0;
#endif
#else
	return dvec4(exp(n.x), exp(n.y), exp(n.z), exp(n.w));
#endif
}
inline dvec4 min(const dvec4& a, const dvec4& b) { return _mm256_min_pd(a, b); }
inline dvec4 max(const dvec4& a, const dvec4& b) { return _mm256_max_pd(a, b); }

inline dvec4 floor(const dvec4 &n) {
	return _mm256_floor_pd(n);
}

struct i64vec4 {
	union {
		__m256i ymm;
		int64_t data[4];
		struct {
			int64_t x, y, z, w;
		};
	};
	
	inline			i64vec4() = default;
	inline			i64vec4(__m256i n)									: ymm(n) {}
	inline			i64vec4(int64_t n)									: ymm(_mm256_set1_epi64x(n)) {}
	inline			i64vec4(int64_t x, int64_t y, int64_t z, int64_t w)	: ymm(_mm256_set_epi64x(w, z, y, x)) {}
	inline explicit	i64vec4(mask64x4 n)									: ymm(_mm256_castpd_si256(n.ymm)) {}

	inline operator __m256i() const { return ymm; }

	inline dvec4 castd() const { return _mm256_castsi256_pd(ymm); }
	inline static __m256d int64_to_double(__m256i x) {
		x = _mm256_add_epi64(x, _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000)));
		return _mm256_sub_pd(_mm256_castsi256_pd(x), _mm256_set1_pd(0x0018000000000000));
	}
	inline dvec4 ToDouble() const {
		__m256i temp = _mm256_add_epi64(ymm, _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000)));
		return _mm256_sub_pd(_mm256_castsi256_pd(temp), _mm256_set1_pd(0x0018000000000000));
	}

	inline int64_t &operator[](const size_t &i) { return data[i]; }

	inline i64vec4 operator+(const i64vec4& n)	const	{ return i64vec4(_mm256_add_epi64(ymm, n.ymm));						}
	inline i64vec4 operator-(const i64vec4& n)	const	{ return i64vec4(_mm256_sub_epi64(ymm, n.ymm));						}
	
	inline i64vec4 srl(const i64vec4& n)		const	{ return i64vec4(_mm256_srlv_epi64(ymm, n.ymm));					}
	inline i64vec4 srl(const int& n)			const	{ return i64vec4(_mm256_srli_epi64(ymm, n));						}

	inline i64vec4 operator<<(const i64vec4& n)	const	{ return i64vec4(_mm256_sllv_epi64(ymm, n.ymm));					}
	inline i64vec4 operator<<(const int& n)		const	{ return i64vec4(_mm256_slli_epi64(ymm, n));						}
	inline i64vec4 operator>>(const i64vec4 &n) const {
		__m256i SignMask = *this < i64vec4(0);
		return srl(n) | (i64vec4(SignMask) << (i64vec4(64) - n));
	}
	inline i64vec4 operator>>(const int &n) const {
		__m256i SignMask = *this < i64vec4(0);
		return srl(n) | (i64vec4(SignMask) << i64vec4(64 - n));
	}

	inline i64vec4 operator&(const i64vec4& n)	const	{ return i64vec4(_mm256_and_si256(ymm, n.ymm));						}
	inline i64vec4 operator|(const i64vec4& n)	const	{ return i64vec4(_mm256_or_si256 (ymm, n.ymm));						}
	inline i64vec4 operator^(const i64vec4& n)	const	{ return i64vec4(_mm256_xor_si256(ymm, n.ymm));						}
	inline i64vec4 operator~()					const	{ return i64vec4(_mm256_xor_si256(ymm, _mm256_set1_epi64x(~0ull)));	}

	inline i64vec4 operator+=(const i64vec4& n)			{ ymm = _mm256_add_epi64(ymm, n.ymm);	return *this;				}
	inline i64vec4 operator-=(const i64vec4& n)			{ ymm = _mm256_sub_epi64(ymm, n.ymm);	return *this;				}

	inline i64vec4 operator<<=(const i64vec4& n)		{ ymm = _mm256_sllv_epi64(ymm, n.ymm);	return *this;				}
	inline i64vec4 operator>>=(const i64vec4& n)		{ ymm = _mm256_srlv_epi64(ymm, n.ymm);	return *this;				}
	inline i64vec4 operator<<=(const int& n)			{ ymm = _mm256_slli_epi64(ymm, n);		return *this;				}
	inline i64vec4 operator>>=(const int& n)			{ ymm = _mm256_srli_epi64(ymm, n);		return *this;				}

	inline i64vec4 operator&=(const i64vec4& n)			{ ymm = _mm256_and_si256(ymm, n.ymm);	return *this;				}
	inline i64vec4 operator|=(const i64vec4& n)			{ ymm = _mm256_or_si256(ymm, n.ymm);	return *this;				}
	inline i64vec4 operator^=(const i64vec4& n)			{ ymm = _mm256_xor_si256(ymm, n.ymm);	return *this;				}

	inline mask64x4 operator< (const i64vec4 &n) const { return mask64x4(_mm256_cmpgt_epi64(n.ymm, ymm)); }
	inline mask64x4 operator> (const i64vec4 &n) const { return mask64x4(_mm256_cmpgt_epi64(ymm, n.ymm)); }
	inline mask64x4 operator<=(const i64vec4 &n) const { return mask64x4(_mm256_xor_si256(_mm256_cmpgt_epi64(ymm, n.ymm), _mm256_set1_epi64x(~0ll))); }
	inline mask64x4 operator>=(const i64vec4 &n) const { return mask64x4(_mm256_xor_si256(_mm256_cmpgt_epi64(n.ymm, ymm), _mm256_set1_epi64x(~0ll))); }
	inline mask64x4 operator==(const i64vec4 &n) const { return mask64x4(_mm256_cmpeq_epi64(ymm, n.ymm)); }
	inline mask64x4 operator!=(const i64vec4 &n) const { return mask64x4(_mm256_xor_si256(_mm256_cmpeq_epi64(ymm, n.ymm), _mm256_set1_epi64x(~0ll))); }

	inline int64_t HorizontalMin() const { return std::min(std::min(x, y), std::min(z, w)); }
};

inline i64vec4 operator&(mask64x4 x, i64vec4 y) { return _mm256_and_si256(x, y); }
inline i64vec4 operator|(mask64x4 x, i64vec4 y) { return _mm256_or_si256 (x, y); }
inline i64vec4 operator^(mask64x4 x, i64vec4 y) { return _mm256_xor_si256(x, y); }
inline i64vec4 operator&(i64vec4 x, mask64x4 y) { return _mm256_and_si256(x, y); }
inline i64vec4 operator|(i64vec4 x, mask64x4 y) { return _mm256_or_si256 (x, y); }
inline i64vec4 operator^(i64vec4 x, mask64x4 y) { return _mm256_xor_si256(x, y); }

inline std::ostream &operator<< (std::ostream &stream, const i64vec4 &vec) {
	stream << '[';
	stream << vec.data[0] << ',';
	stream << vec.data[1] << ',';
	stream << vec.data[2] << ',';
	stream << vec.data[3];
	stream << "]";
	return stream;
}

inline mask64x4 Select(mask64x4 sel, mask64x4 t, mask64x4 f) { return _mm256_blendv_pd(f, t, sel); }
inline dvec4 Select(mask64x4 sel, dvec4 t, dvec4 f) { return _mm256_blendv_pd(f, t, sel); }
inline i64vec4 Select(mask64x4 sel, i64vec4 t, i64vec4 f) { return _mm256_blendv_epi8(f, t, sel); }
inline DExpVec4 Select(mask64x4 sel, DExpVec4 t, DExpVec4 f);
#include "FloatExpVector4.h"

inline i64vec4 dvec4::casti64() const { return _mm256_castpd_si256(ymm); }
inline i64vec4 dvec4::ToInt64() const {
	__m256d temp = _mm256_add_pd(ymm, _mm256_set1_pd(0x0018000000000000));
	return _mm256_sub_epi64(_mm256_castpd_si256(temp), _mm256_castpd_si256(_mm256_set1_pd(0x0018000000000000)));
}

template<typename T>
inline std::complex<T> Select(mask64x4 sel, std::complex<T> t, std::complex<T> f) {
	return std::complex<T>(Select(sel, t.real(), f.real()), Select(sel, t.imag(), f.imag()));
}

inline dvec4 ArrayToVec4(double *x) {
	return dvec4(x);
}

inline bool isinf(const dvec4 &x) { return false; }
inline bool isnan(const dvec4 &x) { return false; }
inline dvec4 copysign(const dvec4 &mag, const dvec4 &sgn) {
	return dvec4(_mm256_or_pd(mag.abs().ymm, _mm256_and_pd(_mm256_set1_pd(-0.0), sgn.ymm)));
}
