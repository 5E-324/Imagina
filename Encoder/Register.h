#pragma once
#include <stdint.h>
#include <concepts>

enum class Register16 : uint8_t {
	ax = 0,
	cx = 1,
	dx = 2,
	bx = 3,
	sp = 4,
	bp = 5,
	si = 6,
	di = 7,
	r8w = 8,
	r9w = 9,
	r10w = 10,
	r11w = 11,
	r12w = 12,
	r13w = 13,
	r14w = 14,
	r15w = 15,
};

enum class Register32 : uint8_t {
	eax = 0,
	ecx = 1,
	edx = 2,
	ebx = 3,
	esp = 4,
	ebp = 5,
	esi = 6,
	edi = 7,
	r8d = 8,
	r9d = 9,
	r10d = 10,
	r11d = 11,
	r12d = 12,
	r13d = 13,
	r14d = 14,
	r15d = 15,
};

enum class Register64 : uint8_t {
	rax = 0,
	rcx = 1,
	rdx = 2,
	rbx = 3,
	rsp = 4,
	rbp = 5,
	rsi = 6,
	rdi = 7,
	r8 = 8,
	r9 = 9,
	r10 = 10,
	r11 = 11,
	r12 = 12,
	r13 = 13,
	r14 = 14,
	r15 = 15,
};

enum class RegisterXMM : uint8_t {
	xmm0 = 0,
	xmm1 = 1,
	xmm2 = 2,
	xmm3 = 3,
	xmm4 = 4,
	xmm5 = 5,
	xmm6 = 6,
	xmm7 = 7,
	xmm8 = 8,
	xmm9 = 9,
	xmm10 = 10,
	xmm11 = 11,
	xmm12 = 12,
	xmm13 = 13,
	xmm14 = 14,
	xmm15 = 15,
};

template <typename T> inline constexpr bool IsRegisterType = false;
template <> inline constexpr bool IsRegisterType<Register16> = true;
template <> inline constexpr bool IsRegisterType<Register32> = true;
template <> inline constexpr bool IsRegisterType<Register64> = true;
template <> inline constexpr bool IsRegisterType<RegisterXMM> = true;

template <typename T> constexpr bool IsRegisterTypeOrUint8 = IsRegisterType<T> || std::convertible_to<T, uint8_t>;

template <typename T> concept Register = IsRegisterType<T>;
template <typename T> concept RegisterOrInt8 = IsRegisterTypeOrUint8<T>;
template <typename T> concept Register64OrInt8 = std::is_same_v<T, Register64> || std::convertible_to<T, uint8_t>;
template <typename T> concept Register32Or64 = std::is_same_v<T, Register64> || std::is_same_v<T, Register32>;

template <typename T>
constexpr size_t RegisterSize = 0;

template <> inline constexpr size_t RegisterSize<Register16> = 16;
template <> inline constexpr size_t RegisterSize<Register32> = 32;
template <> inline constexpr size_t RegisterSize<Register64> = 64;
