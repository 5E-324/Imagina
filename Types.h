#ifndef _TYPES_H_
#define _TYPES_H_

#include <mpirxx.h>
#include "FloatExp.h"
#include "ArbitraryPrecision.h"
#include "Vector4.h"
#include "FloatExpVector4.h"

template <size_t n, typename T> struct VectorSelector { using type = void; };
template <typename T> struct VectorSelector<1, T> { using type = T; };
template <typename T> struct VectorSelector<4, T> { using type = typename Vector4Selector<T>::type; };

template <size_t n, typename T>
using vec = typename VectorSelector<n, T>::type;

template <size_t n, typename T>
using complexvec = std::complex<typename VectorSelector<n, T>::type>;

using HPReal = mpf_class;

using HRReal = FloatExp;
using SRReal = double;

using SRComplex = std::complex<SRReal>;
using HRComplex = std::complex<HRReal>;

template <size_t n> using VHRReal = vec<n, HRReal>;
template <size_t n> using VSRReal = vec<n, SRReal>;

template <size_t n> using VHRComplex = complexvec<n, HRReal>;
template <size_t n> using VSRComplex = complexvec<n, SRReal>;

template <size_t n> using VBool = vec<n, bool>;
template <size_t n> using VInt64 = vec<n, int64_t>;

using Int = intptr_t;
using Uint = uintptr_t;

using IntIter = int_fast64_t;
using UintIter = uint_fast64_t;

using IntPix = int_least32_t;
using UintPix = uint_least32_t;

constexpr size_t VectorSize = 4;

template <size_t n> VHRReal<n> ArrayToVector(HRReal *x);

template <> inline VHRReal<1> ArrayToVector<1>(HRReal *x) { return *x; }
template <> inline VHRReal<4> ArrayToVector<4>(HRReal *x) { return ArrayToVec4(x); }

template <size_t n> VInt64<n> VInt64IntegerSequence();

template <> inline VInt64<1> VInt64IntegerSequence<1>() { return 0; }
template <> inline VInt64<4> VInt64IntegerSequence<4>() { return VInt64<4>(0, 1, 2, 3); }

constexpr struct {
	IntPix x, y;
} VectorSize2D{ 2, 2 };

constexpr SRReal operator""_sr(long double x) {
	return x;
}

constexpr HRReal operator""_hr(long double x) {
	return x;
}

inline HPReal operator""_hp(long double x) {
	return HPReal(double(x));
}

inline HPReal operator""_hp(unsigned long long x) {
	return HPReal(x);
}

inline double ToReal(const mpf_class& n) {
	return n.get_d();
}

struct RelRect {
	HRReal X, Y;
	HRReal HalfW, HalfH;
};

struct RelLocation {
	HRReal X, Y;
	HRReal HalfH;
};

struct Rect {
	SRReal X, Y;
	SRReal Width, Height;
};

struct TextureDescription {
	GLuint Texture;
	RelRect FractalRect;
	Rect TextureRect;
};

struct Coordinate {
	HPReal X, Y;
};

template <typename T>
class ReservedMemoryArray {
	size_t ReservedSize = 0, CommittedSize = 0;
	size_t MaxCapacity = 0, Capacity = 0;
	size_t Size = 0;
	T *Data = nullptr;
	static constexpr size_t CommitSize = 256 * 1024;

	void Commit() {
		if (CommittedSize + CommitSize > ReservedSize) throw "";
		VirtualAlloc(((uint8_t *)Data) + CommittedSize, CommitSize, MEM_COMMIT, PAGE_READWRITE);
		CommittedSize += CommitSize;
		Capacity = CommittedSize / sizeof(T);
	}

public:

	ReservedMemoryArray() {}
	ReservedMemoryArray(size_t ReserveSize) : ReservedSize(ReserveSize), MaxCapacity(ReserveSize / sizeof(T)), Data((T *)VirtualAlloc(nullptr, ReserveSize, MEM_RESERVE, PAGE_READWRITE)) {}
	ReservedMemoryArray(const ReservedMemoryArray &a) = delete;
	ReservedMemoryArray(ReservedMemoryArray &&a) : ReservedSize(a.ReservedSize), CommittedSize(a.CommittedSize), MaxCapacity(a.MaxCapacity), Capacity(a.Capacity), Size(a.Size), Data(a.Data) {
		a.Data = nullptr;
	}
	~ReservedMemoryArray() {
		if (Data) {
			if constexpr (!std::is_trivially_destructible_v<T>) {
				for (size_t i = 0; i < Size; i++) {
					Data[i].~T();
				}
			}
			VirtualFree(Data, 0, MEM_RELEASE);
		}
	}

	ReservedMemoryArray &operator=(const ReservedMemoryArray &a) = delete;
	ReservedMemoryArray &operator=(ReservedMemoryArray &&a) {
		ReservedSize	= a.ReservedSize;
		CommittedSize	= a.CommittedSize;
		MaxCapacity		= a.MaxCapacity;
		Capacity		= a.Capacity;
		Size			= a.Size;
		Data			= a.Data;
		a.Data = nullptr;
		return *this;
	}

	size_t size() const { return Size; }

	size_t GetReservedSize() const { return ReservedSize; }

	void clear() { Size = 0; }

	void SetSize(size_t NewSize) {
		if (NewSize > Capacity) {
			size_t NewCommittedSize = (NewSize * sizeof(T) + CommitSize - 1) & ~(CommitSize - 1); // i * sizeof(complex) round up to nearest multiple of CommitSize
			if (NewCommittedSize > ReservedSize) throw "";
			VirtualAlloc(((uint8_t *)Data) + CommittedSize, NewCommittedSize - CommittedSize, MEM_COMMIT, PAGE_READWRITE);
			CommittedSize = NewCommittedSize;
			Capacity = CommittedSize / sizeof(T);
		}
		Size = NewSize;
	}

	void RecalculateCapacity() {
		MaxCapacity = ReservedSize / sizeof(T);
		Capacity = CommittedSize / sizeof(T);
	}

	void Append(const T &a) {
		if (Size >= Capacity) {
			Commit();
		}
		new (Data + Size) T(a);
		Size++;
	}

	void Pop() {
		if (Size) Size--;
	}

	void Append(T &&a) {
		if (Size >= Capacity) {
			Commit();
		}
		new (Data + Size) T(std::move(a));
		Size++;
	}
	void Concatenate(ReservedMemoryArray &&a) {
		static_assert(std::is_trivially_move_constructible_v<T>, "");
		static_assert(std::is_trivially_destructible_v<T>, "");
		size_t size = Size * sizeof(T);
		size_t ASize = a.Size * sizeof(T);
		uint8_t *data = (uint8_t *)Data;
		uint8_t *AData = (uint8_t *)a.Data;
		size_t Offset = 0;
		if (size < CommittedSize) {
			size_t CopySize = std::min(CommittedSize - size, ASize - Offset);
			memcpy(data + size, AData + Offset, CopySize);
			size += CopySize;
			Offset += CopySize;
		}
		while (Offset < ASize) {
			Commit();
			size_t CopySize = std::min(CommitSize, ASize - Offset);
			memcpy(data + size, AData + Offset, CopySize);
			VirtualFree(((uint8_t *)AData) + (Offset & ~(CommitSize - 1)), CommitSize, MEM_DECOMMIT);
			size += CopySize;
			Offset += CopySize;
		}
		Size += a.Size;
		VirtualFree(AData, 0, MEM_RELEASE);
		a = ReservedMemoryArray();
	}
	void ShrinkToFit() {
		if (Size + CommitSize / sizeof(T) * 2 < Capacity) {
			size_t NewCommittedSize = (Size * sizeof(T) + CommitSize - 1) & ~(CommitSize - 1); // i * sizeof(complex) round up to nearest multiple of CommitSize
			VirtualFree(((uint8_t *)Data) + NewCommittedSize, CommittedSize - NewCommittedSize, MEM_DECOMMIT);
			CommittedSize = NewCommittedSize;
			Capacity = CommittedSize / sizeof(T);
		}
	}
	operator T *() { return Data; }
	operator const T *() const { return Data; }

	void Serialize(std::ostream &stream) const {
		stream.write((const char *)&Size, sizeof(Size));
		stream.write((const char *)Data, Size * sizeof(T));
	}

	void Deserialize(std::istream &stream) {
		size_t NewSize;
		stream.read((char *)&NewSize, sizeof(Size));
		SetSize(NewSize);
		stream.read((char *)Data, Size * sizeof(T));
	}
};

struct RGBA {
	float R, G, B, A;
};

template <typename To, typename From>
inline To convert(const From &x) {
	return To(x);
}

template <>
inline double convert<double, mpf_class>(const mpf_class &x) {
	return x.get_d();
}

template <>
inline SRComplex convert<SRComplex, HRComplex>(const HRComplex &x) {
	return SRComplex(SRReal(x.real()), SRReal(x.imag()));
}

inline complexvec4<double> operator*(const complexvec4<double> &a, const std::complex<double> &b) {
	return a * complexvec4<double>(b);
}

inline complexvec4<double> operator*(const std::complex<double> &a, const complexvec4<double> &b) {
	return complexvec4<double>(a) * b;
}

inline complexvec4<FExpDouble> operator*(const complexvec4<FExpDouble> &a, const std::complex<FExpDouble> &b) {
	return a * complexvec4<FExpDouble>(b);
}

inline complexvec4<FExpDouble> operator*(const std::complex<FExpDouble> &a, const complexvec4<FExpDouble> &b) {
	return complexvec4<FExpDouble>(a) * b;
}

struct EvaluationParameters {
	IntIter MaxIterations;

	Coordinate CenterCoordinate;
	RelLocation Location;
};
#endif