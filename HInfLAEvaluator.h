#pragma once

#include "Evaluator.h"
#include "Computation.h"
#include <complex>

#define PERIOD_DETECTION_METHOD2

template<typename T> inline T ChebyshevNorm(std::complex<T> z) {
	return std::max(std::abs(z.real()), std::abs(z.imag()));
}

template<> inline DExpVec4 ChebyshevNorm<DExpVec4>(std::complex<DExpVec4> z) {
	return DExpVec4::absmax(z.real(), z.imag());
}

template<> inline dvec4 ChebyshevNorm<dvec4>(std::complex<dvec4> z) {
	return max(z.real().abs(), z.imag().abs());
}

struct ATInfo {
	size_t StepLength;
	HRReal ThresholdC;
	SRReal SqrEscapeRadius;
	std::complex<SRReal> RefC;
	std::complex<HRReal> ZCoeff, CCoeff, InvZCoeff;

	bool Usable(HRReal SqrRadius) {
		return SqrRadius * norm(CCoeff) * 0x1.0p32_hr > norm(std::complex<HRReal>(RefC)) && SqrEscapeRadius > 4.0;
	}
};

#ifdef PERIOD_DETECTION_METHOD2
#define IF_PERIOD_DETECTION_METHOD2(...) __VA_ARGS__
#else
#define IF_PERIOD_DETECTION_METHOD2(...)
#endif

template<typename Real = HRReal>
struct LAInfo {
	static constexpr Real Stage0PeriodDetectionThreshold = 0x1.0p-10;
	static constexpr Real PeriodDetectionThreshold = 0x1.0p-10;
	static constexpr Real Stage0PeriodDetectionThreshold2 = 0x1.0p-6;
	static constexpr Real PeriodDetectionThreshold2 = 0x1.0p-3;
	static constexpr Real LAThresholdScale = 0x1.0p-24;
	static constexpr Real LAThresholdCScale = 0x1.0p-24;
	static constexpr size_t Order = 4;

	using complex = std::complex<Real>;
	using vreal = vec4<Real>;
	using vcomplex = std::complex<vreal>;
	complex Ref, ZCoeff, CCoeff;
	Real LAThreshold, LAThresholdC;
	IF_PERIOD_DETECTION_METHOD2(Real MinMag);

	LAInfo() {}
	explicit LAInfo(complex z) :
		Ref(z), ZCoeff(1.0), CCoeff(1.0), LAThreshold(1.0), LAThresholdC(1.0) IF_PERIOD_DETECTION_METHOD2(, MinMag(4.0)) {}

	template <typename T> LAInfo(const LAInfo<T> &LA) :
		Ref(Real(LA.Ref.real()), Real(LA.Ref.imag())),
		ZCoeff(Real(LA.ZCoeff.real()), Real(LA.ZCoeff.imag())),
		CCoeff(Real(LA.CCoeff.real()), Real(LA.CCoeff.imag())),
		LAThreshold(LA.LAThreshold),
		LAThresholdC(LA.LAThresholdC)
		IF_PERIOD_DETECTION_METHOD2(, MinMag(LA.MinMag)) {}

	bool Stage0DetectPeriod(complex z) {
#ifdef PERIOD_DETECTION_METHOD2
		return ChebyshevNorm(z) < MinMag * Stage0PeriodDetectionThreshold2;
#else
		return ChebyshevNorm(z) / ChebyshevNorm(ZCoeff) * LAThresholdScale < LAThreshold *Stage0PeriodDetectionThreshold;
#endif
	}

	bool DetectPeriod(complex z) {
#ifdef PERIOD_DETECTION_METHOD2
		return ChebyshevNorm(z) < MinMag * PeriodDetectionThreshold2;
#else
		return ChebyshevNorm(z) / ChebyshevNorm(ZCoeff) * LAThresholdScale < LAThreshold *PeriodDetectionThreshold;
#endif
	}

	bool Step(LAInfo &out, complex z) const {
		Real ChebyMagz = ChebyshevNorm(z);
		Real ChebyMagZCoeff = ChebyshevNorm(ZCoeff);
		Real ChebyMagCCoeff = ChebyshevNorm(CCoeff);

#ifdef PERIOD_DETECTION_METHOD2
		out.MinMag = std::min(ChebyMagz, MinMag);
#endif

		out.LAThreshold = std::min(LAThreshold, ChebyMagz / ChebyMagZCoeff * LAThresholdScale);
		out.LAThresholdC = std::min(LAThresholdC, ChebyMagz / ChebyMagCCoeff * LAThresholdCScale);
		out.ZCoeff = (Real)2.0 * z * ZCoeff;
		out.CCoeff = (Real)2.0 * z * CCoeff + (Real)1.0;
		out.Ref = Ref;

#ifdef PERIOD_DETECTION_METHOD2
		return out.MinMag < MinMag *Stage0PeriodDetectionThreshold2;
#else
		return out.LAThreshold < LAThreshold *Stage0PeriodDetectionThreshold;
#endif
	}
	LAInfo Step(complex z) const {
		LAInfo Result;

		Step(Result, z);
		return Result;
	}
	bool Composite(LAInfo &out, LAInfo LA) const {
		complex z = LA.Ref;
		Real ChebyMagz = ChebyshevNorm(z);
		Real ChebyMagZCoeff = ChebyshevNorm(ZCoeff);
		Real ChebyMagCCoeff = ChebyshevNorm(CCoeff);

		out.LAThreshold = std::min(LAThreshold, ChebyMagz / ChebyMagZCoeff * LAThresholdScale);
		out.LAThresholdC = std::min(LAThresholdC, ChebyMagz / ChebyMagCCoeff * LAThresholdCScale);
		out.ZCoeff = (Real)2.0 * z * ZCoeff;
		Real RescaleFactor = out.LAThreshold / LAThreshold;
		out.CCoeff = (Real)2.0 * z * CCoeff;
		ChebyMagZCoeff = ChebyshevNorm(out.ZCoeff);
		ChebyMagCCoeff = ChebyshevNorm(out.CCoeff);
		Real temp = out.LAThreshold;
		out.LAThreshold = std::min(out.LAThreshold, LA.LAThreshold / ChebyMagZCoeff);
		out.LAThresholdC = std::min(out.LAThresholdC, LA.LAThreshold / ChebyMagCCoeff);
		out.ZCoeff = out.ZCoeff * LA.ZCoeff;
		RescaleFactor = out.LAThreshold / temp;
		out.CCoeff = out.CCoeff * LA.ZCoeff + LA.CCoeff;
		out.Ref = Ref;

#ifdef PERIOD_DETECTION_METHOD2
		temp = std::min(ChebyMagz, MinMag);
		out.MinMag = std::min(temp, LA.MinMag);

		return temp < MinMag *PeriodDetectionThreshold2;
#else
		return temp < LAThreshold *PeriodDetectionThreshold;
#endif
	}

	LAInfo Composite(LAInfo LA) const {
		LAInfo Result;

		Composite(Result, LA);
		return Result;
	}

	template<typename real>
	struct Temp {
		bool unusable;
		std::complex<real> newdz;
	};

	template<>
	struct Temp<vreal> {
		mask64x4 unusable;
		vcomplex newdz;
	};

	template<typename real>
	Temp<real> Prepare(const std::complex<real> &dz, const std::complex<real> &dc) const {
		Temp<real> temp;
		temp.newdz = dz * ((real(2.0) * std::complex<real>(Ref)) + dz);
		vreal ChebyMagNewdz = ChebyshevNorm(temp.newdz);
		temp.unusable = (ChebyMagNewdz >= real(LAThreshold));
		return temp;
	}

	template<typename real>
	void Evaluate(const Temp<real> &temp, std::complex<real> &dz, const std::complex<real> &dc) const {
		using complex = std::complex<real>;
		dz = temp.newdz * complex(ZCoeff) + dc * complex(CCoeff);
	}
	template<typename real>
	void EvaluateDzdz(const Temp<real> &temp, std::complex<real> &dz, std::complex<real> &dzdz, const std::complex<real> &dc, const real &ScalingFactor) const {
		using complex = std::complex<real>;

		dzdz = dzdz * real(2.0) * (dz + complex(Ref)) * complex(ZCoeff);
	}
	template<typename real>
	void EvaluateDzdc(const Temp<real> &temp, std::complex<real> &dz, std::complex<real> &dzdc, const std::complex<real> &dc, const real &ScalingFactor) const {
		using complex = std::complex<real>;

		dzdc = dzdc * real(2.0) * (dz + complex(Ref)) * complex(ZCoeff) + complex(CCoeff) * ScalingFactor;
	}
	template<typename real>
	void EvaluateDerivatives(const Temp<real> &temp, std::complex<real> &dz, std::complex<real> &dzdz, std::complex<real> &dzdc, const std::complex<real> &dc, const real &ScalingFactor) const {
		using complex = std::complex<real>;

		complex DzdzStep = real(2.0) * (dz + complex(Ref)) * complex(ZCoeff);

		dzdz = dzdz * DzdzStep;
		dzdc = dzdc * DzdzStep + complex(CCoeff) * ScalingFactor;
	}

	ATInfo CreateAT(const LAInfo &Next) {
		ATInfo Result;
		using HRComplex = std::complex<HRReal>;

		Result.ZCoeff = ZCoeff;
		Result.CCoeff = HRComplex(ZCoeff) * HRComplex(CCoeff);
		Result.InvZCoeff = 1.0_hr / HRComplex(Result.ZCoeff);
		Result.RefC = convert<SRComplex>(HRComplex(Next.Ref) * HRComplex(ZCoeff));
		Result.SqrEscapeRadius = SRReal(std::min(norm(HRComplex(ZCoeff)) * LAThreshold, 0x1.0p256_hr));
		Result.ThresholdC = std::min(HRReal(LAThresholdC), 0x1.0p256_hr / ChebyshevNorm(Result.CCoeff));

		return Result;
	}
	ATInfo CreateAT() {
		return CreateAT(*(this + 1));
	}
};

class HInfLAEvaluator : public Evaluator {
	static constexpr size_t MaxLAStages = 256;
	static constexpr size_t VSize = VectorSize;
	static constexpr size_t VSizeX = VectorSize2D.x;
	static constexpr size_t VSizeY = VectorSize2D.y;

	size_t RefReserveSize;
	size_t LAReserveSize;
	static constexpr size_t CommitSize = 256 * 1024;
	using Complex = std::complex<HRReal>;
	struct LAInfoI {
		size_t StepLength, NextStageLAIndex;
	};

	struct LAStageInfo {
		uint64_t LAIndex;
		uint64_t MacroItCount : 63;
		uint64_t UseDoublePrecision : 1;
	};

	struct JuliaPixel {
		double LogInvDE;
		double Iteration;
		double Weight = 1.0;

		JuliaPixel operator*(SRReal x) {
			return JuliaPixel{ LogInvDE * x, Iteration * x, Weight };
		}
		JuliaPixel operator+(JuliaPixel x) {
			return JuliaPixel{ LogInvDE + x.LogInvDE, Iteration + x.Iteration, Weight + x.Weight };
		}
	};

	static JuliaPixel Interpolate(SRReal x, JuliaPixel Pixel0, JuliaPixel Pixel1) {
		return Pixel0 * (1.0 - x) + Pixel1 * x;
	}

	struct JuliaMipmap {
		JuliaPixel *Data;
		IntPix Width, Height;
		IntPix OriginOffset;
		SRReal PixelDensity, PixelSize;
		SRReal FirstPixelWidth = 1.0, InvFirstPixelsCenterDistance;
		SRReal LastPixelWidth = 1.0, InvLastPixelsCenterDistance;

		JuliaPixel FetchTexel(IntPix TexelX, IntPix TexelY) {
			TexelY &= Height - 1; // Works only if Height is power of two

			TexelX = std::clamp<IntPix>(TexelX, 0, Width - 1);

			JuliaPixel &Pixel = Data[TexelX * Height + TexelY];

			return Pixel;
		}

		JuliaPixel Sample(SRReal x, SRReal y) {
			SRReal SampleX = -x * PixelDensity - 0.5;
			SRReal SampleY = y * PixelDensity - 0.5;

			SampleX += OriginOffset;
			SampleX = std::max(SampleX, 1.0);

			IntPix TexelX = floor(SampleX);
			IntPix TexelY = floor(SampleY);

			TexelX = std::clamp(TexelX, 0, Width - 2);

			SRReal FractX = SampleX - TexelX;
			SRReal FractY = SampleY - TexelY;

			if (TexelX == 0) {
				FractX = 1.0 - FractX;
				FractX *= InvFirstPixelsCenterDistance;
				FractX = 1.0 - FractX;
			} else if (TexelX == Width - 2) {
				FractX *= InvLastPixelsCenterDistance;
			}

			JuliaPixel Pixel1 = Interpolate(FractY, FetchTexel(TexelX, TexelY), FetchTexel(TexelX, TexelY + 1));
			TexelX++;
			JuliaPixel Pixel2 = FetchTexel(TexelX, TexelY) * (1 - FractY) + FetchTexel(TexelX, TexelY + 1) * FractY;

			return Interpolate(FractX, Pixel1, Pixel2);
		}
	};

	struct JuliaMipmaps {
		Int MipmapLevel;
		// Stored in column major order
		JuliaMipmap *Maps = nullptr;
		HRReal MinMagZ, MaxMagZ;
		SRReal PixelDensity, PixelSize;
		bool Valid = false;

		JuliaPixel Sample(HRComplex z, HRReal PixelSize);
	};

	struct LAReferenceTrivialContent {
		std::complex<double> Refc;
		size_t RefIt;
		size_t MaxIt;
		bool DoublePrecisionPT;
		bool DirectEvaluate;
		bool IsPeriodic;
		bool UseAT;


		ATInfo AT;

		size_t LAStageCount;
	};

	static_assert(std::is_trivially_copyable_v<LAReferenceTrivialContent>, "");

	template<typename Real = HRReal>
	struct LAReference : Reference, LAReferenceTrivialContent {
		using complex = std::complex<Real>;
		LAStageInfo LAStages[MaxLAStages];
		ReservedMemoryArray<complex> Ref;
		ReservedMemoryArray<LAInfo<Real>> LA;
		ReservedMemoryArray<LAInfoI> LAI;

		ReservedMemoryArray<std::complex<HRReal>> ExtendedRef;
		ReservedMemoryArray<uint64_t> ExtendedIterations;

		HRReal Radius;
		JuliaMipmaps Julia;

		LAReference(size_t RefReserveSize, size_t LAReserveSize) : Ref(RefReserveSize), LA(LAReserveSize), LAI(LAReserveSize) {}

		void ShrinkToFit() {
			Ref.ShrinkToFit();
			LA.ShrinkToFit();
			LAI.ShrinkToFit();
		}
	};

	LAReference<> *ref = nullptr;

	struct Cancellation {};

	template<typename vreal, bool Dzdc, bool Dzdz, bool Scaled>
	struct Derivatives {
		static_assert(!Dzdc && !Dzdz);
		Derivatives() {}

		template <typename vreal2, bool Scaled2>
		Derivatives(const Derivatives<vreal2, false, false, Scaled2> &derivatives) {}
	};

	template<typename vreal, bool Scaled>
	struct Derivatives<vreal, true, false, Scaled> {
		vec4<HRReal> InvScalingFactor, HRScalingFactor;
		vreal ScalingFactor;

		std::complex<vreal> dzdc;

		Derivatives() :
			InvScalingFactor(1.0), HRScalingFactor(1.0),
			ScalingFactor(1.0),
			dzdc(vreal(0.0)) {
		}

		template <typename vreal2>
		Derivatives(const Derivatives<vreal2, true, false, false> &derivatives, vec4<HRReal> InvScalingFactor) :
			InvScalingFactor(InvScalingFactor),
			HRScalingFactor(vec4<HRReal>(1.0) / InvScalingFactor),
			ScalingFactor(HRScalingFactor),
			dzdc(derivatives.dzdc *HRScalingFactor) {
		}

		template <typename vreal2>
		Derivatives(const Derivatives<vreal2, true, false, false> &derivatives) :
			Derivatives(derivatives, max(ChebyshevNorm(derivatives.dzdc), vreal2(1.0)) *vreal2(0x1.0p448)) {
		}
	};

	template<typename vreal, bool Scaled>
	struct Derivatives<vreal, true, true, Scaled> : Derivatives<vreal, true, false, Scaled> {
		using base = Derivatives<vreal, true, false, Scaled>;
		std::complex<vreal> dzdz;

		Derivatives() :
			base(),
			dzdz(vreal(1.0)) {
		}

		template <typename vreal2>
		Derivatives(const Derivatives<vreal2, true, true, false> &derivatives) :
			base(derivatives, max(ChebyshevNorm(derivatives.dzdz), vreal2(1.0)) *vreal2(0x1.0p448)),
			dzdz(derivatives.dzdz *base::HRScalingFactor) {
		}
	};

	template<typename real, bool EvaluateDzdc, bool EvaluateDzdz = false>
	struct EvaluationContext : Derivatives<vec4<real>, EvaluateDzdc, EvaluateDzdz, !std::is_same_v<real, HRReal>> {
		using derivatives = Derivatives<vec4<real>, EvaluateDzdc, EvaluateDzdz, !std::is_same_v<real, HRReal>>;
		using vreal = vec4<real>;
		using HRRealvec = vec4<HRReal>;
		using vcomplex = std::complex<vreal>;
		using ivec = i64vec4;
		using mask = mask64x4;

		volatile bool &Cancelled;

		LAReference<real> &Ref;
		size_t CurrentLAStage;
		size_t MaxIt;
		mask ValidMask;
		mask SamplingEnabledMask;
		mask SamplingDisablingMask = false;
		mask SamplingMask;
		mask MixedSamplingMask;
		mask SamplingUpdatedMask;
		ivec Iteration, RefIteration;
		vcomplex z, c;
		vcomplex dz, dc;

		vreal PixScale;
		vreal SqrPixScale;
		real SqrJuliaPixelSize;
		real SqrJuliaHalfPixelSize;
		real SqrJuliaMinMagZ;
		real SqrJuliaMaxMagZ;

		ivec SamplingIteration;
		vcomplex SamplingZ;
		vcomplex SamplingDzdc;

		EvaluationContext(volatile bool &Cancelled, LAReference<real> &Ref, real PixelScale, mask ValidMask, vcomplex dz, vcomplex dc, uint64_t MaxIter) :
			derivatives(),
			Cancelled(Cancelled), Ref(Ref),
			CurrentLAStage(Ref.LAStageCount), MaxIt(MaxIter),
			ValidMask(ValidMask), SamplingEnabledMask(true), SamplingMask(false), MixedSamplingMask(false), SamplingUpdatedMask(false),
			Iteration(0), RefIteration(0),
			z(dz), c(dc + vcomplex(Ref.Refc)),
			dz(dz), dc(dc),
			PixScale(PixelScale), SqrPixScale(PixScale *PixScale),
			SqrJuliaPixelSize(Ref.Julia.PixelSize *Ref.Julia.PixelSize),
			SqrJuliaHalfPixelSize(SqrJuliaPixelSize * 0.25),
			SqrJuliaMinMagZ(std::max(FContext.EvalLocation.HalfH * 0x1.0p48, Ref.Julia.MinMagZ *Ref.Julia.MinMagZ)),
			SqrJuliaMaxMagZ(Ref.Julia.MaxMagZ *Ref.Julia.MaxMagZ) {
		}

		EvaluationContext(volatile bool &Cancelled, LAReference<real> &Ref, real PixelScale, mask ValidMask, vcomplex dc, uint64_t MaxIter) :
			EvaluationContext(Cancelled, Ref, PixelScale, ValidMask, vreal(0.0), dc, MaxIter) {
		}

		template <typename real2>
		EvaluationContext(const EvaluationContext<real2, EvaluateDzdc, EvaluateDzdz> &Context) :
			derivatives(Context),
			Cancelled(Context.Cancelled),
			Ref(reinterpret_cast<LAReference<real> &>(Context.Ref)),
			CurrentLAStage(Context.CurrentLAStage), MaxIt(Context.MaxIt),
			ValidMask(Context.ValidMask), SamplingEnabledMask(Context.SamplingEnabledMask), SamplingMask(Context.SamplingMask), MixedSamplingMask(Context.MixedSamplingMask), SamplingUpdatedMask(false),
			Iteration(Context.Iteration), RefIteration(Context.RefIteration),
			z(Context.z), c(Context.c),
			dz(Context.dz), dc(Context.dc),
			PixScale(derivatives::InvScalingFactor *Context.PixScale), SqrPixScale(PixScale *PixScale),
			SqrJuliaPixelSize(Ref.Julia.PixelSize *Ref.Julia.PixelSize),
			SqrJuliaHalfPixelSize(SqrJuliaPixelSize * 0.25),
			SqrJuliaMinMagZ(std::max(FContext.EvalLocation.HalfH * 0x1.0p48, Ref.Julia.MinMagZ *Ref.Julia.MinMagZ)),
			SqrJuliaMaxMagZ(Ref.Julia.MaxMagZ *Ref.Julia.MaxMagZ) {
		}

		void TrySampleJulia(mask &ActiveMask);

		void EvaluateLA();

		void EvaluatePT();

		void EvaluateDirectly();

		void EvaluateAT(const real &SqrEscapeRadius);
	};

	template<typename real>
	struct ReferenceGenerationContext {
		size_t MaxIt = 0, Iteration = 0;
		LAReference<real> &reference;
		const volatile bool &Cancelled;

		ReferenceGenerationContext(LAReference<real> &reference, const volatile bool &Cancelled) : reference(reference), Cancelled(Cancelled) {}

		void CalculateOrbit(const EvaluationParameters &parameters, HRReal radius);
		bool CreateLAFromOrbit(HRReal radius);
		bool CreateLAFromOrbit_MT(HRReal radius);
		bool CreateNewLAStage();
		void CreateATFromLA(HRReal SqrRadius);
		void ConvertStageToDouble(size_t Stage);
		void GenerateApproximationData(HRReal radius);
		void GenerateReference(const EvaluationParameters &parameters);
	};

	class HInfLAReferenceTask : public Evaluator::ReferenceTask, public Task::Cancellable, public ProgressTrackable {
		HInfLAEvaluator *evaluator;
		EvaluationParameters parameters;
		volatile bool Cancelled = false;
		ReferenceGenerationContext<HRReal> *ComputationContext = nullptr;
		bool UseHRReal;

	public:
		HInfLAReferenceTask(HInfLAEvaluator *evaluator, const EvaluationParameters &parameters, Reference *ref) :
			evaluator(evaluator), parameters(parameters) {
			reference = ref;
		}
		virtual std::string_view GetDescription() const override;
		virtual void Execute() override;
		virtual void Cancel() override;
		virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const override;
	};

	class HInfLAEvaluationTask : public ParallelTask, public Task::Cancellable, public ProgressTrackable {
		volatile bool Cancelled = false;
		EvaluationParameters parameters;
		PixelManager *rasterizer;
		Reference *reference;

		class JuliaRenderingContext {
			std::mutex Mutex;
			volatile bool Started = false, Finished = false, Failed = false;

			size_t SuperSample = 2;
			IntPix PixelGroupSize = (SuperSample == 1) ? 2 : 1;

			size_t SuperSampleSqr = SuperSample * SuperSample;
			SRReal SubPixelSize = 1.0 / SuperSample;

			IntPix Width;
			IntPix Height = 16; // Must be power of two
			IntPix OriginOffset;
			SRReal PixelSize = 3.14159265358979323846 / Height;

			size_t MipmapLevel = 32;

			std::atomic<IntPix> Column = 0;

			VSRReal<VectorSize> EvaluateJulia(const EvaluationParameters &parameters, LAReference<> &reference, VHRComplex<VectorSize> dc, VHRReal<VectorSize> &InvDE);

		public:
			void RenderJulia(const EvaluationParameters &parameters, LAReference<> &reference, bool IsMainThread);
		};

		JuliaRenderingContext JuliaContext;

		template <bool DistanceEstimation>
		void Evaluate(GroupedRasterizingInterface &RI, Reference *Reference);

	public:
		HInfLAEvaluationTask(const EvaluationParameters &parameters, PixelManager *rasterizer, Reference *reference) : parameters(parameters), rasterizer(rasterizer), reference(reference) {}
		virtual std::string_view GetDescription() const override;
		virtual void Execute(size_t ThreadID) override;
		virtual void Cancel() override;
		virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const override;
	};

public:
	HInfLAEvaluator();
	~HInfLAEvaluator() override;

	virtual ReferenceTask *CreateReferenceTask(const EvaluationParameters &parameters) override;
	virtual void FreeReference(Reference *reference) override;
	virtual Task *CreateEvaluationTask(const EvaluationParameters &parameters, PixelManager *rasterizer, Reference *reference) override;

	template<typename real, typename orbitreal = real> void SaveCompressedOrbit(std::ostream &stream, const LAReference<real> &Ref);
	template<typename real> void LoadCompressedOrbit(std::istream &stream, LAReference<real> &Ref);

	virtual size_t SaveReference(std::ostream &stream, const Reference *reference) override;
	virtual size_t LoadReference(std::istream &stream, Reference *&reference) override;

	class FeatureFinder;
	FeatureFinder *featureFinder;

	virtual Evaluator::FeatureFinder *GetFeatureFinder() override;

};

#include "FeatureFinder.h"