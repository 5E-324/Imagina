#pragma once

#include <vector>
#include <assert.h>

inline void pause() {
	_mm_pause();
}

class SpinBarrier {
	std::atomic<ptrdiff_t> value;
	std::atomic<ptrdiff_t> initValue;

	static constexpr ptrdiff_t PhaseMask = 1;
	static constexpr ptrdiff_t RemainingThreadCountMask = ~PhaseMask;

public:
	explicit SpinBarrier(size_t threadCount) {
		value = initValue = threadCount << 1;
	}

	void Arrive() {
		ptrdiff_t newValue = (value -= 2);
		if (newValue <= 1) {
			assert(newValue >= 0);

			value = initValue | (newValue ^ 1);
		}
	}
	size_t ArriveAndWait() {
		ptrdiff_t newValue = (value -= 2);
		size_t newPhase = (~newValue) & PhaseMask;
		if (newValue <= 1) {
			assert(newValue >= 0);

			value = initValue | newPhase;
		} else {
			while ((value & PhaseMask) != newPhase); // pause();
		}

		return newPhase;
	}
	void ArriveAndDrop() {
		initValue -= 2;
		Arrive();
	}
	size_t WaitWithoutArrive() {
		size_t value2;
		while ((value2 = value) > 3) pause();
		return (~value2) & PhaseMask;
	}
	size_t Phase() const { return value & PhaseMask; }
	size_t NextPhase() const { return (~value) & PhaseMask; }
};

class HInfLAEvaluator::FeatureFinder : public Evaluator::FeatureFinder {
	struct StepLimitReached {};
	static constexpr IntIter StepLimit = 1048576;

	static constexpr SRReal NearLinearRadiusScale = 0.25;
	static constexpr SRReal SqrNearLinearRadiusScale = NearLinearRadiusScale * NearLinearRadiusScale;
	static constexpr SRReal MisiurewiczCaptureRadiusScale = 0.5;
	static constexpr SRReal SqrMisiurewiczCaptureRadiusScale = MisiurewiczCaptureRadiusScale * MisiurewiczCaptureRadiusScale;

	static constexpr SRReal Tolerance = std::numeric_limits<float>::epsilon();
	static constexpr SRReal SqrTolerance = Tolerance * Tolerance;

	class MandelbrotFeature : public Feature {
		enum class FeatureType {
			PeriodicPoint,
			MisiurewiczPoint,
			RealAxis,
		};
		HRReal X, Y, Scale;
		FeatureType type;
		IntIter Preperiod = 0, Period = 0;
		std::wstring information;

		MandelbrotFeature() = default;
		MandelbrotFeature(const HRReal &X, const HRReal &Y, FeatureType type) : X(X), Y(Y), type(type) {}

	public:
		virtual std::wstring_view Name() override;
		virtual std::wstring_view Information() override;
		virtual IntIter ItLimForZoomLevel(HRReal HalfH) override;
		virtual void GetCoordinate(HRReal &X, HRReal &Y) override;
		virtual bool GetRadius(HRReal &targetRadius) override;
		virtual bool CanLocatePrecisely() override;

		friend class HInfLAEvaluator::FeatureFinder;
	};

	template<typename real = HRReal>
	struct EvaluationContext {
		using complex = std::complex<real>;

		const LAReference<real> &reference;

		struct Point {
			complex z, dzdc;
			real Magnitude;
			IntIter Iteration;
		};

		std::vector<Point> MisiurewiczCheckingStack;

		IntIter Iteration, RefIteration, MaxIteration, PrePeriod, Period;

		HRReal InvScalingFactor;
		real ScalingFactor;

		complex c, dc;
		complex z, dz, zcoeff, dzdc;
		//complex PrevZ, PrevDzdc;
		complex ZAfterPrePeriod, DzdcAfterPrePeriod;
		bool PrePeriodReached;

		//real SqrRadius;
		real SqrMisiurewiczCaptureRadius;
		real SqrNearLinearRadius;

		EvaluationContext(const LAReference<real> &ref) : reference(ref) {}

		EvaluationContext(const LAReference<real> &ref, complex dc, complex RefC, real sqrRadius, uint64_t MaxIter, IntIter PrePeriod, IntIter Period) :
			reference(ref),
			Iteration(0), RefIteration(0), MaxIteration(MaxIter), PrePeriod(PrePeriod), Period(Period),
			InvScalingFactor(1.0_hr),
			ScalingFactor(real(1.0)),
			c(dc + RefC), dc(dc),
			z(real(0.0)), dz(real(0.0)), zcoeff(real(0.0)), dzdc(real(0.0)),
			PrePeriodReached(PrePeriod == 0),
			SqrMisiurewiczCaptureRadius(sqrRadius *SqrMisiurewiczCaptureRadiusScale),
			SqrNearLinearRadius(sqrRadius) {
		}

		template<typename T>
		EvaluationContext &operator=(const EvaluationContext<T> &x) {
			Iteration = x.Iteration; RefIteration = x.RefIteration; MaxIteration = x.MaxIteration; PrePeriod = x.PrePeriod; Period = x.Period;
			InvScalingFactor = std::max<HRReal>(ChebyshevNorm(x.dzdc), 1.0_hr) * 0x1.0p448_hr;
			HRReal SqrInvScalingFactor = InvScalingFactor * InvScalingFactor;
			HRReal HRScalingFactor = 1.0_hr / InvScalingFactor;
			ScalingFactor = SRReal(HRScalingFactor);

			c = convert<complex>(x.c); dc = convert<complex>(x.dc);
			z = convert<complex>(x.z); dz = convert<complex>(x.dz); zcoeff = convert<complex>(x.zcoeff * HRScalingFactor); dzdc = convert<complex>(x.dzdc * HRScalingFactor);
			ZAfterPrePeriod = convert<complex>(x.ZAfterPrePeriod); DzdcAfterPrePeriod = convert<complex>(x.DzdcAfterPrePeriod * HRScalingFactor);
			SqrMisiurewiczCaptureRadius = real(x.SqrMisiurewiczCaptureRadius * SqrInvScalingFactor);
			SqrNearLinearRadius = real(x.SqrNearLinearRadius * SqrInvScalingFactor);

			return *this;
		}

		template<typename T>
		EvaluationContext(const EvaluationContext<T> &x) : reference(x.reference) {
			*this = x;
		}

		bool FindMisiurewicz();
		bool CheckPeriodicity();
		template<bool FindPeriod> bool EvaluateDirectly();
		template<bool FindPeriod> bool EvaluatePT();
		template<bool FindPeriod> bool EvaluateLAStage(size_t LAStage);
	};

	template<bool FindPeriod>
	bool Evaluate(const LAReference<> &Ref, std::complex<HRReal> &diff, std::complex<HRReal> &zcoeff, std::complex<HRReal> &dzdc, std::complex<HRReal> dc, HRReal SqrRadius, uint64_t MaxIt, IntIter &PrePeriod, IntIter &Period);

	class MandelbrotFeatureLocatingTask : public PreciseLocatingTask, public CooperativeTask, public Task::Cancellable, public ProgressTrackable {
		MandelbrotFeature &Feature;
		size_t Precision;
		MandelbrotFeatureLocatingTask(MandelbrotFeature *feature, size_t Precision, const Coordinate &centerCoordinate);
		std::string ProgressString;
		size_t NRIteration = 0;
		size_t PredictedIterationCount = 0;
		size_t Iteration = 0;
		volatile bool Cancelled = false;
		volatile bool ShouldStop = false;
		SpinBarrier barrier{ 3 };
		HPReal ZR[2];
		HPReal ZI[2];
		HPReal CR, CI;

	public:
		virtual std::string_view GetDescription() const override;
		virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const override;
		virtual std::string_view GetDetailedProgress() override;
		virtual void Execute(size_t ThreadID) override;
		virtual void Cancel() override;

		friend class FeatureFinder;
	};

public:
	bool FindPeriodicPoint(MandelbrotFeature &Feature, HRReal x, HRReal y, HRReal radius, LAReference<> &reference);
	virtual Feature *FindFeature(HRReal x, HRReal y, HRReal radius, Reference *reference) override;
	virtual PreciseLocatingTask *CreatePreciseLocatingTask(Feature *feature, size_t Precision, const Coordinate &centerCoordinate) override;
};