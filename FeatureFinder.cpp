#include "Includes.h"
#include "FractalContext.h"
#include <vector>
#include <sstream>
#include <barrier>

Evaluator::FeatureFinder *HInfLAEvaluator::GetFeatureFinder() {
	return featureFinder;
}

std::wstring_view HInfLAEvaluator::FeatureFinder::MandelbrotFeature::Name() {
	using namespace std;
	switch (type) {
		case FeatureType::PeriodicPoint: return L"Periodic point"sv;
		case FeatureType::MisiurewiczPoint: return L"Misiurewicz point"sv;
		case FeatureType::RealAxis: return L"Real axis"sv;
	}
	return std::wstring_view();
}

std::wstring_view HInfLAEvaluator::FeatureFinder::MandelbrotFeature::Information() {
	if (!information.empty()) return information;
	if (type == FeatureType::RealAxis) return std::wstring_view();
	std::wstringstream stream;
	switch (type) {
		case FeatureType::PeriodicPoint: {
			stream << L"Period: " << Period << std::endl;
			char Buffer[32];
			gmp_sprintf(Buffer, "%.7Fg", Scale.to_mpf_class(53).get_mpf_t());
			stream << L"Scale: " << Buffer;
			break;
		}
		case FeatureType::MisiurewiczPoint: {
			stream << L"Preperiod: ≤" << Preperiod << std::endl;
			stream << L"Period: " << Period << std::endl;
			break;
		}
	}
	information = stream.str();
	return information;
}

void HInfLAEvaluator::FeatureFinder::MandelbrotFeature::GetCoordinate(HRReal &X, HRReal &Y) {
	X = this->X;
	Y = this->Y;
}

IntIter HInfLAEvaluator::FeatureFinder::MandelbrotFeature::ItLimForZoomLevel(HRReal HalfH) {
	if (type != FeatureType::PeriodicPoint) return 0;
	return Period * 256;
}

bool HInfLAEvaluator::FeatureFinder::MandelbrotFeature::CanLocatePrecisely() {
	return type == FeatureType::PeriodicPoint;
}

bool HInfLAEvaluator::FeatureFinder::MandelbrotFeature::GetRadius(HRReal &targetRadius) {
	if (type != FeatureType::PeriodicPoint) return false;
	targetRadius = Scale * 4.0_hr;
	return true;
}

template<typename real>
inline bool HInfLAEvaluator::FeatureFinder::EvaluationContext<real>::CheckPeriodicity() {
	real Magnitude = norm(z);
	if (Magnitude < SqrNearLinearRadius * norm(dzdc)) {
		PrePeriod = 0;
		Period = Iteration;
		return true;
		//return Iteration;
	}
	if (!MisiurewiczCheckingStack.empty()) {
		if (MisiurewiczCheckingStack.back().Magnitude > Magnitude) {
			Point p;
			do {
				p = MisiurewiczCheckingStack.back();
				MisiurewiczCheckingStack.pop_back();
			} while (!MisiurewiczCheckingStack.empty() && MisiurewiczCheckingStack.back().Magnitude > Magnitude);
			if (norm(z - p.z) < SqrMisiurewiczCaptureRadius * norm(dzdc - p.dzdc) && norm(z - p.z) > 0x1.0p-24 * Magnitude) {
				ZAfterPrePeriod = p.z;
				DzdcAfterPrePeriod = p.dzdc;
				PrePeriod = p.Iteration;
				Period = Iteration;
				return true;
			}
		}
		if (!MisiurewiczCheckingStack.empty()) {
			Point p = MisiurewiczCheckingStack.back();
			if (norm(z - p.z) < SqrMisiurewiczCaptureRadius * norm(dzdc - p.dzdc) && norm(z - p.z) > 0x1.0p-24 * Magnitude) {
				ZAfterPrePeriod = p.z;
				DzdcAfterPrePeriod = p.dzdc;
				PrePeriod = p.Iteration;
				Period = Iteration;
				return true;
			}
		}
	}
	MisiurewiczCheckingStack.push_back({ z, dzdc, Magnitude, IntIter(Iteration) });
	if (Magnitude * SqrNearLinearRadiusScale < SqrNearLinearRadius * norm(dzdc)) {
		SqrNearLinearRadius = Magnitude * SqrNearLinearRadiusScale / norm(dzdc);
	}
	return false;
}

template <typename real> template<bool FindPeriod>
bool HInfLAEvaluator::FeatureFinder::EvaluationContext<real>::EvaluateDirectly() {
	uint64_t IterationLimit = std::min(Iteration + StepLimit, MaxIteration);
	while (Iteration < IterationLimit) {
		if (Iteration == 0) { // TODO: Optimize this
			zcoeff = ScalingFactor;
		} else {
			zcoeff = zcoeff * real(2.0) * z;
		}
		dzdc = dzdc * real(2.0) * z + ScalingFactor;
		z = c + z * z;
		Iteration++;

		if (norm(z) > real(4096.0)) break;

		if constexpr (FindPeriod) {
			bool PeriodFound = CheckPeriodicity();
			if (PeriodFound) return true;
		} else if (Iteration == PrePeriod) {
			PrePeriodReached = true;
			ZAfterPrePeriod = z;
			DzdcAfterPrePeriod = dzdc;
		}
	}
	if (Iteration >= IterationLimit) {
		if (IterationLimit < MaxIteration) {
			throw StepLimitReached();
		} else {
			return true;
		}
	}
	return false;
}

template <typename real> template<bool FindPeriod>
bool HInfLAEvaluator::FeatureFinder::EvaluationContext<real>::EvaluatePT() {
	uint64_t IterationLimit = std::min(Iteration + StepLimit, MaxIteration);
	while (Iteration < IterationLimit) {
		if (Iteration == 0) { // TODO: Optimize this
			zcoeff = ScalingFactor;
		} else {
			zcoeff = zcoeff * real(2.0) * z;
		}
		dzdc = dzdc * real(2.0) * z + ScalingFactor;
		dz = dz * (complex(reference.Ref[RefIteration]) + z) + dc;
		Iteration++;
		RefIteration++;
		z = dz + complex(reference.Ref[RefIteration]);
		if (RefIteration == reference.RefIt || norm(z) < norm(dz)) {
			dz = z;
			RefIteration = 0;
		}

		if (norm(z) > real(4096.0)) break;

		if constexpr (FindPeriod) {
			bool PeriodFound = CheckPeriodicity();
			if (PeriodFound) return true;
		} else if (Iteration == PrePeriod) {
			PrePeriodReached = true;
			ZAfterPrePeriod = z;
			DzdcAfterPrePeriod = dzdc;
		}
	}
	if (Iteration >= IterationLimit) {
		if (IterationLimit < MaxIteration) {
			throw StepLimitReached();
		} else {
			return true;
		}
	}
	return false;
}

template<typename real> template<bool FindPeriod>
bool HInfLAEvaluator::FeatureFinder::EvaluationContext<real>::EvaluateLAStage(size_t LAStage) {
	size_t Steps = 0;
	size_t LAIndex = reference.LAStages[LAStage].LAIndex;
	const LAInfo<real> *LA = reference.LA + LAIndex;
	const LAInfoI *LAI = reference.LAI + LAIndex;
	size_t MacroItCount = reference.LAStages[LAStage].MacroItCount;

	if (ChebyshevNorm(dc) >= LA[0].LAThresholdC) {
		RefIteration = LAI[RefIteration].NextStageLAIndex;
		return 0;
	}

	while (Iteration < MaxIteration) {
		Steps++;
		if (Steps > StepLimit) throw StepLimitReached();
		auto temp = LA[RefIteration].Prepare(dz, dc);
		if (temp.unusable) {
			RefIteration = LAI[RefIteration].NextStageLAIndex;
			break;
		}
		if (Iteration == 0) { // TODO: Optimize this
			zcoeff = ScalingFactor * LA[RefIteration].ZCoeff;
		} else {
			LA[RefIteration].EvaluateDzdz(temp, dz, zcoeff, dc, ScalingFactor);
		}
		LA[RefIteration].EvaluateDzdc(temp, dz, dzdc, dc, ScalingFactor);
		LA[RefIteration].Evaluate(temp, dz, dc);
		Iteration += LAI[RefIteration].StepLength;

		RefIteration++;
		z = dz + complex(LA[RefIteration].Ref);

		if (RefIteration == MacroItCount || ChebyshevNorm(dz) > ChebyshevNorm(z)) {
			RefIteration = 0;
			dz = z;
		}

		if constexpr (FindPeriod) {
			bool PeriodFound = CheckPeriodicity();
			if (PeriodFound) return true;
		} else if (Iteration == PrePeriod) {
			PrePeriodReached = true;
			ZAfterPrePeriod = z;
			DzdcAfterPrePeriod = dzdc;
		}
	}
	if (Iteration >= MaxIteration) {
		return true;
	} else {
		return false;
	}
}

template<bool FindPeriod>
bool HInfLAEvaluator::FeatureFinder::Evaluate(const LAReference<> &Ref, std::complex<HRReal> &diff, std::complex<HRReal> &zcoeff, std::complex<HRReal> &dzdc, std::complex<HRReal> dc, HRReal SqrRadius, uint64_t MaxIt, IntIter &PrePeriod, IntIter &Period) {
	const LAReference<double> &DPRef = reinterpret_cast<const LAReference<double> &>(Ref);

	EvaluationContext<HRReal> HRContext(Ref, dc, Ref.Refc, SqrRadius, MaxIt, PrePeriod, Period);
	EvaluationContext<SRReal> SRContext(DPRef);

	size_t CurrentLAStage = Ref.LAStageCount;

	if (Ref.DirectEvaluate) {
		SRContext = HRContext;
		SRContext.EvaluateDirectly<FindPeriod>();

		goto SRExit;
	}

	while (CurrentLAStage--) {
		if (Ref.LAStages[CurrentLAStage].UseDoublePrecision) {
			CurrentLAStage++;
			SRContext = HRContext;
			goto SRLA;
		}
		bool Finished = HRContext.EvaluateLAStage<FindPeriod>(CurrentLAStage);
		if (Finished) goto HRExit;
	}

	if (Ref.DoublePrecisionPT) {
		SRContext = HRContext;
		goto SRPT;
	}
	HRContext.EvaluatePT<FindPeriod>();

HRExit:
	{
		bool Succeed = FindPeriod ? (HRContext.Period != 0) : (!PrePeriod || HRContext.PrePeriodReached) && HRContext.Iteration == Period;

		if (Succeed) {
			diff = HRContext.z;
			zcoeff = HRContext.zcoeff;
			dzdc = HRContext.dzdc;
			if (HRContext.PrePeriod != 0) {
				diff -= HRContext.ZAfterPrePeriod;
				dzdc -= HRContext.DzdcAfterPrePeriod;
			}
			Period = HRContext.Period;
			PrePeriod = HRContext.PrePeriod;
		}

		return Succeed;
	}
SRLA:
	while (CurrentLAStage--) {
		bool Finished = SRContext.EvaluateLAStage<FindPeriod>(CurrentLAStage);
		if (Finished) goto SRExit;
	}

SRPT:
	SRContext.EvaluatePT<FindPeriod>();

SRExit:
	{
		bool Succeed = FindPeriod ? (SRContext.Period != 0) : (!PrePeriod || SRContext.PrePeriodReached) && SRContext.Iteration == Period;

		if (Succeed) {
			diff = SRContext.z;
			zcoeff = SRContext.zcoeff;
			dzdc = SRContext.dzdc;
			if (SRContext.PrePeriod != 0) {
				diff -= SRContext.ZAfterPrePeriod;
				dzdc -= SRContext.DzdcAfterPrePeriod;
			}
		}
		zcoeff *= SRContext.InvScalingFactor;
		dzdc *= SRContext.InvScalingFactor;

		Period = SRContext.Period;
		PrePeriod = SRContext.PrePeriod;

		return Succeed;
	}
}

HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::
MandelbrotFeatureLocatingTask(MandelbrotFeature *feature, size_t Precision, const Coordinate &centerCoordinate) :
	CooperativeTask(3), Feature(*feature), Precision(Precision),
	ZR{ HPReal(0, Precision), HPReal(0, Precision) }, ZI{ HPReal(0, Precision), HPReal(0, Precision) },
	CR(centerCoordinate.X), CI(centerCoordinate.Y) {
	CR += Feature.X;
	CI += Feature.Y;
	CR.set_prec(Precision);
	CI.set_prec(Precision);
}

std::string_view HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::GetDescription() const {
	using namespace std;
	return "Locating periodic point"sv;
}

bool HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::GetProgress(SRReal &Numerator, SRReal &Denoninator) const {
	if (!PredictedIterationCount) return false;

	Numerator = (SRReal)NRIteration + (SRReal)Iteration / Feature.Period;
	Denoninator = PredictedIterationCount;

	return true;
}

std::string_view HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::GetDetailedProgress() {
	std::stringstream stream;

	stream << "Newton raphson iteration " << NRIteration + 1;
	if (PredictedIterationCount) {
		stream << " / " << PredictedIterationCount;
	}
	stream << std::endl;
	stream << std::fixed << std::setprecision(1) << ((double)Iteration / Feature.Period * 100) << "% (" << Iteration << "/" << Feature.Period << ")" << std::endl;
	ProgressString = stream.str();

	return ProgressString;
}

void HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::Execute(size_t ThreadID) {
	HRComplex z, dzdc, d2zdc2;

	size_t FeatureCoordinateExponent = std::max(Feature.X.Exponent, Feature.Y.Exponent);
	size_t DerivativePrecision = (-Feature.Scale.Exponent + 32) / 4;
	DerivativePrecision = std::min(DerivativePrecision, Precision + FeatureCoordinateExponent + 32);

	HPReal DZDCR = HPReal(0, DerivativePrecision);
	HPReal DZDCI = HPReal(0, DerivativePrecision);

	HPReal Temp = HPReal(0, Precision);
	HPReal Templp = HPReal(0, DerivativePrecision);
	HPReal Templp2 = HPReal(0, DerivativePrecision);

	switch (ThreadID) {
		case 1: {
			while (true) {
				size_t Phase = barrier.ArriveAndWait();
				if (ShouldStop) {
					return;
				}

				size_t NextPhase = Phase ^ 1;

				ZR[NextPhase] = ZR[Phase] + ZI[Phase];
				Temp = ZR[Phase] - ZI[Phase];

				ZR[NextPhase] = ZR[NextPhase] * Temp + CR;

			}
			return;
		}
		case 2: {
			while (true) {
				size_t Phase = barrier.ArriveAndWait();
				if (ShouldStop) {
					return;
				}

				size_t NextPhase = Phase ^ 1;

				ZI[NextPhase] = ZI[Phase] * ZR[Phase];
				ZI[NextPhase] += ZI[NextPhase];
				ZI[NextPhase] += CI;
			}
			return;
		}
	}

	for (NRIteration = 0; NRIteration < 32; NRIteration++) {
		size_t Phase = barrier.NextPhase();
		ZR[Phase] = CR;
		ZI[Phase] = CI;
		DZDCR = 1;
		DZDCI = 0;
		dzdc = 1.0_hr;
		d2zdc2 = 0.0_hr;

		for (Iteration = 1; Iteration < Feature.Period; Iteration++) {
			if (Cancelled) {
				ShouldStop = true;
				barrier.ArriveAndDrop();
				return;
			}

			Phase = barrier.ArriveAndWait();
			size_t NextPhase = Phase ^ 1;

			Templp = DZDCR * ZR[Phase];
			Templp2 = DZDCI * ZI[Phase];
			Templp -= Templp2;
			Templp += Templp;
			Templp += 1;
			DZDCI *= ZR[Phase];
			Templp2 = DZDCR * ZI[Phase];
			DZDCI += Templp2;
			DZDCI += DZDCI;
			DZDCR = Templp;

			z = { ZR[Phase], ZI[Phase] };
			d2zdc2 = 2.0_hr * (dzdc * dzdc + z * d2zdc2);
			dzdc *= 2.0_hr * z;
		}
		Phase = barrier.WaitWithoutArrive();

		// Z / DZDC
		Temp = ZR[Phase] * DZDCR + ZI[Phase] * DZDCI;
		ZI[Phase] = ZI[Phase] * DZDCR - ZR[Phase] * DZDCI;
		ZR[Phase] = Temp;
		Temp = DZDCR * DZDCR + DZDCI * DZDCI; // denominator
		ZR[Phase] /= Temp;
		ZI[Phase] /= Temp;

		CR -= ZR[Phase];
		CI -= ZI[Phase];

		HRReal diffr = ZR[Phase];
		HRReal diffi = ZI[Phase];
		HRReal normdiff = diffr * diffr + diffi * diffi;
		HRReal err = normdiff * normdiff * norm(d2zdc2) / norm(dzdc);
		if (-err.Exponent >= Precision * 2) break;
		SRReal newCorrectBits = log2(normdiff / err);
		int IterRemaining = ceil(log2((Precision * 2 + err.Exponent) / (newCorrectBits * 2) + 1));
		PredictedIterationCount = NRIteration + std::max(IterRemaining, 1) + 1;
	}

	ShouldStop = true;
	barrier.ArriveAndDrop();
	coordinate.X = std::move(CR);
	coordinate.Y = std::move(CI);
}

void HInfLAEvaluator::FeatureFinder::MandelbrotFeatureLocatingTask::Cancel() {
	Cancelled = true;
}

bool HInfLAEvaluator::FeatureFinder::FindPeriodicPoint(MandelbrotFeature &Feature, HRReal x, HRReal y, HRReal radius, LAReference<> &reference) {
	double Tolerance = 0x1p-40;
	double SqrTolerance = Tolerance * Tolerance;
	using complex = std::complex<HRReal>;

	HRReal SqrRadius = radius * radius;

	size_t MaxIt = reference.MaxIt;
	if (reference.IsPeriodic || Global::ItLim < reference.MaxIt) { // FIXME
		MaxIt = Global::ItLim;
	}

	complex InitialDc = { x, y };
	complex dc = InitialDc;
	complex diff;
	complex zcoeff;
	complex dzdc;

	IntIter PrePeriod = 0, Period = 0;
	bool PeriodFound = Evaluate<true>(reference, diff, zcoeff, dzdc, dc, SqrRadius, MaxIt, PrePeriod, Period);
	if (!PeriodFound) return false;
	Feature.Preperiod = PrePeriod;
	Feature.Period = Period - PrePeriod;

	dc -= diff / dzdc;
	Feature.X = dc.real();
	Feature.Y = dc.imag();

	for (size_t NRIteration = 0; NRIteration < 32; NRIteration++) {
		bool Succeed = Evaluate<false>(reference, diff, zcoeff, dzdc, dc, SqrRadius, Period, PrePeriod, Period);
		if (!Succeed) return false;
		diff /= dzdc;
		dc -= diff;
		if (norm(diff) < norm(dc) * SqrTolerance) break;
	}

	{
		bool Succeed = Evaluate<false>(reference, diff, zcoeff, dzdc, dc, SqrRadius, Period, PrePeriod, Period);
		if (!Succeed) return false;
		diff /= dzdc;
		dc -= diff;
	}

	if (norm(InitialDc - dc) > SqrRadius) return false;
	Feature.X = dc.real();
	Feature.Y = dc.imag();
	Feature.Scale = HRReal(1.0) / sqrt(norm(zcoeff * dzdc));
	return true;
}

Evaluator::Feature *HInfLAEvaluator::FeatureFinder::FindFeature(HRReal x, HRReal y, HRReal radius, Reference *reference) {
	LAReference<> &Ref = *static_cast<LAReference<>*>(reference);
	MandelbrotFeature *Feature = new MandelbrotFeature;
	using enum MandelbrotFeature::FeatureType;
	try {
		bool Found = FindPeriodicPoint(*Feature, x, y, radius, Ref);
		Feature->type = Feature->Preperiod ? MisiurewiczPoint : PeriodicPoint;
		if (Found) return Feature;
	} catch (StepLimitReached) {

	}

	if (std::abs(y + Ref.Refc.imag()) < radius) {
		Feature->X = x;
		Feature->Y = -Ref.Refc.imag();
		Feature->type = RealAxis;
		return Feature;
	}
	delete Feature;
	return nullptr;
}

Evaluator::FeatureFinder::PreciseLocatingTask *HInfLAEvaluator::FeatureFinder::CreatePreciseLocatingTask(Feature *feature, size_t Precision, const Coordinate &centerCoordinate) {
	MandelbrotFeature *mandelbrotFeature = dynamic_cast<MandelbrotFeature *>(feature);
	if (!mandelbrotFeature || mandelbrotFeature->type != MandelbrotFeature::FeatureType::PeriodicPoint) return nullptr;
	return new MandelbrotFeatureLocatingTask(mandelbrotFeature, Precision, centerCoordinate);
}