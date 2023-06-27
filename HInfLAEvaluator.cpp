#include "Includes.h"
#include "FractalContext.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <random>

#ifndef M_E
#define M_E 2.71828182845905
#endif

HInfLAEvaluator::JuliaPixel HInfLAEvaluator::JuliaMipmaps::Sample(HRComplex z, HRReal PixelSize) {
	SRComplex SamplePosition = log(z);

	HRReal HRPixelSizeRatio = (PixelSize) / (sqrt(std::norm(z)) * this->PixelSize);
	SRReal PixelSizeRatio = SRReal(HRPixelSizeRatio);
	SRReal MipmapLevel = log2(PixelSizeRatio);
	Int MipmapLevelI = floor(MipmapLevel);
	MipmapLevelI = std::clamp(MipmapLevelI, Int(0), this->MipmapLevel - 1);

	JuliaMipmap &Mipmap = Maps[MipmapLevelI];

	JuliaPixel Pixel = Mipmap.Sample(SamplePosition.real(), SamplePosition.imag());

	if (MipmapLevel > 0.0 && MipmapLevelI < this->MipmapLevel - 1) {
		SRReal FractMipmapLevel = MipmapLevel - MipmapLevelI;
		MipmapLevelI++;
		JuliaMipmap &Mipmap2 = Maps[MipmapLevelI];
		JuliaPixel Pixel2 = Mipmap2.Sample(SamplePosition.real(), SamplePosition.imag());
		Pixel = Interpolate(FractMipmapLevel, Pixel, Pixel2);
	}
	Pixel.Weight = std::clamp(MipmapLevel + 1, 0.0, 1.0);
	return Pixel;
}
HInfLAEvaluator::HInfLAEvaluator() {
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	RefReserveSize = (statex.ullTotalPhys) & ~(CommitSize - 1);
	LAReserveSize = (statex.ullTotalPhys / 2) & ~(CommitSize - 1);
	featureFinder = new FeatureFinder;
}
HInfLAEvaluator::~HInfLAEvaluator() {
	if (ref) delete ref;
	delete featureFinder;
}

template <typename real>
struct WayPoint {
	std::complex<real> Z;
	uint64_t Iteration : 63;
	uint64_t Rebase : 1;
};

template<typename real, typename orbitreal> void HInfLAEvaluator::SaveCompressedOrbit(std::ostream &stream, const LAReference<real> &Ref) {
	static_assert(std::is_same_v<real, orbitreal> || std::is_same_v<real, HRReal>);
	constexpr bool RecoverHRValue = !std::is_same_v<real, orbitreal>;

	using WayPoint = WayPoint<real>;

	auto &Orbit = reinterpret_cast<const ReservedMemoryArray<std::complex<orbitreal>>&>(Ref.Ref);
	const LAInfo<real> *RecoveryLA = nullptr;
	const LAInfoI *RecoveryLAI = nullptr;
	size_t PrevExtendedIteration = 0;
	size_t NextExtendedIteration = 0;
	size_t ExtendedLAIndex = 0;
	double Threshold = 1ULL << Global::ReferenceQuality;

	if constexpr (RecoverHRValue) {
		if (Ref.LAStages[Ref.LAStageCount - 1].UseDoublePrecision) throw ""; // FIXME
		size_t RecoveryLAStage = 0;
		while (Ref.LAStages[RecoveryLAStage].UseDoublePrecision) RecoveryLAStage++;
		size_t Index = Ref.LAStages[RecoveryLAStage].LAIndex;
		RecoveryLA = Ref.LA + Index;
		RecoveryLAI = Ref.LAI + Index;
		NextExtendedIteration = RecoveryLAI[0].StepLength;
	}

	auto GetZ = [&](size_t i) -> std::complex<real> {
		if constexpr (RecoverHRValue) {
			if (i <= PrevExtendedIteration) {
				if (i == PrevExtendedIteration) {
					return RecoveryLA[ExtendedLAIndex].Ref;
				} else {
					PrevExtendedIteration = 0;
					NextExtendedIteration = RecoveryLAI[0].StepLength;
					ExtendedLAIndex = 0;
				}
			}
			while (i >= NextExtendedIteration) {
				PrevExtendedIteration = NextExtendedIteration;
				ExtendedLAIndex++;
				NextExtendedIteration += RecoveryLAI[ExtendedLAIndex].StepLength;
				if (i == PrevExtendedIteration) {
					return RecoveryLA[ExtendedLAIndex].Ref;
				}
			}
		}
		return Orbit[i];
	};

	std::complex<real> z, c = Ref.Refc;
	ReservedMemoryArray<WayPoint> WayPoints(1024 * 1024 * 1024);
	ReservedMemoryArray<size_t> Rebases(1024 * 1024 * 1024);
	size_t i = 1;
	z = c;
	for (; i <= Ref.RefIt; i++) {
		if (ChebyshevNorm(GetZ(i)) < 0x1.0p-4) {
			z = GetZ(i);
			WayPoints.Append(WayPoint{ z, i, true });
			break;
		} else if (ChebyshevNorm(z - GetZ(i)) * Threshold > ChebyshevNorm(GetZ(i))) {
			z = GetZ(i);
			WayPoints.Append(WayPoint{ z, i, false });
		}
		z = z * z + c;
	}
	std::complex<real> dz = z;
	size_t PrevWayPointIteration = i;

	dz = (GetZ(0) * real(2.0) + dz) * dz;
	i++;
	size_t j = 1;
	for (; i <= Ref.RefIt; i++, j++) {
		z = dz + GetZ(j);
		if (j >= PrevWayPointIteration || ChebyshevNorm(z - GetZ(i)) * Threshold > ChebyshevNorm(GetZ(i))) { // TODO: Use a different threshold for uncorrectable parts of the orbit.
			PrevWayPointIteration = i;
			z = GetZ(i);
			dz = z - GetZ(j);
			if (ChebyshevNorm(z) < ChebyshevNorm(dz) || (i - j) * 4 < i) { // TODO: Improve second condition
				dz = z;
				j = 0;
				WayPoints.Append(WayPoint{ dz, i, true });
			} else {
				WayPoints.Append(WayPoint{ dz, i, false });
			}
		} else if (ChebyshevNorm(z) < ChebyshevNorm(dz) * 0x1.000001p0) {
			dz = z;
			j = 0;
			if (Rebases.size() && Rebases[Rebases.size() - 1] > WayPoints[WayPoints.size() - 1].Iteration) {
				Rebases[Rebases.size() - 1] = i;
			} else {
				Rebases.Append(i);
			}
		}
		dz = (GetZ(j) * real(2.0) + dz) * dz;
	}

	WayPoints.Serialize(stream);
	Rebases.Serialize(stream);
}
template<typename real>__declspec(noinline) void HInfLAEvaluator::LoadCompressedOrbit(std::istream &stream, LAReference<real> &Ref) {
	using WayPoint = WayPoint<real>;
	using complex = std::complex<real>;

	Ref.Ref.RecalculateCapacity();
	complex z = real(0.0), c = Ref.Refc;
	ReservedMemoryArray<WayPoint> WayPoints(1024 * 1024 * 1024);
	ReservedMemoryArray<size_t> Rebases(1024 * 1024 * 1024);
	WayPoints.Deserialize(stream);
	Rebases.Deserialize(stream);
	WayPoints.Append(WayPoint{ real(0.0), ~0ull, false }); // Dummy
	Rebases.Append(~0ull); // Dummy
	WayPoint NextWayPoint;
	size_t NextRebase;
	size_t WayPointIndex = 0;
	size_t RebaseIndex = 0;
	NextWayPoint = WayPoints[0];
	NextRebase = Rebases[0];
	auto CorrectOrbit = [&](size_t Begin, size_t End, complex diff) {
		constexpr real ScalingFactor = 0x1.0p-384;
		complex dzdc = ScalingFactor;
		for (size_t i = End; i > Begin; ) {
			i--;
			dzdc *= Ref.Ref[i] * real(2.0);
			Ref.Ref[i] += diff * (conj(dzdc) * (ScalingFactor / std::norm(dzdc)));
		}
	};
	size_t i = 0;
	size_t UncorrectedOrbitBegin = 1;
	for (; i <= Ref.RefIt; i++) {
		if (i == NextWayPoint.Iteration) {
			CorrectOrbit(UncorrectedOrbitBegin, i, NextWayPoint.Z - z);
			UncorrectedOrbitBegin = i + 1;
			z = NextWayPoint.Z;
			bool Rebase = NextWayPoint.Rebase;
			WayPointIndex++;
			NextWayPoint = WayPoints[WayPointIndex];
			if (Rebase) {
				break;
			}
		}
		Ref.Ref.Append(z);
		z = z * z + c;
	}
	size_t j = 0;
	std::complex<real> dz = z;
	for (; i <= Ref.RefIt; i++, j++) {
		z = dz + Ref.Ref[j];
		if (i == NextWayPoint.Iteration) {
			if (NextWayPoint.Rebase) {
				dz = z;
				j = 0;
			}
			CorrectOrbit(UncorrectedOrbitBegin, i, NextWayPoint.Z - dz);
			UncorrectedOrbitBegin = i + 1;
			dz = NextWayPoint.Z;
			z = dz + Ref.Ref[j];
			WayPointIndex++;
			NextWayPoint = WayPoints[WayPointIndex];
		} else if (i == NextRebase) {
			RebaseIndex++;
			NextRebase = Rebases[RebaseIndex];
			dz = z;
			j = 0;
		} else if (ChebyshevNorm(z) < ChebyshevNorm(dz)) {
			dz = z;
			j = 0;
		}
		Ref.Ref.Append(z);
		dz = (Ref.Ref[j] * real(2.0) + dz) * dz;
	}
}

struct ReferenceHeader { // TODO
	bool ExtendedRange;
};

size_t HInfLAEvaluator::SaveReference(std::ostream &stream, const Reference *reference) {
	const LAReference<> &Ref = *static_cast<const LAReference<>*>(reference);
	const LAReference<double> &DPRef = *static_cast<const LAReference<double>*>(reference);

	ReferenceHeader Header;

	Header.ExtendedRange = !Ref.DoublePrecisionPT || (Ref.LAStageCount && !Ref.LAStages[Ref.LAStageCount - 1].UseDoublePrecision);

	stream.write((const char *)&Header, sizeof(ReferenceHeader));
	stream.write((const char *)static_cast<const ReferenceTrivialContent *>(&Ref), sizeof(ReferenceTrivialContent));
	stream.write((const char *)static_cast<const LAReferenceTrivialContent *>(&Ref), sizeof(LAReferenceTrivialContent));

	if (Header.ExtendedRange) {
		if (Ref.DoublePrecisionPT) {
			SaveCompressedOrbit<HRReal, SRReal>(stream, Ref);
		} else {
			SaveCompressedOrbit(stream, Ref);
		}
	} else {
		SaveCompressedOrbit(stream, DPRef);
	}

	return 1;
}

size_t HInfLAEvaluator::LoadReference(std::istream &stream, Reference *&reference) {
	LAReference<> &Ref = *(ref ? ref : new LAReference<>(RefReserveSize, LAReserveSize));
	ref = nullptr;
	LAReference<double> &DPRef = reinterpret_cast<LAReference<double>&>(Ref);

	ReferenceHeader Header;

	stream.read((char *)&Header, sizeof(ReferenceHeader));
	stream.read((char *)static_cast<ReferenceTrivialContent *>(&Ref), sizeof(ReferenceTrivialContent));
	stream.read((char *)static_cast<LAReferenceTrivialContent *>(&Ref), sizeof(LAReferenceTrivialContent));

	Ref.Ref.clear();
	Ref.LA.clear();
	Ref.LAI.clear();

	Ref.Julia.Valid = false;

	bool Cancelled = false;

	HRReal radius = FContext.CurrentLocation.HalfH * 2.0;
	Ref.LAStageCount = 0;
	if (Header.ExtendedRange) {
		LoadCompressedOrbit(stream, Ref);
		Ref.LA.RecalculateCapacity();
		Ref.LAI.RecalculateCapacity();
		ReferenceGenerationContext Context(Ref, Cancelled);
		Context.GenerateApproximationData(radius);
	} else {
		LoadCompressedOrbit(stream, DPRef);
		DPRef.LA.RecalculateCapacity();
		DPRef.LAI.RecalculateCapacity();
		ReferenceGenerationContext Context(DPRef, Cancelled);
		Context.GenerateApproximationData(radius);
	}

	reference = &Ref;

	return 1;
}

template<typename real>
void HInfLAEvaluator::ReferenceGenerationContext<real>::GenerateReference(const EvaluationParameters &parameters) {
	real radius = real(FContext.CurrentLocation.HalfH) * 2.0;

	reference.MaxIt = parameters.MaxIterations;
	reference.Refc = SRComplex{ convert<SRReal>(parameters.CenterCoordinate.X), convert<SRReal>(parameters.CenterCoordinate.Y) };

	reference.IsPeriodic = false;
	reference.UseAT = false;

	reference.AbsolutePrecision = pow2(-double(parameters.CenterCoordinate.X.get_prec())) * 16.0_hr;

	reference.LAStageCount = 0;

	reference.Ref.clear();
	reference.LA.clear();
	reference.LAI.clear();

	reference.Ref.RecalculateCapacity();
	reference.LA.RecalculateCapacity();
	reference.LAI.RecalculateCapacity();

	if (radius > real(0x1.0p-32)) {
		reference.DirectEvaluate = true;
		reference.ShrinkToFit();
		return;
	}

	CalculateOrbit(parameters, radius);

	if (reference.RefIt < 8) {
		reference.DirectEvaluate = true;

		reference.Ref.clear();
		reference.ShrinkToFit();
		return;
	} else {
		reference.DirectEvaluate = false;
	}

	GenerateApproximationData(radius);
}

std::string_view HInfLAEvaluator::HInfLAReferenceTask::GetDescription() const {
	using namespace std;
	return "Computing reference"sv;
}

void HInfLAEvaluator::HInfLAReferenceTask::Execute() {
	try {
		if (FContext.CurrentLocation.HalfH > 0x1.0p-896_hr) {
			ReferenceGenerationContext<SRReal> Context(static_cast<LAReference<SRReal>&>(*reference), Cancelled);
			UseHRReal = false;
			ComputationContext = reinterpret_cast<ReferenceGenerationContext<HRReal> *> (&Context);
			Context.GenerateReference(parameters);
		} else {
			ReferenceGenerationContext<HRReal> Context(static_cast<LAReference<HRReal>&>(*reference), Cancelled);
			UseHRReal = true;
			ComputationContext = &Context;
			Context.GenerateReference(parameters);
		}
	} catch (Cancellation) {
		evaluator->ref = static_cast<LAReference<> *>(reference);
		reference = nullptr;
	}
}

void HInfLAEvaluator::HInfLAReferenceTask::Cancel() {
	Cancelled = true;
}

bool HInfLAEvaluator::HInfLAReferenceTask::GetProgress(SRReal &Numerator, SRReal &Denoninator) const {
	if (!ComputationContext) {
		return false;
	}
	if (ComputationContext->MaxIt == 0) {
		return false;
	}
	Numerator = ComputationContext->Iteration;
	Denoninator = 0.0;
	return true;
}

Evaluator::ReferenceTask *HInfLAEvaluator::CreateReferenceTask(const EvaluationParameters &parameters) {
	LAReference<> *reference = ref ? ref : new LAReference<>(RefReserveSize, LAReserveSize);
	ref = nullptr;
	return new HInfLAReferenceTask(this, parameters, reference);
}

template <typename T> struct IsFloatExp { constexpr operator bool() { return false; }; };
template <> struct IsFloatExp<FExpDouble> { constexpr operator bool() { return true; }; };

template<typename real>
void HInfLAEvaluator::ReferenceGenerationContext<real>::CalculateOrbit(const EvaluationParameters &parameters, HRReal radius) {
	using complex = std::complex<real>;

	reference.DoublePrecisionPT = !IsFloatExp<real>();

	MaxIt = reference.MaxIt = parameters.MaxIterations;
	complex dzdc = real(1.0);

	using HPReal2 = HPReal;
	OMITemp<HPReal2> T;

	HPReal2 CR = parameters.CenterCoordinate.X;
	HPReal2 CI = parameters.CenterCoordinate.Y;
	HPReal2 ZR = CR;
	HPReal2 ZI = CI;
	complex z = { convert<real>(ZR), convert<real>(ZI) };
	reference.Refc = SRComplex{ convert<SRReal>(CR), convert<SRReal>(CI) };
	reference.Ref.Append(complex(0.0));

	size_t i = 1;
	for (; i <= MaxIt; i++) {
		Iteration = i;
		if (Cancelled) throw Cancellation();

		reference.Ref.Append(z);

		if (radius * ChebyshevNorm(dzdc) * 2 > ChebyshevNorm(z)) {
			complex PredCenter = z / dzdc;
			FContext.PredictCenter(RelLocation{ -PredCenter.real(), -PredCenter.imag(), 0 });
			i++;
			reference.IsPeriodic = true;
			break;
		}
		OMI(ZR, ZI, CR, CI, T);

		dzdc = real(2.0) * z * dzdc + real(1.0);
		z = complex{ convert<real>(ZR), convert<real>(ZI) };

		real Magnitude = std::norm(z);

		if (Magnitude > 5.0) {
			i++;
			break;
		}
	}
	reference.Ref.ShrinkToFit();
	reference.RefIt = i - 1;
}

template<typename real>
bool HInfLAEvaluator::ReferenceGenerationContext<real>::CreateLAFromOrbit(HRReal radius) {
	reference.LAStages[0].UseDoublePrecision = !IsFloatExp<real>();

	reference.LAStages[0].LAIndex = 0;

	size_t Period = 0;

	LAInfo<real> LA = LAInfo<real>(real(0.0));
	LA = LA.Step(reference.Ref[1]);
	LAInfoI LAI;
	LAI.NextStageLAIndex = 0;

	size_t i;
	for (i = 2; i < reference.RefIt; i++) {
		if (Cancelled) throw Cancellation();
		LAInfo<real> NewLA;
		bool PeriodDetected = LA.Step(NewLA, reference.Ref[i]);

		if (PeriodDetected) {
			Period = i;
			LAI.StepLength = Period;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			LAI.NextStageLAIndex = i;

			if (i + 1 < reference.RefIt) {
				LA = LAInfo<real>(reference.Ref[i]).Step(reference.Ref[i + 1]);
				i += 2;
			} else {
				LA = LAInfo<real>(reference.Ref[i]);
				i += 1;
			}
			break;
		}
		LA = NewLA;
	}

	reference.LAStageCount = 1;

	size_t PeriodBegin = Period;
	size_t PeriodEnd = PeriodBegin + Period;

	if (!Period) {
		if (reference.RefIt > 64) {
			LA = LAInfo<real>(reference.Ref[0]).Step(reference.Ref[1]);
			LAI.NextStageLAIndex = 0;
			i = 2;

			double NthRoot = round(log2(double(reference.RefIt)) / 4); // log16
			Period = round(pow(double(reference.RefIt), 1.0 / NthRoot));

			PeriodBegin = 0;
			PeriodEnd = Period;
		} else {
			LAI.StepLength = reference.RefIt;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));

			reference.LAStages[0].MacroItCount = 1;
			reference.LA.ShrinkToFit();
			reference.LAI.ShrinkToFit();

			return false;
		}
	} else if (Period > 64) {
		reference.LA.Pop();
		reference.LAI.Pop();

		LA = LAInfo<real>(reference.Ref[0]).Step(reference.Ref[1]);
		LAI.NextStageLAIndex = 0;
		i = 2;

		double NthRoot = round(log2(double(Period)) / 4); // log16
		Period = round(pow(double(Period), 1.0 / NthRoot));

		PeriodBegin = 0;
		PeriodEnd = Period;
	}

	for (; i < reference.RefIt; i++) {
		LAInfo<real> NewLA;
		bool PeriodDetected = LA.Step(NewLA, reference.Ref[i]);

		if (PeriodDetected || i >= PeriodEnd) {
			LAI.StepLength = i - PeriodBegin;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			LAI.NextStageLAIndex = i;
			PeriodBegin = i;
			PeriodEnd = PeriodBegin + Period;
			if (NewLA.DetectPeriod(reference.Ref[i + 1]) || i + 1 >= reference.RefIt) {
				LA = LAInfo<real>(reference.Ref[i]);
			} else {
				LA = LAInfo<real>(reference.Ref[i]).Step(reference.Ref[i + 1]);
				i++;
			}
		} else {
			LA = NewLA;
		}
	}

	LAI.StepLength = i - PeriodBegin;

	reference.LA.Append(LA);
	reference.LAI.Append(LAI);

	reference.LAStages[0].MacroItCount = reference.LA.size();

	reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));
	reference.LAI.Append(LAInfoI());

	return true;
}

template<typename real>
bool HInfLAEvaluator::ReferenceGenerationContext<real>::CreateLAFromOrbit_MT(HRReal radius) {
	static constexpr size_t ThreadCount = 8;
	reference.LAStages[0].UseDoublePrecision = !IsFloatExp<real>();

	reference.LAStages[0].LAIndex = 0;

	size_t Period = 0;

	LAInfo<real> LA = LAInfo<real>(real(0.0));
	LA = LA.Step(reference.Ref[1]);
	LAInfoI LAI;
	LAI.NextStageLAIndex = 0;

	size_t i;
	for (i = 2; i < reference.RefIt; i++) {
		LAInfo<real> NewLA;
		bool PeriodDetected = LA.Step(NewLA, reference.Ref[i]);

		if (PeriodDetected) {
			Period = i;
			LAI.StepLength = Period;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			LAI.NextStageLAIndex = i;

			if (i + 1 < reference.RefIt) {
				LA = LAInfo<real>(reference.Ref[i]).Step(reference.Ref[i + 1]);
				i += 2;
			} else {
				LA = LAInfo<real>(reference.Ref[i]);
				i += 1;
			}
			break;
		}
		LA = NewLA;
	}

	reference.LAStageCount = 1;

	if (!Period) {
		LAI.StepLength = reference.RefIt;

		reference.LA.Append(LA);
		reference.LAI.Append(LAI);

		reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));

		reference.LAStages[0].MacroItCount = 1;
		reference.LA.ShrinkToFit();
		reference.LAI.ShrinkToFit();

		return false;
	}

	size_t PeriodBegin = Period;
	size_t PeriodEnd = PeriodBegin + Period;

	if (Period > 64) {
		reference.LA.Pop();
		reference.LAI.Pop();

		LA = LAInfo<real>(reference.Ref[0]).Step(reference.Ref[1]);
		LAI.NextStageLAIndex = 0;
		i = 2;

		Period = round(sqrt(double(Period)));

		PeriodBegin = 0;
		PeriodEnd = Period;
	}

	size_t Start[ThreadCount] = { 0 };

	ReservedMemoryArray<LAInfo<real>> LAs[ThreadCount];
	for (size_t i = 1; i < ThreadCount; i++) {
		LAs[i] = ReservedMemoryArray<LAInfo<real>>(reference.LA.GetReservedSize());
	}
	ReservedMemoryArray<LAInfoI> LAIs[ThreadCount];
	for (size_t i = 1; i < ThreadCount; i++) {
		LAIs[i] = ReservedMemoryArray<LAInfoI>(reference.LAI.GetReservedSize());
	}

	std::thread Workers[ThreadCount];

	Workers[0] = std::thread([&]() {
		for (; i < reference.RefIt; i++) {
			LAInfo<real> NewLA;
			bool PeriodDetected = LA.Step(NewLA, reference.Ref[i]);

			if (PeriodDetected || i >= PeriodEnd) {
				LAI.StepLength = i - PeriodBegin;

				reference.LA.Append(LA);
				reference.LAI.Append(LAI);

				LAI.NextStageLAIndex = i;
				PeriodBegin = i;
				PeriodEnd = PeriodBegin + Period;
				if (NewLA.DetectPeriod(reference.Ref[i + 1]) || i + 1 >= reference.RefIt) {
					LA = LAInfo<real>(reference.Ref[i]);
				} else {
					LA = LAInfo<real>(reference.Ref[i]).Step(reference.Ref[i + 1]);
					i++;
				}
				if (i > reference.RefIt / ThreadCount) {
					while (!Start[1]);
					if (i == Start[1] - 1) {
						i++;
						break;
					} else if (i >= Start[1]) {
						throw "";
					}
				}
			} else {
				LA = NewLA;
			}
		}
		});
	auto Worker = [Period, &Start, &LAs, &LAIs, this](size_t ThreadID) {
		size_t j = reference.RefIt * ThreadID / ThreadCount;
		size_t End = reference.RefIt * (ThreadID + 1) / ThreadCount;
		LAInfo<real> LA_ = LAInfo<real>(reference.Ref[j]);
		LA_ = LA_.Step(reference.Ref[j + 1]);
		LAInfoI LAI_;
		LAI_.NextStageLAIndex = j;
		j += 2;

		size_t PeriodBegin;
		size_t PeriodEnd;

		for (; j < reference.RefIt; j++) {
			LAInfo<real> NewLA;
			bool PeriodDetected = LA_.Step(NewLA, reference.Ref[j]);

			if (PeriodDetected) {
				LAI_.NextStageLAIndex = j;
				PeriodBegin = j;
				PeriodEnd = PeriodBegin + Period;

				if (j + 1 < reference.RefIt) {
					LA_ = LAInfo<real>(reference.Ref[j]).Step(reference.Ref[j + 1]);
					j += 2;
				} else {
					LA_ = LAInfo<real>(reference.Ref[j]);
					j += 1;
				}
				break;
			}
			LA_ = NewLA;
		}
		Start[ThreadID] = j;
		for (; j < reference.RefIt; j++) {
			LAInfo<real> NewLA;
			bool PeriodDetected = LA_.Step(NewLA, reference.Ref[j]);

			if (PeriodDetected || j >= PeriodEnd) {
				LAI_.StepLength = j - PeriodBegin;

				LAs[ThreadID].Append(LA_);
				LAIs[ThreadID].Append(LAI_);

				LAI_.NextStageLAIndex = j;
				PeriodBegin = j;
				PeriodEnd = PeriodBegin + Period;
				if (NewLA.DetectPeriod(reference.Ref[j + 1]) || j + 1 >= reference.RefIt) {
					LA_ = LAInfo<real>(reference.Ref[j]);
				} else {
					LA_ = LAInfo<real>(reference.Ref[j]).Step(reference.Ref[j + 1]);
					j++;
				}
				if (j > End) {
					if (ThreadID == ThreadCount - 1) throw "";
					while (!Start[ThreadID + 1]);
					if (j == Start[ThreadID + 1] - 1) {
						j++;
						break;
					} else if (j >= Start[ThreadID + 1]) {
						throw "";
					}
				}
			} else {
				LA_ = NewLA;
			}
		}
		if (ThreadID == ThreadCount - 1) {
			LAI_.StepLength = j - PeriodBegin;

			LAs[ThreadID].Append(LA_);
			LAIs[ThreadID].Append(LAI_);
		}
	};

	for (size_t i = 1; i < ThreadCount; i++) {
		Workers[i] = std::thread(Worker, i);
	}

	for (size_t i = 0; i < ThreadCount; i++) {
		Workers[i].join();
	}

	for (size_t i = 1; i < ThreadCount; i++) {
		reference.LA.Concatenate(std::move(LAs[i]));
		reference.LAI.Concatenate(std::move(LAIs[i]));
	}

	reference.LAStages[0].MacroItCount = reference.LA.size();

	reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));
	reference.LAI.Append(LAInfoI());

	return true;
}

template<typename real>
bool HInfLAEvaluator::ReferenceGenerationContext<real>::CreateNewLAStage() {
	LAInfo<real> LA;
	LAInfoI LAI;
	size_t i;
	size_t PeriodBegin;
	size_t PeriodEnd;

	size_t PrevStage = reference.LAStageCount - 1;
	size_t CurrentStage = reference.LAStageCount;
	size_t PrevStageLAIndex = reference.LAStages[PrevStage].LAIndex;
	size_t PrevStageMacroItCount = reference.LAStages[PrevStage].MacroItCount;
	LAInfo<real> *PrevStageLA = reference.LA + PrevStageLAIndex;
	LAInfoI *PrevStageLAI = reference.LAI + PrevStageLAIndex;

	size_t Period = 0;

	reference.LAStages[CurrentStage].UseDoublePrecision = !IsFloatExp<real>();
	reference.LAStages[CurrentStage].LAIndex = reference.LA.size();

	LA = PrevStageLA[0].Composite(PrevStageLA[1]);
	LAI.NextStageLAIndex = 0;
	i = PrevStageLAI[0].StepLength + PrevStageLAI[1].StepLength;
	size_t j;

	for (j = 2; j < PrevStageMacroItCount; j++) {
		LAInfo<real> NewLA;
		bool PeriodDetected = LA.Composite(NewLA, PrevStageLA[j]);

		if (PeriodDetected) {
			if (PrevStageLA[j].LAThreshold == 0.0) break;
			Period = i;

			LAI.StepLength = Period;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			LAI.NextStageLAIndex = j;
			if (NewLA.DetectPeriod(PrevStageLA[j + 1].Ref) || j + 1 >= PrevStageMacroItCount) {
				LA = PrevStageLA[j];
				i += PrevStageLAI[j].StepLength;
				j++;
			} else {
				LA = PrevStageLA[j].Composite(PrevStageLA[j + 1]);
				i += PrevStageLAI[j].StepLength + PrevStageLAI[j + 1].StepLength;
				j += 2;
			}
			break;
		}
		LA = NewLA;
		i += PrevStageLAI[j].StepLength;
	}
	reference.LAStageCount++; if (reference.LAStageCount > MaxLAStages) throw "Too many stages"; // FIXME

	PeriodBegin = Period;
	PeriodEnd = PeriodBegin + Period;

	if (!Period) {
		if (reference.RefIt / PrevStageLAI[0].StepLength > 64) {
			LA = PrevStageLA[0].Composite(PrevStageLA[1]);
			i = PrevStageLAI[0].StepLength + PrevStageLAI[1].StepLength;
			LAI.NextStageLAIndex = 0;

			j = 2;

			double Ratio = double(reference.RefIt) / PrevStageLAI[0].StepLength;
			double NthRoot = round(log2(double(Ratio)) / 4); // log16
			Period = PrevStageLAI[0].StepLength * (size_t)round(pow(double(Ratio), 1.0 / NthRoot));

			PeriodBegin = 0;
			PeriodEnd = Period;
		} else {
			LAI.StepLength = reference.RefIt;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));
			reference.LAStages[CurrentStage].MacroItCount = 1;

			return false;
		}
	} else if (Period > PrevStageLAI[0].StepLength * 64) {
		reference.LA.Pop();
		reference.LAI.Pop();

		LA = PrevStageLA[0].Composite(PrevStageLA[1]);
		i = PrevStageLAI[0].StepLength + PrevStageLAI[1].StepLength;
		LAI.NextStageLAIndex = 0;

		j = 2;

		double Ratio = double(Period) / PrevStageLAI[0].StepLength;
		double NthRoot = round(log2(double(Ratio)) / 4); // log16
		Period = PrevStageLAI[0].StepLength * (size_t)round(pow(double(Ratio), 1.0 / NthRoot));

		PeriodBegin = 0;
		PeriodEnd = Period;
	}

	for (; j < PrevStageMacroItCount; j++) {
		LAInfo<real> NewLA;
		bool PeriodDetected = LA.Composite(NewLA, PrevStageLA[j]);

		if (PeriodDetected || i >= PeriodEnd) {
			LAI.StepLength = i - PeriodBegin;

			reference.LA.Append(LA);
			reference.LAI.Append(LAI);

			LAI.NextStageLAIndex = j;
			PeriodBegin = i;
			PeriodEnd = PeriodBegin + Period;
			if (NewLA.DetectPeriod(PrevStageLA[j + 1].Ref) || j + 1 >= PrevStageMacroItCount) {
				LA = PrevStageLA[j];
			} else {
				LA = PrevStageLA[j].Composite(PrevStageLA[j + 1]);
				i += PrevStageLAI[j].StepLength;
				j++;
			}
		} else {
			LA = NewLA;
		}
		i += PrevStageLAI[j].StepLength;
	}

	LAI.StepLength = i - PeriodBegin;

	reference.LA.Append(LA);
	reference.LAI.Append(LAI);

	reference.LAStages[CurrentStage].MacroItCount = reference.LA.size() - reference.LAStages[CurrentStage].LAIndex;

	reference.LA.Append(LAInfo<real>(reference.Ref[reference.RefIt]));
	reference.LAI.Append(LAInfoI());
	return true;
}

template<typename real>
void HInfLAEvaluator::ReferenceGenerationContext<real>::CreateATFromLA(HRReal SqrRadius) {
	constexpr bool UseFloatExp = IsFloatExp<real>();

	for (size_t Stage = reference.LAStageCount; Stage > 0; ) {
		Stage--;
		size_t LAIndex = reference.LAStages[Stage].LAIndex;
		if (UseFloatExp && reference.LAStages[Stage].UseDoublePrecision) {
			LAReference<double> &DPRef = reinterpret_cast<LAReference<double>&>(reference);

			reference.AT = DPRef.LA[LAIndex].CreateAT();
		} else {
			reference.AT = reference.LA[LAIndex].CreateAT();
		}
		reference.AT.StepLength = reference.LAI[LAIndex].StepLength;
		if (reference.AT.Usable(SqrRadius)) {
			reference.UseAT = true;
			return;
		}
	}
	reference.UseAT = false;
}

template<typename real>
void HInfLAEvaluator::ReferenceGenerationContext<real>::ConvertStageToDouble(size_t Stage) {
	size_t LAIndex = reference.LAStages[Stage].LAIndex;
	size_t MacroItCount = reference.LAStages[Stage].MacroItCount;
	LAInfo<real> *LA = reference.LA + LAIndex;

	reference.LAStages[Stage].UseDoublePrecision = true;
	LAReference<double> &DPRef = reinterpret_cast<LAReference<double>&>(reference);
	LAInfo<double> *DPLA = DPRef.LA + LAIndex;
	for (size_t i = 0; i <= MacroItCount; i++) {
		DPLA[i] = LA[i];
	}
}

template<typename real>
void HInfLAEvaluator::ReferenceGenerationContext<real>::GenerateApproximationData(HRReal radius) {
	constexpr bool UseFloatExp = IsFloatExp<real>();
	HRReal SqrRadius = radius * radius;

	bool PeriodDetected = CreateLAFromOrbit(radius);

	auto DoubleUsableAtPrevStage = [&](size_t Stage) {
		size_t LAIndex = reference.LAStages[Stage].LAIndex;
		return radius > real(0x1.0p-896) || (reference.LA[LAIndex].LAThreshold > real(0x1.0p-768) && reference.LA[LAIndex].LAThresholdC > real(0x1.0p-768));
	};

	if constexpr (UseFloatExp) {
		if (DoubleUsableAtPrevStage(0)) {
			reference.DoublePrecisionPT = true;
			LAReference<double> &DPRef = reinterpret_cast<LAReference<double>&>(reference);
			for (size_t i = 0; i <= reference.RefIt; i++) {
				DPRef.Ref[i] = convert<SRComplex>(reference.Ref[i]);
			}
		}
	}

	if (!PeriodDetected) return;

	while (true) {
		if (Cancelled) throw Cancellation();
		size_t PrevStage = reference.LAStageCount - 1;
		size_t CurrentStage = reference.LAStageCount;

		bool PeriodDetected = CreateNewLAStage();

		if constexpr (UseFloatExp) {
			if (DoubleUsableAtPrevStage(CurrentStage)) {
				ConvertStageToDouble(PrevStage);
			}
		}
		if (!PeriodDetected) break;
	}

	CreateATFromLA(SqrRadius);

	if constexpr (UseFloatExp) {
		if (radius > real(0x1.0p-896)) {
			ConvertStageToDouble(reference.LAStageCount - 1);
		}
	}
	reference.LA.ShrinkToFit();
	reference.LAI.ShrinkToFit();
}

void HInfLAEvaluator::FreeReference(Reference *reference) {
	if (ref) {
		delete reference;
	} else {
		ref = static_cast<LAReference<> *>(reference);
		ref->Julia.Valid = false;
	}
}

#define VecGatherFirst() {\
	j = Select(ValidMask, RefIteration, ivec(0x7FFFFFFFFFFFFFFF)).HorizontalMin();\
	ActiveMask = ValidMask && (RefIteration == ivec(j));\
	GatherMask = AndNot(ActiveMask, ValidMask);\
	NextGatherIteration = Select(GatherMask, RefIteration, ivec(0x7FFFFFFFFFFFFFFF)).HorizontalMin();\
}

#define VecGatherAndSync() [&]() {\
	if (!ActiveMask) {\
		if (GatherMask) {\
			j = NextGatherIteration;\
		} else if (PauseMask) {\
			z = InactiveZ;\
			dz = InactiveDz;\
			if constexpr (EvaluateDzdz) derivatives::dzdz = InactiveDzdz;\
			if constexpr (EvaluateDzdc) derivatives::dzdc = InactiveDzdc;\
			ActiveMask = PauseMask;\
			PauseMask = false;\
			j = 0;\
		} else return false;\
	}\
	if (j == NextGatherIteration) {\
		mask ActivatingMask = ValidMask && (RefIteration == ivec(j)) && GatherMask;\
		ActiveMask = ActiveMask || ActivatingMask;\
		GatherMask = AndNot(ActivatingMask, GatherMask);\
		z = Select(ActivatingMask, InactiveZ, z);\
		dz = Select(ActivatingMask, InactiveDz, dz);\
		if constexpr (EvaluateDzdz) derivatives::dzdz = Select(ActivatingMask, InactiveDzdz, derivatives::dzdz);\
		if constexpr (EvaluateDzdc) derivatives::dzdc = Select(ActivatingMask, InactiveDzdc, derivatives::dzdc);\
		NextGatherIteration = Select(GatherMask, RefIteration, ivec(0x7FFFFFFFFFFFFFFF)).HorizontalMin();\
	}\
	return true;\
}()

#define VecBreakIf(x) {																\
	ActiveMask = AndNot((x), ActiveMask);									\
	if (!ActiveMask) continue;														\
}

#define VecResetIf(x) {																\
	mask __x = (x) && ActiveMask;													\
	if (__x) {																		\
		InactiveZ = Select((__x), z, InactiveZ);										\
		InactiveDz = Select((__x), z, InactiveDz);									\
		if constexpr (EvaluateDzdz) InactiveDzdz = Select((__x), derivatives::dzdz, InactiveDzdz);			\
		if constexpr (EvaluateDzdc) InactiveDzdc = Select((__x), derivatives::dzdc, InactiveDzdc);			\
		PauseMask = (__x) || PauseMask;									\
		ActiveMask = AndNot((__x), ActiveMask);							\
		if (!ActiveMask) continue;													\
	}																				\
}

#define VecResetAndSync() {															\
	z = Select(PauseMask, InactiveZ, z);											\
	dz = Select(PauseMask, InactiveDz, z);											\
	if constexpr (EvaluateDzdz) derivatives::dzdz = Select(PauseMask, InactiveDzdz, derivatives::dzdz);			\
	if constexpr (EvaluateDzdc) derivatives::dzdc = Select(PauseMask, InactiveDzdc, derivatives::dzdc);			\
	ActiveMask = PauseMask || ActiveMask;								\
	PauseMask = false;												\
	j = 0;																			\
}

template<typename real, bool EvaluateDzdc, bool EvaluateDzdz>
void HInfLAEvaluator::EvaluationContext<real, EvaluateDzdc, EvaluateDzdz>::TrySampleJulia(mask &ActiveMask) {
	if (!Global::HighQuality || !Ref.Julia.Maps) return;
	vreal Magnitude = std::norm(z);
	mask cmp = ActiveMask
		&& SamplingEnabledMask
		&& Magnitude * SqrJuliaHalfPixelSize < std::norm(derivatives::dzdc) * SqrPixScale
		&& Magnitude > vreal(SqrJuliaMinMagZ)
		&& Magnitude < vreal(SqrJuliaMaxMagZ);
	mask SamplingShouldDisable = Magnitude * SqrJuliaPixelSize * vreal(8.0) < std::norm(derivatives::dzdc) * SqrPixScale && Magnitude < vreal(SqrJuliaMinMagZ);
	SamplingDisablingMask = SamplingShouldDisable || SamplingDisablingMask;
	if (bool(cmp)) {
		SamplingMask = SamplingMask || cmp;

		mask DirectSamplingMask = cmp && Magnitude * SqrJuliaPixelSize < std::norm(derivatives::dzdc) * SqrPixScale;
		mask DirectSamplingMask2 = cmp && Magnitude * SqrJuliaPixelSize * vreal(16.0) < std::norm(derivatives::dzdc) *SqrPixScale;
		MixedSamplingMask = AndNot(DirectSamplingMask, MixedSamplingMask || cmp);

		SamplingIteration = Select(cmp, Iteration, SamplingIteration);
		SamplingZ = Select(cmp, z, SamplingZ);
		SamplingDzdc = Select(cmp, derivatives::dzdc, SamplingDzdc);
		SamplingUpdatedMask = SamplingUpdatedMask || cmp;
		ValidMask = AndNot(DirectSamplingMask2, ValidMask);
		ActiveMask = AndNot(DirectSamplingMask2, ActiveMask);
	}
}

template<typename real, bool EvaluateDzdc, bool EvaluateDzdz>
void HInfLAEvaluator::EvaluationContext<real, EvaluateDzdc, EvaluateDzdz>::EvaluateLA() {
	size_t LAIndex = Ref.LAStages[CurrentLAStage].LAIndex;
	LAInfo<real> *LA = Ref.LA + LAIndex;
	LAInfoI *LAI = Ref.LAI + LAIndex;
	size_t MacroItCount = Ref.LAStages[CurrentLAStage].MacroItCount;
	mask ActiveMask = ValidMask;
	mask PauseMask = false;
	mask GatherMask = false;
	size_t NextGatherIteration = 0x7FFFFFFFFFFFFFFF;

	vcomplex InactiveZ = z;
	vcomplex InactiveDz = dz;
	vcomplex InactiveDzdz;
	vcomplex InactiveDzdc;
	if constexpr (EvaluateDzdz) InactiveDzdz = derivatives::dzdz;
	if constexpr (EvaluateDzdc) InactiveDzdc = derivatives::dzdc;

	size_t j = 0;
	if ((ChebyshevNorm(dc) >= vreal(LA[j].LAThresholdC)) && ValidMask) {
		j = LAI[j].NextStageLAIndex;
		return;
	}

	VecGatherFirst();

	while (VecGatherAndSync()) {
		if (Cancelled) throw Cancellation();
		auto temp = LA[j].Prepare(dz, dc);

		temp.unusable = temp.unusable && ActiveMask;
		if (temp.unusable) {
			InactiveZ = Select(temp.unusable, z, InactiveZ);
			InactiveDz = Select(temp.unusable, dz, InactiveDz);
			if constexpr (EvaluateDzdz) InactiveDzdz = Select(temp.unusable, derivatives::dzdz, InactiveDzdz);
			if constexpr (EvaluateDzdc) InactiveDzdc = Select(temp.unusable, derivatives::dzdc, InactiveDzdc);
			RefIteration = Select(temp.unusable, ivec(LAI[j].NextStageLAIndex), RefIteration);
			ActiveMask = AndNot(temp.unusable, ActiveMask);
			if (!ActiveMask) {
				continue;
			}
		}

		if constexpr (EvaluateDzdc) {
			if constexpr (EvaluateDzdz) {
				LA[j].EvaluateDerivatives(temp, dz, derivatives::dzdz, derivatives::dzdc, dc, derivatives::ScalingFactor);
			} else {
				LA[j].EvaluateDzdc(temp, dz, derivatives::dzdc, dc, derivatives::ScalingFactor);
			}
		}
		LA[j].Evaluate(temp, dz, dc);
		z = dz + vcomplex(LA[j + 1].Ref);

		Iteration += ActiveMask & ivec(LAI[j].StepLength);
		mask CompareResult = (Iteration < ivec(MaxIt));
		ActiveMask = ActiveMask && CompareResult;
		TrySampleJulia(ActiveMask);
		if (!ActiveMask) continue;

		j++;
		if (j == MacroItCount) {
			VecResetAndSync();
		} else {
			VecResetIf(ChebyshevNorm(dz) > ChebyshevNorm(z));
		}
	}
	z = InactiveZ;
	dz = InactiveDz;
	if constexpr (EvaluateDzdz) derivatives::dzdz = InactiveDzdz;
	if constexpr (EvaluateDzdc) derivatives::dzdc = InactiveDzdc;
}

template<typename real, bool EvaluateDzdc, bool EvaluateDzdz>
void HInfLAEvaluator::EvaluationContext<real, EvaluateDzdc, EvaluateDzdz>::EvaluatePT() {
	size_t j = 0;
	mask ActiveMask = ValidMask;
	mask PauseMask = false;
	mask GatherMask = false;
	size_t NextGatherIteration = 0x7FFFFFFFFFFFFFFF;

	vcomplex InactiveZ = z;
	vcomplex InactiveDz = dz;
	vcomplex InactiveDzdz;
	vcomplex InactiveDzdc;
	if constexpr (EvaluateDzdz) InactiveDzdz = derivatives::dzdz;
	if constexpr (EvaluateDzdc) InactiveDzdc = derivatives::dzdc;

	VecGatherFirst();

	while (VecGatherAndSync()) {
		if (Cancelled) throw Cancellation();
		if constexpr (EvaluateDzdz) derivatives::dzdz = derivatives::dzdz * vreal(2.0) * (dz + vcomplex(Ref.Ref[j]));
		if constexpr (EvaluateDzdc) derivatives::dzdc = derivatives::dzdc * vreal(2.0) * (dz + vcomplex(Ref.Ref[j])) + derivatives::ScalingFactor;

		dz = dz * vreal(2.0) * vcomplex(Ref.Ref[j]) + dz * dz + dc;

		Iteration -= ivec(ActiveMask);

		j++;
		z = dz + vcomplex(Ref.Ref[j]);
		vreal Magnitude = std::norm(z);
		TrySampleJulia(ActiveMask);
		mask CompareResult = (Magnitude > vreal(4096.0) || (Iteration >= ivec(MaxIt))) && ActiveMask;
		if (CompareResult) {
			InactiveZ = Select(CompareResult, z, InactiveZ);
			if constexpr (EvaluateDzdz) InactiveDzdz = Select(CompareResult, derivatives::dzdz, InactiveDzdz);
			if constexpr (EvaluateDzdc) InactiveDzdc = Select(CompareResult, derivatives::dzdc, InactiveDzdc);

			VecBreakIf(CompareResult);
		}
		if (j == Ref.RefIt) {
			VecResetAndSync();
		} else {
			VecResetIf(std::norm(dz) * real(64) > Magnitude);
		}
	}
	z = InactiveZ;
	dz = InactiveDz;
	if constexpr (EvaluateDzdz) derivatives::dzdz = InactiveDzdz;
	if constexpr (EvaluateDzdc) derivatives::dzdc = InactiveDzdc;
}

template<typename real, bool EvaluateDzdc, bool EvaluateDzdz>
void HInfLAEvaluator::EvaluationContext<real, EvaluateDzdc, EvaluateDzdz>::EvaluateDirectly() {
	mask ActiveMask = ValidMask;

	vcomplex InactiveZ;
	vcomplex InactiveDzdz;
	vcomplex InactiveDzdc;

	while (ActiveMask) {
		if (Cancelled) throw Cancellation();
		if constexpr (EvaluateDzdz) derivatives::dzdz = derivatives::dzdz * vreal(2.0) * z;
		if constexpr (EvaluateDzdc) derivatives::dzdc = derivatives::dzdc * vreal(2.0) * z + derivatives::ScalingFactor;

		z = c + z * z;

		Iteration -= ivec(ActiveMask.ymm);

		vreal Magnitude = std::norm(z);

		mask CompareResult = (Magnitude > vreal(4096.0) || (Iteration >= ivec(MaxIt))) && ActiveMask;
		if (CompareResult) {
			InactiveZ = Select(CompareResult, z, InactiveZ);
			if constexpr (EvaluateDzdz) InactiveDzdz = Select(CompareResult, derivatives::dzdz, InactiveDzdz);
			if constexpr (EvaluateDzdc) InactiveDzdc = Select(CompareResult, derivatives::dzdc, InactiveDzdc);

			ActiveMask = AndNot(CompareResult, ActiveMask);
		}
	}
	z = InactiveZ;
	if constexpr (EvaluateDzdz) derivatives::dzdz = InactiveDzdz;
	if constexpr (EvaluateDzdc) derivatives::dzdc = InactiveDzdc;
}

template<typename real, bool EvaluateDzdc, bool EvaluateDzdz>
void HInfLAEvaluator::EvaluationContext<real, EvaluateDzdc, EvaluateDzdz>::EvaluateAT(const real &__restrict SqrEscapeRadius) {
	mask ActiveMask = ValidMask;

	vcomplex InactiveZ;
	vcomplex InactiveDzdz;
	vcomplex InactiveDzdc;

	while (true) {
		if (Cancelled) throw Cancellation();
		if constexpr (EvaluateDzdz) derivatives::dzdz = derivatives::dzdz * vreal(2.0) * z;
		if constexpr (EvaluateDzdc) derivatives::dzdc = derivatives::dzdc * vreal(2.0) * z + derivatives::ScalingFactor;

		z = z * z + c;

		Iteration -= ivec(ActiveMask);
		mask CompareResult = Iteration < ivec(MaxIt);
		ActiveMask = CompareResult && ActiveMask;
		if (!ActiveMask) break;

		vreal Magnitude = std::norm(z);

		CompareResult = (Magnitude > vreal(SqrEscapeRadius)) && ActiveMask;
		if (CompareResult) {
			InactiveZ = Select(CompareResult, z, InactiveZ);
			if constexpr (EvaluateDzdz) InactiveDzdz = Select(CompareResult, derivatives::dzdz, InactiveDzdz);
			if constexpr (EvaluateDzdc) InactiveDzdc = Select(CompareResult, derivatives::dzdc, InactiveDzdc);
			ActiveMask = AndNot(CompareResult, ActiveMask);
			if (!ActiveMask) break;
		}
	}
	z = InactiveZ;
	if constexpr (EvaluateDzdz) derivatives::dzdz = InactiveDzdz;
	if constexpr (EvaluateDzdc) derivatives::dzdc = InactiveDzdc;
}

Task *HInfLAEvaluator::CreateEvaluationTask(const EvaluationParameters &parameters, PixelManager *rasterizer, Reference *reference) {
	return new HInfLAEvaluationTask(parameters, rasterizer, reference);
}

void HInfLAEvaluator::HInfLAEvaluationTask::Execute(size_t ThreadID) {
	GroupedRasterizingInterface &RI = rasterizer->GetGroupedInterface(VSize);
	LAReference<> &Ref = *static_cast<LAReference<>*>(reference);
	if (!Ref.Julia.Valid) {
		JuliaContext.RenderJulia(parameters, static_cast<LAReference<>&>(*reference), ThreadID == 0);
	}
	try {
		if (FContext.UsingDE) {
			Evaluate<true>(RI, reference);
		} else {
			//throw "";
			Evaluate<false>(RI, reference);
		}
	} catch (Cancellation) {

	}
	rasterizer->FreeInterface(RI);
}

void HInfLAEvaluator::HInfLAEvaluationTask::Cancel() {
	Cancelled = true;
}

bool HInfLAEvaluator::HInfLAEvaluationTask::GetProgress(SRReal &Numerator, SRReal &Denoninator) const {
	const ProgressTrackable *progressTrackable = dynamic_cast<const ProgressTrackable *>(rasterizer);
	if (progressTrackable) {
		return progressTrackable->GetProgress(Numerator, Denoninator);
	}
	return false;
}

std::string_view HInfLAEvaluator::HInfLAEvaluationTask::GetDescription() const {
	using namespace std;
	return "Computing image"sv;
}

template <bool DistanceEstimation>
void HInfLAEvaluator::HInfLAEvaluationTask::Evaluate(GroupedRasterizingInterface &RI, Reference *reference) {
	LAReference<> &Ref = *static_cast<LAReference<>*>(reference);
	bool HighQuality = Global::HighQuality;
	constexpr bool EvaluateDzdc = true;

	size_t MaxIt = Ref.MaxIt;
	if (Ref.IsPeriodic || parameters.MaxIterations < Ref.MaxIt) {
		MaxIt = parameters.MaxIterations;
	}

	if constexpr (DistanceEstimation) Global::MaxIt = 1ull << 48;
	else Global::MaxIt = MaxIt;

	using vreal = VHRReal<VSize>;
	using mask = VBool<VSize>;
	mask ValidMask;
	using vcomplex = VHRComplex<VSize>;

	vreal dcr, dci;

	HRReal PixScale = 2.0_hr * FContext.EvalLocation.HalfH / FContext.ImageHeight;

	size_t ValidCount;
	HRReal dcrs[VSize], dcis[VSize];

	while ((ValidCount = RI.GetCoordinate(dcrs, dcis))) {
		ValidMask = VInt64IntegerSequence<VSize>() < VInt64<VSize>(ValidCount);
		dcr = ArrayToVector<VSize>(dcrs);
		dci = ArrayToVector<VSize>(dcis);
		vcomplex dc = { dcr, dci };

		EvaluationContext<HRReal, EvaluateDzdc> HRContext(Cancelled, Ref, PixScale, ValidMask, dc, MaxIt);

		if (Ref.DirectEvaluate) {
			EvaluationContext<SRReal, EvaluateDzdc> SRContext(HRContext);

			SRContext.EvaluateDirectly();
			HRContext.Iteration = SRContext.Iteration;
			HRContext.z = SRContext.z;
			if constexpr (EvaluateDzdc) HRContext.dzdc = vcomplex(SRContext.dzdc) * SRContext.InvScalingFactor;
			goto end;
		}

		if (Ref.UseAT && !(ChebyshevNorm(dc) > vreal(Ref.AT.ThresholdC))) {
			size_t ATMaxIt = (Ref.UseAT) ? (MaxIt - 1) / Ref.AT.StepLength + 1 : 0;

			EvaluationContext<SRReal, EvaluateDzdc> SRContext(HRContext);
			SRContext.c = VSRComplex<VSize>(dc * vcomplex(Ref.AT.CCoeff)) + VSRComplex<VSize>(Ref.AT.RefC);
			SRContext.MaxIt = ATMaxIt;
			SRContext.EvaluateAT(Ref.AT.SqrEscapeRadius);

			HRContext.z = vcomplex(SRContext.z) * vcomplex(Ref.AT.InvZCoeff);
			HRContext.dz = HRContext.z;
			if constexpr (EvaluateDzdc) HRContext.dzdc = vcomplex(SRContext.dzdc) * vcomplex(Ref.AT.CCoeff) * vcomplex(Ref.AT.InvZCoeff) * SRContext.InvScalingFactor;
			HRContext.Iteration = _mm256_mul_epu32(SRContext.Iteration, _mm256_set1_epi64x(Ref.AT.StepLength)); // FIXME

			HRContext.TrySampleJulia(HRContext.ValidMask);
			HRContext.SamplingEnabledMask = AndNot(HRContext.SamplingDisablingMask, HRContext.SamplingEnabledMask);
			if (!HRContext.ValidMask) {
				goto end;
			}
		}

		while (HRContext.CurrentLAStage) {
			if (Ref.LAStages[HRContext.CurrentLAStage - 1].UseDoublePrecision) {
				goto DPEvaluation;
			}
			HRContext.CurrentLAStage--;
			HRContext.EvaluateLA();
			HRContext.SamplingEnabledMask = AndNot(HRContext.SamplingDisablingMask, HRContext.SamplingEnabledMask);

			if (!HRContext.ValidMask) {
				goto end;
			}
		}

		if (Ref.DoublePrecisionPT) {
			goto DPEvaluation;
		}
		HRContext.EvaluatePT();
		goto end;

		DPEvaluation :
		{
			EvaluationContext<SRReal, EvaluateDzdc> SRContext(HRContext);

			while (SRContext.CurrentLAStage) {
				SRContext.CurrentLAStage--;
				SRContext.EvaluateLA();
				SRContext.SamplingEnabledMask = AndNot(SRContext.SamplingDisablingMask, SRContext.SamplingEnabledMask);

				if (!SRContext.ValidMask) {
					goto SREnd;
				}
			}

			SRContext.EvaluatePT();

			SREnd:
			HRContext.Iteration = SRContext.Iteration;
			HRContext.z = Select(HRContext.ValidMask, VHRComplex<VSize>(SRContext.z), HRContext.z);
			if constexpr (EvaluateDzdc) HRContext.dzdc = Select(HRContext.ValidMask, vcomplex(SRContext.dzdc) * SRContext.InvScalingFactor, HRContext.dzdc);
			HRContext.SamplingIteration = Select(SRContext.SamplingUpdatedMask, SRContext.SamplingIteration, HRContext.SamplingIteration);
			HRContext.SamplingZ = Select(SRContext.SamplingUpdatedMask, VHRComplex<VSize>(SRContext.SamplingZ), HRContext.SamplingZ);
			if constexpr (EvaluateDzdc) HRContext.SamplingDzdc = Select(SRContext.SamplingUpdatedMask, vcomplex(SRContext.SamplingDzdc) * SRContext.InvScalingFactor, HRContext.SamplingDzdc);

			HRContext.SamplingMask = SRContext.SamplingMask;
			HRContext.MixedSamplingMask = SRContext.MixedSamplingMask;
		}

	end:
		vreal FinalMagDzdc = vreal(-1.0);
		vreal FinalMagnitude = std::norm(HRContext.z);
		i64vec4 Iteration = HRContext.Iteration;

		if constexpr (EvaluateDzdc) {
			FinalMagDzdc = sqrt(std::norm(HRContext.dzdc));
		}
		VSRReal<VSize> VSamplingLogInvDE;
		VSRReal<VSize> VSamplingIteration;
		VSRReal<VSize> VSamplingAcceptance;
		if (HighQuality) {
			SRReal SamplingIteration[4];
			SRReal SamplingLogInvDE[4];
			SRReal SamplingAcceptance[4];

			VHRReal<VSize> MagSamplingDzdc = sqrt(std::norm(HRContext.SamplingDzdc));
			for (size_t i = 0; i < VectorSize; i++) {
				if (HRContext.SamplingMask[i]) {
					HRComplex z = HRComplex(HRContext.SamplingZ.real()[i], HRContext.SamplingZ.imag()[i]);
					HRReal MagDzdc = HRReal(MagSamplingDzdc.Mantissa[i], MagSamplingDzdc.Exponent[i]);
					JuliaPixel Pixel = Ref.Julia.Sample(z, MagDzdc * PixScale);

					SamplingLogInvDE[i] = Pixel.LogInvDE;
					SamplingIteration[i] = Pixel.Iteration;
					SamplingAcceptance[i] = Pixel.Weight;
				}
			}
			VSamplingLogInvDE = VSRReal<VSize>(SamplingLogInvDE[0], SamplingLogInvDE[1], SamplingLogInvDE[2], SamplingLogInvDE[3]);
			VSamplingIteration = VSRReal<VSize>(SamplingIteration[0], SamplingIteration[1], SamplingIteration[2], SamplingIteration[3]);
			VSamplingAcceptance = VSRReal<VSize>(SamplingAcceptance[0], SamplingAcceptance[1], SamplingAcceptance[2], SamplingAcceptance[3]);

		}
		SRReal Results[VSize];
		VSRReal<VSize> Result;
		if constexpr (DistanceEstimation) {
			FinalMagnitude = sqrt(FinalMagnitude);
			dvec4 LogMagnitude = log(FinalMagnitude);
			vreal InvDE = (FinalMagDzdc + FinalMagDzdc) / (FinalMagnitude * vreal(LogMagnitude));

			VSRReal<VSize> LogVal = log(std::norm(HRContext.SamplingDzdc) / std::norm(HRContext.SamplingZ)) * 0.5;

			dvec4 LogInvDE = log(InvDE);
			LogVal += VSamplingLogInvDE;
			LogVal = Select(HRContext.MixedSamplingMask, LogInvDE * (VSRReal<VSize>(1.0) - VSamplingAcceptance) + LogVal * VSamplingAcceptance, LogVal);
			LogInvDE = Select(HRContext.SamplingMask, LogVal, LogInvDE);

			LogInvDE = log(HRExp<VHRReal<VSize>>(LogInvDE) + VHRReal<VSize>(1));
			LogInvDE *= 8;

			Result = Select(Iteration >= i64vec4(MaxIt), dvec4(std::numeric_limits<double>::infinity()), LogInvDE);
		} else {
			VSRReal<VSize> IterDouble = Iteration.ToDouble();
			IterDouble = Select(FinalMagnitude > vreal(4096.0), IterDouble + VSRReal<VSize>(log2(log2(4096.0))) - log2(log2(FinalMagnitude)), IterDouble);

			if (HighQuality) {
				VSRReal<VSize> SRRSamplingIteration = HRContext.SamplingIteration.ToDouble();
				SRRSamplingIteration += VSamplingIteration;

				SRRSamplingIteration = Select(HRContext.MixedSamplingMask, IterDouble * (VSRReal<VSize>(1.0) - VSamplingAcceptance) + SRRSamplingIteration * VSamplingAcceptance, SRRSamplingIteration);
				IterDouble = Select(HRContext.SamplingMask, SRRSamplingIteration, IterDouble);
			}
			Result = IterDouble;
		}
		for (size_t i = 0; i < VSize; i++) {
			Results[i] = Result[i];
		}
		RI.WriteResults(Results);
	}
}

void HInfLAEvaluator::HInfLAEvaluationTask::JuliaRenderingContext::RenderJulia(const EvaluationParameters &parameters, LAReference<> &reference, bool IsMainThread) {
	if (!Global::HighQuality) return;
	if (IsMainThread) {
		reference.Radius = FContext.CurrentLocation.HalfH;
		if (reference.Julia.Maps) {
			for (Int i = 0; i < reference.Julia.MipmapLevel; i++) {
				delete[]reference.Julia.Maps[i].Data;
			}
			delete[]reference.Julia.Maps;
			reference.Julia.Maps = nullptr;
		}

		SRReal PixelDensity = 1.0 / PixelSize;
		OriginOffset = ceil(PixelDensity * 2 + 0.5);
		reference.Julia.MinMagZ = sqrt(FContext.CurrentLocation.HalfH) * 0x1.0p24; // FIXME
		reference.Julia.MaxMagZ = M_E;
		reference.Julia.PixelDensity = PixelDensity;
		reference.Julia.PixelSize = PixelSize;

		if (reference.Julia.MinMagZ > 0x1.0p-4) {
			Failed = true;
			return;
		}

		SRReal MaxX = -log(reference.Julia.MinMagZ * 0x1.0p-12);
		IntPix MaxXPixel = ceil(MaxX * PixelDensity + 0.5);
		Width = OriginOffset + MaxXPixel;

		reference.Julia.Maps = new JuliaMipmap[MipmapLevel];
		IntPix MipmapMaxXPixel = MaxXPixel, MipmapOriginOffset = OriginOffset, MipmapHeight = Height;
		SRReal MipmapPixelDensity = PixelDensity, MipmapPixelSize = PixelSize;
		for (size_t i = 0; i < MipmapLevel; i++) {
			JuliaMipmap &Mipmap = reference.Julia.Maps[i];
			Mipmap.Width = MipmapOriginOffset + MipmapMaxXPixel;
			Mipmap.Height = MipmapHeight;
			Mipmap.OriginOffset = MipmapOriginOffset;
			Mipmap.PixelDensity = MipmapPixelDensity;
			Mipmap.PixelSize = MipmapPixelSize;

			size_t AllocationWidth = (Mipmap.Width + (PixelGroupSize - 1)) & ~(PixelGroupSize - 1);

			Mipmap.Data = new JuliaPixel[AllocationWidth * Mipmap.Height];

			if (MipmapMaxXPixel <= 2) {
				MipmapLevel = i + 1;
				break;
			}

			MipmapOriginOffset = (MipmapOriginOffset + 1) >> 1;
			MipmapMaxXPixel = (MipmapMaxXPixel + 1) >> 1;
			MipmapHeight = (MipmapHeight + 1) >> 1;
			MipmapPixelDensity = MipmapPixelDensity * 0.5;
			MipmapPixelSize = MipmapPixelSize * 2.0;
		}
		reference.Julia.MipmapLevel = MipmapLevel;

		Started = true;
	} else {
		while (!Started) { // TODO: Improve synchronization mechanism
			if (Failed) return;
		}
	}

	std::random_device random_device;
	std::default_random_engine random(random_device());
	std::uniform_real_distribution<SRReal> distribution(0.0, SubPixelSize);

	for (IntPix PixX = Column.fetch_add(PixelGroupSize); PixX < Width; PixX = Column.fetch_add(PixelGroupSize)) {
		for (IntPix PixY = 0; PixY < Height; PixY += PixelGroupSize) {
			SRReal x = IntPix(PixX - OriginOffset);
			HRReal z0r[VectorSize];
			HRReal z0i[VectorSize];

			SRReal y = PixY;

			SRReal SumLogInvDE = 0.0;
			SRReal SumIteration = 0.0;
			for (size_t SubPixX = 0; SubPixX < SuperSample; SubPixX += 2) {
				for (size_t SubPixY = 0; SubPixY < SuperSample; SubPixY += 2) {
					for (size_t i = 0; i < 4; i++) {
						SRReal SampleX = (x + distribution(random) + (SubPixX + ((i & 1) ? 1 : 0)) * SubPixelSize) * PixelSize;
						SRReal SampleY = (y + distribution(random) + (SubPixY + ((i & 2) ? 1 : 0)) * SubPixelSize) * PixelSize;

						HRReal Magnitude = HRExp<HRReal>(-SampleX);
						SRReal Argument = SampleY;
						z0r[i] = Magnitude * cos(Argument);
						z0i[i] = Magnitude * sin(Argument);
					}

					VHRComplex<VSize> z0 = VHRComplex<VSize>(VHRReal<VSize>(z0r[0], z0r[1], z0r[2], z0r[3]), VHRReal<VSize>(z0i[0], z0i[1], z0i[2], z0i[3]));

					VHRReal<VSize> InvDE;
					VSRReal<VSize> Iteration;

					Iteration = EvaluateJulia(parameters, reference, z0, InvDE);

					VHRReal<VSize> MagnitudeZ = sqrt(std::norm(z0));
					InvDE *= MagnitudeZ;

					if (PixelGroupSize == 1) {
						for (size_t i = 0; i < VectorSize; i++) {
							SumLogInvDE += log(HRReal(InvDE.Mantissa[i], InvDE.Exponent[i]));
							SumIteration += Iteration[i];
						}
					} else {
						for (size_t i = 0; i < 2; i++) {
							for (size_t j = 0; j < 2; j++) {
								JuliaPixel &Pixel = reference.Julia.Maps[0].Data[size_t(Height) * (PixX + i) + (PixY + j)];
								Pixel.LogInvDE = log(HRReal(InvDE.Mantissa[i * 2 + j], InvDE.Exponent[i * 2 + j]));
								Pixel.Iteration = Iteration[i];
							}
						}
					}
				}
			}

			if (PixelGroupSize == 1) {
				JuliaPixel &Pixel = reference.Julia.Maps[0].Data[size_t(Height) * PixX + PixY];
				Pixel.LogInvDE = SumLogInvDE * (1.0 / SuperSampleSqr);
				Pixel.Iteration = SumIteration * (1.0 / SuperSampleSqr);
			}
		}
	}

	if (!IsMainThread) {
		while (!Finished); // TODO: Improve synchronization mechanism
		return;
	}

	size_t Level = 1;
	for (; Level < MipmapLevel; Level++) {
		JuliaMipmap &Mipmap = reference.Julia.Maps[Level];
		JuliaMipmap &PrevMipmap = reference.Julia.Maps[Level - 1];
		if (PrevMipmap.Height == 1) break;

		IntPix x = 0, PrevMipmapX = 0;

		Mipmap.FirstPixelWidth = PrevMipmap.FirstPixelWidth * 0.5;

		if (PrevMipmap.OriginOffset & 1) {
			for (IntPix y = 0; y < Mipmap.Height; y++) {
				JuliaPixel Pixel0 = PrevMipmap.Data[y * 2];
				JuliaPixel Pixel1 = PrevMipmap.Data[y * 2 + 1];
				SRReal SumLogInvDE = Pixel0.LogInvDE * Pixel0.Weight + Pixel1.LogInvDE * Pixel1.Weight;
				SRReal SumIteration = Pixel0.Iteration * Pixel0.Weight + Pixel1.Iteration * Pixel1.Weight;
				SRReal Weight = Pixel0.Weight + Pixel1.Weight;
				Mipmap.Data[y].LogInvDE = SumLogInvDE / Weight;
				Mipmap.Data[y].Iteration = SumIteration / Weight;
				Mipmap.Data[y].Weight = Weight;
			}
			x = 1; PrevMipmapX = 1;
		} else {
			Mipmap.FirstPixelWidth += 0.5;
		}

		for (; PrevMipmapX < PrevMipmap.Width - 1; x++, PrevMipmapX += 2) {
			for (IntPix y = 0; y < Mipmap.Height; y++) {
				JuliaPixel Pixel0 = PrevMipmap.Data[(PrevMipmapX + 0) * PrevMipmap.Height + (y * 2 + 0)];
				JuliaPixel Pixel1 = PrevMipmap.Data[(PrevMipmapX + 0) * PrevMipmap.Height + (y * 2 + 1)];
				JuliaPixel Pixel2 = PrevMipmap.Data[(PrevMipmapX + 1) * PrevMipmap.Height + (y * 2 + 0)];
				JuliaPixel Pixel3 = PrevMipmap.Data[(PrevMipmapX + 1) * PrevMipmap.Height + (y * 2 + 1)];

				SRReal SumLogInvDE = Pixel0.LogInvDE * Pixel0.Weight + Pixel1.LogInvDE * Pixel1.Weight;
				SumLogInvDE += Pixel2.LogInvDE * Pixel2.Weight + Pixel3.LogInvDE * Pixel3.Weight;
				SRReal SumIteration = Pixel0.Iteration * Pixel0.Weight + Pixel1.Iteration * Pixel1.Weight;
				SumIteration += Pixel2.Iteration * Pixel2.Weight + Pixel3.Iteration * Pixel3.Weight;

				SRReal Weight = Pixel0.Weight + Pixel1.Weight + Pixel2.Weight + Pixel3.Weight;
				Mipmap.Data[x * Mipmap.Height + y].LogInvDE = SumLogInvDE / Weight;
				Mipmap.Data[x * Mipmap.Height + y].Iteration = SumIteration / Weight;
				Mipmap.Data[x * Mipmap.Height + y].Weight = Weight;
			}
		}

		Mipmap.LastPixelWidth = PrevMipmap.LastPixelWidth * 0.5;

		if (PrevMipmapX == PrevMipmap.Width - 1) {
			for (IntPix y = 0; y < Mipmap.Height; y++) {
				JuliaPixel Pixel0 = PrevMipmap.Data[(PrevMipmapX + 0) * PrevMipmap.Height + (y * 2 + 0)];
				JuliaPixel Pixel1 = PrevMipmap.Data[(PrevMipmapX + 0) * PrevMipmap.Height + (y * 2 + 1)];
				SRReal SumLogInvDE = Pixel0.LogInvDE * Pixel0.Weight + Pixel1.LogInvDE * Pixel1.Weight;
				SRReal SumIteration = Pixel0.Iteration * Pixel0.Weight + Pixel1.Iteration * Pixel1.Weight;
				SRReal Weight = Pixel0.Weight + Pixel1.Weight;
				Mipmap.Data[x * Mipmap.Height + y].LogInvDE = SumLogInvDE / Weight;
				Mipmap.Data[x * Mipmap.Height + y].Iteration = SumIteration / Weight;
				Mipmap.Data[x * Mipmap.Height + y].Weight = Weight;
			}
		} else {
			Mipmap.LastPixelWidth += 0.5;
		}

		Mipmap.InvFirstPixelsCenterDistance = 1.0 / (Mipmap.FirstPixelWidth * 0.5 + 0.5);
		Mipmap.InvLastPixelsCenterDistance = 1.0 / (Mipmap.LastPixelWidth * 0.5 + 0.5);
	}
	for (; Level < MipmapLevel; Level++) {
		JuliaMipmap &Mipmap = reference.Julia.Maps[Level];
		JuliaMipmap &PrevMipmap = reference.Julia.Maps[Level - 1];

		IntPix x = 0, PrevMipmapX = 0;

		Mipmap.FirstPixelWidth = PrevMipmap.FirstPixelWidth * 0.5;

		if (PrevMipmap.OriginOffset & 1) {
			Mipmap.Data[0] = PrevMipmap.Data[0];
			x = 1; PrevMipmapX = 1;
		} else {
			Mipmap.FirstPixelWidth += 0.5;
		}

		for (; PrevMipmapX < PrevMipmap.Width - 1; x++, PrevMipmapX += 2) {
			JuliaPixel Pixel0 = PrevMipmap.Data[PrevMipmapX];
			JuliaPixel Pixel1 = PrevMipmap.Data[PrevMipmapX + 1];
			SRReal SumLogInvDE = Pixel0.LogInvDE * Pixel0.Weight + Pixel1.LogInvDE * Pixel1.Weight;
			SRReal SumIteration = Pixel0.LogInvDE * Pixel0.Weight + Pixel1.LogInvDE * Pixel1.Weight;
			SRReal Weight = Pixel0.Weight + Pixel1.Weight;
			Mipmap.Data[x].LogInvDE = SumLogInvDE / Weight;
			Mipmap.Data[x].LogInvDE = SumIteration / Weight;
			Mipmap.Data[x].Weight = Weight;
		}

		Mipmap.LastPixelWidth = PrevMipmap.LastPixelWidth * 0.5;

		if (PrevMipmapX == PrevMipmap.Width - 1) {
			Mipmap.Data[x] = PrevMipmap.Data[PrevMipmapX];
		} else {
			Mipmap.LastPixelWidth += 0.5;
		}

		Mipmap.InvFirstPixelsCenterDistance = 1.0 / (Mipmap.FirstPixelWidth * 0.5 + 0.5);
		Mipmap.InvLastPixelsCenterDistance = 1.0 / (Mipmap.LastPixelWidth * 0.5 + 0.5);
	}
	reference.Julia.Valid = true;
	Finished = true;
}

VSRReal<VectorSize> HInfLAEvaluator::HInfLAEvaluationTask::JuliaRenderingContext::EvaluateJulia(const EvaluationParameters &parameters, LAReference<> &reference, VHRComplex<VectorSize> dz, VHRReal<VectorSize> &InvDE) {
	size_t MaxIt = reference.MaxIt;
	if (reference.IsPeriodic || parameters.MaxIterations < reference.MaxIt) {
		MaxIt = parameters.MaxIterations;
	}

	using mask = mask64x4;
	using HREvaluationContext = EvaluationContext<HRReal, true, true>;
	using SREvaluationContext = EvaluationContext<SRReal, true, true>;

	bool Cancelled = false; // FIXME
	HREvaluationContext HRContext(Cancelled, reference, 0.0, mask(true), dz, VHRReal<VSize>(0.0), MaxIt);

	if (reference.DirectEvaluate) {
		SREvaluationContext SRContext(HRContext);

		SRContext.EvaluateDirectly();
		HRContext.Iteration = SRContext.Iteration;

		HRContext.z = SRContext.z;
		HRContext.dzdz = VHRComplex<VSize>(SRContext.dzdz) * SRContext.InvScalingFactor;
		HRContext.dzdc = VHRComplex<VSize>(SRContext.dzdc) * SRContext.InvScalingFactor;
		goto end;
	}

	while (HRContext.CurrentLAStage) {
		if (reference.LAStages[HRContext.CurrentLAStage - 1].UseDoublePrecision) {
			goto SREvaluation;
		}
		HRContext.CurrentLAStage--;
		HRContext.EvaluateLA();
	}

	if (reference.DoublePrecisionPT) {
		goto SREvaluation;
	}
	HRContext.EvaluatePT();
	goto end;

SREvaluation:
	{
		SREvaluationContext SRContext(HRContext);

		while (SRContext.CurrentLAStage) {
			SRContext.CurrentLAStage--;
			SRContext.EvaluateLA();
		}

		SRContext.EvaluatePT();

		HRContext.Iteration = SRContext.Iteration;

		HRContext.z = SRContext.z ;
		HRContext.dzdz = VHRComplex<VSize>(SRContext.dzdz) * SRContext.InvScalingFactor;
		HRContext.dzdc = VHRComplex<VSize>(SRContext.dzdc) * SRContext.InvScalingFactor;
	}

end:
	VHRReal<VSize> FinalMagDzdz = VHRReal<VSize>(-1.0);
	VHRReal<VSize> FinalMagnitude = std::norm(HRContext.z);
	i64vec4 Iteration = HRContext.Iteration;

	FinalMagDzdz = std::norm(HRContext.dzdz);

	VSRReal<VSize> IterDouble = Iteration.ToDouble();
	IterDouble = Select(FinalMagnitude != VSRReal<VSize>(-1.0), IterDouble + VSRReal<VSize>(log2(log2(4096.0))) - log2(log2(FinalMagnitude)), IterDouble);

	FinalMagnitude = sqrt(FinalMagnitude);
	FinalMagDzdz = sqrt(FinalMagDzdz);

	dvec4 LogMagnitude = log(FinalMagnitude);
	InvDE = (FinalMagDzdz + FinalMagDzdz) / (FinalMagnitude * VHRReal<VSize>(LogMagnitude));

	return IterDouble;
}