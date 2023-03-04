#include "Includes.h"

Reference::~Reference() {}
Evaluator::~Evaluator() {}
Evaluator::FeatureFinder::~FeatureFinder() {}

void Evaluator::FreeReference(Reference *reference) {
	delete reference;
}

std::wstring_view Evaluator::Feature::Name() {
	return std::wstring_view();
}

std::wstring_view Evaluator::Feature::Information() {
    return std::wstring_view();
}

IntIter Evaluator::Feature::ItLimForZoomLevel(HRReal HalfH) {
	return 0;
}

bool Evaluator::Feature::GetRadius(HRReal &targetRadius) {
	return false;
}

bool Evaluator::Feature::CanLocatePrecisely() {
	return false;
}

std::string_view Evaluator::FeatureFinder::PreciseLocatingTask::GetDetailedProgress() {
	return std::string_view();
}

HInfLAEvaluator::FeatureFinder::PreciseLocatingTask *Evaluator::FeatureFinder::CreatePreciseLocatingTask(Feature *feature, size_t Precision, const Coordinate &centerCoordinate) {
	return nullptr;
}

bool Evaluator::FeatureFinder::LocatePrecisely(Feature *feature, size_t Precision, Coordinate &featureCoordinate, const Coordinate &centerCoordinate) {
	return false;
}

Evaluator::Feature::~Feature() {}

Evaluator::FeatureFinder *Evaluator::GetFeatureFinder() {
	return nullptr;
}
size_t Evaluator::SaveReference(unused std::ostream &stream, unused const Reference *reference) {
	return 0;
}
size_t Evaluator::LoadReference(unused std::istream &stream, unused Reference *&reference) {
	return 0;
}


std::string_view SimpleEvaluator::StandardReferenceTask::GetDescription() const {
	using namespace std;
	return "Computing reference"sv;
}

void SimpleEvaluator::StandardReferenceTask::Execute() {
	reference = evaluator->GenerateReference(parameters);
}

std::string_view SimpleEvaluator::StandardEvaluationTask::GetDescription() const {
	using namespace std;
	return "Computing image"sv;
}

void SimpleEvaluator::StandardEvaluationTask::Execute() {
	evaluator->Evaluate(parameters, *rasterizer, reference);
}

bool SimpleEvaluator::StandardEvaluationTask::GetProgress(SRReal &Numerator, SRReal &Denoninator) const {
	const ProgressTrackable *progressTrackable = dynamic_cast<const ProgressTrackable *>(rasterizer);
	if (progressTrackable) {
		return progressTrackable->GetProgress(Numerator, Denoninator);
	}
	return false;
}

Reference *SimpleEvaluator::GenerateReference(unused const EvaluationParameters &parameters) {
	return nullptr;
}

Evaluator::ReferenceTask *SimpleEvaluator::CreateReferenceTask(const EvaluationParameters &parameters) {
	return new StandardReferenceTask(this, parameters);
}

Task *SimpleEvaluator::CreateEvaluationTask(const EvaluationParameters &parameters, PixelManager *rasterizer, Reference *reference) {
	return new StandardEvaluationTask(this, parameters, rasterizer, reference);
}

template<FractalTypeEnum FractalType>
PerturbationEvaluator<FractalType>::PerturbationEvaluator() {
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	ReserveSize = (statex.ullTotalPhys) & ~(CommitSize - 1);
}

template<FractalTypeEnum FractalType>
Reference *PerturbationEvaluator<FractalType>::GenerateReference(const EvaluationParameters &parameters) {
	if (FContext.CurrentLocation.HalfH > 0x1.0p-896_hr) {
		return _GenerateReference<double>(parameters);
	} else {
		return _GenerateReference<FExpDouble>(parameters);
	}
}
template<FractalTypeEnum FractalType>
template <typename real2>
Reference *PerturbationEvaluator<FractalType>::_GenerateReference(const EvaluationParameters &parameters) {
	Global::MaxIt = Global::ItLim;

	PTReference<real2, FractalType> &Ref = *new PTReference<real2, FractalType>(ReserveSize);

	Ref.AbsolutePrecision = pow2(-double(parameters.CenterCoordinate.X.get_prec())) * 16.0_hr;

	size_t MaxIt = Ref.MaxIt = Global::ItLim;

	HPReal cr = parameters.CenterCoordinate.X, ci = parameters.CenterCoordinate.Y;
	HPReal zr = 0.0, zi = 0.0;
	size_t i = 0;
	HRReal zrl = zr, zil = zi;

	for (; i <= MaxIt; i++) {
		Ref.Refr.Append(real2(zrl));
		Ref.Refi.Append(real2(zil));
		HPReal newzr = zr * zr - zi * zi + cr;
		HPReal newzi = zr * zi;
		if constexpr (FractalType == FractalTypeEnum::BurningShip) newzi = abs(newzi);
		newzi = newzi + newzi + ci;
		if constexpr (FractalType == FractalTypeEnum::Tricorn) newzi = -newzi;
		zrl = newzr;
		zil = newzi;
		HRReal Magnitude = zrl * zrl + zil * zil;
		if (Magnitude > 64.0) {
			i++;
			break;
		}
		zr = newzr;
		zi = newzi;
	}
	Ref.RefIt = i - 1;
	return &Ref;
}

double diffabs(double X, double x) {
	if (X >= 0) {
		if (X + x >= 0) { return x; } else { return -(2 * X + x); }
	} else {
		if (X + x > 0) { return 2 * X + x; } else { return -x; }
	}
}
inline dvec4 diffabs(dvec4 X, dvec4 x) {
	dvec4 Sum = X + x;
	dvec4 SignSum = Sum & dvec4(-0.0);
	mask64x4 DifferentSign = (Sum ^ X).ymm;
	return Select(DifferentSign, X + Sum, x) ^ SignSum;
}
inline DExpVec4 diffabs(DExpVec4 X, DExpVec4 x) {
	DExpVec4 Sum = X + x;
	dvec4 SignSum = Sum.Mantissa & dvec4(-0.0);
	mask64x4 DifferentSign = (Sum.Mantissa ^ X.Mantissa).ymm;
	DExpVec4 Result = Select(DifferentSign, X + Sum, x);
	Result.Mantissa ^= SignSum;
	return Result;
}

template<FractalTypeEnum FractalType>
void PerturbationEvaluator<FractalType>::Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *reference) {
	_PTReference &Ref = *static_cast<_PTReference *>(reference);
	if (Ref.Type == typeid(double)) {
		_Evaluate<double>(rasterizer, reference);
	} else {
		_Evaluate<FExpDouble>(rasterizer, reference);
	}
}
template<FractalTypeEnum FractalType>
template <typename real2>
void PerturbationEvaluator<FractalType>::_Evaluate(PixelManager &rasterizer, Reference *reference) {
	GroupedRasterizingInterface &RI = rasterizer.GetGroupedInterface(VSize);
	PTReference<real2, FractalType> &Ref = *static_cast<PTReference<real2, FractalType> *>(reference);
	size_t MaxIt = std::min(Ref.MaxIt, Global::ItLim);
	Global::MaxIt = MaxIt;

	using vreal2 = vec4<real2>;
	using mask = mask64x4;

	mask ValidMask;

	VHRReal<VSize> Dcr, Dci;
	HRReal dcrs[VSize], dcis[VSize];
	size_t ValidCount;
	while ((ValidCount = RI.GetCoordinate(dcrs, dcis))) {
		ValidMask = VInt64IntegerSequence<VSize>() < VInt64<VSize>(ValidCount);
		Dcr = ArrayToVector<VSize>(dcrs);
		Dci = ArrayToVector<VSize>(dcis);
		vreal2 dcr = vreal2(Dcr), dci = vreal2(Dci);
		vreal2 dr = 0.0, di = 0.0;
		i64vec4 Iteration = 0;
		size_t j = 0;

		vreal2 FinalMagnitude = -1.0;
		mask ActiveMask = ValidMask;
		mask PauseMask = false;
		vreal2 InactiveDr, InactiveDi;
		while (true) {
			if (!ActiveMask) {
				if (PauseMask) {
					dr = InactiveDr;
					di = InactiveDi;
					ActiveMask = PauseMask;
					PauseMask = false;
					j = 0;
				} else break;
			}
			real2 Zr, Zi;
			Zr = Ref.Refr[j];
			Zi = Ref.Refi[j];
			vreal2 newdr = (dr * Zr - di * Zi) * real2(2.0) + dr * dr - di * di + dcr;
			vreal2 newdi = dr * Zi + di * Zr + dr * di;
			if constexpr (FractalType == FractalTypeEnum::BurningShip) newdi = diffabs(Zr * Zi, newdi);
			newdi = newdi * real2(2.0) + dci;
			if constexpr (FractalType == FractalTypeEnum::Tricorn) newdi = -newdi;
			j++;

			dr = newdr;
			di = newdi;
			vreal2 zr = dr + Ref.Refr[j], zi = di + Ref.Refi[j];

			vreal2 Magnitude = zr * zr + zi * zi;
			vreal2 magd = dr * dr + di * di;
			Iteration -= i64vec4(ActiveMask);
			ActiveMask = ActiveMask && (Iteration < i64vec4(Ref.MaxIt));
			mask CompareResult = (Magnitude > vreal2(4096.0)) && ActiveMask;
			if (CompareResult) {
				FinalMagnitude = Select(CompareResult, Magnitude, FinalMagnitude);
				ActiveMask = AndNot(CompareResult, ActiveMask);
			}
			if (!ActiveMask) continue;
			if (j >= Ref.RefIt) {
				dr = Select(PauseMask, InactiveDr, zr);
				di = Select(PauseMask, InactiveDi, zi);
				ActiveMask = PauseMask || ActiveMask;
				PauseMask = false;
				j = 0;
			} else {
				mask CompareResult = (Magnitude < magd *real2(64.0)) && ActiveMask;
				if (CompareResult) {
					InactiveDr = Select(CompareResult, zr, InactiveDr);
					InactiveDi = Select(CompareResult, zi, InactiveDi);
					PauseMask = CompareResult || PauseMask;
					ActiveMask = AndNot(CompareResult, ActiveMask);
				}
			}
		}

		dvec4 IterDouble = Iteration.ToDouble();
		IterDouble = Select(FinalMagnitude != vreal2(-1.0), IterDouble + dvec4(log2(log2(4096.0))) - log2(log2(FinalMagnitude)), IterDouble);

		SRReal Results[VSize];
		for (size_t i = 0; i < VSize; i++) {
			Results[i] = IterDouble[i];
		}
		RI.WriteResults(Results);
	}
	rasterizer.FreeInterface(RI);
}

template class PerturbationEvaluator<FractalTypeEnum::Mandelbrot>;
template class PerturbationEvaluator<FractalTypeEnum::BurningShip>;
template class PerturbationEvaluator<FractalTypeEnum::Tricorn>;

NovaEvaluator::NovaEvaluator() {
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);
	ReserveSize = (statex.ullTotalPhys) & ~(CommitSize - 1);
}

Reference *NovaEvaluator::GenerateReference(const EvaluationParameters &parameters) {
	if (FContext.CurrentLocation.HalfH > 0x1.0p-896_hr) {
		return _GenerateReference<double>(parameters.CenterCoordinate);
	} else {
		return _GenerateReference<FExpDouble>(parameters.CenterCoordinate);
	}
}

std::complex<HPReal> div(const std::complex<HPReal> &a, const std::complex<HPReal> &b) {
	return a * std::complex<HPReal>(b.real(), -b.imag()) / norm(b);
}

template <typename real2>
Reference *NovaEvaluator::_GenerateReference(const Coordinate &coordinate) {
	using complexh = std::complex<HPReal>;
	using complex = std::complex<real2>;

	PTReference<real2> &Ref = *new PTReference<real2>;
	size_t MaxIt = Ref.MaxIt = Global::ItLim;


	Ref.AbsolutePrecision = pow2(-double(coordinate.X.get_prec())) * 16.0_hr;

	Ref.Zr = ReservedMemoryArray<real2>(ReserveSize);
	Ref.Zsqr = ReservedMemoryArray<complex>(ReserveSize);
	Ref.Zpow3m1 = ReservedMemoryArray<complex>(ReserveSize);
	Ref.Zm1 = ReservedMemoryArray<complex>(ReserveSize);

	complexh C = { coordinate.X, coordinate.Y };
	complexh Z = { 1.0, 0.0 };
	complex z = { 1.0, 0.0 };
	Ref.Refc = std::complex<HRReal>(C);

	size_t i = 0;
	Ref.Zr.Append(z.real());
	Ref.Zsqr.Append(z * z);
	Ref.Zpow3m1.Append(real2(0.0));
	Ref.Zm1.Append(z - real2(1.0));

	complexh ZSQR = Z * Z;
	for (; i <= MaxIt; i++) {
		Z = Z - div((Z * ZSQR - 1_hp), (3_hp * ZSQR)) + C;
		ZSQR = Z * Z;
		z = std::complex<HRReal>(Z);
		complex zpow3m1;
		zpow3m1 = std::complex<HRReal>(ZSQR * Z - 1_hp);
		real2 zrm1 = real2(HRReal(Z.real() - 1));
		complex zm1 = complex(zrm1, z.imag());


		Ref.Zr.Append(z.real());
		Ref.Zsqr.Append(z * z);
		Ref.Zpow3m1.Append(zpow3m1);
		Ref.Zm1.Append(zm1);
		HRReal Magnitude = norm(z);
		if (Magnitude <= real2(0.0)) {
			i++;
			break;
		}
	}
	Ref.RefIt = i - 1;
	return &Ref;
}

void NovaEvaluator::Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *reference) {
	_PTReference &Ref = *static_cast<_PTReference *>(reference);
	if (Ref.Type == typeid(double)) {
		_Evaluate<double>(rasterizer, reference);
	} else {
		_Evaluate<FExpDouble>(rasterizer, reference);
	}
}

template<typename T>
__forceinline std::complex<T> operator^(std::complex<T> a, size_t b) {
	switch (b) {
		case 1: return a;
		case 2: return a * a;
		case 3: return a * a * a;
		case 4: return (a * a) ^ 2;
		default:
			std::complex<T> Result = a;
			for (size_t i = 1; i < b; i++) {
				Result *= a;
			}
			return Result;
	}
}

template<typename real2>
void NovaEvaluator::_Evaluate(PixelManager &rasterizer, Reference *reference) {
	using complex = std::complex<real2>;
	RasterizingInterface &RI = rasterizer.GetInterface();
	//bool intcheck = false;
	PTReference<real2> &Ref = *static_cast<PTReference<real2> *>(reference);
	size_t MaxIt = std::min(Ref.MaxIt, Global::ItLim);
	Global::MaxIt = MaxIt;
	constexpr double EffectiveBailoutRadius = 1.0 / (1ull << 48);

	HRReal Dcr, Dci;
	while (RI.GetCoordinate(Dcr, Dci)) {
		complex dc = { real2(Dcr), real2(Dci) };
		complex c = dc + Ref.Refc;
		complex dz = real2(0.0);
		complex dzdc = real2(0.0);
		complex zm1 = real2(0.0);
		complex prevzm1 = real2(0.0);
		int i = 0, j = 0;

		real2 FinalDistance = -1.0;
		double CurrentRootDerivative;
		SRComplex CurrentRootLocationDerivative;
		double BailoutRadius;
		bool NearRoot = false, CurrentRootIsRepulsive = false;
		SRComplex Rootm1;

		SRComplex Rootsm1[3];
		SRComplex RootDerivatives[3];
		double MagRootDerivatives[3];
		SRComplex RootLocationDerivatives[3];
		double RootBailoutRadius[3];
		double HalfRootDistance;

		{
			SRComplex cd = { SRReal(c.real()), SRReal(c.imag()) };
			SRComplex delta0 = 9.0 * cd * cd;
			SRComplex delta1 = -3.0 * delta0 * cd - 13.5;
			SRComplex CC = std::pow((delta1 - sqrt(delta1 * delta1 - delta0 * delta0 * delta0)), 1.0 / 3.0);
			SRComplex xi = SRComplex(-0.5, sqrt(3.0) / 2.0);
			SRComplex CC1 = xi * CC;
			SRComplex CC2 = xi * CC1;
			Rootsm1[0] = cd - (CC + delta0 / CC) * (1.0 / 3.0);
			Rootsm1[1] = cd - (CC1 + delta0 / CC1) * (1.0 / 3.0);
			Rootsm1[2] = cd - (CC2 + delta0 / CC2) * (1.0 / 3.0);

			RootDerivatives[0] = (2.0 / 3.0) - (2.0 / 3.0) / (Rootsm1[0] ^ 3);
			RootDerivatives[1] = (2.0 / 3.0) - (2.0 / 3.0) / (Rootsm1[1] ^ 3);
			RootDerivatives[2] = (2.0 / 3.0) - (2.0 / 3.0) / (Rootsm1[2] ^ 3);

			MagRootDerivatives[0] = norm(RootDerivatives[0]);
			MagRootDerivatives[1] = norm(RootDerivatives[1]);
			MagRootDerivatives[2] = norm(RootDerivatives[2]);

			RootLocationDerivatives[0] = -1.0 / (RootDerivatives[0] - 1.0);
			RootLocationDerivatives[1] = -1.0 / (RootDerivatives[1] - 1.0);
			RootLocationDerivatives[2] = -1.0 / (RootDerivatives[2] - 1.0);


			for (size_t i = 0; i < 3; i++) {
				double InvRootDerivative2_div2 = norm(Rootsm1[i] ^ 4);
				double OneMinusDerivative = 1.0 - sqrt(MagRootDerivatives[i]);
				OneMinusDerivative *= OneMinusDerivative;
				RootBailoutRadius[i] = std::min(MagRootDerivatives[i], OneMinusDerivative) * InvRootDerivative2_div2 * double(1.0 / 0x10000);
				RootBailoutRadius[i] = std::max(RootBailoutRadius[i], EffectiveBailoutRadius);
			}

			HalfRootDistance = norm(Rootsm1[0] - Rootsm1[1]);
			HalfRootDistance = std::min(HalfRootDistance, norm(Rootsm1[1] - Rootsm1[2]));
			HalfRootDistance = std::min(HalfRootDistance, norm(Rootsm1[0] - Rootsm1[2]));
			HalfRootDistance *= 0.25;

			Rootsm1[0] -= 1.0;
			Rootsm1[1] -= 1.0;
			Rootsm1[2] -= 1.0;
		}

		while (i < Ref.MaxIt) {
			complex Z = complex(Ref.Zr[j], Ref.Zm1[j].imag());
			complex Zsqr = Ref.Zsqr[j];
			complex Zpow3m1 = Ref.Zpow3m1[j];
			complex Zpow4mZ = Zpow3m1 * Z;
			complex Z4 = Zsqr * Zsqr;
			complex subex = dz * (Zsqr * real2(3.0) + dz * (Z * real2(3.0) + dz));
			complex zpow3m1 = Zpow3m1 + subex;
			complex zpow3 = Zsqr * Z + subex;
			dzdc = dzdc * zpow3m1 / (real2(1.5) * zpow3) + real2(1.0);

			subex = (Z * real2(2.0) + dz) * Zsqr * dz;
			dz = dz * ((subex + Zpow4mZ) * real2(2.0) - dz) / (real2(3.0) * (subex + Z4)) + dc;
			j++;
			i++;
			complex Zm1 = Ref.Zm1[j];
			zm1 = dz + Zm1;
			Z = complex(Ref.Zr[j], Ref.Zm1[j].imag());
			real2 magzm1 = norm(zm1);
			real2 magd = norm(dz);
			real2 magZ = norm(Z);
			if (magZ * real2(0.5) < magd) {
				complex z = dz + Z;
				prevzm1 = zm1;
				while (i < Ref.MaxIt) {
					complex zsqr = z * z;
					complex zpow3 = zsqr * z;
					dzdc = dzdc * (zpow3 - real2(1.0)) / (real2(1.5) * zpow3) + real2(1.0);
					z = real2(2.0 / 3.0) * z + real2(1.0 / 3.0) / (zsqr)+c;
					i++;
					zm1 = z - real2(1.0);
					SRComplex zm1d = SRComplex(SRReal(zm1.real()), SRReal(zm1.imag()));
					if (NearRoot) {
						real2 distance = norm(zm1 - complex(Rootm1));
						if (distance > HalfRootDistance) NearRoot = false;
						if (CurrentRootIsRepulsive) {
						} else if (distance < BailoutRadius) {
							FinalDistance = distance;
							break;
						}
					} else if (norm(zm1 - prevzm1) < (1.0 / (1ull << 10))) {
						NearRoot = true;
						Rootm1 = Rootsm1[0];
						CurrentRootDerivative = MagRootDerivatives[0];
						CurrentRootLocationDerivative = RootLocationDerivatives[0];
						BailoutRadius = RootBailoutRadius[0];
						if (norm(Rootsm1[1] - zm1d) < norm(Rootm1 - zm1d)) {
							Rootm1 = Rootsm1[1];
							CurrentRootDerivative = MagRootDerivatives[1];
							CurrentRootLocationDerivative = RootLocationDerivatives[1];
							BailoutRadius = RootBailoutRadius[1];
						}
						if (norm(Rootsm1[2] - zm1d) < norm(Rootm1 - zm1d)) {
							Rootm1 = Rootsm1[2];
							CurrentRootDerivative = MagRootDerivatives[2];
							CurrentRootLocationDerivative = RootLocationDerivatives[2];
							BailoutRadius = RootBailoutRadius[2];
						}
						if (CurrentRootDerivative <= 1.0) {
							double distance = norm(zm1d - Rootm1);
							if (distance < BailoutRadius) {
								FinalDistance = distance;
								break;
							}
							CurrentRootIsRepulsive = false;
						} else {
							CurrentRootIsRepulsive = true;
						}
					}
					prevzm1 = zm1;
				}
				break;
			}
			if (magzm1 < magd || j >= Ref.RefIt) {
				dz = zm1;
				j = 0;
			}

			SRComplex zm1d = SRComplex(SRReal(zm1.real()), SRReal(zm1.imag()));
			if (NearRoot) {
				real2 distance = norm(zm1 - complex(Rootm1));
				if (distance > HalfRootDistance) NearRoot = false;
				if (CurrentRootIsRepulsive) {
				} else if (distance < BailoutRadius) {
					FinalDistance = distance;
					break;
				}
			} else if (norm(zm1 - prevzm1) < (1.0 / (1ull << 10))) {
				NearRoot = true;
				Rootm1 = Rootsm1[0];
				CurrentRootDerivative = MagRootDerivatives[0];
				CurrentRootLocationDerivative = RootLocationDerivatives[0];
				BailoutRadius = RootBailoutRadius[0];
				if (norm(Rootsm1[1] - zm1d) < norm(Rootm1 - zm1d)) {
					Rootm1 = Rootsm1[1];
					CurrentRootDerivative = MagRootDerivatives[1];
					CurrentRootLocationDerivative = RootLocationDerivatives[1];
					BailoutRadius = RootBailoutRadius[1];
				}
				if (norm(Rootsm1[2] - zm1d) < norm(Rootm1 - zm1d)) {
					Rootm1 = Rootsm1[2];
					CurrentRootDerivative = MagRootDerivatives[2];
					CurrentRootLocationDerivative = RootLocationDerivatives[2];
					BailoutRadius = RootBailoutRadius[2];
				}
				if (CurrentRootDerivative <= 1.0) {
					double distance = norm(zm1d - Rootm1);
					if (distance < BailoutRadius) {
						FinalDistance = distance;
						break;
					}
					CurrentRootIsRepulsive = false;
				} else {
					CurrentRootIsRepulsive = true;
				}
			}
			prevzm1 = zm1;
		}

		SRReal fracIter = ((FinalDistance != -1.0) ? (log2(FinalDistance) - log2(EffectiveBailoutRadius)) / log2(CurrentRootDerivative) : 0);
		HRReal density = -sqrt(norm(dzdc - complex(CurrentRootLocationDerivative))) / (sqrt(FinalDistance) * log(CurrentRootDerivative));
		SRReal LogDensity = log(density + 1) * 8.0;
		if (FContext.UsingDE) {
			RI.WriteResults((FinalDistance != -1.0) ? LogDensity : std::numeric_limits<double>::infinity());
		} else {
			RI.WriteResults(double(i) - fracIter);
		}
	}
	rasterizer.FreeInterface(RI);
}