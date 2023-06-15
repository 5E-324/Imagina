#include "Includes.h"

void FractalContext::UpdateRelativeCoordinate(HRReal OffsetX, HRReal OffsetY) {
	pixelManager.UpdateRelativeCoordinate(OffsetX, OffsetY);
	RenderLocation.X += OffsetX;
	RenderLocation.Y += OffsetY;
	EvalLocation.X += OffsetX;
	EvalLocation.Y += OffsetY;
}
std::chrono::high_resolution_clock::time_point PrevFrameTime;

void FractalContext::FindFeature(SRReal CenterX, SRReal CenterY) {
	delete Feature;
	Feature = nullptr;
	Evaluator::FeatureFinder *featureFinder = evaluator->GetFeatureFinder();
	if (featureFinder && reference) {
		SRReal AttractiveRadius = 1.0 / 12.0;
		SRReal MinSearchRadius = AttractiveRadius / 2.0;
		SRReal MaxSearchRadius = AttractiveRadius;
		SRReal Distance, PredCenterSearchRadius;
		HRReal x = CenterX * RenderLocation.HalfH + RenderLocation.X;
		HRReal y = CenterY * RenderLocation.HalfH + RenderLocation.Y;
		for (double SearchRadius = MinSearchRadius; SearchRadius <= MaxSearchRadius; SearchRadius *= 2.0) {
			HRReal CandidateCenterX, CandidateCenterY, CandidateDistance;
			Evaluator::Feature *candidateFeature = featureFinder->FindFeature(x, y, RenderLocation.HalfH * SearchRadius, reference);
			if (candidateFeature) {
				candidateFeature->GetCoordinate(CandidateCenterX, CandidateCenterY);
				CandidateCenterX = (CandidateCenterX - RenderLocation.X) / RenderLocation.HalfH;
				CandidateCenterY = (CandidateCenterY - RenderLocation.Y) / RenderLocation.HalfH;

				HRReal DiffX = CandidateCenterX - CenterX;
				HRReal DiffY = CandidateCenterY - CenterY;

				CandidateDistance = sqrt(SRReal(DiffX * DiffX + DiffY * DiffY));

				if (!Feature || (CandidateDistance * (PredCenterSearchRadius / SearchRadius) * 0.5 < Distance)) {
					delete Feature;
					Feature = candidateFeature;
					Distance = SRReal(CandidateDistance);
					PredCenterSearchRadius = SearchRadius;
				} else {
					delete candidateFeature;
				}
			}
			if (Feature) {
				if (Distance < AttractiveRadius) {
					Attraction = 1.0 - (Distance / AttractiveRadius);

				} else {
					delete Feature;
					Feature = nullptr;
				}
			}
		}
	}
}

void FractalContext::ZoomIn(SRReal CenterX, SRReal CenterY) {
	if (Feature) {
		HRReal FeatureX, FeatureY;
		Feature->GetCoordinate(FeatureX, FeatureY);

		FeatureX = (FeatureX - RenderLocation.X) / RenderLocation.HalfH;
		FeatureY = (FeatureY - RenderLocation.Y) / RenderLocation.HalfH;

		SRReal DiffX = SRReal(FeatureX) - CenterX;
		SRReal DiffY = SRReal(FeatureY) - CenterY;

		Attraction = 1.0 + sqrt(Attraction) * 0.5 * SRReal(CurrentLocation.HalfH / RenderLocation.HalfH);

		CenterX += DiffX * Attraction;
		CenterY += DiffY * Attraction;

		CenterX = std::clamp<double>(CenterX, -double(ImageWidth) / ImageHeight, double(ImageWidth) / ImageHeight);
		CenterY = std::clamp<double>(CenterY, -1.0, 1.0);
	}

	RelLocation NewLocation;
	NewLocation.HalfH = CurrentLocation.HalfH * 0.5_hr;
	HRReal fac = RenderLocation.HalfH - NewLocation.HalfH;
	NewLocation.X = RenderLocation.X + CenterX * fac;
	NewLocation.Y = RenderLocation.Y + CenterY * fac;

	ChangeLocation(NewLocation);
}

void FractalContext::ZoomOut(HRReal CenterX, HRReal CenterY) {
	RelLocation NewLocation = CurrentLocation;
	NewLocation.X -= CenterX * NewLocation.HalfH;
	NewLocation.Y -= CenterY * NewLocation.HalfH;
	NewLocation.HalfH *= 2.0_hr;

	ChangeLocation(NewLocation);
}

void FractalContext::Move(HRReal X, HRReal Y) {
	RelLocation NewLocation = CurrentLocation;
	NewLocation.X += X * NewLocation.HalfH;
	NewLocation.Y += Y * NewLocation.HalfH;

	ChangeLocation(NewLocation);
}

void FractalContext::ChangeLocation(RelLocation NewLocation) {
	if (Zooming && NewLocation.HalfH == CurrentLocation.HalfH) return;
	if (NewLocation.X == CurrentLocation.X && NewLocation.Y == CurrentLocation.Y && NewLocation.HalfH == CurrentLocation.HalfH) {
		return;
	}
	if (LockReference && (abs(NewLocation.X) > NewLocation.HalfH * 0x1p20_hr || abs(NewLocation.Y) > NewLocation.HalfH * HRReal(1ull << 20))) {
		return;
	}
	if (NewLocation.HalfH <= CurrentLocation.HalfH * 0.5_hr || NewLocation.HalfH >= CurrentLocation.HalfH * 2.0_hr) {
		Zooming = true;
		RemainingZoomTime = 0.25;
		PrevFrameTime = std::chrono::high_resolution_clock::now();
	} else {
		RenderLocation = NewLocation;
		EvalLocation = NewLocation;
		ComputePixel = true;
	}
	CurrentLocation = NewLocation;
}

void FractalContext::ZoomToAnimated(Coordinate newCenterCoordinate, size_t precision, HRReal halfH) {
	HRReal Ratio = CurrentLocation.HalfH / halfH;
	SRReal DepthDiff = log2(Ratio);

	UpdateRelativeCoordinate(HRReal(CenterCoordinate.X - newCenterCoordinate.X), HRReal(CenterCoordinate.Y - newCenterCoordinate.Y));
	CenterCoordinate = std::move(newCenterCoordinate);

	CurrentLocation.HalfH = halfH;
	CurrentLocation.X = 0.0_hr;
	CurrentLocation.Y = 0.0_hr;
	CenterChanged = true;

	if (abs(DepthDiff) > 30) {
		if (DepthDiff < 0) {
			SetLocation(newCenterCoordinate, precision, halfH);
			return;
		}
		DepthDiff = 5;
		RenderLocation = FContext.CurrentLocation;
		EvalLocation = FContext.CurrentLocation;
		RenderLocation.HalfH *= 0x1p5;
		EvalLocation.HalfH *= 0x1p5;
		InvalidatePixel();
		Zooming = true;
		RemainingZoomTime = 0.25;
		PrevFrameTime = std::chrono::high_resolution_clock::now();
	} else if (Ratio != 1.0) {
		Zooming = true;
		RemainingZoomTime = std::max(0.35, DepthDiff * 0.05);
		PrevFrameTime = std::chrono::high_resolution_clock::now();
	} else {
		RenderLocation = CurrentLocation;
		EvalLocation = CurrentLocation;
		ComputePixel = true;
	}

	SetDefaultPrecision(precision);
	InvalidateReference();
}

void FractalContext::SetLocation(Coordinate newCenterCoordinate, size_t precision, HRReal halfH) {
	InvalidateAll();
	CenterCoordinate = std::move(newCenterCoordinate);

	CurrentLocation.HalfH = halfH;
	CurrentLocation.X = 0.0_hr;
	CurrentLocation.Y = 0.0_hr;
	RenderLocation = FContext.CurrentLocation;
	EvalLocation = FContext.CurrentLocation;

	SetDefaultPrecision(precision);
}

void FractalContext::ChangeCenter() {
	uint64_t Precision = -std::min(0ll, CurrentLocation.HalfH.Exponent) + 64;
	SetDefaultPrecision(Precision);
	CenterCoordinate.X.set_prec(Precision);
	CenterCoordinate.Y.set_prec(Precision);

	RelLocation PrevLocation = CurrentLocation;
	if (HasPredictedCenter
		&& abs(CurrentLocation.X - PredictedCenter.X) < CurrentLocation.HalfH * 2
		&& abs(CurrentLocation.Y - PredictedCenter.Y) < CurrentLocation.HalfH * 2) {
		CenterCoordinate.X += PredictedCenter.X;
		CenterCoordinate.Y += PredictedCenter.Y;
		CurrentLocation.X -= PredictedCenter.X;
		CurrentLocation.Y -= PredictedCenter.Y;
	} else {
		CenterCoordinate.X += CurrentLocation.X;
		CenterCoordinate.Y += CurrentLocation.Y;
		CurrentLocation.X = 0.0;
		CurrentLocation.Y = 0.0;
	}
	UpdateRelativeCoordinate(CurrentLocation.X - PrevLocation.X, CurrentLocation.Y - PrevLocation.Y);
	CenterChanged = true;
	HasPredictedCenter = false;
}

void FractalContext::AsyncEvaluate(Evaluator *evaluator, PixelManager *pixelManager, Reference *reference) {
	EvaluationParameters parameters;
	parameters.MaxIterations = Global::ItLim;
	//Global::MaxIt = parameters.MaxIterations;
	parameters.CenterCoordinate = CenterCoordinate;
	parameters.Location = EvalLocation;
	EvaluationTask = evaluator->CreateEvaluationTask(parameters, pixelManager, reference);
	EvaluationTaskContext = Computation::AddTask(EvaluationTask);
}

bool FractalContext::IsIdle() {
	if (!pixelManager.Completed()) {
		return false;
	}
	std::lock_guard<std::shared_mutex> lock(TaskMutex);
	return true;
}

void FractalContext::CancelReferenceComputation() {
	if (ReferenceTaskContext) {
		if (!ReferenceTaskContext->Terminated()) {
			ReferenceTaskContext->Cancel();
		}
		ReferenceTask = nullptr;
		ReferenceTaskContext->Wait();
		ReferenceTaskContext->Release();
		ReferenceTaskContext = nullptr;
	}
}

void FractalContext::CancelPixelComputation() {
	if (!IsIdle()) {
		pixelManager.Abort();
	}
	if (EvaluationTaskContext) {
		if (!EvaluationTaskContext->Terminated()) {
			EvaluationTaskContext->Cancel();
		}
		EvaluationTask = nullptr;
		EvaluationTaskContext->Wait();
		EvaluationTaskContext->Release();
		EvaluationTaskContext = nullptr;
	}
	while (!IsIdle());
}

void FractalContext::CancelAll() { // Cancel both, then wait.
	if (!IsIdle()) {
		pixelManager.Abort();
	}
	if (ReferenceTaskContext && !ReferenceTaskContext->Terminated()) {
		ReferenceTaskContext->Cancel();
	}
	if (EvaluationTaskContext && !EvaluationTaskContext->Terminated()) {
		EvaluationTaskContext->Cancel();

		EvaluationTask = nullptr;
		EvaluationTaskContext->Wait();
		EvaluationTaskContext->Release();
		EvaluationTaskContext = nullptr;
	}
	if (ReferenceTaskContext) {
		ReferenceTask = nullptr;
		ReferenceTaskContext->Wait();
		ReferenceTaskContext->Release();
		ReferenceTaskContext = nullptr;
	}
}

FractalContext::FractalContext() : evaluator(new HInfLAEvaluator) {
	RelLocation PrevLocation = CurrentLocation;
	CenterCoordinate.X += CurrentLocation.X;
	CenterCoordinate.Y += CurrentLocation.Y;
	CurrentLocation.X = 0.0;
	CurrentLocation.Y = 0.0;
	if (CurrentLocation.X != PrevLocation.X || CurrentLocation.Y != PrevLocation.Y) {
		UpdateRelativeCoordinate(-CurrentLocation.X, -CurrentLocation.Y);
	}
	NeedReference = true;
}

static size_t OldMaxIt = 0;
void FractalContext::Update() {
	auto CurrentFrameTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> Duration = CurrentFrameTime - PrevFrameTime;
	PrevFrameTime = CurrentFrameTime;

	double DeltaTime = Duration.count();
	if (DeltaTime > 0.1) DeltaTime = 0.1;

	Global::ColoringValueOffset += Global::ColorCyclingSpeed * DeltaTime;
	if (Global::ColoringValueOffset >= 1.0) Global::ColoringValueOffset -= 1.0;
	if (Global::ColoringValueOffset < 1.0) Global::ColoringValueOffset += 1.0;

	if (Zooming) {
		SRReal x = DeltaTime / RemainingZoomTime;
		x = 1.0 - x;
		x = 1.0 - pow(x, 1.6666);
		RemainingZoomTime -= DeltaTime;

		if (RemainingZoomTime > 0.0) {
			HRReal r = RenderLocation.HalfH / CurrentLocation.HalfH;
			HRReal fac = HRReal(pow(double(CurrentLocation.HalfH / RenderLocation.HalfH), x));
			HRReal a = (r * (1.0_hr - fac)) / (r - 1.0_hr);
			RenderLocation.HalfH *= fac;
			RenderLocation.X += (CurrentLocation.X - RenderLocation.X) * a;
			RenderLocation.Y += (CurrentLocation.Y - RenderLocation.Y) * a;

			if (CurrentLocation.HalfH < RenderLocation.HalfH && EvalLocation.HalfH >= RenderLocation.HalfH) {
				EvalLocation.HalfH *= 0.5_hr;
				if (EvalLocation.HalfH <= CurrentLocation.HalfH) {
					EvalLocation = CurrentLocation;
				} else {
					HRReal r2 = RenderLocation.HalfH / CurrentLocation.HalfH;
					HRReal fac2 = EvalLocation.HalfH / RenderLocation.HalfH;
					HRReal a2 = (r2 * (1.0_hr - fac2)) / (r2 - 1.0_hr);
					EvalLocation.X = RenderLocation.X + (CurrentLocation.X - RenderLocation.X) * a2;
					EvalLocation.Y = RenderLocation.Y + (CurrentLocation.Y - RenderLocation.Y) * a2;
				}
				ComputePixel = true;
			} else if (CurrentLocation.HalfH > RenderLocation.HalfH && EvalLocation.HalfH <= RenderLocation.HalfH) {
				EvalLocation.HalfH *= 2.0_hr;
				if (EvalLocation.HalfH >= CurrentLocation.HalfH) {
					EvalLocation = CurrentLocation;
				} else {
					HRReal r2 = RenderLocation.HalfH / CurrentLocation.HalfH;
					HRReal fac2 = EvalLocation.HalfH / RenderLocation.HalfH;
					HRReal a2 = (r2 * (1.0_hr - fac2)) / (r2 - 1.0_hr);
					EvalLocation.X = RenderLocation.X + (CurrentLocation.X - RenderLocation.X) * a2;
					EvalLocation.Y = RenderLocation.Y + (CurrentLocation.Y - RenderLocation.Y) * a2;
				}
				ComputePixel = true;
			}
		} else {
			RenderLocation = CurrentLocation;
			RemainingZoomTime = 0.0;
			Zooming = false;
#ifndef USE_BASIC_PIXEL_MANAGER
			pixelManager.RemovePrevTextures();
#endif
			if (EvalLocation.HalfH != CurrentLocation.HalfH) {
				EvalLocation = CurrentLocation;
				ComputePixel = true;
			}
		}
	}

	if (ComputePixel || ParameterChanged || RecomputeReference) {
		CancelPixelComputation();
		if (LockReference && !RecomputeReference) {
			if (!reference && !ReferenceTask) {
				EvaluationParameters parameters;
				parameters.MaxIterations = Global::ItLim;
				parameters.CenterCoordinate = CenterCoordinate;
				parameters.Location = EvalLocation;

				ReferenceTask = evaluator->CreateReferenceTask(parameters);
				ReferenceTaskContext = Computation::AddTask(ReferenceTask);
			}
			ParameterChanged = false;
		} else {
			HRReal Distance = std::max(abs(EvalLocation.X), abs(EvalLocation.Y));
			if ((!reference || Distance > EvalLocation.HalfH * 64 || EvalLocation.HalfH * HRReal(1.0 / 512.0) < reference->AbsolutePrecision || ParameterChanged || RecomputeReference) && /*!ReferenceWorker.joinable() &&*/ !ReferenceTask) {
				ChangeCenter();
				CenterChanged = false;
				RecomputeReference = false;
				ParameterChanged = false;
				if (reference) evaluator->FreeReference(reference);
				reference = nullptr;
				if (NeedReference) {
					EvaluationParameters parameters;
					parameters.MaxIterations = Global::ItLim;
					parameters.CenterCoordinate = CenterCoordinate;
					parameters.Location = EvalLocation;

					ReferenceTask = evaluator->CreateReferenceTask(parameters);
					ReferenceTaskContext = Computation::AddTask(ReferenceTask);
				}
			}
		}
		LocationChanged = ComputePixel;
		ComputePixel = false;
		EvaluationPending = true;
	}

	if (ReferenceTaskContext && ReferenceTaskContext->Finished()) {
		reference = ReferenceTask->reference;
		ReferenceTask = nullptr;
		ReferenceTaskContext->Release();
		ReferenceTaskContext = nullptr;
	}
	if (EvaluationTaskContext && EvaluationTaskContext->Finished()) {
		EvaluationTask = nullptr;
		EvaluationTaskContext->Release();
		EvaluationTaskContext = nullptr;
	}
	if (EvaluationPending && (!NeedReference || reference)) {
		if (Global::SizeChanged) {
			pixelManager.SetResolution(ImageWidth, ImageHeight);
			Global::SizeChanged = false;
		}
		if (LocationChanged) {
			pixelManager.SetLocation(EvalLocation);
			if (RenderLocation.HalfH == EvalLocation.HalfH) {
				RenderLocation.X = EvalLocation.X;
				RenderLocation.Y = EvalLocation.Y;
			}
		}
#ifndef USE_BASIC_PIXEL_MANAGER
		if (!Zooming) pixelManager.RemovePrevTextures();
#endif
		LocationChanged = false;

		if (Global::ItLim != OldMaxIt) {
#ifndef USE_BASIC_PIXEL_MANAGER
			if (Global::ItLim > OldMaxIt) pixelManager.ChangeMaxit(OldMaxIt);
			else
#endif
			if ((CurrentFractalType == FractalTypeEnum::Mandelbrot && UsingDE && UsingLA)
				|| (CurrentFractalType == FractalTypeEnum::Nova && UsingDE)
				|| Global::PreModulo) pixelManager.Clear();
			//pixelManager.Clear();
			OldMaxIt = Global::ItLim;
		}
		pixelManager.Begin();
		AsyncEvaluate(evaluator, &pixelManager, reference ? reference : reinterpret_cast<Reference *>(&CenterCoordinate));
		EvaluationPending = false;
		using namespace std;
		this_thread::sleep_for(4ms);
	}
}

void FractalContext::GetTextures(TextureDescription *TD, size_t NumDesired, size_t &NumObtained) {
	if (!NumDesired) return;

	pixelManager.GetTextures(TD, NumDesired, NumObtained);
}

void FractalContext::InvalidatePixel() {
	CancelPixelComputation();
	ComputePixel = true;
	pixelManager.Clear();
}

void FractalContext::InvalidateAll() {
	CancelAll();
	RecomputeReference = true;
	ComputePixel = true;
	pixelManager.Clear();
}

void FractalContext::SetResolution(int32_t width, int32_t height) {
	FContext.ImageWidth = width;
	FContext.ImageHeight = height;
	FContext.ComputePixel = true;
	Global::SizeChanged = true;
}

void FractalContext::ChangeFractalType(FractalTypeEnum FractalType, bool UseLinearApproximation) {
	CancelAll();
	if (reference) {
		evaluator->FreeReference(reference);
		reference = nullptr;
	}

	if (FractalType != CurrentFractalType) {
		CurrentLocation.HalfH = 2.0_hr;
		CurrentLocation.X = 0.0_hr;
		CurrentLocation.Y = 0.0_hr;
		RenderLocation = FContext.CurrentLocation;
		EvalLocation = FContext.CurrentLocation;

		CenterCoordinate.X = 0.0;
		CenterCoordinate.Y = 0.0;

		SetDefaultPrecision(64);
		FContext.CenterCoordinate.X.set_prec(64);
		FContext.CenterCoordinate.Y.set_prec(64);
		Global::ItLim = 1024;

		CurrentFractalType = FractalType;
	}


	InvalidateAll();

	delete evaluator;
	switch (FractalType) {
		case FractalTypeEnum::Mandelbrot: {
			UsingLA = UseLinearApproximation;
			if (UseLinearApproximation) {
				evaluator = new HInfLAEvaluator;
			} else {
				evaluator = new PerturbationEvaluator<FractalTypeEnum::Mandelbrot>;
			}
			break;
		}
		case FractalTypeEnum::Tricorn: {
			evaluator = new PerturbationEvaluator<FractalTypeEnum::Tricorn>;
			break;
		}
		case FractalTypeEnum::BurningShip: {
			evaluator = new PerturbationEvaluator<FractalTypeEnum::BurningShip>;
			break;
		}
		case FractalTypeEnum::Nova: {
			evaluator = new NovaEvaluator;
			break;
		}
#ifdef ENABLE_CUSTOM_FORMULA
		case FractalTypeEnum::Custom: {
			evaluator = new JITEvaluator;
			break;
		}
#endif
	}
}
void FractalContext::SetUseDE(bool UseDE) {
	if (UsingDE == UseDE) return;
	CancelPixelComputation();
	pixelManager.Clear();
	UsingDE = UseDE;
	ComputePixel = true;
}

void FractalContext::ResetLocation() {
	CancelAll();
	if (reference) {
		evaluator->FreeReference(reference);
		reference = nullptr;
	}

	pixelManager.Clear();

	CurrentLocation.HalfH = 2.0_hr;
	CurrentLocation.X = 0.0_hr;
	CurrentLocation.Y = 0.0_hr;
	RenderLocation = FContext.CurrentLocation;
	EvalLocation = FContext.CurrentLocation;

	CenterCoordinate.X = 0.0;
	CenterCoordinate.Y = 0.0;

	SetDefaultPrecision(64);
	FContext.CenterCoordinate.X.set_prec(64);
	FContext.CenterCoordinate.Y.set_prec(64);
	Global::ItLim = 1024;

	InvalidateAll();
}

void FractalContext::InvalidateReference() {
	CancelReferenceComputation();
	RecomputeReference = true;
}
