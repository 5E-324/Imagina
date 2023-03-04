#ifndef _FRACTAL_CONTEXT_H_
#define _FRACTAL_CONTEXT_H_

#include "Evaluator.h"
#include "FractalContext.h"
#include <thread>
#include <condition_variable>
#include <shared_mutex>

class FractalContext {
public:
	RelLocation CurrentLocation;
	size_t ImageWidth, ImageHeight;

	bool Zooming = false;
	RelLocation RenderLocation;
	double RemainingZoomTime;

	size_t ThreadCount;

	Evaluator *evaluator;
	FractalTypeEnum CurrentFractalType = FractalTypeEnum::Mandelbrot;
	bool UsingDE = true;
	bool UsingLA = true;

	StandardPixelManager pixelManager;
	bool ComputePixel = true, CenterChanged = true, EvaluationPending = true, LocationChanged = true;
	bool RecomputeReference, ParameterChanged;
	bool NeedReference = false;
	Reference *reference = nullptr;

	Evaluator *TaskEvaluator;
	PixelManager *TaskRasterizer;
	Reference *TaskReference;
	std::shared_mutex TaskMutex;
	std::condition_variable_any TaskConditionVariable;
	bool TaskValid = false;
	std::thread *worker = nullptr;
	Coordinate CenterCoordinate;
	std::thread ReferenceWorker;
	Evaluator::ReferenceTask *ReferenceTask = nullptr;
	ExecutionContext *ReferenceTaskContext = nullptr;

	Task *EvaluationTask = nullptr;
	ExecutionContext *EvaluationTaskContext = nullptr;

	RelLocation PredictedCenter;
	bool HasPredictedCenter = false;

	Evaluator::Feature *Feature = nullptr;
	SRReal Attraction;

	RelLocation EvalLocation;

	bool LockReference = false;

	void UpdateRelativeCoordinate(HRReal OffsetX, HRReal OffsetY);
	void FindFeature(SRReal CenterX, SRReal CenterY);
	void ZoomIn(SRReal CenterX, SRReal CenterY);
	void ZoomOut(HRReal CenterX, HRReal CenterY);
	void Move(HRReal X, HRReal Y);
	void ChangeLocation(RelLocation NewLocation);
	void ZoomToAnimated(Coordinate newCenterCoordinate, size_t precision, HRReal halfH);
	void SetLocation(Coordinate newCenterCoordinate, size_t precision, HRReal halfH);
	void ChangeCenter();

	void AsyncEvaluate(Evaluator *e, PixelManager *rasterizer, Reference *Reference);
	bool IsIdle();
	void CancelReferenceComputation();
	void CancelPixelComputation();
	void CancelAll();

	FractalContext();
	void Update();
	void GetTextures(TextureDescription *TD, size_t NumDesired, size_t &NumObtained);

	void PredictCenter(RelLocation NewCenter) {
		PredictedCenter = NewCenter;
		HasPredictedCenter = true;
	}

	void InvalidatePixel();
	void InvalidateAll();

	void SetResolution(int32_t width, int32_t height);

	void ChangeFractalType(FractalTypeEnum FractalType, bool UseLinearApproximation);
	void SetUseDE(bool UseDE);
	void ResetLocation();
	void InvalidateReference();
};

#endif