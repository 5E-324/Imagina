#ifndef _MAIN_WINDOW_H_
#define _MAIN_WINDOW_H_

#include "FractalContext.h"

extern HGLRC GLContext;
extern HWND HWnd;
extern HDC MainDC;
extern int32_t MouseX, MouseY;

void CreateMainWindow();

enum class MenuID : UINT_PTR {
	Open,
	Save,
	SaveImage,
	SaveRawPixelData,
	ReferenceSaving,

	DistanceEstimation,
	HighQuality,
	Location,
	IterationLimit,
	Transformation,

	IncreaseColorDensity,
	DecreaseColorDensity,
	IncreaseIterationLimit,
	DecreaseIterationLimit,

	ResetLocation,

	LockReference,
	Tasks,

	ImageSize,
	PaletteMipmaps,
	BilinearFilter,
	PreModulo,

	FractalTypeBegin,
	FractalTypeMandelbrot = FractalTypeBegin,
	FractalTypeTricorn,
	FractalTypeBurningShip,
	FractalTypeNova,
	FractalTypeCustom,
	FractalTypeEnd = FractalTypeCustom,

	AlgorithmBegin,
	AlgorithmPerturbation = AlgorithmBegin,
	AlgorithmLinearApproximation,
	AlgorithmEnd = AlgorithmLinearApproximation,

	RecomputeReference,
	RecomputePixel,
	RecomputeAll,
};


extern size_t WindowWidth;
extern size_t WindowHeight;
extern FractalTypeEnum CurrentFractalType;

void SetFractalType(FractalTypeEnum Type);

LRESULT CALLBACK WindowProcess(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam);

#endif