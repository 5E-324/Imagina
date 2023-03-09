#include "Includes.h"

FractalTypeEnum CurrentFractalType = FractalTypeEnum::Mandelbrot;

FractalContext FContext;

namespace Global {
	SRReal ItDiv = 256.0;
	SRReal ItMod = 256.0 * DefaultPreModMultiplier;
	size_t ItLim = 1024;
	uint64_t MaxIt = 1024;
	bool UseBilinearFilter = false;
	bool FlipVertically = false;
	bool PreModulo = false;
	bool SizeChanged = false;
	bool Redraw = true;
	bool Initialized = false;
	bool Transform = false;

	std::vector<RGBA> Palette = {
		{ 0.0,  0.0,  0.75, 1.0 },
		{ 0.5,  0.0,  1.0,  1.0 },
		{ 1.0,  1.0,  0.0,  1.0 },
		{ 0.75, 0.0,  0.0,  1.0 },
		{ 0.75, 0.0,  1.0,  1.0 },
		{ 0.0,  1.0,  1.0,  1.0 },
	};
	bool GammaCorrection = true;

	SRReal ColoringValueOffset = 0.0;
	SRReal ColorCyclingSpeed = 0.0;

	SRReal Rotation = 0.0;
	SRReal StretchAngle = 0.0;
	SRReal StretchRatio = 0.0;
	glm::dmat2 TransformMatrix;
	glm::dmat2 InvTransformMatrix;

	std::string CustomFormula = "z^2 + c";

	bool HighQuality = true;
}