#pragma once

#include <vector>
enum class FractalTypeEnum;
class FractalContext;

extern FractalTypeEnum CurrentFractalType;

extern FractalContext FContext;

extern size_t InitialWindowWidth;
extern size_t InitialWindowHeight;

namespace Global {
	constexpr SRReal MaxPreModMultiplier = 4096.0_sr;
	constexpr SRReal DefaultPreModMultiplier = 64.0_sr;

	extern SRReal ItDiv;
	extern SRReal ItMod;
	extern size_t ItLim;
	extern uint64_t MaxIt;
	extern bool UseBilinearFilter;
	extern bool FlipImaginaryAxis;
	extern bool PreModulo;
	extern bool SizeChanged;
	extern bool Redraw;
	extern bool Initialized;
	extern bool Transform;

	extern std::vector<RGBA> Palette;
	extern bool GammaCorrection;

	extern SRReal ColoringValueOffset;
	extern SRReal ColorCyclingSpeed;

	extern SRReal Rotation;
	extern SRReal StretchAngle;
	extern SRReal StretchRatio;
	extern glm::dmat2 TransformMatrix;
	extern glm::dmat2 InvTransformMatrix;

	extern std::string CustomFormula;

	extern bool HighQuality;

	extern bool SaveReference;
	extern int ReferenceQuality;
}
