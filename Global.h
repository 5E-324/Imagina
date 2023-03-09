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
	extern bool FlipVertically;
	extern bool PreModulo;
	extern bool SizeChanged;
	extern bool Redraw;
	extern bool Initialized;

	extern std::vector<RGBA> Palette;
	extern bool GammaCorrection;

	extern SRReal ColoringValueOffset;
	extern SRReal ColorCyclingSpeed;

	extern std::string CustomFormula;

	extern bool HighQuality;
}
