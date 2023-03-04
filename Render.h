#ifndef _RENDER_H_
#define _RENDER_H_

#include "FractalContext.h"

void InitRenderResources();
void ResizeViewport(int32_t Width, int32_t Height);
void BeginRender();
void EndRender();
void RenderFractal(FractalContext &Context);
void EnablePaletteMipmap(bool Enable);

void SetPalette(const std::vector<RGBA> &Palette);

void RenderLineToFeature(double mouseX, double mouseY, double featureX, double fetureY);

#endif