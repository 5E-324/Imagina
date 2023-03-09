#include "Includes.h"
#include <GL/glu.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <fenv.h>
#include <gmp-impl.h>

void OpenFile(wchar_t *FileName, size_t ExtensionOffset = 0);

int WINAPI wWinMain(HINSTANCE /*hInstance*/, HINSTANCE /*hPrevInstance*/, LPWSTR /*lpCmdLine*/, int /*nCmdShow*/) {
	CreateMainWindow();
	InitRenderResources();
	Computation::Init();

	FContext.CurrentLocation = { 0.0, 0.0, 2.0 };
	FContext.RenderLocation = { 0.0, 0.0, 2.0 };
	FContext.EvalLocation = { 0.0, 0.0, 2.0 };

	FContext.ImageWidth = InitialWindowWidth;
	FContext.ImageHeight = InitialWindowHeight;

	Global::Initialized = true;
	bool UseGetMessage = false;

	int argc;
	LPWSTR *argv = CommandLineToArgvW(GetCommandLineW(), &argc);
	if (argc == 2) {
		OpenFile(argv[1]);
	}

	LocalFree(argv);

	while (true) {
		MSG Message;
		if (UseGetMessage) {
			while (true) {
				BOOL ret = GetMessage(&Message, nullptr, 0, 0);
				if (!ret) {
					FContext.pixelManager.Abort();
					while (!FContext.pixelManager.Completed());
					if (FContext.ReferenceWorker.joinable()) FContext.ReferenceWorker.detach();
					return 0;
				}
				TranslateMessage(&Message);
				DispatchMessage(&Message);
				if (FContext.ReferenceTaskContext || FContext.ParameterChanged || FContext.RecomputeReference || FContext.ComputePixel || FContext.Zooming || Global::Redraw) {
					UseGetMessage = false;
					break;
				}
			}
		} else {
			while (PeekMessage(&Message, nullptr, 0, 0, PM_REMOVE)) {
				TranslateMessage(&Message);
				DispatchMessage(&Message);
				if (Message.message == WM_QUIT) {
					FContext.pixelManager.Abort();
					while (!FContext.pixelManager.Completed());
					if (FContext.ReferenceWorker.joinable()) FContext.ReferenceWorker.detach();
					return 0;
				}
			}
		}

		FContext.Update();

		BeginRender();

		if (FContext.pixelManager.Completed() && !(FContext.ReferenceTaskContext || FContext.ParameterChanged || FContext.RecomputeReference || FContext.ComputePixel || FContext.Zooming || Global::Redraw)) UseGetMessage = true;
		else UseGetMessage = false;
		RenderFractal(FContext);

		POINT MousePoint{ MouseX, MouseY };
		SRReal X = (SRReal((int)MousePoint.x) - WindowWidth * 0.5) / WindowHeight * 2.0;
		SRReal Y = -(SRReal((int)MousePoint.y) / WindowHeight * 2.0 - 1.0);

		SRReal FractalX = Global::InvTransformMatrix[0][0] * X + Global::InvTransformMatrix[1][0] * Y;
		SRReal FractalY = Global::InvTransformMatrix[0][1] * X + Global::InvTransformMatrix[1][1] * Y;

		FContext.FindFeature(FractalX, FractalY);
		if (FContext.Feature) {
			HRReal FeatureX, FeatureY;
			FContext.Feature->GetCoordinate(FeatureX, FeatureY);

			FeatureX = (FeatureX - FContext.RenderLocation.X) / FContext.RenderLocation.HalfH; // To screen coordinate
			FeatureY = (FeatureY - FContext.RenderLocation.Y) / FContext.RenderLocation.HalfH;

			HRReal FeatureScreenX = Global::TransformMatrix[0][0] * FeatureX + Global::TransformMatrix[1][0] * FeatureY;
			HRReal FeatureScreenY = Global::TransformMatrix[0][1] * FeatureX + Global::TransformMatrix[1][1] * FeatureY;

			RenderLineToFeature(X, Y, (double)FeatureScreenX, (double)FeatureScreenY);
		}

		Global::Redraw = false;
		if (Global::ColorCyclingSpeed != 0.0) Global::Redraw = true;

		if (FContext.pixelManager.Completed() && !(FContext.ReferenceTaskContext || FContext.ParameterChanged || FContext.RecomputeReference || FContext.ComputePixel || FContext.Zooming || Global::Redraw)) {
		} else UseGetMessage = false;

		EndRender();
	}
}