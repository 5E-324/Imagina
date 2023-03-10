#include "Includes.h"
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <png.h>
#include "Resource.h"
#include "windowsx.h"
#include "gmp-impl.h"

#include "File.h"

size_t InitialWindowWidth = 2048;
size_t InitialWindowHeight = 1024;

HGLRC GLContext;
HWND HWnd;
HDC MainDC;

size_t WindowWidth = InitialWindowWidth;
size_t WindowHeight = InitialWindowHeight;

bool FullScreen = false;

HMENU MainMenu;
HMENU File;
HMENU Fractal;
HMENU Computation;
HMENU Image;
HMENU FractalType;
HMENU AlgorithmMenu;
HMENU Recompute;

HWND TasksDialog = nullptr;

bool UseLinearApproximation = true;
bool UseDE = true;
bool UsePalleteMipmaps = false;

int32_t MouseX = 0, MouseY = 0;
bool Panning = false;
static bool SizeChanged = false;
static bool SizeMove = false;

namespace Global {
	bool moved = false;
	bool ImageSizeFollowWindowSize = true;
}

bool InitPixelFormat(HDC hDC) {
	PIXELFORMATDESCRIPTOR pfd = {}, *ppfd;
	int pixelformat;

	ppfd = &pfd;

	ppfd->nSize = sizeof(PIXELFORMATDESCRIPTOR);
	ppfd->nVersion = 1;
	ppfd->dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
	ppfd->dwLayerMask = PFD_MAIN_PLANE;
	ppfd->iPixelType = PFD_TYPE_RGBA;
	ppfd->cColorBits = 8;
	ppfd->cDepthBits = 16;
	ppfd->cAccumBits = 0;
	ppfd->cStencilBits = 0;

	pixelformat = ChoosePixelFormat(hDC, ppfd);

	if (pixelformat == 0) {
		MessageBox(nullptr, L"ChoosePixelFormat failed", L"Error", MB_OK);
		return false;
	}

	if (SetPixelFormat(hDC, pixelformat, ppfd) == FALSE) {
		MessageBox(nullptr, L"SetPixelFormat failed", L"Error", MB_OK);
		return false;
	}

	return true;
}

void CreateMainWindow() {
	MainMenu = CreateMenu();
	File = CreateMenu();
	Fractal = CreateMenu();
	Computation = CreateMenu();
	Image = CreateMenu();
	FractalType = CreateMenu();
	AlgorithmMenu = CreateMenu();
	Recompute = CreateMenu();

	AppendMenuW(MainMenu, MF_POPUP, (UINT_PTR)File, L"File");
	AppendMenuW(MainMenu, MF_POPUP, (UINT_PTR)Fractal, L"Fractal");
	AppendMenuW(MainMenu, MF_POPUP, (UINT_PTR)Computation, L"Computation");
	AppendMenuW(MainMenu, MF_POPUP, (UINT_PTR)Image, L"Image");

	AppendMenuW(File, MF_STRING, (UINT_PTR)MenuID::Open, L"Open\tCtrl+O");
	AppendMenuW(File, MF_STRING, (UINT_PTR)MenuID::Save, L"Save\tCtrl+S");
	AppendMenuW(File, MF_STRING, (UINT_PTR)MenuID::SaveImage, L"Save image\tCtrl+Shift+S");
	AppendMenuW(File, MF_STRING, (UINT_PTR)MenuID::SaveRawPixelData, L"Save raw pixel data");

	AppendMenuW(Fractal, MF_POPUP, (UINT_PTR)FractalType, L"Formula");
	AppendMenuW(Fractal, MF_STRING | MF_CHECKED, (UINT_PTR)MenuID::DistanceEstimation, L"Distance estimation");
	AppendMenuW(Fractal, MF_STRING | MF_CHECKED, (UINT_PTR)MenuID::HighQuality, L"High quality (Experimental)");

	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::Location, L"Location");
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::IterationLimit, L"Iteration limit");
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::Transformation, L"Transformation");

	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::IncreaseColorDensity, L"Increase color density\tA");
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::DecreaseColorDensity, L"Decrease color density\tS");
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::IncreaseIterationLimit, L"Increase iteration limit\tD");
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::DecreaseIterationLimit, L"Decrease iteration limit\tF");
	AppendMenuW(Fractal, MF_SEPARATOR, 0, nullptr);
	AppendMenuW(Fractal, MF_STRING, (UINT_PTR)MenuID::ResetLocation, L"Reset location");

	AppendMenuW(Computation, MF_POPUP, (UINT_PTR)AlgorithmMenu, L"Algorithm");
	AppendMenuW(Computation, MF_POPUP, (UINT_PTR)Recompute, L"Recompute");
	AppendMenuW(Computation, MF_STRING, (UINT_PTR)MenuID::LockReference, L"Lock reference");
	AppendMenuW(Computation, MF_STRING, (UINT_PTR)MenuID::Tasks, L"Tasks");

	AppendMenuW(Image, MF_STRING, (UINT_PTR)MenuID::ImageSize, L"Image size");
	AppendMenuW(Image, MF_STRING, (UINT_PTR)MenuID::BilinearFilter, L"Bilinear filter");
	AppendMenuW(Image, MF_STRING, (UINT_PTR)MenuID::FlipVertically, L"Flip imaginary axis");
	AppendMenuW(Image, MF_STRING, (UINT_PTR)MenuID::PreModulo, L"Pre-modulo (workaround for precision loss)");

	AppendMenuW(FractalType, MF_STRING, (UINT_PTR)MenuID::FractalTypeMandelbrot, L"Mandelbrot");
	AppendMenuW(FractalType, MF_STRING, (UINT_PTR)MenuID::FractalTypeTricorn, L"Tricorn");
	AppendMenuW(FractalType, MF_STRING, (UINT_PTR)MenuID::FractalTypeBurningShip, L"Burning ship");
	AppendMenuW(FractalType, MF_STRING, (UINT_PTR)MenuID::FractalTypeNova, L"Nova");
#ifdef ENABLE_CUSTOM_FORMULA
	AppendMenuW(FractalType, MF_STRING, (UINT_PTR)MenuID::FractalTypeCustom, L"Custom (Experimental)");
#endif
	CheckMenuRadioItem(FractalType, (UINT_PTR)MenuID::FractalTypeBegin, (UINT_PTR)MenuID::FractalTypeEnd, (UINT_PTR)MenuID::FractalTypeMandelbrot, MF_BYCOMMAND);

	AppendMenuW(AlgorithmMenu, MF_STRING, (UINT_PTR)MenuID::AlgorithmPerturbation, L"Perturbation");
	AppendMenuW(AlgorithmMenu, MF_STRING, (UINT_PTR)MenuID::AlgorithmLinearApproximation, L"Linear approximation");
	CheckMenuRadioItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmBegin, (UINT_PTR)MenuID::AlgorithmEnd, (UINT_PTR)MenuID::AlgorithmLinearApproximation, MF_BYCOMMAND);

	AppendMenuW(Recompute, MF_STRING, (UINT_PTR)MenuID::RecomputeReference, L"Reference");
	AppendMenuW(Recompute, MF_STRING, (UINT_PTR)MenuID::RecomputePixel, L"Pixel");
	AppendMenuW(Recompute, MF_STRING, (UINT_PTR)MenuID::RecomputeAll, L"All");





	WNDCLASSEX WindowClass;
	WindowClass.cbSize = sizeof(WNDCLASSEX);

	WindowClass.style = CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS | CS_OWNDC;
	WindowClass.lpfnWndProc = WindowProcess;
	WindowClass.cbClsExtra = 0;
	WindowClass.cbWndExtra = 0;
	WindowClass.hInstance = nullptr;
	WindowClass.hIcon = nullptr;
	WindowClass.hCursor = LoadCursor(nullptr, IDC_ARROW);
	WindowClass.hbrBackground = nullptr;
	WindowClass.lpszMenuName = nullptr;
	WindowClass.lpszClassName = L"ImaginaMainWindow";
	WindowClass.hIconSm = nullptr;

	RegisterClassEx(&WindowClass);

	MONITORINFO MonitorInfo = {};
	MonitorInfo.cbSize = sizeof(MONITORINFO);

	RECT rect;
	rect.left = 0;
	rect.top = 0;
	rect.right = InitialWindowWidth;
	rect.bottom = InitialWindowHeight;

	AdjustWindowRectEx(&rect, WS_OVERLAPPEDWINDOW, TRUE, NULL);
	size_t BorderWidth = rect.right - rect.left - InitialWindowWidth;
	size_t BorderHeight = rect.bottom - rect.top - InitialWindowHeight;
	InitialWindowWidth += BorderWidth;
	InitialWindowHeight += BorderHeight;

	GetMonitorInfo(MonitorFromWindow(GetDesktopWindow(), MONITOR_DEFAULTTOPRIMARY), &MonitorInfo);
	InitialWindowWidth = std::min(long(InitialWindowWidth), MonitorInfo.rcWork.right - MonitorInfo.rcWork.left - 50);
	InitialWindowHeight = std::min(long(InitialWindowHeight), MonitorInfo.rcWork.bottom - MonitorInfo.rcWork.top - 50);

	if ((InitialWindowHeight - BorderHeight) * 2 < InitialWindowWidth - BorderWidth) {
		InitialWindowWidth = (InitialWindowHeight - BorderHeight) * 2 + BorderWidth;
	}

	HWnd = CreateWindowEx(
		NULL,
		L"ImaginaMainWindow",
		L"Imagina",
		WS_OVERLAPPEDWINDOW,
		(MonitorInfo.rcWork.left + MonitorInfo.rcWork.right - InitialWindowWidth) / 2,
		(MonitorInfo.rcWork.top + MonitorInfo.rcWork.bottom - InitialWindowHeight) / 2,
		InitialWindowWidth,
		InitialWindowHeight,
		nullptr,
		MainMenu,
		nullptr,
		nullptr);

	InitialWindowWidth -= BorderWidth;
	InitialWindowHeight -= BorderHeight;

	WindowWidth = InitialWindowWidth;
	WindowHeight = InitialWindowHeight;
	Global::SizeChanged = true;

	MainDC = GetDC(HWnd);



	InitPixelFormat(MainDC);
	GLContext = wglCreateContext(MainDC);
	wglMakeCurrent(MainDC, GLContext);

	ShowWindow(HWnd, SW_NORMAL);
}

struct FractalTypeInfo {
	MenuID MenuID;
	bool SupportPerturbation : 1;
	bool SupportLinearApproximation : 1;
	bool PTSupportDistanceEstimation : 1;
	bool LASupportDistanceEstimation : 1;
	bool SupportHighQuality : 1;
};

const std::unordered_map<FractalTypeEnum, FractalTypeInfo> FractalInfoMap{
	{ FractalTypeEnum::Mandelbrot,	{ MenuID::FractalTypeMandelbrot,	true,  true,  false, true , true  } },
	{ FractalTypeEnum::Tricorn,		{ MenuID::FractalTypeTricorn,		true,  false, false, false, false } },
	{ FractalTypeEnum::BurningShip, { MenuID::FractalTypeBurningShip,	true,  false, false, false, false } },
	{ FractalTypeEnum::Nova,		{ MenuID::FractalTypeNova,			true,  false, true,  false, false } },
	{ FractalTypeEnum::Custom,		{ MenuID::FractalTypeCustom,		false, false, false, false, false } },
};

void SetFractalType(FractalTypeEnum Type) {
	FractalTypeInfo Info = FractalInfoMap.find(Type)->second;

	CheckMenuRadioItem(FractalType, (UINT_PTR)MenuID::FractalTypeBegin, (UINT_PTR)MenuID::FractalTypeEnd, (UINT_PTR)Info.MenuID, MF_BYCOMMAND);
	EnableMenuItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmPerturbation, MF_BYCOMMAND | (Info.SupportPerturbation ? MF_ENABLED : MF_GRAYED));
	EnableMenuItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmLinearApproximation, MF_BYCOMMAND | (Info.SupportLinearApproximation ? MF_ENABLED : MF_GRAYED));
	if (UseLinearApproximation && Info.SupportLinearApproximation) {
		CheckMenuRadioItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmBegin, (UINT_PTR)MenuID::AlgorithmEnd, (UINT_PTR)MenuID::AlgorithmLinearApproximation, MF_BYCOMMAND);
		EnableMenuItem(Fractal, (UINT_PTR)MenuID::DistanceEstimation, MF_BYCOMMAND | (Info.LASupportDistanceEstimation ? MF_ENABLED : MF_GRAYED));
		EnableMenuItem(Fractal, (UINT_PTR)MenuID::HighQuality, MF_BYCOMMAND | (Info.SupportHighQuality ? MF_ENABLED : MF_GRAYED));
	} else if (Info.SupportPerturbation) {
		CheckMenuRadioItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmBegin, (UINT_PTR)MenuID::AlgorithmEnd, (UINT_PTR)MenuID::AlgorithmPerturbation, MF_BYCOMMAND);
		EnableMenuItem(Fractal, (UINT_PTR)MenuID::DistanceEstimation, MF_BYCOMMAND | (Info.PTSupportDistanceEstimation ? MF_ENABLED : MF_GRAYED));
		EnableMenuItem(Fractal, (UINT_PTR)MenuID::HighQuality, MF_BYCOMMAND | MF_GRAYED);
	}
	FContext.ChangeFractalType(Type, UseLinearApproximation && Info.SupportLinearApproximation);
	CurrentFractalType = Type;
}

void OpenLocation() {
	OPENFILENAMEW OpenFileName;
	wchar_t FileName[MAX_PATH] = L"";

	OpenFileName.lStructSize = sizeof(OPENFILENAMEW);
	OpenFileName.hwndOwner = HWnd;
	OpenFileName.hInstance = nullptr;
	OpenFileName.lpstrFilter = L""
		"All supported formats(*.im;*.imt;*.kfr;*.kfp)\0*.im;*.imt;*.kfr;*.kfp\0"
		"Imagina(*.im)\0*.im\0"
		"Imagina text(*.imt)\0*.imt\0"
		"Kalle's fraktaler(*.kfr)\0*.kfr\0"
		"Kalle's fraktaler palette(*.kfp)\0*.kfp\0"
		"All Files(*.*)\0*.*\0";
	OpenFileName.lpstrCustomFilter = nullptr;
	OpenFileName.nMaxCustFilter = 0;
	OpenFileName.nFilterIndex = 1;
	OpenFileName.lpstrFile = FileName;
	OpenFileName.nMaxFile = sizeof(FileName) / sizeof(wchar_t);
	OpenFileName.lpstrFileTitle = nullptr;
	OpenFileName.nMaxFileTitle = 0;
	OpenFileName.lpstrInitialDir = nullptr;
	OpenFileName.lpstrTitle = nullptr;
	OpenFileName.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST;
	OpenFileName.lpstrDefExt = L"im";

	if (!GetOpenFileNameW(&OpenFileName)) return;

	OpenFile(FileName, OpenFileName.nFileExtension);
}
void SaveLocation() {
	using enum FileType;

	OPENFILENAMEW OpenFileName;
	wchar_t FileName[MAX_PATH] = L"";

	OpenFileName.lStructSize = sizeof(OPENFILENAMEW);
	OpenFileName.hwndOwner = HWnd;
	OpenFileName.hInstance = nullptr;
	OpenFileName.lpstrFilter = L""
		"Imagina(*.im)\0*.im\0"
		"Imagina text(*.imt)\0*.imt\0"
		"Kalle's fraktaler(*.kfr)\0*.kfr\0";
	OpenFileName.lpstrCustomFilter = nullptr;
	OpenFileName.nMaxCustFilter = 0;
	OpenFileName.nFilterIndex = 1;
	OpenFileName.lpstrFile = FileName;
	OpenFileName.nMaxFile = sizeof(FileName) / sizeof(wchar_t);
	OpenFileName.lpstrFileTitle = nullptr;
	OpenFileName.nMaxFileTitle = 0;
	OpenFileName.lpstrInitialDir = nullptr;
	OpenFileName.lpstrTitle = nullptr;
	OpenFileName.Flags = OFN_EXPLORER | OFN_OVERWRITEPROMPT;
	OpenFileName.lpstrDefExt = L"im";

	if (!GetSaveFileNameW(&OpenFileName)) return;

	FileType Type;
	switch (OpenFileName.nFilterIndex) {
		case 1: Type = Imagina; break;
		case 2: Type = ImaginaText; break;
		case 3: Type = Kfr; break;
	}
	SaveFile(FileName, Type);
}
void SaveImage(bool raw = false) {
	if (!FContext.pixelManager.Completed()) {
		MessageBoxA(nullptr, "Please wait for computations to finish.", "Save image", MB_OK | MB_ICONASTERISK | MB_TASKMODAL);
		return;
	}

	OPENFILENAMEW OpenFileName;
	wchar_t FileName[MAX_PATH] = L"";

	OpenFileName.lStructSize = sizeof(OPENFILENAMEW);
	OpenFileName.hwndOwner = HWnd;
	OpenFileName.hInstance = nullptr;
	OpenFileName.lpstrFilter = raw ? L"All Files(*.*)\0*.*\0" : L"PNG(*.png)\0*.png\0All Files(*.*)\0*.*\0";
	OpenFileName.lpstrCustomFilter = nullptr;
	OpenFileName.nMaxCustFilter = 0;
	OpenFileName.nFilterIndex = 1;
	OpenFileName.lpstrFile = FileName;
	OpenFileName.nMaxFile = sizeof(FileName) / sizeof(wchar_t);
	OpenFileName.lpstrFileTitle = nullptr;
	OpenFileName.nMaxFileTitle = 0;
	OpenFileName.lpstrInitialDir = nullptr;
	OpenFileName.lpstrTitle = nullptr;
	OpenFileName.Flags = OFN_EXPLORER | OFN_OVERWRITEPROMPT;
	OpenFileName.lpstrDefExt = raw ? L"" : L"png";

	if (!GetSaveFileNameW(&OpenFileName)) return;

	if (raw) {
		SaveRawPixelData(FileName);
	} else {
		SaveImage(FileName);
	}
}

INT_PTR SetFormulaProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM /*lParam*/) {
	switch (message) {
		case WM_INITDIALOG: {
			char Buffer[32];
			sprintf_s(Buffer, "%zu", Global::ItLim);

			SetDlgItemTextA(hWndDlg, IDC_FORMULA, Global::CustomFormula.c_str());

			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDOK: {
					int Length = GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_FORMULA));
					Length++;
					char *Buffer = new char[Length];

					if (GetDlgItemTextA(hWndDlg, IDC_FORMULA, Buffer, Length)) {
						Global::CustomFormula = Buffer;
					}
					delete[]Buffer;
					SetFractalType(FractalTypeEnum::Custom);
				}
				[[fallthrough]];
				case IDCANCEL: {
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
			}
			break;
		}
	}

	return 0;
}

enum class ZoomType {
	Height = 0,
	HalfHeight = 1,
	Magnification = 2,
	Depth = 3,
};

std::string ZoomToString(ZoomType type, HRReal halfH, HRReal OriginalHalfH) {
	char Buffer[32];
	switch (type) {
		case ZoomType::Height: {
			gmp_sprintf(Buffer, "%.12Fg", (halfH * 2.0_hr).to_mpf_class(53).get_mpf_t());
			break;
		}
		case ZoomType::HalfHeight: {
			gmp_sprintf(Buffer, "%.12Fg", halfH.to_mpf_class(53).get_mpf_t());
			break;
		}
		case ZoomType::Magnification: {
			gmp_sprintf(Buffer, "%.12Fg", (OriginalHalfH / halfH).to_mpf_class(53).get_mpf_t());
			break;
		}
		case ZoomType::Depth: {
			sprintf_s(Buffer, "%.15g", log2(OriginalHalfH / halfH));
			break;
		}
	}
	return Buffer;
}

HRReal ParseZoom(ZoomType type, std::string_view string, HRReal OriginalHalfH) {
	HRReal Zoom = mpf_class(string.data(), 52);

	switch (type) {
		case ZoomType::Height: {
			Zoom *= 0.5;
			break;
		}
		case ZoomType::HalfHeight: {
			break;
		}
		case ZoomType::Magnification: {
			Zoom = OriginalHalfH / Zoom;
			break;
		}
		case ZoomType::Depth: {
			Zoom = OriginalHalfH / pow2(SRReal(Zoom));
			break;
		}
	}
	return Zoom;
}

INT_PTR SetLocationProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam) {
	static int ItemIndex = 2;
	switch (message) {
		case WM_INITDIALOG: {
			HWND ComboBox = GetDlgItem(hWndDlg, IDC_COMBO1);

			SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Height");
			SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Half height");
			SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Magnification");
			SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Depth");

			SendMessage(ComboBox, CB_SETCURSEL, (WPARAM)ItemIndex, (LPARAM)0);

			size_t NumberOfDigits = (size_t)std::max(0.0, -log10(FContext.CurrentLocation.HalfH)) + 10;
			size_t Length = NumberOfDigits + 32;
			char *Buffer = new char[Length];

			gmp_sprintf(Buffer, "%.*Ff", NumberOfDigits, mpf_class(FContext.CenterCoordinate.X + HPReal(FContext.CurrentLocation.X)).get_mpf_t());
			SetDlgItemTextA(hWndDlg, IDC_REAL, Buffer);

			gmp_sprintf(Buffer, "%.*Ff", NumberOfDigits, mpf_class(FContext.CenterCoordinate.Y + HPReal(FContext.CurrentLocation.Y)).get_mpf_t());
			SetDlgItemTextA(hWndDlg, IDC_IMAGINARY, Buffer);

			gmp_sprintf(Buffer, "%.12Fg", (2.0_hr / FContext.CurrentLocation.HalfH).to_mpf_class(53).get_mpf_t());
			SetDlgItemTextA(hWndDlg, IDC_ZOOM, Buffer);

			sprintf_s(Buffer, 32, "%zu", Global::ItLim);
			SetDlgItemTextA(hWndDlg, IDC_ITERATIONS, Buffer);

			delete[]Buffer;
			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDC_COMBO1: {
					if (HIWORD(wParam) == CBN_SELCHANGE) {
						ItemIndex = SendMessage((HWND)lParam, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
						SetDlgItemTextA(hWndDlg, IDC_ZOOM, ZoomToString(ZoomType(ItemIndex), FContext.CurrentLocation.HalfH, 2.0_hr).data());
					}
					break;
				}
				case IDOK: {
					int Length = GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_REAL));
					Length = std::max(Length, GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_IMAGINARY)));
					Length = std::max(Length, GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_ZOOM)));
					Length = std::max(Length, GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_ITERATIONS)));
					Length++;
					char *Buffer = new char[Length];

					if (GetDlgItemTextA(hWndDlg, IDC_ZOOM, Buffer, Length)) {
						mpf_class PrevX(std::move(FContext.CenterCoordinate.X)), PrevY(std::move(FContext.CenterCoordinate.Y));
						try {
							HRReal HalfH = ParseZoom(ZoomType(ItemIndex), Buffer, 2.0_hr);
							if (HalfH > 16.0_hr) HalfH = 16.0_hr;
							uint64_t Precision = -std::min(0ll, HalfH.Exponent) + 64;

							GetDlgItemTextA(hWndDlg, IDC_REAL, Buffer, Length);
							FContext.CenterCoordinate.X = mpf_class(Buffer, Precision);
							GetDlgItemTextA(hWndDlg, IDC_IMAGINARY, Buffer, Length);
							FContext.CenterCoordinate.Y = mpf_class(Buffer, Precision);

							if (GetDlgItemTextA(hWndDlg, IDC_ITERATIONS, Buffer, Length)) {
								size_t NewItLim;
								char Dummy;
								if (sscanf_s(Buffer, "%zu%c", &NewItLim, &Dummy) == 1) {
									Global::ItLim = std::max(UINT64_C(2), NewItLim);
								} else {
									throw std::invalid_argument("Iterations");
								}
							}

							FContext.CurrentLocation.HalfH = HalfH;
							FContext.CurrentLocation.X = 0.0_hr;
							FContext.CurrentLocation.Y = 0.0_hr;
							FContext.RenderLocation = FContext.CurrentLocation;
							FContext.EvalLocation = FContext.CurrentLocation;
							FContext.InvalidateAll();

							SetDefaultPrecision(Precision);

						} catch (std::invalid_argument) {
							FContext.CenterCoordinate.X = std::move(PrevX);
							FContext.CenterCoordinate.Y = std::move(PrevY);
							delete[]Buffer;
							MessageBoxA(hWndDlg, "Value invalid.", nullptr, MB_OK | MB_ICONERROR);
							break;
						}
					}
					delete[]Buffer;
				}
				[[fallthrough]];
				case IDCANCEL: {
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
			}
			break;
		}
	}

	return 0;
}

INT_PTR SetIterationLimitProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM /*lParam*/) {
	switch (message) {
		case WM_INITDIALOG: {
			char Buffer[32];
			sprintf_s(Buffer, "%zu", Global::ItLim);

			SetDlgItemTextA(hWndDlg, IDC_EDIT1, Buffer);

			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDOK: {
					int Length = GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_EDIT1));
					Length++;
					char *Buffer = new char[Length];

					if (GetDlgItemTextA(hWndDlg, IDC_EDIT1, Buffer, Length)) {
						size_t NewItLim;
						char Dummy;
						if (sscanf_s(Buffer, "%zu%c", &NewItLim, &Dummy) == 1) {
							Global::ItLim = std::max(UINT64_C(2), NewItLim);
							FContext.ParameterChanged = true;
						}
					}
					delete[]Buffer;
				}
				[[fallthrough]];
				case IDCANCEL: {
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
			}
			break;
		}
	}

	return 0;
}
INT_PTR TransformProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM /*lParam*/) {
	switch (message) {
		case WM_INITDIALOG: {
			char buffer[32];
			
			SendMessage(GetDlgItem(hWndDlg, IDC_FLIP_IMAGINARY), BM_SETCHECK, (WPARAM)(Global::FlipVertically ? BST_CHECKED : BST_UNCHECKED), (LPARAM)0);

			sprintf_s(buffer, "%g", Global::Rotation);
			SetDlgItemTextA(hWndDlg, IDC_ROTATION, buffer);

			sprintf_s(buffer, "%g", Global::StretchAngle);
			SetDlgItemTextA(hWndDlg, IDC_STRETCH_ANGLE, buffer);

			sprintf_s(buffer, "%g", Global::StretchRatio);
			SetDlgItemTextA(hWndDlg, IDC_STRETCH_RATIO, buffer);

			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDOK: {
					char buffer[32];
					buffer[32] = 0;
					SRReal newRotation, newStretchAngle, newStretchRatio;
					try {
						GetDlgItemTextA(hWndDlg, IDC_ROTATION, buffer, 31);
						newRotation = std::stof(buffer);

						GetDlgItemTextA(hWndDlg, IDC_STRETCH_ANGLE, buffer, 31);
						newStretchAngle = std::stof(buffer);

						GetDlgItemTextA(hWndDlg, IDC_STRETCH_RATIO, buffer, 31);
						newStretchRatio = std::stof(buffer);

						if (newStretchRatio == 0.0) throw std::invalid_argument("StretchRatio");
					} catch (std::invalid_argument) {
						MessageBoxA(hWndDlg, "Value invalid.", nullptr, MB_OK | MB_ICONERROR);
						break;
					} catch (std::out_of_range) {
						MessageBoxA(hWndDlg, "Value out of range.", nullptr, MB_OK | MB_ICONERROR);
						break;
					}

					Global::FlipVertically = SendMessage(GetDlgItem(hWndDlg, IDC_FLIP_IMAGINARY), BM_GETCHECK, (WPARAM)0, (LPARAM)0) == BST_CHECKED;

					Global::Rotation = newRotation;
					Global::StretchAngle = newStretchAngle;
					Global::StretchRatio = newStretchRatio;

					if (!Global::FlipVertically && newRotation == 0.0 && newStretchRatio == 1.0) {
						Global::TransformMatrix = glm::dmat2(1.0);
						Global::InvTransformMatrix = glm::dmat2(1.0);
						Global::Transform = false;
					} else {
						Global::Transform = true;

						SRReal sinRotation = sin(glm::radians(Global::Rotation));
						SRReal cosRotation = cos(glm::radians(Global::Rotation));

						SRReal sinStretchAngle = sin(glm::radians(Global::StretchAngle));
						SRReal cosStretchAngle = cos(glm::radians(Global::StretchAngle));

						Global::TransformMatrix = glm::dmat2(cosStretchAngle, sinStretchAngle, -sinStretchAngle, cosStretchAngle);
						Global::TransformMatrix *= glm::dmat2(Global::StretchRatio, 0.0, 0.0, 1.0);
						Global::TransformMatrix *= glm::dmat2(cosStretchAngle, -sinStretchAngle, sinStretchAngle, cosStretchAngle);
						Global::TransformMatrix *= glm::dmat2(cosRotation, sinRotation, -sinRotation, cosRotation);
						Global::TransformMatrix *= glm::dmat2(1.0, 0.0, 0.0, Global::FlipVertically ? -1.0 : 1.0);

						Global::InvTransformMatrix = glm::inverse(Global::TransformMatrix);
					}

					FContext.InvalidatePixel();
				}
				[[fallthrough]];
				case IDCANCEL: {
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
			}
			break;
		}
	}

	return 0;
}

INT_PTR SetImageSizeProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM /*lParam*/) {
	static bool SetWindowSize = true;
	char Buffer[32];
	switch (message) {
		case WM_INITDIALOG: {
			CheckDlgButton(hWndDlg, IDC_FOLLOW_WINDOW_SIZE, Global::ImageSizeFollowWindowSize ? BST_CHECKED : BST_UNCHECKED);

			sprintf_s(Buffer, "%zu", FContext.ImageWidth);
			SetDlgItemTextA(hWndDlg, IDC_WIDTH, Buffer);

			sprintf_s(Buffer, "%zu", FContext.ImageHeight);
			SetDlgItemTextA(hWndDlg, IDC_HEIGHT, Buffer);

			CheckDlgButton(hWndDlg, IDC_SET_WINDOW_SIZE, SetWindowSize ? BST_CHECKED : BST_UNCHECKED);

			EnableWindow(GetDlgItem(hWndDlg, IDC_WIDTH), !Global::ImageSizeFollowWindowSize);
			EnableWindow(GetDlgItem(hWndDlg, IDC_HEIGHT), !Global::ImageSizeFollowWindowSize);
			EnableWindow(GetDlgItem(hWndDlg, IDC_SET_WINDOW_SIZE), !Global::ImageSizeFollowWindowSize);

			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDC_FOLLOW_WINDOW_SIZE: {
					if (HIWORD(wParam) == BN_CLICKED) {
						bool Enabled = IsDlgButtonChecked(hWndDlg, IDC_FOLLOW_WINDOW_SIZE) != BST_CHECKED;
						EnableWindow(GetDlgItem(hWndDlg, IDC_WIDTH), Enabled);
						EnableWindow(GetDlgItem(hWndDlg, IDC_HEIGHT), Enabled);
						EnableWindow(GetDlgItem(hWndDlg, IDC_SET_WINDOW_SIZE), Enabled);
					}
					break;
				}
				case IDOK: {
					bool ImageSizeFollowWindowSize = IsDlgButtonChecked(hWndDlg, IDC_FOLLOW_WINDOW_SIZE) == BST_CHECKED;

					if (ImageSizeFollowWindowSize && !Global::ImageSizeFollowWindowSize && (FContext.ImageWidth != WindowWidth || FContext.ImageHeight != WindowHeight)) {
						FContext.SetResolution(WindowWidth, WindowHeight);
					}

					Global::ImageSizeFollowWindowSize = ImageSizeFollowWindowSize;

					if (Global::ImageSizeFollowWindowSize) goto Cancel;

					size_t NewWidth, NewHeight;
					char Dummy;

					if (!GetDlgItemTextA(hWndDlg, IDC_WIDTH, Buffer, 32)) goto Cancel;
					if (sscanf_s(Buffer, "%zu%c", &NewWidth, &Dummy) != 1) goto Cancel;

					if (!GetDlgItemTextA(hWndDlg, IDC_HEIGHT, Buffer, 32)) goto Cancel;
					if (sscanf_s(Buffer, "%zu%c", &NewHeight, &Dummy) != 1) goto Cancel;
					
					FContext.SetResolution(NewWidth, NewHeight);

					SetWindowSize = IsDlgButtonChecked(hWndDlg, IDC_SET_WINDOW_SIZE) == BST_CHECKED;

					if (SetWindowSize) {
						RECT rect;
						rect.left = 0;
						rect.top = 0;
						rect.right = NewWidth;
						rect.bottom = NewHeight;

						AdjustWindowRectEx(&rect, WS_OVERLAPPEDWINDOW, TRUE, NULL);
						NewWidth = rect.right - rect.left;
						NewHeight = rect.bottom - rect.top;

						SetWindowPos(HWnd, nullptr, 0, 0, NewWidth, NewHeight, SWP_NOMOVE | SWP_NOZORDER | SWP_FRAMECHANGED);
					}
					
				}
				Cancel:
				[[fallthrough]];
				case IDCANCEL: {
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
			}
			break;
		}
	}

	return 0;
}
int n = 0;
INT_PTR TasksProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam) {
	switch (message) {
		case WM_INITDIALOG: {
			SetTimer(hWndDlg, 1, 100, nullptr);
		}
		[[fallthrough]];
		case WM_TIMER: {
			char Buffer[256] = "(no task in progress)";
			Buffer[255] = 0;
			auto Contexts = Computation::GetExecutionContexts();

			size_t BufferIndex = 0;

			for (auto &Context : Contexts) {
				Task *task = Context->GetTask();
				std::string_view Description = task->GetDescription();

				ProgressTrackable *TaskProgress = dynamic_cast<ProgressTrackable *>(task);

				size_t CopySize = std::min(255 - BufferIndex, Description.size());
				memcpy(Buffer + BufferIndex, Description.data(), CopySize);
				BufferIndex += CopySize;
				Buffer[BufferIndex] = 0;

				std::chrono::duration<SRReal> duration = Context->GetDuration();

				if (BufferIndex < 254) {
					Buffer[BufferIndex++] = ' ';
					Buffer[BufferIndex++] = '(';
					Buffer[BufferIndex] = 0;
				}

				if (TaskProgress) {
					SRReal Numerator, Denominator;
					if (TaskProgress->GetProgress(Numerator, Denominator)) {
						if (Denominator == 0.0_sr) {
							BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%g ", Numerator);
						} else {
							BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%.1f%% ", ((Numerator * 100.0) / Denominator));
						}
					} else {
						size_t CopySize = std::min(255 - BufferIndex, sizeof("-.-% ") - 1);
						memcpy(Buffer + BufferIndex, "-.-% ", CopySize);
						BufferIndex += CopySize;
						Buffer[BufferIndex] = 0;
					}
				}
				if (duration.count() == 0.0) {
					size_t CopySize = std::min(255 - BufferIndex, sizeof("Waiting") - 1);
					memcpy(Buffer + BufferIndex, "Waiting", CopySize);
					BufferIndex += CopySize;
					Buffer[BufferIndex] = 0;
				} if (duration.count() < 1) {
					BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%gms", round(duration.count() * 1000));
				} else if (duration.count() < 10) {
					BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%.2fs", duration.count());
				} else if (duration.count() < 60) {
					BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%.1fs", duration.count());
				} else {
					Uint Seconds = round(duration.count());
					Uint Minutes = Seconds / 60;
					Seconds %= 60;

					if (Minutes < 60) {
						BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%d:%02d", Minutes, Seconds);
					} else {
						BufferIndex += sprintf_s(Buffer + BufferIndex, 256 - BufferIndex, "%d:%02d:%02d", Minutes / 60, Minutes % 60, Seconds);
					}
				}

				if (BufferIndex < 254) {
					Buffer[BufferIndex++] = ')';
					Buffer[BufferIndex++] = '\n';
					Buffer[BufferIndex] = 0;
				}

				Context->Release();
			}

			SetDlgItemTextA(hWndDlg, IDC_TASKS, Buffer);

			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDCANCEL: {
					KillTimer(hWndDlg, 1);
					DestroyWindow(TasksDialog);
					TasksDialog = nullptr;
					return TRUE;
				}
			}
			break;
		}
	}
	return 0;
}

std::string DurationToString(std::chrono::duration<SRReal> duration) {
	char Buffer[32];
	size_t size;
	if (duration.count() < 1) {
		size = sprintf_s(Buffer, 32, "%gms", round(duration.count() * 1000));
	} else if (duration.count() < 10) {
		size = sprintf_s(Buffer, 32, "%.2fs", duration.count());
	} else if (duration.count() < 60) {
		size = sprintf_s(Buffer, 32, "%.1fs", duration.count());
	} else {
		Uint Seconds = round(duration.count());
		Uint Minutes = Seconds / 60;
		Seconds %= 60;

		if (Minutes < 60) {
			size = sprintf_s(Buffer, 32, "%d:%02d", Minutes, Seconds);
		} else {
			size = sprintf_s(Buffer, 32, "%d:%02d:%02d", Minutes / 60, Minutes % 60, Seconds);
		}
	}
	return std::string(Buffer, size);
}

INT_PTR FeatureFinderProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam) {
	static Evaluator::FeatureFinder::PreciseLocatingTask *preciseLocatingTask = nullptr;
	static ExecutionContext *preciseLocatingTaskContext = nullptr;
	static HRReal TargetRadius;
	static uint64_t Precision;
	static IntIter newItLim;
	static bool Cancelling;
	static int TargetTypeIndex = 4;
	static bool Relative = false;
	static HRReal FeatureRadius;
	static std::string TargetValueText;

	auto InitTargetValue = [](HWND hWndDlg) {
		if (TargetTypeIndex == 4) {
			SetDlgItemTextA(hWndDlg, IDC_ZOOM_TARGET_VALUE, "1.0");
		} else if (Relative) {
			SetDlgItemTextA(hWndDlg, IDC_ZOOM_TARGET_VALUE, ZoomToString(ZoomType(TargetTypeIndex), FeatureRadius / FContext.CurrentLocation.HalfH, 1.0_hr).data());
		} else {
			SetDlgItemTextA(hWndDlg, IDC_ZOOM_TARGET_VALUE, ZoomToString(ZoomType(TargetTypeIndex), FeatureRadius, 2.0_hr).data());
		}
	};

	auto CalculateTargetRadius = [](HWND hWndDlg) {
		size_t Length = GetWindowTextLengthA(GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_VALUE)) + 1;
		char *Buffer = new char[Length];
		HRReal TargetRadius;

		if (GetDlgItemTextA(hWndDlg, IDC_ZOOM_TARGET_VALUE, Buffer, Length)) {
			TargetValueText = Buffer;
			if (TargetTypeIndex == 4) {
				double Depth = std::stof(TargetValueText);
				double FeatureDepth = log2(FeatureRadius);

				if (Relative) {
					double CurrentDepth = log2(FContext.CurrentLocation.HalfH);
					FeatureDepth -= CurrentDepth;
					Depth *= FeatureDepth;
					Depth += CurrentDepth;
				} else {
					Depth *= FeatureDepth;
				}
					
				TargetRadius = pow2(Depth);
			} else {
				TargetRadius = ParseZoom(ZoomType(TargetTypeIndex), TargetValueText, Relative ? 1.0_hr : 2.0_hr);
				if (Relative) {
					TargetRadius *= FContext.CurrentLocation.HalfH;
				}
			}
		}
		delete[]Buffer;
		return TargetRadius;
	};

	Evaluator::Feature *Feature = FContext.Feature;
	switch (message) {
		case WM_INITDIALOG: {
			Cancelling = false;
			std::wstring_view name = Feature->Name();
			std::wstring_view information = Feature->Information();

			if (!name.empty()) SetWindowTextW(hWndDlg, Feature->Name().data());
			if (!information.empty()) SetDlgItemTextW(hWndDlg, IDC_INFORMATION, information.data());

			if (Feature->CanLocatePrecisely()) {
				HWND ComboBox = GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_TYPE);

				SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Height");
				SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Half height");
				SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Magnification");
				SendMessageA(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)"Depth");
				SendMessageW(ComboBox, (UINT)CB_ADDSTRING, 0, (LPARAM)L"Feature depth \xD7");

				SendMessage(ComboBox, CB_SETCURSEL, (WPARAM)TargetTypeIndex, (LPARAM)0);
				SendMessage(GetDlgItem(hWndDlg, IDC_RELATIVE), BM_SETCHECK, (WPARAM)(Relative ? BST_CHECKED : BST_UNCHECKED), (LPARAM)0);
				if (!TargetValueText.empty() && (TargetTypeIndex == 4 || Relative)) {
					SetDlgItemTextA(hWndDlg, IDC_ZOOM_TARGET_VALUE, TargetValueText.data());
				} else InitTargetValue(hWndDlg);
			} else {
				EnableWindow(GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_TYPE), FALSE);
				EnableWindow(GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_VALUE), FALSE);
				EnableWindow(GetDlgItem(hWndDlg, IDC_RELATIVE), FALSE);
				EnableWindow(GetDlgItem(hWndDlg, IDC_BUTTON_ZOOM), FALSE);
			}

			FeatureRadius = FContext.CurrentLocation.HalfH;
			FContext.Feature->GetRadius(FeatureRadius);
			break;
		}
		case WM_TIMER: {
			if (preciseLocatingTaskContext->Terminated()) {
				if (preciseLocatingTaskContext->Finished()) {
					FContext.ZoomToAnimated(preciseLocatingTask->coordinate, Precision, TargetRadius);
					if (newItLim) Global::ItLim = newItLim;
				}
				preciseLocatingTaskContext->Release();
				preciseLocatingTask = nullptr;
				preciseLocatingTaskContext = nullptr;
				EndDialog(hWndDlg, wParam);
				return TRUE;
			}
			std::stringstream stream;
			ProgressTrackable *TaskProgress = dynamic_cast<ProgressTrackable *>(preciseLocatingTask);

			SRReal Progress = 0.0;

			stream << preciseLocatingTask->GetDescription();
			if (TaskProgress) {
				SRReal Numerator, Denominator;
				if (TaskProgress->GetProgress(Numerator, Denominator) && Denominator != 0.0_sr) {
					Progress = Numerator / Denominator;
					stream << std::fixed << std::setprecision(1) << " (" << Progress * 100.0 << "%)";
				}
			}
			stream << std::endl;
			auto Duration = preciseLocatingTaskContext->GetDuration();
			stream << preciseLocatingTask->GetDetailedProgress();
			stream << "Time: " << DurationToString(Duration);
			if (Progress != 0.0_sr) {
				auto TotalTime = Duration / Progress;
				stream << " / " << DurationToString(TotalTime) << ", ";
				stream << DurationToString(TotalTime * (1.0_sr - Progress)) << " remaining.";
			}
			stream << std::endl;
			std::string text = stream.str();
			if (!text.empty()) SetDlgItemTextA(hWndDlg, IDC_INFORMATION, text.data());
			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case IDC_ZOOM_TARGET_TYPE: {
					if (HIWORD(wParam) == CBN_SELCHANGE) {
						TargetTypeIndex = SendMessage((HWND)lParam, (UINT)CB_GETCURSEL, (WPARAM)0, (LPARAM)0);
						InitTargetValue(hWndDlg);
					}
					break;
				}
				case IDC_RELATIVE: {
					if (HIWORD(wParam) == BN_CLICKED) {
						Relative = SendMessage(GetDlgItem(hWndDlg, IDC_RELATIVE), BM_GETCHECK, (WPARAM)0, (LPARAM)0) == BST_CHECKED;
						InitTargetValue(hWndDlg);
					}
					break;
				}
				case IDC_BUTTON_CENTER: {
					RelLocation NewLocation = FContext.CurrentLocation;
					FContext.Feature->GetCoordinate(NewLocation.X, NewLocation.Y);
					FContext.ChangeLocation(NewLocation);
					EndDialog(hWndDlg, wParam);
					return TRUE;
				}
				case IDC_BUTTON_ZOOM: {
					try {
						TargetRadius = CalculateTargetRadius(hWndDlg);
					} catch (std::invalid_argument) {
						MessageBoxA(hWndDlg, "Zoom target value invalid.", nullptr, MB_OK | MB_ICONERROR);
						break;
					} catch (std::out_of_range) {
						MessageBoxA(hWndDlg, "Zoom target value out of range.", nullptr, MB_OK | MB_ICONERROR);
						break;
					}
					Precision = -std::min(0ll, TargetRadius.Exponent) + 64;

					Evaluator::FeatureFinder *featureFinder = FContext.evaluator->GetFeatureFinder();
					newItLim = FContext.Feature->ItLimForZoomLevel(TargetRadius);
					if (!featureFinder) goto failed;

					using PreciseLocatingTask = Evaluator::FeatureFinder::PreciseLocatingTask;
					preciseLocatingTask = featureFinder->CreatePreciseLocatingTask(FContext.Feature, Precision, FContext.CenterCoordinate);

					if (!preciseLocatingTask) goto failed;
					FContext.CancelAll();
					preciseLocatingTaskContext = Computation::AddTask(preciseLocatingTask);

					using namespace std;

					preciseLocatingTaskContext->Wait(100ms);

					if (preciseLocatingTaskContext->Finished()) {
						FContext.ZoomToAnimated(preciseLocatingTask->coordinate, Precision, TargetRadius);
						preciseLocatingTaskContext->Release();
						preciseLocatingTask = nullptr;
						preciseLocatingTaskContext = nullptr;
						if (newItLim) Global::ItLim = newItLim;
					} else {
						EnableWindow(GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_TYPE), FALSE);
						EnableWindow(GetDlgItem(hWndDlg, IDC_ZOOM_TARGET_VALUE), FALSE);
						EnableWindow(GetDlgItem(hWndDlg, IDC_RELATIVE), FALSE);
						EnableWindow(GetDlgItem(hWndDlg, IDC_BUTTON_CENTER), FALSE);
						EnableWindow(GetDlgItem(hWndDlg, IDC_BUTTON_ZOOM), FALSE);
						SetTimer(hWndDlg, 0, 100, nullptr);
						break;
					}
				}
				failed:
				[[fallthrough]];
				case IDCANCEL: {
					if (preciseLocatingTaskContext) {
						if (preciseLocatingTaskContext->Terminated()) {
							preciseLocatingTaskContext->Release();
							preciseLocatingTaskContext = nullptr;
							EndDialog(hWndDlg, wParam);
							return TRUE;
						} else {
							preciseLocatingTaskContext->Cancel();
							EnableWindow(GetDlgItem(hWndDlg, IDCANCEL), FALSE);
							Cancelling = true;
						}
					} else {
						EndDialog(hWndDlg, wParam);
						return TRUE;
					}
					break;
				}
			}
			break;
		}
	}
	return 0;
}

LRESULT CALLBACK WindowProcess(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam) {
	switch (Message) {
		case WM_LBUTTONDOWN: {
			MouseX = (int32_t)LOWORD(lParam);
			MouseY = (int32_t)HIWORD(lParam);
			Panning = true;
			break;
		}
		case WM_LBUTTONUP: {
			MouseX = (int32_t)LOWORD(lParam);
			MouseY = (int32_t)HIWORD(lParam);
			Panning = false;
			break;
		}
		case WM_LBUTTONDBLCLK: {
			if (FContext.Feature) {
				DialogBox(nullptr, MAKEINTRESOURCE(IDD_FEATURE_FINDER), HWnd, FeatureFinderProc);
			}
			break;
		}
		case WM_MOUSEMOVE: {
			if (!(wParam & MK_LBUTTON)) Panning = false;
			if (Panning) {
				double MovementX = -double((int32_t)LOWORD(lParam) - MouseX) / WindowHeight * 2.0;
				double MovementY = double((int32_t)HIWORD(lParam) - MouseY) / WindowHeight * 2.0;

				double FractalMovementX = Global::InvTransformMatrix[0][0] * MovementX + Global::InvTransformMatrix[1][0] * MovementY;
				double FractalMovementY = Global::InvTransformMatrix[0][1] * MovementX + Global::InvTransformMatrix[1][1] * MovementY;

				Global::moved = true;
				FContext.Move(FractalMovementX, FractalMovementY);
			}
			MouseX = (int32_t)LOWORD(lParam);
			MouseY = (int32_t)HIWORD(lParam);
			Global::Redraw = true;
			break;
		}
		case WM_MOUSEWHEEL: {
			POINT MousePoint{ (short)LOWORD(lParam), (short)HIWORD(lParam) };
			ScreenToClient(hWnd, &MousePoint);
			SRReal X = (SRReal((int)MousePoint.x) - WindowWidth * 0.5) / WindowHeight * 2.0;
			SRReal Y = -(SRReal((int)MousePoint.y) / WindowHeight * 2.0 - 1.0);

			SRReal FractalX = Global::InvTransformMatrix[0][0] * X + Global::InvTransformMatrix[1][0] * Y;
			SRReal FractalY = Global::InvTransformMatrix[0][1] * X + Global::InvTransformMatrix[1][1] * Y;

			int delta = GET_WHEEL_DELTA_WPARAM(wParam);

			if (delta > 0) {
				FContext.ZoomIn(FractalX, FractalY);
			} else if (delta < 0) {
				FContext.ZoomOut(FractalX, FractalY);
			}

			break;
		}
		case WM_KEYDOWN: {
			if (GetKeyState(VK_CONTROL) & 0x8000) {
				switch (wParam) {
					case 'S': {
						if (GetKeyState(VK_SHIFT) & 0x8000) {
							SaveImage();
						} else {
							SaveLocation();
						}
						break;
					}
					case 'O': {
						OpenLocation();
						break;
					}
				}
			} else {
				switch (wParam) {
					case 'A': Global::ItDiv = std::max(8.0_sr, Global::ItDiv / 2); Global::Redraw = true; break;
					case 'S': Global::ItDiv = std::min(0x1p48_sr, Global::ItDiv * 2); Global::Redraw = true; break;
					case 'E': Global::ColorCyclingSpeed = -0.1; Global::Redraw = true; break;
					case 'R': Global::ColorCyclingSpeed = 0.1; Global::Redraw = true; break;
					case 'D': Global::ItLim = std::min(UINT64_C(1) << 48, Global::ItLim * 2); FContext.ParameterChanged = true; break;
					case 'F': Global::ItLim = std::max(UINT64_C(8), Global::ItLim / 2); FContext.ParameterChanged = true; break;
					default:
						break;
				}
			}
			break;
		}
		case WM_KEYUP: {
			switch (wParam) {
				case 'E':
				case 'R': Global::ColorCyclingSpeed = 0.0; Global::Redraw = true; break;
			}
			break;
		}
		case WM_COMMAND: {
			switch (LOWORD(wParam)) {
				case (UINT_PTR)MenuID::Open: {
					OpenLocation();
					break;
				}
				case (UINT_PTR)MenuID::Save: {
					SaveLocation();
					break;
				}
				case (UINT_PTR)MenuID::SaveImage: {
					SaveImage();
					break;
				}
				case (UINT_PTR)MenuID::SaveRawPixelData: {
					SaveImage(true);
					break;
				}
				case (UINT_PTR)MenuID::FractalTypeMandelbrot: {
					SetFractalType(FractalTypeEnum::Mandelbrot);
					break;
				}
				case (UINT_PTR)MenuID::FractalTypeTricorn: {
					SetFractalType(FractalTypeEnum::Tricorn);
					break;
				}
				case (UINT_PTR)MenuID::FractalTypeBurningShip: {
					SetFractalType(FractalTypeEnum::BurningShip);
					break;
				}
				case (UINT_PTR)MenuID::FractalTypeNova: {
					SetFractalType(FractalTypeEnum::Nova);
					break;
				}
				case (UINT_PTR)MenuID::FractalTypeCustom: {
					DialogBox(nullptr, MAKEINTRESOURCE(IDD_CUSTOM_FORMULA), HWnd, SetFormulaProc);
					break;
				}
				case (UINT_PTR)MenuID::DistanceEstimation: {
					if (CurrentFractalType != FractalTypeEnum::Nova && (CurrentFractalType != FractalTypeEnum::Mandelbrot || !UseLinearApproximation)) break;
					UseDE = !UseDE;
					CheckMenuItem(Fractal, (UINT_PTR)MenuID::DistanceEstimation, MF_BYCOMMAND | (UseDE ? MF_CHECKED : MF_UNCHECKED));
					FContext.SetUseDE(UseDE);
					break;
				}
				case (UINT_PTR)MenuID::HighQuality: {
					if (CurrentFractalType != FractalTypeEnum::Mandelbrot || !UseLinearApproximation) break;
					Global::HighQuality = !Global::HighQuality;
					CheckMenuItem(Fractal, (UINT_PTR)MenuID::HighQuality, MF_BYCOMMAND | (Global::HighQuality ? MF_CHECKED : MF_UNCHECKED));
					FContext.InvalidatePixel();
					break;
				}
				case (UINT_PTR)MenuID::Location: {
					DialogBox(nullptr, MAKEINTRESOURCE(IDD_LOCATION), HWnd, SetLocationProc);
					break;
				}
				case (UINT_PTR)MenuID::IterationLimit: {
					DialogBox(nullptr, MAKEINTRESOURCE(IDD_ITERATION_LIMIT), HWnd, SetIterationLimitProc);
					break;
				}
				case (UINT_PTR)MenuID::Transformation: {
					DialogBox(nullptr, MAKEINTRESOURCE(IDD_TRANSFORMATION), HWnd, TransformProc);
					break;
				}
				case (UINT_PTR)MenuID::IncreaseColorDensity: {
					Global::ItDiv = std::max(8.0_sr, Global::ItDiv / 2);
					if (Global::PreModulo && Global::ItDiv * Global::MaxPreModMultiplier < Global::ItMod) {
						Global::ItMod = Global::ItDiv * Global::DefaultPreModMultiplier;
						FContext.InvalidatePixel();
					}
					break;
				}
				case (UINT_PTR)MenuID::DecreaseColorDensity: {
					Global::ItDiv = std::min(0x1p48_sr, Global::ItDiv * 2);
					if (Global::PreModulo && Global::ItDiv > Global::ItMod) {
						Global::ItMod = Global::ItDiv * Global::DefaultPreModMultiplier;
						FContext.InvalidatePixel();
					}
					break;
				}
				case (UINT_PTR)MenuID::IncreaseIterationLimit: {
					Global::ItLim = std::min(UINT64_C(1) << 48, Global::ItLim * 2); FContext.ParameterChanged = true;
					break;
				}
				case (UINT_PTR)MenuID::DecreaseIterationLimit: {
					Global::ItLim = std::max(UINT64_C(8), Global::ItLim / 2); FContext.ParameterChanged = true;
					break;
				}
				case (UINT_PTR)MenuID::ResetLocation: {
					FContext.ResetLocation();
					break;
				}
				case (UINT_PTR)MenuID::AlgorithmPerturbation: {
					if (CurrentFractalType != FractalTypeEnum::Mandelbrot) break;
					EnableMenuItem(Fractal, (UINT_PTR)MenuID::DistanceEstimation, MF_BYCOMMAND | MF_GRAYED);
					EnableMenuItem(Fractal, (UINT_PTR)MenuID::HighQuality, MF_BYCOMMAND | MF_GRAYED);
					CheckMenuRadioItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmBegin, (UINT_PTR)MenuID::AlgorithmEnd, (UINT_PTR)MenuID::AlgorithmPerturbation, MF_BYCOMMAND);
					FContext.ChangeFractalType(CurrentFractalType, false);
					UseLinearApproximation = false;
					break;
				}
				case (UINT_PTR)MenuID::AlgorithmLinearApproximation: {
					if (CurrentFractalType != FractalTypeEnum::Mandelbrot) break;
					EnableMenuItem(Fractal, (UINT_PTR)MenuID::DistanceEstimation, MF_BYCOMMAND | MF_ENABLED);
					EnableMenuItem(Fractal, (UINT_PTR)MenuID::HighQuality, MF_BYCOMMAND | MF_ENABLED);
					CheckMenuRadioItem(AlgorithmMenu, (UINT_PTR)MenuID::AlgorithmBegin, (UINT_PTR)MenuID::AlgorithmEnd, (UINT_PTR)MenuID::AlgorithmLinearApproximation, MF_BYCOMMAND);
					FContext.ChangeFractalType(CurrentFractalType, true);
					UseLinearApproximation = true;
					break;
				}
				case (UINT_PTR)MenuID::LockReference: {
					FContext.LockReference = !FContext.LockReference;
					CheckMenuItem(Computation, (UINT_PTR)MenuID::LockReference, MF_BYCOMMAND | (FContext.LockReference ? MF_CHECKED : MF_UNCHECKED));
					break;
				}
				case (UINT_PTR)MenuID::Tasks: {
					if (!TasksDialog) {
						TasksDialog = CreateDialogW(nullptr, MAKEINTRESOURCE(IDD_TASKS), HWnd, TasksProc);
					}
					ShowWindow(TasksDialog, SW_SHOW);
					break;
				}
				case (UINT_PTR)MenuID::RecomputeReference: {
					FContext.InvalidateReference();
					break;
				}
				case (UINT_PTR)MenuID::RecomputePixel: {
					FContext.InvalidatePixel();
					break;
				}
				case (UINT_PTR)MenuID::RecomputeAll: {
					FContext.InvalidateAll();
					break;
				}
				case (UINT_PTR)MenuID::ImageSize: {
					DialogBox(nullptr, MAKEINTRESOURCE(IDD_IMAGE_SIZE), HWnd, SetImageSizeProc);
					break;
				}
				case (UINT_PTR)MenuID::PaletteMipmaps: {
					UsePalleteMipmaps = !UsePalleteMipmaps;
					EnablePaletteMipmap(UsePalleteMipmaps);
					CheckMenuItem(Image, (UINT_PTR)MenuID::PaletteMipmaps, MF_BYCOMMAND | (UsePalleteMipmaps ? MF_CHECKED : MF_UNCHECKED));
					break;
				}
				case (UINT_PTR)MenuID::BilinearFilter: {
					Global::UseBilinearFilter = !Global::UseBilinearFilter;
					CheckMenuItem(Image, (UINT_PTR)MenuID::BilinearFilter, MF_BYCOMMAND | (Global::UseBilinearFilter ? MF_CHECKED : MF_UNCHECKED));
					Global::Redraw = true;
					break;
				}
				case (UINT_PTR)MenuID::FlipVertically: {
					Global::FlipVertically = !Global::FlipVertically;
					CheckMenuItem(Image, (UINT_PTR)MenuID::FlipVertically, MF_BYCOMMAND | (Global::FlipVertically ? MF_CHECKED : MF_UNCHECKED));
					Global::Redraw = true;
					break;
				}
				case (UINT_PTR)MenuID::PreModulo: {
					Global::PreModulo = !Global::PreModulo;
					CheckMenuItem(Image, (UINT_PTR)MenuID::PreModulo, MF_BYCOMMAND | (Global::PreModulo ? MF_CHECKED : MF_UNCHECKED));
					FContext.InvalidatePixel();
					break;
				}
			}
			break;
		}
		case WM_ENTERSIZEMOVE: {
			SizeMove = true;
			return 0;
		}
		case WM_SIZE: {
			if (wParam == SIZE_MINIMIZED) break;
			if (WindowWidth == LOWORD(lParam) && WindowHeight == HIWORD(lParam)) break;
			WindowWidth = LOWORD(lParam);
			WindowHeight = HIWORD(lParam);
			if (SizeMove) {
				SizeChanged = true;
			} else {
				ResizeViewport(WindowWidth, WindowHeight);
				if (Global::ImageSizeFollowWindowSize) {
					FContext.SetResolution(WindowWidth, WindowHeight);
				}
			}
			break;
		}
		case WM_PAINT: {
			if (SizeMove && SizeChanged && Global::Initialized) {
				ResizeViewport(WindowWidth, WindowHeight);
				BeginRender();
				RenderFractal(FContext);
				EndRender();
			} else {
				Global::Redraw = true;
			}
			return DefWindowProc(hWnd, Message, wParam, lParam);
		}
		case WM_EXITSIZEMOVE: {
			if (SizeChanged && Global::ImageSizeFollowWindowSize) {
				ResizeViewport(WindowWidth, WindowHeight);
				FContext.SetResolution(WindowWidth, WindowHeight);
				SizeChanged = false;
			}
			SizeMove = false;
			return 0;
		}
		case WM_DESTROY: {
			PostQuitMessage(0);
			break;
		}
		default: {
			return DefWindowProc(hWnd, Message, wParam, lParam);
		}
	}
	return 0;
}