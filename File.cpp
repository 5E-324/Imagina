#include "Includes.h"
#include "File.h"

#include <fstream>
#include <png.h>
#include <gmp-impl.h>
#include <stack>

// Temporary file format
constexpr uint64_t IMMagicNumber = 0x000A0D56504D49FF;

struct IMFileHeader {
	uint64_t Magic;
	uint64_t Reserved;
	uint64_t LocationOffset;
	uint64_t ReferenceOffset;
};

extern bool UseDE;

size_t mpz_inp_raw_stream(mpz_ptr x, std::istream &fp) {
	unsigned char  csize_bytes[4];
	mpir_out_struct out;

	/* 4 bytes for size */
	fp.read((char *)csize_bytes, sizeof(csize_bytes));

	mpz_inp_raw_p(x, csize_bytes, out);

	if (out->writtenSize != 0) {
		fp.read(out->written, out->writtenSize);

		mpz_inp_raw_m(x, out);
	}
	return out->writtenSize + 4;
}
size_t mpz_out_raw_stream(std::ostream &fp, mpz_srcptr x) {
	mpir_out_struct out;

	mpz_out_raw_m(out, x);

	fp.write(out->written, out->writtenSize);

	void (*gmp_free_func)(void *, size_t);

	mp_get_memory_functions(nullptr, nullptr, &gmp_free_func);
	(*gmp_free_func) (out->allocated, out->allocatedSize);
	return out->writtenSize;
}
size_t mpf_out_raw_stream(std::ostream &f, mpf_srcptr X) {
	long int expt; mpz_t Z; int nz;
	expt = X->_mp_exp;
	f.write((char *)&expt, sizeof(long int));
	nz = X->_mp_size;
	Z->_mp_alloc = std::abs(nz);
	Z->_mp_size = nz;
	Z->_mp_d = X->_mp_d;
	return mpz_out_raw_stream(f, Z) + sizeof(int);
}
void mpf_inp_raw_stream(std::istream &f, mpf_ptr X) {
	long int expt; mpz_t Z;
	mpz_init(Z);
	f.read((char *)&expt, sizeof(long int));
	mpz_inp_raw_stream(Z, f);
	mpf_set_z(X, Z);
	X->_mp_exp = expt;
	mpz_clear(Z);
}
int mpf_out_raw_stream(std::ostream &f, const mpf_class &X) {
	return mpf_out_raw_stream(f, X.get_mpf_t());
}
void mpf_inp_raw_stream(std::istream &f, mpf_class &X) {
	return mpf_inp_raw_stream(f, X.get_mpf_t());
}

FileType DetermineFileType(std::istream &File, std::wstring_view Extension) {
	using enum FileType;
	char temp[64];
	File.seekg(0, std::ios::beg);
	File.get(temp, 64);
	File.seekg(0, std::ios::beg);

	uint64_t Magic;
	File.read((char *)&Magic, sizeof(uint64_t));
	File.seekg(0, std::ios::beg);

	FileType Type = Other;

	if (Magic == IMMagicNumber) {
		Type = Imagina;
	} else {
		for (char *p = temp; *p; p++) if (*p == ':') Type = ImaginaText;
		if (Type == ImaginaText) {
			if (Extension == L"kfp") Type = Kfp;
			if (Extension == L"kfr") Type = Kfr;
		}
	}

	return Type;
}

void LoadImFile(std::istream &File) {
	FContext.InvalidatePixel();
	IMFileHeader Header;
	File.read((char *)&Header, sizeof(IMFileHeader));

	File.seekg(Header.LocationOffset);

	File.read((char *)&FContext.CurrentLocation.HalfH, sizeof(HRReal));
	File.read((char *)&Global::ItLim, sizeof(size_t));

	FContext.CurrentLocation.X = 0.0_hr;
	FContext.CurrentLocation.Y = 0.0_hr;
	FContext.RenderLocation = FContext.CurrentLocation;
	FContext.EvalLocation = FContext.CurrentLocation;

	uint64_t Precision = -std::min(0ll, FContext.CurrentLocation.HalfH.Exponent) + 64;
	SetDefaultPrecision(Precision);
	FContext.CenterCoordinate.X.set_prec(Precision);
	FContext.CenterCoordinate.Y.set_prec(Precision);
	FContext.pixelManager.Clear();
	FContext.ComputePixel = true;

	mpf_inp_raw_stream(File, FContext.CenterCoordinate.X);
	mpf_inp_raw_stream(File, FContext.CenterCoordinate.Y);

	if (Header.ReferenceOffset) {
		File.seekg(Header.ReferenceOffset);

		//FContext.evaluator->LoadReference(OpenFile, FContext.reference);
		if (FContext.reference) FContext.evaluator->FreeReference(FContext.reference);
		FContext.reference = nullptr;

		FContext.evaluator->LoadReference(File, FContext.reference);
	} else {
		FContext.InvalidateReference();
	}
}

const std::unordered_map<std::string, FractalTypeEnum> FormulaNameMap = {
	{ "Mandelbrot", FractalTypeEnum::Mandelbrot },
	{ "Tricorn", FractalTypeEnum::Tricorn },
	{ "Burning ship", FractalTypeEnum::BurningShip },
	{ "Nova", FractalTypeEnum::Nova },
};

// Workaround for a bug of MPIR
void CorrectMpf(mpf_ptr mp) {
	bool sign = mp->_mp_size < 0;
	mp->_mp_size = abs(mp->_mp_size);
	while (mp->_mp_size && !mp->_mp_d[mp->_mp_size - 1]) { mp->_mp_size--; mp->_mp_exp--; }
	mp->_mp_size = sign ? -mp->_mp_size : mp->_mp_size;
}

void LoadImtFile(std::istream &File) {
	FContext.InvalidateAll();
	std::unordered_map<std::string, std::string> Parameters;
	std::string KeyAndIndentation, Key, Value;
	struct BlockInfo {
		std::string Indentation;
		std::string Prefix;
	};
	std::stack<BlockInfo> Blocks;
	while (std::getline(File, KeyAndIndentation, ':').good() && std::getline(File, Value).good()) {
		size_t IndentationSize = 0;
		while (IndentationSize < KeyAndIndentation.size() && (KeyAndIndentation[IndentationSize] == ' ' || KeyAndIndentation[IndentationSize] == '\t')) {
			IndentationSize++;
		}
		std::string_view Indentation = std::string_view(KeyAndIndentation).substr(0, IndentationSize);
		std::string_view BlockIndentation;
		if (Blocks.empty()) {
			BlockIndentation = std::string_view();
		} else {
			BlockIndentation = std::string_view(Blocks.top().Indentation);
		}
		size_t BlockIndentationSize = BlockIndentation.size();
		if (IndentationSize > BlockIndentationSize) {
			if (Indentation.substr(0, BlockIndentationSize) != BlockIndentation) throw;
			Blocks.push(BlockInfo{ std::string(Indentation), Key + '.' });
		} else if (IndentationSize < BlockIndentationSize) {
			while (!Blocks.empty() && Blocks.top().Indentation.size() > IndentationSize) Blocks.pop();
			if (Blocks.empty()) {
				if (Indentation.size() != 0) throw;
			} else if (Indentation != Blocks.top().Indentation) throw;
		} else if (Indentation != BlockIndentation) {
			throw;
		}

		size_t WhitespaceCount = 0;
		while (WhitespaceCount < Value.size() && (Value[WhitespaceCount] == ' ' || Value[WhitespaceCount] == '\t')) {
			WhitespaceCount++;
		}
		Key = Blocks.empty() ? std::string() : Blocks.top().Prefix;
		Key.append(std::string_view(KeyAndIndentation).substr(IndentationSize));
		Parameters.try_emplace(Key, Value.substr(WhitespaceCount));
	}
	if (auto FormulaText = Parameters.find("Formula"); FormulaText != Parameters.end()) {
		if (auto Formula = FormulaNameMap.find(FormulaText->second); Formula != FormulaNameMap.end()) {
		} else {
			Global::CustomFormula = FormulaText->second;
			SetFractalType(FractalTypeEnum::Custom);
		}
	}
	// TODO: Restore previous values in case of exception
	if (Parameters.contains("Location")) {
		try {
			HRReal HalfH = HRReal(mpf_class(Parameters.at("Location.Size"), 52, -16));

			if (HalfH > 16.0_hr) HalfH = 16.0_hr;
			uint64_t Precision = -std::min(0ll, HalfH.Exponent) + 64;

			Coordinate NewCoordinate{ mpf_class(Parameters.at("Location.Re"), Precision, 16), mpf_class(Parameters.at("Location.Im"), Precision, -16) };

			CorrectMpf(NewCoordinate.X.get_mpf_t()); // Workaround for a bug of MPIR
			CorrectMpf(NewCoordinate.Y.get_mpf_t());

			//FContext.ZoomToAnimated(NewCoordinate, Precision, HalfH);
			FContext.SetLocation(NewCoordinate, Precision, HalfH);

			size_t NewItLim;
			char Dummy;
			if (sscanf_s(Parameters.at("Location.Iterations").c_str(), "%zu%c", &NewItLim, &Dummy) == 1) {
				Global::ItLim = std::max(UINT64_C(2), std::min(UINT64_C(1) << 32, NewItLim));
			}

		} catch (std::out_of_range) {
			ErrorMessage("Error while loading location", "Argument missing");
		} catch (std::invalid_argument) {
			ErrorMessage("Error while loading location", "Invalid argument");
		}
	}
}

void LoadKfFile(std::istream &File, bool IsKfr) {
	FContext.InvalidateAll();
	std::unordered_map<std::string, std::string> Parameters;
	std::string Key, Value;
	while (std::getline(File, Key, ':').good() && std::getline(File, Value).good()) {
		Parameters.try_emplace(Key, Value);
	}
	// TODO: Restore previous values in case of exception
	if (IsKfr) {
		if (auto Zoom = Parameters.find("Zoom"); Zoom != Parameters.end()) {
			try {
				HRReal HalfH = 2.0_hr / HRReal(mpf_class(Zoom->second, 52));

				if (HalfH > 16.0_hr) HalfH = 16.0_hr;
				uint64_t Precision = -std::min(0ll, HalfH.Exponent) + 64;

				Coordinate NewCoordinate{ mpf_class(Parameters.at("Re"), Precision), mpf_class(Parameters.at("Im"), Precision) };

				CorrectMpf(NewCoordinate.X.get_mpf_t()); // Workaround for a bug of MPIR
				CorrectMpf(NewCoordinate.Y.get_mpf_t());

				FContext.SetLocation(NewCoordinate, Precision, HalfH);

				size_t NewItLim;
				char Dummy;
				if (sscanf_s(Parameters.at("Iterations").c_str(), "%zu%c", &NewItLim, &Dummy) == 1) {
					Global::ItLim = std::max(UINT64_C(2), std::min(UINT64_C(1) << 32, NewItLim));
				}

			} catch (std::out_of_range) { // Parameter not found
				// TODO
			} catch (std::invalid_argument) {

			}
		}
	}

	if (auto Colors = Parameters.find("Colors"); Colors != Parameters.end()) {
		const char *ColorsString = Colors->second.c_str();
		int NumberOfArguments, NumberOfCharacters;
		RGBA Color;
		Color.A = 1.0;
		Global::Palette.clear();
		while ((NumberOfArguments = sscanf_s(ColorsString, "%f,%f,%f,%n", &Color.B, &Color.G, &Color.R, &NumberOfCharacters)) == 3) {
			Color.R *= 1.0f / 255.0f;
			Color.G *= 1.0f / 255.0f;
			Color.B *= 1.0f / 255.0f;
			Global::Palette.push_back(Color);
			ColorsString += NumberOfCharacters;
		}
		if (Global::Palette.size() == 0) { // TODO: Restore default palette
			Global::Palette.push_back(RGBA{ 0.0, 0.0, 0.0, 1.0 });
		}
		Global::GammaCorrection = false;
		SetPalette(Global::Palette);
	}
	if (auto IterDiv = Parameters.find("IterDiv"); IterDiv != Parameters.end()) {
		double ItDiv;
		sscanf_s(IterDiv->second.c_str(), "%lf", &ItDiv);
		Global::ItDiv = ItDiv * (UseDE ? 8192.0 : 1024.0);
	}
}

void OpenFile(wchar_t *FileName, size_t ExtensionOffset) {
	using enum FileType;
	std::ifstream OpenFile(FileName, std::ifstream::in | std::ifstream::binary);

	FileType Type = DetermineFileType(OpenFile, ExtensionOffset ? FileName + ExtensionOffset : L"");

	if (Type != Imagina) {
		OpenFile.close();
		OpenFile.open(FileName, std::ifstream::in);
	}

	switch (Type) {
		case Imagina:	  LoadImFile(OpenFile);		   break;
		case ImaginaText: LoadImtFile(OpenFile);	   break;
		case Kfr:		  LoadKfFile(OpenFile, true);  break;
		case Kfp:		  LoadKfFile(OpenFile, false); break;
		default: {
			MessageBoxW(HWnd, L"File format not recognized.", nullptr, MB_OK | MB_ICONERROR);
			break;
		}
	}

	OpenFile.close();
}

void SaveKfrFile(std::ostream &File) {

	size_t NumberOfDigits = (size_t)std::max(0.0, -log10(FContext.CurrentLocation.HalfH)) + 10;

	File.precision(NumberOfDigits);
	File << std::fixed;

	File << "Re: " << (FContext.CenterCoordinate.X + HPReal(FContext.CurrentLocation.X)) << std::endl;
	File << "Im: " << (FContext.CenterCoordinate.Y + HPReal(FContext.CurrentLocation.Y)) << std::endl;

	File.precision(12);
	File << std::defaultfloat;

	File << "Zoom: " << (2.0_hr / FContext.CurrentLocation.HalfH).to_mpf_class(53) << std::endl;
	File << "Iterations: " << Global::ItLim << std::endl;

	File << "IterDiv: " << Global::ItDiv / (UseDE ? 8192.0 : 1024.0) << std::endl;

	if (UseDE) {
		File << "ColorMethod: 7" << std::endl;
		File << "Differences: 7" << std::endl;
	} else {
		File << "ColorMethod: 0" << std::endl;
		File << "Differences: 0" << std::endl;
	}
}

// Temporary file format
void SaveImFile(std::ostream &File) {
	IMFileHeader Header = {};
	Header.Magic = IMMagicNumber;
	File.write((char *)&Header, sizeof(IMFileHeader));
	Header.LocationOffset = File.tellp();

	File.write((char *)&FContext.CurrentLocation.HalfH, sizeof(HRReal));
	File.write((char *)&Global::ItLim, sizeof(size_t));
	mpf_out_raw_stream(File, FContext.CenterCoordinate.X);
	mpf_out_raw_stream(File, FContext.CenterCoordinate.Y);

	uint64_t ReferenceOffset = File.tellp();
	if (FContext.evaluator->SaveReference(File, FContext.reference)) { // FIXME
		Header.ReferenceOffset = ReferenceOffset;
	} else {
		Header.ReferenceOffset = 0;
	}

	File.seekp(0);
	File.write((char *)&Header, sizeof(IMFileHeader));
}

void SaveImtFile(std::ostream &File) {
	size_t NumberOfDigits = (size_t)std::max<int64_t>(0, -FContext.CurrentLocation.HalfH.Exponent / 4) + 10;

	File << "Formula: ";
	switch (CurrentFractalType) {
		case FractalTypeEnum::Mandelbrot: {
			File << "Mandelbrot";
			break;
		}
		case FractalTypeEnum::Tricorn: {
			File << "Tricorn";
			break;
		}
		case FractalTypeEnum::BurningShip: {
			File << "Burning ship";
			break;
		}
		case FractalTypeEnum::Nova: {
			File << "Nova";
			break;
		}
		case FractalTypeEnum::Custom: {
			File << Global::CustomFormula;
			break;
		}
		default: {
			break;
		}
	}
	File << std::endl;

	File.precision(NumberOfDigits);
	//File << std::hexfloat;
	File.setf(std::ios_base::hex | std::ios_base::fixed, std::ios_base::basefield | std::ios_base::floatfield);
	File << "Location:\n";
	File << "\tRe: " << (FContext.CenterCoordinate.X + HPReal(FContext.CurrentLocation.X)) << std::endl;
	File << "\tIm: " << (FContext.CenterCoordinate.Y + HPReal(FContext.CurrentLocation.Y)) << std::endl;

	File.precision(12);
	File.unsetf(std::ios_base::floatfield);
	File << "\tSize: " << (FContext.CurrentLocation.HalfH).to_mpf_class(53) << std::endl;

	File.unsetf(std::ios_base::basefield);
	File << "\tIterations: " << Global::ItLim << std::endl;
}

void SaveFile(wchar_t *FileName, FileType Type) {
	using enum FileType;
	switch (Type) {
		case Imagina: {
			std::ofstream File(FileName, std::ofstream::out | std::ofstream::binary);
			SaveImFile(File);
			File.close();
			break;
		}
		case ImaginaText: {
			std::ofstream File(FileName, std::ofstream::out);
			SaveImtFile(File);
			File.close();
			break;
		}
		case Kfr: {
			std::ofstream File(FileName, std::ofstream::out);
			SaveKfrFile(File);
			File.close();
			break;
		}
	}
}

void SaveImage(wchar_t *FileName) {
	int32_t width = FContext.pixelManager.Width, height = FContext.pixelManager.Height;

	FILE *fp;
	_wfopen_s(&fp, FileName, L"wb");
	if (!fp) abort();

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) abort();

	png_infop info = png_create_info_struct(png);
	if (!info) abort();

	if (setjmp(png_jmpbuf(png))) abort();

	png_init_io(png, fp);

	// Output is 8bit depth, RGBA format.
	png_set_IHDR(
		png,
		info,
		width, height,
		8,
		PNG_COLOR_TYPE_RGBA,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT
	);
	png_write_info(png, info);

	png_bytep *row_pointers = new png_bytep[height];
	png_bytep pixData = new png_byte[width * height * 4]{ 128 };
	memset(pixData, 128, width * height * 4);

	for (int32_t i = 0; i < height; i++) {
		row_pointers[height - 1 - i] = pixData + i * width * 4;
	}

	struct RGB {
		float r, g, b;
	};
	float ItMul = 1.0 / Global::ItDiv;
	float FMaxIt = Global::MaxIt;
	float PaletteSize = Global::Palette.size();
	size_t PaletteSizeI = Global::Palette.size();
	auto ValToRGBA = [&](float val) {
		if (val >= FMaxIt) return RGBA{ 0.0, 0.0, 0.0, 1.0 };
		val *= ItMul;
		val += Global::ColoringValueOffset;
		val -= floor(val);
		val *= PaletteSize;
		size_t index1 = floor(val);
		size_t index2 = index1 + 1;
		if (index2 == PaletteSizeI) index2 = 0;
		val -= floor(val);
		float weight1 = 1.0 - val;
		float weight2 = val;
		RGBA result;
		result.R = Global::Palette[index1].R * weight1 + Global::Palette[index2].R * weight2;
		result.G = Global::Palette[index1].G * weight1 + Global::Palette[index2].G * weight2;
		result.B = Global::Palette[index1].B * weight1 + Global::Palette[index2].B * weight2;
		result.A = Global::Palette[index1].A * weight1 + Global::Palette[index2].A * weight2;
		return result;
	};
	size_t X = FContext.pixelManager.Passes[FContext.pixelManager.PassCount - 1].X;
	size_t Y = FContext.pixelManager.Passes[FContext.pixelManager.PassCount - 1].Y;
	size_t TextureWidth = FContext.pixelManager.TextureWidth;
	for (size_t y = Y, i = 0; y < Y + height; y++) {
		for (size_t x = X; x < X + width; x++, i++) {
			RGBA rgba = ValToRGBA(FContext.pixelManager.Data[y * TextureWidth + x]);
			if (Global::GammaCorrection) {
				rgba.R = sqrt(rgba.R);
				rgba.G = sqrt(rgba.G);
				rgba.B = sqrt(rgba.B);
			}
			pixData[i * 4] = round(rgba.R * 255.0);
			pixData[i * 4 + 1] = round(rgba.G * 255.0);
			pixData[i * 4 + 2] = round(rgba.B * 255.0);
			pixData[i * 4 + 3] = round(rgba.A * 255.0);
		}
	}

	png_write_image(png, row_pointers);
	png_write_end(png, nullptr);

	delete[] row_pointers;
	delete[] pixData;

	fclose(fp);

	png_destroy_write_struct(&png, &info);
}

void SaveRawPixelData(wchar_t *FileName) {
	std::ofstream File(FileName, std::ofstream::out | std::ofstream::binary);

	int32_t width = FContext.pixelManager.Width, height = FContext.pixelManager.Height;

	size_t X = FContext.pixelManager.Passes[FContext.pixelManager.PassCount - 1].X;
	size_t Y = FContext.pixelManager.Passes[FContext.pixelManager.PassCount - 1].Y;
	size_t TextureWidth = FContext.pixelManager.TextureWidth;

	for (size_t y = Y; y < Y + height; y++) {
		File.write((const char *)&FContext.pixelManager.Data[y * TextureWidth + X], width * sizeof(float));
	}

	File.close();
}