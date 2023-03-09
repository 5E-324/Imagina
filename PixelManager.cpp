#include "Includes.h"

RasterizingInterface::~RasterizingInterface() {}
GroupedRasterizingInterface::~GroupedRasterizingInterface() {}

DefaultGroupedRasterizerInterface::~DefaultGroupedRasterizerInterface() {
	delete[]ValidList;
}

size_t DefaultGroupedRasterizerInterface::GetCoordinate(HRReal *Real, HRReal *Imaginary) {
	ValidCount = 0;
	for (size_t i = 0; i < GroupSize; i++) {
		if (Interfaces[i]->GetCoordinate(Real[ValidCount], Imaginary[ValidCount])) {
			ValidList[ValidCount] = i;
			ValidCount++;
		}
	}
	return ValidCount;
}

void DefaultGroupedRasterizerInterface::WriteResults(SRReal *Values) {
	for (size_t i = 0; i < ValidCount; i++) {
		Interfaces[ValidList[i]]->WriteResults(Values[i]);
	}
}

PixelManager::~PixelManager() {}
void PixelManager::SetResolution(int32_t /*Width*/, int32_t /*Height*/) {

}

GroupedRasterizingInterface &PixelManager::GetGroupedInterface(size_t GroupSize) {
	RasterizingInterface **Interfaces = new RasterizingInterface * [GroupSize];
	for (size_t i = 0; i < GroupSize; i++) {
		Interfaces[i] = &GetInterface();
	}
	return *new DefaultGroupedRasterizerInterface(GroupSize, Interfaces);
}

void PixelManager::FreeInterface(RasterizingInterface &Interface) {
	delete &Interface;
}

void PixelManager::FreeInterface(GroupedRasterizingInterface &Interface) {
	DefaultGroupedRasterizerInterface &I = static_cast<DefaultGroupedRasterizerInterface &>(Interface);
	for (size_t i = 0; i < I.GroupSize; i++) {
		FreeInterface(*I.Interfaces[i]);
	}
	delete &I;
}

void PixelManager::Abort() {}

template<typename T> T RoundUpToEven(T x) {
	return (x + T(1)) & ~T(1);
}

template<typename T> T RoundUpToMultipleOfFour(T x) {
	return (x + T(3)) & ~T(3);
}

void StandardPixelManager::Init() {
	AspectRatio = double(Width) / double(Height);
	OriginOffsetY = -CurrentLocation.HalfH + PixelSize * 0.5;
	OriginOffsetX = -CurrentLocation.HalfH * HRReal(Width) / HRReal(Height) + PixelSize * 0.5;
	Origin.Y = CurrentLocation.Y + OriginOffsetY;
	Origin.X = CurrentLocation.X + OriginOffsetX;

	uint32_t PassWidth = RoundUpToMultipleOfFour(Width); // TODO: Round according to pixel group size
	uint32_t PassHeight = RoundUpToEven(Height);
	int32_t PixelSizeExpX = 0, PixelSizeExpY = 0;

	for (size_t i = PassCount; i > 1;) {
		i--;
		Passes[i].Type = PassType::XDoublingPos;
		Passes[i].Width = PassWidth;
		Passes[i].Height = PassHeight;
		Passes[i].PixelSizeExpX = PixelSizeExpX;
		Passes[i].PixelSizeExpY = PixelSizeExpY;
		Passes[i].PixelGroupSize = std::max(1u, PassWidth / 16);
		Passes[i].CoordOffsetX = 0;
		Passes[i].CoordOffsetY = 0;
		PassWidth = RoundUpToEven(PassWidth / 2 + 1);
		PassHeight = RoundUpToMultipleOfFour(PassHeight);
		PixelSizeExpX++;

		i--;
		Passes[i].Type = PassType::YDoublingPos;
		Passes[i].Width = PassWidth;
		Passes[i].Height = PassHeight;
		Passes[i].PixelSizeExpX = PixelSizeExpX;
		Passes[i].PixelSizeExpY = PixelSizeExpY;
		Passes[i].PixelGroupSize = std::max(1u, PassWidth / 16);
		Passes[i].CoordOffsetX = 0;
		Passes[i].CoordOffsetY = 0;
		PassWidth = RoundUpToMultipleOfFour(PassWidth);
		PassHeight = RoundUpToEven(PassHeight / 2 + 1);
		PixelSizeExpY++;
	}
	Passes[0].Type = PassType::Full;
	Passes[0].Width = PassWidth;
	Passes[0].Height = PassHeight;
	Passes[0].PixelSizeExpX = PixelSizeExpX;
	Passes[0].PixelSizeExpY = PixelSizeExpY;
	Passes[0].PixelGroupSize = std::max(1u, PassWidth / 16);
	Passes[0].CoordOffsetX = 0;
	Passes[0].CoordOffsetY = 0;

	PixelCount = Passes[0].Width * Passes[0].Height;
	PixelGroupCount = PixelCount / PixelGroupSize;
	TextureWidth = Passes[0].Width;
	TextureHeight = Passes[0].Height;

	Passes[0].X = 0;
	Passes[0].Y = 0;
	Passes[0].Begin = 0;
	Passes[0].End = PixelGroupCount; //PixelCount;

	for (size_t i = 1; i < PassCount;) {
		Passes[i].X = TextureWidth; // YDoubling
		Passes[i].Y = 0;
		Passes[i].Begin = PixelGroupCount; //PixelCount;
		PixelCount += Passes[i].Width * Passes[i].Height / 2; // Half of the pixels are new
		PixelGroupCount = PixelCount / PixelGroupSize;
		Passes[i].End = PixelGroupCount; //PixelCount;

		TextureWidth += Passes[i].Width;
		TextureHeight = std::max(TextureHeight, Passes[i].Height);
		i++;

		Passes[i].X = 0; // XDoubling
		Passes[i].Y = TextureHeight;
		Passes[i].Begin = PixelGroupCount; //PixelCount;
		PixelCount += Passes[i].Width * Passes[i].Height / 2;
		PixelGroupCount = PixelCount / PixelGroupSize;
		Passes[i].End = PixelGroupCount; //PixelCount;

		TextureWidth = std::max(TextureWidth, Passes[i].Width);
		TextureHeight += Passes[i].Height;
		i++;
	}

	Coordinates = new Coordinate[PixelGroupCount];
	Data = new float[TextureWidth * TextureHeight];

	size_t CoordIndex = 0;
	size_t Begin = CoordIndex;
	for (int32_t Y = 0; Y < int32_t(Passes[0].Height); Y += PixelGroupHeight) {
		for (int32_t X = 0; X < int32_t(Passes[0].Width); X += PixelGroupWidth) {
			Coordinates[CoordIndex] = { X, Y };
			CoordIndex++;
		}
	}
	size_t End = CoordIndex;

	std::sort(Coordinates + Begin, Coordinates + End,
		[Width = Passes[0].Width, Height = Passes[0].Height](Coordinate a, Coordinate b) { // TODO: Take pixel group size into account
		a.X = a.X * 2 + 1 - Width;
		a.Y = a.Y * 2 + 1 - Height;
		b.X = b.X * 2 + 1 - Width;
		b.Y = b.Y * 2 + 1 - Height;
		return a.X * a.X + a.Y * a.Y < b.X *b.X + b.Y * b.Y;
	});

	for (size_t i = 1; i < PassCount;) {
		Begin = CoordIndex;
		for (int32_t Y = 0; Y < int32_t(Passes[i].Height); Y += 2 * PixelGroupHeight) { //YDoubling
			for (int32_t X = 0; X < int32_t(Passes[i].Width); X += PixelGroupWidth) {
				Coordinates[CoordIndex] = { X, Y };
				CoordIndex++;
			}
		}
		End = CoordIndex;
		std::sort(Coordinates + Begin, Coordinates + End,
			[Width2 = Passes[i].Width * 2, Height = Passes[i].Height](Coordinate a, Coordinate b) {
			a.X = a.X * 4 + 2 - Width2;
			a.Y = a.Y * 2 + 1 - Height;
			b.X = b.X * 4 + 2 - Width2;
			b.Y = b.Y * 2 + 1 - Height;
			return a.X * a.X + a.Y * a.Y < b.X *b.X + b.Y * b.Y;
		});
		i++;

		Begin = CoordIndex;
		for (int32_t Y = 0; Y < int32_t(Passes[i].Height); Y += PixelGroupHeight) { //XDoubling
			for (int32_t X = 0; X < int32_t(Passes[i].Width); X += 2 * PixelGroupWidth) {
				Coordinates[CoordIndex] = { X, Y };
				CoordIndex++;
			}
		}
		End = CoordIndex;
		std::sort(Coordinates + Begin, Coordinates + End,
			[Width = Passes[i].Width, Height = Passes[i].Height](Coordinate a, Coordinate b) {
			a.X = a.X * 2 + 1 - Width;
			a.Y = a.Y * 2 + 1 - Height;
			b.X = b.X * 2 + 1 - Width;
			b.Y = b.Y * 2 + 1 - Height;
			return a.X * a.X + a.Y * a.Y < b.X *b.X + b.Y * b.Y;
		});
		i++;
	}

	memset(Data, 0xFF, sizeof(float) * TextureWidth * TextureHeight);
}

StandardPixelManager::StandardPixelManager() {
	Init();
}

StandardPixelManager::~StandardPixelManager() {
	if (Data) delete[] Data;
	if (Coordinates) delete[] Coordinates;
	if (TextureID) glDeleteTextures(1, &TextureID);
	if (PrevTextureID) glDeleteTextures(1, &PrevTextureID);
}

void StandardPixelManager::ChangeMaxit(uint64_t OldMaxIt) {
	ReusedPass = std::min(std::max(CurrentPass, ReusedPass), PassCount - 1);
	for (size_t i = 0; i < size_t(TextureWidth) * size_t(TextureHeight); i++) {
		if (Data[i] >= float(OldMaxIt)) {
			Data[i] = std::numeric_limits<float>::quiet_NaN();
		}
	}
	OverwriteReused = false;
}

void StandardPixelManager::Clear() {
	ReusedPass = 0;
	memset(Data, 0xFF, sizeof(float) * TextureWidth * TextureHeight);
}

void StandardPixelManager::SetResolution(int32_t width, int32_t height) {
	if (!Completed()) throw "";
	Width = width;
	Height = height;

	if (Data) delete[] Data;
	if (Coordinates) delete[] Coordinates;
	if (TextureID) glDeleteTextures(1, &TextureID);
	if (PrevTextureID) glDeleteTextures(1, &PrevTextureID);
	TextureID = 0;
	PrevTextureID = 0;
	NumPrevTDs = 0;
	ReusedPass = 0;
	PixelReused = false;
	PixelSize = CurrentLocation.HalfH * 2.0_hr / HRReal(Height);

	Init();
}

void StandardPixelManager::UpdateRelativeCoordinate(HRReal OffsetX, HRReal OffsetY) {
	TexRect.X += OffsetX;
	TexRect.Y += OffsetY;
	CurrentLocation.X += OffsetX;
	CurrentLocation.Y += OffsetY;
	Origin.X += OffsetX;
	Origin.Y += OffsetY;

	for (size_t i = 0; i < NumPrevTDs; i++) {
		PrevTDs[i].FractalRect.X += OffsetX;
		PrevTDs[i].FractalRect.Y += OffsetY;
	}
}

void StandardPixelManager::RemovePrevTextures() {
	NumPrevTDs = 0;
}

void MoveImage(float *Data, int32_t RowLength, int32_t Width, int32_t Height, int32_t DiffX, int32_t DiffY) {
	auto CopyRow = [&](int32_t DstRow, int32_t SrcRow) {
		if (DiffX >= 0) {
			memcpy(Data + DstRow * RowLength + DiffX, Data + SrcRow * RowLength, sizeof(float) * (Width - DiffX));
			memset(Data + DstRow * RowLength, 0xFF, sizeof(float) * DiffX);
		} else {
			memcpy(Data + DstRow * RowLength, Data + SrcRow * RowLength + (-DiffX), sizeof(float) * (Width - (-DiffX)));
			memset(Data + DstRow * RowLength + Width - (-DiffX), 0xFF, sizeof(float) * (-DiffX));
		}
	};

	if (DiffY > 0) {
		int32_t y = Height;
		while (y > DiffY) {
			y--;
			CopyRow(y, y - DiffY);
		}
		while (y > 0) {
			y--;
			memset(Data + y * RowLength, 0xFF, sizeof(float) * Width);
		}
	} else if (DiffY < 0) {
		int32_t y = 0;
		while (y < Height - (-DiffY)) {
			CopyRow(y, y - DiffY);
			y++;
		}
		while (y < Height) {
			memset(Data + y * RowLength, 0xFF, sizeof(float) * Width);
			y++;
		}
	} else {
		if (DiffX >= 0) {
			for (int32_t y = 0; y < Height; y++) {
				memmove(Data + y * RowLength + DiffX, Data + y * RowLength, sizeof(float) * (Width - DiffX));
				memset(Data + y * RowLength, 0xFF, sizeof(float) * DiffX);
			}
		} else {
			for (int32_t y = 0; y < Height; y++) {
				memmove(Data + y * RowLength, Data + y * RowLength + (-DiffX), sizeof(float) * (Width - (-DiffX)));
				memset(Data + y * RowLength + Width - (-DiffX), 0xFF, sizeof(float) * (-DiffX));
			}
		}
	}
}

void StandardPixelManager::InvalidatePass(size_t Pass) {
	float *data = Data + Passes[Pass].Y * TextureWidth + Passes[Pass].X;
	size_t Width = Passes[Pass].Width;
	size_t Height = Passes[Pass].Height;
	for (size_t i = 0; i < Height; i++) {
		memset(data + i * TextureWidth, 0xFF, Width * sizeof(float));
	}
}

void StandardPixelManager::AlignPass(size_t Pass, int32_t &DiffX, int32_t &DiffY) {
	int32_t PixelSizeX = 1 << Passes[Pass].PixelSizeExpX;
	int32_t PixelSizeY = 1 << Passes[Pass].PixelSizeExpY;
	DiffX += PixelSizeX - 1;
	DiffY += PixelSizeY - 1;
	Passes[Pass].CoordOffsetX = (DiffX & (PixelSizeX - 1)) - (PixelSizeX - 1);
	Passes[Pass].CoordOffsetY = (DiffY & (PixelSizeY - 1)) - (PixelSizeY - 1);
	DiffX >>= Passes[Pass].PixelSizeExpX;
	DiffY >>= Passes[Pass].PixelSizeExpY;
}
void StandardPixelManager::MovePass(size_t Pass, int32_t DiffX, int32_t DiffY) {
	int32_t DiffX2 = DiffX + Passes[Pass].CoordOffsetX;
	int32_t DiffY2 = DiffY + Passes[Pass].CoordOffsetY;
	int32_t Width = Passes[Pass].Width;
	int32_t Height = Passes[Pass].Height;
	int32_t X = Passes[Pass].X;
	int32_t Y = Passes[Pass].Y;
	AlignPass(Pass, DiffX2, DiffY2);
	MoveImage(Data + Y * TextureWidth + X, TextureWidth, Width, Height, DiffX2, DiffY2);
}

void CopyImage3(float *Dst, int32_t DstWidth, int32_t DstHeight, int32_t DstX, int32_t DstY,
	float *Src, int32_t SrcWidth, int32_t SrcHeight, int32_t SrcX, int32_t SrcY,
	int32_t Width, int32_t Height) {
	int32_t NegXAdjustment = -std::min({ DstX, SrcX, 0 });
	int32_t NegYAdjustment = -std::min({ DstY, SrcY, 0 });
	int32_t PosXAdjustment = std::max({ SrcX + Width - SrcWidth, DstX + Width - DstWidth, 0 });
	int32_t PosYAdjustment = std::max({ SrcY + Height - SrcHeight, DstY + Height - DstHeight, 0 });

	Width -= NegXAdjustment + PosXAdjustment;
	SrcX += NegXAdjustment;
	DstX += NegXAdjustment;

	Height -= NegYAdjustment + PosYAdjustment;
	SrcY += NegYAdjustment;
	DstY += NegYAdjustment;

	for (int32_t y = 0; y < Height; y++) {
		memcpy(Dst + uintptr_t(DstY + y) * DstWidth + DstX, Src + uintptr_t(SrcY + y) * SrcWidth + SrcX, sizeof(float) * Width);
	}
}

void StandardPixelManager::ShiftDownPass(size_t DstPass, int32_t DiffX, int32_t DiffY) {
	size_t SrcPass = DstPass + 2;
	int32_t DiffX2 = DiffX + Passes[SrcPass].CoordOffsetX;
	int32_t DiffY2 = DiffY + Passes[SrcPass].CoordOffsetY;
	int32_t Width = Passes[DstPass].Width;
	int32_t Height = Passes[DstPass].Height;
	int32_t SrcX = Passes[SrcPass].X;
	int32_t SrcY = Passes[SrcPass].Y;
	int32_t DstX = Passes[DstPass].X;
	int32_t DstY = Passes[DstPass].Y;
	AlignPass(SrcPass, DiffX2, DiffY2);
	Passes[DstPass].CoordOffsetX = Passes[SrcPass].CoordOffsetX * 2;
	Passes[DstPass].CoordOffsetY = Passes[SrcPass].CoordOffsetY * 2;
	CopyImage3(Data, TextureWidth, TextureHeight, DstX, DstY, // FIXME: Clamp to edges of each pass, not the whole texture
		Data, TextureWidth, TextureHeight, SrcX - DiffX2, SrcY - DiffY2,
		Width, Height);
}

// TODO: Refactorize
void CopyImage3(float *Dst, int32_t DstWidth, int32_t DstHeight, int32_t DstX, int32_t DstY,
	float *Src, int32_t SrcX, int32_t SrcY,
	int32_t Width, int32_t Height, int32_t RowLength) {

	int32_t NegXAdjustment = -std::min({ DstX, SrcX, 0 });
	int32_t NegYAdjustment = -std::min({ DstY, SrcY, 0 });
	int32_t PosXAdjustment = std::max(DstX + Width - DstWidth, 0);
	int32_t PosYAdjustment = std::max(DstY + Height - DstHeight, 0);

	Width -= NegXAdjustment + PosXAdjustment;
	SrcX += NegXAdjustment;
	DstX += NegXAdjustment;

	Height -= NegYAdjustment + PosYAdjustment;
	SrcY += NegYAdjustment;
	DstY += NegYAdjustment;

	for (int32_t y = 0; y < DstY; y++) {
		memset(Dst + uintptr_t(y) * RowLength, 0xFF, sizeof(float) * DstWidth);
	}
	for (int32_t y = 0; y < Height; y++) {
		memset(Dst + uintptr_t(DstY + y) * RowLength, 0xFF, sizeof(float) * DstX);
		memcpy(Dst + uintptr_t(DstY + y) * RowLength + DstX, Src + uintptr_t(SrcY + y) * RowLength + SrcX, sizeof(float) * Width);
		memset(Dst + uintptr_t(DstY + y) * RowLength + DstX + Width, 0xFF, sizeof(float) * (DstWidth - (DstX + Width)));
	}
	for (int32_t y = DstY + Height; y < DstHeight; y++) {
		memset(Dst + uintptr_t(y) * RowLength, 0xFF, sizeof(float) * DstWidth);
	}
}
void StandardPixelManager::ShiftUpPass(size_t DstPass, int32_t DiffX, int32_t DiffY) {
	size_t SrcPass = DstPass - 2;
	Passes[DstPass].CoordOffsetX = Passes[SrcPass].CoordOffsetX >> 1;
	Passes[DstPass].CoordOffsetY = Passes[SrcPass].CoordOffsetY >> 1;
	int32_t DiffX2 = DiffX + Passes[DstPass].CoordOffsetX;
	int32_t DiffY2 = DiffY + Passes[DstPass].CoordOffsetY;
	int32_t Width = Passes[DstPass].Width;
	int32_t Height = Passes[DstPass].Height;
	int32_t SrcWidth = Passes[SrcPass].Width;
	int32_t SrcHeight = Passes[SrcPass].Height;
	int32_t SrcX = Passes[SrcPass].X;
	int32_t SrcY = Passes[SrcPass].Y;
	int32_t DstX = Passes[DstPass].X;
	int32_t DstY = Passes[DstPass].Y;
	AlignPass(DstPass, DiffX2, DiffY2);

	CopyImage3(Data + DstY * TextureWidth + DstX, Width, Height, DiffX2, DiffY2,
		Data + SrcY * TextureWidth + SrcX, 0, 0,
		SrcWidth, SrcHeight, TextureWidth);
}

void StandardPixelManager::MoveImageStack(int32_t DiffX, int32_t DiffY) {
	size_t Pass = PassCount - 1;
	MovePass(Pass, DiffX, DiffY);

	while (Pass > 0) {
		Pass--;
		MovePass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetX > -(1 << (Passes[Pass].PixelSizeExpX - 1))) {
			Passes[Pass + 1].Type = PassType::XDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::XDoublingNeg;
		}

		Pass--;
		MovePass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetY > -(1 << (Passes[Pass].PixelSizeExpY - 1))) {
			Passes[Pass + 1].Type = PassType::YDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::YDoublingNeg;
		}
	}
}

void StandardPixelManager::ShiftDownImageStack(int32_t DiffX, int32_t DiffY) {
	size_t Pass = 0;
	while (Pass < PassCount - 3) {
		ShiftDownPass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetY > -(1 << (Passes[Pass].PixelSizeExpY - 1))) {
			Passes[Pass + 1].Type = PassType::YDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::YDoublingNeg;
		}
		Pass++;

		ShiftDownPass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetX > -(1 << (Passes[Pass].PixelSizeExpX - 1))) {
			Passes[Pass + 1].Type = PassType::XDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::XDoublingNeg;
		}
		Pass++;
	}
	ShiftDownPass(PassCount - 3, DiffX, DiffY);

	InvalidatePass(PassCount - 2);
	Passes[PassCount - 2].Type = PassType::YDoublingPos;
	Passes[PassCount - 2].CoordOffsetX = 0;
	Passes[PassCount - 2].CoordOffsetY = 0;
	InvalidatePass(PassCount - 1);
	Passes[PassCount - 1].Type = PassType::XDoublingPos;
	Passes[PassCount - 1].CoordOffsetX = 0;
	Passes[PassCount - 1].CoordOffsetY = 0;
}

void StandardPixelManager::ShiftUpImageStack(int32_t DiffX, int32_t DiffY) {
	size_t Pass = PassCount - 1;
	ShiftUpPass(Pass, DiffX, DiffY);

	while (Pass > 2) {
		Pass--;
		ShiftUpPass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetX > -(1 << (Passes[Pass].PixelSizeExpX - 1))) {
			Passes[Pass + 1].Type = PassType::XDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::XDoublingNeg;
		}

		Pass--;
		ShiftUpPass(Pass, DiffX, DiffY);
		if (Passes[Pass].CoordOffsetY > -(1 << (Passes[Pass].PixelSizeExpY - 1))) {
			Passes[Pass + 1].Type = PassType::YDoublingPos;
		} else {
			Passes[Pass + 1].Type = PassType::YDoublingNeg;
		}
	}

	Passes[2].Type = PassType::XDoublingPos; // TODO: Reuse from higher passes
	InvalidatePass(1);
	Passes[1].Type = PassType::YDoublingPos;
	Passes[1].CoordOffsetX = Passes[2].CoordOffsetX;
	Passes[1].CoordOffsetY = Passes[2].CoordOffsetY;
	InvalidatePass(0);
	Passes[0].CoordOffsetX = Passes[2].CoordOffsetX;
	Passes[0].CoordOffsetY = Passes[2].CoordOffsetY;
}

void StandardPixelManager::SetLocation(RelLocation &Location) {
	if (CurrentLocation.HalfH != Location.HalfH) {
		NumPrevTDs = 0;
		size_t NumObtained;
		GetTextures(PrevTDs, 16, NumObtained);
		NumPrevTDs = NumObtained;
		std::swap(TextureID, PrevTextureID);
	} else {
		NumPrevTDs = 0;
	}

	HRReal NewPixelSize = Location.HalfH * 2.0_hr / HRReal(Height);

	RelLocation NewOrigin;
	HRReal HalfW = Location.HalfH * HRReal(Width) / HRReal(Height);
	HRReal HalfH = Location.HalfH;
	HRReal NewOriginOffsetX = -HalfW + NewPixelSize * 0.5;
	HRReal NewOriginOffsetY = -HalfH + NewPixelSize * 0.5;
	NewOrigin.X = Location.X + NewOriginOffsetX;
	NewOrigin.Y = Location.Y + NewOriginOffsetY;

	goto ReuseFailed;
	if (Location.HalfH == 0.5 * CurrentLocation.HalfH && PixelValid) {
		int32_t DiffX = round(SRReal((Origin.X - NewOrigin.X) / PixelSize));
		int32_t DiffY = round(SRReal((Origin.Y - NewOrigin.Y) / PixelSize));
		if (DiffX >= -(Width / 2) && DiffX <= 0 && DiffY >= -(Height / 2) && DiffY <= 0) { // FIXME: Condition
			ShiftDownImageStack(DiffX, DiffY);

			ReusedPass = std::min(std::max(CurrentPass, ReusedPass) - 2, PassCount - 1);
			NewOrigin.X = Origin.X - DiffX * PixelSize;
			NewOrigin.Y = Origin.Y - DiffY * PixelSize;
			PixelReused = true;
			OverwriteReused = Global::HighQuality;
		} else {
			goto ReuseFailed;
		}
	} else if (Location.HalfH == CurrentLocation.HalfH) {
		int32_t DiffX = round(SRReal((Origin.X - NewOrigin.X) / PixelSize));
		int32_t DiffY = round(SRReal((Origin.Y - NewOrigin.Y) / PixelSize));
		if (DiffX > -Width && DiffX < Width && DiffY > -Height && DiffY < Height) { // FIXME: Condition
			if (DiffX || DiffY) {
				MoveImageStack(DiffX, DiffY);
				//ReusedPass = PassCount - 1;
				ReusedPass = std::min(std::max(CurrentPass, ReusedPass), PassCount - 1);
				PixelReused = true;
				OverwriteReused = false;
			}
			NewOrigin.X = Origin.X - DiffX * PixelSize;
			NewOrigin.Y = Origin.Y - DiffY * PixelSize;
		} else {
			goto ReuseFailed;
		}
	} else if (Location.HalfH == 2.0 * CurrentLocation.HalfH) {
		if (Passes[0].CoordOffsetX & 1) Origin.X += PixelSize; // The LSB of CoordOffset will be discarded by rounding, so add it to origin first
		if (Passes[0].CoordOffsetY & 1) Origin.Y += PixelSize;
		int32_t DiffX = round(SRReal((Origin.X - NewOrigin.X) / NewPixelSize));
		int32_t DiffY = round(SRReal((Origin.Y - NewOrigin.Y) / NewPixelSize));
		if (DiffX >= 0 && DiffX <= Width / 2 && DiffY >= 0 && DiffY <= Height / 2) { // FIXME: Condition
			ShiftUpImageStack(DiffX, DiffY);

			ReusedPass = std::min(std::max(CurrentPass, ReusedPass) + 2, PassCount - 1);
			NewOrigin.X = Origin.X - DiffX * NewPixelSize;
			NewOrigin.Y = Origin.Y - DiffY * NewPixelSize;
			PixelReused = true;
			OverwriteReused = Global::HighQuality;
		} else {
			goto ReuseFailed;
		}
	} else {
	ReuseFailed:
		PixelReused = false;
		ReusedPass = 0;
		memset(Data, 0xFF, sizeof(float) * TextureWidth * TextureHeight);
	}
	NeedUpdateTexture = true;

	PixelSize = NewPixelSize;
	Origin = NewOrigin;
	OriginOffsetX = NewOriginOffsetX;
	OriginOffsetY = NewOriginOffsetY;
	Location.X = Origin.X + HalfW - PixelSize * 0.5;
	Location.Y = Origin.Y + HalfH - PixelSize * 0.5;
	CurrentLocation = Location;
}

void StandardPixelManager::Begin() {
	i = 0;
	PixelValid = true;
	idle = false;
	abort = false;
	CurrentPass = 0;
	NeedUpdateTexture = true;
}

RasterizingInterface &StandardPixelManager::GetInterface() {
	StandardRI *Interface = new StandardRI(*this);
	InterfaceListMutex.lock();
	if (FirstInterface) {
		FirstInterface->PrevInterface = Interface;
		Interface->NextInterface = FirstInterface;
	}
	FirstInterface = Interface;
	InterfaceListMutex.unlock();
	return *Interface;
}

GroupedRasterizingInterface &StandardPixelManager::GetGroupedInterface(size_t GroupSize) {
	if (GroupSize != 4) throw "Not implemented";
	StandardGRI *Interface = new StandardGRI(*this);
	InterfaceListMutex.lock();
	if (FirstVInterface) {
		FirstVInterface->PrevInterface = Interface;
		Interface->NextInterface = FirstVInterface;
	}
	FirstVInterface = Interface;
	InterfaceListMutex.unlock();
	return *Interface;
}

void StandardPixelManager::FreeInterface(RasterizingInterface &Interface) {
	StandardRI *RI = static_cast<StandardRI *>(&Interface);
	InterfaceListMutex.lock();
	RI->TryCopyData();
	if (FirstCopyPendingInterface) { // TODO: Improve copying
		bool Copied = false;
		do {
			Copied = false;
			for (StandardRI *CPRI = FirstCopyPendingInterface; CPRI; ) {
				Copied = Copied || CPRI->TryCopyData();
				if (CPRI->CopyQueue.size() == 0) {
					StandardRI *NextInterface = CPRI->NextInterface;
					if (CPRI->NextInterface) {
						CPRI->NextInterface->PrevInterface = CPRI->PrevInterface;
					}
					if (CPRI->PrevInterface) {
						CPRI->PrevInterface->NextInterface = CPRI->NextInterface;
					} else {
						FirstCopyPendingInterface = CPRI->NextInterface;
					}
					delete CPRI;
					CPRI = NextInterface;
				} else {
					CPRI = CPRI->NextInterface;
				}
			}
		} while (FirstCopyPendingInterface && Copied);
	}

	if (RI->NextInterface) {
		RI->NextInterface->PrevInterface = RI->PrevInterface;
	}
	if (RI->PrevInterface) {
		RI->PrevInterface->NextInterface = RI->NextInterface;
	} else {
		FirstInterface = RI->NextInterface;
	}
	if (RI->CopyQueue.size() != 0) {
		if (FirstCopyPendingInterface) {
			FirstCopyPendingInterface->PrevInterface = RI;
			RI->NextInterface = FirstCopyPendingInterface;
		} else {
			RI->NextInterface = nullptr;
		}
		RI->PrevInterface = nullptr;
		FirstCopyPendingInterface = RI;
	} else {
		delete &Interface;
	}
	if (!FirstInterface && FirstCopyPendingInterface) {
		bool Copied = false;
		do {
			Copied = false;
			for (StandardRI *CPRI = FirstCopyPendingInterface; CPRI; ) {
				Copied = Copied || CPRI->TryCopyData();
				if (CPRI->CopyQueue.size() == 0) {
					StandardRI *NextInterface = CPRI->NextInterface;
					if (CPRI->NextInterface) {
						CPRI->NextInterface->PrevInterface = CPRI->PrevInterface;
					}
					if (CPRI->PrevInterface) {
						CPRI->PrevInterface->NextInterface = CPRI->NextInterface;
					} else {
						FirstCopyPendingInterface = CPRI->NextInterface;
					}
					delete CPRI;
					CPRI = NextInterface;
				} else {
					CPRI = CPRI->NextInterface;
				}
			}
		} while (FirstCopyPendingInterface && Copied);
	}
	InterfaceListMutex.unlock();
}

void StandardPixelManager::FreeInterface(GroupedRasterizingInterface &Interface) {
	StandardGRI *RI = static_cast<StandardGRI *>(&Interface);
	InterfaceListMutex.lock();
	RI->TryCopyData();
	if (FirstCopyPendingVInterface) { // TODO: Improve copying
		bool Copied = false;
		do {
			Copied = false;
			for (StandardGRI *CPRI = FirstCopyPendingVInterface; CPRI; ) {
				Copied = Copied || CPRI->TryCopyData();
				if (CPRI->CopyQueue.size() == 0) {
					StandardGRI *NextInterface = CPRI->NextInterface;
					if (CPRI->NextInterface) {
						CPRI->NextInterface->PrevInterface = CPRI->PrevInterface;
					}
					if (CPRI->PrevInterface) {
						CPRI->PrevInterface->NextInterface = CPRI->NextInterface;
					} else {
						FirstCopyPendingVInterface = CPRI->NextInterface;
					}
					delete CPRI;
					CPRI = NextInterface;
				} else {
					CPRI = CPRI->NextInterface;
				}
			}
		} while (FirstCopyPendingVInterface && Copied);
	}

	if (RI->NextInterface) {
		RI->NextInterface->PrevInterface = RI->PrevInterface;
	}
	if (RI->PrevInterface) {
		RI->PrevInterface->NextInterface = RI->NextInterface;
	} else {
		FirstVInterface = RI->NextInterface;
	}
	if (RI->CopyQueue.size() != 0) {
		if (FirstCopyPendingVInterface) {
			FirstCopyPendingVInterface->PrevInterface = RI;
			RI->NextInterface = FirstCopyPendingVInterface;
		} else {
			RI->NextInterface = nullptr;
		}
		RI->PrevInterface = nullptr;
		FirstCopyPendingVInterface = RI;
	} else {
		delete &Interface;
	}
	if (!FirstVInterface && FirstCopyPendingVInterface) {
		bool Copied = false;
		do {
			Copied = false;
			for (StandardGRI *CPRI = FirstCopyPendingVInterface; CPRI;) {
				Copied = Copied || CPRI->TryCopyData();
				if (CPRI->CopyQueue.size() == 0) {
					StandardGRI *NextInterface = CPRI->NextInterface;
					if (CPRI->NextInterface) {
						CPRI->NextInterface->PrevInterface = CPRI->PrevInterface;
					}
					if (CPRI->PrevInterface) {
						CPRI->PrevInterface->NextInterface = CPRI->NextInterface;
					} else {
						FirstCopyPendingVInterface = CPRI->NextInterface;
					}
					delete CPRI;
					CPRI = NextInterface;
				} else {
					CPRI = CPRI->NextInterface;
				}
			}
		} while (FirstCopyPendingVInterface && Copied);
	}
	InterfaceListMutex.unlock();
}

void StandardPixelManager::GetTextures(TextureDescription *TD, size_t NumDesired, size_t &NumObtained) {
	if (!NumDesired)
		return;

	size_t Pass = std::min(CurrentPass, PassCount - 1);
	Pass = std::max(Pass, ReusedPass);
	size_t MinPass = Pass;

	if (abort) {
		MinPass = MinPassBeforeAbort;
	} else {
		InterfaceListMutex.lock();
		for (StandardRI *RI = FirstInterface; RI; RI = RI->NextInterface) {
			MinPass = std::min(MinPass, RI->Pass);
		}
		for (StandardGRI *RI = FirstVInterface; RI; RI = RI->NextInterface) {
			MinPass = std::min(MinPass, RI->Pass);
		}
		MinPass = std::max(MinPass, size_t(1)) - 1;
		InterfaceListMutex.unlock();
	}
	if (!TextureID) {
		glGenTextures(1, &TextureID);
		glBindTexture(TARGET_TEXTRUE, TextureID);
		glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_MAG_FILTER, Global::UseBilinearFilter ? GL_LINEAR : GL_NEAREST);
		glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_MIN_FILTER, Global::UseBilinearFilter ? GL_LINEAR : GL_NEAREST);
		glTexImage2D(TARGET_TEXTRUE, 0, GL_R32F, TextureWidth, TextureHeight, 0, GL_RED, GL_FLOAT, Data);
	} else if (NeedUpdateTexture) {
		glBindTexture(TARGET_TEXTRUE, TextureID);
		if (Completed()) NeedUpdateTexture = false;
		glTexSubImage2D(TARGET_TEXTRUE, 0, 0, 0, TextureWidth, TextureHeight, GL_RED, GL_FLOAT, Data);
	}
	glBindTexture(TARGET_TEXTRUE, TextureID);
	glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_MAG_FILTER, Global::UseBilinearFilter ? GL_LINEAR : GL_NEAREST);
	glTexParameteri(TARGET_TEXTRUE, GL_TEXTURE_MIN_FILTER, Global::UseBilinearFilter ? GL_LINEAR : GL_NEAREST);
	size_t i = 0;
	for (; i < NumDesired; i++, Pass--) {
		size_t X = Passes[Pass].X, Y = Passes[Pass].Y;
		size_t Width = Passes[Pass].Width, Height = Passes[Pass].Height;
		size_t ExpX = Passes[Pass].PixelSizeExpX, ExpY = Passes[Pass].PixelSizeExpY;
		TD[i].Texture = TextureID;
		if (Global::UseBilinearFilter) {
			TD[i].TextureRect.X = SRReal(X) / 2048 + (0.5 / 2048.0);
			TD[i].TextureRect.Y = SRReal(Y) / 2048 + (0.5 / 2048.0);
			TD[i].TextureRect.Width = SRReal(Width) / 2048 - (1.0 / 2048.0);
			TD[i].TextureRect.Height = SRReal(Height) / 2048 - (1.0 / 2048.0);
			TD[i].TextureRect.X = SRReal(X) + 0.5;
			TD[i].TextureRect.Y = SRReal(Y) + 0.5;
			TD[i].TextureRect.Width = SRReal(Width - 1);
			TD[i].TextureRect.Height = SRReal(Height - 1);
		} else {
			TD[i].TextureRect.X = SRReal(X);
			TD[i].TextureRect.Y = SRReal(Y);
			TD[i].TextureRect.Width = SRReal(Width);
			TD[i].TextureRect.Height = SRReal(Height);
		}
#ifndef USE_GL_TEXTURE_RECTANGLE
		TD[i].TextureRect.X /= HRReal(TextureWidth);
		TD[i].TextureRect.Y /= HRReal(TextureHeight);
		TD[i].TextureRect.Width /= HRReal(TextureWidth);
		TD[i].TextureRect.Height /= HRReal(TextureHeight);
#endif
		int32_t OffsetX = Passes[Pass].CoordOffsetX;
		int32_t OffsetY = Passes[Pass].CoordOffsetY;
		TD[i].FractalRect.X = Origin.X + HRReal(OffsetX * 2 + ((Width - 1) << ExpX)) * 0.5 * PixelSize;
		TD[i].FractalRect.Y = Origin.Y + HRReal(OffsetY * 2 + ((Height - 1) << ExpY)) * 0.5 * PixelSize;
		if (Global::UseBilinearFilter) {
			TD[i].FractalRect.HalfW = HRReal((Width << ExpX) - 1) * PixelSize * 0.5;
			TD[i].FractalRect.HalfH = HRReal((Height << ExpY) - 1) * PixelSize * 0.5;
		} else {
			TD[i].FractalRect.HalfW = HRReal(Width << ExpX) * PixelSize * 0.5;
			TD[i].FractalRect.HalfH = HRReal(Height << ExpY) * PixelSize * 0.5;
		}
		NumObtained = i + 1;
		if (Pass <= MinPass) {
			i++;
			break;
		}
	}
	for (size_t j = 0; i < NumDesired && j < NumPrevTDs; i++, j++) {
		TD[i] = PrevTDs[j];
		NumObtained = i + 1;
	}

	std::sort(TD, TD + NumObtained, [](const TextureDescription &a, const TextureDescription &b) {
		return (a.FractalRect.HalfH * a.FractalRect.HalfW) / (a.TextureRect.Height * a.TextureRect.Width)
			< (b.FractalRect.HalfH * b.FractalRect.HalfW) / (b.TextureRect.Height * b.TextureRect.Width);
		});
}

bool StandardPixelManager::GetProgress(SRReal &Numerator, SRReal &Denoninator) const {
	Numerator = i;
	Denoninator = PixelGroupCount;
	return true;
}

bool StandardPixelManager::Completed() {
	return (idle || (i >= PixelGroupCount) || abort) && !FirstInterface && !FirstVInterface;
}

void StandardPixelManager::Abort() {
	size_t MinPass = std::min(CurrentPass, PassCount - 1);
	MinPass = std::max(MinPass, ReusedPass);

	InterfaceListMutex.lock();
	for (StandardRI *RI = FirstInterface; RI; RI = RI->NextInterface) {
		MinPass = std::min(MinPass, RI->Pass);
	}
	for (StandardGRI *RI = FirstVInterface; RI; RI = RI->NextInterface) {
		MinPass = std::min(MinPass, RI->Pass);
	}
	InterfaceListMutex.unlock();
	MinPass = std::max(MinPass, size_t(1)) - 1;
	MinPassBeforeAbort = MinPass;
	abort = true;
}

bool StandardRI::TryCopyData() {
	bool Copied = false;
	while (!CopyQueue.empty()) {
		Copy &cp = CopyQueue.front();
		float Data = R.Data[cp.SrcIndex];
		if (Data != Data) break;
		R.Data[cp.DstIndex] = Data;
		CopyQueue.pop();
		Copied = true;
	}
	return Copied;
}

void StandardRI::CopyData() {
	while (!CopyQueue.empty()) {
		Copy &cp = CopyQueue.front();
		float Data = R.Data[cp.SrcIndex];
		while (Data != Data) {
			if (R.abort) break;
			std::this_thread::yield();
			Data = R.Data[cp.SrcIndex];
		}
		R.Data[cp.DstIndex] = Data;
		CopyQueue.pop();
	}
}

StandardRI::~StandardRI() {
}

bool StandardRI::GetCoordinate(HRReal &Real, HRReal &Imaginary) {
	if (R.abort) return false;
start:
	HRReal X = R.Origin.X;
	HRReal Y = R.Origin.Y;
	if (YInGroup >= StandardPixelManager::PixelGroupHeight) {
		YInGroup = 0;
		j++;
		if (j >= End) {
			Pass = R.CurrentPass;
			if (Pass >= R.PassCount) return false;
			uint32_t Increment = R.Passes[Pass].PixelGroupSize;
			j = R.i.fetch_add(Increment);
			if (j >= R.PixelGroupCount) return false;
			if (j >= R.Passes[Pass].End) {
				do {
					Pass++;
					if (Pass >= R.PassCount) {
						R.CurrentPass = R.PassCount;
						return false;
					}
				} while (j >= R.Passes[Pass].End);
				R.CurrentPass = Pass;
			}
			End = std::min(j + Increment, R.PixelGroupCount);
		}
		if (j >= R.Passes[Pass].End) {
			Pass++;
		}
	} else if (j >= R.PixelGroupCount) return false;

	using PassType = StandardPixelManager::PassType;

	PassType passType = R.Passes[Pass].Type;
	StandardPixelManager::Coordinate coord = R.Coordinates[j];

	if (passType == PassType::XDoublingPos || passType == PassType::XDoublingNeg) {
		coord.X += XInGroup * 2;
	} else {
		coord.X += XInGroup;
	}
	if (passType == PassType::YDoublingPos || passType == PassType::YDoublingNeg) {
		coord.Y += YInGroup * 2;
	} else {
		coord.Y += YInGroup;
	}
	XInGroup++;
	if (XInGroup >= StandardPixelManager::PixelGroupWidth) {
		YInGroup++;
		XInGroup = 0;
	}
	Coord = coord;

	if (passType == PassType::XDoublingPos) coord.X += 1;
	if (passType == PassType::YDoublingPos) coord.Y += 1;

	size_t Index = int64_t(R.TextureWidth) * int64_t(coord.Y + R.Passes[Pass].Y) + int64_t(coord.X + R.Passes[Pass].X);
	float Value = R.Data[Index];
	if (!R.OverwriteReused && Value == Value) {
		WriteResults(Value);
		goto start;
	}

	Real = X + ((coord.X << R.Passes[Pass].PixelSizeExpX) + R.Passes[Pass].CoordOffsetX) * R.PixelSize;
	Imaginary = Y + ((coord.Y << R.Passes[Pass].PixelSizeExpY) + R.Passes[Pass].CoordOffsetY) * R.PixelSize;



	return true;
}
void StandardRI::WriteResults(SRReal Value) {
	TryCopyData();

	if (Global::PreModulo && Value < Global::MaxIt) {
		Value = fmod(Value, Global::ItMod);
	}

	StandardPixelManager::PassInfo &CurrentPass = R.Passes[Pass], &PreviousPass = R.Passes[Pass - 1];

	size_t Index = int64_t(R.TextureWidth) * int64_t(Coord.Y + CurrentPass.Y) + int64_t(Coord.X + CurrentPass.X);

	if (R.Passes[Pass].Type == StandardPixelManager::PassType::Full) {
		R.Data[Index] = Value;
	} else if (R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingPos || R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingNeg) {
		size_t PrevPassIndex = int64_t(R.TextureWidth) * int64_t(Coord.Y + PreviousPass.Y) + int64_t((Coord.X >> 1) + PreviousPass.X);
		size_t CopyDst;
		if (R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingPos) {
			R.Data[Index + 1] = Value;
			CopyDst = Index;
		} else {
			R.Data[Index] = Value;
			CopyDst = Index + 1;
			PrevPassIndex += 1;
		}
		float PrevPassData = R.Data[PrevPassIndex];
		if (PrevPassData != PrevPassData) {
			CopyQueue.emplace(CopyDst, PrevPassIndex);
		} else {
			R.Data[CopyDst] = PrevPassData;
		}
	} else if (R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingPos || R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingNeg) {
		size_t PrevPassIndex = int64_t(R.TextureWidth) * int64_t((Coord.Y >> 1) + PreviousPass.Y) + int64_t(Coord.X + PreviousPass.X);
		size_t CopyDst;
		if (R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingPos) {
			R.Data[Index + R.TextureWidth] = Value;
			CopyDst = Index;
		} else {
			R.Data[Index] = Value;
			CopyDst = Index + R.TextureWidth;
			PrevPassIndex += R.TextureWidth;
		}
		float PrevPassData = R.Data[PrevPassIndex];
		if (PrevPassData != PrevPassData) {
			CopyQueue.emplace(CopyDst, PrevPassIndex);
		} else {
			R.Data[CopyDst] = PrevPassData;
		}
	}
}

bool StandardGRI::TryCopyData() {
	bool Copied = false;
	while (!CopyQueue.empty()) {
		Copy &cp = CopyQueue.front();
		float Data = R.Data[cp.SrcIndex];
		if (Data != Data) break;
		R.Data[cp.DstIndex] = Data;
		CopyQueue.pop();
		Copied = true;
	}
	return Copied;
}

void StandardGRI::CopyData() {
	while (!CopyQueue.empty()) {
		Copy &cp = CopyQueue.front();
		float Data = R.Data[cp.SrcIndex];
		while (Data != Data) {
			if (R.abort) break;
			std::this_thread::yield();
			Data = R.Data[cp.SrcIndex];
		}
		R.Data[cp.DstIndex] = Data;
		CopyQueue.pop();
	}
}

StandardGRI::~StandardGRI() {
	TryCopyData();
}

size_t StandardGRI::GetCoordinate(HRReal *Real, HRReal *Imaginary) {
	if (R.abort) return 0;
start:
	//HRReal X = R.Origin.X;
	//HRReal Y = R.Origin.Y;
	HRReal X = R.CurrentLocation.X;
	HRReal Y = R.CurrentLocation.Y;
	if (j >= i) {
		Pass = R.CurrentPass;
		if (Pass >= R.PassCount) return 0;
		uint32_t Increment = R.Passes[Pass].PixelGroupSize;
		j = R.i.fetch_add(Increment);
		if (j >= R.PixelGroupCount) return 0;
		if (j >= R.Passes[Pass].End) {
			do {
				Pass++;
				if (Pass >= R.PassCount) {
					R.CurrentPass = R.PassCount;
					return 0;
				}
			} while (j >= R.Passes[Pass].End);
			R.CurrentPass = Pass;
		}
		i = std::min(j + Increment, R.PixelGroupCount);
	} else if (j >= R.PixelGroupCount) return 0;
	size_t i = j;
	j++;
	if (i >= R.Passes[Pass].End) {
		Pass++;
	}

	using PassType = StandardPixelManager::PassType;

	PassType passType = R.Passes[Pass].Type;
	SRReal Value[4];

	StandardPixelManager::Coordinate GroupCoordinate = R.Coordinates[i];

	for (int k = 0; k < 4; k++) {
		IntPix XInGroup = k & 1;
		IntPix YInGroup = k >> 1;
		StandardPixelManager::Coordinate coord = GroupCoordinate;
		if (passType == PassType::XDoublingPos || passType == PassType::XDoublingNeg) {
			coord.X += XInGroup * 2;
		} else {
			coord.X += XInGroup;
		}
		if (passType == PassType::YDoublingPos || passType == PassType::YDoublingNeg) {
			coord.Y += YInGroup * 2;
		} else {
			coord.Y += YInGroup;
		}
		Coord[k] = coord;
		if (passType == PassType::XDoublingPos) coord.X += 1;
		if (passType == PassType::YDoublingPos) coord.Y += 1;

		size_t Index = int64_t(R.TextureWidth) * int64_t(coord.Y + R.Passes[Pass].Y) + int64_t(coord.X + R.Passes[Pass].X);
		Value[k] = R.Data[Index];

		//Real[k] = X + ((coord.X << R.Passes[Pass].PixelSizeExpX) + R.Passes[Pass].CoordOffsetX) * R.PixelSize;
		//Imaginary[k] = Y + ((coord.Y << R.Passes[Pass].PixelSizeExpY) + R.Passes[Pass].CoordOffsetY) * R.PixelSize;

		HRReal x = R.OriginOffsetX + ((coord.X << R.Passes[Pass].PixelSizeExpX) + R.Passes[Pass].CoordOffsetX) * R.PixelSize;
		HRReal y = R.OriginOffsetY + ((coord.Y << R.Passes[Pass].PixelSizeExpY) + R.Passes[Pass].CoordOffsetY) * R.PixelSize;

		Real[k] = X + Global::InvTransformMatrix[0][0] * x + Global::InvTransformMatrix[1][0] * y;
		Imaginary[k] = Y + Global::InvTransformMatrix[0][1] * x + Global::InvTransformMatrix[1][1] * y;
	}
	if (!R.OverwriteReused && Value[0] == Value[0] && Value[1] == Value[1] && Value[2] == Value[2] && Value[3] == Value[3]) {
		WriteResults(Value);
		goto start;
	}

	return 4;
}

void StandardGRI::WriteResults(SRReal *Values) {
	TryCopyData();
	StandardPixelManager::PassInfo &CurrentPass = R.Passes[Pass], &PreviousPass = R.Passes[Pass - 1];

	size_t Index[4];
	for (size_t i = 0; i < 4; i++) {
		if (Global::PreModulo && Values[i] < Global::MaxIt) {
			Values[i] = fmod(Values[i], Global::ItMod);
		}

		Index[i] = int64_t(R.TextureWidth) * int64_t(Coord[i].Y + CurrentPass.Y) + int64_t(Coord[i].X + CurrentPass.X);
	}

	if (R.Passes[Pass].Type == StandardPixelManager::PassType::Full) {
		for (int i = 0; i < 4; i++) {
			float Value = Values[i];
			R.Data[Index[i]] = Value;
		}
	} else if (R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingPos || R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingNeg) {
		for (int i = 0; i < 4; i++) {
			float Value = Values[i];
			size_t PrevPassIndex = int64_t(R.TextureWidth) * int64_t(Coord[i].Y + PreviousPass.Y) + int64_t((Coord[i].X >> 1) + PreviousPass.X);
			size_t CopyDst;
			if (R.Passes[Pass].Type == StandardPixelManager::PassType::XDoublingPos) {
				R.Data[Index[i] + 1] = Value;
				CopyDst = Index[i];
			} else {
				R.Data[Index[i]] = Value;
				CopyDst = Index[i] + 1;
				PrevPassIndex += 1;
			}
			float PrevPassData = R.Data[PrevPassIndex];
			if (PrevPassData != PrevPassData) {
				CopyQueue.emplace(CopyDst, PrevPassIndex);
			} else {
				R.Data[CopyDst] = PrevPassData;
			}
		}
	} else if (R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingPos || R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingNeg) {
		for (int i = 0; i < 4; i++) {
			float Value = Values[i];
			size_t PrevPassIndex = int64_t(R.TextureWidth) * int64_t((Coord[i].Y >> 1) + PreviousPass.Y) + int64_t(Coord[i].X + PreviousPass.X);
			size_t CopyDst;
			if (R.Passes[Pass].Type == StandardPixelManager::PassType::YDoublingPos) {
				R.Data[Index[i] + R.TextureWidth] = Value;
				CopyDst = Index[i];
			} else {
				R.Data[Index[i]] = Value;
				CopyDst = Index[i] + R.TextureWidth;
				PrevPassIndex += R.TextureWidth;
			}
			float PrevPassData = R.Data[PrevPassIndex];
			if (PrevPassData != PrevPassData) {
				CopyQueue.emplace(CopyDst, PrevPassIndex);
			} else {
				R.Data[CopyDst] = PrevPassData;
			}
		}
	}
}