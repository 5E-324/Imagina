#pragma once

#include <queue>

class RasterizingInterface {
public:
	virtual ~RasterizingInterface();

	virtual bool GetCoordinate(HRReal &Real, HRReal &Imaginary) = 0;
	virtual void WriteResults(SRReal Value) = 0;
};

class GroupedRasterizingInterface {
public:
	virtual ~GroupedRasterizingInterface();

	virtual size_t GetCoordinate(HRReal *Real, HRReal *Imaginary) = 0;
	virtual void WriteResults(SRReal *Values) = 0;
};

class DefaultGroupedRasterizerInterface : public GroupedRasterizingInterface {
	friend class PixelManager;
	size_t GroupSize;
	size_t ValidCount;
	size_t *ValidList;
	RasterizingInterface **Interfaces;
public:
	inline DefaultGroupedRasterizerInterface(size_t GroupSize, RasterizingInterface **Interfaces) : GroupSize(GroupSize), Interfaces(Interfaces) {
		ValidList = new size_t[GroupSize];
	}
	virtual ~DefaultGroupedRasterizerInterface() override;
	virtual size_t GetCoordinate(HRReal *Real, HRReal *Imaginary) override;
	virtual void WriteResults(SRReal *Values) override;
};

class PixelManager {
public:
	virtual ~PixelManager();

	virtual void Clear() = 0;
	virtual void SetResolution(int32_t Width, int32_t Height);
	virtual void UpdateRelativeCoordinate(HRReal OffsetX, HRReal OffsetY) = 0;
	virtual void SetLocation(RelLocation &Location) = 0;
	virtual void Begin() = 0;

	virtual RasterizingInterface &GetInterface() = 0;
	virtual GroupedRasterizingInterface &GetGroupedInterface(size_t GroupSize);
	virtual void FreeInterface(RasterizingInterface &Interface);
	virtual void FreeInterface(GroupedRasterizingInterface &Interface);
	virtual bool Completed() = 0;
	virtual void Abort();

	virtual void GetTextures(TextureDescription *TD, size_t NumDesired, size_t &NumObtained) = 0;
};

class StandardRI;
class StandardGRI;

class StandardPixelManager : public PixelManager, public ProgressTrackable {
public:
	int32_t Width = 2048, Height = 1024;
	int32_t TextureWidth, TextureHeight;
	static constexpr size_t DoublingCount = 8;
	static constexpr size_t PassCount = DoublingCount * 2 + 1;

	static constexpr IntPix PixelGroupWidth = 2;
	static constexpr IntPix PixelGroupHeight = 2;
	static constexpr size_t PixelGroupSize = size_t(PixelGroupWidth) * PixelGroupHeight;

	HRReal AspectRatio;

	GLuint TextureID = 0;
	GLuint PrevTextureID = 0;

	size_t NumPrevTDs = 0;
	TextureDescription PrevTDs[16];

	RelRect TexRect = {};
	RelLocation CurrentLocation = {};
	RelLocation Origin;
	float *Data = nullptr;
	std::atomic_uint i = 0;
	bool PixelReused = false;
	bool OverwriteReused = false;
	bool PixelValid = false;
	enum class PassType {
		Full = 0,
		//XDoubling,
		//YDoubling,

		XDoublingPos = 0b100,
		XDoublingNeg = 0b101,
		YDoublingPos = 0b110,
		YDoublingNeg = 0b111,
	};
	struct PassInfo {
		PassType Type;
		int32_t Width, Height;
		int32_t PixelSizeExpX, PixelSizeExpY;

		uint32_t X, Y;
		size_t Begin, End;
		uint32_t PixelGroupSize;
		int32_t CoordOffsetX, CoordOffsetY;
	};
	struct Coordinate {
		int32_t X, Y;
	};
	HRReal PixelSize;
	Coordinate *Coordinates = nullptr;
	PassInfo Passes[PassCount];
	size_t CurrentPass;
	size_t MinPassBeforeAbort;
	size_t PixelCount, PixelGroupCount;
	bool idle = true;
	bool abort = false;
	bool NeedUpdateTexture = false;

	size_t ReusedPass = 0;

	std::mutex InterfaceListMutex;
	StandardRI *FirstInterface = nullptr;
	StandardRI *FirstCopyPendingInterface = nullptr;
	StandardGRI *FirstVInterface = nullptr;
	StandardGRI *FirstCopyPendingVInterface = nullptr;

private:
	void Init();
	void InvalidatePass(size_t Pass);
	void AlignPass(size_t Pass, int32_t &DiffX, int32_t &DiffY);
	void MovePass(size_t Pass, int32_t DiffX, int32_t DiffY);
	void ShiftDownPass(size_t DstPass, int32_t DiffX, int32_t DiffY);
	void ShiftUpPass(size_t DstPass, int32_t DiffX, int32_t DiffY);
	void MoveImageStack(int32_t DiffX, int32_t DiffY);
	void ShiftDownImageStack(int32_t DiffX, int32_t DiffY);
	void ShiftUpImageStack(int32_t DiffX, int32_t DiffY);

public:
	StandardPixelManager();
	virtual ~StandardPixelManager() override;

	void ChangeMaxit(uint64_t OldMaxIt);
	virtual void Clear() override;
	virtual void SetResolution(int32_t Width, int32_t Height) override;
	virtual void UpdateRelativeCoordinate(HRReal OffsetX, HRReal OffsetY) override;
	void RemovePrevTextures();
	virtual void SetLocation(RelLocation &Location) override;
	virtual void Begin() override;

	virtual RasterizingInterface &GetInterface() override;
	virtual GroupedRasterizingInterface &GetGroupedInterface(size_t GroupSize) override;
	virtual void FreeInterface(RasterizingInterface &Interface) override;
	virtual void FreeInterface(GroupedRasterizingInterface &Interface) override;
	virtual bool Completed() override;
	virtual void Abort() override;

	virtual void GetTextures(TextureDescription *TD, size_t NumDesired, size_t &NumObtained) override;

	virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const override;
};

class StandardRI : public RasterizingInterface {
	friend class StandardPixelManager;
	StandardPixelManager &R;
	size_t Pass = 0;
	size_t End = 0, j = 0;
	StandardPixelManager::Coordinate Coord;
	StandardRI *PrevInterface = nullptr;
	StandardRI *NextInterface = nullptr;

	IntPix XInGroup = 0, YInGroup = StandardPixelManager::PixelGroupHeight;

	struct Copy {
		size_t DstIndex, SrcIndex;
		Copy() = default;
		Copy(size_t dst, size_t src) :DstIndex(dst), SrcIndex(src) {}
	};

	std::queue<Copy> CopyQueue;

	bool TryCopyData();
	void CopyData();

public:
	inline StandardRI(StandardPixelManager &r) : R(r) {}
	~StandardRI() override;
	virtual bool GetCoordinate(HRReal &Real, HRReal &Imaginary) override;
	virtual void WriteResults(SRReal Value) override;
};

class StandardGRI : public GroupedRasterizingInterface {
	friend class StandardPixelManager;
	StandardPixelManager &R;
	//int I;
	size_t Pass = 0;
	size_t i = 0, j = 0;
	StandardPixelManager::Coordinate Coord[4];
	StandardGRI *PrevInterface = nullptr;
	StandardGRI *NextInterface = nullptr;

	struct Copy {
		size_t DstIndex, SrcIndex;
		Copy() = default;
		Copy(size_t dst, size_t src) :DstIndex(dst), SrcIndex(src) {}
	};

	std::queue<Copy> CopyQueue;

	bool TryCopyData();
	void CopyData();

public:
	inline StandardGRI(StandardPixelManager &r) : R(r) {}
	~StandardGRI() override;
	virtual size_t GetCoordinate(HRReal *Real, HRReal *Imaginary) override;
	virtual void WriteResults(SRReal *Values) override;
};
