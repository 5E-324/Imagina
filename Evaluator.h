#ifndef _EVALUATOR_H_
#define _EVALUATOR_H_

#include "PixelManager.h"
#include "Computation.h"

#include <complex>

using real2 = HRReal;

enum class FractalTypeEnum {
	Mandelbrot,
	Tricorn,
	BurningShip,
	Nova,
	Custom,
};

struct ReferenceTrivialContent {
	HRReal AbsolutePrecision;
	HRReal RelativePrecision;
	HRReal ValidRadius;
};

struct Reference : ReferenceTrivialContent {
	Coordinate CenterCoordinate;

	virtual ~Reference();
};

class Evaluator {
public:
	class Feature {
	public:
		virtual ~Feature();
		virtual std::wstring_view Name();
		virtual std::wstring_view Information();
		virtual IntIter ItLimForZoomLevel(HRReal HalfH);
		virtual void GetCoordinate(HRReal &X, HRReal &Y) = 0;
		virtual bool GetRadius(HRReal &targetRadius);
		virtual bool CanLocatePrecisely();
	};

	class FeatureFinder {
	public:
		class PreciseLocatingTask : public virtual Task {
		public:
			Coordinate coordinate;

			virtual std::string_view GetDetailedProgress();
		};

		virtual ~FeatureFinder();
		virtual Feature *FindFeature(HRReal x, HRReal y, HRReal radius, Reference *reference) = 0;
		virtual PreciseLocatingTask *CreatePreciseLocatingTask(Feature *feature, size_t Precision, const Coordinate &centerCoordinate);
		virtual bool LocatePrecisely(Feature *feature, size_t Precision, Coordinate &featureCoordinate, const Coordinate &centerCoordinate);
	};

	class ReferenceTask : public virtual Task {
	public:
		Reference *reference;
	};

	virtual ~Evaluator();

	virtual ReferenceTask *CreateReferenceTask(const EvaluationParameters &parameters) = 0;
	virtual void FreeReference(Reference *reference);
	virtual Task *CreateEvaluationTask(const EvaluationParameters &parameters, PixelManager *pixelManager, Reference *reference) = 0;

	virtual FeatureFinder *GetFeatureFinder();

	virtual size_t SaveReference(std::ostream &stream, const Reference *reference);
	virtual size_t LoadReference(std::istream &stream, Reference *&reference);
};

class SimpleEvaluator : public Evaluator {
protected:
	class StandardReferenceTask : public ReferenceTask {
		SimpleEvaluator *evaluator;
		EvaluationParameters parameters;
	public:
		StandardReferenceTask(SimpleEvaluator *evaluator, const EvaluationParameters &parameters) : evaluator(evaluator), parameters(parameters) {}
		virtual std::string_view GetDescription() const override;
		virtual void Execute() override;
	};

	class StandardEvaluationTask : public ParallelTask, public ProgressTrackable {
		SimpleEvaluator *evaluator;
		EvaluationParameters parameters;
		PixelManager *rasterizer;
		Reference *reference;
	public:
		StandardEvaluationTask(SimpleEvaluator *evaluator, const EvaluationParameters &parameters, PixelManager *rasterizer, Reference *reference) : evaluator(evaluator), parameters(parameters), rasterizer(rasterizer), reference(reference) {}
		virtual std::string_view GetDescription() const override;
		virtual void Execute() override;
		virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const override;
	};

	virtual Reference *GenerateReference(const EvaluationParameters &parameters);
	virtual void Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *Reference) = 0;

public:
	virtual ReferenceTask *CreateReferenceTask(const EvaluationParameters &parameters);
	virtual Task *CreateEvaluationTask(const EvaluationParameters &parameters, PixelManager *pixelManager, Reference *reference);
};

template<FractalTypeEnum FractalType>
class PerturbationEvaluator : public SimpleEvaluator {
	static constexpr size_t VSize = VectorSize;

	static constexpr size_t CommitSize = 256 * 1024;
	size_t ReserveSize;

	struct _PTReference : Reference {
		const std::type_info &Type;
		size_t RefIt;
		size_t MaxIt;
		_PTReference(const std::type_info &type) : Type(type) {}
	};

	template <typename real2, FractalTypeEnum>
	struct PTReference : _PTReference {
		ReservedMemoryArray<real2> Refr, Refi;
		PTReference(size_t ReserveSize) : _PTReference(typeid(real2)), Refr(ReserveSize), Refi(ReserveSize) {}
	};

public:
	PerturbationEvaluator();

	template <typename real2>
	Reference *_GenerateReference(const EvaluationParameters &parameters);
	template <typename real2>
	void _Evaluate(PixelManager &rasterizer, Reference *Reference);
	virtual Reference *GenerateReference(const EvaluationParameters &parameters) override;
	virtual void Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *Reference) override;
};

class NovaEvaluator : public SimpleEvaluator {
	static constexpr size_t CommitSize = 256 * 1024;
	size_t ReserveSize;

	struct _PTReference : Reference {
		const std::type_info &Type;
		_PTReference(const std::type_info &type) : Type(type) {}
	};

	template <typename real2>
	struct PTReference : _PTReference {
		ReservedMemoryArray<real2> Zr;
		ReservedMemoryArray<std::complex<real2>> Zsqr;
		ReservedMemoryArray<std::complex<real2>> Zpow3m1;
		ReservedMemoryArray<std::complex<real2>> Zm1;
		size_t RefIt;
		size_t MaxIt;
		std::complex<real2> Refc;
		PTReference() : _PTReference(typeid(real2)) {}
	};

public:
	NovaEvaluator();

	template <typename real2>
	Reference *_GenerateReference(const Coordinate &coordinate);
	template <typename real2>
	void _Evaluate(PixelManager &rasterizer, Reference *Reference);
	virtual Reference *GenerateReference(const EvaluationParameters &parameters) override;
	virtual void Evaluate(const EvaluationParameters &parameters, PixelManager &rasterizer, Reference *Reference) override;
};

#endif