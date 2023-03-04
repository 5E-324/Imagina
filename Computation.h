#ifndef _COMPUTATION_H_
#define _COMPUTATION_H_

//#include "StandardIncludes.h"
#include "Types.h"

#include <thread>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <vector>
#include <chrono>
#include <functional>
using std::chrono::steady_clock;
class Task {
	friend class Computation;
private:
	virtual bool IsSynchronous();
public:
	virtual ~Task();
	virtual std::string_view GetDescription() const;
	virtual void Execute() = 0;

	class Cancellable;

};

class Task::Cancellable : public virtual Task {
public:
	virtual ~Cancellable() {}
	virtual void Cancel() = 0;
};

class ProgressTrackable {
public:
	virtual ~ProgressTrackable() {}
	virtual bool GetProgress(SRReal &Numerator, SRReal &Denoninator) const = 0;
};

class ParallelTask : public virtual Task {
public:
	virtual void Execute() override;
	virtual void Execute(size_t ThreadID);
};

class CooperativeTask : public ParallelTask {
	friend class Computation;
protected:
	const size_t ThreadCount;
	CooperativeTask(size_t ThreadCount) : ThreadCount(ThreadCount) {}
};

class SynchronousTask : public virtual Task {
private:
	virtual bool IsSynchronous() override final;
};

class ExecutionContext {
	friend class Computation;
private:
	std::atomic_uint32_t ReferenceCount = { 1 };
	Task *task;
	std::mutex Mutex;
	std::condition_variable ConditionVariable;

	steady_clock::time_point StartTime;
	steady_clock::time_point EndTime;

	size_t ExecutingThreadCount = 0;
	bool AcceptNewThreads = true;
	bool finished = false;
	bool Cancelled = false;
	bool InTaskQueue = true;

	void RemoveFromQueue();

public:
	using duration = steady_clock::duration;

	ExecutionContext(Task *task) : task(task) {}
	~ExecutionContext() { delete task; }

	void Cancel();
	bool Finished();
	bool Terminated();
	void Wait();
	template <class Representation, class Period>
	void Wait(std::chrono::duration<Representation, Period> duration) {
		std::unique_lock lock(Mutex);

		ConditionVariable.wait_for(lock, duration, [this]() { return ExecutingThreadCount == 0 && !AcceptNewThreads; });
	}
	Task *GetTask();
	duration GetDuration();

	size_t AddReference();
	size_t Release();
	size_t WaitAndRelease();
};

class Computation {
	friend class ExecutionContext;
private:
	static size_t WorkerCount;
	static std::thread *Workers;
	static std::deque<ExecutionContext> TaskQueue;

	static std::mutex TaskQueueMutex;
	static std::condition_variable TaskConditionVariable;

	static size_t HardwareConcurrency;

	static void WorkerFunction();

public:
	static void Init();

	static ExecutionContext *AddTask(Task *task);
	static ExecutionContext *AddTask(const std::function<void()> &function);
	static std::vector<ExecutionContext *> GetExecutionContexts();
};

#endif