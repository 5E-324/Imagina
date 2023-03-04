#include "Includes.h"

bool Task::IsSynchronous() {
	return false;
}

Task::~Task() {}

std::string_view Task::GetDescription() const {
	using namespace std;
	return "(no description)"sv;
}

void ParallelTask::Execute() {
	throw std::exception("Not implemented");
}

void ParallelTask::Execute(unused size_t ThreadID) {
	Execute();
}

bool SynchronousTask::IsSynchronous() {
	return true;
}

void ExecutionContext::RemoveFromQueue() {
	while (!Computation::TaskQueue.empty() && Computation::TaskQueue.front().ReferenceCount == 0) {
		Computation::TaskQueue.pop_front();
	}
}

void ExecutionContext::Cancel() {
	Task::Cancellable *Cancellable = dynamic_cast<Task::Cancellable *>(task);
	std::unique_lock<std::mutex> lock(Mutex);
	AcceptNewThreads = false;
	Cancelled = true;
	if (ExecutingThreadCount == 0) {
		lock.unlock();
		lock.release();
		Release();
	} else if (Cancellable) {
		lock.unlock();
		Cancellable->Cancel();
	}
}

bool ExecutionContext::Finished() {
	return finished;
}

bool ExecutionContext::Terminated() {
	std::lock_guard<std::mutex> lock(Mutex);
	return ExecutingThreadCount == 0 && !AcceptNewThreads;
}

void ExecutionContext::Wait() {
	std::unique_lock lock(Mutex);

	ConditionVariable.wait(lock, [this]() { return ExecutingThreadCount == 0 && !AcceptNewThreads; });
}

Task *ExecutionContext::GetTask() {
	return task;
}

ExecutionContext::duration ExecutionContext::GetDuration() {
	std::lock_guard lock(Mutex);

	if (AcceptNewThreads && ExecutingThreadCount == 0) {
		return duration::zero();
	} else if (ExecutingThreadCount != 0) {
		return steady_clock::now() - StartTime;
	} else {
		return EndTime - StartTime;
	}
}

size_t ExecutionContext::AddReference() {
	return ++ReferenceCount;
}

size_t ExecutionContext::Release() {
	int32_t NewReferenceCount = --ReferenceCount;
	if (!NewReferenceCount) {
		if (!InTaskQueue) {
			delete this;
		} else {
			std::lock_guard<std::mutex> lock(Computation::TaskQueueMutex);
			RemoveFromQueue();
		}
	}
	return NewReferenceCount;
}

size_t ExecutionContext::WaitAndRelease() {
	Wait();
	return Release();
}

size_t Computation::WorkerCount = 0;
std::thread *Computation::Workers = nullptr;
std::deque<ExecutionContext> Computation::TaskQueue;

std::mutex Computation::TaskQueueMutex;
std::condition_variable Computation::TaskConditionVariable;

size_t Computation::HardwareConcurrency;

void Computation::WorkerFunction() {
	SetWorkerPriority();

	std::unique_lock<std::mutex> TaskQueueLock(TaskQueueMutex);

	while (true) {
		ExecutionContext *Context = nullptr;
		std::unique_lock<std::mutex> ContextLock;
		TaskConditionVariable.wait(TaskQueueLock, [&Context, &ContextLock, &TaskQueueLock]() { // Wait for a task to be available
			for (ExecutionContext &context : TaskQueue) {
				if (context.AcceptNewThreads) {
					ContextLock = std::unique_lock(context.Mutex, std::try_to_lock);
					if (!ContextLock.owns_lock()) { // Prevent dead lock
						TaskQueueLock.unlock();
						ContextLock.lock();
					}
					if (context.AcceptNewThreads) { // Still accepting new threads after locking?
						Context = &context;
						return true;
					}
					ContextLock.unlock(); // Try another one
					ContextLock.release();
					if (!TaskQueueLock.owns_lock()) TaskQueueLock.lock(); // The task queue may be unlocked. If so, lock it again.
				}
			}
			return false;
		});
		if (TaskQueueLock.owns_lock()) TaskQueueLock.unlock();

		size_t ThreadID = Context->ExecutingThreadCount;
		if (ThreadID == 0) {
			Context->StartTime = steady_clock::now();
		}
		Context->ExecutingThreadCount++;

		Task *task = Context->task;

		ParallelTask *parallelTask = dynamic_cast<ParallelTask *>(task);

		if (!parallelTask) Context->AcceptNewThreads = false;
		CooperativeTask *cooperativeTask = dynamic_cast<CooperativeTask *>(parallelTask);
		if (cooperativeTask && Context->ExecutingThreadCount >= cooperativeTask->ThreadCount) Context->AcceptNewThreads = false;

		ContextLock.unlock();

		DWORD_PTR originalAffinityMask = 0;
		if (cooperativeTask) { // PLATFORM DEPENDENT
			if (cooperativeTask->ThreadCount <= HardwareConcurrency / 2) {
				SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
				DWORD_PTR affinityMask = DWORD_PTR(0b11) << (ThreadID * 2);
				if (affinityMask) originalAffinityMask = SetThreadAffinityMask(GetCurrentThread(), affinityMask);
			} else {
				SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_NORMAL);
			}
		}

		if (parallelTask) {
			parallelTask->Execute(ThreadID);
		} else {
			task->Execute();
		}

		if (cooperativeTask) {
			if (originalAffinityMask) {
				SetThreadAffinityMask(GetCurrentThread(), originalAffinityMask);
			}
			SetWorkerPriority();
		}

		ContextLock.lock();

		Context->AcceptNewThreads = false; // At least one thread executing the task is terminated, stop accepting new threads
		Context->ExecutingThreadCount--;
		bool Notify = false;
		if (Context->ExecutingThreadCount == 0) {
			if (!Context->Cancelled) Context->finished = true;
			Context->EndTime = steady_clock::now();
			ContextLock.unlock(); // Unlock the context before locking the task queue
			ContextLock.release();
			TaskQueueLock.lock();
			size_t NewReferenceCount = --Context->ReferenceCount;
			if (NewReferenceCount == 0) {
				Context->RemoveFromQueue();
			} else {
				Notify = true; // Must lock the task queue befor notifying
			}
		} else {
			ContextLock.unlock();
			ContextLock.release();
			TaskQueueLock.lock();
		}
		if (Notify) { // The context will not be destructed while the task queue is locked
			Context->ConditionVariable.notify_all();
		}
	}
}

void Computation::Init() {
	HardwareConcurrency = std::max(std::thread::hardware_concurrency(), 1u);
#ifdef _DEBUG
	WorkerCount = 1;
#else
	WorkerCount = HardwareConcurrency;
#endif

	Workers = new std::thread[WorkerCount];
	for (size_t i = 0; i < WorkerCount; i++) {
		Workers[i] = std::thread(Computation::WorkerFunction);
	}
}

ExecutionContext *Computation::AddTask(Task *task) {
	ExecutionContext *Context = nullptr;
	if (task->IsSynchronous()) {
		Context = new ExecutionContext(task);
		Context->InTaskQueue = false;

		task->Execute();

		Context->AcceptNewThreads = false;
		Context->finished = true;
		return Context;
	}
	{
		std::lock_guard<std::mutex> lock(Computation::TaskQueueMutex);
		Context = &TaskQueue.emplace_back(task);
		Context->ReferenceCount = 2;
	}
	if (dynamic_cast<ParallelTask *>(task)) {
		TaskConditionVariable.notify_all();
	} else {
		TaskConditionVariable.notify_one();
	}
	return Context;
}

class FunctionTask : public Task {
	std::function<void()> function;

public:
	FunctionTask(std::function<void()> function) : function(function) {}

	virtual void Execute() override {
		function();
	}
};

ExecutionContext *Computation::AddTask(const std::function<void()> &function) {
	return AddTask(new FunctionTask(function));
}

std::vector<ExecutionContext *> Computation::GetExecutionContexts() {
	std::lock_guard<std::mutex> TaskQueueLock(TaskQueueMutex);
	std::vector<ExecutionContext *> Contexts;
	Contexts.reserve(TaskQueue.size());
	for (auto &Context : TaskQueue) {
		if (Context.ReferenceCount == 0) continue;
		Context.AddReference();
		Contexts.push_back(&Context);
	}
	return Contexts;
}