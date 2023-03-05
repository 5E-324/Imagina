#include "PlatformDependent.h"

#if __has_include(<Windows.h>)
#include <Windows.h>
#else
#include <windows.h>
#endif

void SetWorkerPriority() {
	SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_LOWEST);
}

void ErrorMessage(const char *Message) {
	MessageBoxA(nullptr, Message, nullptr, MB_TASKMODAL | MB_OK);
}

void ErrorMessage(const char *Title, const char *Message) {
	MessageBoxA(nullptr, Message, Title, MB_TASKMODAL | MB_OK);
}