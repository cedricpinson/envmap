#pragma once

#include <chrono>

typedef std::chrono::time_point<std::chrono::steady_clock> Clock;

Clock logStart(const char* str);
void logEnd(Clock start);
