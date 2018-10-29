#include "log.h"

#include <stdio.h>

Clock logStart(const char* str)
{
    printf("%s", str);
    auto start = std::chrono::steady_clock::now();
    return start;
}
void logEnd(Clock start)
{
    auto t = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = t - start;
    printf(" : %fs\n", duration.count());
}
