#pragma once

#include <math.h>
#include <thread>

//#define USE_THREAD

template <class T, class ArgObject> void threadLines(T* func, ArgObject context, int nbLines, int nbThread)
{

#ifdef USE_THREAD
    std::thread threadList[64];

    if (nbThread > 64)
        nbThread = 64;

    if (nbLines < nbThread)
    {
        nbThread = nbLines;
    }

    float linesPerThread = nbLines / nbThread;
    int startY = 0;
    int stopY;
    for (int i = 0; i < nbThread; i++)
    {
        stopY = startY + ceil(linesPerThread);
        if (stopY > nbLines - 1)
            stopY = nbLines - 1;

        threadList[i] = std::thread(func, context, startY, stopY);

        startY = stopY + 1;
    }

    for (int i = 0; i < nbThread; i++)
    {
        threadList[i].join();
    }
#else
    for (int i = 0; i < nbThread; i++)
        break;
    (*func)(context, 0, nbLines - 1);
#endif
}
