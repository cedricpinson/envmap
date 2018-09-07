#pragma once

#include <math.h>
#include <thread>

template <class T, class ArgObject> void threadLines(T* func, ArgObject context, int nbLines, int nbThread)
{
    if (nbThread > 1)
    {
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
            stopY = startY + (int)ceilf(linesPerThread);
            if (stopY > nbLines - 1)
                stopY = nbLines - 1;

            threadList[i] = std::thread(func, context, startY, stopY);

            startY = stopY + 1;
        }

        for (int i = 0; i < nbThread; i++)
        {
            threadList[i].join();
        }
    }
    else
    {
        (*func)(context, 0, nbLines - 1);
    }
}
