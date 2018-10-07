#pragma once

struct Light
{
    double2 _centroidPosition;

    // light Area
    double2 _size;
    double2 _position;

    // if one values is out of bound
    // flag the light as error
    bool _error = false;

    bool _merged = false;
    unsigned short _padding;
    int _mergedNum = 0;

    double _areaSize = 0;
    double _sum = 0;

    // average
    double _lumAverage = 0;
    double3 _colorAverage;
};
