#pragma once
#include <string.h>

typedef char Path[1024];

inline const char* get_filename_ext(const char* filename)
{
    const char* dot = strrchr(filename, '.');
    if (!dot || dot == filename)
        return "";
    return dot + 1;
}

inline void create_path(Path& path, const char* dir, const char* name)
{
    memset(path, 0, 1024);
    int l = strlen(dir);
    strncpy(path, dir, 512);
    if (dir[l - 1] != '/')
    {
        path[l] = '/';
    }
    strncat(path, name, 64);
}
