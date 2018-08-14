#include "utils.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>

const char* get_filename_ext(const char* filename)
{
    const char* dot = strrchr(filename, '.');
    if (!dot || dot == filename)
        return "";
    return dot + 1;
}

void create_path(Path& path, const char* dir, const char* name)
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

int make_directory(const char* dir)
{
    struct stat st = {0};

    if (stat(dir, &st) == -1)
    {
        mkdir(dir, 0700);
    }
    return 0;
}
