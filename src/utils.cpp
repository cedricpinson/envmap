#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <sys/stat.h>

int compressGZ(const char* path)
{
    Path gzFile;
    memcpy(gzFile, path, sizeof(Path));
    strcat(gzFile, ".gz");

    char command[1024];
    // remove the gz file if already exist
    remove(gzFile);
    snprintf(command, 1024, "7z a -y -tgzip %s %s >/dev/null", gzFile, path);
    int ret = system(command);
    if (ret != 0)
    {
        printf("error when running %s\n", command);
        return 0;
    }

    return getFileSize(gzFile);
}

const char* getFilenameExtension(const char* filename)
{
    const char* dot = strrchr(filename, '.');
    if (!dot || dot == filename)
        return "";
    return dot + 1;
}

void createPath(Path& path, const char* dir, const char* name)
{
    memset(path, 0, sizeof(Path));
    int l = strlen(dir);
    strncpy(path, dir, 512);
    if (dir[l - 1] != '/')
    {
        path[l] = '/';
    }
    strncat(path, name, 64);
}

int getFileSize(const char* path)
{
    struct stat st;
    stat(path, &st);
    int size = st.st_size;
    return size;
}

int makeDirectory(const char* dir)
{
    struct stat st
    {};

    if (stat(dir, &st) == -1)
    {
        mkdir(dir, 0700);
    }
    return 0;
}
