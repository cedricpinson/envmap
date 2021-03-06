#pragma once

typedef char Path[1024];

const char* getFilenameExtension(const char* filename);
int getFileSize(const char* filename);
void createPath(Path& path, const char* dir, const char* name);
int makeDirectory(const char* dir);
int compressGZ(const char* path);
