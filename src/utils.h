#pragma once

typedef char Path[1024];

const char* get_filename_ext(const char* filename);
void create_path(Path& path, const char* dir, const char* name);
int make_directory(const char* dir);
