#pragma once

struct Light;
struct Image;

int extractMainLight(Light& light, const Image& equirect);
