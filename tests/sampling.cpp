#include "Cubemap.h"

void printVector(const char* name, float3 vec)
{
    printf("%s: %f, %f, %f\n", name, (double)vec[0], (double)vec[1], (double)vec[2]);
}

void printAddress(const char* name, Cubemap::Address address)
{
    printf("%s: face %d - s(%f) t(%f)\n", name, (int)address.face, (double)address.s, (double)address.t);
}

// original code
void regular(int x, int y, int size)
{
    // transform from [0..res - 1] to [- (1 - 1 / res) .. (1 - 1 / res)]
    // + 0.5f is for texel center addressing
    float nvcU = (2.0f * ((float)x + 0.5f) / (float)size) - 1.0f;
    float nvcV = (2.0f * ((float)y + 0.5f) / (float)size) - 1.0f;
    printf("regular %f,%f\n", (double)nvcU, (double)nvcV);
}

// https://seblagarde.wordpress.com/2012/06/10/amd-cubemapgen-for-physically-based-rendering/
// The Stretch edge fixup method of ModifiedCubemapgen is based on NVTT implementation [8]
void fixUp(int x, int y, int size)
{
    // transform from [0..res - 1] to [-1 .. 1], match up edges exactly.
    float nvcU = (2.0f * (float)x / ((float)size - 1.0f)) - 1.0f;
    float nvcV = (2.0f * (float)y / ((float)size - 1.0f)) - 1.0f;
    printf("fixup %f,%f\n", (double)nvcU, (double)nvcV);
}

int main( // int argc, char** argv
)
{

    Cubemap cm;

    int textureSize = 8;
    cm.size = textureSize;

    int y = 0;
    for (int x = 0; x < textureSize; x++)
    {
        float centerX = x + 0.5;
        float centerY = y + 0.5;

        printf("pixel %d,%d - center %f,%f\n", x, y, (double)centerX, (double)centerY);

        float3 dir, dirCenter;
        cm.getDirectionFor(dir.ptr(), Cubemap::Face::PX, x, y);
        printVector("dir", dir);
        cm.getDirectionFor(dirCenter.ptr(), Cubemap::Face::PX, centerX, centerY);
        printVector("dirCenter", dir);

        Cubemap::Address address;
        cm.getAddressFor(address, dirCenter.ptr());
        printAddress("address", address);
        printf("2d vector %f, %f\n", (double)address.s * textureSize, (double)address.t * textureSize);
        cm.getDirectionFor(dir.ptr(), Cubemap::Face::PX, textureSize + 1, y);
        printVector("dirCenter+1", dir);

        // data for fixup vector
        printf("data for fixup vector\n");
        cm.getDirectionFixUpFor(dir.ptr(), Cubemap::Face::PX, x, y);
        printVector("fixup vector", dir);
        Cubemap::getAddressFor(address, dir.ptr());
        printAddress("address fixup", address);
        float x0 = address.s * textureSize;
        float y0 = address.t * textureSize;
        printf("fetch %f,%f\n", (double)x0, (double)y0);


        printf("\n");
    }

    return 0;
}
