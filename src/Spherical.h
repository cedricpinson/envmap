#pragma once

#define MAX_SH_ORDER 5
#define NUM_SH_COEFFICIENT (MAX_SH_ORDER * MAX_SH_ORDER)

typedef double Spherical[NUM_SH_COEFFICIENT * 3];
