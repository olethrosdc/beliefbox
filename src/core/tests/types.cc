#include "real.h"
#include <cstdio>

// examine type sizes and conversions

int main()
{
    printf("Real: %d\n", sizeof(real));
    printf("Float: %d\n", sizeof(float));
    printf("Double: %d\n", sizeof(double));
    //printf("Long Real: %d\n", sizeof(long real));
    //printf("Long Float: %d\n", sizeof(long float));
    printf("Long Double: %d\n", sizeof(long double));
    printf("Int: %d\n", sizeof(int));
    printf("Short: %d\n", sizeof(short));
    printf("Long: %d\n", sizeof(long));
    //printf("Long long: %d\n", sizeof(long long));
    
    return 0;
}
