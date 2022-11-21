
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "../include/function.h"

int main()
{
    clock_t start1, end1, start0, end0, start2, end2;
    start1 = clock();
    struct Matrix *matrix3 = createRam(8000, 8000);
    end1 = clock();
    printf("improve spend time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);
    return 0;
}
