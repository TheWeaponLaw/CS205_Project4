#include <cblas.h>
#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include <time.h>
#include "../include/function.h"
#define SIZE 1000
int main()
{
    int range = 2;
    clock_t start1, end1, start0, end0, start2, end2;
    start1 = clock();
    Matrix *matrix3 = createRam(SIZE, SIZE, range);
    sleep(1);
    Matrix *matrix4 = createRam(SIZE, SIZE, range);
    Matrix *matrix5 = createZero(SIZE, SIZE);
    Matrix *matrix6 = createZero(SIZE, SIZE);
    end1 = clock();
    printf("Create time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);

    // start1 = clock();
    // matmul_plain(matrix3, matrix4, matrix5);
    // end1 = clock();
    // printf("before time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);
    // showMatrix(matrix5);

    start1 = clock();
    matmul_improved(matrix3, matrix4, matrix5);
    end1 = clock();
    printf("improved time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);

    start1 = clock();
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SIZE, SIZE, SIZE, 1, matrix3->data, SIZE, matrix4->data, SIZE, 0, matrix6->data, SIZE);
    end1 = clock();
    printf("Cblas time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);

    printf("error: %f\n", test(matrix5->data, matrix6->data, SIZE));
    return 0;
}
