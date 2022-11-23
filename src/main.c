#include <cblas.h>
#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include <time.h>
#include <omp.h>
#include "../include/function.h"
#define SIZE 2000
#define RANGE 100

#define START start = omp_get_wtime();
#define END(NAME)          \
    end = omp_get_wtime(); \
    printf("%s time is %lf\n", NAME, end - start);
int main()
{
    double start = 0;
    double end = 0;
    START
    Matrix *matrix3 = createRam(SIZE, SIZE, RANGE);
    sleep(1);
    Matrix *matrix4 = createRam(SIZE, SIZE, RANGE);
    Matrix *matrix5 = createZero(SIZE, SIZE);
    Matrix *matrix6 = createZero(SIZE, SIZE);
    END("Create")

    // start1 = clock();
    // matmul_plain(matrix3, matrix4, matrix5);
    // end1 = clock();
    // printf("before time: %lf\n", (double)(end1 - start1) / CLOCKS_PER_SEC);
    // showMatrix(matrix5);

    START
    matmul_improved(matrix3, matrix4, matrix5);
    END("Improve")

    START
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SIZE, SIZE, SIZE, 1, matrix3->data, SIZE, matrix4->data, SIZE, 0, matrix6->data, SIZE);
    END("Cblas")
    printf("error: %f\n", test(matrix5->data, matrix6->data, SIZE));
    return 0;
}
