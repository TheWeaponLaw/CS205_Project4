#include <cblas.h>
#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
#include <time.h>
#include <omp.h>
#include "../include/function.h"
#define SIZE_ROW 1000
#define SIZE_COL 1000
#define RANGE 1000

#define START start = omp_get_wtime();
#define END(NAME)          \
    end = omp_get_wtime(); \
    printf("%s time is %lf\n", NAME, end - start);
int main()
{
    double start = 0;
    double end = 0;
    printf("Compute matrix of %d * %d\n", SIZE_ROW, SIZE_COL);
    START
    Matrix *matrix3 = createRam(SIZE_ROW, SIZE_COL, RANGE);
    sleep(1);
    Matrix *matrix4 = createRam(SIZE_COL, SIZE_ROW, RANGE);
    Matrix *matrix5 = createZero(SIZE_ROW, SIZE_ROW);
    Matrix *matrix6 = createZero(SIZE_ROW, SIZE_ROW);
    END("Create")

    // START
    // matmul_plain(matrix3, matrix4, matrix5);
    // END("Plain")

    START
    matmul_improved(matrix3, matrix4, matrix5);
    END("Improve")

    START
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, SIZE_ROW, SIZE_ROW, SIZE_COL, 1, matrix3->data, SIZE_COL, matrix4->data, SIZE_ROW, 0, matrix6->data, SIZE_ROW);
    END("Cblas")

    printf("error: %f\n", test_2(matrix5->data, matrix6->data, SIZE_ROW));
    printf("average: %f\n", average(matrix6->data, SIZE_ROW));
    // showMatrix(matrix5);
    // showMatrix(matrix6);
    // test_time();
    return 0;
}

double test_time()
{
    double start, end;
    double sum;
    for (int i = 0; i < 10; ++i)
    {
        Matrix *matrix1 = createRam(SIZE_ROW, SIZE_COL, RANGE);
        sleep(1);
        Matrix *matrix2 = createRam(SIZE_COL, SIZE_ROW, RANGE);
        Matrix *matrix3 = createZero(SIZE_ROW, SIZE_ROW);
        start = omp_get_wtime();
        matmul_improved(matrix1, matrix2, matrix3);
        end = omp_get_wtime();
        printf("%d improved time is %lf\n", i, end - start);
        sum += (end - start);
    }
    printf("The average time is %lf\n", sum / 10);
}
