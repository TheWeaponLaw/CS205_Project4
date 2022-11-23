#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <immintrin.h>
#include "../include/function.h"

char judgeValid(const struct Matrix *matrix)
{
    if (matrix == NULL)
    {
        return 1;
    }
    else if (matrix->data == NULL)
    {
        return 2;
    }
    else if (matrix->row <= 0 || matrix->column <= 0)
    {
        return 4;
    }
    else
    {
        return 0;
    }
}

struct Matrix *createZero(int row, int column)
{
    //判断两个数据是否合法
    if (row <= 0 || column <= 0)
    {
        printf("The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        struct Matrix *matrix = (struct Matrix *)malloc(sizeof(int) * 2 + sizeof(float *));
        if (matrix == NULL)
        {
            printf("The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc((unsigned long long)(sizeof(float) * (unsigned long long)row * (unsigned long long)column));
        //判断申请空间是否成功
        if (matrix->data == NULL)
        {
            printf("The memory allocated failed\n");
            if (matrix != NULL)
                free(matrix);
            return NULL;
        }
        for (size_t i = 0; i < (unsigned long long)row * column; i++)
        {
            matrix->data[i] = 0;
        }
        return matrix;
    }
}

struct Matrix *createSpe(int row, int column, const float *data)
{
    if (row <= 0 || column <= 0)
    { //判断行列是否合法
        printf("The row or column are invalid!\n");
        return NULL;
    }
    else if (data == NULL)
    { //判断输入数据是否为空
        printf("Wrong data pointer!\n");
        return NULL;
    }
    else
    {
        struct Matrix *matrix = (struct Matrix *)malloc(sizeof(struct Matrix));
        if (matrix == NULL)
        { //判断内存是否申请成功
            printf("The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc(sizeof(float) * row * column);
        memcpy(matrix->data, data, sizeof(float) * row * column);
        return matrix;
    }
}

struct Matrix *createRam(int row, int column, int range)
{
    //判断两个数据是否合法
    if (row <= 0 || column <= 0)
    {
        printf("The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        struct Matrix *matrix = (struct Matrix *)malloc(sizeof(int) * 2 + sizeof(float *));
        if (matrix == NULL)
        {
            printf("The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc((sizeof(float) * (unsigned long long)row * (unsigned long long)column));
        //判断申请空间是否成功
        if (matrix->data == NULL)
        {
            printf("The memory allocated failed\n");
            if (matrix != NULL)
                free(matrix);
            return NULL;
        }
        srandom(time(NULL));
        for (size_t i = 0; i < (unsigned long long)row * column; i++)
        {
            matrix->data[i] = 1.0 * rand() / RAND_MAX * range;
        }
        return matrix;
    }
}

void deleteMatrix(struct Matrix **matrix)
{
    if (*matrix == NULL)
    { //是否struct结构体指针是否为空
        printf("The matrix doesn't exist for deleting!\n");
    }
    else
    { //数据指针是否为空
        if ((*matrix)->data == NULL)
        {
            printf("There's don't have data in matrix.\n");
        }
        else
        {
            free((*matrix)->data);
            (*matrix)->data = NULL;
        }
        //释放结构体其它申请空间
        free(*matrix);
        *matrix = NULL; //清除指针
    }
}

void showMatrix(const struct Matrix *matrix)
{
    if (judgeValid(matrix))
    {
        printf("The matrix isn't valid for showing!\n");
    }
    else
    {
        printf("The address of the matrix is: %p\n", matrix);
        for (size_t i = 0; i < matrix->row; i++)
        {
            for (size_t j = 0; j < matrix->column; j++)
            {
                if (fabsf(1.f - fabsf(matrix->data[i * (matrix->column) + j])) < 0.0001f)
                {
                    if (matrix->data[i * (matrix->column) + j] < 0)
                        matrix->data[i * (matrix->column) + j] = -1.f;
                    else
                        matrix->data[i * (matrix->column) + j] = 1.f;
                }
                printf("%-10g", matrix->data[i * (matrix->column) + j]);
            }
            printf("\n");
        }
        printf("\n");
    }
}

void matmul_plain(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix *matrix3)
{
    //判断两个相乘矩阵是否合法
    if (judgeValid(matrix1) || judgeValid(matrix2))
    {
        printf("The matrix aren't valid for multiplying!\n");
    }
    else if (matrix1->column != matrix2->row)
    { //判断矩阵1的列是否等于矩阵2的行
        printf("The two matrix aren't match for multiplying!\n");
    }
    else
    {
        //若matrix3的指针不存在或matrix3的数据指针申请空间与要求空间不一样
        if (matrix3 == NULL)
        {
            printf("The matrix3 doesn't exist!");
        }
        if (matrix1->row * matrix2->column != matrix3->row * matrix3->column) //_msize((*matrix3)->data))
        {                                                                     //若空间大小不同
            if (matrix3->data != NULL)
            { //释放之前matrix数据空间
                free(matrix3->data);
            }
            matrix3->data = (float *)malloc(sizeof(float) * matrix1->row * matrix2->column);
            if (matrix3->data == NULL)
            {
                printf("Memory allocated failed!\n");
                return;
            }
        }

        for (size_t i = 0; i < matrix1->row * matrix2->column; i++)
        {
            matrix3->data[i] = 0;
            for (size_t j = 0; j < matrix1->column; j++)
            {
                matrix3->data[i] += matrix1->data[i / matrix2->column * matrix1->column + j] *
                                    matrix2->data[i % matrix2->column + j * matrix2->column];
            }
        }
        matrix3->row = matrix1->row;
        matrix3->column = matrix2->column;
    }
}

void matmul_improved(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix *matrix3)
{
    //判断两个相乘矩阵是否合法
    if (judgeValid(matrix1) || judgeValid(matrix2))
    {
        printf("The matrix aren't valid for multiplying!\n");
    }
    else if (matrix1->column != matrix2->row)
    { //判断矩阵1的列是否等于矩阵2的行
        printf("The two matrix aren't match for multiplying!\n");
    }
    else
    {
        //若matrix3的指针不存在或matrix3的数据指针申请空间与要求空间不一样
        if (matrix3 == NULL)
        {
            printf("The matrix3 doesn't exist!");
        }
        if (matrix1->row * matrix2->column != matrix3->row * matrix3->column) //_msize((*matrix3)->data))
        {                                                                     //若空间大小不同
            if (matrix3->data != NULL)
            { //释放之前matrix数据空间
                free(matrix3->data);
            }
            matrix3->data = (float *)malloc(sizeof(float) * matrix1->row * matrix2->column);
            if (matrix3->data == NULL)
            {
                printf("Memory allocated failed!\n");
                return;
            }
        }
        matrix3->row = matrix1->row;
        matrix3->column = matrix2->column;

        float *temp = (float *)malloc(sizeof(float) * matrix1->row * matrix2->column);
        if (temp == NULL)
        {
            printf("Memory allocated failed!\n");
            return;
        }
        for (size_t i = 0; i < matrix1->row * matrix1->column; i += 8)
        {
            temp[i] = matrix2->data[(i % matrix1->row) * matrix1->row + i / matrix1->row];
            temp[i + 1] = matrix2->data[((i + 1) % matrix1->row) * matrix1->row + (i + 1) / matrix1->row];
            temp[i + 2] = matrix2->data[((i + 2) % matrix1->row) * matrix1->row + (i + 2) / matrix1->row];
            temp[i + 3] = matrix2->data[((i + 3) % matrix1->row) * matrix1->row + (i + 3) / matrix1->row];
            temp[i + 4] = matrix2->data[((i + 4) % matrix1->row) * matrix1->row + (i + 4) / matrix1->row];
            temp[i + 5] = matrix2->data[((i + 5) % matrix1->row) * matrix1->row + (i + 5) / matrix1->row];
            temp[i + 6] = matrix2->data[((i + 6) % matrix1->row) * matrix1->row + (i + 6) / matrix1->row];
            temp[i + 7] = matrix2->data[((i + 7) % matrix1->row) * matrix1->row + (i + 7) / matrix1->row];
        }
        __m256 a, b;
        __m256 c = _mm256_setzero_ps();
        float sum[8] = {0};
        float rest = 0;
        for (size_t i = 0; i < matrix1->row; i++)
        {
            for (size_t j = 0; j < matrix1->row; j++)
            {
                size_t lim_k = matrix1->row / 8 * 8;
                for (size_t k = 0; k < lim_k; k += 8)
                {
                    a = _mm256_loadu_ps(matrix1->data + i * matrix1->row + k);
                    b = _mm256_loadu_ps(temp + j * matrix1->row + k);
                    c = _mm256_add_ps(c, _mm256_mul_ps(a, b));
                }
                _mm256_storeu_ps(sum, c);
                c = _mm256_setzero_ps();

                matrix3->data[i * matrix1->row + j] = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7];
                //剩余部分
                for (int k = lim_k; k < matrix1->column; k++)
                {
                    rest += matrix1->data[i * matrix1->row + k] * matrix2->data[k * matrix1->column + j];
                }
                matrix3->data[i * matrix1->row + j] += rest;
                rest = 0;
            }
        }
        free(temp);
    }
}

float test(float *sample, float *test, size_t size)
{
    float sum[8] = {0};
    __m256 a, b;
    __m256 c = _mm256_setzero_ps();
    for (size_t i = 0; i < size * size; i += 8)
    {
        a = _mm256_loadu_ps(sample + i);
        b = _mm256_loadu_ps(test + i);
        c = _mm256_add_ps(c, _mm256_mul_ps(_mm256_sub_ps(a, b), _mm256_sub_ps(a, b)));
    }
    _mm256_storeu_ps(sum, c);

    return (sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7]);
}
