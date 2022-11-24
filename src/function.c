#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <immintrin.h>
#include <malloc.h>
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

struct Matrix *createZero(size_t row, size_t column)
{
    //判断两个数据是否合法
    if (row <= 0 || column <= 0)
    {
        printf("The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        struct Matrix *matrix = (struct Matrix *)malloc(sizeof(size_t) * 2 + sizeof(float *));
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
        for (size_t i = 0; i < (unsigned long long)row * column; ++i)
        {
            matrix->data[i] = 0;
        }
        return matrix;
    }
}

struct Matrix *createSpe(size_t row, size_t column, const float *data)
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

struct Matrix *createRam(size_t row, size_t column, size_t range)
{
    //判断两个数据是否合法
    if (row <= 0 || column <= 0)
    {
        printf("The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        struct Matrix *matrix = (struct Matrix *)malloc(sizeof(size_t) * 2 + sizeof(float *));
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
        srand(time(NULL));
        for (size_t i = 0; i < (unsigned long long)row * column; ++i)
        {
            // matrix->data[i] = 1.0 * rand() / RAND_MAX * range;
            matrix->data[i] = rand() % range;
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
        for (size_t i = 0; i < matrix->row; ++i)
        {
            for (size_t j = 0; j < matrix->column; ++j)
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

        for (size_t i = 0; i < matrix1->row * matrix2->column; ++i)
        {
            matrix3->data[i] = 0;
            for (size_t j = 0; j < matrix1->column; ++j)
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
        // mul_matrix(matrix1->row, matrix1->column, matrix1->data, matrix2->row, matrix2->column, matrix2->data, matrix3->data);
        float *temp = (float *)malloc(sizeof(float) * matrix2->row * matrix2->column);
        if (temp == NULL)
        {
            printf("Memory allocated failed!\n");
            return;
        }
        for (size_t j = 0; j < matrix2->row; ++j)
        {
            for (size_t i = 0; i < matrix2->column; ++i)
            {
                temp[i * matrix2->row + j] = matrix2->data[j * matrix2->column + i];
            }
        }
        if (matrix1->row < 128)
        {
            mul_matrix(matrix1->row, matrix1->column, matrix1->data, matrix2->row, matrix2->column, temp, matrix3->data);
        }
        else
        {
#pragma omp parallel for
            for (int n = 0; n < 8; ++n)
            {
                mul_matrix(matrix1->row / 8, matrix1->column, matrix1->data + (matrix1->row / 8 * n) * matrix1->column, matrix2->row, matrix2->column, temp, matrix3->data + (matrix1->row / 8 * n) * matrix2->column);
            }
            mul_matrix(matrix1->row % 8, matrix1->column, matrix1->data + (matrix1->row / 8 * 8) * matrix1->column, matrix2->row, matrix2->column, temp, matrix3->data + (matrix1->row / 8 * 8) * matrix2->column);
        }
        free(temp);
    }
}

void mul_matrix(const size_t matrix1_row, const size_t matrix1_col, const float *matrix1, const size_t matrix2_row, const size_t matrix2_col, const float *temp, float *matrix3)
{
    __m256 a, b;
    __m256 c = _mm256_setzero_ps();
    float sum[8] = {0};
    for (size_t i = 0; i < matrix1_row; ++i)
    {
        for (size_t j = 0; j < matrix2_col; ++j)
        {
            size_t lim_k = matrix1_col / 8 * 8;
            for (size_t k = 0; k < lim_k; k += 8)
            {
                a = _mm256_loadu_ps(matrix1 + i * matrix1_col + k);
                b = _mm256_loadu_ps(temp + j * matrix2_row + k);
                c = _mm256_add_ps(c, _mm256_mul_ps(a, b));
            }
            _mm256_storeu_ps(sum, c);
            c = _mm256_setzero_ps();

            matrix3[i * matrix2_col + j] = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7];
            //剩余部分
            for (size_t k = lim_k; k < matrix1_col; ++k)
            {
                matrix3[i * matrix2_col + j] += matrix1[i * matrix1_col + k] * temp[j * matrix2_row + k];
            }
        }
    }
}

float average(float *test, size_t size)
{
    float ave = test[0];
    for (size_t i = 1; i < size * size; ++i)
    {
        ave = ave + (test[i] - ave) / (i + 1);
    }
    return (ave);
}

float test(float *sample, float *test, size_t size)
{
    float sum[8] = {0};
    __m256 a, b;
    __m256 c = _mm256_setzero_ps();
    for (size_t i = 0; i < (size * size) / 8 * 8; i += 8)
    {
        a = _mm256_loadu_ps(sample + i);
        b = _mm256_loadu_ps(test + i);
        c = _mm256_add_ps(c, _mm256_mul_ps(_mm256_sub_ps(a, b), _mm256_sub_ps(a, b)));
    }
    _mm256_storeu_ps(sum, c);
    float dif_square = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7];
    for (size_t i = (size * size) / 8 * 8; i < size * size; i += 8)
    {
        dif_square += (test[i] - sample[i]) * (test[i] - sample[i]);
    }

    return (dif_square);
}

float test_2(float *sample, float *test, size_t size)
{
    float dif = 0;
    for (size_t i = 0; i < size * size; ++i)
    {
        dif += fabsf(sample[i] - test[i]);
    }
    return dif;
}