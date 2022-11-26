#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <immintrin.h>
#include <malloc.h>
#include <stdbool.h>
#include "../include/function.h"

char judgeValid(const Matrix *matrix)
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

Matrix *createZero(size_t row, size_t column)
{
    Matrix *matrix = NULL;
    //判断两个数据是否合法
    if (row == 0 || column == 0)
    {
        fprintf(stderr, "The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        matrix = (Matrix *)malloc(sizeof(Matrix));
        if (matrix == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc((sizeof(float) * row * column));
        //判断申请空间是否成功
        if (matrix->data == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            if (matrix != NULL)
                free(matrix);
            return NULL;
        }
        for (size_t i = 0; i < row * column; ++i)
        {
            matrix->data[i] = 0;
        }
        return matrix;
    }
}

Matrix *createSpe(size_t row, size_t column, const float *data)
{
    Matrix *matrix = NULL;
    //判断两个数据是否合法
    if (row == 0 || column == 0)
    {
        fprintf(stderr, "The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        matrix = (Matrix *)malloc(sizeof(Matrix));
        if (matrix == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc((sizeof(float) * row * column));
        //判断申请空间是否成功
        if (matrix->data == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            if (matrix != NULL)
                free(matrix);
            return NULL;
        }
        memcpy(matrix->data, data, sizeof(float) * row * column);
        return matrix;
    }
}

Matrix *createRam(size_t row, size_t column, size_t range)
{
    Matrix *matrix = NULL;
    //判断两个数据是否合法
    if (row == 0 || column == 0)
    {
        fprintf(stderr, "The row or column are invalid!\n");
        return NULL;
    }
    else
    {
        matrix = (Matrix *)malloc(sizeof(Matrix));
        if (matrix == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            return NULL;
        }
        matrix->row = row;
        matrix->column = column;
        matrix->data = (float *)malloc((sizeof(float) * row * column));
        //判断申请空间是否成功
        if (matrix->data == NULL)
        {
            fprintf(stderr, "The memory allocated failed\n");
            if (matrix != NULL)
                free(matrix);
            return NULL;
        }
        for (size_t i = 0; i < row * column; ++i)
        {
            matrix->data[i] = 1.0 * rand() / RAND_MAX * range;
        }
        return matrix;
    }
}

bool deleteMatrix(Matrix **matrix)
{
    if (*matrix == NULL)
    { //是否struct结构体指针是否为空
        fprintf(stderr, "The matrix doesn't exist for deleting!\n");
        return false;
    }
    else
    { //数据指针是否为空
        if ((*matrix)->data == NULL)
        {
            fprintf(stderr, "There's don't have data in matrix.\n");
        }
        else
        {
            free((*matrix)->data);
            (*matrix)->data = NULL;
        }
        //释放结构体其它申请空间
        free(*matrix);
        *matrix = NULL; //清除指针
        return true;
    }
}

void showMatrix(const Matrix *matrix)
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

bool matmul_plain(const Matrix *matrix1, const Matrix *matrix2, Matrix *matrix3)
{
    //判断两个相乘矩阵是否合法
    if (judgeValid(matrix1) || judgeValid(matrix2))
    {
        fprintf(stderr, "The matrix aren't valid for multiplying!\n");
        return false;
    }
    else if (matrix1->column != matrix2->row)
    { //判断矩阵1的列是否等于矩阵2的行
        fprintf(stderr, "The two matrix aren't match for multiplying!\n");
        return false;
    }
    else
    {
        //若matrix3的指针不存在或matrix3的数据指针申请空间与要求空间不一样
        if (matrix3 == NULL)
        {
            fprintf(stderr, "The matrix3 doesn't exist!");
            return false;
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
                fprintf(stderr, "Memory allocated failed!\n");
                return false;
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
        return true;
    }
}

bool matmul_improved(const Matrix *matrix1, const Matrix *matrix2, Matrix *matrix3)
{
    //判断两个相乘矩阵是否合法
    if (judgeValid(matrix1) || judgeValid(matrix2))
    {
        fprintf(stderr, "The matrix aren't valid for multiplying!\n");
        return false;
    }
    else if (matrix1->column != matrix2->row)
    { //判断矩阵1的列是否等于矩阵2的行
        fprintf(stderr, "The two matrix aren't match for multiplying!\n");
        return false;
    }
    else
    {
        //若matrix3的指针不存在或matrix3的数据指针申请空间与要求空间不一样
        if (matrix3 == NULL)
        {
            fprintf(stderr, "The matrix3 doesn't exist!");
            return false;
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
                fprintf(stderr, "Memory allocated failed!\n");
                return false;
            }
        }
        matrix3->row = matrix1->row;
        matrix3->column = matrix2->column;
        float *temp = (float *)malloc(sizeof(float) * matrix2->row * matrix2->column);
        if (temp == NULL)
        {
            fprintf(stderr, "Memory allocated failed!\n");
            return false;
        }
        float *data_pre = matrix2->data; //修改处
        for (size_t j = 0; j < matrix2->row; ++j)
        {
            for (size_t i = 0; i < matrix2->column; ++i)
            {
                temp[i * matrix2->row + j] = *(data_pre++); //转置
            }
        }
        if (matrix1->row < 128)
        {
            mul_matrix(matrix1->row, matrix1->column, matrix1->data, matrix2->row, matrix2->column, temp, matrix3->data);
        }
        else
        {
#pragma omp parallel for
            for (size_t n = 0; n < 8; ++n)
            {
                mul_matrix(matrix1->row / 8, matrix1->column, matrix1->data + (matrix1->row / 8 * n) * matrix1->column, matrix2->row, matrix2->column, temp, matrix3->data + (matrix1->row / 8 * n) * matrix2->column);
            }
            mul_matrix(matrix1->row % 8, matrix1->column, matrix1->data + (matrix1->row / 8 * 8) * matrix1->column, matrix2->row, matrix2->column, temp, matrix3->data + (matrix1->row / 8 * 8) * matrix2->column);
        }
        free(temp);
        return true;
    }
}

void mul_matrix(const size_t matrix1_row, const size_t matrix1_col, const float *matrix1, const size_t matrix2_row, const size_t matrix2_col, const float *temp, float *matrix3)
{
    __m256 a, b;
    __m256 c = _mm256_setzero_ps();
    float sum[8] = {0};
    for (size_t i = 0; i < matrix1_row; ++i)
    {
        float *matrix1_add = matrix1 + i * matrix1_col; // matrix1的读取位置

        for (size_t j = 0; j < matrix2_col; ++j)
        {
            size_t lim_k = matrix1_col / 8 * 8;
            float *matrix2_add = temp + j * matrix2_row; // matrix2的读取位置

            for (size_t k = 0; k < lim_k; k += 8)
            {
                a = _mm256_loadu_ps(matrix1_add + k);
                b = _mm256_loadu_ps(matrix2_add + k);
                c = _mm256_add_ps(c, _mm256_mul_ps(a, b));
            }
            _mm256_storeu_ps(sum, c);
            c = _mm256_setzero_ps();

            matrix3[i * matrix2_col + j] = sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7];
            //剩余部分
            float *matrix1_data = matrix1 + i * matrix1_col;
            float *matrix2_data = temp + j * matrix2_row;
            for (size_t k = lim_k; k < matrix1_col; ++k)
            {
                matrix3[i * matrix2_col + j] += *(matrix1_data + k) * *(matrix2_data + k);
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