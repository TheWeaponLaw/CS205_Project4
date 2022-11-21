#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "../include/function.h"

struct Matrix
{
    int row;
    int column;
    float *data;
};

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
    // else if (_msize(matrix->data) != sizeof(float) * matrix->column * matrix->row)
    // {
    //     return 3;
    // }
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

#define RANGE 3;

struct Matrix *createRam(int row, int column)
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
        srand((unsigned)time(NULL));
        for (size_t i = 0; i < (unsigned long long)row * column; i++)
        {
            matrix->data[i] = rand() % RANGE + 1;
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

void copyMatrix(struct Matrix *matrix, const struct Matrix *matrixCopy)
{
    if (judgeValid(matrixCopy))
    { //被复制的矩阵是否合法
        printf("The copy matrix is invalid!\n");
    }
    else if (matrix == NULL)
    { //指针为空
        printf("Don't exist such matrix!\n");
    }
    else
    {
        if (matrix->data == NULL || matrix->row * matrix->column != matrixCopy->column *
                                                                        matrix->row) //_msize(matrix->data) != _msize(matrixCopy->data))
        {                                                                            //若空间大小不同
            if (matrix->data != NULL)
            { //释放之前matrix数据空间
                free(matrix->data);
                matrix->data = NULL;
            }
            matrix->data = (float *)malloc(sizeof(float) * matrixCopy->row * matrixCopy->column);
        }
        memcpy(matrix->data, matrixCopy->data, sizeof(float) * matrixCopy->row * matrixCopy->column);
        matrix->row = matrixCopy->row;
        matrix->column = matrixCopy->column;
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

void mulMatrix(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix **matrix3)
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
        if (*matrix3 == NULL)
        {
            *matrix3 = createZero(matrix1->row, matrix2->column);
        }
        if (matrix1->row * matrix2->column != (*matrix3)->row * (*matrix3)->column) //_msize((*matrix3)->data))
        {                                                                           //若空间大小不同
            if ((*matrix3)->data != NULL)
            { //释放之前matrix数据空间
                free((*matrix3)->data);
            }
            (*matrix3)->data = (float *)malloc(sizeof(float) * matrix1->row * matrix2->column);
            if ((*matrix3)->data == NULL)
            {
                printf("Memory allocated failed!\n");
                return;
            }
        }
        for (size_t i = 0; i < matrix1->row * matrix2->column; i++)
        {
            (*matrix3)->data[i] = 0;
            for (size_t j = 0; j < matrix1->column; j++)
            {
                (*matrix3)->data[i] += matrix1->data[i / matrix2->column * matrix1->column + j] *
                                       matrix2->data[i % matrix2->column + j * matrix2->column];
            }
        }
        (*matrix3)->row = matrix1->row;
        (*matrix3)->column = matrix2->column;
    }
}