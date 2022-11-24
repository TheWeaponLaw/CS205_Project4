#ifndef MATRIX_FUNCTION_H
#define MATRIX_FUNCTION_H

typedef struct Matrix
{
    size_t row;
    size_t column;
    float *data;
} Matrix;

struct Matrix *createZero(size_t row, size_t column);
struct Matrix *createSpe(size_t row, size_t column, const float *data);
struct Matrix *createRam(size_t row, size_t column, int range);
void deleteMatrix(struct Matrix **matrix); //删除矩阵

void matmul_plain(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix *matrix3); //矩阵相乘
void matmul_improved(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix *matrix3);
void showMatrix(const struct Matrix *matrix); //展示矩阵
float test(float *sample, float *test, size_t size);
float test_2(float *sample, float *test, size_t size);
void mul_matrix(size_t matrix1_row, size_t matrix1_col, float *matrix1, size_t matrix2_row, size_t matrix2_col, float *matrix2, float *matrix3);
float average(float *test, size_t size);
#endif