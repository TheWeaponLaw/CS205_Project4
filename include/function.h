#ifndef MATRIX_FUNCTION_H
#define MATRIX_FUNCTION_H

struct Matrix;                                                                                            //矩阵结构体
struct Matrix *createZero(int row, int column);                                                           //创建一个全为0的矩阵
struct Matrix *createSpe(int row, int column, const float *data);                                          //根据一个string创建一个新的矩阵
struct Matrix *createRam(int row, int column);
void setMatrix(struct Matrix *matrix);                                                                    //设置矩阵里面的值
void deleteMatrix(struct Matrix **matrix);                                                                //删除矩阵
void copyMatrix(struct Matrix *matrix, const struct Matrix *matrixCopy);                                  //复制矩阵
void addMatrix(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix **matrix3);      //矩阵相加
void subtractMatrix(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix **matrix3); //矩阵相减
void addScalarMat(struct Matrix *matrix, const float scalar);                                             //矩阵加上一个值
void subScalarMat(struct Matrix *matrix, const float scalar);                                             //矩阵减去一个值
void mulScalarMat(struct Matrix *matrix, const float scalar);                                             //矩阵乘上一个值
void mulMatrix(const struct Matrix *matrix1, const struct Matrix *matrix2, struct Matrix **matrix3);      //矩阵相乘
void showMatrix(const struct Matrix *matrix);                                                             //展示矩阵
float findMin(struct Matrix *matrix);                                                                     //找最小值
float findMax(struct Matrix *matrix);                                                                     //找最大值
void tranMatrix(struct Matrix *matrix);                                                                   //矩阵转置
struct Matrix *invMatrix(struct Matrix *matrix);                                                          //矩阵求逆
#endif                                                                                                    // MATRIX_FUNCTION_H