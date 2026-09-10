#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <math.h>
#include <string.h>

// ==================== 矩阵运算函数 ====================

// 矩阵加法：C = A + B
void matrix_add(const double* A, const double* B, double* C, int rows, int cols)
{
    int i = 0;
    int j = 0;
    int idx = 0;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            idx = i * cols + j;//存储映射
            //^二维逻辑平铺成一维连续数组，GNSS 里矩阵都这么存，方便内存连续、计算更快。
            C[idx] = A[idx] + B[idx];
        }
    }
}


// 矩阵减法：C = A - B
void matrix_subtract(const double* A, const double* B, double* C, int rows, int cols)
{
    int i = 0;
    int j = 0;
    int idx = 0;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            idx = i * cols + j;
            C[idx] = A[idx] - B[idx];
        }
    }
}

// 矩阵乘法：C = A * B  (A: m×n, B: n×p, C: m×p)
void matrix_multiply(const double* A, const double* B, double* C, int m, int n, int p)
{
    int i = 0;   // C的行号（A的行号）
    int j = 0;   // C的列号（B的列号）
    int k = 0;   // 求和循环下标
    double sum = 0.0; // 存放内积累加结果

    // 外层循环：遍历 C 的每一行 i（A 的行）
    for (i = 0; i < m; i++)
    {
        // 中层循环：遍历 C 的每一列 j（B 的列）
        for (j = 0; j < p; j++)
        {
            sum = 0.0;  // 计算新一个C[i,j]前，累加器清零

            // 内层循环：A一行、B一列逐元素相乘求和
            for (k = 0; k < n; k++)
            {
                sum = sum + A[i * n + k] * B[k * p + j];// A[i * n + k]  A 每行有n个元素；i是行号，k是列号。
                /*A[i*n + k] = A[0*3 + 0] = A[0] → a00
                  B[k*p + j] = B[0*2 + 0] = B[0] → b00
                  sum = 0 + a00*b00*/
            }
            // 把求和结果存入 C[i][j]
            // C一维下标：C第i行，第j列 → i * p + j
            C[i * p + j] = sum;
        }
    }
}


// 矩阵×向量：v_out = A * v  (A: m×n, v: n×1, v_out: m×1)
void matrix_vector_multiply(const double* A, const double* v, double* v_out, int m, int n)
{
    int i = 0;   // 输出向量v_out下标，对应A的行号
    int k = 0;   // 内层遍历下标
    double sum = 0.0;

    // 外层循环：遍历A每一行（也就是v_out每一个元素）
    for (i = 0; i < m; i++)
    {
        sum = 0.0;    // 计算新元素前，累加器清零
        // 内层循环：取出A第i行所有元素，和v逐元素相乘累加
        for (k = 0; k < n; k++)
        {
            // A[i][k]一维下标：i*n + k
            sum = sum + A[i * n + k] * v[k];
        }
        v_out[i] = sum; // 点乘结果存入输出向量
    }
}

// 矩阵×常数：C = A * num
void matrix_num_multiply(const double* A, double num, double* C, int rows, int cols)
{
    int i = 0;
    int j = 0;
    int idx = 0;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            idx = i * cols + j;//二维位置->一维数组
            C[idx] = A[idx] * num;
        }
    }
}

// 矩阵转置：C = A^T  (A: m×n, C: n×m)
void matrix_trans(const double* A, double* C, int m, int n)
{
    int i = 0;
    int j = 0;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C[j * m + i] = A[i * n + j];
        }
    }
}

// 矩阵初始化：置为全 0
void matrix_initialize(double* mat, int rows, int cols)
{
    int total = 0;
    int i = 0;

    total = rows * cols;   // 计算矩阵总元素个数
    for (i = 0; i < total; i++)
    {
        mat[i] = 0.0;      // 每个位置赋值0,用的是一维数组的方式存储矩阵
    }
}

// 4x4矩阵求逆（高斯消元）
bool matrix_invert4x4(const double m[4][4], double inv[4][4])
{
    // ========== 步骤1：初始化4×8增广矩阵 mat[4][8] ==========
    double mat[4][8];
    for (int i = 0; i < 4; i++) {
        // 左4列：复制原矩阵m
        for (int j = 0; j < 4; j++)
            mat[i][j] = m[i][j];
        // 右4列：单位矩阵I，i == j-4 代表对角线位置赋值1
        for (int j = 4; j < 8; j++)
            mat[i][j] = (i == j - 4) ? 1.0 : 0.0;//条件 ? 条件成立取值 : 不成立取值
    }

    const double eps = 1e-12; // 极小值，用来判断元素是否等于0（浮点不能直接==0）

    // 逐列处理（每一列对应一个主元）
    for (int col = 0; col < 4; col++)
    {
        // ========== 步骤2：列主元选择（提升数值稳定性！重点） ==========
        int pivot = col;
        // 在当前列【col行往下】寻找绝对值最大元素作为主元
        for (int r = col; r < 4; r++)
        {
            if (fabs(mat[r][col]) > fabs(mat[pivot][col]))//pivot：保存当前找到最大元素所在行号，初始值 pivot = col
                pivot = r;//如果当前行的绝对值大于之前找到的最大值，则更新pivot为当前行号
        }

        // 主元接近0 → 矩阵奇异，无法求逆
        if (fabs(mat[pivot][col]) < eps)
            return false;

        // 如果最大元素不在当前行，交换【pivot行】与【col行】
        if (pivot != col) {
            double temp[8];
            //memcpy(目标地址, 源地址, 拷贝字节数)memcpy 直接整块内存批量复制，简洁高效。
            memcpy(temp, mat[col], sizeof(temp));
            memcpy(mat[col], mat[pivot], sizeof(temp));
            memcpy(mat[pivot], temp, sizeof(temp));
        }

        // ========== 步骤3：主元行归一化，让 mat[col][col] = 1 ==========
        double div = mat[col][col];//col：当前主元行、主元列编号；
        // 从主元列开始向右全部除以div  div是主元的值，归一化后主元位置变为1
        for (int j = col; j < 8; j++)
            mat[col][j] /= div;//a /= b等价a = a / b;

        // ========== 步骤4：消去本列其他所有行的元素，全部变成0 ==========
        for (int r = 0; r < 4; r++)
        {
            if (r != col && fabs(mat[r][col]) > eps)
            {
                double fac = mat[r][col];
                // 当前行 = 当前行 − fac × 主元行
                for (int j = col; j < 8; j++)
                    //第r行的第j列元素 = 第r行的第j列元素 - fac * 第col行的第j列元素
                    mat[r][j] -= fac * mat[col][j];//a -= b 等价于 a = a - b
            }
        }
    }

    // ========== 步骤5：提取逆矩阵，增广矩阵右侧4列就是A⁻¹ ==========
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            inv[i][j] = mat[i][j + 4];

    return true;
}

// 5x5矩阵求逆（高斯消元）
bool matrix_invert5x5(const double m[5][5], double inv[5][5])
{
    double mat[5][10];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) mat[i][j] = m[i][j];
        for (int j = 5; j < 10; j++) mat[i][j] = (i == j - 5) ? 1.0 : 0.0;
    }

    const double eps = 1e-12;
    for (int col = 0; col < 5; col++) {
        // 初始主元行为当前行 col
        int pivot = col;
        // 在当前列 [col ~ 4] 中寻找绝对值最大的元素作为主元
        for (int r = col; r < 5; r++) {
            if (fabs(mat[r][col]) > fabs(mat[pivot][col]))
                pivot = r;
        }
        // 若主元接近 0，说明矩阵奇异，无法求逆
        if (fabs(mat[pivot][col]) < eps)
            return false;

        // 主元行与当前行交换
        if (pivot != col) {
            double temp[10];
            memcpy(temp, mat[col], sizeof(temp));
            memcpy(mat[col], mat[pivot], sizeof(temp));
            memcpy(mat[pivot], temp, sizeof(temp));
        }

        // -------------------- 主元行归一化（将主元位置变为 1） --------------------
        double div = mat[col][col]; // 主元值
        for (int j = col; j < 10; j++)
            mat[col][j] /= div;

        // -------------------- 消去当前列其他所有行（变为 0） --------------------
        for (int r = 0; r < 5; r++) {
            // 跳过主元行，只处理其他行
            if (r != col && fabs(mat[r][col]) > eps) {
                double fac = mat[r][col]; // 消元系数
                // 整行消元，右侧单位矩阵同步变换
                for (int j = col; j < 10; j++)
                    mat[r][j] -= fac * mat[col][j];
            }
        }
    }

    // 3. 提取逆矩阵
    // 增广矩阵右侧 5 列即为原矩阵的逆矩阵 A⁻¹
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            inv[i][j] = mat[i][j + 5];

    return true;
}

// ==================== 向量运算函数 ====================

// 向量初始化：置为全 0
void vector_initialize(double* vec, int size)
{
    int i = 0;
    for (i = 0; i < size; i++)
    {
        vec[i] = 0.0;
    }
}

// 向量加法：c = a + b
void vector_add(const double* a, const double* b, double* c, int size)
{
    int i = 0;
    for (i = 0; i < size; i++)
    {
        c[i] = a[i] + b[i];
    }
}

// 向量减法：c = a - b
void vector_subtract(const double* a, const double* b, double* c, int size)
{
    int i = 0;
    for (i = 0; i < size; i++)
    {
        c[i] = a[i] - b[i];
    }
}

// 向量×常数：c = a * num
void vector_num_multiply(const double* a, double num, double* c, int size)
{
    int i = 0;
    for (i = 0; i < size; i++)
    {
        c[i] = a[i] * num;
    }
}

// 向量点乘：返回 a·b
double vector_multiply(const double* a, const double* b, int size)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < size; i++)
    {
        sum = sum + a[i] * b[i];
    }
    return sum;
}

// 3维向量叉乘：c = a × b
void vector_cross_multiply(const double a[3], const double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

// 向量模长计算
double vector_magnitude(const double* vec, int size)
{
    double sum = 0.0;
    int i = 0;
    for (i = 0; i < size; i++)
    {
        sum = sum + vec[i] * vec[i];
    }
    double result = sqrt(sum);
    return result;
}
