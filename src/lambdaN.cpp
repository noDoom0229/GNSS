/*-----------------------------------------------------------------------------
    lambdaN.cpp  —— Lambda 整数模糊度搜索（对应报告 2.6、3.2.4 lambdaN.cpp）

    实现 LAMBDA(MLAMBDA) 方法：
      1. Q = L D Lᵀ 分解
      2. 整数高斯变换 + 置换  ->  降相关 (Z 变换，报告 (2-39)(2-40))
      3. 在变换后的空间中做整数最小二乘搜索，找出 m 组最优解
      4. 反变换回原模糊度空间

    矩阵存储约定（与常见 LAMBDA 实现一致）：按列展开，A[i + j*n] 是第 i 行第 j 列。
    Q 是对称阵，所以按行 / 按列展开结果相同，调用方直接传按行展开的 Qnn 即可。
    输出 F 的第 k 组解存放在 F[k*n] ~ F[k*n + n - 1]。
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cmath>
#include <cstring>

#define LOOPMAX     1000000                    // 搜索循环上限（报告规定）
#define SGN(x)      ((x) <= 0.0 ? -1.0 : 1.0)  // 符号函数
#define ROUND(x)    (floor((x) + 0.5))         // 四舍五入
#define SWAP(x,y)   do { double tmp_ = x; x = y; y = tmp_; } while (0)

// ----------------------------------------------------------
//   (2) LD 分解：Q = Lᵀ D L（L 为单位下三角，D 为对角）
//   返回 0 成功，-1 失败（Q 非正定）
// ----------------------------------------------------------
static int LD(int n, const double* Q, double* L, double* D)
{
    int info = 0;
    vector<double> A(Q, Q + n * n);
    memset(L, 0, sizeof(double) * n * n);

    for (int i = n - 1; i >= 0; i--)
    {
        D[i] = A[i + i * n];
        if (D[i] <= 0.0) { info = -1; break; }
        double a = sqrt(D[i]);
        for (int j = 0; j <= i; j++) L[i + j * n] = A[i + j * n] / a;
        for (int j = 0; j <= i - 1; j++)
            for (int k = 0; k <= j; k++)
                A[j + k * n] -= L[i + k * n] * L[i + j * n];
        for (int j = 0; j <= i; j++) L[i + j * n] /= L[i + i * n];
    }
    if (info) printf("lambda: LD factorization error.\n");
    return info;
}

// ----------------------------------------------------------
//   (3) 整数高斯变换：用 L[i][j] 的四舍五入整数消去，使 |L[i][j]| <= 0.5
// ----------------------------------------------------------
static void gauss(int n, double* L, double* Z, int i, int j)
{
    int mu = (int)ROUND(L[i + j * n]);
    if (mu != 0)
    {
        for (int k = i; k < n; k++) L[k + n * j] -= (double)mu * L[k + i * n];
        for (int k = 0; k < n; k++) Z[k + n * j] -= (double)mu * Z[k + i * n];
    }
}

// ----------------------------------------------------------
//   (4) 置换：交换第 j 和 j+1 个模糊度，并更新 L、D、Z
// ----------------------------------------------------------
static void perm(int n, double* L, double* D, int j, double del, double* Z)
{
    double eta = D[j] / del;
    double lam = D[j + 1] * L[j + 1 + j * n] / del;
    D[j] = eta * D[j + 1];
    D[j + 1] = del;
    for (int k = 0; k <= j - 1; k++)
    {
        double a0 = L[j + k * n];
        double a1 = L[j + 1 + k * n];
        L[j + k * n] = -L[j + 1 + j * n] * a0 + a1;
        L[j + 1 + k * n] = eta * a0 + lam * a1;
    }
    L[j + 1 + j * n] = lam;
    for (int k = j + 2; k < n; k++) SWAP(L[k + j * n], L[k + (j + 1) * n]);
    for (int k = 0; k < n; k++)     SWAP(Z[k + j * n], Z[k + (j + 1) * n]);
}

// ----------------------------------------------------------
//   (5) Lambda 降相关：反复做高斯变换和置换，直到 D 按升序排列且 L 元素都很小
// ----------------------------------------------------------
static void reduction(int n, double* L, double* D, double* Z)
{
    int j = n - 2, k = n - 2;
    while (j >= 0)
    {
        if (j <= k)
            for (int i = j + 1; i < n; i++) gauss(n, L, Z, i, j);

        double del = D[j] + L[j + 1 + j * n] * L[j + 1 + j * n] * D[j + 1];
        if (del + 1E-6 < D[j + 1])
        {
            perm(n, L, D, j, del, Z);
            k = j;
            j = n - 2;
        }
        else
        {
            j--;
        }
    }
}

// ----------------------------------------------------------
//   (6) 整数最小二乘搜索：在椭球内找 m 组最优整数解 zn，残差 s（按升序）
//   返回 0 成功，-1 循环次数超过 LOOPMAX
// ----------------------------------------------------------
static int search(int n, int m, const double* L, const double* D,
    const double* zs, double* zn, double* s)
{
    int nn = 0, imax = 0, c;
    double newdist, maxdist = 1E99, y;
    vector<double> S(n * n, 0.0), dist(n, 0.0), zb(n, 0.0), z(n, 0.0), step(n, 0.0);

    int k = n - 1;
    dist[k] = 0.0;
    zb[k] = zs[k];
    z[k] = ROUND(zb[k]);
    y = zb[k] - z[k];
    step[k] = SGN(y);

    for (c = 0; c < LOOPMAX; c++)
    {
        newdist = dist[k] + y * y / D[k];
        if (newdist < maxdist)
        {
            if (k != 0)
            {
                // 向下一层继续搜索
                dist[--k] = newdist;
                for (int i = 0; i <= k; i++)
                    S[k + i * n] = S[k + 1 + i * n] + (z[k + 1] - zb[k + 1]) * L[k + 1 + i * n];
                zb[k] = zs[k] + S[k + k * n];
                z[k] = ROUND(zb[k]);
                y = zb[k] - z[k];
                step[k] = SGN(y);
            }
            else
            {
                // 到达最底层，得到一组候选解
                if (nn < m)
                {
                    if (nn == 0 || newdist > s[imax]) imax = nn;
                    for (int i = 0; i < n; i++) zn[i + nn * n] = z[i];
                    s[nn++] = newdist;
                }
                else
                {
                    if (newdist < s[imax])
                    {
                        for (int i = 0; i < n; i++) zn[i + imax * n] = z[i];
                        s[imax] = newdist;
                        for (int i = imax = 0; i < m; i++) if (s[imax] < s[i]) imax = i;
                    }
                    maxdist = s[imax];
                }
                z[0] += step[0];
                y = zb[0] - z[0];
                step[0] = -step[0] - SGN(step[0]);
            }
        }
        else
        {
            // 超出椭球，回到上一层
            if (k == n - 1) break;
            k++;
            z[k] += step[k];
            y = zb[k] - z[k];
            step[k] = -step[k] - SGN(step[k]);
        }
    }

    // 按残差升序排列
    for (int i = 0; i < m - 1; i++)
    {
        for (int j = i + 1; j < m; j++)
        {
            if (s[i] < s[j]) continue;
            SWAP(s[i], s[j]);
            for (int kk = 0; kk < n; kk++) SWAP(zn[kk + i * n], zn[kk + j * n]);
        }
    }

    if (c >= LOOPMAX)
    {
        printf("lambda: search loop count overflow.\n");
        return -1;
    }
    return 0;
}

// ----------------------------------------------------------
//   (7) 整数最小二乘估计
//   n 浮点模糊度个数，m 需要的候选解个数（Ratio 检验用 2）
//   a 浮点模糊度，Q 协因数阵，F 输出 m 组整数解，s 输出对应残差
//   返回 0 成功，-1 失败
// ----------------------------------------------------------
int lambda(int n, int m, const double* a, const double* Q, double* F, double* s)
{
    if (n <= 0 || m <= 0) return -1;

    vector<double> L(n * n, 0.0), D(n, 0.0), Z(n * n, 0.0), z(n, 0.0), E(n * m, 0.0);
    for (int i = 0; i < n; i++) Z[i + i * n] = 1.0;   // Z 初始化为单位阵

    // 1. LD 分解
    if (LD(n, Q, L.data(), D.data()) != 0) return -1;

    // 2. 降相关
    reduction(n, L.data(), D.data(), Z.data());

    // 3. 变换后的浮点模糊度 z = Zᵀ a
    for (int i = 0; i < n; i++)
    {
        z[i] = 0.0;
        for (int k = 0; k < n; k++) z[i] += Z[k + i * n] * a[k];
    }

    // 4. 搜索 m 组最优整数解
    if (search(n, m, L.data(), D.data(), z.data(), E.data(), s) != 0) return -1;

    // 5. 反变换 F = Z⁻ᵀ E（Zᵀ F = E）
    Mat Zt = mat_zero(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Zt[i][j] = Z[j + i * n];
    Mat ZtInv;
    if (!mat_inv(Zt, ZtInv))
    {
        printf("lambda: Z matrix is singular.\n");
        return -1;
    }
    for (int k = 0; k < m; k++)
    {
        for (int i = 0; i < n; i++)
        {
            double v = 0.0;
            for (int j = 0; j < n; j++) v += ZtInv[i][j] * E[j + k * n];
            F[i + k * n] = ROUND(v);   // Z 是整数幺模矩阵，结果应为整数，四舍五入消除浮点误差
        }
    }
    return 0;
}
