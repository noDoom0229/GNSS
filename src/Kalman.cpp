/*-----------------------------------------------------------------------------
    Kalman.cpp  —— 卡尔曼滤波 RTK 浮点解（对应报告 2.5、3.2.3 Kalman.cpp）

    状态向量 X = [X_R, Y_R, Z_R, N_dd(1) ... N_dd(n)]
      模糊度排列与 RTK.cpp 相同：[GPS f1 | GPS f2 | BDS f1 | BDS f2]

    状态方程   X(k) = Φ X(k-1) + w        (2-27)
    观测方程   L(k) = H X(k) + v          (2-28)
    预测       X(k|k-1) = Φ X(k-1|k-1)                     (2-32)
               P(k|k-1) = Φ P(k-1|k-1) Φᵀ + Q               (2-33)
    更新       K = P Hᵀ (H P Hᵀ + R)⁻¹                      (2-34)
               X(k|k) = X(k|k-1) + K (L - H X(k|k-1))       (2-35)
               P(k|k) = (I-KH) P (I-KH)ᵀ + K R Kᵀ            (2-36)
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cmath>

// ==========================================================
//   辅助：按 DDObs 填写卡尔曼结构体中的卫星信息 / 标签
// ==========================================================
static void set_kalman_struct(SDEpochObs& SDObs, DDCObs& DDObs, KalmanFilter& kf)
{
    kf.clear();
    kf.Time = SDObs.Time;
    kf.RefPRN[0] = DDObs.RefPrn[0];
    kf.RefPRN[1] = DDObs.RefPrn[1];

    for (int s = 0; s < 2; s++)
    {
        vector<int>& idxList = (s == 0) ? DDObs.GPSidx : DDObs.BDSidx;
        for (size_t k = 0; k < idxList.size(); k++)
        {
            SDSatObs& sd = SDObs.SdSatObs[idxList[k]];
            kf.PRN.push_back(sd.Prn);
            kf.System.push_back((int)sd.System);
        }
        // 模糊度标签：{系统, PRN, 频率}，顺序 f1 全部卫星，再 f2 全部卫星
        for (int f = 0; f < 2; f++)
        {
            for (size_t k = 0; k < idxList.size(); k++)
            {
                SDSatObs& sd = SDObs.SdSatObs[idxList[k]];
                vector<int> label;
                label.push_back((int)sd.System);
                label.push_back(sd.Prn);
                label.push_back(f);
                kf.Label.push_back(label);
            }
        }
    }
    kf.Dim = 3 + (int)kf.Label.size();
}

// 在 kf 的模糊度标签中查找 {sys, prn, f}，返回模糊度序号（不含前 3 个位置参数），找不到返回 -1
static int find_label(KalmanFilter& kf, int sys, int prn, int f)
{
    for (size_t i = 0; i < kf.Label.size(); i++)
    {
        if (kf.Label[i][0] == sys && kf.Label[i][1] == prn && kf.Label[i][2] == f)
            return (int)i;
    }
    return -1;
}

// 当前历元第 i 个模糊度能否由上一历元状态转移得到。
// 返回 1 表示可以（并给出转移系数：coefSelf 对应 pre 中同一卫星，coefRef 对应 pre 中新参考星），
// 返回 0 表示这是新升起的卫星（或参考星换成了新卫星），需要重新初始化。
static int amb_transition(KalmanFilter& pre, KalmanFilter& now, int i,
    int& jSelf, double& coefSelf, int& jRef, double& coefRef)
{
    int sys = now.Label[i][0];
    int prn = now.Label[i][1];
    int f = now.Label[i][2];
    int s = (sys == (int)GPS) ? 0 : 1;

    jSelf = jRef = -1;
    coefSelf = coefRef = 0.0;

    if (pre.RefPRN[s] == now.RefPRN[s])
    {
        // 参考星不变：N_now(prn) = N_pre(prn)
        jSelf = find_label(pre, sys, prn, f);
        if (jSelf < 0) return 0;
        coefSelf = 1.0;
        return 1;
    }

    // 参考星改变：N_now(prn, 相对新参考星) = N_pre(prn) - N_pre(新参考星)
    //   其中 N_pre(旧参考星) = 0
    jRef = find_label(pre, sys, now.RefPRN[s], f);
    if (jRef < 0) return 0;          // 新参考星上一历元没有被跟踪，无法转换
    coefRef = -1.0;

    if (prn == pre.RefPRN[s])
    {
        return 1;                    // 当前星就是旧参考星，只有 -N_pre(新参考星) 一项
    }
    jSelf = find_label(pre, sys, prn, f);
    if (jSelf < 0) return 0;
    coefSelf = 1.0;
    return 1;
}

// ==========================================================
//   (1) 卡尔曼滤波初始化  报告 3.2.3 (1)
// ==========================================================
// 位置用 SPP 结果，模糊度用 (ΔΔL - ΔΔP)/λ；
// 位置协方差 5 m、模糊度协方差 50 周（取平方作为方差）。
int kalman_initialize(SDEpochObs& SDObs, DDCObs& DDObs, PosRes& rovPosRes, KalmanFilter& nowKalman)
{
    set_kalman_struct(SDObs, DDObs, nowKalman);

    vector<double> amb0;
    if (!init_dd_ambiguity(SDObs, DDObs, amb0)) return 0;
    if ((int)amb0.size() != nowKalman.Dim - 3) return 0;

    nowKalman.X.assign(nowKalman.Dim, 0.0);
    nowKalman.X[0] = rovPosRes.Position[0];
    nowKalman.X[1] = rovPosRes.Position[1];
    nowKalman.X[2] = rovPosRes.Position[2];
    for (int i = 0; i < nowKalman.Dim - 3; i++) nowKalman.X[3 + i] = amb0[i];

    nowKalman.Xcov = mat_zero(nowKalman.Dim, nowKalman.Dim);
    for (int i = 0; i < 3; i++) nowKalman.Xcov[i][i] = KF_POS_SIGMA * KF_POS_SIGMA;
    for (int i = 3; i < nowKalman.Dim; i++) nowKalman.Xcov[i][i] = KF_AMB_SIGMA * KF_AMB_SIGMA;

    nowKalman.IsInit = true;
    return 1;
}

// ==========================================================
//   (2) 状态转移矩阵 Φ (nowDim x preDim)  报告 3.2.3 (2)
// ==========================================================
int get_phi(KalmanFilter& preKalman, KalmanFilter& nowKalman, vector<vector<double>>& phi)
{
    phi = mat_zero(nowKalman.Dim, preKalman.Dim);
    for (int i = 0; i < 3; i++) phi[i][i] = 1.0;     // 位置：单位阵

    for (int i = 0; i < nowKalman.Dim - 3; i++)
    {
        int jSelf, jRef;
        double cSelf, cRef;
        if (!amb_transition(preKalman, nowKalman, i, jSelf, cSelf, jRef, cRef))
            continue;                                  // 新卫星：整行为 0
        if (jSelf >= 0) phi[3 + i][3 + jSelf] += cSelf;
        if (jRef >= 0)  phi[3 + i][3 + jRef] += cRef;
    }
    return 1;
}

// ==========================================================
//   (3) 过程噪声 Q (nowDim x nowDim)  报告 3.2.3 (3)
// ==========================================================
// 位置过程噪声 5 m；连续跟踪卫星模糊度噪声取极小值；新升起卫星 5 周。
int get_Q(KalmanFilter& preKalman, KalmanFilter& nowKalman, vector<vector<double>>& Q)
{
    Q = mat_zero(nowKalman.Dim, nowKalman.Dim);
    for (int i = 0; i < 3; i++) Q[i][i] = KF_POS_NOISE * KF_POS_NOISE;

    for (int i = 0; i < nowKalman.Dim - 3; i++)
    {
        int jSelf, jRef;
        double cSelf, cRef;
        if (amb_transition(preKalman, nowKalman, i, jSelf, cSelf, jRef, cRef))
            Q[3 + i][3 + i] = KF_AMB_NOISE;
        else
            Q[3 + i][3 + i] = KF_NEW_AMB_NOISE * KF_NEW_AMB_NOISE;
    }
    return 1;
}

// ==========================================================
//   (4) 观测矩阵 H  报告 (2-16)(2-17)(2-21)
// ==========================================================
// 每颗非参考星 j（参考星 i）的方向系数：
//   l = (X_R - X^j)/ρ^j - (X_R - X^i)/ρ^i   （m、n 同理）           (2-18)~(2-20)
// 伪距行：[l m n | 0 ... 0]
// 相位行：[l m n | 0 ... λ ... 0]  λ 在对应模糊度列
// 流动站位置取 Xk_k1 的前 3 个分量（最小二乘调用时为当前迭代值）。
int get_H(EpochData& rovEpkObs, SDEpochObs& SDObs, DDCObs& DDObs,
    vector<double>& Xk_k1, vector<vector<double>>& H)
{
    int nG = (int)DDObs.GPSidx.size();
    int nB = (int)DDObs.BDSidx.size();
    int rows = 4 * nG + 4 * nB;
    int cols = 3 + 2 * nG + 2 * nB;
    if (rows == 0 || (int)Xk_k1.size() != cols) return 0;

    H = mat_zero(rows, cols);

    XYZCoord rovPos;
    rovPos.xyz[0] = Xk_k1[0];
    rovPos.xyz[1] = Xk_k1[1];
    rovPos.xyz[2] = Xk_k1[2];

    int row0 = 0;     // 当前系统观测值行起点
    int amb0 = 0;     // 当前系统模糊度列起点（不含前 3 列）
    for (int s = 0; s < 2; s++)
    {
        vector<int>& idxList = (s == 0) ? DDObs.GPSidx : DDObs.BDSidx;
        int n = (int)idxList.size();
        if (n == 0) continue;

        // 参考星方向余弦
        int refIdx = DDObs.RefPos[s];
        SDSatObs& ref = SDObs.SdSatObs[refIdx];
        int refRov = ref.nRov;
        DisRecSat dRef;
        if (!calc_rec_sat_dis(rovPos, rovEpkObs, refRov, dRef)) return 0;
        double eRef[3];
        for (int k = 0; k < 3; k++)
            eRef[k] = (rovPos.xyz[k] - rovEpkObs.SatPVT[refRov].SatPos[k]) / dRef.Dis;

        for (int j = 0; j < n; j++)
        {
            SDSatObs& sd = SDObs.SdSatObs[idxList[j]];
            int rovIdx = sd.nRov;
            DisRecSat dSat;
            if (!calc_rec_sat_dis(rovPos, rovEpkObs, rovIdx, dSat)) return 0;

            double lmn[3];
            for (int k = 0; k < 3; k++)
                lmn[k] = (rovPos.xyz[k] - rovEpkObs.SatPVT[rovIdx].SatPos[k]) / dSat.Dis - eRef[k];

            int rP1 = row0 + j;            // 伪距 f1
            int rP2 = row0 + n + j;        // 伪距 f2
            int rL1 = row0 + 2 * n + j;    // 相位 f1
            int rL2 = row0 + 3 * n + j;    // 相位 f2
            for (int k = 0; k < 3; k++)
            {
                H[rP1][k] = lmn[k];
                H[rP2][k] = lmn[k];
                H[rL1][k] = lmn[k];
                H[rL2][k] = lmn[k];
            }
            H[rL1][3 + amb0 + j] = get_wavelength(sd.System, 0);
            H[rL2][3 + amb0 + n + j] = get_wavelength(sd.System, 1);
        }
        row0 += 4 * n;
        amb0 += 2 * n;
    }
    return 1;
}

// ==========================================================
//   (5) 观测噪声 R  报告 3.2.3 (5)：与 P 矩阵结构对应，R = P⁻¹
// ==========================================================
// 非差等方差 σ² 时，n 个双差观测值的协方差阵（课件 II-3 P43）：
//   cov(DD) = 2σ² (I + 1·1ᵀ)   对角线 4σ²，非对角线 2σ²
static void fill_R_block(vector<vector<double>>& R, int start, int n, double sigma)
{
    double s2 = sigma * sigma;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            R[start + i][start + j] = (i == j) ? 4.0 * s2 : 2.0 * s2;
}

int get_R(int GPSDDnum, int BDSDDnum, vector<vector<double>>& R)
{
    int rows = 4 * GPSDDnum + 4 * BDSDDnum;
    if (rows <= 0) return 0;
    R = mat_zero(rows, rows);

    int start = 0;
    int n[2] = { GPSDDnum, BDSDDnum };
    for (int s = 0; s < 2; s++)
    {
        fill_R_block(R, start, n[s], RTK_SIGMA_CODE);  start += n[s];
        fill_R_block(R, start, n[s], RTK_SIGMA_CODE);  start += n[s];
        fill_R_block(R, start, n[s], RTK_SIGMA_PHASE); start += n[s];
        fill_R_block(R, start, n[s], RTK_SIGMA_PHASE); start += n[s];
    }
    return 1;
}

// ==========================================================
//   (6) 残差向量 V = 观测值 - 理论值  报告 (2-22)
// ==========================================================
// 伪距：ΔΔP - ΔΔρ
// 相位：ΔΔL - ΔΔρ - λ N
// 其中 ΔΔρ = (ρ_R^j - ρ_B^j) - (ρ_R^i - ρ_B^i)，基准站坐标已知不变。
// ΔΔP、ΔΔL、ΔΔρ 单位均为米（解码时相位已乘波长），N 为周，故用 λN 换算。
int get_V(EpochData& basEpkObs, EpochData& rovEpkObs, PosRes& basPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs, vector<double>& Xk_k1, vector<double>& V)
{
    int nG = (int)DDObs.GPSidx.size();
    int nB = (int)DDObs.BDSidx.size();
    int rows = 4 * nG + 4 * nB;
    int cols = 3 + 2 * nG + 2 * nB;
    if (rows == 0 || (int)Xk_k1.size() != cols) return 0;

    V.assign(rows, 0.0);

    XYZCoord rovPos, basPos;
    for (int k = 0; k < 3; k++)
    {
        rovPos.xyz[k] = Xk_k1[k];
        basPos.xyz[k] = basPosRes.Position[k];
    }

    int row0 = 0;
    int amb0 = 0;
    for (int s = 0; s < 2; s++)
    {
        vector<int>& idxList = (s == 0) ? DDObs.GPSidx : DDObs.BDSidx;
        int n = (int)idxList.size();
        if (n == 0) continue;

        SDSatObs& ref = SDObs.SdSatObs[DDObs.RefPos[s]];
        int refRov = ref.nRov, refBas = ref.nBas;
        DisRecSat dRefRov, dRefBas;
        if (!calc_rec_sat_dis(rovPos, rovEpkObs, refRov, dRefRov)) return 0;
        if (!calc_rec_sat_dis(basPos, basEpkObs, refBas, dRefBas)) return 0;
        double sdRhoRef = dRefRov.Dis - dRefBas.Dis;   // 参考星单差距离

        for (int j = 0; j < n; j++)
        {
            SDSatObs& sd = SDObs.SdSatObs[idxList[j]];
            int rovIdx = sd.nRov, basIdx = sd.nBas;
            DisRecSat dRov, dBas;
            if (!calc_rec_sat_dis(rovPos, rovEpkObs, rovIdx, dRov)) return 0;
            if (!calc_rec_sat_dis(basPos, basEpkObs, basIdx, dBas)) return 0;
            double ddRho = (dRov.Dis - dBas.Dis) - sdRhoRef;

            for (int f = 0; f < 2; f++)
            {
                double ddP = sd.dP[f] - ref.dP[f];
                double ddL = sd.dL[f] - ref.dL[f];
                double wl = get_wavelength(sd.System, f);
                double N = Xk_k1[3 + amb0 + f * n + j];

                V[row0 + f * n + j] = ddP - ddRho;                 // 伪距行
                V[row0 + 2 * n + f * n + j] = ddL - ddRho - wl * N; // 相位行
            }
        }
        row0 += 4 * n;
        amb0 += 2 * n;
    }
    return 1;
}

// ==========================================================
//   (7) 状态预测  (2-32)(2-33)
// ==========================================================
int kalman_predict(KalmanFilter& preKalman, KalmanFilter& nowKalman,
    vector<vector<double>>& Pk_k1, vector<double>& Xk_k1)
{
    if ((int)preKalman.X.size() != preKalman.Dim) return 0;
    if ((int)nowKalman.X.size() != nowKalman.Dim) return 0;

    Mat phi, Q;
    get_phi(preKalman, nowKalman, phi);
    get_Q(preKalman, nowKalman, Q);

    Xk_k1 = mat_mul_vec(phi, preKalman.X);
    Pk_k1 = mat_add(mat_mul(mat_mul(phi, preKalman.Xcov), mat_trans(phi)), Q);

    // 新升起的卫星：Φ 对应行全 0，预测值为 0，用初始化值 (ΔΔL-ΔΔP)/λ 赋值
    for (int i = 0; i < nowKalman.Dim - 3; i++)
    {
        int jSelf, jRef;
        double cSelf, cRef;
        if (!amb_transition(preKalman, nowKalman, i, jSelf, cSelf, jRef, cRef))
            Xk_k1[3 + i] = nowKalman.X[3 + i];
    }
    return 1;
}

// ==========================================================
//   (8) 测量更新  (2-34)(2-35)(2-36)
// ==========================================================
int kalman_update(EpochData& basEpkObs, EpochData& rovEpkObs, PosRes& basPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs, KalmanFilter& nowKalman,
    vector<vector<double>>& Pk_k1, vector<double>& Xk_k1,
    vector<vector<double>>& Pk, vector<double>& Xk)
{
    int nG = (int)DDObs.GPSidx.size();
    int nB = (int)DDObs.BDSidx.size();

    Mat H, R;
    Vec V;
    if (!get_H(rovEpkObs, SDObs, DDObs, Xk_k1, H)) return 0;
    if (!get_V(basEpkObs, rovEpkObs, basPosRes, SDObs, DDObs, Xk_k1, V)) return 0;
    if (!get_R(nG, nB, R)) return 0;

    Mat Ht = mat_trans(H);
    Mat PHt = mat_mul(Pk_k1, Ht);
    Mat S = mat_add(mat_mul(H, PHt), R);
    Mat Sinv;
    if (!mat_inv(S, Sinv))
    {
        printf("kalman_update: innovation covariance is singular.\n");
        return 0;
    }
    Mat K = mat_mul(PHt, Sinv);

    // X = X_pred + K V
    Vec KV = mat_mul_vec(K, V);
    Xk = Xk_k1;
    for (int i = 0; i < nowKalman.Dim; i++) Xk[i] += KV[i];

    // P = (I-KH) P (I-KH)ᵀ + K R Kᵀ
    Mat IKH = mat_sub(mat_eye(nowKalman.Dim), mat_mul(K, H));
    Pk = mat_add(mat_mul(mat_mul(IKH, Pk_k1), mat_trans(IKH)),
        mat_mul(mat_mul(K, R), mat_trans(K)));

    nowKalman.X = Xk;
    nowKalman.Xcov = Pk;
    return 1;
}

// ==========================================================
//   (9) 卡尔曼滤波浮点解  报告 3.2.3 (9)
// ==========================================================
int float_kalman(EpochData& basEpkObs, EpochData& rovEpkObs,
    PosRes& basPosRes, PosRes& rovPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs,
    KalmanFilter& preKalman, KalmanFilter& nowKalman)
{
    int nDD = (int)DDObs.GPSidx.size() + (int)DDObs.BDSidx.size();
    if (nDD < RTK_MIN_DD_SAT)
    {
        printf("float_kalman: not enough double-difference satellites (%d).\n", nDD);
        nowKalman.IsInit = false;
        return 0;
    }

    // 1. 用当前历元信息初始化结构体（也给新卫星提供初值）
    if (!kalman_initialize(SDObs, DDObs, rovPosRes, nowKalman)) return 0;

    // 2. 初始历元：只初始化，不滤波
    if (!preKalman.IsInit)
    {
        printf("float_kalman: filter initialized.\n");
        return 0;
    }

    // 3. 历元不连续：重新初始化
    double dt = (nowKalman.Time.Week - preKalman.Time.Week) * 604800.0
        + (nowKalman.Time.SecOfWeek - preKalman.Time.SecOfWeek);
    if (dt <= 0.0 || dt > KF_MAX_GAP)
    {
        printf("float_kalman: epoch gap %.1f s, filter re-initialized.\n", dt);
        return 0;
    }

    // 4. 预测 + 更新
    Mat Pk_k1, Pk;
    Vec Xk_k1, Xk;
    if (!kalman_predict(preKalman, nowKalman, Pk_k1, Xk_k1))
    {
        nowKalman.IsInit = false;
        return 0;
    }
    if (!kalman_update(basEpkObs, rovEpkObs, basPosRes, SDObs, DDObs, nowKalman, Pk_k1, Xk_k1, Pk, Xk))
    {
        nowKalman.IsInit = false;
        return 0;
    }

    // 5. 保存浮点解，为 Lambda 准备数据
    int n = DDObs.AmbNum;
    DDObs.dPos[0] = Xk[0] - basPosRes.Position[0];
    DDObs.dPos[1] = Xk[1] - basPosRes.Position[1];
    DDObs.dPos[2] = Xk[2] - basPosRes.Position[2];
    DDObs.Q = Pk;
    DDObs.FloatAmb.assign(Xk.begin() + 3, Xk.end());
    DDObs.Qnn.clear();
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            DDObs.Qnn.push_back(Pk[3 + i][3 + j]);

    rovPosRes.Position[0] = Xk[0];
    rovPosRes.Position[1] = Xk[1];
    rovPosRes.Position[2] = Xk[2];
    rovPosRes.Type = RTKFloat;
    return 1;
}
