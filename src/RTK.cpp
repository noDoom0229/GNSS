/*-----------------------------------------------------------------------------
    RTK.cpp  —— RTK 相对定位核心（对应报告 3.2.2 RTK.cpp）

    包含：时间同步、站间单差、单差周跳探测、参考星选取、
          站星距离、权阵、最小二乘浮点解、Lambda 固定解。

    观测值 / 参数排列约定（所有函数共用，见报告 (2-21)~(2-23)）：
      非参考星顺序：DDObs.GPSidx 在前，DDObs.BDSidx 在后
      模糊度参数顺序：[GPS f1 (nG) | GPS f2 (nG) | BDS f1 (nB) | BDS f2 (nB)]
      观测值行顺序：  GPS: P f1 (nG), P f2 (nG), L f1 (nG), L f2 (nG)
                      BDS: P f1 (nB), P f2 (nB), L f1 (nB), L f2 (nB)
      注意：原有工程在解码时已把载波相位乘以波长转换为“米”，
            因此 dL 单位是米，模糊度 N 单位是周，相位方程中系数为波长 λ。
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cmath>
#include <cstring>

// ==========================================================
//   辅助函数
// ==========================================================
double get_wavelength(GNSSSys sys, int freq)
{
    if (sys == GPS) return (freq == 0) ? WL1_GPS : WL2_GPS;
    if (sys == BDS) return (freq == 0) ? WL1_BDS : WL3_BDS;
    return 0.0;
}

// BDS GEO 卫星：C01~C05, C59~C63
bool is_bds_geo(GNSSSys sys, int prn)
{
    if (sys != BDS) return false;
    if (prn >= 1 && prn <= 5) return true;
    if (prn >= 59 && prn <= 63) return true;
    return false;
}

static double time_diff(const GPSTime& t1, const GPSTime& t2)
{
    return (t1.Week - t2.Week) * 604800.0 + (t1.SecOfWeek - t2.SecOfWeek);
}

// ==========================================================
//   (0) 从文件 / 网络解码直到得到一个新的观测历元
// ==========================================================
int process_binary_file(ifstream& binFile, unsigned char* buff, size_t& bytesRead,
    EpochData& epk, Ephem* gpsEph, Ephem* bdsEph, BestPos& bestPos)
{
    while (true)
    {
        // 先尝试解码缓冲区中已有的数据（mode = 1：解出一个 RANGE 就返回）
        int len = (int)bytesRead;
        int ret = DecodeNovOem7Dat(buff, len, &epk, gpsEph, bdsEph, &bestPos, 1);
        bytesRead = (size_t)len;
        if (ret == 1) return 1;

        // 缓冲区已满却解不出完整消息：丢掉 1 个字节重新同步，避免死循环
        if (bytesRead >= MAXRAWLEN)
        {
            memmove(buff, buff + 1, bytesRead - 1);
            bytesRead--;
            continue;
        }

        // 读入更多数据
        if (!binFile.good()) return 0;
        binFile.read((char*)buff + bytesRead, MAXRAWLEN - bytesRead);
        streamsize n = binFile.gcount();
        if (n <= 0) return 0;        // 文件结束
        bytesRead += (size_t)n;
    }
}

int process_socket_data(SOCKET sock, unsigned char* buff, size_t& bytesRead,
    EpochData& epk, Ephem* gpsEph, Ephem* bdsEph, BestPos& bestPos)
{
    while (true)
    {
        int len = (int)bytesRead;
        int ret = DecodeNovOem7Dat(buff, len, &epk, gpsEph, bdsEph, &bestPos, 1);
        bytesRead = (size_t)len;
        if (ret == 1) return 1;

        if (bytesRead >= MAXRAWLEN)
        {
            memmove(buff, buff + 1, bytesRead - 1);
            bytesRead--;
            continue;
        }

        int n = recv(sock, (char*)buff + bytesRead, (int)(MAXRAWLEN - bytesRead), 0);
        if (n <= 0)
        {
            printf("Socket receive failed or connection closed.\n");
            return 0;
        }
        bytesRead += (size_t)n;
    }
}

// ==========================================================
//   (1) 文件流时间同步
// ==========================================================
// 思路（报告 3.2.2 (1) 与心得(1)）：两站分别独立解码，不嵌套。
//   进入函数时上一历元两站数据都已处理完，所以各自先读一个新历元；
//   然后比较时间：谁的时间小，就继续读谁，直到 |dt| < dT。
int get_synch_obs_file(
    ifstream& basBinFile, unsigned char* basBuff, size_t& basBytesRead,
    ifstream& rovBinFile, unsigned char* rovBuff, size_t& rovBytesRead,
    RTKData& rtkData, double dT)
{
    if (!process_binary_file(basBinFile, basBuff, basBytesRead, rtkData.BasEpkData,
        rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.BasBestPos))
        return 0;
    if (!process_binary_file(rovBinFile, rovBuff, rovBytesRead, rtkData.RovEpkData,
        rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.RovBestPos))
        return 0;

    while (true)
    {
        double dt = time_diff(rtkData.RovEpkData.Time, rtkData.BasEpkData.Time);
        if (fabs(dt) < dT)
        {
            rtkData.SdObs.Time = rtkData.RovEpkData.Time;
            return 1;
        }

        if (dt > 0)
        {
            // 基准站时间落后，继续读基准站
            if (!process_binary_file(basBinFile, basBuff, basBytesRead, rtkData.BasEpkData,
                rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.BasBestPos))
                return 0;
        }
        else
        {
            // 流动站时间落后，继续读流动站
            if (!process_binary_file(rovBinFile, rovBuff, rovBytesRead, rtkData.RovEpkData,
                rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.RovBestPos))
                return 0;
        }
    }
}

// ==========================================================
//   (2) 实时流时间同步：逻辑同上，只是数据来自 socket
// ==========================================================
int get_synch_obs_net(
    SOCKET NetGpsBas, unsigned char* basBuff, size_t& basBytesRead,
    SOCKET NetGpsRov, unsigned char* rovBuff, size_t& rovBytesRead,
    RTKData& rtkData, double dT)
{
    if (!process_socket_data(NetGpsBas, basBuff, basBytesRead, rtkData.BasEpkData,
        rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.BasBestPos))
        return 0;
    if (!process_socket_data(NetGpsRov, rovBuff, rovBytesRead, rtkData.RovEpkData,
        rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.RovBestPos))
        return 0;

    while (true)
    {
        double dt = time_diff(rtkData.RovEpkData.Time, rtkData.BasEpkData.Time);
        if (fabs(dt) < dT)
        {
            rtkData.SdObs.Time = rtkData.RovEpkData.Time;
            return 1;
        }

        if (dt > 0)
        {
            if (!process_socket_data(NetGpsBas, basBuff, basBytesRead, rtkData.BasEpkData,
                rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.BasBestPos))
                return 0;
        }
        else
        {
            if (!process_socket_data(NetGpsRov, rovBuff, rovBytesRead, rtkData.RovEpkData,
                rtkData.GPSEphemList.data(), rtkData.BDSEphemList.data(), rtkData.RovBestPos))
                return 0;
        }
    }
}

// ==========================================================
//   (3) 站间单差  报告 (2-5) (2-6)：ΔP = P_R - P_B，ΔL = L_R - L_B
// ==========================================================
int form_sd_obs(EpochData& basEpkData, EpochData& rovEpkData, SDEpochObs& SDObs, ConfigInfo& cfg)
{
    SDObs.Time = rovEpkData.Time;
    SDObs.SatNum = 0;

    for (int i = 0; i < rovEpkData.SatNum; i++)
    {
        SATOBSDATA& rov = rovEpkData.SatObs[i];
        if (rov.Prn == 0) continue;
        if (rov.System != GPS && rov.System != BDS) continue;

        // 在基准站中寻找同一颗卫星
        int j = -1;
        for (int k = 0; k < basEpkData.SatNum; k++)
        {
            if (basEpkData.SatObs[k].System == rov.System && basEpkData.SatObs[k].Prn == rov.Prn)
            {
                j = k;
                break;
            }
        }
        if (j < 0) continue;
        SATOBSDATA& bas = basEpkData.SatObs[j];

        // 观测值是否有效（粗差探测结果）
        if (!rov.Valid || !bas.Valid) continue;

        // 双频观测是否完整
        bool complete = true;
        for (int f = 0; f < 2; f++)
        {
            if (fabs(rov.P[f]) < 1e-5 || fabs(rov.L[f]) < 1e-5 ||
                fabs(bas.P[f]) < 1e-5 || fabs(bas.L[f]) < 1e-5)
                complete = false;
        }
        if (!complete) continue;

        // 双频伪距差过大（电离层 + 噪声正常只有几米）
        if (fabs(rov.P[0] - rov.P[1]) > cfg.PseuThreshold) continue;
        if (fabs(bas.P[0] - bas.P[1]) > cfg.PseuThreshold) continue;

        // 高度角（SatPVT.Elevation 单位为度，由 SPP 计算得到）
        if (!rovEpkData.SatPVT[i].Valid || !basEpkData.SatPVT[j].Valid) continue;
        if (rovEpkData.SatPVT[i].Elevation < cfg.ElevThreshold) continue;
        if (basEpkData.SatPVT[j].Elevation < cfg.ElevThreshold) continue;

        if (SDObs.SatNum >= MAXCHANNUM) break;
        SDSatObs& sd = SDObs.SdSatObs[SDObs.SatNum];
        sd.Prn = rov.Prn;
        sd.System = rov.System;
        sd.nBas = (short)j;
        sd.nRov = (short)i;

        for (int f = 0; f < 2; f++)
        {
            // ParityFlag = 0 表示可能存在半周，标记该频率不可用
            if (rov.half[f] == 0 || bas.half[f] == 0)
                sd.Valid[f] = 0;
            else
                sd.Valid[f] = 1;

            sd.dP[f] = rov.P[f] - bas.P[f];
            sd.dL[f] = rov.L[f] - bas.L[f];
        }
        SDObs.SatNum++;
    }
    return SDObs.SatNum;
}

// ==========================================================
//   (4) 单差观测值周跳探测  报告 (2-7) (2-8)
// ==========================================================
// MW (Melbourne-Wubbena) 组合：
//     MW = (f1*L1 - f2*L2)/(f1 - f2) - (f1*P1 + f2*P2)/(f1 + f2)
//   宽巷相位减窄巷伪距，消去了几何距离、电离层、钟差，只剩宽巷模糊度和噪声，
//   所以相邻历元 MW 应该基本不变；若变化超限，说明 L1 或 L2 发生了周跳。
// GF (Geometry-Free) 组合：
//     GF = L1 - L2
//   消去了几何距离和钟差，只剩电离层和模糊度项。站间单差后电离层几乎抵消，
//   所以相邻历元 GF 变化很小；变化超限同样说明有周跳。
// 两个组合对不同类型的周跳敏感度不同（MW 对 L1、L2 同时跳等量周不敏感，
// GF 可以补上），所以两者一起使用。阈值沿用原工程 GFThres / MWThres。
int detect_cycle_slip(SDEpochObs& sdEpk)
{
    ComObs cur[MAXCHANNUM];   // 当前历元的单差组合值
    int validNum = 0;

    for (int i = 0; i < sdEpk.SatNum; i++)
    {
        SDSatObs& sd = sdEpk.SdSatObs[i];
        if (sd.Valid[0] != 1 || sd.Valid[1] != 1) continue;   // 有半周，不参与

        double f1, f2;
        if (sd.System == GPS) { f1 = FG1_GPS; f2 = FG2_GPS; }
        else if (sd.System == BDS) { f1 = FG1_BDS; f2 = FG3_BDS; }
        else continue;

        double gf = sd.dL[0] - sd.dL[1];
        double mw = (f1 * sd.dL[0] - f2 * sd.dL[1]) / (f1 - f2)
            - (f1 * sd.dP[0] + f2 * sd.dP[1]) / (f1 + f2);

        cur[i].Prn = sd.Prn;
        cur[i].Sys = sd.System;
        cur[i].GF = gf;
        cur[i].MW = mw;
        cur[i].n = 1;

        // 查找上一历元同一颗卫星的组合值；第一次出现的卫星没有历史值，无法探测，按可用处理
        for (int j = 0; j < MAXCHANNUM; j++)
        {
            ComObs& pre = sdEpk.SdComObs[j];
            if (pre.n <= 0) continue;
            if (pre.Prn != sd.Prn || pre.Sys != sd.System) continue;

            double dGF = fabs(gf - pre.GF);
            double dMW = fabs(mw - pre.MW);
            if (dGF < GFThres && dMW < MWThres)
            {
                // 无周跳：MW 取历元平均值，抑制噪声
                cur[i].MW = (pre.MW * pre.n + mw) / (pre.n + 1);
                cur[i].n = pre.n + 1;
            }
            else
            {
                // 有周跳：该卫星本历元两个频率都标记为不可用，MW 重新开始平均
                sd.Valid[0] = sd.Valid[1] = 0;
            }
            break;
        }

        if (sd.Valid[0] == 1 && sd.Valid[1] == 1) validNum++;
    }

    // 保存当前历元组合值，供下一历元比较
    for (int i = 0; i < MAXCHANNUM; i++) sdEpk.SdComObs[i] = cur[i];
    return validNum;
}

// ==========================================================
//   (5) 参考星选取  报告 3.2.2 (5)
// ==========================================================
int DetRefSat(EpochData& basEpkData, EpochData& rovEpkData, SDEpochObs& SDObs, DDCObs& DDObs)
{
    DDObs.clear();

    double maxElev[2] = { -1.0, -1.0 };
    vector<int> usable[2];   // [0]=GPS [1]=BDS 可用单差卫星索引

    for (int i = 0; i < SDObs.SatNum; i++)
    {
        SDSatObs& sd = SDObs.SdSatObs[i];
        int s;
        if (sd.System == GPS) s = 0;
        else if (sd.System == BDS) s = 1;
        else continue;

        if (is_bds_geo(sd.System, sd.Prn)) continue;                          // 跳过 GEO
        if (!basEpkData.SatPVT[sd.nBas].Valid || !rovEpkData.SatPVT[sd.nRov].Valid) continue; // 位置计算失败
        if (sd.Valid[0] != 1 || sd.Valid[1] != 1) continue;                  // 半周 / 周跳

        usable[s].push_back(i);

        double elev = rovEpkData.SatPVT[sd.nRov].Elevation;
        if (elev > maxElev[s])
        {
            maxElev[s] = elev;
            DDObs.RefPrn[s] = sd.Prn;
            DDObs.RefPos[s] = i;
        }
    }

    for (int s = 0; s < 2; s++)
    {
        DDObs.SDNum[s] = (int)usable[s].size();
        if (DDObs.SDNum[s] < 2)
        {
            // 该系统只有参考星或没有卫星，无法组成双差
            DDObs.RefPrn[s] = DDObs.RefPos[s] = -1;
            DDObs.SDNum[s] = 0;
            continue;
        }
        for (size_t k = 0; k < usable[s].size(); k++)
        {
            int idx = usable[s][k];
            if (idx == DDObs.RefPos[s]) continue;
            if (s == 0) DDObs.GPSidx.push_back(idx);
            else        DDObs.BDSidx.push_back(idx);
        }
    }

    int nDD = (int)DDObs.GPSidx.size() + (int)DDObs.BDSidx.size();
    DDObs.AmbNum = 2 * nDD;   // 每颗非参考星两个频率各一个双差模糊度
    return nDD;
}

// ==========================================================
//   (9) 站星几何距离
// ==========================================================
int calc_rec_sat_dis(XYZCoord& recPos, EpochData& epkObs, int& idx, DisRecSat& disRecSat)
{
    if (idx < 0 || idx >= MAXCHANNUM) return 0;
    if (!epkObs.SatPVT[idx].Valid) return 0;

    double dx = epkObs.SatPVT[idx].SatPos[0] - recPos.xyz[0];
    double dy = epkObs.SatPVT[idx].SatPos[1] - recPos.xyz[1];
    double dz = epkObs.SatPVT[idx].SatPos[2] - recPos.xyz[2];

    disRecSat.Prn = epkObs.SatObs[idx].Prn;
    disRecSat.System = epkObs.SatObs[idx].System;
    disRecSat.Dis = sqrt(dx * dx + dy * dy + dz * dz);
    return disRecSat.Dis > 1.0 ? 1 : 0;
}

// ==========================================================
//   双差模糊度初值：N0 = (ΔΔL - ΔΔP) / λ   (周)
// ==========================================================
int init_dd_ambiguity(SDEpochObs& SDObs, DDCObs& DDObs, vector<double>& amb)
{
    amb.clear();
    for (int s = 0; s < 2; s++)
    {
        vector<int>& idxList = (s == 0) ? DDObs.GPSidx : DDObs.BDSidx;
        if (idxList.empty()) continue;
        int r = DDObs.RefPos[s];
        if (r < 0) return 0;
        SDSatObs& ref = SDObs.SdSatObs[r];

        for (int f = 0; f < 2; f++)
        {
            for (size_t k = 0; k < idxList.size(); k++)
            {
                SDSatObs& sd = SDObs.SdSatObs[idxList[k]];
                double ddL = sd.dL[f] - ref.dL[f];
                double ddP = sd.dP[f] - ref.dP[f];
                double wl = get_wavelength(sd.System, f);
                amb.push_back((ddL - ddP) / wl);
            }
        }
    }
    return (int)amb.size() == DDObs.AmbNum ? 1 : 0;
}

// ==========================================================
//   (8) 权阵  报告 (2-24)
// ==========================================================
// 同一系统、同一频率、同一类型的 n 个双差观测值来自 n+1 个单差观测值，
// 互相之间数学相关。若每个单差观测值方差为 σ²，则双差协方差阵为
//     D = σ² (I + 1·1ᵀ)     对角线 2σ²，非对角线 σ²
// 其逆（权阵）为
//     P = 1/σ² · 1/(n+1) · [ n  -1 ... ; -1  n ... ; ... ]
// 即对角线 n/(n+1)，非对角线 -1/(n+1)，再乘 1/σ²。
// 报告 (2-24) 中 σ = 1；这里伪距与相位分别用 RTK_SIGMA_CODE / RTK_SIGMA_PHASE。
static void fill_block(vector<vector<double>>& P, int start, int n, double sigma)
{
    if (n <= 0) return;
    double w = 1.0 / (sigma * sigma);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j) P[start + i][start + j] = w * (double)n / (n + 1.0);
            else        P[start + i][start + j] = -w / (n + 1.0);
        }
    }
}

int create_P_matrix(int nGPS, int nBDS, vector<vector<double>>& P)
{
    int rows = 4 * nGPS + 4 * nBDS;
    if (rows <= 0) return 0;
    P = mat_zero(rows, rows);

    int start = 0;
    int n[2] = { nGPS, nBDS };
    for (int s = 0; s < 2; s++)
    {
        // P f1, P f2, L f1, L f2
        fill_block(P, start, n[s], RTK_SIGMA_CODE);  start += n[s];
        fill_block(P, start, n[s], RTK_SIGMA_CODE);  start += n[s];
        fill_block(P, start, n[s], RTK_SIGMA_PHASE); start += n[s];
        fill_block(P, start, n[s], RTK_SIGMA_PHASE); start += n[s];
    }
    return 1;
}

// ==========================================================
//   (6) RTK 最小二乘浮点解  报告 2.4、3.2.2 (6)
// ==========================================================
// 参数向量 X = [X_R, Y_R, Z_R, N_dd ...]，
//   B 由 get_H() 构建（与卡尔曼观测矩阵相同，只是初值不同），
//   W 由 get_V() 构建（观测值 - 由当前参数计算的理论值），
//   P 由 create_P_matrix() 构建。
//   dX = (BᵀPB)⁻¹ BᵀPW,  Qxx = (BᵀPB)⁻¹                        (2-25)(2-26)
int LSS_RTK_float(RTKData& rtkData, PosRes& basPosRes, PosRes& rovPosRes)
{
    DDCObs& dd = rtkData.DDObs;
    SDEpochObs& sd = rtkData.SdObs;

    int nG = (int)dd.GPSidx.size();
    int nB = (int)dd.BDSidx.size();
    int nDD = nG + nB;
    if (nDD < RTK_MIN_DD_SAT)
    {
        printf("LSS_RTK_float: not enough double-difference satellites (%d).\n", nDD);
        return 0;
    }

    // 1. 参数初值：流动站位置用 SPP 结果，模糊度用 (ΔΔL - ΔΔP)/λ
    vector<double> amb0;
    if (!init_dd_ambiguity(sd, dd, amb0)) return 0;

    int dim = 3 + dd.AmbNum;
    vector<double> X(dim, 0.0);
    X[0] = rovPosRes.Position[0];
    X[1] = rovPosRes.Position[1];
    X[2] = rovPosRes.Position[2];
    for (int i = 0; i < dd.AmbNum; i++) X[3 + i] = amb0[i];

    // 2. 权阵
    Mat P;
    if (!create_P_matrix(nG, nB, P)) return 0;

    // 3. 迭代
    Mat Qxx;
    bool converged = false;
    for (int iter = 0; iter < RTK_LS_MAX_ITER; iter++)
    {
        Mat B;
        Vec W;
        if (!get_H(rtkData.RovEpkData, sd, dd, X, B)) return 0;
        if (!get_V(rtkData.BasEpkData, rtkData.RovEpkData, basPosRes, sd, dd, X, W)) return 0;

        Mat Bt = mat_trans(B);
        Mat BtP = mat_mul(Bt, P);
        Mat N = mat_mul(BtP, B);
        if (!mat_inv(N, Qxx))
        {
            printf("LSS_RTK_float: normal matrix is singular.\n");
            return 0;
        }
        Vec dX = mat_mul_vec(Qxx, mat_mul_vec(BtP, W));

        for (int i = 0; i < dim; i++) X[i] += dX[i];

        double norm = sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
        if (norm < RTK_LS_CONV_THRES)
        {
            converged = true;
            break;
        }
    }
    if (!converged)
    {
        printf("LSS_RTK_float: iteration did not converge.\n");
        return 0;
    }

    // 4. 保存浮点解，为 Lambda 准备数据
    dd.dPos[0] = X[0] - basPosRes.Position[0];
    dd.dPos[1] = X[1] - basPosRes.Position[1];
    dd.dPos[2] = X[2] - basPosRes.Position[2];
    dd.Q = Qxx;
    dd.FloatAmb.assign(X.begin() + 3, X.end());
    dd.Qnn.clear();
    for (int i = 0; i < dd.AmbNum; i++)
        for (int j = 0; j < dd.AmbNum; j++)
            dd.Qnn.push_back(Qxx[3 + i][3 + j]);

    rovPosRes.Position[0] = X[0];
    rovPosRes.Position[1] = X[1];
    rovPosRes.Position[2] = X[2];
    rovPosRes.Type = RTKFloat;
    return 1;
}

// ==========================================================
//   (7) RTK 固定解  报告 2.6、3.2.2 (7)
// ==========================================================
// 1) Lambda 得到两组最优整数解及残差 s[0] <= s[1]
// 2) Ratio = s[1] / s[0]，小于阈值则保持浮点解
// 3) 固定解基线 (2-41)：b_fix = b_float - Q_ba · Q_aa⁻¹ · (a_float - a_fix)
int RTK_fixed(DDCObs& ddObs, PosRes& rovPosRes, PosRes& basPosRes, ConfigInfo& cfg)
{
    int n = ddObs.AmbNum;
    if (n <= 0 || (int)ddObs.FloatAmb.size() != n || (int)ddObs.Qnn.size() != n * n)
        return 0;
    if ((int)ddObs.Q.size() != 3 + n) return 0;

    // 1. Lambda
    vector<double> F(n * 2, 0.0);
    double s[2] = { 0.0, 0.0 };
    if (lambda(n, 2, ddObs.FloatAmb.data(), ddObs.Qnn.data(), F.data(), s) != 0)
    {
        printf("RTK_fixed: lambda search failed.\n");
        return 0;
    }

    ddObs.FixedAmb.assign(F.begin(), F.begin() + n);
    ddObs.FixRMS.clear();
    ddObs.FixRMS.push_back(s[0]);
    ddObs.FixRMS.push_back(s[1]);

    // 2. Ratio 检验
    ddObs.Ratio = (s[0] > 0.0) ? (s[1] / s[0]) : 0.0;
    if (ddObs.Ratio < cfg.RatioThreshold)
    {
        return 0;   // 固定失败，保持 RTK Float
    }

    // 3. 固定解基线
    Mat Qaa = mat_zero(n, n);
    Mat Qba = mat_zero(3, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++) Qaa[i][j] = ddObs.Q[3 + i][3 + j];
        for (int k = 0; k < 3; k++) Qba[k][i] = ddObs.Q[k][3 + i];
    }
    Mat QaaInv;
    if (!mat_inv(Qaa, QaaInv))
    {
        printf("RTK_fixed: ambiguity covariance is singular.\n");
        return 0;
    }

    Vec da(n, 0.0);
    for (int i = 0; i < n; i++) da[i] = ddObs.FloatAmb[i] - ddObs.FixedAmb[i];

    Vec corr = mat_mul_vec(mat_mul(Qba, QaaInv), da);
    for (int k = 0; k < 3; k++) ddObs.dPos[k] -= corr[k];

    // 4. 流动站坐标 = 基准站坐标 + 固定解基线
    for (int k = 0; k < 3; k++)
        rovPosRes.Position[k] = basPosRes.Position[k] + ddObs.dPos[k];
    rovPosRes.Type = RTKFixed;
    return 1;
}
