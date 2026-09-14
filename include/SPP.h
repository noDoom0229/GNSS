#pragma once
/*-----------------------------------------------------------------------------
    SPP.h  —— 程序总头文件（对应报告 3.1.1 SPP.h）

    原有 SPP 项目的常量 / 结构体 / 函数声明分别放在
        GnssConsts.h、GnssStructs.h、GnssFuncDeclare.h
    这里统一包含它们，并在其基础上新增 RTK 相对定位所需的
        类型别名、结构体、常量和函数声明。
-----------------------------------------------------------------------------*/
#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"

#include <string>
#include <vector>
#include <fstream>
#include <cstdio>

using namespace std;

// ==========================================================
//   类型别名：报告中的名字  ->  原有工程中的结构体
// ==========================================================
typedef GPSTIME    GPSTime;    // GPS 时（周 + 周内秒）
typedef XYZ        XYZCoord;   // 地心直角坐标
typedef EPOCHOBS   EpochData;  // 一个历元的观测数据 + 卫星 PVT + 组合观测值
typedef GPSEPHREC  Ephem;      // 广播星历
typedef POSRES     BestPos;    // 解码得到的 BESTPOS/PSRPOS 结果
typedef PPRESULT   PosRes;     // 定位结果（SPP / RTK Float / RTK Fixed）
typedef MWGF       ComObs;     // MW / GF 组合观测值

// ==========================================================
//   RTK 相关常量（报告未给出具体数值的项集中放在这里，便于修改）
// ==========================================================
#define RTK_SYNC_DT        0.001   // 基准站/流动站时间同步阈值 (s)  课件 II-3 P19：对齐误差 < 0.001 s
#define RTK_MIN_DD_SAT     2       // 最少双差卫星数（4n 个观测 >= 3 + 2n 个未知数 => n >= 2）
#define RTK_LS_MAX_ITER    10      // 最小二乘最大迭代次数
#define RTK_LS_CONV_THRES  1e-4    // 最小二乘位置改正数收敛阈值 (m)
// 非差观测值中误差（等方差模型，课件 II-3 P42~P44 的权阵推导以此为前提）
#define RTK_SIGMA_CODE     0.3     // 非差伪距中误差 (m)
#define RTK_SIGMA_PHASE    0.003   // 非差相位中误差 (m)
#define KF_POS_SIGMA       5.0     // 卡尔曼初始化位置中误差 (m)        报告：5 m
#define KF_AMB_SIGMA       50.0    // 卡尔曼初始化模糊度中误差 (周)     报告：50 周
#define KF_POS_NOISE       5.0     // 位置过程噪声中误差 (m)            报告：5 m
#define KF_AMB_NOISE       1e-6    // 连续跟踪卫星模糊度过程噪声（极小值）
#define KF_NEW_AMB_NOISE   5.0     // 新升起卫星模糊度过程噪声中误差 (周) 报告：5 周
#define KF_MAX_GAP         1.5     // 相邻历元间隔超过该值 (s) 认为历元不连续，重新初始化滤波
#define GPS_UTC_LEAPSEC    18      // GPS 时与 UTC 之差 (s)，NMEA 输出 UTC 时使用

// ==========================================================
//   (2) 每颗卫星的单差观测值
// ==========================================================
struct SDSatObs
{
    short    Prn;
    GNSSSys  System;
    short    Valid[2];      // 每个频率单差观测值的质量：-1 未定, 0 有半周/周跳, 1 可用
    double   dP[2], dL[2];  // 单差伪距 / 单差相位 (m)   流动站 - 基准站
    short    nBas, nRov;    // 该卫星在基准站 / 流动站观测数组中的索引

    SDSatObs()
    {
        Prn = nBas = nRov = 0;
        System = UnknownSys;
        dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
        Valid[0] = Valid[1] = -1;
    }
};

// ==========================================================
//   (3) 一个历元的单差观测值
// ==========================================================
struct SDEpochObs
{
    GPSTime    Time;
    short      SatNum;
    SDSatObs   SdSatObs[MAXCHANNUM];
    ComObs     SdComObs[MAXCHANNUM];   // 上一历元的单差 MW/GF 组合值（周跳探测用，跨历元保留）

    SDEpochObs()
    {
        SatNum = 0;
    }

    void clear_except_comObs()
    {
        Time = GPSTime();
        SatNum = 0;
        for (int i = 0; i < MAXCHANNUM; i++)
        {
            SdSatObs[i] = SDSatObs();
        }
    }
};

// ==========================================================
//   (4) 双差信息
// ==========================================================
struct DDCObs
{
    int RefPrn[2], RefPos[2];          // [0]=GPS [1]=BDS 参考星 PRN 及其在 SdSatObs 中的索引
    int AmbNum, SDNum[2];              // 待估模糊度个数；GPS/BDS 可用单差卫星数（含参考星）
    vector<int> GPSidx, BDSidx;        // 参与双差的非参考星在 SdSatObs 中的索引
    double dPos[3];                    // 基线向量 (流动站 - 基准站)
    vector<vector<double>> Q;          // 参数协因数阵 (3+AmbNum)x(3+AmbNum)
    vector<double> FloatAmb;           // 浮点模糊度 (周)
    vector<double> Qnn;                // 模糊度协因数阵，按行展开 AmbNum*AmbNum
    vector<double> FixedAmb;           // 固定模糊度
    vector<double> FixRMS;             // Lambda 返回的两组最优解残差
    double Ratio;

    DDCObs()
    {
        for (int i = 0; i < 2; i++)
        {
            SDNum[i] = 0;
            RefPos[i] = RefPrn[i] = -1;
        }
        AmbNum = 0;
        dPos[0] = dPos[1] = dPos[2] = 0.0;
        Ratio = 0.0;
    }

    void clear()
    {
        for (int i = 0; i < 2; i++)
        {
            SDNum[i] = 0;
            RefPos[i] = RefPrn[i] = -1;
        }
        AmbNum = 0;
        dPos[0] = dPos[1] = dPos[2] = 0.0;
        Ratio = 0.0;
        GPSidx.clear(); BDSidx.clear();
        Q.clear(); FloatAmb.clear(); Qnn.clear(); FixedAmb.clear(); FixRMS.clear();
    }
};

// ==========================================================
//   (5) 卡尔曼滤波结构体
// ==========================================================
struct KalmanFilter
{
    GPSTime Time;
    int Dim;                      // 状态维数 = 3 + 模糊度个数
    int RefPRN[2];                // [0]=GPS [1]=BDS 参考星 PRN
    vector<int> PRN;              // 非参考星 PRN（与 Label 一一对应）
    vector<int> System;           // 非参考星系统
    vector<vector<int>> Label;    // 每个模糊度状态的标签 {系统, PRN, 频率}
    vector<double> X;             // 状态向量 [X Y Z N1 ... Nn]
    vector<vector<double>> Xcov;  // 状态协方差
    bool IsInit;

    KalmanFilter()
    {
        RefPRN[0] = RefPRN[1] = 0;
        Dim = 3;
        IsInit = false;
    }

    void clear()
    {
        RefPRN[0] = RefPRN[1] = 0;
        Dim = 3;
        PRN.clear(); System.clear(); Label.clear(); X.clear(); Xcov.clear();
    }
};

// ==========================================================
//   (6) RTK 相关信息
// ==========================================================
struct RTKData
{
    EpochData BasEpkData;
    EpochData RovEpkData;
    SDEpochObs SdObs;
    DDCObs DDObs;
    vector<Ephem> GPSEphemList{ MAXGPSNUM };
    vector<Ephem> BDSEphemList{ MAXBDSNUM };
    BestPos BasBestPos;
    BestPos RovBestPos;
    KalmanFilter PreKalman;
    KalmanFilter NowKalman;
};

// ==========================================================
//   (7) 配置信息
// ==========================================================
struct ConfigInfo
{
    string BasBinFile, RovBinFile;
    string BasIP, RovIP;
    int BasPort, RovPort;
    string NMEAOutputFile;
    string OutputFile;
    int PosMode;          // 数据来源：0 = 二进制文件，1 = 网络实时流
    int CalcMode;         // 浮点解方式：0 = 最小二乘，1 = 卡尔曼滤波
    double BasX, BasY, BasZ;   // 基准站已知坐标 (全为 0 时用基准站 SPP 结果)
    double RovX, RovY, RovZ;   // 流动站参考坐标，用于计算 dE/dN/dU (全为 0 时用流动站 BESTPOS)
    double ElevThreshold;      // 高度角阈值 (deg)
    double PseuThreshold;      // 双频伪距差阈值 (m)
    double RatioThreshold;     // Ratio 检验阈值

    ConfigInfo()
    {
        BasPort = RovPort = 0;
        PosMode = CalcMode = 0;
        BasX = BasY = BasZ = RovX = RovY = RovZ = 0.0;
        ElevThreshold = 10.0;
        PseuThreshold = 10.0;
        RatioThreshold = 3.0;
    }
};

// ==========================================================
//   (9) 精度评定信息（课件 II-2 P19「精度评定与结果输出模块」，
//        II-3 P54~P56 精度评定公式）
// ==========================================================
struct QualityInfo
{
    double Sigma0;   // 单位权中误差 sqrt(VᵀPV / r)
    double RMS;      // 观测值残差 RMS sqrt(VᵀV / n)
    double mENU[3];  // E / N / U 方向中误差 (m)
    int    nObs;     // 双差观测值个数
    int    nPar;     // 待估参数个数
    bool   Valid;

    QualityInfo()
    {
        Sigma0 = RMS = 0.0;
        mENU[0] = mENU[1] = mENU[2] = 0.0;
        nObs = nPar = 0;
        Valid = false;
    }
};

// ==========================================================
//   (8) 接收机到卫星距离
// ==========================================================
struct DisRecSat
{
    short    Prn;
    GNSSSys  System;
    double   Dis;

    DisRecSat()
    {
        Prn = 0;
        System = UnknownSys;
        Dis = 0.0;
    }
};

// ==========================================================
//   config.cpp
// ==========================================================
bool load_config(const string& filename, ConfigInfo& config);

// ==========================================================
//   RTK.cpp
// ==========================================================
// 读取二进制文件直到解出一个新的观测历元；返回 1 成功，0 文件结束
int process_binary_file(ifstream& binFile, unsigned char* buff, size_t& bytesRead,
    EpochData& epk, Ephem* gpsEph, Ephem* bdsEph, BestPos& bestPos);
// 从网络接收数据直到解出一个新的观测历元；返回 1 成功，0 连接断开
int process_socket_data(SOCKET sock, unsigned char* buff, size_t& bytesRead,
    EpochData& epk, Ephem* gpsEph, Ephem* bdsEph, BestPos& bestPos);

int get_synch_obs_file(
    ifstream& basBinFile, unsigned char* basBuff, size_t& basBytesRead,
    ifstream& rovBinFile, unsigned char* rovBuff, size_t& rovBytesRead,
    RTKData& rtkData, double dT);

int get_synch_obs_net(
    SOCKET NetGpsBas, unsigned char* basBuff, size_t& basBytesRead,
    SOCKET NetGpsRov, unsigned char* rovBuff, size_t& rovBytesRead,
    RTKData& rtkData, double dT);

int form_sd_obs(EpochData& basEpkData, EpochData& rovEpkData, SDEpochObs& SDObs, ConfigInfo& cfg);
int detect_cycle_slip(SDEpochObs& sdEpk);
int DetRefSat(EpochData& basEpkData, EpochData& rovEpkData, SDEpochObs& SDObs, DDCObs& DDObs);
int calc_rec_sat_dis(XYZCoord& recPos, EpochData& epkObs, int& idx, DisRecSat& disRecSat);
int create_P_matrix(int nGPS, int nBDS, vector<vector<double>>& P);
int LSS_RTK_float(RTKData& rtkData, PosRes& basPosRes, PosRes& rovPosRes);
int RTK_fixed(DDCObs& ddObs, PosRes& rovPosRes, PosRes& basPosRes, ConfigInfo& cfg);
// 精度评定：在当前解上重建 B / P / V，计算单位权中误差、RMS 与 ENU 中误差
int calc_rtk_quality(RTKData& rtkData, PosRes& basPosRes, PosRes& rovPosRes,
    int calcMode, QualityInfo& quality);

// 辅助：某频率的载波波长 (m)
double get_wavelength(GNSSSys sys, int freq);
// 辅助：是否为 BDS GEO 卫星
bool is_bds_geo(GNSSSys sys, int prn);
// 辅助：双差模糊度初值 (L_dd - P_dd) / lambda，按 [GPS f1, GPS f2, BDS f1, BDS f2] 顺序
int init_dd_ambiguity(SDEpochObs& SDObs, DDCObs& DDObs, vector<double>& amb);

// ==========================================================
//   Kalman.cpp
// ==========================================================
int kalman_initialize(SDEpochObs& SDObs, DDCObs& DDObs, PosRes& rovPosRes, KalmanFilter& nowKalman);
int get_phi(KalmanFilter& preKalman, KalmanFilter& nowKalman, vector<vector<double>>& phi);
int get_Q(KalmanFilter& preKalman, KalmanFilter& nowKalman, vector<vector<double>>& Q);
int get_H(EpochData& rovEpkObs, SDEpochObs& SDObs, DDCObs& DDObs,
    vector<double>& Xk_k1, vector<vector<double>>& H);
int get_R(int GPSDDnum, int BDSDDnum, vector<vector<double>>& R);
int get_V(EpochData& basEpkObs, EpochData& rovEpkObs, PosRes& basPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs, vector<double>& Xk_k1, vector<double>& V);
int kalman_predict(KalmanFilter& preKalman, KalmanFilter& nowKalman,
    vector<vector<double>>& Pk_k1, vector<double>& Xk_k1);
int kalman_update(EpochData& basEpkObs, EpochData& rovEpkObs, PosRes& basPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs, KalmanFilter& nowKalman,
    vector<vector<double>>& Pk_k1, vector<double>& Xk_k1,
    vector<vector<double>>& Pk, vector<double>& Xk);
int float_kalman(EpochData& basEpkObs, EpochData& rovEpkObs,
    PosRes& basPosRes, PosRes& rovPosRes,
    SDEpochObs& SDObs, DDCObs& DDObs,
    KalmanFilter& preKalman, KalmanFilter& nowKalman);

// ==========================================================
//   lambdaN.cpp
// ==========================================================
// n 维浮点模糊度 a、协因数阵 Q(n*n, 按行展开)，搜索 m 组最优整数解 F(n*m)，s 为对应残差
// 返回 0 成功，-1 失败
int lambda(int n, int m, const double* a, const double* Q, double* F, double* s);

// ==========================================================
//   writeToFile.cpp
// ==========================================================
void write_header_RTK(FILE* fp);
void write_2_screen_RTK(PosRes& rovPosRes, PosRes& basPosRes, DDCObs& ddObs, SDEpochObs& sdObs,
    QualityInfo& quality, const double refXYZ[3]);
void write_2_file_RTK(FILE* fp, PosRes& rovPosRes, PosRes& basPosRes, DDCObs& ddObs, SDEpochObs& sdObs,
    QualityInfo& quality, const double refXYZ[3]);
