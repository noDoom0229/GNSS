#pragma once
#include "GnssConsts.h"

/*-----------------------------------------------
    数据结构
------------------------------------------------*/

// ==========================================================
//                     时间系统
// ==========================================================
// 通用年月日时分秒时间
struct COMMONTIME
{
    unsigned short Year;       // 年
    unsigned short Month;      // 月
    unsigned short Day;        // 日
    unsigned short Hour;       // 时
    unsigned short Minute;     // 分
    double Second;             // 秒（支持小数）
    COMMONTIME() : Year(0), Month(0), Day(0), Hour(0), Minute(0), Second(0.0) {}
};

// 简化儒略日时间
struct MJDTIME
{
    int Days;                  // 儒略日整数天
    double FracDay;            // 天内小数部分
    MJDTIME() : Days(0), FracDay(0.0) {}
};

// GPS时间（周+周内秒）
struct GPSTIME
{
    unsigned short Week;       // GPS周数
    double SecOfWeek;          // 周内秒（支持小数）
    GPSTIME() : Week(0), SecOfWeek(0.0) {}
};

// ==========================================================
//                     坐标系统
// ==========================================================
// 地心直角坐标 XYZ
union XYZ
{
    struct {
        double x;              // X 分量（ECEF）
        double y;              // Y 分量
        double z;              // Z 分量
    };
    double xyz[3];             // 数组形式，方便遍历
};

// 大地经纬度高程 BLH
union GEOCOOR
{
    struct {
        double latitude;       // 纬度（rad）
        double longitude;      // 经度（rad）
        double height;         // 高程（m）
    };
    double Blh[3];             // 数组形式
};

// 站心地平 NEU 坐标
union NEU
{
    struct {
        double dE;             // 东分量
        double dN;             // 北分量
        double dU;             // 天分量
    };
    double Neu[3];             // 数组形式
};

// ==========================================================
//                     观测值
// ==========================================================
// 单颗卫星原始观测数据
struct SATOBSDATA
{
    short Prn;                  // 卫星PRN编号
    GNSSSys System;            // 卫星系统类型（GPS/BDS/GLONASS等）
    double P[2];                // 伪距观测值（m），双频
    double L[2];                // 载波相位观测值，双频。注意：解码时已乘波长，单位是【米】不是周
    double D[2];                // 多普勒观测值（已乘波长，m/s）
    double cn0[2];              // 载噪比 C/N0 (dB-Hz)
    double LockTime[2];         // 相位锁定时间 (s)
    unsigned char half[2];      // 半周标记（= NovAtel Parity known flag，0 表示可能存在半周）
    unsigned char LLI[2];       // RINEX 失锁标记：bit0 = 失锁重捕，bit1 = 半周模糊
    bool Valid;                 // 观测值是否有效

    SATOBSDATA()
    {
        Prn = 0;
        System = UnknownSys;
        for (int i = 0; i < 2; i++)
        {
            P[i] = L[i] = D[i] = cn0[i] = LockTime[i] = 0.0;
            half[i] = LLI[i] = 0;
        }
        Valid = false;
    }
};

// 广播星历结构体（GPS/BDS通用）
struct GPSEPHREC
{
    short PRN;                  // 卫星PRN
    GNSSSys System;            // 卫星系统
    GPSTIME TOC;               // 钟差参考时刻
    GPSTIME TOE;               // 轨道参考时刻
    double ClkBias;            // 卫星钟差
    double ClkDrift;           // 钟漂
    double ClkDriftRate;       // 钟漂率
    double IODE;               // 轨道数据编号
    double IODC;               // 钟差数据编号
    double SqrtA;              // 长半轴平方根
    double M0;                 // 平近点角
    double e;                  // 偏心率
    double OMEGA;              // 升交点赤经
    double i0;                 // 轨道倾角
    double omega;              // 近地点角距
    double Crs, Cuc, Cus;      // 轨道修正项
    double Cic, Cis, Crc;      // 轨道修正项
    double DeltaN;            // 平均角速度修正
    double OMEGADot;          // 升交点赤经变化率
    double iDot;              // 倾角变化率
    int SVHealth;              // 卫星健康状态
    double TGD1, TGD2;        // 硬件群延迟

    GPSEPHREC()
    {
        PRN = 0;
        System = UnknownSys;
        TOC.Week = TOE.Week = 0;
        TOC.SecOfWeek = TOE.SecOfWeek = 0.0;
        ClkBias = ClkDrift = ClkDriftRate = 0.0;
        IODE = IODC = SqrtA = M0 = e = OMEGA = i0 = omega = 0.0;
        Crs = Cuc = Cus = Cic = Cis = Crc = 0.0;
        DeltaN = OMEGADot = iDot = 0.0;
        SVHealth = 0;
        TGD1 = TGD2 = 0.0;
    }
};

// NovAtel原始单点定位输出结果
struct POSRES
{
    GPSTIME Time;              // 定位时间
    double Pos[3];             // 定位坐标 BLH / XYZ
    double SigmaPos;          // 位置精度
    double SigmaVel;          // 速度精度
    double PDOP;               // 几何精度因子
    unsigned char SatNum_tracked;  // 跟踪卫星数
    unsigned char SatNum_used;     // 解算卫星数
    float undulation;         // 大地水准面差距
    float lat_sigma, lon_sigma, hgt_sigma;  // 各方向误差
    unsigned int sol_status;  // 解状态
    unsigned int pos_type;    // 定位类型

    POSRES()
    {
        Pos[0] = Pos[1] = Pos[2] = 0.0;
        SigmaPos = SigmaVel = PDOP = 0.0;
        SatNum_tracked = SatNum_used = 0;
        undulation = lat_sigma = lon_sigma = hgt_sigma = 0.0f;
        sol_status = pos_type = 0;
    }
};

// 卫星位置、速度、钟差、改正量打包
struct SATPVT
{
    double SatPos[3];         // 卫星三维位置（ECEF）
    double SatVel[3];         // 卫星三维速度
    double SatClkOft;         // 卫星钟差
    double SatClkSft;         // 卫星钟漂
    double Elevation;         // 高度角（rad）
    double Azimuth;           // 方位角（rad）
    double TropCorr;          // 对流层改正
    double Tgd1, Tgd2;        // 硬件延迟改正
    bool Valid;               // 计算是否有效

    SATPVT()
    {
        SatPos[0] = SatPos[1] = SatPos[2] = 0.0;
        SatVel[0] = SatVel[1] = SatVel[2] = 0.0;
        SatClkOft = SatClkSft = Tgd1 = Tgd2 = TropCorr = 0.0;
        Elevation = PAI / 2.0;
        Azimuth = 0.0;
        Valid = false;
    }
};

// MW、GF线性组合（周跳探测专用）
struct MWGF
{
    short Prn;                // 卫星PRN
    GNSSSys Sys;              // 卫星系统
    double MW, GF, PIF;       // 组合观测值
    int n;                    // 历元计数
    double D_prev;            // 上一历元多普勒
    double P_prev;            // 上一历元伪距
    double L_prev;            // 上一历元相位

    MWGF()
    {
        Prn = n = 0;
        Sys = UnknownSys;
        MW = GF = PIF = D_prev = P_prev = L_prev = 0.0;
    }
};

// 单台接收机完整历元观测包
struct EPOCHOBS
{
    GPSTIME Time;             // 历元时间
    short SatNum;             // 卫星数量
    SATOBSDATA SatObs[MAXCHANNUM];  // 原始观测值
    SATPVT     SatPVT[MAXCHANNUM];  // 卫星PVT
    MWGF       ComObs[MAXCHANNUM];  // 组合观测值
    double Pos[3];            // 接收机位置

    // 上一历元的相位锁定时间，用于判断是否失锁重捕（RINEX LLI 的 bit0）。
    // 观测值数组每历元都会清零，所以锁定时间要单独按 PRN 存一份。
    // 第一维 0 = GPS，1 = BDS；第二维 PRN-1；第三维频率
    double LockPrev[2][MAXBDSNUM][2];

    EPOCHOBS()
    {
        SatNum = 0;
        Pos[0] = Pos[1] = Pos[2] = 0.0;
        for (int s = 0; s < 2; s++)
            for (int i = 0; i < MAXBDSNUM; i++)
                LockPrev[s][i][0] = LockPrev[s][i][1] = -1.0;   // -1 表示还没见过这颗卫星
    }
};

// ==========================================================
//                     差分观测结构体
// ==========================================================
// 单差（基准站-流动站）单卫星观测
struct SDSATOBS
{
    short Prn;                // 卫星PRN
    GNSSSys System;          // 卫星系统
    short Valid;             // 是否有效
    double dP[2], dL[2];     // 单差伪距 / 相位
    short nBas, nRov;        // 基准站 / 流动站索引

    SDSATOBS()
    {
        Prn = nBas = nRov = 0;
        System = UnknownSys;
        dP[0] = dP[1] = dL[0] = dL[1] = 0.0;
        Valid = -1;
    }
};

// 整历元单差数据包
struct SDEPOCHOBS
{
    GPSTIME Time;             // 时间
    short SatNum;             // 卫星数
    SDSATOBS SdSatObs[MAXCHANNUM];  // 单差观测
    MWGF       SdCObs[MAXCHANNUM];  // 单差组合值

    SDEPOCHOBS() : SatNum(0) {}
};

// 双差模糊度固定解算结果
struct DDCOBS
{
    int RefPrn[2], RefPos[2];        // 参考星PRN与位置
    int Sats, DDSatNum[2];           // 双差卫星数
    double FixedAmb[MAXCHANNUM * 4]; // 固定模糊度
    double ResAmb[2], Ratio;         // 模糊度残差 & 检验比
    float FixRMS[2];                 // 固定解RMS
    double dPos[3];                  // 基线向量
    bool bFixed;                     // 是否固定成功

    DDCOBS()
    {
        for (int i = 0; i < 2; i++)
        {
            DDSatNum[i] = 0;
            RefPrn[i] = RefPos[i] = -1;
        }
        Sats = 0;
        dPos[0] = dPos[1] = dPos[2] = 0.0;
        ResAmb[0] = ResAmb[1] = FixRMS[0] = FixRMS[1] = 0.0;
        Ratio = 0.0;
        bFixed = false;
        for (int i = 0; i < MAXCHANNUM * 2; i++)
            FixedAmb[2 * i] = FixedAmb[2 * i + 1] = 0.0;
    }
};

// ==========================================================
//                     单点定位SPP/SPV结果
// ==========================================================
// 解算结果类型（报告 3.1.1 SPP.h）：编号不可改动
enum resultType { None = 0, SPP_SPV = 1, RTKFloat = 2, RTKFixed = 3, BESTpos = 4 };

struct PPRESULT
{
    GPSTIME Time;             // 解算时间
    double Position[3];       // 定位结果 XYZ
    double Velocity[3];       // 测速结果
    double RcvClkOft[2];      // 接收机钟差
    double RcvClkSft;         // 接收机钟漂
    double PDOP, SigmaPos, SigmaVel;  // 精度指标
    short GPSSatNum, BDSSatNum, AllSatNum;  // 卫星数
    bool IsSuccess;           // 解算是否成功
    resultType Type;          // 结果类型（None/SPP_SPV/RTKFloat/RTKFixed）

    PPRESULT()
    {
        Type = None;
        Time.Week = 0; Time.SecOfWeek = 0.0;
        Position[0] = -2267810.173;
        Position[1] = 5009324.109;
        Position[2] = 3221016.632;
        Velocity[0] = Velocity[1] = Velocity[2] = 0.0;
        RcvClkOft[0] = RcvClkOft[1] = 0.0;
        RcvClkSft = 0.0;
        PDOP = SigmaPos = SigmaVel = 999.9;
        GPSSatNum = BDSSatNum = AllSatNum = 0;
        IsSuccess = false;
    }
};

// ==========================================================
//                     全局数据总容器
// ==========================================================
// 整个系统的原始数据、观测、星历、差分结果总管理
struct RAWDAT
{
    EPOCHOBS   BasEpk;        // 基准站历元观测
    EPOCHOBS   RovEpk;        // 流动站历元观测
    SDEPOCHOBS SdObs;         // 单差观测
    DDCOBS     DDObs;         // 双差解算结果
    GPSEPHREC  GpsEph[MAXGPSNUM];  // GPS星历
    GPSEPHREC  BdsEph[MAXBDSNUM];  // BDS星历
};