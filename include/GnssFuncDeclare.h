#pragma once
#include "GnssConsts.h"
#include "GnssStructs.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <windows.h>

using namespace std;


// ==========================================================
//                     矩阵向量运算
// ==========================================================
// 矩阵加法
void matrix_add(const double* A, const double* B, double* C, int rows, int cols);
// 矩阵减法
void matrix_subtract(const double* A, const double* B, double* C, int rows, int cols);
// 矩阵乘法
void matrix_multiply(const double* A, const double* B, double* C, int m, int n, int p);
// 矩阵 × 向量
void matrix_vector_multiply(const double* A, const double* v, double* v_out, int m, int n);
// 矩阵 × 常数
void matrix_num_multiply(const double* A, double num, double* C, int rows, int cols);
// 矩阵转置
void matrix_trans(const double* A, double* C, int m, int n);
// 矩阵初始化为 0
void matrix_initialize(double* mat, int rows, int cols);
// 4x4矩阵求逆
bool matrix_invert4x4(const double m[4][4], double inv[4][4]);
// 5x5矩阵求逆
bool matrix_invert5x5(const double m[5][5], double inv[5][5]);

// 向量初始化为 0
void vector_initialize(double* vec, int size);
// 向量加法
void vector_add(const double* a, const double* b, double* c, int size);
// 向量减法
void vector_subtract(const double* a, const double* b, double* c, int size);
// 向量 × 常数
void vector_num_multiply(const double* a, double num, double* c, int size);
// 向量点乘
double vector_multiply(const double* a, const double* b, int size);
// 3 维向量叉乘
void vector_cross_multiply(const double a[3], const double b[3], double c[3]);
// 向量模长
double vector_magnitude(const double* vec, int size);


// ==========================================================
//                     时间系统转换
// ==========================================================
// 通用时转简化儒略日
void CommonTime_2_MjdTime(const COMMONTIME* ct, MJDTIME* mjd);
// 通用时转GPS时
void CommonTime_2_GpsTime(const COMMONTIME* ct, GPSTIME* gps);
// 简化儒略日转通用时
void MjdTime_2_CommonTime(const MJDTIME* mjd, COMMONTIME* ct);
// 简化儒略日转GPS时
void MjdTime_2_GpsTime(const MJDTIME* mjd, GPSTIME* gps);
// GPS时转通用时
void GpsTime_2_CommonTime(const GPSTIME* gps, COMMONTIME* ct);
// GPS时转简化儒略日
void GpsTime_2_MjdTime(const GPSTIME* gps, MJDTIME* mjd);


// ==========================================================
//                     坐标系统转换
// ==========================================================
// 坐标旋转函数
void RotX(double x, double y, double z, double alpha, double& ox, double& oy, double& oz);
void RotZ(double x, double y, double z, double theta, double& ox, double& oy, double& oz);
void RotZ_Dot(double x, double y, double z, double theta, double omega, double& ox, double& oy, double& oz);
// Z轴标量旋转
void RotZScalar(double x, double y, double z, double alpha, double& xo, double& yo, double& zo);

// BLH坐标 ⇌ 地心直角坐标
void BLH_2_XYZ(const double BLH[3], double XYZ[3], const double R, const double F);
void XYZ_2_BLH(const double XYZ[3], double BLH[3], const double R, const double F);

// BLH → 站心NEU坐标
void BLH_2_Neu(const GEOCOOR* stationBlh, double satX, double satY, double satZ,
    double staX, double staY, double staZ, double& N, double& E, double& U);

// 计算位置差对应的NEU偏差
void CompEnudPos(const double X0[], const double Xr[], const GEOCOOR* Blh, double dNeu[]);

// 计算卫星高度角与方位角
void CompSatElAz(const double Xr[], const double Xs[], const GEOCOOR* Blh, double* Elev, double* Azim);


// ==========================================================
//                    网络与串口通信
// ==========================================================
// 打开Socket
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);
// 关闭Socket
void CloseSocket(SOCKET& sock);
// 保存1分钟Socket数据流到文件
bool SaveSocketStreamToFile_1min(const char* ip, unsigned short port, const char* binFilePath, double time);


// ==========================================================
//                     NovAtel 数据解码
// ==========================================================
// NovAtel OEM7 数据解码入口
int  __cdecl DecodeNovOem7Dat(unsigned char* buf, int& len, EPOCHOBS* obs,
    GPSEPHREC* gpsEph, GPSEPHREC* bdsEph, POSRES* pos, int mode);

// 解析RANGE观测值
void DecodeRange(const unsigned char* msg, EPOCHOBS* obs);
// 解析GPS星历
void DecodeGpsEphem(const unsigned char* msg, GPSEPHREC geph[]);
// 解析BDS星历
void DecodeBdsEphem(const unsigned char* msg, GPSEPHREC beph[]);
// 解析PSRPOS单点定位结果
void DecodePos(const unsigned char* msg, POSRES* pos);


// ==========================================================
//                   卫星轨道与钟差计算
// ==========================================================
// 计算GPS卫星PVT
int  CalculateGPSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid);
// 计算BDS卫星PVT
int  CalculateBDSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid);
// 计算卫星钟差并检查星历有效性
bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t,
    GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATPVT* Mid);


// ==========================================================
//                    误差改正与质量控制
// ==========================================================
// Hopfield对流层延迟
double Hopfield(const double H, const double Elev);
// 观测值粗差检测
void DetectOutlier(EPOCHOBS* Obs);


// ==========================================================
//                    核心解算接口
// ==========================================================
// 信号发射时刻卫星PVT计算
void ComputeSatPVTAtSignalTrans(EPOCHOBS* Epk, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, double UserPos[3]);
// 单点定位解算
bool SPP(EPOCHOBS* Epoch, RAWDAT* Raw, PPRESULT* Result);
// 单点测速解算
void SPV(EPOCHOBS* Epoch, PPRESULT* Result);


// ==========================================================
//                      结果输出
// ==========================================================
// 输出综合定位结果到文件
void OutputCombinedResult(const EPOCHOBS* Obs, const POSRES* pos, const PPRESULT* Result, const char* outputFile);
// 实时输出解算结果到控制台
void OutputResult_RealTime(const PPRESULT* Result, const POSRES* Pos, const EPOCHOBS* Obs);
// 实时输出ENU误差到文件
void OutputENUResult_RealTime(const EPOCHOBS* Obs, const POSRES* pos, const PPRESULT* Result, const char* outputFile);


