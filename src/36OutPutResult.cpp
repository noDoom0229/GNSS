#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"

#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cmath>

// ===================== 内部私有工具 =====================
namespace {
    const double TRUE_XYZ[3] = { -2267810.173, 5009324.109, 3221016.632 };
}

// ===================== 对外接口实现 =====================

/*-----------------------------------------------
    事后定位卫星与接收机的位置输出
    功能：将每一历元 SPP 定位、SPV 测速结果格式化，追加写入文本结果文件；
    同时完成两类坐标转换：XYZ → BLH，XYZ → ENU
------------------------------------------------*/
void OutputCombinedResult(const EPOCHOBS* Obs, const POSRES* pos, const PPRESULT* Result, const char* outputFile)
{
    std::ofstream fout_PVT(outputFile, std::ios::app);

    if (!fout_PVT.is_open()) {
        printf("输出文件打开失败！\n");
        return;
    }

    GEOCOOR blh;
    XYZ_2_BLH(Result->Position, blh.Blh, R_WGS84, F_WGS84);

    double dENU[3] = { 0.0 };
    CompEnudPos(Result->Position, TRUE_XYZ, &blh, dENU);
    //format: Week SecOfWeek X Y Z B L H Vx Vy Vz VN VE VU dN dE dU PDOP RcvClkOft_GPS RcvClkOft_BDS GPSSatNum BDSSatNum AllSatNum
    fout_PVT << "--------------------------------------------------------" << std::endl;
    fout_PVT << "Epoch Time: Week " << std::setw(4) << Obs->Time.Week                                //历元头部
        << " Sec " << std::fixed << std::setprecision(3) << Obs->Time.SecOfWeek << std::endl;   //打印 GPS 时：周号 + 周内秒，所有结果按同一时间戳归类。
    fout_PVT << "--------------------------------------------------------" << std::endl;

    fout_PVT << "SPP Position:" << std::endl;
    // ECEF XYZ
    fout_PVT << "  X: " << std::fixed << std::setprecision(4) << std::setw(15) << Result->Position[0]
        << "  Y: " << std::setw(15) << Result->Position[1]
        << "  Z: " << std::setw(15) << Result->Position[2] << std::endl;
    // BLH 纬度/经度/大地高
    fout_PVT << "  B: " << std::fixed << std::setprecision(4) << std::setw(15) << blh.Blh[0]
        << "  L: " << std::setw(15) << blh.Blh[1]
        << "  H: " << std::setw(15) << blh.Blh[2] << std::endl;
    // XYZ To ENU
    double B = blh.Blh[0] * Rad;
    double L = blh.Blh[1] * Rad;
    double sinB = sin(B), cosB = cos(B);
    double sinL = sin(L), cosL = cos(L);

    double vx = Result->Velocity[0];
    double vy = Result->Velocity[1];
    double vz = Result->Velocity[2];

    double dV_ENU[3];
    dV_ENU[0] = -sinL * vx + cosL * vy;
    dV_ENU[1] = -sinB * cosL * vx - sinB * sinL * vy + cosB * vz;
    dV_ENU[2] = cosB * cosL * vx + cosB * sinL * vy + sinB * vz;

    //format: Vx Vy Vz VN VE VU
    fout_PVT << "SPV Velocity:" << std::endl;

    fout_PVT << "  Vx: " << std::fixed << std::setprecision(4) << std::setw(10) << Result->Velocity[0]
        << " Vy: " << std::setw(10) << Result->Velocity[1]
        << " Vz: " << std::setw(10) << Result->Velocity[2] << std::endl;

    fout_PVT << "  VN: " << std::fixed << std::setprecision(4) << std::setw(10) << dV_ENU[1]
        << " VE: " << std::setw(10) << dV_ENU[0]
        << " VU: " << std::setw(10) << dV_ENU[2] << std::endl;

    fout_PVT << "Position Error (m):" << std::endl;
    fout_PVT << "  dN: " << std::fixed << std::setprecision(4) << std::setw(10) << dENU[1]
        << "  dE: " << std::setw(10) << dENU[0]
        << "  dU: " << std::setw(10) << dENU[2] << std::endl;

    fout_PVT << "Solution Metrics:" << std::endl;
    fout_PVT << "  PDOP: " << std::fixed << std::setprecision(3) << Result->PDOP << std::endl;

    fout_PVT << "  Clock: GPS " << std::setw(8) << Result->RcvClkOft[0]
        << "s, BDS " << std::setw(8) << Result->RcvClkOft[1] << "s" << std::endl;

    fout_PVT << "  Satellites: Tracked " << std::setw(2) << pos->SatNum_tracked
        << ", Used " << std::setw(2) << pos->SatNum_used << std::endl;
    fout_PVT << "--------------------------------------------------------" << std::endl << std::endl;

    fout_PVT.close();
}
/*    ===========example============

--------------------------------------------------------
Epoch Time: Week 2260 Sec 43210.123
--------------------------------------------------------
SPP Position:
  X:   -2267810.1730  Y:    5009324.1090  Z:    3221016.6320
  B:        30.5900  L:       114.3100  H:        23.5600
SPV Velocity:
  Vx:   0.0040 Vy:  -0.0120 Vz:   0.0020
  VN:   0.0010 VE:  -0.0130 VU:   0.0005
Position Error (m):
  dN:   0.3420  dE:   0.5120  dU:   0.8850
Solution Metrics:
  PDOP: 1.842
  Clock: GPS -0.00000123s, BDS -0.00000211s
  Satellites: Tracked 12, Used 9
--------------------------------------------------------*/


/*-----------------------------------------------
    实时数据结果输出
------------------------------------------------*/
void OutputResult_RealTime(const PPRESULT* Result, const POSRES* Pos, const EPOCHOBS* Obs)
{
    GEOCOOR blh;
    XYZ_2_BLH(Result->Position, blh.Blh, R_WGS84, F_WGS84);
	// format: SPP: Week SecOfWeek GPSSats BDSSats Sats B L H X Y Z SigmaPos Vx Vy Vz SigmaVel PDOP
    printf("SPP: %4d %9.3f ", Obs->Time.Week, Obs->Time.SecOfWeek);
    printf("GPSSats:%2d BDSSats:%2d Sats:%2d ", Result->GPSSatNum, Result->BDSSatNum, Result->AllSatNum);
    printf("B:%13.8f L:%13.8f H:%8.3f ", blh.Blh[0], blh.Blh[1], blh.Blh[2]);
    printf("X:%13.4f Y:%13.4f Z:%13.4f SigmaPos:%7.4f ", Result->Position[0], Result->Position[1], Result->Position[2], Result->SigmaPos);

    printf("Vx:%+8.4f Vy:%+8.4f Vz:%+8.4f SigmaVel:%7.4f ", Result->Velocity[0], Result->Velocity[1], Result->Velocity[2], Result->SigmaVel);
    printf("PDOP:%7.4f\n", Result->PDOP);
}

/*-----------------------------------------------
    实时数据结果输出（包含dENU和其他指标）
------------------------------------------------*/
void OutputENUResult_RealTime(const EPOCHOBS* Obs, const POSRES* pos, const PPRESULT* Result, const char* outputFile)
{
    std::ofstream fout_PVT(outputFile, std::ios::app);

    if (!fout_PVT.is_open()) {
        printf("ENU误差文件打开失败！\n");
        return;
    }
	//To Blh   
    GEOCOOR refBlh;
    XYZ_2_BLH(TRUE_XYZ, refBlh.Blh, R_WGS84, F_WGS84);
    // ENU error 
    double dENU[3] = { 0.0 };
    CompEnudPos(Result->Position, TRUE_XYZ, &refBlh, dENU);
    //output 格式
    fout_PVT << std::fixed << std::setprecision(5);
	// dENU: Week SecOfWeek dE dN dU PDOP SigmaPos SigmaVel GPS_RcvClkOft BDS_RcvClkOft GPSSats BDSSats Sats
    fout_PVT << "dENU: "
        << Result->Time.Week << " "
        << Result->Time.SecOfWeek << "  "
        << "dE：" << dENU[0] << " "
        << "dN：" << dENU[1] << " "
        << "dU：" << dENU[2] << " "
        << "PDOP：" << Result->PDOP << " "
        << "SigmaPos：" << Result->SigmaPos << " "
        << "SigmaVel：" << Result->SigmaVel << " "
        << "GPS_RcvClkOft：" << Result->RcvClkOft[0] << " "
        << "BDS_RcvClkOft：" << Result->RcvClkOft[1] << " "
        << "GPS_Sats：" << Result->GPSSatNum << " "
        << "BDS_Sats：" << Result->BDSSatNum << " "
        << "Sum_Sats：" << Result->AllSatNum
        << std::endl;

    fout_PVT.close();
}