/*-----------------------------------------------------------------------------
    writeToFile.cpp  —— RTK 结果输出（对应报告 3.2.6 writeToFile.cpp，表 3）

    第 1~14 列严格按报告表 3：
      WEEK SOW X Y Z dX dY dZ Ratio SatNum Type dE dN dU
      X Y Z    流动站 WGS84 ECEF 坐标
      dX dY dZ 基线向量（流动站 - 基准站）
      SatNum   单差卫星数
      Type     结果类型 1=SPP 2=RTK Float 3=RTK Fixed
      dE dN dU 流动站相对参考坐标的 ENU 误差；没有参考坐标时输出 nan

    第 15~20 列为课件 II-2 P19「精度评定与结果输出模块」要求的精度信息：
      Sigma0 RMS PDOP mE mN mU
      （只有 RTK 解算成功的历元才有意义，SPP 历元输出 0）
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cmath>
#include <limits>

// 计算 ENU 误差（流动站 - 参考坐标）
static void calc_enu_error(PosRes& rovPosRes, const double refXYZ[3], double dENU[3])
{
    if (fabs(refXYZ[0]) < 1.0 && fabs(refXYZ[1]) < 1.0 && fabs(refXYZ[2]) < 1.0)
    {
        double nan = numeric_limits<double>::quiet_NaN();
        dENU[0] = dENU[1] = dENU[2] = nan;
        return;
    }
    GEOCOOR blh;
    XYZ_2_BLH(refXYZ, blh.Blh, R_WGS84, F_WGS84);
    CompEnudPos(rovPosRes.Position, refXYZ, &blh, dENU);   // dENU[0]=E [1]=N [2]=U
}

// 生成一行结果文本
static void format_line(char* buf, int size, PosRes& rovPosRes, PosRes& basPosRes,
    DDCObs& ddObs, SDEpochObs& sdObs, QualityInfo& quality, const double refXYZ[3])
{
    double dXYZ[3];
    for (int k = 0; k < 3; k++) dXYZ[k] = rovPosRes.Position[k] - basPosRes.Position[k];

    double dENU[3];
    calc_enu_error(rovPosRes, refXYZ, dENU);

    snprintf(buf, size,
        "%4d %10.3f %14.4f %14.4f %14.4f %10.4f %10.4f %10.4f %8.2f %3d %2d %9.4f %9.4f %9.4f"
        " %8.4f %8.4f %6.2f %8.4f %8.4f %8.4f",
        rovPosRes.Time.Week, rovPosRes.Time.SecOfWeek,
        rovPosRes.Position[0], rovPosRes.Position[1], rovPosRes.Position[2],
        dXYZ[0], dXYZ[1], dXYZ[2],
        ddObs.Ratio, (int)sdObs.SatNum, (int)rovPosRes.Type,
        dENU[0], dENU[1], dENU[2],
        quality.Sigma0, quality.RMS, rovPosRes.PDOP,
        quality.mENU[0], quality.mENU[1], quality.mENU[2]);
}

void write_header_RTK(FILE* fp)
{
    const char* header =
        "%WEEK        SOW              X              Y              Z         dX         dY         dZ"
        "    Ratio Sat Ty        dE        dN        dU   Sigma0      RMS  PDOP       mE       mN       mU";
    if (fp != NULL)
    {
        fprintf(fp, "%s\n", header);
        fflush(fp);
    }
    printf("%s\n", header);
}

void write_2_screen_RTK(PosRes& rovPosRes, PosRes& basPosRes, DDCObs& ddObs, SDEpochObs& sdObs,
    QualityInfo& quality, const double refXYZ[3])
{
    char buf[512];
    format_line(buf, sizeof(buf), rovPosRes, basPosRes, ddObs, sdObs, quality, refXYZ);
    printf("%s\n", buf);
}

void write_2_file_RTK(FILE* fp, PosRes& rovPosRes, PosRes& basPosRes, DDCObs& ddObs, SDEpochObs& sdObs,
    QualityInfo& quality, const double refXYZ[3])
{
    if (fp == NULL) return;
    char buf[512];
    format_line(buf, sizeof(buf), rovPosRes, basPosRes, ddObs, sdObs, quality, refXYZ);
    fprintf(fp, "%s\n", buf);
    fflush(fp);
}
