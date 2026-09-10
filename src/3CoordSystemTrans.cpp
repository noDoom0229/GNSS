#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <math.h>
/*-----------------------------------------------
大地坐标(BLH) → 直角坐标(XYZ)
直角坐标(XYZ) → 大地坐标(BLH)
计算 NEU 旋转矩阵（定位误差用）
计算卫星高度角、方位角
------------------------------------------------*/

// 三维坐标绕 X 轴旋转
/**
 * @brief 三维坐标绕 X 轴旋转
 * @param x,y,z 输入坐标
 * @param alpha 旋转角（弧度）
 * @param ox,oy,oz 输出坐标
 */
void RotX(double x, double y, double z, double alpha, double& ox, double& oy, double& oz)
{
    double c = cos(alpha);
    double s = sin(alpha);
    ox = x;
    oy = y * c - z * s;
    oz = y * s + z * c;
}

// 三维坐标绕 Z 轴旋转
void RotZ(double x, double y, double z, double theta, double& ox, double& oy, double& oz)
{
    double c = cos(theta);
    double s = sin(theta);
    ox = x * c - y * s;
    oy = x * s + y * c;
    oz = z;
}

// Z 轴旋转速度牵连项（地球自转改正） 旋转角theta，地球自转角速度omega
void RotZ_Dot(double x, double y, double z, double theta, double omega, double& ox, double& oy, double& oz)
{
    double c = cos(theta);
    double s = sin(theta);
    ox = omega * (-x * s - y * c);// 自转牵连速度 X,Y 分量
    oy = omega * (x * c - y * s);
    oz = 0.0;
}

// Z轴标量旋转,另一种旋转顺序，用于坐标变换 旋转角alpha
void RotZScalar(double x, double y, double z, double alpha, double& xo, double& yo, double& zo)
{
    double c = cos(alpha);
    double s = sin(alpha);
    xo = c * x + s * y;
    yo = -s * x + c * y;
    zo = z;
}



/*-----------------------------------------------
    大地坐标 -> 笛卡尔坐标 (BLH -> XYZ)
------------------------------------------------*/
void BLH_2_XYZ(const double BLH[3], double XYZ[3], const double R, const double F)
{
    double a = R;                           // 长半轴
    double f = F;                           // 扁率
    double b = a * (1.0 - f);               // 短半轴
    double e2 = (a * a - b * b) / (a * a);  // 第一偏心率平方

    double L = BLH[0] * Rad;                // 经度 (弧度)
    double B = BLH[1] * Rad;                // 纬度 (弧度)
    double H = BLH[2];                      // 高程 (米)

    double sinB = sin(B);
    double cosB = cos(B);
    double sinL = sin(L);
    double cosL = cos(L);

    double N = a / sqrt(1.0 - e2 * sinB * sinB); // 卯酉圈半径
    // 计算 XYZ 坐标
    XYZ[0] = (N + H) * cosB * cosL;
    XYZ[1] = (N + H) * cosB * sinL;
    XYZ[2] = (N * (1.0 - e2) + H) * sinB;
}

/*-----------------------------------------------
    笛卡尔坐标 -> 大地坐标 (XYZ -> BLH)
------------------------------------------------*/
void XYZ_2_BLH(const double XYZ[3], double BLH[3], const double R, const double F)
{
    double a = R;
    double f = F;
    double b = a * (1 - f);
    double e2 = (a * a - b * b) / (a * a);
    double x = XYZ[0], y = XYZ[1], z = XYZ[2];
    double p = sqrt(x * x + y * y);  //p = 地面投影半径

    double B = atan2(z, p * (1 - e2)); // 初始纬度 (弧度)
    double H = 0.0;                    //初始高度

    int count = 0;
    double delta, B_new;
    //迭代逼近真实值
    do
    {
        double sinB = sin(B);
        double N = a / sqrt(1.0 - e2 * sinB * sinB);
        H = p / cos(B) - N;
        B_new = atan2(z + e2 * N * sinB, p);
        delta = B - B_new;
        B = B_new;
        count++;
    } while (fabs(delta) > 1e-12 && count < 1000);

    BLH[0] = B * Deg;           // 纬度 度
    BLH[1] = atan2(y, x) * Deg; // 经度 度
    BLH[2] = H;                 // 高程 米

}

/*-----------------------------------------------
    生成 BLH -> NEU 的旋转矩阵
------------------------------------------------*/
void BLH_2_Neu(const GEOCOOR* stationBlh,
    double satX, double satY, double satZ,
    double staX, double staY, double staZ,
    double& N, double& E, double& U)
{
    //经纬度转弧度
    double B = stationBlh->Blh[0] * Rad;
    double L = stationBlh->Blh[1] * Rad;

    double sinB = sin(B);
    double cosB = cos(B);
    double sinL = sin(L);
    double cosL = cos(L);
    // XYZ 坐标差
    double dx = satX - staX;
    double dy = satY - staY;
    double dz = satZ - staZ;

    E = -sinL * dx + cosL * dy;
    N = -sinB * cosL * dx - sinB * sinL * dy + cosB * dz;
    U = cosB * cosL * dx + cosB * sinL * dy + sinB * dz;
}



/*-----------------------------------------------
    计算 ENU 定位误差
------------------------------------------------*/
/**
 * @brief 计算解算坐标与参考真值的 ENU 误差
 * @param X0 解算坐标 XYZ
 * @param Xr 参考真值 XYZ
 * @param Blh 测站BLH
 * @param dNeu 输出 dN dE dU 误差
 */
void CompEnudPos(const double X0[], const double Xr[], const GEOCOOR* Blh, double dNeu[])
{
    double dB = Blh->Blh[0] * Rad;
    double dL = Blh->Blh[1] * Rad;
    double sinB = sin(dB), cosB = cos(dB);
    double sinL = sin(dL), cosL = cos(dL);

    double dx = X0[0] - Xr[0];
    double dy = X0[1] - Xr[1];
    double dz = X0[2] - Xr[2];

    dNeu[0] = -sinL * dx + cosL * dy;
    dNeu[1] = -sinB * cosL * dx - sinB * sinL * dy + cosB * dz;
    dNeu[2] = cosB * cosL * dx + cosB * sinL * dy + sinB * dz;
}




/*-----------------------------------------------
    计算卫星高度角与方位角
------------------------------------------------*/
void CompSatElAz(const double Xr[], const double Xs[], const GEOCOOR* Blh, double* Elev, double* Azim)
{
    double dENU[3];
    CompEnudPos(Xs, Xr, Blh, dENU);
    double EN = sqrt(dENU[0] * dENU[0] + dENU[1] * dENU[1]);
    *Elev = atan2(dENU[2], EN) * Deg;
    *Azim = atan2(dENU[0], dENU[1]) * Deg;  // E在前、N在后
    if (*Azim < 0.0) *Azim += 360.0;
}/*

              卫星
               ↑ dU
               |
       N ↑     |
         |     · 视线
         |    /
         |   /
         |  / E
接收机原点——————→ E
水平投影 EN = √(dE?+dN?)
高度角Elev：视线与水平面夹角
方位角Azim：从正北顺时针转到水平投影线

*/