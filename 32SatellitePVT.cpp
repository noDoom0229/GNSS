#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <cmath>

/*-----------------------------------------------
    卫星导航算法与程序设计 - 6 卫星位置计算V4
------------------------------------------------*/

/*-----------------------------------------------
    判断星历是否过期、计算卫星钟差
    输入：卫星PRN、系统、时间、星历
    输出：卫星钟差存入结构体 Mid->SatClkOft
    返回：true=星历有效，false=星历无效
------------------------------------------------*/

bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATPVT* Mid)
{
    double dt, LT = 7500.0; // GPS时间有效性阈值
    GPSTIME T;
    GPSEPHREC* EPH = nullptr;

    T = *t;// 观测时刻

    // 根据卫星系统类型选择对应星历数据
    if (Sys == GPS)
    {
        EPH = GPSEph + Prn - 1; // GPS星历数组索引
    }
    else if (Sys == BDS)
    {
        EPH = BDSEph + Prn - 1; // BDS星历数组索引
        T.Week -= 1356;         //北斗时 → GPS 时转换     GPST等价时刻=BDT-1356周-14 s
        T.SecOfWeek -= 14;
        LT = 3900.0;            // BDS时间有效性阈值
    }
    else
    {
        return false;           // 非GPS或BDS
    }

    // 计算钟差参考时间差    dt = t_obs + t_TOC
    dt = (T.Week - EPH->TOC.Week) * 604800.0 + T.SecOfWeek - EPH->TOC.SecOfWeek;//前后两者都是GPS周内秒，单位秒
    if (fabs(dt) > LT || EPH->SVHealth != 0)
    {
        return false;           // 星历过期或不健康
    }

	// 计算卫星钟差 GNSS 标准卫星钟差二阶多项式；钟偏、钟漂、钟漂率
    Mid->SatClkOft = EPH->ClkBias + EPH->ClkDrift * dt + EPH->ClkDriftRate * dt * dt;
    return true;
}


/*-----------------------------------------------
    GPS卫星位置、速度、钟差与钟速计算

 GPS 卫星广播星历 + 观测时刻 GPS 时，计算卫星【ECEF 地心地固坐标系位置 XYZ、速度 VxVyVz、卫星钟差、钟速】，
 存入结构体 SATPVT* Mid结果直接供给伪距单点定位 SPP构建观测方程。

 输入：卫星编号 Prn、观测GPS时刻 t = 周+周内秒、卫星星历指针 Eph
 输出：卫星状态结构体，存放位置、速度、钟差、TGD


    DecodeGpsEphem() → GPS星历数组geph[]
            ↓
    SPP单点定位遍历每颗GPS卫星
            ↓
    1. CompSatClkOff()：校验星历时效性、卫星健康
    2. CalculateGPSSatPVT()：计算卫星XYZ位置、速度、精确钟差
            ↓
    构建伪距观测方程：
    ρ_obs = c·(t_r - t_s) + |r_s - r_r| + ε

------------------------------------------------*/
int CalculateGPSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid)
{
    Mid->Valid = true;

    /*  一、计算卫星位置  */
    double A, n0, deltaT, n, Mk, Ek, Et, vk, phi_k, delta_u, delta_r, delta_i,
        u_k, r_k, i_k, Xk, Yk, Omega_k, X, Y, Z;

    // 1、计算轨道长半轴
    A = Eph->SqrtA * Eph->SqrtA;

    // 2、计算平均运动角速度
    n0 = sqrt(GM_WGS84 / (A * A * A));

    // 3、计算相对于星历参考历元的时间         deltaT = t_obs ? t_TOE
    deltaT = (t->Week - Eph->TOE.Week) * 604800.0 + t->SecOfWeek - Eph->TOE.SecOfWeek;

    // 4、对平均运动角速度进行改正 （加上摄动Δn）
    n = n0 + Eph->DeltaN;

    // 5、计算平近点角
    Mk = Eph->M0 + n * deltaT;

    // 6、迭代求解偏近点角 Ek      开普勒方程 M = E ? e sinE
    Ek = 0;
    Et = Mk;// 初始值设为平近点角
    while (fabs(Et - Ek) > 1e-12)
    {
        Et = Ek;
        Ek = Mk + Eph->e * sin(Ek);
    }

    // 7、计算真近点角
    double sin_vk = sqrt(1.0 - Eph->e * Eph->e) * sin(Ek) / (1.0 - Eph->e * cos(Ek));
    double cos_vk = (cos(Ek) - Eph->e) / (1.0 - Eph->e * cos(Ek));
    vk = atan2(sin_vk, cos_vk);// 会自动处理象限

    // 8、计算升交角距
    phi_k = vk + Eph->omega;// 升交角距 = 真近点角 + 近地点角距

    // 9、计算二阶调和改正数
    delta_u = Eph->Cuc * cos(2.0 * phi_k) + Eph->Cus * sin(2.0 * phi_k);
    delta_r = Eph->Crc * cos(2.0 * phi_k) + Eph->Crs * sin(2.0 * phi_k);
    delta_i = Eph->Cic * cos(2.0 * phi_k) + Eph->Cis * sin(2.0 * phi_k);

    // 10-12、 计算改正后的升交角距、轨道半径和轨道倾角
    u_k = phi_k + delta_u;
    r_k = A * (1.0 - Eph->e * cos(Ek)) + delta_r;// 轨道半径 = 长半轴 * （1 - 偏心率 * cos(偏近点角)）+ 轨道修正数
    i_k = Eph->i0 + Eph->iDot * deltaT + delta_i;// 轨道倾角 = 初始倾角 + 倾角变化率 * 时间 + 轨道修正数

    // 13、计算卫星在轨道上面的位置
    Xk = r_k * cos(u_k);
    Yk = r_k * sin(u_k);

    // 14、计算改正后的升交点经度
    Omega_k = Eph->OMEGA + (Eph->OMEGADot - W_WGS84) * deltaT - W_WGS84 * Eph->TOE.SecOfWeek;
    double cosO = cos(Omega_k);
    double sinO = sin(Omega_k);
    double cosi = cos(i_k);
    double sini = sin(i_k);

    // 15、轨道平面 -> ECEF 标量展开
    X = Xk * cosO - Yk * cosi * sinO;
    Y = Xk * sinO + Yk * cosi * cosO;
    Z = Yk * sini;

    /*  二、计算卫星速度  （链式求导） */
    double E_k_dot, phi_k_dot, u_k_dot, r_k_dot, i_k_dot, omega_k_dot,
        Xk_dot, Yk_dot, Vx, Vy, Vz;

    // 偏近点角变化率
    E_k_dot = n / (1.0 - Eph->e * cos(Ek));
    // 真近点角变化率
    phi_k_dot = E_k_dot * sqrt(1.0 - Eph->e * Eph->e) / (1.0 - Eph->e * cos(Ek));
    // 升交角距变化率
    u_k_dot = phi_k_dot + 2.0 * (Eph->Cus * cos(2.0 * phi_k) - Eph->Cuc * sin(2.0 * phi_k)) * phi_k_dot;
    // 轨道半径变化率
    r_k_dot = A * Eph->e * sin(Ek) * E_k_dot + 2.0 * (Eph->Crs * cos(2.0 * phi_k) - Eph->Crc * sin(2.0 * phi_k)) * phi_k_dot;
    // 轨道倾角变化率
    i_k_dot = Eph->iDot + 2.0 * (Eph->Cis * cos(2.0 * phi_k) - Eph->Cic * sin(2.0 * phi_k)) * phi_k_dot;
    // 升交点经度变化率
    omega_k_dot = Eph->OMEGADot - W_WGS84;

    // 轨道平面内速度分量
    // 轨道半径变化率 * cos(升交角距) - 轨道半径 * sin(升交角距) * 升交角距变化率
    Xk_dot = r_k_dot * cos(u_k) - r_k * sin(u_k) * u_k_dot;
    Yk_dot = r_k_dot * sin(u_k) + r_k * cos(u_k) * u_k_dot;

    // 速度分量标量公式
    // 轨道平面内速度分量 * 旋转矩阵分量 + 升交点经度变化率 * 其他项
    Vx = Xk_dot * cosO - Xk * sinO * omega_k_dot
        - Yk_dot * cosi * sinO + Yk * sini * sinO * i_k_dot
        - Yk * cosi * cosO * omega_k_dot;
    Vy = Xk_dot * sinO + Xk * cosO * omega_k_dot
        + Yk_dot * cosi * cosO - Yk * sini * cosO * i_k_dot
        - Yk * cosi * sinO * omega_k_dot;//  omega_k_dot符号与 Vx 相反
    Vz = Yk_dot * sini + Yk * cosi * i_k_dot;

    /*  三、计算卫星钟差与钟速  */
    double T_sv, t_rel, delta_T_sv, delta_t_rel, F;
    // 相对论改正系数
    F = -4.442807633e-10;//GPS相对论常数 IS-GPS-200

    // 计算相对论效应改正数（椭圆轨道带来的附加钟差）
    t_rel = F * Eph->e * Eph->SqrtA * sin(Ek);
    // 计算卫星钟差       总卫星钟差 = 二阶多项式 + 相对论改正
    T_sv = Eph->ClkBias + Eph->ClkDrift * deltaT + Eph->ClkDriftRate * deltaT * deltaT + t_rel;

    // 相对论效应改正数的时间导数
    delta_t_rel = F * Eph->e * Eph->SqrtA * cos(Ek) * E_k_dot;
    // 计算卫星钟速
    delta_T_sv = Eph->ClkDrift + 2.0 * Eph->ClkDriftRate * deltaT + delta_t_rel;

    /*  四、保存结果  */
    Mid->SatPos[0] = X;
    Mid->SatPos[1] = Y;
    Mid->SatPos[2] = Z;

    Mid->SatVel[0] = Vx;
    Mid->SatVel[1] = Vy;
    Mid->SatVel[2] = Vz;

    Mid->SatClkOft = T_sv;      //卫星钟差 (s)
    Mid->SatClkSft = delta_T_sv;//卫星钟速 (s/s)

    Mid->Tgd1 = Eph->TGD1;
    Mid->Tgd2 = Eph->TGD2;

    return 0;
}

/*-----------------------------------------------
    BDS卫星位置、速度、钟差与钟速计算
------------------------------------------------*/
int CalculateBDSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATPVT* Mid)
{
    Mid->Valid = true;

    /*  一、计算卫星位置  */
    double A, n0, deltaT, n, Mk, Ek, Et, vk, phi_k, delta_u, delta_r, delta_i,
        u_k, r_k, i_k, Xk, Yk, X, Y, Z;

    // 1、计算轨道长半轴
    A = Eph->SqrtA * Eph->SqrtA;

    // 2、计算平均角速度
    n0 = sqrt(GM_CGCS2000 / (A * A * A));

    // 3、BDS时间系统修正，计算相对于星历参考历元的时间   GPST → BDT时间转换，计算tk = t_obs(BDT) ? t_TOE(BDT
    deltaT = (t->Week - 1356 - Eph->TOE.Week) * 604800.0 + t->SecOfWeek - 14.0 - Eph->TOE.SecOfWeek;//BDT=GPST?1356周?14秒

    if (deltaT > 302400.0) deltaT -= 604800.0;// 周折叠，把时间差约束在 ±0.5周以内（±302400s）
    if (deltaT < -302400.0) deltaT += 604800.0;

    // 4、对平均运动角速度进行改正
    n = n0 + Eph->DeltaN;

    // 5、计算平近点角
    Mk = Eph->M0 + n * deltaT;

    // 6、迭代解偏近点角
    Ek = 0;
    Et = Mk;
    while (fabs(Et - Ek) > 1e-12)
    {
        Et = Ek;
        Ek = Mk + Eph->e * sin(Ek);
    }

    // 7、计算真近点角
    double sin_vk = sqrt(1.0 - Eph->e * Eph->e) * sin(Ek) / (1.0 - Eph->e * cos(Ek));
    double cos_vk = (cos(Ek) - Eph->e) / (1.0 - Eph->e * cos(Ek));
    vk = atan2(sin_vk, cos_vk);

    // 8、计算升交角距
    phi_k = vk + Eph->omega;

    // 9、轨道摄动改正
    delta_u = Eph->Cuc * cos(2.0 * phi_k) + Eph->Cus * sin(2.0 * phi_k);
    delta_r = Eph->Crc * cos(2.0 * phi_k) + Eph->Crs * sin(2.0 * phi_k);
    delta_i = Eph->Cic * cos(2.0 * phi_k) + Eph->Cis * sin(2.0 * phi_k);

    // 10-12、计算改正后的升交角距、轨道半径和轨道倾角
    u_k = phi_k + delta_u;
    r_k = A * (1.0 - Eph->e * cos(Ek)) + delta_r;
    i_k = Eph->i0 + Eph->iDot * deltaT + delta_i;

    // 13、计算卫星在轨道上面的位置
    Xk = r_k * cos(u_k);
    Yk = r_k * sin(u_k);

    // BDS GEO 判定（倾角小于30°）
    bool isBDS_GEO = (Eph->i0 * 180.0 / PAI < 30.0);//Eph->i0 * 180.0 / PAI是弧度转角度
    double Vx, Vy, Vz;
    double E_k_dot = n / (1.0 - Eph->e * cos(Ek));

    // GEO卫星
    if (isBDS_GEO)
    {
        double Omega_k, cosO, sinO, cosi, sini, Xg, Yg, Zg;
        // 14、计算改正后的升交点经度【不带地球自转修正项】
        Omega_k = Eph->OMEGA + Eph->OMEGADot * deltaT - W_CGCS2000 * Eph->TOE.SecOfWeek;//只有OMEGADot，没有减去地球自转角速度；
        cosO = cos(Omega_k);
        sinO = sin(Omega_k);
        cosi = cos(i_k);
        sini = sin(i_k);

        // 1. 轨道坐标系 -> 临时地固系
        Xg = Xk * cosO - Yk * cosi * sinO;
        Yg = Xk * sinO + Yk * cosi * cosO;
        Zg = Yk * sini;

        // 2. GEO  Rx(5°) 倾斜改正
        const double tilt = 5.0 * Rad;//tilt改正数
        double X1, Y1, Z1;
        RotX(Xg, Yg, Zg, tilt, X1, Y1, Z1);

        // 3. 绕Z轴旋转，地球自转 Rz(-ωe*ΔT)
        double theta = -W_CGCS2000 * deltaT;
        RotZ(X1, Y1, Z1, theta, X, Y, Z);

        /*  二、GEO卫星速度计算  */
        double phi_k_dot, u_k_dot, r_k_dot, i_k_dot, Omega_k_dot,
            Xk_dot, Yk_dot, Xg_dot, Yg_dot, Zg_dot, X1d, Y1d, Z1d, X2d, Y2d, Z2d, rot_dx, rot_dy, rot_dz;

        // 各参数变化率
        phi_k_dot = E_k_dot * sqrt(1.0 - Eph->e * Eph->e) / (1.0 - Eph->e * cos(Ek));                           //真近点角变化率
        u_k_dot = phi_k_dot + 2.0 * (Eph->Cus * cos(2.0 * phi_k) - Eph->Cuc * sin(2.0 * phi_k)) * phi_k_dot;    //升交角距变化率
        r_k_dot = A * Eph->e * sin(Ek) * E_k_dot + 2.0 * (Eph->Crs * cos(2.0 * phi_k)                           //轨道半径变化率
            - Eph->Crc * sin(2.0 * phi_k)) * phi_k_dot;
        i_k_dot = Eph->iDot + 2.0 * (Eph->Cis * cos(2.0 * phi_k) - Eph->Cic * sin(2.0 * phi_k)) * phi_k_dot;    //轨道倾角变化率
        Omega_k_dot = Eph->OMEGADot;                                                                            //升交点经度变化率

        // 轨道平面内速度
        Xk_dot = r_k_dot * cos(u_k) - r_k * sin(u_k) * u_k_dot;
        Yk_dot = r_k_dot * sin(u_k) + r_k * cos(u_k) * u_k_dot;

        // 轨道系速度
        Xg_dot = Xk_dot * cosO - Xk * sinO * Omega_k_dot
            - Yk_dot * cosi * sinO + Yk * sini * sinO * i_k_dot - Yk * cosi * cosO * Omega_k_dot;
        Yg_dot = Xk_dot * sinO + Xk * cosO * Omega_k_dot
            + Yk_dot * cosi * cosO - Yk * sini * cosO * i_k_dot - Yk * cosi * sinO * Omega_k_dot;
        Zg_dot = Yk_dot * sini + Yk * cosi * i_k_dot;

        // Rx 变换速度
        RotX(Xg_dot, Yg_dot, Zg_dot, tilt, X1d, Y1d, Z1d);

        // Rz 变换速度
        RotZ(X1d, Y1d, Z1d, theta, X2d, Y2d, Z2d);

        // Rz 牵连速度（旋转矩阵导数）
        RotZ_Dot(X1, Y1, Z1, theta, -W_CGCS2000, rot_dx, rot_dy, rot_dz);

        // 总速度 = 旋转后速度 + 牵连速度
        Vx = X2d + rot_dx;
        Vy = Y2d + rot_dy;
        Vz = Z2d + rot_dz;
    }
    // 非GEO卫星，同GPS计算
    else
    {
        double Omega_k, cosO, sinO, cosi, sini;
        // 14、计算改正后的升交点经度
        Omega_k = Eph->OMEGA + (Eph->OMEGADot - W_CGCS2000) * deltaT - W_CGCS2000 * Eph->TOE.SecOfWeek;
        cosO = cos(Omega_k);
        sinO = sin(Omega_k);
        cosi = cos(i_k);
        sini = sin(i_k);

        // 15、计算在地固坐标系下的位置
        X = Xk * cosO - Yk * cosi * sinO;
        Y = Xk * sinO + Yk * cosi * cosO;
        Z = Yk * sini;

        /*  二、计算卫星速度  */
        double phi_k_dot, u_k_dot, r_k_dot, i_k_dot, omega_k_dot,
            Xk_dot, Yk_dot;

        // 各参数变化率
        phi_k_dot = E_k_dot * sqrt(1.0 - Eph->e * Eph->e) / (1.0 - Eph->e * cos(Ek));
        u_k_dot = phi_k_dot + 2.0 * (Eph->Cus * cos(2.0 * phi_k) - Eph->Cuc * sin(2.0 * phi_k)) * phi_k_dot;
        r_k_dot = A * Eph->e * sin(Ek) * E_k_dot + 2.0 * (Eph->Crs * cos(2.0 * phi_k) - Eph->Crc * sin(2.0 * phi_k)) * phi_k_dot;
        i_k_dot = Eph->iDot + 2.0 * (Eph->Cis * cos(2.0 * phi_k) - Eph->Cic * sin(2.0 * phi_k)) * phi_k_dot;
        omega_k_dot = Eph->OMEGADot - W_CGCS2000;

        // 轨道平面内速度
        Xk_dot = r_k_dot * cos(u_k) - r_k * sin(u_k) * u_k_dot;
        Yk_dot = r_k_dot * sin(u_k) + r_k * cos(u_k) * u_k_dot;

        // 速度分量标量公式
        Vx = Xk_dot * cosO - Xk * sinO * omega_k_dot
            - Yk_dot * cosi * sinO + Yk * sini * sinO * i_k_dot - Yk * cosi * cosO * omega_k_dot;
        Vy = Xk_dot * sinO + Xk * cosO * omega_k_dot
            + Yk_dot * cosi * cosO - Yk * sini * cosO * i_k_dot - Yk * cosi * sinO * omega_k_dot;
        Vz = Yk_dot * sini + Yk * cosi * i_k_dot;
    }

    /*  三、BDS卫星钟差与钟速计算  */
    double T_sv, t_rel, delta_T_sv, delta_t_rel, F;
    // 相对论改正系数
    F = -4.442807633e-10;

    // 相对论效应改正数
    t_rel = F * Eph->e * Eph->SqrtA * sin(Ek);
    // 卫星钟差
    T_sv = Eph->ClkBias + Eph->ClkDrift * deltaT + Eph->ClkDriftRate * deltaT * deltaT + t_rel;

    // 相对论效应改正数导数
    delta_t_rel = F * Eph->e * Eph->SqrtA * cos(Ek) * E_k_dot;
    // 卫星钟速
    delta_T_sv = Eph->ClkDrift + 2.0 * Eph->ClkDriftRate * deltaT + delta_t_rel;

    /*  四、保存结果  */
    Mid->SatPos[0] = X;
    Mid->SatPos[1] = Y;
    Mid->SatPos[2] = Z;

    Mid->SatVel[0] = Vx;
    Mid->SatVel[1] = Vy;
    Mid->SatVel[2] = Vz;

    Mid->SatClkOft = T_sv;
    Mid->SatClkSft = delta_T_sv;

    Mid->Tgd1 = Eph->TGD1;
    Mid->Tgd2 = Eph->TGD2;

    return 0;
}


/*-----------------------------------------------
    信号发射时刻卫星位置的计算
接收机观测时刻是【信号接收时刻 t_r?】；卫星轨道函数必须使用【信号发射时刻 t_tr】计算卫星位置

流程：
1.迭代求解卫星信号发射时刻与卫星钟差；
2.使用发射时刻计算 GPS/BDS 卫星 ECEF 位置、速度；
3.地球自转效应改正（信号传播期间地球转动，卫星坐标需要旋转修正）；
4.计算卫星高度角、方位角；
5.Hopfield 模型对流层延迟；所有结果存入 Epk->SatPVT[i]，供给伪距单点定位构观测方程。


传入：EPOCHOBS* Epk ： 前历元观测数据（所有卫星伪距、系统、PRN、接收时刻）
    GPSEPHREC* GPSEph	GPSEPHREC* BDSEph ：卫星星历数组首地址（GPS、BDS）
    double UserPos[3]	用户位置 ECEF XYZ
------------------------------------------------*/
void ComputeSatPVTAtSignalTrans(EPOCHOBS* Epk, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, double UserPos[3])
{
    GPSTIME T = Epk->Time;	// 定义当前历元的接收机接收GPS时间T

    // 循环对每个星历求卫星位置
    for (int i = 0; i < Epk->SatNum; i++)
    {
        GPSTIME t_AtSig;			// 定义卫星发射时刻t_tr
        double T_Threshold = 0.0;	// 设置迭代计算卫星钟差的阈值

        // 卫星钟差初步计算，失败则标记无效
        if (!CompSatClkOff(Epk->SatObs[i].Prn, Epk->SatObs[i].System, &Epk->Time, GPSEph, BDSEph, Epk->SatPVT + i))
        {
            Epk->SatPVT[i].Valid = false;
            continue;	// 星历过期则跳过计算，并标记为false
        }

        // 当系统是GPS系统时
        if (Epk->SatObs[i].System == GPS)
        {
            GPSEPHREC* eph = GPSEph + Epk->SatObs[i].Prn - 1;	// 取本组观测值对应的卫星的星历
            Epk->SatPVT[i].SatClkOft = 0.0;						// 初始钟差设置为零

            // 迭代计算卫星钟差
            do
            {
                // 计算卫星发射时刻 t_tr = t_r - P/c - delta*t_j
                t_AtSig.Week = T.Week;
                t_AtSig.SecOfWeek = T.SecOfWeek - Epk->SatObs[i].P[0] / C_Light - Epk->SatPVT[i].SatClkOft;

                // 计算卫星钟差 delta_t = a0 + a1*(t-toc) + a2*(t-toc)^2
                double t_toc = (t_AtSig.Week - eph->TOC.Week) * 604800.0 + t_AtSig.SecOfWeek - eph->TOC.SecOfWeek;
                double dt = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;

                T_Threshold = fabs(dt - Epk->SatPVT[i].SatClkOft);	// 更新迭代判断阈值
                Epk->SatPVT[i].SatClkOft = dt;						// 更新卫星钟差

            } while (T_Threshold > 1e-12);	// 收敛判断

            // 计算GPS卫星位置、速度（含相对论效应）
            CalculateGPSSatPVT(Epk->SatObs[i].Prn, &t_AtSig, eph, Epk->SatPVT + i);

            // 计算信号传播时间
            double Dx = Epk->SatPVT[i].SatPos[0] - UserPos[0];
            double Dy = Epk->SatPVT[i].SatPos[1] - UserPos[1];
            double Dz = Epk->SatPVT[i].SatPos[2] - UserPos[2];
            double t_trans = sqrt(Dx * Dx + Dy * Dy + Dz * Dz) / C_Light;
            double alpha = W_WGS84 * t_trans;	// 地球自转角度

            // 卫星位置、速度原始值
            double px = Epk->SatPVT[i].SatPos[0];
            double py = Epk->SatPVT[i].SatPos[1];
            double pz = Epk->SatPVT[i].SatPos[2];
            double vx = Epk->SatPVT[i].SatVel[0];
            double vy = Epk->SatPVT[i].SatVel[1];
            double vz = Epk->SatPVT[i].SatVel[2];

            double p_new[3], v_new[3];
            RotZScalar(px, py, pz, alpha, p_new[0], p_new[1], p_new[2]);	// 地球自转改正（位置）
            RotZScalar(vx, vy, vz, alpha, v_new[0], v_new[1], v_new[2]);	// 地球自转改正（速度）

            // 更新改正后的卫星位置、速度
            Epk->SatPVT[i].SatPos[0] = p_new[0];
            Epk->SatPVT[i].SatPos[1] = p_new[1];
            Epk->SatPVT[i].SatPos[2] = p_new[2];
            Epk->SatPVT[i].SatVel[0] = v_new[0];
            Epk->SatPVT[i].SatVel[1] = v_new[1];
            Epk->SatPVT[i].SatVel[2] = v_new[2];

            // 计算卫星高度角、方位角
            GEOCOOR blh;
            XYZ_2_BLH(UserPos, blh.Blh, R_WGS84, F_WGS84);
            //根据接收机位置 + 改正后卫星坐标计算高度角 Elevation、方位角 Azimuth（弧度）；
            CompSatElAz(UserPos, Epk->SatPVT[i].SatPos, &blh, &Epk->SatPVT[i].Elevation, &Epk->SatPVT[i].Azimuth);

            // 对流层延迟改正Hopfield 模型
            Epk->SatPVT[i].TropCorr = Hopfield(blh.Blh[2], Epk->SatPVT[i].Elevation);
        }
        // 当系统是BDS系统时
        else if (Epk->SatObs[i].System == BDS)
        {
            GPSEPHREC* eph = BDSEph + Epk->SatObs[i].Prn - 1;	// 取本组观测值对应的卫星的星历
            GPSTIME t_BDS;										// BDS时间
            Epk->SatPVT[i].SatClkOft = 0.0;						// 初始钟差设置为零

            // 迭代计算卫星钟差
            do
            {
                // 计算卫星发射时刻
                t_AtSig.Week = T.Week;
                t_AtSig.SecOfWeek = T.SecOfWeek - Epk->SatObs[i].P[0] / C_Light - Epk->SatPVT[i].SatClkOft;
                t_BDS.Week = t_AtSig.Week - 1356;				// BDT与GPST周数差1356
                t_BDS.SecOfWeek = t_AtSig.SecOfWeek - 14.0;		// BDT与GPST秒数差14s

                // 计算卫星钟差
                double t_toc = (t_BDS.Week - eph->TOC.Week) * 604800.0 + t_BDS.SecOfWeek - eph->TOC.SecOfWeek;
                double dt = eph->ClkBias + eph->ClkDrift * t_toc + eph->ClkDriftRate * t_toc * t_toc;

                T_Threshold = fabs(dt - Epk->SatPVT[i].SatClkOft);	// 更新迭代判断阈值
                Epk->SatPVT[i].SatClkOft = dt;						// 更新卫星钟差

            } while (T_Threshold > 1e-12);	// 收敛判断

            // 计算BDS卫星位置、速度
            CalculateBDSSatPVT(Epk->SatObs[i].Prn, &t_AtSig, eph, Epk->SatPVT + i);

            // 计算信号传播时间
            double Dx = Epk->SatPVT[i].SatPos[0] - UserPos[0];
            double Dy = Epk->SatPVT[i].SatPos[1] - UserPos[1];
            double Dz = Epk->SatPVT[i].SatPos[2] - UserPos[2];
            double t_trans = sqrt(Dx * Dx + Dy * Dy + Dz * Dz) / C_Light;//  距离 / 光速
            double alpha = W_CGCS2000 * t_trans;	// 地球自转角度

            // 卫星位置、速度原始值
            double px = Epk->SatPVT[i].SatPos[0];
            double py = Epk->SatPVT[i].SatPos[1];
            double pz = Epk->SatPVT[i].SatPos[2];
            double vx = Epk->SatPVT[i].SatVel[0];
            double vy = Epk->SatPVT[i].SatVel[1];
            double vz = Epk->SatPVT[i].SatVel[2];

            double p_new[3], v_new[3];
            RotZScalar(px, py, pz, alpha, p_new[0], p_new[1], p_new[2]);	// 地球自转改正（位置）
            RotZScalar(vx, vy, vz, alpha, v_new[0], v_new[1], v_new[2]);	// 地球自转改正（速度）

            // 更新改正后的卫星位置、速度
            Epk->SatPVT[i].SatPos[0] = p_new[0];
            Epk->SatPVT[i].SatPos[1] = p_new[1];
            Epk->SatPVT[i].SatPos[2] = p_new[2];
            Epk->SatPVT[i].SatVel[0] = v_new[0];
            Epk->SatPVT[i].SatVel[1] = v_new[1];
            Epk->SatPVT[i].SatVel[2] = v_new[2];

            // 计算卫星高度角、方位角
            GEOCOOR blh;
            XYZ_2_BLH(UserPos, blh.Blh, R_CGCS2000, F_CGCS2000);
            CompSatElAz(UserPos, Epk->SatPVT[i].SatPos, &blh, &Epk->SatPVT[i].Elevation, &Epk->SatPVT[i].Azimuth);

            // 对流层延迟改正
            Epk->SatPVT[i].TropCorr = Hopfield(blh.Blh[2], Epk->SatPVT[i].Elevation);
        }
    }
}

