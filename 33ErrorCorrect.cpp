#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <math.h>
#include <string.h>

/*-----------------------------------------------
	对流层误差改正、电离层延迟改正算法
	卫星导航算法与程序设计-7 误差改正 - v2
------------------------------------------------*/

/*-----------------------------------------------
	对流层误差改正 Hopfield模型

	ppt倒着计算

	输入：
H：接收机大地高（m）
Elev：卫星高度角（弧度）
	输出：
trop_delay：对流层总延迟（单位：米），用于伪距观测改正
------------------------------------------------*/
double Hopfield(const double H, const double Elev)
{
	double trop_delay = 0.0;

	// 合法性判断
	if (fabs(H) > 100000.0 || Elev <= 0.0)
	{
		return trop_delay;
	}

	// 1.气象参数计算
	// 相对湿度（随高度衰减）
	double RH = RelHum * exp(-0.0006396 * (H - SeaLevel));
	// 大气压随高度衰减
	double p = AtmosPre * pow(1.0 - 0.0000226 * (H - SeaLevel), 5.225);
	// 气温：标准大气垂直递减率 6.5K/km
	double T = TempT0 - 0.0065 * (H - SeaLevel);
	// 水汽压 e (mbar)
	double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);

	// 2.湿/干分量高度
	double h_w = 11000.0;
	double h_d = 40136.0 + 148.72 * (TempT0 - 273.16);

	// 3.湿/干延迟系数计算
	double K_w = 155.2e-7 * 4810.0 * e * (h_w - H) / (T * T);
	double K_d = 155.2e-7 * p * (h_d - H) / T;

	// 4.高度角映射
	double arg1 = sqrt(Elev * Elev + 6.25);
	double arg2 = sqrt(Elev * Elev + 2.25);
	double s1 = sin(arg1 * Rad);
	double s2 = sin(arg2 * Rad);

	trop_delay = K_d / s1 + K_w / s2;
	return trop_delay;
}



/*-----------------------------------------------
	粗差探测
	卫星导航算法与程序设计-8 单点定位与测速V5

功能：基于载波相位 GF 组合、MW 组合历元间差分，进行观测值粗差检测（周跳 + 伪距粗差筛查），给单颗卫星观测设置 SatObs[i].Valid 有效标记
------------------------------------------------*/
void DetectOutlier(EPOCHOBS* Obs)
{
	//CurComObs：保存当前历元所有卫星的 GF、MW、IF 伪距组合，临时容器；全部清零。
	MWGF CurComObs[MAXCHANNUM];
	memset(CurComObs, 0, sizeof(CurComObs));

	for (int i = 0; i < Obs->SatNum; i++)
	{
		Obs->SatObs[i].Valid = false;//先默认标记观测无效，校验通过再置 true。

		// 1、只要有一个观测值异常，该卫星直接判定为无效P：伪距 (m)，L：载波相位 (周)
		if (fabs(Obs->SatObs[i].P[0]) < 1e-5 || fabs(Obs->SatObs[i].P[1]) < 1e-5 ||
			fabs(Obs->SatObs[i].L[0]) < 1e-5 || fabs(Obs->SatObs[i].L[1]) < 1e-5)
			continue;

		// 2、计算当前历元该卫星的GF和MW组合值
		double Lgf, Lmw;

		if (Obs->SatObs[i].System == GPS)	// GPS系统
		{
			Lgf = Obs->SatObs[i].L[0] - Obs->SatObs[i].L[1];												//L_GF=L_1-L_2
			Lmw = (FG1_GPS * Obs->SatObs[i].L[0] - FG2_GPS * Obs->SatObs[i].L[1]) / (FG1_GPS - FG2_GPS)		//L_MW=(f1*L1-f2*L2)/(f1-f2)
				- (FG1_GPS * Obs->SatObs[i].P[0] + FG2_GPS * Obs->SatObs[i].P[1]) / (FG1_GPS + FG2_GPS);
		}
		else if (Obs->SatObs[i].System == BDS)	// BDS系统
		{
			Lgf = Obs->SatObs[i].L[0] - Obs->SatObs[i].L[1];												//L_GF=L_1-L_2
			Lmw = (FG1_BDS * Obs->SatObs[i].L[0] - FG3_BDS * Obs->SatObs[i].L[1]) / (FG1_BDS - FG3_BDS)		//L_MW=(f1*L1-f2*L2)/(f1-f2)
				- (FG1_BDS * Obs->SatObs[i].P[0] + FG3_BDS * Obs->SatObs[i].P[1]) / (FG1_BDS + FG3_BDS);
		}
		else
			continue;

		// 存储当前MWGF计算值
		CurComObs[i].Prn = Obs->SatObs[i].Prn;
		CurComObs[i].Sys = Obs->SatObs[i].System;
		CurComObs[i].GF = Lgf;
		CurComObs[i].MW = Lmw;
		CurComObs[i].n = 1;

		// 3、伪距的IF组合观测值   P_IF = (f1^2*P1-f2^2*P2)/(f1^2-f2^2)
		if (Obs->SatObs[i].System == GPS)
			CurComObs[i].PIF = (FG1_GPS * FG1_GPS * Obs->SatObs[i].P[0] - FG2_GPS * FG2_GPS * Obs->SatObs[i].P[1]) / (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS);
		else
			CurComObs[i].PIF = (FG1_BDS * FG1_BDS * Obs->SatObs[i].P[0] - FG3_BDS * FG3_BDS * Obs->SatObs[i].P[1]) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);

		// 与上一历元组合值对比做历元间差分
		for (int j = 0; j < MAXCHANNUM; j++)
		{
			// 从上个历元的MWGF数据中查找该卫星
			if (Obs->ComObs[j].Prn == CurComObs[i].Prn && Obs->ComObs[j].Sys == CurComObs[i].Sys)
			{
				// 4、差值 dGF / dMW
				double dGF = fabs(CurComObs[i].GF - Obs->ComObs[j].GF);
				double dMW = fabs(CurComObs[i].MW - Obs->ComObs[j].MW);

				// 5、检查是否超限，超限标记为粗差
				if (dGF < GFThres && dMW < MWThres)
				{
					Obs->SatObs[i].Valid = true;
					CurComObs[i].MW = (Obs->ComObs[j].MW * (Obs->ComObs[j].n - 1) + CurComObs[i].MW) / Obs->ComObs[j].n;
					CurComObs[i].n = Obs->ComObs[j].n + 1;
				}
				break;
			}
		}
	}

	// 更新历元组合值数组
	memcpy(Obs->ComObs, CurComObs, sizeof(CurComObs));
}
