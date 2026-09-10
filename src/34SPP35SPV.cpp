#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <cmath>
#include <cstring>

/*-----------------------------------------------
* 卫星导航算法与程序设计-8 单点定位与测速V5
SPP 根据卫星位置 + 伪距观测值，解算接收机三维坐标（X/Y/Z）+ 接收机钟差。
SPV 根据卫星速度 + 多普勒观测值，解算接收机三维速度（Vx/Vy/Vz）+ 接收机钟漂。
------------------------------------------------*/

/*-----------------------------------------------
	单点定位SPP计算  伪距观测方程 + 最小二乘迭代

输入：
	Epoch：当前历元观测、卫星 PVT（位置 / 钟差 / 对流层延迟、粗差标记）
	Raw：存放 GPS/BDS 广播星历数组
	Result【输入输出】：
		输入：上一历元定位坐标（作为本次迭代初值）
		输出：最终接收机 XYZ、接收机钟差、PDOP、单位权中误差、参与解算卫星数量

返回：true 迭代收敛解算成功；false 解算失败（卫星不足 / 矩阵奇异 / 迭代发散）。
------------------------------------------------*/

bool SPP(EPOCHOBS* Epoch, RAWDAT* Raw, PPRESULT* Result)
{
	// 1. 保存当前历元时间
	Result->Time = Epoch->Time;

	// 2. 接收机初始位置初始化
	//初始坐标,不能是0！！！发散
	if (fabs(Result->Position[0]) < 1.0 &&
		fabs(Result->Position[1]) < 1.0 &&
		fabs(Result->Position[2]) < 1.0)
	{
		Result->Position[0] = -2267810.173;
		Result->Position[1] = 5009324.109;
		Result->Position[2] = 3221016.632;
	}

	// 3. 迭代初始值
	double X_R0[3] = { 0 };                // 接收机坐标迭代初值
	memcpy(X_R0, Result->Position, sizeof(X_R0));

	double dt_R0[2] = { 0, 0 };            // 接收机钟差 0G1B
	bool flag = true;                      // 迭代控制标志
	int count = 0;                         // 迭代次数计数器

	// 4. 最小二乘迭代定位（非线性方程线性化迭代）
	do {
		// 4.1 计算信号发射时刻的卫星精确位置、钟差、对流层延迟
		// 迭代重算卫星发射时刻PVT，固定接收机位置算有误差    迭代更新：发射时刻 → 卫星坐标 → Sagnac 自转改正 → 对流层延迟。
		ComputeSatPVTAtSignalTrans(Epoch, Raw->GpsEph, Raw->BdsEph, X_R0);

		// 4.2 定义观测方程相关变量
		double B_arr[128][5] = { 0 };      // 设计矩阵 B（行：卫星，列：X/Y/Z/GPS钟差/BDS钟差）最多128颗卫星，5列对应5个待估参数
		double w_arr[128] = { 0 };          // 观测残差向量 w = 观测值 - 理论计算值
		int satnum[2] = { 0, 0 };           // 有效卫星数 [GPS, BDS]   satnum[0]=GPS有效星；satnum[1]=BDS有效星
		int validSat = 0;                   // 总有效卫星数

		// 5. 遍历所有卫星，构建观测方程
		for (int i = 0; i < Epoch->SatNum; i++) {
			// 5.1 跳过无效卫星（观测无效/卫星位置无效）
			if (!Epoch->SatObs[i].Valid || !Epoch->SatPVT[i].Valid)
				continue;

			GNSSSys sys = Epoch->SatObs[i].System;

			// 5.2 根据解算模式过滤卫星（仅GPS / 仅BDS / 联合解）
			if (g_solveMode == ONLY_GPS)
			{
				if (sys != GPS) continue;
			}
			else if (g_solveMode == ONLY_BDS)
			{
				if (sys != BDS) continue;
			}

			// 5.3 计算卫星到接收机的几何距离 ρ
			double Dx = Epoch->SatPVT[i].SatPos[0] - X_R0[0];
			double Dy = Epoch->SatPVT[i].SatPos[1] - X_R0[1];
			double Dz = Epoch->SatPVT[i].SatPos[2] - X_R0[2];
			double Rou = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
			if (Rou < 1e-3) continue;        // 距离异常，跳过

			// 5.4 方向余弦（线性化观测方程系数）
			double cx = -Dx / Rou;//-!!
			double cy = -Dy / Rou;
			double cz = -Dz / Rou;

			double w_i; // 单颗卫星残差

			// 5.5 GPS 卫星观测方程构建
			if (Epoch->SatObs[i].System == GPS) {
				// 设计矩阵 B 赋值
				B_arr[validSat][0] = cx;
				B_arr[validSat][1] = cy;
				B_arr[validSat][2] = cz;
				B_arr[validSat][3] = 1.0;    // GPS钟差系数
				B_arr[validSat][4] = 0.0;    // BDS钟差系数
				satnum[0]++;

				// 伪距残差 = 观测伪距 - （几何距离 + 接收机钟差 - 卫星钟差 + 对流层延迟）
				w_i = Epoch->ComObs[i].PIF - (Rou + dt_R0[0] - C_Light * Epoch->SatPVT[i].SatClkOft + Epoch->SatPVT[i].TropCorr);
			}
			// 5.6 BDS 卫星观测方程构建（含硬件延迟TGD改正）
			else if (Epoch->SatObs[i].System == BDS) {
				// 设计矩阵 B 赋值
				B_arr[validSat][0] = cx;
				B_arr[validSat][1] = cy;
				B_arr[validSat][2] = cz;
				B_arr[validSat][3] = 0.0;    // GPS钟差系数
				B_arr[validSat][4] = 1.0;    // BDS钟差系数
				satnum[1]++;

				// BDS 硬件延迟（群延迟）改正   IF 组合伪距需要引入卫星端 TGD 硬件延迟改正项tgd。
				double tgd = (FG1_BDS * FG1_BDS * C_Light * Epoch->SatPVT[i].Tgd1) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);

				// 伪距残差（含TGD改正）
				w_i = Epoch->ComObs[i].PIF - (Rou + dt_R0[1] - C_Light * Epoch->SatPVT[i].SatClkOft + Epoch->SatPVT[i].TropCorr + tgd);
			}
			else continue;

			// 5.7 保存残差，有效卫星数+1
			w_arr[validSat] = w_i;
			validSat++;
			if (validSat >= 128) break;
		}

		// 6. 解算参数判断：参数=3个坐标 + 可用系统数（GPS/BDS钟差）
		int Para = 3 + (satnum[0] > 0) + (satnum[1] > 0);
		if (validSat < Para)            // 卫星数不足，无法解算
			return false;

		// 7. 构建最小二乘法方程
		double x[5] = { 0 };             // 参数改正量 X/Y/Z/GPS钟差/BDS钟差
		double N[5][5] = { 0 };
		double W[5] = { 0 };
		double N_inv[5][5] = { 0 };      // 法方程逆矩阵
		double PDOP_N_inv[5][5] = { 0 }; // 用于计算PDOP的协方差矩阵

		// 7.1 计算法矩阵 N 和右端项 W
		for (int i = 0; i < validSat; i++) {
			for (int r = 0; r < 5; r++) {
				for (int c = 0; c < 5; c++) {
					N[r][c] += B_arr[i][r] * B_arr[i][c];	// 法方程矩阵 N = B^T * B
				}
				W[r] += B_arr[i][r] * w_arr[i];				// 法方程右端项 W = B^T * w
			}
		}

		bool ok = false;

		// 8. system数量不同，解算参数维数不同，分别求逆解算
		// 8.1 仅GPS解算（4参数：X/Y/Z/GPS钟差）截取N[0~3][0~3]4阶矩阵
		if (satnum[0] > 0 && satnum[1] == 0) {
			double N4[4][4] = { 0 }, W4[4] = { 0 }, inv4[4][4] = { 0 };
			for (int r = 0; r < 4; r++) {
				for (int c = 0; c < 4; c++) N4[r][c] = N[r][c];
				W4[r] = W[r];
			}
			ok = matrix_invert4x4(N4, inv4);
			if (ok) {
				for (int r = 0; r < 4; r++) {
					x[r] = inv4[r][0] * W4[0] + inv4[r][1] * W4[1] + inv4[r][2] * W4[2] + inv4[r][3] * W4[3];
				}
				x[4] = 0;
				memcpy(PDOP_N_inv, inv4, sizeof(inv4));// PDOP计算用
			}
		}
		// 8.2 仅BDS解算（4参数：X/Y/Z/BDS钟差）
		else if (satnum[1] > 0 && satnum[0] == 0) {
			double N4[4][4] = { 0 }, W4[4] = { 0 }, inv4[4][4] = { 0 };
			for (int r = 0; r < 3; r++) {
				for (int c = 0; c < 3; c++) N4[r][c] = N[r][c];
				N4[r][3] = N[r][4];
			}
			// 最后一行：X/Y/Z列来自N的第4行，BDS钟差列来自N的第5行
			for (int c = 0; c < 3; c++) N4[3][c] = N[4][c];
			N4[3][3] = N[4][4];
			W4[0] = W[0]; W4[1] = W[1]; W4[2] = W[2]; W4[3] = W[4];
			// 求逆解算
			ok = matrix_invert4x4(N4, inv4);
			if (ok) {
				x[0] = inv4[0][0] * W4[0] + inv4[0][1] * W4[1] + inv4[0][2] * W4[2] + inv4[0][3] * W4[3];
				x[1] = inv4[1][0] * W4[0] + inv4[1][1] * W4[1] + inv4[1][2] * W4[2] + inv4[1][3] * W4[3];
				x[2] = inv4[2][0] * W4[0] + inv4[2][1] * W4[1] + inv4[2][2] * W4[2] + inv4[2][3] * W4[3];
				x[4] = inv4[3][0] * W4[0] + inv4[3][1] * W4[1] + inv4[3][2] * W4[2] + inv4[3][3] * W4[3];
				x[3] = 0;
				memcpy(PDOP_N_inv, inv4, sizeof(inv4));
			}
		}
		// 8.3 GPS+BDS联合解算（5参数：X/Y/Z/GPS钟差/BDS钟差）
		else {
			ok = matrix_invert5x5(N, N_inv);
			if (ok) {
				for (int r = 0; r < 5; r++) {
					x[r] = 0;
					for (int c = 0; c < 5; c++)
						x[r] += N_inv[r][c] * W[c];
				}
				memcpy(PDOP_N_inv, N_inv, sizeof(N_inv));
			}
		}

		if (!ok) return false;   // 矩阵求逆失败，解算失败

		// 9. 判断迭代收敛（坐标改正量 < 0.0000001 米）最大迭代上限 15 次，防止死循环。
		double norm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		if (norm < 1e-7)
			flag = false;

		// 10. 更新迭代参数：坐标 + 接收机钟差
		X_R0[0] += x[0];
		X_R0[1] += x[1];
		X_R0[2] += x[2];
		dt_R0[0] += x[3];
		dt_R0[1] += x[4];

		// 11. 计算 PDOP 位置精度因子
		Result->PDOP = sqrt(PDOP_N_inv[0][0] + PDOP_N_inv[1][1] + PDOP_N_inv[2][2]);

		// 12. 计算验后单位权中误差（定位精度）
		double vTv = 0;
		for (int i = 0; i < validSat; i++) {
			double res = B_arr[i][0] * x[0] + B_arr[i][1] * x[1] + B_arr[i][2] * x[2]
				+ B_arr[i][3] * x[3] + B_arr[i][4] * x[4] - w_arr[i];
			vTv += res * res;
		}
		Result->SigmaPos = sqrt(vTv / (validSat - Para));// 位置中误差Sigma^=sqrt(vTv/(n-m))，n=观测数，m=参数个数

		// 13. 记录卫星数，迭代次数保护（最多15次）
		count++;
		Result->BDSSatNum = satnum[1];
		Result->GPSSatNum = satnum[0];
		Result->AllSatNum = validSat;
		if (count > 15)
			flag = false;

	} while (flag);    // 迭代循环：未收敛则继续，收敛则退出

	// 14. 保存最终解算结果
	Result->RcvClkOft[0] = dt_R0[0];   // GPS接收机钟差
	Result->RcvClkOft[1] = dt_R0[1];   // BDS接收机钟差
	memcpy(Result->Position, X_R0, sizeof(X_R0)); // 最终定位坐标

	return true;
}


/*-----------------------------------------------
	单点测速SPV计算
	输入：当前历元观测值 + 卫星PVT
	、利用接收机输出的多普勒观测值，先算轨道平面内瞬时速度分量 → 再乘轨道旋转导数矩阵
	采用线性最小二乘直接求解接收机 ECEF 三维速度 + 接收机钟漂（钟差变化率）。
------------------------------------------------*/
void SPV(EPOCHOBS* Epoch, PPRESULT* Result)
{
	// 1. 时间同步：把当前历元时间赋值给定位结果
	Result->Time = Epoch->Time;

	// 2. 定义最小二乘变量
	// 观测系数矩阵 B (行：卫星，列：dx/dy/dz + 接收机钟漂)
	double B[128][4] = { 0 };

	// 观测残差向量 w (多普勒观测值 - 几何距离变化率)
	double w[128] = { 0 };

	// 有效卫星数量（满足观测值+星历有效）
	int valid_satnum = 0;

	// 3. 遍历所有卫星，构建观测方程
	for (int i = 0; i < Epoch->SatNum; i++)
	{
		// 过滤无效卫星：观测值无效 / 卫星位置无效 → 跳过
		if (!Epoch->SatObs[i].Valid || !Epoch->SatPVT[i].Valid)
			continue;

		// 4. 计算接收机到卫星的几何距离 ρ
		// 站星坐标差 = 接收机位置 - 卫星位置
		double xsr = Result->Position[0] - Epoch->SatPVT[i].SatPos[0];
		double ysr = Result->Position[1] - Epoch->SatPVT[i].SatPos[1];
		double zsr = Result->Position[2] - Epoch->SatPVT[i].SatPos[2];

		// 几何距离
		double rou = sqrt(xsr * xsr + ysr * ysr + zsr * zsr);
		if (rou < 1e-3)  // 距离异常，直接跳过
			continue;

		// 5. 计算几何距离变化率 ρ_dot (站星相对速度)
		double roudot = -(xsr * Epoch->SatPVT[i].SatVel[0]
			+ ysr * Epoch->SatPVT[i].SatVel[1]
			+ zsr * Epoch->SatPVT[i].SatVel[2]) / rou;

		// 6. 构建测速观测矩阵 B   行数=有效卫星数，列数=4（dx/dy/dz + 接收机钟漂）
		// B = [ -dx/ρ  -dy/ρ  -dz/ρ  1 ]
		B[valid_satnum][0] = xsr / rou;
		B[valid_satnum][1] = ysr / rou;
		B[valid_satnum][2] = zsr / rou;
		B[valid_satnum][3] = 1;  // 接收机钟漂系数

		// 7. 构建多普勒残差 w
		// 残差 = 观测多普勒 - (几何距离变化率 - 光速×卫星钟漂)
		w[valid_satnum] = Epoch->SatObs[i].D[0] - (roudot - C_Light * Epoch->SatPVT[i].SatClkSft);

		// 有效卫星数 +1
		valid_satnum++;
	}

	// 8. 卫星数不足4颗 → 无法解算4个参数（dx,dy,dz,钟漂）
	if (valid_satnum < 4)
		return;

	// 9. 构建最小二乘法方程
	// 法方程系数矩阵 N = B^T * B
	double N[4][4] = { 0 };
	// 法方程右端项 W = B^T * w
	double W[4] = { 0 };

	// 计算 N 和 W
	for (int i = 0; i < valid_satnum; i++) {
		for (int r = 0; r < 4; r++) {
			for (int c = 0; c < 4; c++) {
				N[r][c] += B[i][r] * B[i][c];
			}
			W[r] += B[i][r] * w[i];
		}
	}

	// 10. 法方程矩阵求逆 N^(-1)
	double N_inv[4][4];
	if (!matrix_invert4x4(N, N_inv))
		return; // 求逆失败 → 退出

	// 11. 最小二乘解：
	// X = [dVx, dVy, dVz, dtr] (接收机速度 + 钟漂)
	double X[4] = { 0 };
	for (int r = 0; r < 4; r++) {
		for (int c = 0; c < 4; c++) {
			X[r] += N_inv[r][c] * W[c];//X = N^(-1) * W
		}
	}

	// 12. 保存解算的接收机速度
	Result->Velocity[0] = X[0];
	Result->Velocity[1] = X[1];
	Result->Velocity[2] = X[2];

	// 13. 验后单位权中误差（测速精度评估）
	double vTv = 0;
	for (int i = 0; i < valid_satnum; i++) {
		double res = B[i][0] * X[0] + B[i][1] * X[1]
			+ B[i][2] * X[2] + B[i][3] * X[3]
			- w[i];
		vTv += res * res;
	}

	// 单位权中误差 = sqrt(残差平方和 / 自由度)
	Result->SigmaVel = sqrt(vTv / (valid_satnum - 4));
}