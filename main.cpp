#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <windows.h>

using namespace std;

SolveMode g_solveMode;
int mode;


int main()
{
    // ===================== 模式选择 =====================
    int mode = 1;
    while (true)
    {
        printf("请输入模式选择（1为事后定位 2为实时定位）：");
        int input_ret = scanf_s("%d", &mode);

        if (input_ret == 1 && (mode == 1 || mode == 2))
        {
            break;
        }
        else
        {
            printf("无效的模式选择，请输入 1 或 2 ！\n");
            while (getchar() != '\n');        //清除缓冲区，防止死循环
        }
    }

    if (mode == 1)
    {
        printf("你选择了事后定位模式。\n");
        int solMode = 0;
        while (mode == 1)
        {
            printf("\n请选择解算系统：\n");
            printf("1 - 仅 GPS 解算\n");
            printf("2 - 仅 BDS 解算\n");
            printf("3 - GPS + BDS 联合解算\n");
            printf("请输入序号：");
            int ret = scanf_s("%d", &solMode);
            if (ret == 1 && solMode >= 1 && solMode <= 3)
            {
                g_solveMode = (SolveMode)solMode;
                break;
            }
            else
            {
                printf("输入无效，请重新输入！\n");
                while (getchar() != '\n');
            }
        }
    }
    else
        printf("你选择了实时定位模式。\n");


    POSRES recv_pos;        // 存储NovAtel接收机原生输出的定位结果
    EPOCHOBS current_obs;   
    RAWDAT nav_eph;         // 星历
	PPRESULT solve_result;  //最小二乘迭代解算结果
    int error_code = 0;     

    // ===================== 模式 1：事后定位 =====================
    if (mode == 1)
    {
        FILE* obs_file = nullptr;
        error_code = fopen_s(&obs_file, "NovatelOEM20211114-01.log", "rb");

        if (error_code != 0 || obs_file == nullptr)
        {
            printf("打开观测文件失败！\n");
            return -1;
        }

        unsigned char data_buf[MAXRAWLEN] = { 0 };
        int data_len = 0;
        int read_len = 0;

        while (1)
        {
            read_len = fread(data_buf + data_len, sizeof(unsigned char), MAXRAWLEN - data_len, obs_file);
            if (read_len <= 0) break;
            data_len += read_len;

            if (DecodeNovOem7Dat(data_buf, data_len, &current_obs, nav_eph.GpsEph, nav_eph.BdsEph, &recv_pos, mode) != 0) continue;
            DetectOutlier(&current_obs);

            bool spp_ok = SPP(&current_obs, &nav_eph, &solve_result);
            if (spp_ok)
            {
                SPV(&current_obs, &solve_result);

                const char* outFile;
                switch (g_solveMode)
                {
                case ONLY_GPS:
                    outFile = "Result/Result_OnlyGPS.txt";
                    break;
                case ONLY_BDS:
                    outFile = "Result/Result_OnlyBDS.txt";
                    break;
                case GPS_BDS:
                    outFile = "Result/Result_GPSBDS.txt";
                    break;
                default:
                    outFile = "Result/AfterResult.txt";
                    break;
                }
                OutputCombinedResult(&current_obs, &recv_pos, &solve_result, outFile);
            }
        }
        fclose(obs_file);
        printf("验证完成，结果已按模式输出到不同文件！\n");
    }

    // ===================== 模式 2：实时定位 =====================
    else if (mode == 2)
    {
        double recv_hour = 0.0;
        cout << "请输入接收数据的时间长度（小时）：" << endl;
        cin >> recv_hour;

        unsigned char net_buf[MAXRAWLEN] = { 0 };       //net_buf：网络总缓冲区，拼接分片 TCP 数据包；
        unsigned char temp_buf[MAXRAWLEN] = { 0 };      //temp_buf：单次 recv 临时接收缓存；
        SOCKET gps_socket = INVALID_SOCKET;             //gps_socket：TCP 套接字。

        bool socket_ok = OpenSocket(gps_socket, "8.148.22.229", 7003);
        if (!socket_ok)
        {
            printf("Socket 连接失败！\n");
            return -1;
        }

        int recv_len = 0;
        int cache_len = 0;
        DWORD time_begin = GetTickCount64();

        while (1)
        {
            DWORD time_now = GetTickCount64();
            ULONGLONG time_cost = time_now - time_begin;
            ULONGLONG time_limit = (ULONGLONG)(recv_hour * 60 * 60 * 1000);

            if (time_cost >= time_limit)
            {
                printf("接收已满 %.2f 小时，保存完毕。\n", recv_hour);
                break;
            }

            Sleep(980);//1Hz
            recv_len = recv(gps_socket, (char*)temp_buf, MAXRAWLEN, 0);

            if (recv_len > 0)
            {
                memcpy(net_buf + cache_len, temp_buf, recv_len);
                cache_len += recv_len;
                memset(temp_buf, 0, MAXRAWLEN);
                //二进制解码
                if (DecodeNovOem7Dat(net_buf, cache_len, &current_obs, nav_eph.GpsEph, nav_eph.BdsEph, &recv_pos, mode) != 0)
                {
                    continue;
                }
				// 观测值粗差检测
                DetectOutlier(&current_obs);
				// 单点定位迭代解算
                bool calc_ok = SPP(&current_obs, &nav_eph, &solve_result);
                if (calc_ok)
                {
					// 单点测速解算
                    SPV(&current_obs, &solve_result);
					// 实时输出结果
                    OutputResult_RealTime(&solve_result, &recv_pos, &current_obs);
                    OutputCombinedResult(&current_obs, &recv_pos, &solve_result, "Result/Net_SPP_SPV.txt");
                    OutputENUResult_RealTime(&current_obs, &recv_pos, &solve_result, "Result/Error_Calculation.txt");
                }
            }
        }

        closesocket(gps_socket);
    }

    printf("接受完成。\n");
    system("pause");
    return 0;
}