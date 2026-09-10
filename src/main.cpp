/*-----------------------------------------------------------------------------
    main.cpp  —— GPS/BDS 双频 RTK 相对定位主程序（对应报告 3.2.5 main.cpp）

    用法：GNSS_RTK [config.ini]   默认读取 config/config.ini

    流程：
      1. 读取配置文件
      2. 根据 PosMode 打开基准站 / 流动站的二进制文件或 TCP 连接
      3. 打开输出文件并写表头
      4. 循环：时间同步 -> 基准站 SPP -> 流动站 SPP -> 单差 -> 周跳探测
               -> 参考星选取 -> 浮点解(最小二乘 / 卡尔曼) -> Lambda 固定
               -> 输出屏幕 / 结果文件 / NMEA -> 清空历元数据
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include "OutputNMEA.h"
#include <cstdio>
#include <cmath>

SolveMode g_solveMode;   // 原 SPP 工程的全局变量
int mode;

// 流动站参考坐标：优先配置文件，其次 BESTPOS；都没有时返回全 0
static void get_reference_xyz(ConfigInfo& cfg, BestPos& rovBestPos, double refXYZ[3])
{
    if (fabs(cfg.RovX) > 1.0 || fabs(cfg.RovY) > 1.0 || fabs(cfg.RovZ) > 1.0)
    {
        refXYZ[0] = cfg.RovX;
        refXYZ[1] = cfg.RovY;
        refXYZ[2] = cfg.RovZ;
        return;
    }
    if (rovBestPos.Time.Week > 0 && fabs(rovBestPos.Pos[1]) > 1e-6)
    {
        // DecodePos 存储顺序为 Pos[0]=经度 Pos[1]=纬度 Pos[2]=椭球高，与 BLH_2_XYZ 的输入一致
        BLH_2_XYZ(rovBestPos.Pos, refXYZ, R_WGS84, F_WGS84);
        return;
    }
    refXYZ[0] = refXYZ[1] = refXYZ[2] = 0.0;
}

// 一个历元结束：清空单差 / 双差信息，滚动卡尔曼滤波结构体
static void end_epoch(RTKData& rtk, ConfigInfo& cfg)
{
    rtk.SdObs.clear_except_comObs();
    rtk.DDObs.clear();
    if (cfg.CalcMode == 1)
    {
        rtk.PreKalman = rtk.NowKalman;
        rtk.NowKalman.clear();
        rtk.NowKalman.IsInit = false;
    }
}

int main(int argc, char* argv[])
{
    // ===================== 1. 读取配置 =====================
    string cfgFile = (argc > 1) ? argv[1] : "config/config.ini";
    ConfigInfo cfg;
    if (!load_config(cfgFile, cfg))
    {
        printf("Load config failed: %s\n", cfgFile.c_str());
        return -1;
    }
    printf("Config loaded: %s\n", cfgFile.c_str());
    printf("PosMode = %d (%s), CalcMode = %d (%s)\n",
        cfg.PosMode, cfg.PosMode == 0 ? "file" : "network",
        cfg.CalcMode, cfg.CalcMode == 0 ? "least squares" : "kalman filter");

    g_solveMode = GPS_BDS;
    mode = 1;

    // ===================== 2. 打开数据源 =====================
    ifstream basFile, rovFile;
    SOCKET basSock = INVALID_SOCKET, rovSock = INVALID_SOCKET;

    if (cfg.PosMode == 0)
    {
        basFile.open(cfg.BasBinFile.c_str(), ios::binary);
        if (!basFile.is_open())
        {
            printf("Cannot open base file: %s\n", cfg.BasBinFile.c_str());
            return -1;
        }
        rovFile.open(cfg.RovBinFile.c_str(), ios::binary);
        if (!rovFile.is_open())
        {
            printf("Cannot open rover file: %s\n", cfg.RovBinFile.c_str());
            return -1;
        }
        printf("Base file : %s\nRover file: %s\n", cfg.BasBinFile.c_str(), cfg.RovBinFile.c_str());
    }
    else
    {
        if (!OpenSocket(basSock, cfg.BasIP.c_str(), (unsigned short)cfg.BasPort))
        {
            printf("Cannot connect base station %s:%d\n", cfg.BasIP.c_str(), cfg.BasPort);
            WSACleanup();
            return -1;
        }
        if (!OpenSocket(rovSock, cfg.RovIP.c_str(), (unsigned short)cfg.RovPort))
        {
            printf("Cannot connect rover station %s:%d\n", cfg.RovIP.c_str(), cfg.RovPort);
            CloseSocket(basSock);
            WSACleanup();
            return -1;
        }
        printf("Base : %s:%d\nRover: %s:%d\n", cfg.BasIP.c_str(), cfg.BasPort, cfg.RovIP.c_str(), cfg.RovPort);
    }

    // ===================== 3. 输出文件 =====================
    FILE* fpOut = fopen(cfg.OutputFile.c_str(), "w");
    if (fpOut == NULL)
    {
        printf("Cannot open output file: %s\n", cfg.OutputFile.c_str());
        return -1;
    }
    FILE* fpNMEA = fopen(cfg.NMEAOutputFile.c_str(), "w");
    if (fpNMEA == NULL)
    {
        printf("Cannot open NMEA file: %s\n", cfg.NMEAOutputFile.c_str());
        fclose(fpOut);
        return -1;
    }
    write_header_RTK(fpOut);

    // ===================== 4. 解算所需结构体 =====================
    RTKData rtk;
    PosRes basPosRes, rovPosRes;
    vector<unsigned char> basBuff(MAXRAWLEN, 0), rovBuff(MAXRAWLEN, 0);
    size_t basBytes = 0, rovBytes = 0;

    bool basKnown = (fabs(cfg.BasX) > 1.0 || fabs(cfg.BasY) > 1.0 || fabs(cfg.BasZ) > 1.0);
    double basKnownPos[3] = { cfg.BasX, cfg.BasY, cfg.BasZ };

    int epochCount = 0, floatCount = 0, fixedCount = 0;

    // ===================== 5. 逐历元处理 =====================
    while (true)
    {
        // 5.1 时间同步
        int ret;
        if (cfg.PosMode == 0)
            ret = get_synch_obs_file(basFile, basBuff.data(), basBytes, rovFile, rovBuff.data(), rovBytes, rtk, RTK_SYNC_DT);
        else
            ret = get_synch_obs_net(basSock, basBuff.data(), basBytes, rovSock, rovBuff.data(), rovBytes, rtk, RTK_SYNC_DT);
        if (!ret)
        {
            printf("No more synchronized data.\n");
            break;
        }
        epochCount++;

        // 5.2 基准站 SPP（同时得到卫星位置、高度角）
        DetectOutlier(&rtk.BasEpkData);
        bool basOK = SPP(&rtk.BasEpkData, rtk.GPSEphemList.data(), rtk.BDSEphemList.data(), &basPosRes);
        if (basKnown)
        {
            // 基准站坐标已知：位置用已知值，SPP 失败时仍用已知坐标计算卫星位置
            if (!basOK)
                ComputeSatPVTAtSignalTrans(&rtk.BasEpkData, rtk.GPSEphemList.data(), rtk.BDSEphemList.data(), basKnownPos);
            for (int k = 0; k < 3; k++) basPosRes.Position[k] = basKnownPos[k];
        }
        else if (!basOK)
        {
            printf("%4d %10.3f  Base SPP failed.\n", rtk.BasEpkData.Time.Week, rtk.BasEpkData.Time.SecOfWeek);
            end_epoch(rtk, cfg);
            continue;
        }

        // 5.3 流动站 SPP + SPV
        DetectOutlier(&rtk.RovEpkData);
        bool rovOK = SPP(&rtk.RovEpkData, rtk.GPSEphemList.data(), rtk.BDSEphemList.data(), &rovPosRes);
        if (!rovOK)
        {
            rovPosRes.Type = None;
            printf("%4d %10.3f  Rover SPP failed.\n", rtk.RovEpkData.Time.Week, rtk.RovEpkData.Time.SecOfWeek);
            end_epoch(rtk, cfg);
            continue;
        }
        SPV(&rtk.RovEpkData, &rovPosRes);
        rovPosRes.Type = SPP_SPV;

        // 5.4 单差 -> 周跳探测 -> 参考星（双差信息）
        form_sd_obs(rtk.BasEpkData, rtk.RovEpkData, rtk.SdObs, cfg);
        detect_cycle_slip(rtk.SdObs);
        int nDD = DetRefSat(rtk.BasEpkData, rtk.RovEpkData, rtk.SdObs, rtk.DDObs);

        // 5.5 浮点解
        int floatOK = 0;
        if (nDD >= RTK_MIN_DD_SAT)
        {
            if (cfg.CalcMode == 0)
                floatOK = LSS_RTK_float(rtk, basPosRes, rovPosRes);
            else
                floatOK = float_kalman(rtk.BasEpkData, rtk.RovEpkData, basPosRes, rovPosRes,
                    rtk.SdObs, rtk.DDObs, rtk.PreKalman, rtk.NowKalman);
        }
        else
        {
            printf("%4d %10.3f  Not enough double-difference satellites (%d).\n",
                rtk.RovEpkData.Time.Week, rtk.RovEpkData.Time.SecOfWeek, nDD);
        }

        // 5.6 Lambda 固定 + Ratio 检验
        if (floatOK)
        {
            floatCount++;
            if (RTK_fixed(rtk.DDObs, rovPosRes, basPosRes, cfg)) fixedCount++;
        }
        else
        {
            printf("%4d %10.3f  RTK solution failed, SPP result only.\n",
                rtk.RovEpkData.Time.Week, rtk.RovEpkData.Time.SecOfWeek);
        }

        // 5.7 输出
        double refXYZ[3];
        get_reference_xyz(cfg, rtk.RovBestPos, refXYZ);
        write_2_screen_RTK(rovPosRes, basPosRes, rtk.DDObs, rtk.SdObs, refXYZ);
        write_2_file_RTK(fpOut, rovPosRes, basPosRes, rtk.DDObs, rtk.SdObs, refXYZ);
        print_NMEA(rovPosRes, fpNMEA);

        // 5.8 清空当前历元
        end_epoch(rtk, cfg);
    }

    // ===================== 6. 收尾 =====================
    printf("Epochs: %d, RTK float: %d, RTK fixed: %d\n", epochCount, floatCount, fixedCount);

    fclose(fpOut);
    fclose(fpNMEA);
    if (cfg.PosMode == 0)
    {
        basFile.close();
        rovFile.close();
    }
    else
    {
        CloseSocket(basSock);
        CloseSocket(rovSock);
        WSACleanup();
    }
    return 0;
}
