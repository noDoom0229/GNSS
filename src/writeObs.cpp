/*-----------------------------------------------------------------------------
    writeObs.cpp —— 解码结果按自定义格式输出（讲义「数据解码」模块要求）

    基准站与流动站各自输出两个文件，用于和 RINEX 转换结果逐历元对照检查：
      观测值文件：历元头一行 + 每颗卫星一行，含 GPS L1/L2、BDS B1I/B3I 的
                  伪距、载波相位、多普勒、信噪比、LLI
      星历文件：  每颗卫星一条最新广播星历，字段顺序与 RINEX NAV 一致

    单位换算：程序内部相位与多普勒都乘过波长（米、米/秒），而 RINEX 里相位是
    周、多普勒是 Hz，所以输出时都除以波长还原，方便直接对比。
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cmath>

// 系统字母：与 RINEX 3 一致，G = GPS，C = BDS
static char sys_char(GNSSSys sys)
{
    if (sys == GPS) return 'G';
    if (sys == BDS) return 'C';
    return '?';
}

// ==========================================================
//   观测值文件
// ==========================================================
void write_obs_header(FILE* fp, const char* station)
{
    if (fp == NULL) return;
    fprintf(fp, "%% GNSS_RTK decoded observation file    station: %s\n", station);
    fprintf(fp, "%% 历元行: > YYYY MM DD HH MM SS.SSS  WEEK  SOW  SatNum\n");
    fprintf(fp, "%% 卫星行: SYS PRN  P1 L1 D1 S1 LLI1   P2 L2 D2 S2 LLI2\n");
    fprintf(fp, "%%   SYS  G = GPS (L1/L2)   C = BDS (B1I/B3I)\n");
    fprintf(fp, "%%   P 伪距 (m)   L 相位 (cycle)   D 多普勒 (Hz)   S 载噪比 (dB-Hz)\n");
    fprintf(fp, "%%   LLI  bit0 = 失锁重捕, bit1 = 半周模糊；观测值为 0 表示该频点未观测到\n");
    fprintf(fp, "%%   时间为 GPS 时；BDS 观测值的时间标签同样是 GPS 时\n");
    fflush(fp);
}

void write_obs_epoch(FILE* fp, EpochData& obs)
{
    if (fp == NULL || obs.SatNum <= 0) return;

    COMMONTIME ct;
    GpsTime_2_CommonTime(&obs.Time, &ct);

    // 只统计真正解出观测值的卫星
    int nSat = 0;
    for (int i = 0; i < obs.SatNum; i++)
        if (obs.SatObs[i].Prn != 0) nSat++;

    fprintf(fp, "> %04d %02d %02d %02d %02d %6.3f %5d %10.3f %3d\n",
        ct.Year, ct.Month, ct.Day, ct.Hour, ct.Minute, ct.Second,
        obs.Time.Week, obs.Time.SecOfWeek, nSat);

    for (int i = 0; i < obs.SatNum; i++)
    {
        SATOBSDATA& o = obs.SatObs[i];
        if (o.Prn == 0) continue;

        fprintf(fp, " %c%02d", sys_char(o.System), o.Prn);
        for (int f = 0; f < 2; f++)
        {
            double wl = get_wavelength(o.System, f);
            // 内部存的是米，还原成 RINEX 的周与 Hz
            double Lcyc = (wl > 0.0) ? o.L[f] / wl : 0.0;
            double Dhz = (wl > 0.0) ? o.D[f] / wl : 0.0;
            fprintf(fp, " %14.3f %15.3f %11.3f %5.1f %2d",
                o.P[f], Lcyc, Dhz, o.cn0[f], (int)o.LLI[f]);
        }
        fprintf(fp, "\n");
    }
    fflush(fp);
}

// ==========================================================
//   星历文件
// ==========================================================
// 按 RINEX NAV 的排列输出：第一行为 PRN + 参考时刻 + 三个钟参数，
// 其后每行 4 个轨道参数。本工程没有解码的字段（如 L2 码指示、拟合区间）留空为 0。
static void write_one_eph(FILE* fp, const GPSEPHREC& e)
{
    COMMONTIME ct;
    GPSTIME toc = e.TOC;
    // BDS 星历的 TOC/TOE/Week 存的都是北斗时。BDS 周数加 1356 后周内秒不变，
    // 按 GPS 时的算法转出来的年月日时分秒正好就是北斗时的钟面时，
    // 与 RINEX 中 BDS 导航电文用北斗时标注历元的做法一致。
    if (e.System == BDS) toc.Week = (unsigned short)(toc.Week + 1356);
    GpsTime_2_CommonTime(&toc, &ct);

    fprintf(fp, "%c%02d %04d %02d %02d %02d %02d %5.1f %20.12e %20.12e %20.12e\n",
        sys_char(e.System), e.PRN,
        ct.Year, ct.Month, ct.Day, ct.Hour, ct.Minute, ct.Second,
        e.ClkBias, e.ClkDrift, e.ClkDriftRate);

    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n", e.IODE, e.Crs, e.DeltaN, e.M0);
    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n", e.Cuc, e.e, e.Cus, e.SqrtA);
    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n",
        e.TOE.SecOfWeek, e.Cic, e.OMEGA, e.Cis);
    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n", e.i0, e.Crc, e.omega, e.OMEGADot);
    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n",
        e.iDot, 0.0, (double)e.TOE.Week, 0.0);
    fprintf(fp, "    %20.12e %20.12e %20.12e %20.12e\n",
        0.0, (double)e.SVHealth, e.TGD1, e.System == BDS ? e.TGD2 : e.IODC);
}

int write_nav_file(const string& path, GPSEPHREC* gpsEph, GPSEPHREC* bdsEph, const char* station)
{
    if (path.empty()) return 0;

    FILE* fp = fopen(path.c_str(), "w");
    if (fp == NULL)
    {
        printf("Cannot open nav output file: %s\n", path.c_str());
        return 0;
    }

    fprintf(fp, "%% GNSS_RTK decoded broadcast ephemeris   station: %s\n", station);
    fprintf(fp, "%% 每颗卫星一条最新星历，字段顺序与 RINEX NAV 一致：\n");
    fprintf(fp, "%%   SYS PRN  TOC(年 月 日 时 分 秒)  ClkBias ClkDrift ClkDriftRate\n");
    fprintf(fp, "%%   IODE  Crs     DeltaN  M0\n");
    fprintf(fp, "%%   Cuc   e       Cus     SqrtA\n");
    fprintf(fp, "%%   TOE   Cic     OMEGA   Cis\n");
    fprintf(fp, "%%   i0    Crc     omega   OMEGADot\n");
    fprintf(fp, "%%   IDOT  -       Week    -\n");
    fprintf(fp, "%%   -     SVHealth TGD1   IODC(GPS) / TGD2(BDS)\n");
    fprintf(fp, "%% 时间系统与 RINEX 一致：GPS 星历用 GPS 时，BDS 星历用北斗时\n");
    fprintf(fp, "%% 首行的年月日时分秒由 TOC 换算，BDS 即为北斗时钟面时\n");

    int n = 0;
    for (int i = 0; i < MAXGPSNUM; i++)
    {
        if (gpsEph[i].PRN == 0 || gpsEph[i].System != GPS) continue;
        write_one_eph(fp, gpsEph[i]);
        n++;
    }
    for (int i = 0; i < MAXBDSNUM; i++)
    {
        if (bdsEph[i].PRN == 0 || bdsEph[i].System != BDS) continue;
        write_one_eph(fp, bdsEph[i]);
        n++;
    }

    fclose(fp);
    printf("Ephemeris written: %s (%d satellites)\n", path.c_str(), n);
    return n;
}
