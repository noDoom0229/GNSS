#pragma once
/*-----------------------------------------------------------------------------
    OutputNMEA.h  —— 将定位结果以 NMEA 0183 $GPGGA 格式写入文件
                     （对应报告 3.1.1 OutputNMEA.h）

    $GPGGA,hhmmss.ss,ddmm.mmmmmm,N,dddmm.mmmmmm,E,fix,nsat,hdop,alt,M,geo,M,,*cs
      fix : 1 = 单点 SPP，4 = RTK 固定解，5 = RTK 浮点解
      alt : 高度 (m)；geo : 大地水准面差距 (m)
      最后两个逗号之间为差分数据龄期、差分站 ID，按报告留空
-----------------------------------------------------------------------------*/
#include "SPP.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// NMEA 校验和：'$' 与 '*' 之间所有字符按位异或
static int NMEACRC(char* buff, int len)
{
    int crc = 0;
    for (int i = 0; i < len; i++)
    {
        crc ^= (unsigned char)buff[i];
    }
    return crc;
}

static void print_NMEA(PosRes& rovPosRes, FILE* fpNMEA)
{
    if (fpNMEA == NULL) return;
    if (rovPosRes.Type == None) return;   // 没有定位结果，不输出

    // 1. GPS 时 -> UTC（减去跳秒）-> 年月日时分秒
    GPSTIME utc = rovPosRes.Time;
    utc.SecOfWeek -= GPS_UTC_LEAPSEC;
    if (utc.SecOfWeek < 0.0)
    {
        utc.SecOfWeek += 604800.0;
        utc.Week -= 1;
    }
    COMMONTIME ct;
    GpsTime_2_CommonTime(&utc, &ct);

    // 2. ECEF -> BLH（度）
    double blh[3];
    XYZ_2_BLH(rovPosRes.Position, blh, R_WGS84, F_WGS84);
    double lat = blh[0], lon = blh[1], hgt = blh[2];

    char ns = (lat >= 0.0) ? 'N' : 'S';
    char ew = (lon >= 0.0) ? 'E' : 'W';
    lat = fabs(lat);
    lon = fabs(lon);
    int latDeg = (int)lat;
    int lonDeg = (int)lon;
    double latMin = (lat - latDeg) * 60.0;
    double lonMin = (lon - lonDeg) * 60.0;

    // 3. 定位状态
    int fix = 1;
    if (rovPosRes.Type == RTKFixed) fix = 4;
    else if (rovPosRes.Type == RTKFloat) fix = 5;

    // 4. GPGGA 报文
    //    HDOP 字段填写原 SPP 计算的 PDOP；大地水准面差距无模型可用，填 0.0，
    //    此时 alt 直接为椭球高（alt + geo 仍等于椭球高，可被 RTKPLOT 正确还原）。
    char buff[256];
    snprintf(buff, sizeof(buff),
        "$GPGGA,%02d%02d%05.2f,%02d%09.6f,%c,%03d%09.6f,%c,%d,%02d,%.1f,%.3f,M,%.3f,M,,",
        ct.Hour, ct.Minute, ct.Second,
        latDeg, latMin, ns,
        lonDeg, lonMin, ew,
        fix, (int)rovPosRes.AllSatNum, rovPosRes.PDOP,
        hgt, 0.0);

    // 5. 校验和：跳过开头的 '$'
    int len = (int)strlen(buff);
    int cs = NMEACRC(buff + 1, len - 1);

    // 6. 写入文件
    fprintf(fpNMEA, "%s*%02X\r\n", buff, cs);
    fflush(fpNMEA);
}
