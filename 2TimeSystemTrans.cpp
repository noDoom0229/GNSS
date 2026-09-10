#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"
#include <math.h>
/*-----------------------------------------------
	时间转换算法 -> 卫星导航算法与程序设计-4 时间与坐标算法
	通用时 < - > 儒略日 < - > GPS时间
                  ||
                年积日
 -----------------------------------------------*/

 /*-----------------------------------------------
     通用时到简化儒略日的转换算法
 ------------------------------------------------*/
void CommonTime_2_MjdTime(const COMMONTIME* ct, MJDTIME* mjd)
{
    // 处理月份小于3的情况,小于3则月+12，年-1
    double JD = (int)(365.25 * ((ct->Month < 2) ? (ct->Year - 1) : ct->Year))
        + (int)(30.6001 * (((ct->Month < 2) ? (ct->Month + 12) : ct->Month) + 1))
        + ct->Day
        + ((ct->Hour * 3600 + ct->Minute * 60 + ct->Second) / 86400.0)
        + 1720981.5;

    double MJD = JD - 2400000.5;

    mjd->Days = (int)(MJD);//地址指针指向结构体里的Days成员，存储整数天数
    mjd->FracDay = ((ct->Hour * 3600.0 + ct->Minute * 60.0 + ct->Second) / 86400.0); // 直接用通用时的 时分秒计算小数天
}

/*-----------------------------------------------
    通用时转化为GPS时间
------------------------------------------------*/
void CommonTime_2_GpsTime(const COMMONTIME* ct, GPSTIME* gps)
{
    MJDTIME mjd1;
    CommonTime_2_MjdTime(ct, &mjd1); //通用时 → 简化儒略日  &mjd1：取mjd1的地址传入函数
    MjdTime_2_GpsTime(&mjd1, gps); //简化儒略日 → GPS时间
}

/*-----------------------------------------------
    简化儒略日时间到通用时的转换算法
------------------------------------------------*/
void MjdTime_2_CommonTime(const MJDTIME* mjd, COMMONTIME* ct)
{
    // 简化儒略日转儒略日
    double JD = mjd->Days + mjd->FracDay + 2400000.5;

    int a = int(JD + 0.5);                // 对JD四舍五入，取整天数
    int b = a + 1537;                    // 历法修正常数
    int c = int((b - 122.1) / 365.25);   // 年份计数
    int d = int(365.25 * c);             // 年积日
    int e = int((b - d) / 30.6001);       // 月份计数

    // 计算年月日
    ct->Day = b - d - int(30.6001 * e) + (JD + 0.5 - int(JD + 0.5));// 小数部分用于恢复日期内的时分秒
    ct->Month = e - 1 - 12 * int(e / 14);
    ct->Year = c - 4715 - int((7 + ct->Month) / 10);

    // 计算时分秒
    double seconds_of_day = mjd->FracDay * 86400.0;
    ct->Hour = int(seconds_of_day / 3600);
    ct->Minute = int((seconds_of_day - ct->Hour * 3600) / 60);//分钟 =（总秒数 - 小时 ×3600）÷ 60
    ct->Second = seconds_of_day - ct->Hour * 3600 - ct->Minute * 60;//秒 = 总秒数 - 小时 ×3600 - 分钟 ×60
}

/*-----------------------------------------------
    简化儒略日时间到GPS时间的转换算法
------------------------------------------------*/
void MjdTime_2_GpsTime(const MJDTIME* mjd, GPSTIME* gps)
{
    gps->Week = (int)((mjd->Days - 44244.0) / 7.0 + mjd->FracDay / 7.0);
    gps->SecOfWeek = (mjd->Days + mjd->FracDay - 44244.0 - gps->Week * 7.0) * 86400.0;
}

/*-----------------------------------------------
    GPS时间转通用时间
------------------------------------------------*/
void GpsTime_2_CommonTime(const GPSTIME* gps, COMMONTIME* ct)
{
    MJDTIME mjd1;
    GpsTime_2_MjdTime(gps, &mjd1);
    MjdTime_2_CommonTime(&mjd1, ct);
}

/*-----------------------------------------------
    GPS时间到简化儒略日时间的转换算法
------------------------------------------------*/
void GpsTime_2_MjdTime(const GPSTIME* gps, MJDTIME* mjd)
{
    // MJD = GPS起始MJD + 周数×7天 + 周内秒换算的天数
    double MJD = 44244.0 + gps->Week * 7 + gps->SecOfWeek / 86400.0;
    mjd->Days = (int)MJD;
    //fmod(SecOfWeek, 86400) / 86400.0
    mjd->FracDay = (gps->SecOfWeek - (int)(gps->SecOfWeek / 86400.0) * 86400.0) / 86400.0;
}
