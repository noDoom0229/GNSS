#pragma once
#ifndef GNSS_CONST_H
#define GNSS_CONST_H
#include <cmath>

// ============ 1.卫星系统枚举  =============
enum GNSSSys { UnknownSys = 255, GPS = 0, GLONASS = 1, GALILEO = 2, QZSS = 3, BDS = 4 };

// ============ 2.信号频段枚举 ============
enum SignalType { L1_CA = 0, L2_PY = 9, B1I_D1 = 0, B3I_D1 = 2, B1I_D2 = 4, B3I_D2 = 6, UnknownSig = 255 };

// ============ 3.解算模式枚举 ============
enum SolveMode { ONLY_GPS = 1, ONLY_BDS = 2, GPS_BDS = 3 };

// 全局变量声明（防止多重定义）
extern SolveMode g_solveMode;
extern int mode;

// 数学常量
#define PAI 3.14159265358979323846264338338279  // 圆周率
#define Rad (PAI / 180.0)                       // 角度转弧度系数
#define Deg (180.0 / PAI)                       // 弧度转角度系数

// WGS84椭球参数（GPS）
#define R_WGS84     6378137.0                    // 长半轴(m)
#define F_WGS84     (1.0 / 298.257223563)        // 扁率
#define W_WGS84     7.2921151467e-5              // 自转角速度(rad/s)
#define GM_WGS84    3986005e+8                   // 地心引力常数(m³/s²)

// CGCS2000椭球参数（BDS）
#define R_CGCS2000  6378137.0                    // 长半轴(m)
#define F_CGCS2000  (1.0 / 298.257222101)        // 扁率
#define W_CGCS2000  7.292115e-5                  // 自转角速度(rad/s)
#define GM_CGCS2000 3986004.418e+8               // 地心引力常数(m³/s²)

// 物理常数
#define C_Light 299792458.0                      // 真空中光速(m/s)

// 工作模式
#define FILEMODE 1                                // 1=文件解析，0=实时接收

// NovAtel协议帧头
#define NOVATEL_SYNC1           0xAA             // 同步字1
#define NOVATEL_SYNC2           0x44             // 同步字2
#define NOVATEL_SYNC3           0x12             // 同步字3
#define NOVATEL_HEADER_LEN      28               // 消息头长度(Byte)
#define NOVATEL_CRC_LEN         4                // CRC校验长度(Byte)

// 时间系统
#define GPST_BDT  14                             // GPS与北斗时差(s)

// 通道与卫星参数
#define MAXCHANNUM 36                            // 最大跟踪通道数
#define MAXGPSNUM  32                            // GPS最大PRN号
#define MAXBDSNUM 63                             // BDS最大PRN号
#define MAXRAWLEN 40960                          // 数据缓冲区长度(Byte)

// CRC校验
#define POLYCRC32   0xEDB88320u                  // CRC32多项式

// RANGE观测值消息
#define MSGID_RANGE             43               // RANGE消息ID
#define RANGE_OBS_LEN           44               // 单观测值长度(Byte)
#define OFFSET_OBS_NUM          0                // 观测值数量偏移
#define OFFSET_PRIMARY_PRN      0                // PRN编号偏移
#define OFFSET_PRIMARY_P        4                // 伪距偏移
#define OFFSET_PRIMARY_L        16               // 载波相位偏移
#define OFFSET_PRIMARY_D        28               // 多普勒偏移
#define OFFSET_PRIMARY_CN0      32               // 载噪比偏移
#define OFFSET_PRIMARY_LOCK     36               // 相位锁定时间偏移
#define OFFSET_CHAN_STATUS      40               // 通道状态字偏移

// 星历与定位消息ID
#define MSGID_GPS_EPHEM         7                // GPS星历ID
#define MSGID_BDS_EPHEM         1696             // BDS星历ID
#define MSGID_PSRPOS            42               // 单点定位ID
#define MSGID_BESTPOS           47               // 最优融合定位ID

// 通道状态位
#define BIT_PARITY              11               // 奇偶校验位
#define BIT_PHASE_LOCK          10               // 相位锁定位
#define BIT_CODE_LOCK           12               // 码锁定位
#define BIT_SAT_SYS             16               // 卫星系统位
#define BIT_SIG_TYPE            21               // 信号类型位

// 状态掩码
#define MASK_PARITY             0x01             // 奇偶校验掩码
#define MASK_PHASE_LOCK         0x01             // 相位锁定掩码
#define MASK_CODE_LOCK          0x01             // 码锁定掩码
#define MASK_SAT_SYS            0x07             // 卫星系统掩码
#define MASK_SIG_TYPE           0x1F             // 信号类型掩码

// 载波频率(Hz)
#define  FG1_GPS  1575.42E6                      // GPS L1频率
#define  FG2_GPS  1227.60E6                      // GPS L2频率
#define  FG1_BDS  1561.098E6                     // BDS B1I频率
#define  FG3_BDS  1268.520E6                     // BDS B3I频率

// 载波波长(m)
#define  WL1_GPS  (C_Light/FG1_GPS)              // GPS L1波长
#define  WL2_GPS  (C_Light/FG2_GPS)              // GPS L2波长
#define  WL1_BDS  (C_Light/FG1_BDS)              // BDS B1I波长
#define  WL3_BDS  (C_Light/FG3_BDS)              // BDS B3I波长

//  周跳/粗差探测阈值 
#define GFThres         0.5                   // GF 组合阈值
#define MWThres         1.0                   // MW 组合阈值

//对流层 Hopfield 模型
#define  Htrop        15000.0                // 对流层有效高度 (m)
#define  SeaLevel     0.0                    // 海平面/基准海拔
#define  TempT0       293.15                 // 标准地面温度 (20℃，常用标准)
#define  AtmosPre     1013.0                 // 标准大气压 (hPa，工程常用值)
#define  RelHum       0.75                   // 相对湿度 (0.75 为气象观测/导航解算默认值)

// 周跳检测阈值
#define  GFThres    0.05                          // GF检测阈值
#define  MWThres    3.0                           // MW检测阈值


#endif