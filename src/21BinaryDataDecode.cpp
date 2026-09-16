#include "GnssConsts.h"
#include "GnssStructs.h"
#include "GnssFuncDeclare.h"

/*----------------------------------------------
字节提取函数
------------------------------------------------*/
double R8(unsigned char* p) // 提取 8 字节 double
{
    double r;
	memcpy(&r, p, 8);// 从字节数组 p 中复制 8 字节到 r      memcpy(目标地址，源地址，拷贝字节数量)  &取地址
    return r;
}
float R4(unsigned char* p)  // 提取 4 字节 float
{
    float r;
    memcpy(&r, p, 4);
    return r;
}
int I4(unsigned char* p)    // 提取 4 字节有符号整数
{
    int r;
    memcpy(&r, p, 4);
    return r;
}
unsigned int UI4(unsigned char* p)  // 提取 4 字节无符号整数
{
    unsigned int r;
    memcpy(&r, p, 4);
    return r;
}
short I2(unsigned char* p)  // 提取 2 字节有符号整数
{
    short r;  //短整型
    memcpy(&r, p, 2);
    return r;
}
unsigned short UI2(unsigned char* p)    // 提取 2 字节无符号整数
{
    unsigned short r;
    memcpy(&r, p, 2);
    return r;
}


// CRC32校验
#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */

/* crc-32 parity ---------------------------------------------------------------
* compute crc-32 parity for novatel raw
* args   : unsigned char *buff I data
*          int    len    I      data length (bytes)
* return : crc-32 parity
* notes  : see NovAtel OEMV firmware manual 1.7 32-bit CRC
*-----------------------------------------------------------------------------*/
unsigned int crc32(const unsigned char* buff, int len)
{
    int i, j;
    unsigned int crc = 0;

    for (i = 0; i < len; i++) {
        crc ^= buff[i];
        for (j = 0; j < 8; j++) {
            if (crc & 1) crc = (crc >> 1) ^ POLYCRC32;
            else crc >>= 1;
        }
    }
    return crc;
}

/*-----------------------------------------------
    观测值数据解码
    时间 → 观测值数量 → 观测值 → 存数据 → 卫星数
------------------------------------------------*/
void DecodeRange(unsigned char* data, EPOCHOBS* obs)
{
    GNSSSys sys;                             // 卫星系统类型     0 GPS，4 北斗
    int i, n = -1, j, k = -1;                 // 循环变量卫星索引  
    //i =第几条观测    同一颗卫星多条频点观测，记录在同一个 n 的位置     k记录最后一颗有效卫星所在数组下标
    // k 初值必须为 -1：本历元一颗卫星都没解出来时，SatNum 才会是 0 而不是 1
    int Prn, Freq;                            // 卫星PRN号，频率索引    系统+PRN号唯一确定一颗卫星，Freq = 0/1 代表第一/第二频点
    double wl;                                // 载波波长   
    unsigned int ChanStatus;                  // 通道状态字 
    int PhaseLockFlag, CodeLockedFlag;        // 相位/码锁定标志  1有效，0无效清零
    int ParityFlag, SatSystem, SigType;       // 奇偶校验、卫星系统、信号类型
    // p：数据指针，跳过报文头部，指向RANGE消息载荷起始位置
    unsigned char* p = data + NOVATEL_HEADER_LEN; // 指针指向消息数据体

    // 1解析观测历元时间（GPS周+周内秒）GPS 周、周内毫秒在报文头部固定偏移,不是在载荷里面
    obs->Time.Week = UI2(data + 14);
    obs->Time.SecOfWeek = UI4(data + 16) * 1E-3;  // ms -> s

    // 2获取当前历元观测值总数，初始化观测值存储区
    int ObsNum = UI4(p + OFFSET_OBS_NUM);

    // 清空观测值之前，先把上一历元的锁定时间按 PRN 存下来，用于本历元判断失锁（LLI bit0）
    double lockPrev[2][MAXBDSNUM][2];
    for (int s = 0; s < 2; s++)
        for (int q = 0; q < MAXBDSNUM; q++)
            for (int f = 0; f < 2; f++)
                lockPrev[s][q][f] = obs->LockPrev[s][q][f];

    memset(obs->SatObs, 0, MAXCHANNUM * sizeof(SATOBSDATA));//memset 将卫星观测数组全部清零，清除上一历元残留脏数据。

    // 本历元没再出现的卫星，锁定时间置回 -1，下次出现时按“重新捕获”处理
    for (int s = 0; s < 2; s++)
        for (int q = 0; q < MAXBDSNUM; q++)
            obs->LockPrev[s][q][0] = obs->LockPrev[s][q][1] = -1.0;

    // 3循环解析每一个观测值
    //for(【初始化语句】;【循环条件】;【每次循环末尾执行】)
    for (i = 0, p += 4; i < ObsNum; i++, p += RANGE_OBS_LEN)//p += 4 先偏移4字节后再每解析完一条，指针 p 向后跳跃一条记录长度RANGE_OBS_LEN。
    {
        // 解码得到跟踪状态标记，从中取出Phase lock flag / Code locked / flag / Parity known flag / Satellite system / signal type
        ChanStatus = UI4(p + OFFSET_CHAN_STATUS);
        //(变量 >> 移位位数) & 掩码      >> 右移：把目标 bit 移动到最低位    &掩码：只保留有效 bit，屏蔽其余位，提取标志位数值
        ParityFlag = (ChanStatus >> BIT_PARITY) & MASK_PARITY;
        PhaseLockFlag = (ChanStatus >> BIT_PHASE_LOCK) & MASK_PHASE_LOCK;
        CodeLockedFlag = (ChanStatus >> BIT_CODE_LOCK) & MASK_CODE_LOCK;
        SatSystem = (ChanStatus >> BIT_SAT_SYS) & MASK_SAT_SYS;
        SigType = (ChanStatus >> BIT_SIG_TYPE) & MASK_SIG_TYPE;

        // 判断卫星系统与信号类型， 如果GPS卫星的信号类型不是L1 C / A或者L2P（Y），BDS卫星不是B1I或B3I
        // continue至①，并记录信号频率类型，第一频率s = 0，第二频率s = 1；
        if (SatSystem == 0)//GPS
        {
            sys = GPS;
            if (SigType == 0) { Freq = 0; wl = WL1_GPS; }      // L1 C/A     wl载入对应载波波长，用于载波相位（周）→距离转换。
            else if (SigType == 9) { Freq = 1; wl = WL2_GPS; } // L2P(Y)
            else continue; // 其余信号直接丢弃，不解析
        }
        else if (SatSystem == 4)//BDS
        {
            sys = BDS;
            if (SigType == 0 || SigType == 4) { Freq = 0; wl = WL1_BDS; } // B1I
            else if (SigType == 2 || SigType == 6) { Freq = 1; wl = WL3_BDS; } // B3I
            else continue;
        }
        else continue; // 其他系统（GLONASS/Galileo）直接跳过

        // 获取卫星PRN号，查找对应存储位置（已存在/空槽位）
        Prn = UI2(p + OFFSET_PRIMARY_PRN);
        n = -1;                 // 每条观测都要重新定位槽位，找不到就丢弃，不能沿用上一条的 n
        for (j = 0; j < MAXCHANNUM; j++)
        {
            //一条观测对应一个频点。同一颗卫星多条频点观测，要存到同一个卫星结构体内部不同 Freq 下标，不能分开存多条卫星。
            // full overwrite 
            if (obs->SatObs[j].System == sys && obs->SatObs[j].Prn == Prn)
            {
                n = j; break;
            }
            //empty record
            if (obs->SatObs[j].Prn == 0)
            {
                k = n = j; break;//k记录当前历元有效卫星数量，n记录当前观测值存储位置
            }
        }
        if (n < 0) continue;    // 通道数已满，丢弃这条观测

        // 填充观测值数据（伪距相位多普勒载噪比锁定时间校验标志）
        obs->SatObs[n].Prn = Prn;
        obs->SatObs[n].System = sys;
        // 码环锁定才读取伪距，否则赋值0
        obs->SatObs[n].P[Freq] = CodeLockedFlag == 1 ? R8(p + OFFSET_PRIMARY_P) : 0.0;
        // 载波相位单位是【周】，转为距离：距离 = -波长 × 相位周数；相位未锁定置0
        obs->SatObs[n].L[Freq] = -wl * (PhaseLockFlag == 1 ? R8(p + OFFSET_PRIMARY_L) : 0.0);
        // 多普勒：多普勒距离变化率 = -波长 × 多普勒(Hz)
        obs->SatObs[n].D[Freq] = -wl * R4(p + OFFSET_PRIMARY_D);
        obs->SatObs[n].cn0[Freq] = R4(p + OFFSET_PRIMARY_CN0);
        double lockTime = R4(p + OFFSET_PRIMARY_LOCK);
        obs->SatObs[n].LockTime[Freq] = lockTime;
        obs->SatObs[n].half[Freq] = ParityFlag;

        // RINEX 失锁标记 LLI：bit0 = 与上一历元之间失过锁，bit1 = 存在半周模糊
        // NovAtel 没有直接给 LLI，按锁定时间和奇偶校验标志推：
        //   相位环没锁 或 锁定时间比上一历元变小（计数被清零重新累积）=> 失锁
        //   Parity known flag = 0 => 相位可能差半周
        unsigned char lli = 0;
        int sysIdx = (sys == GPS) ? 0 : 1;
        double prev = (Prn >= 1 && Prn <= MAXBDSNUM) ? lockPrev[sysIdx][Prn - 1][Freq] : -1.0;
        if (PhaseLockFlag != 1) lli |= 1;
        else if (prev >= 0.0 && lockTime < prev) lli |= 1;
        if (ParityFlag == 0) lli |= 2;
        obs->SatObs[n].LLI[Freq] = lli;

        if (Prn >= 1 && Prn <= MAXBDSNUM)
            obs->LockPrev[sysIdx][Prn - 1][Freq] = lockTime;
    }

    // 4统计当前历元有效卫星数量
    obs->SatNum = k + 1;
}

/*-----------------------------------------------
GPS星历解码
输入：整条二进制报文缓冲区 unsigned char* data
输出：解析后的 GPS 星历存入数组 GPSEPHREC geph[]
作用：把接收机下发的二进制星历 → 提取出 GPS 标准 16 参数星历，供后续卫星位置计算（SPP 定位核心前置步骤）。
------------------------------------------------*/
void DecodeGpsEphem(unsigned char* data, GPSEPHREC geph[])
{
    unsigned char* p = data + NOVATEL_HEADER_LEN;  // 指针跳过报文头部，指向星历消息载荷起始位置
    //PRN
    int prn = UI4(p);
    if (prn < 1 || prn > MAXGPSNUM) return;
    GPSEPHREC* eph = geph + (prn - 1);//使用指针eph，不用反复写geph[prn-1]，
    //卫星号 + 系统
    eph->PRN = prn;
    eph->System = GPS;
    //星历参考历元    TOE：Time of Ephemeris，星历参考时间,GPS 时（周 + 周内秒）    TOC：Time of Clock，钟差参考时间
    eph->TOE.Week = UI4(p + 24);
    eph->TOE.SecOfWeek = R8(p + 32);
    eph->TOC.Week = eph->TOE.Week;
    eph->TOC.SecOfWeek = R8(p + 164);
    //星历标识     IODE：星历数据龄期（轨道参数版本号       IODC：钟差数据龄期（钟参数版本号）
    eph->IODE = UI4(p + 16);    // IODE1
    eph->IODC = UI4(p + 160);   // GPSEPHEM 报文里 IODC 在偏移 160，偏移 20 是 IODE2
    //卫星健康状态
    eph->SVHealth = UI4(p + 12);   //非 0 代表卫星故障，定位时剔除该卫星
    //轨道根数改正项
    eph->SqrtA = sqrt(R8(p + 40));    // √A
    eph->M0 = R8(p + 56);          // 平近点角
    eph->e = R8(p + 64);          // 偏心率
    eph->omega = R8(p + 72);          // 近地点角距

    eph->Cuc = R8(p + 80);      // cuc
    eph->Cus = R8(p + 88);      // cus
    eph->Crc = R8(p + 96);      // crc
    eph->Crs = R8(p + 104);     // crs
    eph->Cic = R8(p + 112);     // cic
    eph->Cis = R8(p + 120);     // cis

    eph->DeltaN = R8(p + 48);   // Δn
    eph->OMEGA = R8(p + 144);   // Ω0 升交点经度
    eph->OMEGADot = R8(p + 152);
    eph->i0 = R8(p + 128);      // i0
    eph->iDot = R8(p + 136);    // IDOT

    //钟差
    eph->ClkBias = R8(p + 180);     // 钟偏
    eph->ClkDrift = R8(p + 188);        // 钟速
    eph->ClkDriftRate = R8(p + 196);   // 钟漂
    //群延时
    // GPSEPHEM 报文只播发一个 TGD（对应 L1），BDS 才有 TGD1/TGD2 两个。
    // 本工程 GPS 用双频无电离层组合定位，TGD 已被组合消去，故两个字段取同值不影响解算。
    eph->TGD1 = R8(p + 172);
    eph->TGD2 = R8(p + 172);

}

/*-----------------------------------------------
BDS星历解码（同GPS
------------------------------------------------*/
void DecodeBdsEphem(unsigned char* data, GPSEPHREC beph[])
{
    unsigned char* p = data + NOVATEL_HEADER_LEN;// skip header

    int prn = UI4(p);
    if (prn < 1 || prn > MAXBDSNUM) return;

    GPSEPHREC* eph = beph + (prn - 1);

    eph->PRN = prn;
    eph->System = BDS;

    eph->TOE.Week = UI4(p + 4);
    eph->TOE.SecOfWeek = UI4(p + 72);
    eph->TOC.Week = eph->TOE.Week;
    eph->TOC.SecOfWeek = UI4(p + 40);

    eph->SVHealth = UI4(p + 16);

    eph->TGD1 = R8(p + 20);
    eph->TGD2 = R8(p + 28);

    eph->ClkBias = R8(p + 44);
    eph->ClkDrift = R8(p + 52);
    eph->ClkDriftRate = R8(p + 60);

    eph->IODC = UI4(p + 36);
    eph->IODE = UI4(p + 68);
    // 轨道根数改正项
    eph->SqrtA = R8(p + 76);
    eph->e = R8(p + 84);
    eph->omega = R8(p + 92);
    eph->DeltaN = R8(p + 100);
    eph->M0 = R8(p + 108);

    eph->OMEGA = R8(p + 116);
    eph->OMEGADot = R8(p + 124);
    eph->i0 = R8(p + 132);
    eph->iDot = R8(p + 140);

    eph->Cuc = R8(p + 148);//周期扰动
    eph->Cus = R8(p + 156);
    eph->Crc = R8(p + 164);
    eph->Crs = R8(p + 172);
    eph->Cic = R8(p + 180);
    eph->Cis = R8(p + 188);
}

/*-----------------------------------------------
PSRPOS解码
------------------------------------------------*/
void DecodePos(unsigned char* data, POSRES* pos)
{
    unsigned char* p = data + NOVATEL_HEADER_LEN;// 跳过28字节消息头

    //解状态 & 位置类型
    pos->sol_status = UI4(p);       // Solution status解算状态（有效解、残差超限、缺少观测等）
    pos->pos_type = UI4(p + 4);     // Position type定位类型，NovAtel 定义：单点定位、伪距差分、RTK 固定、RTK 浮动等。

    //位置
    pos->Pos[1] = R8(p + 8);        // Latitude
    pos->Pos[0] = R8(p + 16);       // Longitude
    pos->Pos[2] = R8(p + 24) + R4(p + 32);       // 高程（m）大地高（海拔高）= 椭球高 − N

    // 大地水准面高度差 N
    pos->undulation = R4(p + 32);

    // 4标准差 
    pos->lat_sigma = R4(p + 40);// 纬度标准差（m）
    pos->lon_sigma = R4(p + 44);// 经度标准差（m）
    pos->hgt_sigma = R4(p + 48);/// 高度标准差（m）

    // 组合成三维位置标准差（平方和开根号）
    pos->SigmaPos = sqrt(
        pos->lat_sigma * pos->lat_sigma +
        pos->lon_sigma * pos->lon_sigma +
        pos->hgt_sigma * pos->hgt_sigma
    );// 位置标准差（m）

    // 卫星数量 
    pos->SatNum_tracked = *(p + 64); // #SVs 接收机当前跟踪卫星总数
    pos->SatNum_used = *(p + 65);    // #solnSVs 真正参与本次定位解算的卫星数

    // 未提供字段 报文不提供的字段，手动初始化
    pos->SigmaVel = 0.0;
    pos->PDOP = 0.0;

    // 时间 
    pos->Time.Week = UI2(data + 14);
    pos->Time.SecOfWeek = UI4(data + 16) * 1E-3;//SecOfWeek 原始单位毫秒，×1e-3 转为秒
}

/*-----------------------------------------------
主解码函数
------------------------------------------------*/

int DecodeNovOem7Dat(
    unsigned char buff[],  // 接收环形缓冲区
    int& len,              // buff里面当前有效字节长度（引用，函数内会修改）
    EPOCHOBS* obs,         // 观测结果输出结构体
    GPSEPHREC geph[],      // GPS星历数组  为什么是 GPSEPHREC geph[] 而不是 GPSEPHREC* geph？
    GPSEPHREC beph[],      // BDS星历数组  为什么是 GPSEPHREC beph[] 而不是 GPSEPHREC* beph？
    POSRES* pos,           // POSRES定位结果结构体
    int mode)              // 工作模式
{
    int i = 0;
    int GotRange = 0;   // 本次调用是否解出了新的 RANGE 历元
    while (i < len)
    {
        // 1、查找同步字符 AA 44 12
        if (i + 3 > len) break;
        if (!(buff[i] == NOVATEL_SYNC1 && buff[i + 1] == NOVATEL_SYNC2 && buff[i + 2] == NOVATEL_SYNC3))
        {
            i++;
            continue;
        }//最后更新得到 i = 当前报文起始地址

        // 2、检查是否有足够的字节读取完整的消息头
        if (i + NOVATEL_HEADER_LEN > len) break;

        int MsgID = UI2(buff + i + 4);//buff+i = 当前报文起始地址
        int MsgLen = UI2(buff + i + 8);


        // 3、检查整条消息是否完整
        //NovAtel 一条完整消息总结构： 【报文头NOVATEL_HEADER_LEN】 + 【载荷MsgLen】 + 【4字节CRC校验码】     
        if (i + NOVATEL_HEADER_LEN + MsgLen + NOVATEL_CRC_LEN > len) break;

        // 4、CRC检验  ？怎么计算的  crc32函数实时计算
        // 计算当前报文的 CRC32 校验码，并与报文尾部的 CRC32 校验码进行比较，如果不一致，跳过当前报文继续查找下一个报文。
        if (crc32(buff + i, NOVATEL_HEADER_LEN + MsgLen) != UI4(buff + i + NOVATEL_HEADER_LEN + MsgLen))
        {
            i += 3; continue;
        }

        // 5. 根据消息ID分发解码
        int Status = 0;
        switch (MsgID)
        {
            //调用解码函数Decode，把解析后的数据存入对应结构体obs
        case MSGID_RANGE:      DecodeRange(buff + i, obs);      Status = 1; GotRange = 1; break;//status = 1 表示已经成功解码了一条观测数据
        case MSGID_GPS_EPHEM:   DecodeGpsEphem(buff + i, geph);  break;//星历数据
        case MSGID_BDS_EPHEM:   DecodeBdsEphem(buff + i, beph);  break;
        case MSGID_PSRPOS:      DecodePos(buff + i, pos);        break;//定位数据
        case MSGID_BESTPOS:     DecodePos(buff + i, pos);        break;//BESTPOS 与 PSRPOS 消息体结构相同
        default: break;
        }

        // 移动指针，丢弃已处理数据
        i += NOVATEL_HEADER_LEN + MsgLen + NOVATEL_CRC_LEN;
        if (Status == 1 && mode == 1) break;
    }

    // 剩余数据前移
    memmove(buff, buff + i, len - i);//把buff+i 开始剩余数据前移到缓冲区buff开头， len - i 更新缓冲区有效长度
    len -= i;//len = len - i
    return GotRange;
}

/*-----------------------------------------------
    NetData Receive
------------------------------------------------*/

#pragma comment(lib,"WS2_32.lib")
#pragma warning(disable:4996)

/*-----------------------------------------------
    打开串口通信
------------------------------------------------*/
bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
    WSADATA wsaData;
    SOCKADDR_IN addrSrv;
    //初始化 Winsock 库
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    {
        printf("WSAStartup failed.\n");
        return false;
    }
    // 创建套接字
    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET)
    {
        printf("socket() failed.\n");
        return false;
    }

    memset(&addrSrv, 0, sizeof(addrSrv));
    addrSrv.sin_addr.s_addr = inet_addr(IP);
    addrSrv.sin_family = AF_INET;
    addrSrv.sin_port = htons(Port);

    if (connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR)) == SOCKET_ERROR)
    {
        printf("connect() to %s:%d failed.\n", IP, Port);
        closesocket(sock);
        sock = INVALID_SOCKET;
        return false;
    }

    return true;
}

/*-----------------------------------------------
    关闭串口通信
------------------------------------------------*/
void CloseSocket(SOCKET& sock)
{
    closesocket(sock);
    sock = INVALID_SOCKET;
}

/*-----------------------------------------------
    采集实时数据为二进制格式
------------------------------------------------*/
bool SaveSocketStreamToFile_1min(const char* ip, unsigned short port, const char* binFilePath, double hour)
{
    SOCKET NetGps;
    OpenSocket(NetGps, ip, port);//链接服务器

    FILE* fout = fopen(binFilePath, "wb");

    unsigned char buffer[MAXRAWLEN];
    int lenR = 0;

    DWORD startTick = GetTickCount64();// 开始时间

    while (true)
    {
        // 判断是否到达设定的接收时间
        DWORD nowTick = GetTickCount64();
        if (nowTick - startTick >= (DWORD)(hour * 3600.0 * 1000.0))
        {
            break;
        }


        lenR = recv(NetGps, (char*)buffer, MAXRAWLEN, 0);
        if (lenR > 0)
        {
            fwrite(buffer, 1, lenR, fout);
            fflush(fout);
        }
    }
    //关文件，断网，套接字
    fclose(fout);
    CloseSocket(NetGps);
    return true;
}
