# GNSS_RTK —— GPS/BDS 双频 RTK 相对定位程序

基于 NovAtel OEM719 接收机二进制数据的 C++ 卫星导航定位程序。
在原有 SPP（单点定位/测速）工程基础上，按《卫星导航算法与程序设计 2 成果总结报告》
新增 RTK 相对定位模块：时间同步、站间单差、周跳探测、参考星选取、站星双差、
最小二乘 / 卡尔曼滤波浮点解、Lambda 模糊度固定、Ratio 检验、固定解、NMEA 与自定义结果输出。

## 1. 功能

| 类别 | 功能 |
|------|------|
| 数据 | 读取二进制文件或 TCP 实时流；解码广播星历（GPS/BDS）、RANGE 观测值、PSRPOS/BESTPOS |
| SPP | 卫星位置/速度/钟差/钟速、Hopfield 对流层、地球自转改正、MW/GF 粗差探测、最小二乘定位与测速 |
| RTK | 基准站/流动站时间同步、站间单差、单差 MW/GF 周跳探测、参考星选取（跳过 GEO）、站星双差 |
| 浮点解 | 最小二乘（相关双差权阵）或卡尔曼滤波（支持参考星变化、新星升起、卫星消失） |
| 固定解 | Lambda（LD 分解、降相关、整数搜索）、Ratio 检验（阈值 3）、固定解基线与协方差 |
| 精度评定 | 单位权中误差 σ0、残差 RMS、PDOP、E/N/U 方向中误差 |
| 输出 | 屏幕、自定义结果文件、NMEA 0183 `$GPGGA`（WGS84） |

系统与频率：GPS L1 + L2，BDS B1I + B3I。采样率最高 1 Hz。

## 2. 环境

- Windows 10/11，Visual Studio 2022（C++17），Winsock2（`ws2_32.lib`，CMake 已自动链接）
- CMake ≥ 3.15
- 源文件统一为 UTF-8（带 BOM），CMake 已为 MSVC 加上 `/utf-8`

## 3. 目录结构

```
GNSS_RTK/
├── CMakeLists.txt
├── README.md
├── include/
│   ├── GnssConsts.h        原有：常量、枚举
│   ├── GnssStructs.h       原有：数据结构（新增 resultType、PPRESULT::Type）
│   ├── GnssFuncDeclare.h   原有：SPP 函数声明（新增 vector 矩阵运算声明）
│   ├── SPP.h               新增：RTK 结构体、常量、函数声明（总头文件）
│   └── OutputNMEA.h        新增：NMEACRC()、print_NMEA()
├── src/
│   ├── 1MatrixAndVector.cpp    原有 + 任意维 vector 矩阵运算
│   ├── 2TimeSystemTrans.cpp    原有：时间转换
│   ├── 3CoordSystemTrans.cpp   原有：坐标转换
│   ├── 21BinaryDataDecode.cpp  原有：NovAtel 解码、Socket
│   ├── 32SatellitePVT.cpp      原有：卫星位置速度钟差
│   ├── 33ErrorCorrect.cpp      原有：对流层、粗差探测
│   ├── 34SPP35SPV.cpp          原有：SPP / SPV
│   ├── 36OutPutResult.cpp      原有：SPP 结果输出
│   ├── config.cpp              新增：load_config()
│   ├── RTK.cpp                 新增：时间同步、单差、周跳、参考星、LS 浮点解、固定解、精度评定
│   ├── Kalman.cpp              新增：卡尔曼滤波浮点解
│   ├── lambdaN.cpp             新增：Lambda 模糊度搜索
│   ├── writeToFile.cpp         新增：RTK 结果输出
│   └── main.cpp                RTK 主流程
├── config/config.ini
├── data/                       放置 base/rover 二进制文件
└── output/                     结果输出目录
```

## 4. 编译

### Visual Studio 2022

1. 打开 VS2022 → “打开本地文件夹”，选择工程根目录（含 `CMakeLists.txt`）。
2. 选择 x64-Debug 或 x64-Release 配置，生成 → 全部生成。
3. 调试工作目录已设为工程根目录，`config/`、`data/`、`output/` 相对路径可直接使用。

### 命令行（Windows）

```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
build\Release\GNSS_RTK.exe config\config.ini
```

## 5. 配置文件 `config/config.ini`

```ini
PosMode = 0            ; 0 = 二进制文件，1 = 网络实时流
CalcMode = 0           ; 0 = 最小二乘浮点解，1 = 卡尔曼滤波浮点解

BasBinFile = data/base.bin
RovBinFile = data/rover.bin

BasIP = 8.148.22.229
BasPort = 7002
RovIP = 8.148.22.229
RovPort = 4002

NMEAOutputFile = output/result.nmea
OutputFile = output/result.txt

BasX = 0   BasY = 0   BasZ = 0   ; 基准站已知坐标（全 0 时用基准站 SPP 结果）
RovX = 0   RovY = 0   RovZ = 0   ; 流动站参考坐标，只用于 dE/dN/dU（全 0 时用流动站 BESTPOS）

ElevThreshold = 10     ; 高度角阈值 (deg)
PseuThreshold = 10     ; 双频伪距差阈值 (m)
RatioThreshold = 3     ; Ratio 检验阈值
```

`#` 或 `;` 开头为注释，空行跳过。IP/端口不写在源码中，全部由配置文件给出。

其余报告未给数值的参数集中在 `include/SPP.h` 顶部（时间同步阈值 `RTK_SYNC_DT`、
最小二乘收敛条件、伪距/相位先验中误差 `RTK_SIGMA_CODE/PHASE`、卡尔曼初始/过程噪声等）。

## 6. 数据文件

- 文件模式：基准站与流动站各一个 NovAtel OEM7 二进制文件（如报告中的
  `zeroline_oem719-202202021500-base.bin` / `-rover.bin`），放在 `data/` 并在配置中指定。
- 文件需包含 RANGE(43)、GPS 星历(7)、BDS 星历(1696)、PSRPOS(42)/BESTPOS(47) 消息。

## 7. 文件流运行

```bat
GNSS_RTK.exe config\config.ini      # PosMode = 0
```

程序逐历元处理直到两个文件任一结束，结束时打印历元数、浮点解数、固定解数。

## 8. 网络实时运行

将 `PosMode = 1`，填写 `BasIP/BasPort/RovIP/RovPort` 后运行。程序通过 Winsock2 建立两条 TCP 连接，
接收实时数据流并逐历元解算；连接断开或接收失败时退出。

## 9. 输出格式

### 自定义结果文件（`OutputFile`），每历元一行

```
WEEK SOW X Y Z dX dY dZ Ratio SatNum Type dE dN dU  Sigma0 RMS PDOP mE mN mU
```

第 1~14 列严格按报告表 3；第 15~20 列是课件「精度评定与结果输出模块」要求的精度信息。

| 字段 | 含义 |
|------|------|
| WEEK / SOW | GPS 周 / 周内秒 |
| X Y Z | 流动站 WGS84 ECEF 坐标 (m) |
| dX dY dZ | 基线向量 流动站 − 基准站 (m) |
| Ratio | Lambda Ratio 值（未做固定时为 0） |
| SatNum | 单差卫星数 |
| Type | 1 = SPP，2 = RTK Float，3 = RTK Fixed |
| dE dN dU | 相对参考坐标的 ENU 误差 (m)，无参考坐标时为 `nan` |
| Sigma0 | 单位权中误差 √(VᵀPV/r)，理想情况接近 1；明显大于 1 说明先验中误差设小了或存在粗差 |
| RMS | 双差残差均方根 √(VᵀV/n) (m) |
| PDOP | 流动站 SPP 的 PDOP |
| mE mN mU | E/N/U 方向坐标中误差 (m)，由 Q_NEU = K Q_XYZ Kᵀ 得到 |

只有 RTK 解算成功的历元才有精度信息，纯 SPP 历元这 6 列输出 0。

### NMEA 文件（`NMEAOutputFile`）

```
$GPGGA,hhmmss.ss,ddmm.mmmmmm,N,dddmm.mmmmmm,E,fix,nsat,hdop,alt,M,geo,M,,*cs
```

fix：1 = SPP，4 = RTK Fixed，5 = RTK Float。UTC = GPS 时 − 18 s。
HDOP 字段填写 SPP 的 PDOP；大地水准面差距填 0.0（alt 为椭球高），可直接载入 RTKPLOT 查看。

## 10. RTK 算法流程

```
读取配置 → 打开基/流数据 → 输出表头
循环每个历元：
  时间同步 get_synch_obs_file / get_synch_obs_net（|t_rov − t_bas| < 0.001 s）
  基准站 SPP（基准站坐标已知时使用配置坐标）
  流动站 SPP + SPV
  form_sd_obs        站间单差 ΔP = P_R − P_B，ΔL = L_R − L_B（剔除无效/双频不全/伪距差大/低高度角，标记半周）
  detect_cycle_slip  单差 MW、GF 组合历元间差分探测周跳
  DetRefSat          每系统选高度角最大且无 GEO/半周/周跳/位置无效的卫星为参考星
  浮点解：CalcMode 0 → LSS_RTK_float   X = (BᵀPB)⁻¹BᵀPW，Q = (BᵀPB)⁻¹
          CalcMode 1 → float_kalman    预测 Φ、Q → 更新 K、V、R（Joseph 形式协方差）
  RTK_fixed          lambda() 得两组整数解 → Ratio = s2/s1 ≥ 阈值 → 固定解基线
                     b_fix = b_float − Q_ba Q_aa⁻¹ (a_float − a_fix)
                     Q_fix = Q_bb − Q_ba Q_aa⁻¹ Q_ab
  calc_rtk_quality   在最终解上重建 B/P/V，算 σ0、RMS、mE/mN/mU
  输出屏幕 / 结果文件 / NMEA → 清空历元数据
```

观测值与参数排列（RTK.cpp 顶部注释）：模糊度 `[GPS f1 | GPS f2 | BDS f1 | BDS f2]`；
观测行 `GPS: P f1, P f2, L f1, L f2; BDS 同`。

**单位约定**：原工程解码时已把载波相位折算成米（`L = -λ·ADR`），所以 `SATOBSDATA::L`、
`SDSatObs::dL` 的单位都是米而不是周（原代码注释写的是“周”，实际是米）。
只有模糊度参数是周，观测方程里用 `λ·N` 换算。

**权阵**（报告 (2-24)，课件 II-3 P42~P44）：非差等方差 σ²、互不相关时，
站间单差 `cov(SD) = 2σ²I`，站星双差 `cov(DD) = 2σ²(I + 1·1ᵀ)`，因此

```
P = 1/(2σ²) · 1/(n+1) · [ n  −1 … ; −1  n … ]
```

伪距、相位分别取 `RTK_SIGMA_CODE`、`RTK_SIGMA_PHASE`（非差中误差，默认 0.3 m / 0.003 m），
GPS 与 BDS、伪距与相位、f1 与 f2 各自成块，整体块对角。

## 11. 常见错误

| 现象 | 原因 / 处理 |
|------|-------------|
| `Cannot open config file` | 工作目录不对，或路径写错；VS 调试工作目录应为工程根目录 |
| `Cannot open base file / rover file` | `data/` 中没有对应二进制文件 |
| `connect() to x.x.x.x:port failed` | 网络不通或服务未开放；确认配置中的 IP/端口 |
| 开头若干历元 `Rover SPP failed` | 粗差探测需要前一历元组合值，且星历尚未解到，属正常现象 |
| `Not enough double-difference satellites` | 高度角/伪距差阈值过严，或两站共视卫星少 |
| `RTK solution failed, SPP result only.` | 浮点解失败（法方程奇异、迭代不收敛），本历元只输出 SPP |
| Ratio 一直小于 3 | 观测质量差、电离层活跃；可检查 `RTK_SIGMA_CODE/PHASE`，或数据本身问题 |
| 一个历元都算不出来 | 先看时间同步：`RTK_SYNC_DT` 默认 0.001 s（课件要求），若两站采样时刻本身有偏差需放宽 |
| Sigma0 明显大于 1 | 先验中误差设小了或存在未探测到的粗差/周跳 |
| 卡尔曼频繁 `filter re-initialized` | 历元间隔大于 `KF_MAX_GAP`（默认 1.5 s），数据缺失导致 |
| MSVC 报 C4819 编码警告 | 确认已用 CMake 生成（含 `/utf-8`） |

## 12. 测试方法

1. 将基准站/流动站二进制文件放入 `data/`，配置 `PosMode = 0`。
2. 运行程序，观察屏幕每历元一行输出；`Type` 列应逐渐出现 2（浮点）和 3（固定）。
3. 检查 `output/result.txt` 与 `output/result.nmea`；NMEA 可载入 RTKPLOT 查看轨迹与 fix 状态。
4. 已知流动站精确坐标时填入 `RovX/Y/Z`，`dE/dN/dU` 列即为定位误差；固定解应在厘米级以内。
5. 将 `CalcMode` 改为 1 重复测试卡尔曼滤波；对比两者 Ratio 与固定率。
6. 网络测试：`PosMode = 1`，填写 IP/端口后运行。

**说明**：本仓库不含实测数据。代码已在 Linux(g++) 下完成编译（零警告）、单元测试
（任意维矩阵求逆、权阵数值与课件 P43 例子一致、P 与 R 互逆、Lambda 搜索、NMEA 校验和、配置读取）
以及仿真几何下的端到端验证：已知基线 + 已知整周模糊度，14 颗星（8 GPS + 6 BDS）单历元双频解算，
最小二乘与卡尔曼均正确固定 24 个双差模糊度，固定解基线误差 1~5 mm。
**尚未用真实 OEM719 数据验证**，真实数据下的精度与固定率需要在 Windows 下运行确认。
