# GNSS

卫星导航程序设计



\# GNSS 卫星导航程序设计

本项目实现GNSS单点定位(SPP)与载波相位平滑伪距定位(SPV)，基于C++开发，可处理NovAtel OEM系列接收机日志，完成GPS/BDS单系统与组合系统定位解算。



\## ✨ 功能说明

1\. 解析 NovAtel OEM 二进制日志，提取伪距、载波相位、星历等观测数据

2\. 卫星位置、地球坐标系、站心坐标系之间的坐标转换

3\. 卫星钟差、电离层、对流层等定位误差模型修正

4\. \*\*SPP 伪距单点定位\*\*：支持仅GPS、GPS+BDS组合定位

5\. \*\*SPV 载波相位平滑伪距定位\*\*：提升定位平滑度与精度

6\. 输出定位结果文本，保存到Result文件夹，便于绘图与精度评估



\## 📁 项目结构

SPP\_SPV/

├── main.cpp                 # 程序入口

├── 3CoordSystemTrans.cpp    # 坐标转换模块

├── 36OutPutResult.cpp        # 结果输出模块

├── GnssConsts.h              # 常数定义

├── GnssStructs.h             # 数据结构体

├── GnssFuncDeclare.h         # 函数声明

├── NovatelOEM\*.log           # NovAtel 接收机原始观测数据

├── Result/

│   ├── Result\_OnlyGPS.txt    # 仅 GPS 定位结果

│   └── Result\_GPSBDS.txt     # GPS+BDS 组合定位结果

├── SPP\_SPV.sln               # VS 工程解决方案

└── SPP\_SPV.vcxproj           # VS 项目文件





\## 🛠 编译环境

\- 开发语言：C++

\- 开发工具：Visual Studio 2019 / 2022

\- 依赖：标准C++库，无第三方GNSS库



\## 🚀 使用方法

1\. 使用Visual Studio打开 `SPP\_SPV.sln`

2\. 编译生成可执行文件

3\. 程序读取根目录下 NovatelOEM\*.log 观测文件

4\. 运行后定位结果自动写入`Result/`目录下txt文件



\## 📌 说明

> 本项目为卫星导航原理课程设计，用于学习GNSS单点定位、相位平滑伪距算法，仅用于教学学习。





