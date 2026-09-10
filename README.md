# GNSS 卫星导航程序设计

本项目包含两部分：

1. **C++ SPP/SPV**：基于 NovAtel OEM 日志的单点定位与载波相位平滑伪距定位
2. **Python GNSS 数据质量分析 (`gnss_qc`)**：按 BD 420022-2019《北斗/GNSS 测量型接收机观测数据质量评估方法》对 RINEX 3.x（`.26o`）观测文件做质量评估

## 功能说明

### SPP / SPV（C++）

1. 解析 NovAtel OEM 二进制日志，提取伪距、载波相位、星历等观测数据
2. 坐标转换与误差模型修正
3. SPP 伪距单点定位（GPS / GPS+BDS）
4. SPV 载波相位平滑伪距定位
5. 结果输出到 `Result/`

### 观测数据质量分析（Python，`gnss_qc`）

依据 BD 420022-2019 实现：

| 模块 | 内容 |
|------|------|
| RINEX 读取 | RINEX 3.x 多系统观测文件（`.26o` / `.rnx`） |
| 完整率 | 双频有效历元 / 理论历元 |
| 周跳探测 | MW 递推检验 + GF 补充探测，统计 o/slps |
| 多路径 | 双频 MP1/MP2 组合 + 滑动窗口去模糊度 |
| 电离层残差 | 双频相位电离层延迟及变化率 IOD |
| 噪声 | 伪距/相位三次差噪声（式 21/23） |
| 载噪比 | 各频点 CNR 均值统计 |
| 输出 | 文本报告、CSV、PNG 图 |

## 项目结构

```
├── main.cpp / *.cpp / *.h     # SPP/SPV C++ 源码
├── NovatelOEM*.log            # NovAtel 原始日志
├── Result/                    # SPP/SPV 定位结果
├── gnss_qc/                   # RINEX 质量分析工具包
│   ├── rinex_reader.py
│   ├── quality.py
│   ├── report.py
│   ├── generate_sample.py
│   └── __main__.py
├── data/                      # 示例 / 用户 RINEX（.26o）
├── Result_QC/                 # 质量分析结果
└── requirements-qc.txt
```

## SPP/SPV 编译与运行

- 语言：C++；工具：Visual Studio 2019/2022
- 打开 `SPP_SPV.sln` 编译，读取 `NovatelOEM*.log`，结果写入 `Result/`

## 质量分析使用方法

```bash
pip install -r requirements-qc.txt

# 分析自己的 RINEX 观测文件（如 yygc1.26o）
python -m gnss_qc /path/to/yygc1.26o --out Result_QC

# 未提供文件时自动生成与用户格式一致的示例 .26o 并分析
python -m gnss_qc --out Result_QC

# 可选：只分析部分系统 / 限制历元
python -m gnss_qc data/yygc1.26o --systems G,C --max-epochs 3600
```

输出：

- `Result_QC/qc_report.txt` — 完整率、周跳比、MP、IOD、噪声、CNR
- `Result_QC/qc_summary.csv` — 逐卫星指标
- `Result_QC/qc_overview.png` 等 — 可视化图

将实测 `.26o` 放到 `data/` 后直接指定路径即可。

## 说明

本项目为卫星导航原理课程相关实现，用于学习 GNSS 定位与观测数据质量评估，仅供教学学习。
