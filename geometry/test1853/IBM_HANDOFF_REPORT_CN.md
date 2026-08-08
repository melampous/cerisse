# NASA LaRC UPWT Test 1853 几何与 IBM 测试交接报告

版本：2026-08-03，公开资料可审计重建、quad 几何扩展与 Cerisse IBM/AWS 启动审计版 1.3  
适用对象：Run 165 中心 single-nozzle，Runs 247/262/263 peripheral tri-nozzle，以及 Runs 307/315 quad-nozzle  
CAD/STL 主坐标单位：mm；IBM 预处理副本单位：m

## 0. 交接结论

这套交付包含五个闭合 CAD/STL 实体：Run 165 single、由于 NASA 未公开 Nozzle 3/4 方位对应关系而各自保留两种编号版本的 tri 与 quad。quad 是中心 N1 加外围 N2/N3/N4 同时启用，不是四个外围喷管。五套都通过了水密、单实体、外法向、退化面、重复面、非流形边和 GTS 自相交检查。

必须同时保留以下两个结论：

1. 这些模型是依据 NASA 官方报告、as-built Table A-2 和工程图建立的“公开资料最高可审计外部 CFD 重建”，不是未公开 NASA 原始 STEP 的逐点复制。
2. 当前 STL 的拓扑和三角形形状质量很好，适合可视化、IBM 接入和 bring-up；但全局最大边长为 `2.863651 mm`，不满足仓库 fixed-STL 文档采用的推荐 production audit gate `max edge / finest dx <= 0.5`。这个数值不是求解器运行时硬限制，也不是充分精度证明；正式亚毫米 IBM 网格序列仍应先按目标最细 `dx` 重新三角化和验证。

关于“圆角是否为 NASA 原始数据”的直接回答：

- 肩部 `R0.100 in = 2.540 mm` 是 NASA DWG 1168296 的直接图纸尺寸。
- 喷管图上的 `R0.125/R0.188/R0.250 in` 是 NASA 名义零件图尺寸，但公开图无法唯一恢复各圆弧、切点与 as-built 修正后的完整轴向轮廓，因此没有把它们擅自加进当前主 CFD 实体。
- 鼻球 `R=1.000 in` 由公开测点和相切几何强烈支持，是名义解析重建，不是公开表中明确给出的 as-built 半径。
- 当前 CAD 没有另外“凭视觉加圆角”。外围喷口斜唇是本重建中喷管锥壁与解析前体的 Boolean 交线，不是人为圆滑；它不等同于未公开 NASA CAD 的已知逐点交线。

### 数据状态标记

本报告采用以下标记，后续 IBM 对话不得混用：

- **[NASA]**：NASA 报告、工程图或运行表直接给出。
- **[AS-BUILT]**：NASA Table A-2 的实测加工喷管参数。
- **[DERIVED]**：由所列公开值或显式名义重建输入解析推导；若依赖名义输入，会在表中注明。
- **[ASSUMPTION]**：为闭合 CFD 几何而明确采用的替代假设。
- **[STL]**：对本次输出文件直接测量或审计所得。
- **[UNKNOWN]**：公开资料不能唯一确定。

## 1. 一手来源与范围

唯一主来源是 NASA 官方 [NTRS 记录](https://ntrs.nasa.gov/citations/20140006403)及其 [NASA/TP-2014-218256 官方 PDF](https://ntrs.nasa.gov/api/citations/20140006403/downloads/20140006403.pdf)。没有用仓库旧 `trinozzle.stl`、后续论文示意图或二手 CFD 模型反推几何。

直接使用的工程图包括：

| 报告打印页 | 图号 | 内容 |
|---:|---|---|
| 92–94 | DWG 1168288 | aftbody cover、剖面和测压孔 |
| 97 | DWG 1168291 | sting |
| 98–99 | DWG 1168292/1168293 | front plenum、interface plate |
| 100–101 | DWG 1168294 | center plug/nozzle 与名义流道 |
| 102–104 | DWG 1168295 | peripheral plug/nozzle 与仪器化 Nozzle 2 |
| 105–106 | DWG 1168296 | forebody OML、坐标和肩圆角 |
| 108–109 | DWG 1168298 | manifold 最终加工与焊接装配 |

原页无损摘图和来源索引见 [`source_extract/SOURCE.md`](source_extract/SOURCE.md)。最终装配图 DWG 1168299 被报告引用，但没有在公开报告中复制，这是若干全局站位仍未知的根本原因。

## 2. Run、硬件构型与风洞条件

不存在“Run 165 tri-nozzle”。正确对应关系如下。

| Run | 实际硬件 | 模型滚转 `phi` | 攻角序列 `alpha`, deg | `M_inf` | `Re_inf`, 1/ft | `C_T` | `P_TJ`, psia | `T_TJ`, R | 总流量, lbm/s |
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
| 165 | 中心 Nozzle 1；外围 2/3/4 为 flush plug | 179.93° | B4-A1: −8, −4, 0, 4, 8, 12, 16, 20, 0 | 4.6 | 1.50e6 | 1.968 | 598.93 | 622.08 | 0.620 |
| 247 | 中心 plug；外围 Nozzles 2/3/4 | −0.03° | B4-A1: −8, −4, 0, 4, 8, 12, 16, 20, 0 | 4.6 | 1.50e6 | 1.919 | 201.19 | 623.56 | 0.599 |
| 262 | 与 247 相同 tri 硬件 | 0.01° | B4-A2: −8, −4, 0, 4, 8, 12, 16, 0 | 4.6 | 1.50e6 | 2.937 | 302.81 | 624.18 | 0.919 |
| 263 | 与 247 相同 tri 硬件 | −179.97° | B4-A2: −8, −4, 0, 4, 8, 12, 16, 0 | 4.6 | 1.50e6 | 2.935 | 303.07 | 630.93 | 0.913 |
| 307 | 中心 N1 + 外围 N2/N3/N4，全四管 | 0.05° | B4-A1: −8, −4, 0, 4, 8, 12, 16, 20, 0 | 4.6 | 1.50e6 | 1.923 | 151.93 | 618.38 | 0.602 |
| 315 | 与 307 相同 quad 硬件，低推力 | 0.05° | B4-A1: −8, −4, 0, 4, 8, 12, 16, 20, 0 | 4.6 | 1.50e6 | 0.912 | 76.39 | 614.81 | 0.285 |

**[NASA]** Runs 247/262/263 的实体硬件相同；Runs 307/315 的 quad 实体硬件相同。它们不是各自不同的 CAD。Run 号包含整段攻角扫描，不能把任一 run 当成单一 `alpha=0` 工况。主同推力构型对比应使用 Run 247 tri (`C_T=1.919`) 与 Run 307 quad (`C_T=1.923`)；Run 315 是低推力 quad。

本次从一手运行表新增锁定的总状态原值如下；单位与 NASA 表保持一致，尚未在本几何任务中转换为求解器静态边界：

| Run | `PTINF`, psf | `TTINF`, R | `PTJ`, psia | `TTJ`, R | `WJ`, lbm/s |
|---:|---:|---:|---:|---:|---:|
| 247 | 3658.01 | 610.37 | 201.19 | 623.56 | 0.599 |
| 307 | 3657.79 | 609.68 | 151.93 | 618.38 | 0.602 |
| 315 | 3658.20 | 609.35 | 76.39 | 614.81 | 0.285 |

**[NASA]** Run 165 的 Table B-3 原始运行平均值还包括 `PTINF=3657.58 psf`、`TTINF=609.60 R`、打印后的静态 `PINF=11.2 psf`、`TINF=116.5 R`、`RHOINF=5.58e-5 slug/ft^3` 和 `QINF=165.38 psf`。**[DERIVED]** 当前 AWS case 没有沿用仓库另一个 Mach-4.6 工况的旧常数，而是采用数值模型假定的 `gamma=1.4`、`MW=0.02896 kg/mol`，从 Run 165 总状态一致推导出 `p_inf=534.5806702 Pa`、`T_inf=64.72986748 K`、`rho_inf=0.02876556480 kg/m^3`、`a_inf=161.2999343 m/s`、`u_inf=741.9796979 m/s`；NASA 表中的静态量只保留到较少有效位，所以与这些推导值有正常的舍入差异。

Nozzles 3/4 没有独立质量流量测量。若只做零阶面积比例核算，三喷管面积分数为 N2/N3/N4 = `0.329616/0.339449/0.330935`；它对应 Run 262 约 `0.302917/0.311954/0.304129 lbm/s`、Run 263 约 `0.300940/0.309917/0.302143 lbm/s`。这些是审计参考值，不是独立实验测量，也不应直接变成强制逐管流量边界；更合理的候选是给三管一致的实验 reservoir/plenum 条件，让不同 `At` 自然产生流量，再核对总量与面积比例。

## 3. 主坐标系、滚转与平移

### 3.1 CAD body frame

所有主 CAD、BREP、STEP 和 `outputs/stl/` 文件均使用 mm。STL 格式本身不携带单位；“mm”来自本生成流程的坐标约定，不是文件头里的可读元数据：

- 原点是 70° 理论尖锥点 TSP，`x_TSP=0`。
- `+x` 沿模型轴指向下游；风洞来流沿下游方向，反推喷流排向 `−x`。
- `r=sqrt(y^2+z^2)`。
- `y=r sin(theta)`，`z=r cos(theta)`。
- `theta=0` 位于模型顶部 `+z`，`theta=90°` 为 `+y`。按 NASA 从上游向下游看的 front view，`+y` 显示在画面左侧。
- master CAD 的滚转为 `phi=0`。

不要分别用各 STL 的 `xmin` 归一化位置。名义物理鼻尖统一位于 `x=1.630115421 mm`；Run 165 的中心喷管唇使实体 `xmin=2.441268921 mm`。若把两个构型各自按 `xmin` 移到零点，会给 Run 165 引入约 `0.811153 mm` 的相对错位。

### 3.2 SI/风洞坐标变换

提供的预处理工具使用

```text
p_world[m] = t_TSP[m] + R_x(phi) * (1.0e-3 * p_master[mm])
```

其中 `R_x(phi)` 是绕 `+x` 的模型滚转。运行表矩阵见 [`outputs/reference/run_roll_transforms.json`](outputs/reference/run_roll_transforms.json)。该文件只有滚转，没有攻角、平移或完整风洞姿态。

推荐保持模型轴沿计算域 `x`，通过来流向量实现每个攻角；这样整个攻角序列可固定同一 STL。若选择旋转实体而不是来流，必须把攻角变换也显式加入并重新导出、审计，不能误以为现有工具已经包含攻角。

## 4. 主体全部尺寸与解析 OML

### 4.1 尺寸表

| 项目 | inch | mm / SI | 状态与解释 |
|---|---:|---:|---|
| 最大直径 | 5.000 | 127.000 mm | **[NASA]** |
| 最大半径 | 2.500 | 63.500 mm | **[NASA]** |
| sphere-cone 半锥角 | 70° | 70° | **[NASA]**，相对模型轴 |
| 名义鼻球半径 | 1.000 | 25.400 mm | **[ASSUMPTION]** 名义重建，非公开 as-built 标注 |
| 肩部圆角 | 0.100 | 2.540 mm | **[NASA]** DWG 1168296 直接尺寸 |
| aft cover 前缘 | `x=1.0066` | 25.56764 mm | **[NASA]**，不是 OML 切点 |
| aft cover 长度 | 9.550 | 242.570 mm | **[NASA]** |
| aft end TSP 站位 | `x=10.5566` | 268.13764 mm | **[DERIVED]** |
| 物理鼻尖至 aft end | 10.492422228 | 266.507524579 mm | **[DERIVED]** |
| aft cover 壁厚 | 0.075 | 1.905 mm | **[NASA]**；主 CFD 实体按实心输出 |
| aft cover 名义内半径 | 2.425 | 61.595 mm | **[NASA]** |
| forebody 后结构面 | `x=1.177` | 29.8958 mm | **[NASA]**，非 OML 切点 |
| 参考面积 | 0.13635 ft² | 0.012667329504 m² | **[NASA]** |
| 参考长度 | 5.000 | 127.000 mm | **[NASA]**，不是模型总长 |
| 表面光洁度 | 63 microinch | 1.6002 micrometre | **[NASA]** 图纸 callout |

解析关键站位：

| 站位 | `x`, mm | `r`, mm | 状态 |
|---|---:|---:|---|
| 名义鼻球球心 | 27.030115421 | 0 | **[DERIVED]**，依赖名义 R1.000 |
| 名义物理鼻尖 | 1.630115421 | 0 | **[DERIVED]**，依赖名义 R1.000 |
| 球—锥相切点 | 3.161922853 | 8.687311640 | **[DERIVED]**，依赖名义 R1.000 |
| 锥—肩圆角相切点 | 22.503817749 | 61.828731164 | **[DERIVED]** |
| 肩圆弧圆心 | 24.890637015 | 60.960000000 | **[DERIVED]** |
| 肩—圆柱相切点 | 24.890637015 | 63.500000000 | **[DERIVED]** |

图纸中的 `x=0.877 in` 是 B–B 剖切/结构站，不是解析锥—圆角切点；`diameter 4.851 in` 是前后体配合区尺寸，不是肩前 OML 外径；`x=1.0066 in` 是 aft-cover 前缘，不是圆角进入圆柱的位置。

### 4.2 OML 方程

用 inch、TSP 为原点，当前公开资料重建的外形为

```text
sphere:   r = sqrt(1 - (x - 1.064177772476)^2)
          0.064177772476 <= x <= 0.124485151690

cone:     r = x tan(70 deg)
          0.124485151690 <= x <= 0.885977076739

shoulder: r = 2.4 + sqrt(0.1^2 - (x - 0.979946339176)^2)
          0.885977076739 <= x <= 0.979946339176

cylinder: r = 2.5
          0.979946339176 <= x <= 10.5566
```

**[ASSUMPTION]** aft end 在 `x=10.5566 in` 用平面封闭，以满足三维 IBM 水密要求。这个平 aft cap 不是完整实验 base/sting。若积分整个 STL 的表面力，它会带入非实验基底贡献，必须排除、单独积分或在报告中解释。

![主体子午面与关键站位](outputs/figures/body_meridional_section.png)

## 5. 喷管 as-built 尺寸

以下是 NASA Table A-2 的实际加工值，不是统一的名义 `0.25/0.50 in`。

| Nozzle | `Dt`, in / mm | virtual `De`, in / mm | `At`, in² / mm² | virtual `Ae`, in² / mm² | `Ae/At` | 扩张半角 | ideal `Me` | `lambda` |
|---:|---|---|---|---|---:|---:|---:|---:|
| 1 | 0.2491 / 6.32714 | 0.5014 / 12.73556 | 0.048735 / 31.44187260 | 0.197451 / 127.38748716 | 4.051550 | 14.981° | 2.9571 | 0.98301 |
| 2 | 0.2498 / 6.34492 | 0.4965 / 12.61110 | 0.049009 / 31.61864644 | 0.193610 / 124.90942760 | 3.950514 | 14.977° | 2.9305 | 0.98302 |
| 3 | 0.2535 / 6.43890 | 0.4989 / 12.67206 | 0.050471 / 32.56187036 | 0.195487 / 126.12039292 | 3.873210 | 14.994° | 2.9097 | 0.98298 |
| 4 | 0.2503 / 6.35762 | 0.4996 / 12.68984 | 0.049205 / 31.74509780 | 0.196036 / 126.47458576 | 3.984035 | 14.985° | 2.9394 | 0.98300 |

三喷管合计 `At=0.148685 in²=95.9256146 mm²`，virtual `Ae=0.585133 in²=377.50440628 mm²`。

所有喷管轴线都平行于全局 `x` 轴，喷流沿 `−x`。外围喷管不是沿 70° 前体局部法向安装。

## 6. 喷管全局位置、编号歧义与真实斜唇

### 6.1 轴线和喉面

外围节圆直径为 `2.500 in=63.500 mm`，轴线半径为 `1.250 in=31.750 mm`，三轴方位为 `0/120/240°`。NASA 只明确 Nozzle 2 位于 `theta=0°`；Nozzles 3/4 与 `120/240°` 的对应关系未公开。

master body frame 中的关键位置为：

| Nozzle | `theta` | 轴心 `(y,z)`, mm，N3@120 版 | virtual-plane `x`, mm | throat `x`, mm | `Rt`, mm | 喉面固体到流体法向 |
|---:|---:|---|---:|---:|---:|---|
| 1 | center | (0, 0) | 2.441268806 | 14.415425035 | 3.163570 | `(-1,0,0)` |
| 2 | 0° | (0, 31.750000000) | 11.556054938 | 23.267709419 | 3.172460 | `(-1,0,0)` |
| 3 | 120° | (27.496306570, −15.875000000) | 11.556054938 | 23.192163829 | 3.219450 | `(-1,0,0)` |
| 4 | 240° | (−27.496306570, −15.875000000) | 11.556054938 | 23.384524196 | 3.178810 | `(-1,0,0)` |

`N3@240` 版本把 N3 设为 `(−27.496306570,−15.875)`，N4 设为 `(+27.496306570,−15.875)`；N2 不变。

位置状态必须这样理解：

- N1 virtual plane 是“as-built virtual-exit 圆与名义球面齐平”的 **[ASSUMPTION]**。
- N2–N4 virtual plane 通过喷管轴线与理想 70° 锥的交点，是显式 **[ASSUMPTION]**，不是 NASA 给出的全局 `x`。
- 上表 throat `x` 由该 virtual plane、as-built `Dt/De` 和壁角推导。因此 throat 直径是 **[AS-BUILT]**，全局轴向站位仍受插入深度假设控制。

### 6.2 外围物理开口

外围可见开口是轴向扩张锥壁与 70° 前体的空间交线，不是 `De/2` 的平面圆。当前 Boolean 重建得到：

| Nozzle | 物理 lip `x` 范围, mm | 相对各轴线局部 `rho` 范围, mm |
|---:|---:|---:|
| 2 | 9.013452 – 13.647450 | 5.746062 – 6.985745 |
| 3 | 9.000834 – 13.657338 | 5.773229 – 7.020413 |
| 4 | 8.997422 – 13.660404 | 5.781652 – 7.029786 |

这条 lip 一般不是圆，也不是严格椭圆。virtual-exit disk 只作为面积/Mach 参考面输出，绝不能与真实物理开口混用。

![喷管前视布置](outputs/figures/nozzle_layout_front_view.png)

![外围物理斜唇变化](outputs/figures/peripheral_physical_lip_variation.png)

## 7. 公开但没有进入主 STL 的硬件尺寸

这些尺寸保留在 [`geometry_parameters.yaml`](geometry_parameters.yaml) 供装配审计，不应误以为当前外部 CFD 实体已经包含它们。

### 7.1 名义喷管流道和 insert

| 项目 | inch | mm | 当前处理 |
|---|---:|---:|---|
| 名义上游孔径 | 0.500 | 12.700 | 保留，未建完整收敛段 |
| 名义收敛半角 | 20° | 20° | 保留 |
| 名义喉径 | 0.250 | 6.350 | 主模型改用各枚 as-built `Dt` |
| 名义扩张半角 | 15° | 15° | 主模型改用各枚 as-built 角度 |
| 名义 virtual-exit 直径 | 0.500 | 12.700 | 主模型改用各枚 as-built `De` |
| 图纸圆弧 | R0.125 / R0.188 / R0.250 | R3.175 / R4.7752 / R6.350 | 引线分别属于不同连接，未擅自拼成 as-built contour |

| insert 项目 | center | peripheral | 状态 |
|---|---|---|---|
| 前部包络直径 | 0.705 in / 17.907 mm | 0.705 in / 17.907 mm | **[NASA]** |
| nozzle 法兰直径 | 0.925 in / 23.495 mm | 0.945 in / 24.003 mm | **[NASA]** |
| plug 法兰直径 | 0.925 in / 23.495 mm | 0.941 in / 23.9014 mm | **[NASA]** |
| 后部 stem 直径 | 0.562 in / 14.2748 mm | 0.562 in / 14.2748 mm | **[NASA]** |
| 法兰前/后局部 `s` | −1.000 / −0.500 in | −1.000 / −0.500 in | `s=0` 取后端面，向前为负 |
| 恒径入口开始 | `s=−1.317 in` / −33.4518 mm | `s=−1.049 in` / −26.6446 mm | **[NASA]** 局部参考 |
| nozzle CAD 面至法兰参考 | 1.218 in / 30.9372 mm | 0.987 in / 25.0698 mm | 不能替代 virtual-plane 尺寸链 |
| plug CAD 面至法兰参考 | 1.186 in / 30.1244 mm | 0.990 in / 25.1460 mm | 同上 |

公开图将暴露前表面标为 “THIS SURFACE DEFINED BY CAD FILE”。完整收敛段、真实 throat blend、plenum 与制造修正后的轴向轮廓均未进入当前主实体。

### 7.2 Plenum、manifold 和 interface plate

| 项目 | inch / in² | SI/mm | 状态 |
|---|---:|---:|---|
| plenum station A 面积 | 0.78540 in² | 506.708664 mm² | **[NASA]** |
| plenum station B 面积 | 0.89185 in² | 575.385946 mm² | **[NASA]**；不保证截面为圆 |
| N1/N2 局部上游面积 | 0.19635 in² | 126.677166 mm² | **[NASA]** |
| 内部 P/T 测点上游距离 | 0.417 in | 10.5918 mm | **[NASA]**，仅 N1/N2 有独立测量 |
| interface plate 外径 | 4.750 in | 120.650 mm | **[NASA]** |
| interface plate 厚度 | 0.550 in | 13.970 mm | **[NASA]** |
| peripheral PCD | 2.500 in | 63.500 mm | **[NASA]** |
| front plenum 原始长度 | 3.850 in | 97.790 mm | **[NASA]** |
| manifold 最终加工长度 | 3.750 in | 95.250 mm | **[NASA]** |
| 后法兰直径 | 1.960 in | 49.784 mm | **[NASA]** |
| 分支参考角 | 17° | 17° | **[NASA]** |
| 工作/水压试验压力 | 2500 / 3750 psi | — | **[NASA]** |

内部弯曲支路和汇合曲面是 **[UNKNOWN]**，未进入外部 CFD solid。

### 7.3 Sting

| 项目 | inch | mm |
|---|---:|---:|
| 零件图参考总长 | 18.500 | 469.900 |
| 前端法兰厚度 | 0.725 | 18.415 |
| 前段外径 | 2.000 | 50.800 |
| 前法兰外径 | 1.960 | 49.784 |
| 后支杆外径 | 1.500 | 38.100 |
| 主供气孔径 | 1.000 | 25.400 |
| 前端内部扩张半角 | 5° | 5° |
| 外部过渡角/圆角 | 10° / R1.000 | 10° / R25.400 |

最终装配 DWG 1168299、外部供气管、线束、胶带及 TSP 到 sting 的尺寸链未公开，因此 sting 未进入当前主 STL。

### 7.4 表面测点

- 167 个 ESP：前体 118，后体 49。
- 9 个 Kulite 安装位；7 个前体传感器正常，后体 P126/P127 装配损坏。
- Table A-1 共 176 行。
- 典型前体测压孔径 `0.020 in=0.508 mm`；柔性压力管 `0.040 in=1.016 mm`；Kulite 外径 `0.0625 in=1.5875 mm`。

全部坐标见 [`surface_ports_table_a1.csv`](surface_ports_table_a1.csv)。CSV 的 `r/theta/x` 来自 NASA 表中已舍入值，`y/z` 是由这些舍入值重新计算的派生列。微孔和传感器槽默认不切入 OML，避免 IBM 出现与主流场无关的小特征。提取 CFD 压力时应沿流体侧法向做一致的小偏移，不能直接用固体表面的最近单元值。

## 8. 交付文件及各自用途

### 8.1 CAD 母版

```text
outputs/cad/test1853_run165_single.{brep,step}
outputs/cad/test1853_run262_263_tri_n3at120.{brep,step}
outputs/cad/test1853_run262_263_tri_n3at240.{brep,step}
outputs/cad/test1853_quad_n3at120.{brep,step}
outputs/cad/test1853_quad_n3at240.{brep,step}
```

STEP/BREP 是重新三角化的母版。正式细网格 STL 应从这些实体或参数生成器重新输出，而不是对旧 STL 反复 subdivide 后冒充 CAD 精度。

### 8.2 mm 闭合 STL

```text
outputs/stl/test1853_run165_single.stl
outputs/stl/test1853_run262_263_tri_n3at120.stl
outputs/stl/test1853_run262_263_tri_n3at240.stl
outputs/stl/test1853_quad_n3at120.stl
outputs/stl/test1853_quad_n3at240.stl
```

这些是 `outputs/stl/` 下的 mm 闭合母版。mm 制计算可直接读它们；SI 制计算应直接读第 8.4 节的闭合米制副本。

### 8.3 开放 patches 与参考面

```text
outputs/patches/<model>/body_wall_excluding_throat_caps.stl
outputs/patches/<model>/nozzle_N_throat_cap.stl
outputs/reference/<model>/nozzle_N_virtual_exit_REFERENCE_ONLY.stl
outputs/reference/<model>/physical_lip_reconstruction.csv
```

开放 patch 只用于离线标识、面积核验、可视化或支持 named-patch 的其他网格器。当前 Cerisse `ib.filename` 不能把它们作为独立几何加载，也不能在完整 STL 之外再重复加载。

### 8.4 SI master STL

`outputs/ibm_prepared/*_master_SI/` 包含 `phi=0`、TSP 位于 `(0,0,0) m` 的米制闭合副本、开放辅助 patches、IBM audit 和精确 throat metadata。正式算例若 TSP 不在原点或要应用 run roll，应使用第 12 节工具生成新的文件。

![Run 165 single 模型](outputs/figures/run165_single_isometric.png)

![Runs 247/262/263 tri 模型](outputs/figures/run262_263_tri_isometric.png)

![Runs 307/315 quad 模型](outputs/figures/quad_isometric.png)

## 9. 正式 STL 面片数、拓扑与文件锁定

### 9.1 总体统计

| 模型 | 三角面 | 焊接后顶点 | 唯一边 | 文件字节 | 面积, mm² | 体积, mm³ | `x` bounds, mm |
|---|---:|---:|---:|---:|---:|---:|---|
| Run 165 single | 353,042 | 176,523 | 529,563 | 17,652,184 | 123,985.403422 | 3,199,916.464716 | 2.441268921 – 268.137634277 |
| Tri N3@120 | 360,834 | 180,419 | 541,251 | 18,041,784 | 124,499.895386 | 3,198,241.721384 | 1.630115390 – 268.137634277 |
| Tri N3@240 | 360,872 | 180,438 | 541,308 | 18,043,684 | 124,499.895826 | 3,198,241.718716 | 1.630115390 – 268.137634277 |
| Quad N3@120 | 364,826 | 182,415 | 547,239 | 18,241,384 | 124,772.973909 | 3,197,304.082031 | 2.441268921 – 268.137634277 |
| Quad N3@240 | 364,880 | 182,442 | 547,320 | 18,244,084 | 124,772.974917 | 3,197,304.074322 | 2.441268921 – 268.137634277 |

三者的 `y/z` bounds 均为 `−63.5…+63.5 mm`。STL 体积相对 OCC 母版少约 `0.0109%`，与曲面内接三角化一致。

mm master SHA-256：

```text
Run165      f01949508523ee76fc034ed43b6fe4e99ee8626f1e69bf8066a38e9def815ca4
Tri N3@120  3c29d3706ab5488ca2740649c22e74ff02c760eca8ab439146928b2c6a8b0999
Tri N3@240  197264bae31a9c9e9caf29e96eba5e2ccee3c16100603b29b5dc144310a72a3d
Quad N3@120 6ff4a44aeddf2c59d0302eac6622600f41a55508f192164876a06cc2d8fd7d52
Quad N3@240 08fed36ce01022b09f01d1667a620a299af0b747cf814562f6ed8c086b964969
```

SI master SHA-256：

```text
Run165      f8d8095fa4d093fe421c1c6cc71712a1709839ddfeefd116cb8b91d91dd9df82
Tri N3@120  2a320b2c3348ef89115459dbb8d9be10e96ce685d57cb46da22fbd7bcbd1b7d4
Tri N3@240  0a9aeef093385a99848e59a94a4bf59813e55d6071eca98ffddc743c0c80b0ae
Quad N3@120 f56ad55ec5f8ff0c5c528cbb10e0779ee45bf74a75142874e3d2f198d6acf6b8
Quad N3@240 3dd8ed82b694b1fa87b7d6f0013ab5d7fa89da0812a4d61e7589f50b9c54541d
```

缩放、滚转或平移后字节和 hash 必然改变；每个计算序列必须记录最终实际输入文件的 hash。

正式生成器已把 `reproducible_num_threads=1` 写入参数并固定 Gmsh 单线程表面网格。旧三套保持原 SHA；新增两套 quad 也以正式全分辨率独立复生核对 SHA。当前参数文件/生成器 SHA-256 分别为：

```text
geometry_parameters.yaml  d2d388e330ecd7123a227c1f446f474c844d7aa03bc291600b3f2cc89da3fabe
build_test1853.py          2af592846540aafbb3eb8d8a209b564fe1b7486df42ef92f00497398a27bcc63
```

### 9.2 拓扑检查

五套模型全部满足：

- 二进制一阶三角形 STL；一个闭合连通体；Euler 特征数 2；正体积。
- watertight、winding consistent、整体外法向一致。
- boundary edge、non-manifold edge、inconsistent edge 均为 0。
- degenerate、weld-collapse、duplicate face 均为 0。
- stored normal 的 zero/mismatch/reversed 均为 0。
- 相对法向闭合残差分别约 `4.62e-16/9.28e-16/4.89e-16`。
- IBM `1e-7` 输入单位坐标量化审计全部 PASS。
- GTS `gtscheck` 三者返回 0：可定向 manifold 且未发现自相交。

IBM JSON 中的 `self_intersection_checked=false` 只说明 Python 审计器本身不做该项；不能解读成“存在自交”。独立 GTS 已完成并通过。

## 10. 三角面形状、边长与逐曲面分布

### 10.1 边长分位数

| 模型 | min | p1 | p5 | median | p95 | p99 | max | mean | 单位 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Run165 | 0.230564 | 0.317215 | 0.385693 | 0.400001 | 2.000003 | 2.000569 | 2.863651 | 0.683415 | mm |
| Tri N3@120 | 0.240679 | 0.315021 | 0.381665 | 0.400001 | 2.000003 | 2.000414 | 2.863651 | 0.677134 | mm |
| Tri N3@240 | 0.228252 | 0.314956 | 0.381328 | 0.400001 | 2.000003 | 2.000412 | 2.863651 | 0.677080 | mm |
| Quad N3@120 | 0.230564 | 0.314789 | 0.381284 | 0.400001 | 2.000003 | 2.000351 | 2.863651 | 0.674083 | mm |
| Quad N3@240 | 0.228252 | 0.314577 | 0.380794 | 0.400001 | 2.000003 | 2.000350 | 2.863651 | 0.674010 | mm |

前部 `x<=35 mm` 的最大边分别为 `0.613127/0.687050/0.678245 mm`，p95 约 `0.407–0.409 mm`；后体 `x>=60 mm` 的 median 约 `2.0 mm`，最大 `2.863651 mm`。因此 Gmsh 的 `front_target=0.4 mm` 和 `aft_target=2.0 mm` 是目标尺寸，不是硬最大边。

### 10.2 透明质量指标

采用

```text
q = 4 sqrt(3) A / (l1^2+l2^2+l3^2)
```

其中等边三角形 `q=1`：

| 模型 | q min | q p1 | q p5 | q mean | q<0.90 | q<0.95 | 最小内角 | 最大内角 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Run165 | 0.659655 | 0.897748 | 0.954481 | 0.993886 | 1.070% | 4.467% | 30.905° | 113.767° |
| Tri N3@120 | 0.649353 | 0.896226 | 0.952722 | 0.993576 | 1.127% | 4.658% | 30.905° | 114.852° |
| Tri N3@240 | 0.659655 | 0.896527 | 0.952878 | 0.993595 | 1.112% | 4.637% | 30.905° | 113.767° |
| Quad N3@120 | 0.649353 | 0.896443 | 0.952533 | 0.993555 | 1.116% | 4.675% | 30.905° | 114.852° |
| Quad N3@240 | 0.659655 | 0.896792 | 0.952582 | 0.993561 | 1.104% | 4.667% | 30.905° | 113.767° |

没有最小角小于 30° 或最大角大于 120° 的面；99% 面的最长边/最短边不超过约 1.44，最大约 1.935。结论是三角形没有针状或退化质量问题，精度瓶颈是面片相对目标流体网格太大。

GTS 自带另一种 face-quality 定义，其 min/mean 为 `0.839152/0.997606`、`0.833233/0.997485`、`0.839152/0.997493`。不要把它与上面的 `q` 数值直接混算。

### 10.3 逐曲面面片数

| 曲面 | Run165 faces | Tri N3@120 faces | Tri N3@240 faces |
|---|---:|---:|---:|
| 鼻球 | 1,766 | 3,671 | 3,671 |
| 主体 70° 锥面 | 181,407 | 175,679 | 175,693 |
| 肩 R2.54 torus | 18,520 | 18,520 | 18,520 |
| aft cylinder | 138,036 | 138,036 | 138,036 |
| aft 平 cap | 7,396 | 7,396 | 7,396 |
| N1 divergent cone / cap | 5,441 / 476 | — | — |
| N2 divergent cone / cap | — | 5,313 / 484 | 5,313 / 484 |
| N3 divergent cone / cap | — | 5,334 / 503 | 5,326 / 509 |
| N4 divergent cone / cap | — | 5,414 / 484 | 5,434 / 490 |

quad 的各曲面就是 single 的 N1 与相应 tri mapping 的 N2/N3/N4 同时存在；正式总面数分别为 364,826/364,880，四个独立 throat cap 面数为 N1 `476`、N2 `484`、N3 `503/509`、N4 `484/490`。完整逐曲面归属以 [`build_audit.json`](outputs/reports/build_audit.json) 和独立 patch 审计为准。

## 11. 曲率、CAD–STL 位置误差与法向误差

### 11.1 解析 CAD 曲率

下表给曲率绝对值；符号取决于法向约定。

| 曲面 | 主曲率，1/mm | 最小有限曲率半径 |
|---|---|---:|
| 鼻球 R25.4 | `k1=k2=0.0393701` | 25.4 mm |
| 70° 锥 | meridional `0`；circumferential 从 `0.0393701` 降至 `0.0055317` | 25.4 mm |
| 肩 R2.54 toroidal fillet | meridional `0.3937008`；circumferential `0.0055317–0.0157480` | 2.54 mm |
| aft cylinder R63.5 | `0` 与 `0.0157480` | 63.5 mm |
| N1 divergent cone 喉部 | `0` 与 `0.305355` | 3.27488 mm |
| N2 divergent cone 喉部 | `0` 与 `0.304505` | 3.28402 mm |
| N3 divergent cone 喉部 | `0` 与 `0.300037` | 3.33293 mm |
| N4 divergent cone 喉部 | `0` 与 `0.303885` | 3.29072 mm |
| throat/aft 平面 | `0,0` | 无穷 |
| lip、throat edge、aft edge | 无唯一曲率 | 锐边，不得填写 0 或有限值 |

STL 每个 facet 本身是平面，面内曲率为 0；上述曲率属于 CAD/解析参考面。当前通用 STL IBM 不从 faceted STL 自动恢复生产级曲率，不能把该表直接当成 solver 已使用的 shape operator。

### 11.2 CAD–STL 弦差

对每个三角形三条边的中点和质心投影到源 OCC 面：

| 模型 | 样本数 | max sampled error, mm | RMS, mm |
|---|---:|---:|---:|
| Run165 | 1,412,168 | 0.014863965 | 0.00220618 |
| Tri N3@120 | 1,443,336 | 0.014863965 | 0.00222741 |
| Tri N3@240 | 1,443,488 | 0.014863965 | 0.00222720 |
| Quad N3@120 | 1,459,304 | 0.014863965 | 0.00223693 |
| Quad N3@240 | 1,459,520 | 0.014863965 | 0.00223643 |

这是高密度抽样误差，不是严格 Hausdorff 上界；它只衡量“当前重建 CAD 到当前 STL”，不衡量“当前重建 CAD 到未公开 NASA 原始 CAD”，也不是 IBM 解误差。

### 11.3 平滑面 facet-normal 误差

排除锐边邻域后，facet normal 相对 OCC 面质心法向的 p99/max 角误差约为：

| 曲面 | Run165 p99/max | Tri N3@120 p99/max | Tri N3@240 p99/max |
|---|---|---|---|
| 鼻球 | 0.149° / 0.328° | 0.144° / 0.427° | 0.144° / 0.427° |
| 主体锥面 | 0.105° / 0.268° | 0.106° / 0.245° | 0.105° / 0.340° |
| 肩圆角 | 2.038° / 2.974° | 2.038° / 2.974° | 2.038° / 2.974° |
| aft cylinder | 0.288° / 0.681° | 0.288° / 0.681° | 0.288° / 0.681° |
| nozzle cones | 1.114° / 1.672° | p99 1.108–1.153° / max 1.686–1.926° | p99 1.121–1.142° / max 1.683–1.806° |

肩 R2.54 和喷管喉部附近是当前最紧曲率区，法向离散误差也最大。锐唇/喉缘没有唯一解析法向，不能包含在平滑面误差上界中。

### 11.4 喉面面积离散

| Nozzle | STL cap area, mm² | NASA as-built `At`, mm² | STL 偏小 |
|---:|---:|---:|---:|
| 1 | 31.358923 | 31.441873 | 0.2638% |
| 2 | 31.535418 | 31.618646 | 0.2632% |
| 3 | 32.479854 | 32.561870 | 0.2519% |
| 4 | 31.661784 | 31.745098 | 0.2624% |
| tri total | 95.677057 | 95.925615 | 0.2591% |

这个 `0.26%` 是 STL 喉面几何控制面积误差，必须进入误差预算，但不能直接宣称当前 IBM 的质量流量/推力也恰好低 `0.26%`。当前 GPIBM callback 规定的是表面 primitive state，守恒通量实际发生在笛卡尔 fluid/solid crossing faces；最终流量还受 crossing 离散、ghost state、Riemann flux、网格相位和 AMR 覆盖影响。代码目前没有“按 STL cap 总面积自动归一化到目标总流量”的现成接口。若必须精确控制总流量，需要实现状态调节/反馈或专用通量边界，并增加按 first-hit throat 分类的 crossing mass/momentum-flux 诊断。正式重网格仍建议把 cap 几何面积误差压到 `<0.05%`。

## 12. 在 Cerisse IBM 中的正确使用方法

### 12.1 本 Test 1853 构型只加载一个闭合实体

三维读取器直接使用 `ib.filename` 的原始坐标，默认强制每个输入几何 closed/watertight。接口本身接受文件数组，但本 Test 1853 构型应只加载所选构型的一个完整闭合 STL，不能再叠加其开放 patches。正确形式为：

```text
ib.filename = /absolute/path/to/test1853_..._SI.stl
```

并在算例 `ibmparm_t` 中使用：

```cpp
static constexpr bool interior_is_solid = true;
```

`interior_is_solid` 是 `prob.h` 的编译期常量，不是可靠的 inputs 级缩放/变换开关。

禁止：

- 把 `body_wall_excluding_throat_caps.stl` 作为 `ib.filename`；它是开放面。
- 把各 `nozzle_N_throat_cap.stl` 分别加入 `ib.filename`；它们也是开放面。
- 同时加载完整 STL 和 patches；这会重叠/重复几何。
- 用 `ib.skip_validation=1` 掩盖开放或损坏模型。

STL 不保留 CAD physical-group 名，当前 SurfElem 只有整体 geometry index。喷注 patch 必须在同一个闭合 STL 上通过 body/world 坐标和法向识别。

### 12.2 单位、缩放和生成 SI 文件

Cerisse 当前静态 STL 读取路径没有 `ib.stl_scale`，transform 默认 identity。若域、气体常数和状态使用 SI，不能把 mm master 直接加载。

示例：将 Run 262 的 N3@120 版本滚转并把 TSP 放到 `(0.275,0.275,0.275) m`：

```bash
cd /home/qiaoj/testcerisse/cerisse
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration tri_n3at120 --run 262 \
  --tsp-m 0.275 0.275 0.275
```

其他有效组合：

```bash
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration single --run 165 --tsp-m X Y Z

.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration tri_n3at240 --run 263 --tsp-m X Y Z

.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration tri_n3at120 --run 247 --tsp-m X Y Z

.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration quad_n3at120 --run 307 --tsp-m X Y Z

.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration quad_n3at240 --run 315 --tsp-m X Y Z
```

工具会同步变换闭合 STL、辅助 patches 和 throat 中心，写出 metadata 与新 hash。metadata 记录参数文件/准备器 SHA；指定 `--run` 时嵌入该 run 的完整核定条件。Run247 只接受 tri，Runs307/315 只接受 quad，错误组合在写文件前失败关闭。IBM 的 `1e-7` 顶点量化在米制下为 `1e-7 m=0.1 micrometre`；当前缩放后最短边约 `2.28e-4 m`，没有 weld-collapse 风险。

正式高分辨率重三角化后，必须用新增的 `--input-stl` 显式指定新文件；否则默认路径仍是当前约 36 万面的标准 master。自定义高分辨率 STL 必须已经是 mm/TSP body frame，并与所选 `--configuration` 的 throat 位置一致：

```bash
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration tri_n3at120 --run 262 \
  --input-stl /absolute/path/to/highres_tri_n3at120_mm.stl \
  --tsp-m X Y Z --name highres_tri_run262_SI
```

若还有与高分辨率实体匹配的开放辅助 patches，可另传 `--patch-dir`；否则工具会有意不复制旧 patches。运行后必须核对 metadata 中的 `closed_ibm_stl.source_sha256`、facets 和 bounds，防止误转换旧低面数文件。

### 12.3 throat patch 识别

在 master frame，对喷管 `j` 使用：

```text
abs(x - x_t[j]) <= eps_x
sqrt((y-y_c[j])^2 + (z-z_c[j])^2) <= Rt[j] + eps_r
n_x approximately -1
```

滚转/平移后必须读取生成目录里的 metadata world coordinates。不要只用三维球距离 `sqrt((x-xc)^2+(y-yc)^2+(z-zc)^2)<R`，否则可能误选相邻 divergent wall。

`eps_x/eps_r` 应略大于 STL/交点数值舍入但远小于相邻几何和流体 `dx`；选定后必须离线统计被识别的面积，并与第 11.4 节相符。不能仅凭“看起来有射流”接受 patch。

默认外表面、喷管锥壁和人工 aft cap 是 wall。throat cap 的固体到流体法向为 `−x`，回调中的 `QU` 是局部法向速度，所以向上游喷注应给正的局部 `QU`，而不是再次手写负全局 `u_x`；`QV=QW=0`。

### 12.4 第一阶段 throat 候选是由总状态推导的 M=1 静态状态

当前闭合面位于 as-built 喉面，而不是 virtual exit。第一阶段可用实验 reservoir 总压/总温构造理想气体 sonic 静态 throat 候选；`gamma=1.4` 时：

```text
T* = T0 * 2/(gamma+1)
p* = p0 * [2/(gamma+1)]^[gamma/(gamma-1)]
rho* = p*/(R T*)
un* = +sqrt(gamma R T*)      # 沿固体到流体局部法向，即全局 -x
```

应使用算例实际 EOS/组分的一致公式。不能直接复用旧二维 `srp_run165_axisym` 或旧 tri case 在 virtual exit 设置 `M≈2.94` 的做法；几何边界站位不同。

固定全部 primitive 的 sonic Dirichlet 只是待验证候选，并非当前已认证的 characteristic/reservoir inlet。更完整的候选应参考 canonical nozzle 的 reservoir-to-boundary 思路，根据流体侧法向状态和 `p0/T0` 构造边界静态状态。两者都必须通过实际 Cartesian crossing mass flux、面积—Mach 关系、总压/总温保持和固定几何网格序列验证。

若实验总质量流量是验收量，应先对三管施加一致 reservoir 条件，再测量 crossing-face 总流量。代码没有现成的 cap 总流量归一化开关；若需强制目标总量，要实现调节/反馈或专用通量边界。不要在没有一致性检查时同时强制总压、总温、单位面积通量和总质量流量四个互相可能过约束的量。

### 12.5 推荐 fixed-STL 配置

仓库 `docs/ibm_fixed_stl_pressure.md` 对一般固定 STL 的推荐是：

```cpp
interp_order = 2;
extrap_order = 2;
interp_order_surf = 2;
extrap_order_surf = 2;
ghost_layers = 1;
alpha = 1.0;
alpha_surf = 1.0;
pressure_closure = shock_aware_fluid_extrapolation;
pressure_sensor_low = 0.03;
pressure_sensor_high = 0.10;
interior_is_solid = true;
```

这些多为 `prob.h` 编译期设置。Euler baseline 使用 stationary slip wall；Navier–Stokes 必须使用物理问题需要的 wall model，但当前通用 3-D faceted STL 的 NS 热流/摩擦精度未认证。

本问题需要自定义 `ibm_user_t` 在同一个闭合面上区分 wall 与 blowing throat。仅在 wall class 中列出 `pressure_closure` 不会自动修正用户自行写入的 `q(1,QPRES)`。实现时应采用当前优先级最高、带 `std::integral_constant<EO>`、`disIM`、`n_valid` 和 `pressure_compat` 的完整 callback：非 throat 分支显式调用 `ibm_wall_pressure_surface<...>()` 取得 shock-aware wall pressure，throat 分支再覆盖 reservoir/sonic 候选状态。

同一个 wall class 含有非零 blowing patch 时，不应为整个类全局声明 `stationary_slip=true` 或 `stationary_no_slip=true`；这些是整个 wall model 的编译期能力声明，会启用只适用于完全静止、不透流壁面的路径。应按当前 `src/ibm/ibm_containers.h`、`ibm_walltypes.h` 和 solver dispatch 接口实现并写专门单元测试。

AMR production audit inputs 建议明确启用：

```text
ib.amr_support_buffer = 1
ib.amr_support_audit = 1
ib.amr_support_strict = 1
ib.boundary_face_audit_strict = 1
```

验证构建的强激波阶段审计：

```text
cns.stage_positivity = 1
cns.stage_positivity_verbose = 1
cns.strict_positivity = 1
cns.soft_positivity = 0
```

不要为让计算继续而打开 `soft_positivity`；它会修改守恒解。不要在 production SRP 中启用 `cns.afd_ibm_conservative_crossing`，仓库文档明确将其列为诊断模式。

## 13. 当前 STL 对 IBM 能支持到什么精度

### 13.1 必须区分四层误差

不能把“IBM 精度”写成单一的 `0.014864 mm`：

1. **实验重建误差**：原始 CAD、外围插入深度等未公开，当前没有可量化上界。
2. **CAD 到 STL 位置/法向误差**：当前最大 sampled chord error `0.014864 mm`，以及第 11.3 节 facet-normal 误差。
3. **边界面积误差**：当前 throat cap 约 `0.25–0.26%`。
4. **IBM 空间/时间/激波误差**：由最细 `dx`、AMR 支撑、闭合、通量、时间窗口和实际收敛序列决定。

### 13.2 当前流体网格对比

最小喉径为 `6.32714 mm`，全局 STL 最大边为 `2.863651 mm`，最坏前体最大边约 `0.68705 mm`。

| `dx_min`, mm | 直径 127 mm 的格数 | 最小喉径格数 | global max-edge/dx | front max-edge/dx | chord/dx | 推荐 production audit gate |
|---:|---:|---:|---:|---:|---:|---|
| 0.500 | 254 | 12.65 | 5.73 | 1.37 | 0.0297 | FAIL |
| 0.333 | 381 | 19.00 | 8.60 | 2.06 | 0.0446 | FAIL |
| 0.250 | 508 | 25.31 | 11.45 | 2.75 | 0.0595 | FAIL |
| 0.200 | 635 | 31.64 | 14.32 | 3.44 | 0.0743 | FAIL |
| 0.125 | 1016 | 50.62 | 22.91 | 5.50 | 0.1189 | FAIL |

仓库 fixed-STL 文档当前采用的推荐 production audit gate 是 `global max edge / finest dx <= 0.5`。它不是 solver 运行时硬限制，也不是满足后就足以证明解的精度。当前模型只有在 `dx>=5.727302 mm` 才通过这一推荐 gate，但此时最小喉径仅约 1.1 个单元，没有喷管解析意义。因此：

- 当前 STL 可以用于文件读取、inside/outside、边界回调、稳定性和初步流程测试。
- 当前 STL 不能作为亚毫米正式 fixed-STL 网格收敛几何，也不能据一次结果声称定量复现实验。
- `chord/dx<5%` 在约 `dx>=0.30 mm` 单项成立，不会覆盖更严格的推荐 edge gate，不能拿它替代完整 production audit。

### 13.3 正式重三角化目标与成本

| 目标最细 `dx`, mm | 推荐 gate 对应 `hmax<=0.5dx`, mm | 等边满铺理论最少面数，tri 表面积估算 |
|---:|---:|---:|
| 0.500 | 0.2500 | 约 4.60 million |
| 0.333 | 0.1665 | 约 10.37 million |
| 0.250 | 0.1250 | 约 18.40 million |
| 0.200 | 0.1000 | 约 28.75 million |
| 0.125 | 0.0625 | 约 73.61 million |

这是所有三角形刚好达到推荐最大边且为等边形时的理论下限；实际网格会更多。只局部细化前体/喷管可显著降成本，但它偏离当前全局推荐 gate，必须另外证明局部法向、交点、载荷和网格收敛。

位置误差目标建议：

- `dx_min=0.25 mm`：chord error 目标 `<=0.005–0.010 mm`，cap area error `<0.05%`。
- `dx_min=0.20 mm`：chord error 目标 `<=0.005 mm`。
- `dx_min=0.125 mm`：chord error 目标 `<=0.002–0.003 mm`。

对一个流体网格收敛序列，只生成一次满足最细 `dx` 的 STL，然后在所有粗/中/细流体分辨率使用完全相同的几何文件、点/面连接和 SHA-256。不能让 STL 随 fluid `dx` 同时变化，否则无法区分几何收敛和流体收敛。

### 13.4 当前 solver 认证边界

仓库 `docs/IBM_CURRENT_STATUS.md` 当前明确列出：

- 一般 faceted STL 曲率：未认证。
- 3-D 曲率、锐边、窄缝、NS 热流/摩擦和 AMR 物理解收敛：未认证。
- smooth isentropic nozzle：ready for verification，尚未 certified。
- 含激波 nozzle：NO-GO，尚缺保守且 positivity-preserving 的生产更新。

因此这套几何可以用于代码开发、稳定性测试、边界实现和网格研究；在 solver 认证矩阵完成以前，不能仅凭 3-D IBM 结果宣布实验定量复现成功。

## 14. 推荐 IBM 验证顺序与验收量

1. **几何审计**：锁定最终 SI STL hash，确认 facets、bounds、闭合、外向、GTS 自交、`hmax/dx`、cap area 和 throat metadata。
2. **无流几何 bring-up**：检查 inside/outside、BVH、所有 throat 命中面积、wall/cap 分类和 AMR support audit；禁止用 skip validation。
3. **无喷流外流**：先关闭 throat blowing，检查 Mach 4.6 外流、模型滚转、攻角符号、人工 aft cap 影响和表面压力提取。
4. **Run 165 单喷管**：第一阶段施加由实验总状态推导的 sonic 静态候选，随后验证 reservoir-to-throat 边界；核对实际 crossing mass/momentum flux、总量、对称性和通量记账。
5. **Run 247/262/263 三喷管**：分别测试 N3@120/N3@240，对三管施加一致 reservoir 条件，让不同 as-built `At` 自然产生流量；与面积比例和实验总量核验，并确认三个 patch 没有误选 divergent wall。
6. **Run 307/315 四喷管**：确认 N1/N2/N3/N4 四个 throat 均被唯一分类；主构型对比先做 Run247/Run307 的近同 `C_T` 条件，再把 Run315 作为低推力 quad，不得将同一总压或总流量从 tri 无条件复制给 quad。
7. **固定几何网格序列**：同一 hash 上做至少三档 `dx`，比较同一物理时间或已证明统计平稳的同一时间窗，不能比较相同步数。
8. **实验对比**：按 176 个 Table A-1 坐标提取 `Cp`，并按实际 run 的完整攻角序列、roll 和构型比较。

每档至少报告：

- 完整 `Cp` 分布、最大/迎风停滞压力、非正压力计数和最小压力。
- 无喷注阶段报告 pressure drag、横向对称载荷和不透流 wall 的 surface/CV force closure。
- 有喷注阶段分别报告：不透流实体壁面的压力/黏性载荷、throat 控制面的压力项与喷注动量流、全部 crossing-face 数值反力，以及包含发动机供气源项的外控制体平衡。此时普通 impermeable-wall 的 surface/CV 差只能作为离散通量记账诊断，不能直接当成物理 body-load 一致性门槛。
- image-point 质量、shock sensor、fallback measure fraction。
- 按 first-hit throat 分类的 Cartesian crossing mass/momentum flux、总温/总压一致性和 N2/N3/N4 分配；不能用 STL cap 面积乘 primitive state 代替实际数值通量。
- 相同物理时间、几何 SHA-256、facets、bounds、`hmax/dx` 和 AMR level 分布。
- 人工 aft cap 的力是排除、单独报告还是包含；不得静默混入实验模型载荷。

气动归一化使用 NASA `Sref=0.012667329504 m²`、`Lref=0.127 m`。外围真实 lip 面积和 virtual `Ae` 不得混作同一归一化面。

## 15. 公开资料无法消除的限制

以下内容仍是 **[UNKNOWN]**，后续对话不能用无来源精确数字填空：

- 原始 forebody/insert STEP、IGES 或逐点 CAD-defined surfaces。
- 四枚 nozzle virtual-exit 的真实全局 `x` 与外围精确插入深度。
- 各枚 as-built 喷管完整轴向 contour、所有 blend 切点和制造公差。
- Nozzle 3/4 与 `120/240°` 的明确物理编号对应。
- 最终装配 DWG 1168299、sting 相对 TSP 站位、外部管线/线束/胶带。
- manifold 三维弯曲支路、接缝、局部台阶和装配间隙。
- 喷管出口边界层、速度剖面与湍流强度；这些本来就不是几何表格参数。

因此无法给“当前 CAD 到真实实验 OML”的严格毫米误差上界。STL 的 `0.014864 mm` 只是在当前解析/假设 CAD 内部的三角化抽样误差。

## 16. 复现与审计命令

CAD/mm-STL 全量生成，不要带单配置 `--config`：

```bash
cd /home/qiaoj/testcerisse/cerisse
.venv/bin/python geometry/test1853/build_test1853.py
MPLCONFIGDIR=/tmp/matplotlib-test1853 \
  .venv/bin/python geometry/test1853/render_test1853.py
.venv/bin/python geometry/test1853/validate_test1853.py
.venv/bin/python geometry/test1853/analyze_test1853_mesh_quality.py
```

生成/刷新报告所列的五个 `phi=0`、TSP 原点 SI master：

```bash
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py --configuration single
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py --configuration tri_n3at120
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py --configuration tri_n3at240
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py --configuration quad_n3at120
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py --configuration quad_n3at240
```

`prepare_test1853_ibm.py` 会写 metadata；随后应对每个 SI master 运行 `tools/audit_ibm_geometry.py` 并把 JSON 写回对应目录。`outputs/reports/ibm_audit_*.json` 则是对 mm masters 的独立审计，两组都不能靠生成器自动假定为最新。

针对实际 SI 文件和最细网格的正式 IBM 审计：

```bash
.venv/bin/python tools/audit_ibm_geometry.py FINAL_SI.stl \
  --kind stl --require-outward-input --weld-tolerance 1e-7 \
  --dx 0.00025 --max-edge-over-dx 0.5 \
  --json FINAL_SI_ibm_audit.json
```

上例的 SI `dx=0.00025 m`，即 `0.25 mm`。审计脚本直接比较原始坐标，不做单位换算；若误填 `--dx 0.25`，它会按 `0.25 m` 给出虚假宽松结果。

`audit_ibm_geometry.py` 不检查三维 STL 自相交。每次新重网格还必须用 GTS/CGAL 等独立路径检查；当前环境可用的规范化 GTS 命令为：

```bash
.venv/bin/python -c 'import trimesh,sys; m=trimesh.load_mesh(sys.argv[1],process=False); sys.stdout.write(m.export(file_type="stl_ascii"))' FINAL_SI.stl \
  | stl2gts | gtscheck -v >/dev/null
```

返回码 0 才表示可定向 manifold 且没有检测到自相交。

详细现有证据：

- [`outputs/reports/build_audit.json`](outputs/reports/build_audit.json)
- [`outputs/reports/independent_validation.md`](outputs/reports/independent_validation.md)
- [`outputs/reports/independent_validation.json`](outputs/reports/independent_validation.json)
- [`outputs/reports/mesh_quality_detailed.json`](outputs/reports/mesh_quality_detailed.json)
- `outputs/reports/ibm_audit_*.json`

## 17. Cerisse `read_geom`/BVH 兼容性、容量与通信审计

### 17.1 审计边界和直接结论

本节针对 2026-07-19 工作区中实际存在的 IBM 源码进行只读审查，重点是：

- [`ibm_backend_bvh.h`](../../src/ibm/ibm_backend_bvh.h)：STL/OFF 读取、焊点、inside/outside、几何 cache；
- [`ibm_bvh.h`](../../src/ibm/ibm_bvh.h) 与 [`ibm_bvh_defs.h`](../../src/ibm/ibm_bvh_defs.h)：BVH2/BVH4、Morton、查询和数据布局；
- [`ibm_solver_io.h`](../../src/ibm/ibm_solver_io.h)：`read_geom()`、surface gather 和 VTP；
- [`ibm_solver.h`](../../src/ibm/ibm_solver.h) 与 [`ibm_containers.h`](../../src/ibm/ibm_containers.h)：逐面 cache、surface SoA、所有权和 CPU/GPU 存储。

审计时仓库 HEAD 为 `dc35646d8eb90aaa362a3292a9d702e30c9d5af5`，但 IBM 文件含用户尚未提交的修改，所以不能只用 HEAD 标识代码。六个关键文件的 SHA-256、测试可执行文件 hash、冒烟测试结果和内存公式已固化在 [`outputs/reports/ibm_code_compatibility_audit.json`](outputs/reports/ibm_code_compatibility_audit.json)。本次没有修改 IBM 源码。

结论分两层：

1. **当前约 35–37 万面模型满足静态几何契约并通过短时集成门禁。** 五个 SI master 都满足 reader 所需的闭合三角网格、外向 winding、单实体、有限坐标和 `1e-7 m` 焊接约定；tri N3@120 与 quad N3@120 已分别在 `IBM/cases/test1853_tri_3d`、`IBM/cases/test1853_quad_3d` 中通过 reader、闭合校验、BVH、marker、严格同层 GP/surface 支撑审计和一个极短 coarse step。它们仍只是 bring-up/capacity 结果，不能写成 SRP 物理认证；另一 mapping 与其余 Run 条件目前是同构配置和静态门禁，不得冒充已逐一长跑。
2. **第 13.3 节的千万级正式 STL 暂不具备容量兼容性。** 18.4M/28.75M 面模型会被“每 rank 全几何复制 + 全面片 surface SoA + rank-0 gather”首先卡住；这不是 STL 格式错误，而是当前数据布局和输出架构的容量上限。

一次实际 CPU 冒烟测试使用

```text
model: test1853_run262_263_tri_n3at120_master_SI.stl
ib.skip_validation=0
ib.plot_surf=0
MPI ranks=1，极短 1 step，plot/checkpoint 均关闭
```

结果为 exit code 0：reader 得到 `180419 vertices / 360834 faces`，BVH2 为 `721667` nodes，BVH4 为 `186071` nodes；BVH build `0.177329 s`，geometry cache `0.023136 s`，完整 `read_geom` `0.893544 s`。这证明该 STL 与测试二进制的实际接入兼容，但不是流动物理认证；测试二进制也不能代替用上述源码 hashes 重编译后的正式回归。

### 17.2 当前真实数据流

当前非 CGAL 路径不是“一次读入、全程序共享”，而是：

```text
每一个 MPI rank
  └─ 独立打开同一个 STL
      └─ 解析 3F 个 STL 顶点记录并焊成 TriMesh
          ├─ TriMesh vertices/faces（ManagedVector）
          ├─ BVH2（长期保留，inside/outside 使用）
          ├─ BVH4（长期保留，closest/segment 使用）
          ├─ InsideTester 再复制一次 vertices/faces
          ├─ LocalFrame + SurfElem（每个 facet 一份）
          └─ 若 ib.plot_surf=1
              ├─ surfImp + surfPhys 仍按全局 F 在每个 rank 分配
              ├─ Phase 2 每 rank 扫描全部 F，重计算 face owner
              ├─ 重计算只在本 rank owned faces 上继续
              └─ 输出时 MPI_Gatherv 全部记录到 rank 0
                  └─ rank 0 串行写 ASCII VTP
```

所以“多 rank”目前会分摊本地 surface 重计算和流体网格，却不会分摊基础 TriMesh、两棵 BVH、frame/cache 以及 surface SoA 的容量。一个节点上若多个 MPI rank 共用一块 GPU，这些逐 rank 数据还会在同一 GPU 地址空间竞争 UVM 容量和迁移带宽。

### 17.3 Test 1853 与 reader 的具体接口合同

- `ib.filename` 对本构型只加载一个完整闭合 SI STL；reader 接口虽支持多个文件，但开放的 `body_wall` 和 `throat_cap` 辅助 STL 不能重复加载。
- 必须用 `outputs/ibm_prepared/` 的 metre 文件或经 `prepare_test1853_ibm.py` 预变换的新文件。reader 没有输入级单位/scale；固定焊接量 `1e-7` 也是按原始坐标单位解释。
- binary STL normal 和 16-bit attribute 被忽略，法向由 triangle winding 重算；STL 内没有可供 wall callback 使用的 patch/material 名称。
- 当前 throat 必须用“平面 + 横向圆”判定：`abs(x-x_t)<=eps` 且 `sqrt((y-y_t)^2+(z-z_t)^2)<=R_t+eps`，再辅以法向朝 `-x`。旧 `IBM/cases/srp_tri_NS/prob.h` 的三维球距离会选入相邻 divergent wall，与本模型不兼容。
- 更稳健的长期接口是在预处理几何格式中保存 `face_patch_id`，并把 sidecar 与 STL SHA-256、weld tolerance 和 face ordering 一起锁定；单独按面号写一个不带 hash 的文本文件不安全。
- 当前内部 `is_closed()` 只统计无向边出现两次，并不真正检查注释所说的相反有向边；全局 signed-volume flip 也不能修复局部 winding 错误。五个交付模型已经由仓库审计、Trimesh 和 GTS 独立证明 winding/闭合/自交合格，因此不会触发这个缺口；任意新 STL 仍必须执行第 16 节外部审计，不能只靠 `read_geom`。

### 17.4 当前逐 rank 内存：实测节点数下的可复核预算

以下按 3-D LP64、`Real=double=8 B`、`int=4 B` 计算；是 payload 下界，不包括 vector capacity、UVM 元数据、流场/AMR/GPStore、MPI buffer 和建树临时量。`ManagedVector` 是一份统一虚拟地址存储，并不意味着同一 payload 必然永久各占一份 host 与 device 物理内存；但每个 MPI 进程都会建立自己的逻辑分配，页面迁移/超额订阅也不能把它当成免费容量。当前 tri N3@120 的实际计数为 `F=360834, V=180419, BVH2=721667, BVH4=186071`。

| 常驻项目 | 数据布局 | 当前 tri / rank |
|---|---:|---:|
| `TriMesh` | `24V + 12F` | 8.259 MiB |
| `InsideTester` 重复网格 | `24V + 12F` | 8.259 MiB |
| BVH2 | 64 B/node | 44.047 MiB |
| BVH4 | 216 B/node | 38.329 MiB |
| `LocalFrame` | 72 B/face | 24.777 MiB |
| `SurfElem` | 40 B/face | 13.765 MiB |
| **几何常驻小计，`plot_surf=0`** | 实际节点数 | **137.435 MiB = 0.144112 GB/rank** |

若启用表面重建/输出：

- 旧 SRP tri 的 `E1/I1`：`surfImp=208 B/face`；加 production `surfPhys≈268 B/face` 和 7-Real conservative impulse `56 B/face` 后，当前 tri 合计约 **0.336 GB/rank**。
- 仓库 fixed-STL 指南的 `E2/I2`：`surfImp=1180 B/face`；同样计算，当前 tri 合计约 **0.687 GB/rank**。
- `CNS_IBM_VALIDATION` 的 `surfPhys` 字段再加约 `176 B/face`；若还持久化 E2/I2 visibility audit，再加约 `288 B/face`。其他 boundary-face/GP validation 容器尚未计入。

把当前实际 BVH4/face 比例外推到第 13.3 节“理论最少面数”，得到：

| `dx` 目标, mm | 理论最少 facets | `plot_surf=0` 几何 GB/rank | surface E1/I1 GB/rank | surface E2/I2 GB/rank |
|---:|---:|---:|---:|---:|
| 0.500 | 4.60M | 1.84 | 4.28 | 8.76 |
| 0.333 | 10.37M | 4.14 | 9.66 | 19.74 |
| 0.250 | 18.40M | 7.35 | 17.14 | 35.02 |
| 0.200 | 28.75M | 11.48 | 26.78 | 54.72 |
| 0.125 | 73.61M | 29.40 | 68.56 | 140.11 |

这些面数本身还是理想等边满铺下界，实际会更高。因而：

- A100 40 GB 上，18.4M + E2/I2 的逐 rank IBM payload 已约 35 GB，尚未放流场、AMR、GP 和临时量，判定 **NO-GO**。
- 28.75M + E2/I2 约 55 GB，即使 80 GB GPU 也不能在没有实测 peak budget 的情况下认为安全；加载/建树峰值可能先于稳态常驻量失败。
- `plot_surf=0` 会显著降内存，适合 reader/BVH bring-up，但它不能替代需要表面压力、载荷和诊断的正式实验比较。

加载峰值还高于上表：binary reader 预留至 `3F` vertices 并维护 `2F` 量级 hash；`is_closed()` 随后建立至多 `3F` edge hash，而且当前被调用两次；GPU build 同时存在 host primitive AABB、device AABB/Morton/parent/counter，BVH4 collapse 又把约 `2F` BVH2 拷回 host 并建立临时 BVH4。千万级运行必须先记录 host RSS、device used/free、UVM faults 和每个容器 capacity，不能只看最终 vector size。

### 17.5 `gatherSurfData` 的确定硬上限

当前 `gatherSurfData()` 把 POD 记录复制为 bytes 后做一次 `MPI_Gatherv`；count/displacement 仍是 32-bit `int`。DOUBLE 常规 layout 下：

| build | `SurfOut` | `INT_MAX/sizeof(SurfOut)` | 当前 tri 数据量 | 18.4M 数据量 |
|---|---:|---:|---:|---:|
| production | 约 192 B | 11,184,810 faces | 69.3 MB | 3.53 GB，必定超限 |
| validation | 约 368 B | 5,835,553 faces | 132.8 MB | 6.77 GB，必定超限 |

10.37M production 已约 1.99 GB，虽略低于单次 count 上限，rank 0 仍要承担接近 2 GB receive buffer、已有完整 SoA 和串行 ASCII 转换；不应视为可生产使用。分块 gather 只能绕开 `int`，不能解决 rank-0 内存和串行文件瓶颈。

### 17.6 已确认的代码热点与正确性风险

1. **闭合验证重复且内存不可控。** `read_mesh()` 已调用 `is_closed()`，`read_geom()` 随后再调用一次；两次都创建巨型 `unordered_map`。应只验证一次并返回结构化结果。大模型可使用 packed directed-edge 的 sort/scan，或读取与源 STL hash 绑定的离线审计 manifest；仍应抽样/全量核验 manifest 与当前 weld 后连接一致。
2. **`InsideTester` 不必要地拥有第二套网格。** 其 view 只需稳定指针，而 `geom_a` 在 tester 生命周期内稳定。改为明确的 non-owning view 可立即省约 `24 B/face`，即 18.4M 约 0.44 GB/rank、28.75M 约 0.69 GB/rank。
3. **BVH2/BVH4 同时常驻。** closest-point 与 segment visibility 已走 BVH4，inside parity 仍走 BVH2。若把 parity ray 移植到 BVH4并通过严格回归，即可释放 BVH2，约省 `128 B/face`：18M 约 2.30 GB/rank、29M 约 3.71 GB/rank。
4. **inside BVH2 栈可能静默漏节点。** 每线程固定 `128 int`，满时只是停止 push child，不报错。必须记录 build-time 最大 frontier/depth；超过容量时应显式 abort 或使用不会漏检的 overflow/stackless 路径。BVH4 的 64-entry fallback 也应在生产中打印并禁止意外退化为全 primitive scan。
5. **30-bit Morton 对千万级表面过粗。** 3-D 每方向只有 10 bit；大量相同 code 不会自动破坏正确性，但会影响树质量、CPU/GPU树形一致性和查询带宽。至少先用 `(Morton, primitive_id)` 确定性排序并记录 visits/depth；随后基准 64-bit Morton + radix sort/HLBVH treelet SAH。
6. **建树有明显 CPU↔GPU 往返。** primitive AABB 先在 CPU 形成再上传，LBVH 在 GPU 生成，collapse 又把整棵 BVH2 拉回 CPU，再把 BVH4 复制回 managed memory。大网格应直接在 device 从 face/vertex 生成 AABB/Morton，构建 wide BVH，或加载与 STL hash、代码版本、精度和 endian 绑定的预建 cache。
7. **surface ownership 的轻量阶段仍每 rank 全扫。** Phase 2/3 对全部 `F` 做 host centroid/FAB 查找；Phase 4 只处理 owned faces，但读取 forced-managed marker/SoA，容易造成 GPU→CPU→GPU 页面迁移。应把 locate、owner compaction 和可并行的 image-point 工作移到 device，并以 compact owned face list 驱动后续 kernel。
8. **surface 字段没有按功能 lazy allocate。** production `surfPhys` 的 conservative/load-consistent 数组即使某功能不用也逐面分配。先根据编译期/runtime 功能开关分配，再进一步把 heavy payload 改为 owned-face-only；保持 `global_face_id` 用于重组。
9. **CPU/GPU 退化面处理不一致。** CPU cache 路径在零面积时替换 centroid/area 后仍可能用零法向长度归一化，GPU 路径则写零 frame。正式 loader 应在建树前统一拒绝 non-finite、重复 index 和退化 triangle；本交付当前退化面数为零。
10. **ASCII VTP 是独立瓶颈。** rank 0 逐数写所有 vertices/connectivity/fields；千万级应改为每 rank binary appended VTP + `.pvtp`，或并行 HDF5/ADIOS2。device owned list 可直接 pack 到 pinned host staging，避免当前 `vector<SurfOut>` 再复制到 `vector<char>`。

### 17.7 分级修改建议

#### P0：任何 10M+ facet 生产尝试前必须完成

1. 启动时打印并门控：`sizeof`、F/V/BVH node 数、每容器 size/capacity、预测 steady/peak host/device bytes、device free memory、surface gather 上限；预算超限应在分配前 abort。
2. 合并闭合/方向校验，补 directed-edge、non-finite、index、degenerate 检查；保存 source/hash/weld 审计 manifest。
3. 去掉 `InsideTester` 的 vertices/faces 副本。
4. 对 `surfPhys` 可选字段 lazy allocate；随后把 `surfImp/surfPhys/visibility/CSR` heavy payload 改为 owned-face-only。
5. 以并行 binary surface output 取代单次 gather + rank-0 ASCII；在完成前，10M 级模型必须保持 `ib.plot_surf=0`，且不得假装已有表面载荷交付能力。
6. 明确 inside/BVH4 栈门槛和 overflow 行为，禁止静默漏节点或线性退化。

第 2、3 项预期不改变数学结果，但仍要回归 pointer lifetime 与 marker；第 4、5 项改变分布式索引/输出顺序，必须按 `global_face_id` 重排后逐字段比较。

#### P1：主要速度、CPU/GPU迁移和中期容量

1. GPU 化 `computeSurfIndices` 的 locate/compact/Phase 4，使用 device owned list；避免 host 全扫 managed 页面。
2. device-direct AABB/Morton/wide-BVH build，或引入版本化预建 cache；构建完成后释放所有 scratch。
3. 实现 BVH4 parity ray并释放 BVH2；closest/segment/inside 都使用一致、可审计的 tie-break。
4. 让 `InsideTester`、TriMesh、BVH 和只读 cache 使用明确的 device storage；CPU I/O 使用普通 host/pinned staging。若暂留 UVM，至少做 device preferred-location/read-mostly/prefetch，并禁止同 GPU 多 rank 重复驻留。
5. 压缩 cache：可先从 `LocalFrame` 只存 `normal+tangent1`、查询时叉乘 `tangent2`；评估从 triangle 重算 centroid/normal 与保留带宽的权衡。任何 float/压缩都必须做法向和载荷误差审计。
6. STL 只在一个 I/O rank 或每节点一个 rank 解析；以分块 broadcast/shared cache 分发已焊 vertices/faces，避免并行文件系统风暴。它只能降低启动 I/O，不能替代 owned-only 的逐 rank 容量修复。

#### P2：需要单独设计和物理认证

- 64-bit Morton、4–8 triangle leaf packet、treelet SAH、stackless/rope BVH4/BVH8。
- child AABB 使用 outward-inflated FP32 或 parent-relative quantization；必须用 `nextafter` 外扩并由 ray/segment/closest oracle 证明零漏检。
- 每 GPU 只保留 subdomain+halo 的细 BVH，配合全局 coarse sign/SDF/TLAS 做 inside/outside；要处理 AMR regrid、移动体和跨 rank 查询。
- 为 Test 1853 的球头、70°锥、R0.100肩、圆柱和喷管锥建立 analytic/CSG backend，以避免用千万面 STL 表示解析曲面。这会改变几何查询后端，必须保持当前假设、斜唇拓扑和 throat 分类，并重新做流体收敛，不是简单内存补丁。

### 17.8 通信优化不能只靠打开 GPU-aware MPI

仓库已有 A100 实测证据：

- [`T2_VERDICT.md`](../../IBM/gpu_a100_test/perf_ibm/T2_surfindices_scaling_A100/T2_VERDICT.md)：运动 E2/I2 案例中 `computeSurfIndices` 为 `12.01 s/40 steps`（1 GPU），2/4/8 GPU 相对为 `0.57/0.33/0.20`；它会扩展，但 8 GPU 时仍是 eflux 的 `4.8x`，约占步时 `15–19%`。因此小网格下不应为一个已不存在的“完全复制平台”盲目重写；千万级时 owned-only 的理由是容量，而 device 化的理由是绝对耗时与 UVM 迁移。
- [`REPORT_FILLED_20260707.md`](../../IBM/gpu_a100_test/REPORT_FILLED_20260707.md)：8xA100 案例 `ParallelCopy_finish + nowait` 约占 `65.31%`，1→8 GPU speedup `2.74x`；这是 AMR halo/同步主导，不全是 IBM BVH。Amazon OpenMPI 4.1.7 没有 CUDA support，`amrex.use_gpu_aware_mpi=1` 直接 MPI abort。

短期建议：

- 在已经测试的 AWS 栈保持 `amrex.use_gpu_aware_mpi=0`；换成确认 CUDA-aware 的 MPI/UCX 后，再做 1/2/8 rank bitwise parity 和性能 A/B，不能仅凭开关名称启用。
- 每 GPU 一个 MPI rank；surface pack 用 device compact + asynchronous D2H 到 pinned host，通信与计算重叠。
- 调整每 GPU cells、BoxArray 和 blocking/max-grid-size，减少 halo 面积/体积比；用 profiler 分开 `ParallelCopy_nowait`、`finish`、实际网络和等待。
- `amrex.the_arena_is_managed=0` 不会自动改变 IBM 显式 `ManagedVector` 和强制 managed marker；必须改具体容器/arena 才会减少 UVM。

### 17.9 每类优化的最低回归门

| 修改 | 必须保持/验证 |
|---|---|
| non-owning `InsideTester`、单次 validation | 同一 STL 的 V/F、bounds、closed/outward、inside markers、closest face ID 逐位一致；重建/移动生命周期无悬空指针 |
| BVH4 inside、stackless、Morton/tie-break | 随机点和 near-edge/near-vertex oracle；CPU/GPU、1/2/8 rank solid/ghost masks；segment first-hit `(fraction, prim_id)`；无 overflow/fallback |
| owned-only surface SoA | 按 global face ID 重排后 pressure/tau/IP/quality 一致；面积、法向闭合、逐几何力和 throat 分类一致；无漏 owner/双 owner |
| distributed binary output | 1/2/8 rank 合并后字段、face count、integrated load 与旧小模型 VTP 一致；文件可被 ParaView/后处理读取 |
| GPU `computeSurfIndices` | CPU/GPU mask 和 owner parity；AMR level/FAB mapping；regrid/restart；TinyProfiler、UVM faults、H2D/D2H bytes |
| cache 压缩/FP32 bounds | CAD/STL 法向与位置误差；零 BVH 漏检；`Cp`、压力/黏性力、first-hit flux 和固定时间网格序列 |

发布优化前至少保留三档测试：小型 CGAL/BVH oracle（若有 CGAL build）、当前 360k Test 1853、以及不能装入旧路径但可验证新容量路径的 synthetic multi-million mesh。对当前模型先做 reader/BVH bring-up 时可设 `ib.plot_surf=0`；需要表面数据时 360k 规模仍可用现路径。18M 以上在完成 P0 前不得进入长时间生产队列。

## 18. 可直接复制给负责 IBM 的对话

```text
请接手 /home/qiaoj/testcerisse/cerisse/geometry/test1853 的 NASA Test 1853 IBM 测试。
首先完整阅读 IBM_HANDOFF_REPORT_CN.md 和 geometry_parameters.yaml。

硬规则：
1. Run165 是中心 N1 single；tri 是 Runs247/262/263，中心 plug、外围 N2/N3/N4；quad 是 Runs307/315，中心 N1 与外围 N2/N3/N4 全启用。主同推力对比是 Run247 tri 与 Run307 quad。
2. NASA 未公开 N3/N4 对应120/240度，tri 与 quad 都必须保留两个 mapping 或做敏感性测试。
3. 本Test1853构型在ib.filename中只加载一个完整闭合SI STL；接口虽支持多个闭合几何，但开放body/throat patches只作离线标识，不能重复加载。
4. STL reader没有输入级scale。用prepare_test1853_ibm.py同步做mm->m、roll、TSP平移；正式高分辨率STL必须显式传--input-stl并核对metadata源hash。
5. throat patch 用平面x + 轴心横向半径 + nx约-1识别，不得用三维球距离。
6. 当前cap是throat，不是M约2.9 virtual exit；由p0/T0推导的M=1静态全状态只是第一阶段候选，必须验证reservoir边界和实际crossing flux。
7. 自定义mixed wall/throat callback的wall分支必须显式走shock-aware pressure helper；不能假定仅列pressure_closure就会自动修改用户写入的压力，也不能把含blowing的整个wall class声明成全局stationary wall。
8. 当前约360k-face STL只用于bring-up。仓库推荐production audit gate是全局hmax/dx_min<=0.5，并非solver硬限制或充分证明；按目标dx从STEP/BREP重新三角化、重审计，在整个fluid序列锁定同一SHA。
9. 报告按first-hit throat分类的Cartesian crossing mass/momentum flux、Cp、fallback、positivity和同一物理时间。有喷注时分开wall载荷、throat压力/动量和含供气源项的CV平衡；人工aft cap载荷排除或单列。
10. 当前仓库对一般3-D faceted STL curvature、shock-containing nozzle和NS载荷尚未认证，不得把一次稳定运行称为实验定量复现。
11. 先读第17节和outputs/reports/ibm_code_compatibility_audit.json。当前read_geom每rank完整读/建，InsideTester重复网格，BVH2/BVH4同时保留；18.4M面的E2/I2逐rank IBM payload已估约35GB，尚未包含流场，A100-40GB为NO-GO。
12. 现有surface MPI_Gatherv在production约11.18M faces、validation约5.84M faces达到32-bit byte-count硬上限；千万级前必须完成owned-only surface payload和分布式binary output，不能只做chunked gather。
13. 当前360k模型bring-up可先ib.plot_surf=0。每GPU只用一个MPI rank；Amazon OpenMPI 4.1.7实测保持amrex.use_gpu_aware_mpi=0，换通信栈后必须重新做1/2/8-rank parity再启用。
14. 第19节把表面h_STL与流体dx分开。目标体网格是L3前体/喷口/肩、L2前圆柱、L1中圆柱、L0后圆柱/人工cap；解析图不是实际AMReX hierarchy，必须由Test1853 plotfile复核。
15. 当前ib.amr_support_buffer会在每一级seed整个IBM界面；仅改user_tagging不能让后体保持粗网格。不要简单关buffer。所有user/flow tags与baseline support seed应服从统一body-frame local_max_level；若其他合法tag在粗表面附近创建child grid，support closure仍须局部提升补足halo。gradient只在完整valid-fluid stencil计算，并保持runtime dilation、amr_support_audit=1、amr_support_strict=1。
16. Run165 外流初值必须是几何体前方整张 yz 启动面：当前 x_front=-0.052 m；左侧为 Run165 M4.6，右侧为同p/T/rho静止气体。只有x-low持续注入来流，不能把几何体周围直接填成M4.6。喷管真实内通道的quasi-1D初值是单独的数值启动保护，不是预置外流。
17. AWS job862 的两步结果只是初始化/IBM/容量/性能门禁，最终时间8.8566e-8 s；启动波尚未到达模型，严禁称为弓形激波、喷流相互作用或Run165稳态已经建立。
18. job862的`regrid_int=1`会在AMR subcycling的eligible level/substep重复regrid；实测10次Amr::regrid、44次IBM rebuild。逐调用MPI-rank最大值之和占求解器`Run Time total`约29.33%；原作业旧analyzer仍打印OK，更新analyzer的离线复算才触发累计性能报警。正式计算必须降低regrid频率并继续记录GP/in-out/image-point开销；该集群关闭Slurm accounting，当前没有可信MaxRSS，不能把“未OOM”写成峰值内存认证。

在修改旧 IBM/cases/srp_tri 之前，先指出其中旧 trinozzle.stl、统一R=0.0028、旧中心坐标和三维球patch选择与本报告不兼容；不要直接复用这些硬编码。
```

## 19. 网格密度可视化与局部 IBM `dx` 分区建议

### 19.1 先区分两个完全不同的尺度

- `h_STL`：固定几何表面每个三角形的最大边长；影响几何近似、BVH 面数和 surface SoA。
- `dx`：Cerisse/AMReX 流体笛卡尔单元宽度；影响流场分辨率、CFL、AMR 单元数和通信。

两者不能用一张图冒充。当前 STL 的逐面 `h_STL` 热图、局部真三角线框和可交互 VTP 已生成：

- [`tri_surface_mesh_density_overview.png`](outputs/figures/tri_surface_mesh_density_overview.png)：全身、前体、肩部和底面同一色标；
- [`tri_surface_mesh_density_forebody_nozzles.png`](outputs/figures/tri_surface_mesh_density_forebody_nozzles.png)：三喷管/前体放大；
- [`tri_surface_mesh_density_local_details.png`](outputs/figures/tri_surface_mesh_density_local_details.png)：N2、R2.54 肩、后圆柱和 aft cap 的真实三角片；
- [`single_tri_surface_mesh_density_comparison.png`](outputs/figures/single_tri_surface_mesh_density_comparison.png)：Run165 single 与 Runs262/263 tri 前视同色标对比；
- [`tri_surface_mesh_density_axial.png`](outputs/figures/tri_surface_mesh_density_axial.png)：沿 x 的 1-mm 分箱 `p50/p95/max` 与单位面积面数；
- [`test1853_run262_263_tri_n3at120_surface_mesh_metrics.vtp`](outputs/visualization/test1853_run262_263_tri_n3at120_surface_mesh_metrics.vtp)：ParaView 可交互逐面查询 `hmax_mm/hmean_mm/area_mm2/centroid_x_mm/visualization_zone_id`；
- [`mesh_distribution_visualization.json`](outputs/reports/mesh_distribution_visualization.json)：上述分区的机器可读计数和统计。

当前 tri 的逐面最大边长实测如下。喷管分区按轴距、喉面和向内壁法向作可视化分类；正式边界语义仍以 throat sidecar patch/metadata 为准。

| 表面位置 | faces | 占比 | `h_STL` p50, mm | p95, mm | max, mm |
|---|---:|---:|---:|---:|---:|
| active nozzles | 17,828 | 4.94% | 0.4003 | 0.4467 | 0.5455 |
| forebody | 179,054 | 49.62% | 0.4000 | 0.4070 | 0.6871 |
| NASA R2.54 shoulder | 18,520 | 5.13% | 0.4017 | 0.4851 | 0.6131 |
| forward cylinder, 至 x=35 mm | 59,514 | 16.49% | 0.4000 | 0.4568 | 0.6131 |
| 35–60 mm size transition | 30,278 | 8.39% | 0.6982 | 1.7429 | 2.5166 |
| aft cylinder, x>=60 mm | 48,244 | 13.37% | 2.0000 | 2.0656 | 2.8637 |
| artificial flat aft cap | 7,396 | 2.05% | 2.0000 | 2.1485 | 2.7262 |

所以当前表面网格已经是前密后疏：Gmsh 请求值为 `x<=35 mm: 0.4 mm`、`35<x<60 mm` 线性过渡、`x>=60 mm: 2.0 mm`。人工底面只占约 2.1% 的 STL faces，单独进一步粗化它对 STL 总面数收益很小；但在三维流体 AMR 中降低整个后体/底面邻域的 level，体单元节省会很显著。

### 19.2 建议的 IBM 体网格分区

[`proposed_ibm_dx_distribution.png`](outputs/figures/proposed_ibm_dx_distribution.png) 是**解析设计图，不是运行产生的 AMReX hierarchy**。数值采用仓库现有 H100 基线的 `0.64 m / 384` 基础间距作为尺度参考，但旧案例的 `trinozzle.stl`、坐标和喷口 tagging 不能用于本 Test 1853。以 TSP/body frame 的 `x` 为准：

| 目标区 | AMR level | 示例 `dx`, mm | 相对 L3 单位体积 cell 数 |
|---|---:|---:|---:|
| nozzle、forebody、R2.54 shoulder、`x<=40 mm` 表面邻域 | L3 | 0.2083 | 1 |
| `40<x<=70 mm` 前圆柱 | L2 | 0.4167 | 1/8 |
| `70<x<=130 mm` 中圆柱 | L1 | 0.8333 | 1/64 |
| `x>130 mm` 后圆柱与人工底面主体 | L0 | 1.6667 | 1/512 |

补充流场区域：

- tri 的每一条离散外围轴周围 `d<=10 mm, -20<=x<=x_throat+5 mm` 保持 L3，并以 `d<=12 mm` 的 L2 parent envelope 包围；`d<=12 mm, -80<=x<-20 mm` 的喷流延伸区保持 L2。Run165 把同一规则移到中心 N1 轴。图中的 `x-r` 色带只是三条轴的方位投影，绝不能在三维代码中实现成连续环形/torus tag；实际必须逐喷管计算横向轴距 `d_j`。
- 弓形激波/剪切层用几何固定 L1 启动包络加密度梯度跟踪到 L2；先不要把整个激波体积永久锁为 L3。
- aft cap 内部留 L0；如圆柱末缘的锐边/尾迹敏感，可只把外圆约 8–10 mm 环带保留 L1。
- 若只做 bring-up，可先停在 L2，最细 `dx=0.4167 mm`；喉径约 15.2 cells、肩半径约 6.1 cells。局部 L3 候选使喉径约 30.4–30.9 cells、肩半径约 12.2 cells。

这里把“后体密度低”解释成 `dx` **更大**。`dx` 越小反而越密、越贵；L0 的 `dx` 是前体 L3 的 8 倍，同体积 cell 数是 1/512。人工 cap 不是 NASA 完整 base/sting；若目标是 forebody `Cp` 和喷流/激波相互作用，可把 cap 载荷排除并保持粗网格。若研究总阻力、底压、尾迹或圆柱壁摩擦，则不能把它视为无关。

### 19.3 当前代码不能仅靠 `user_tagging` 实现上述前密后疏

当前 [`CNS.cpp`](../../src/CNS.cpp) 的 `ib.amr_support_buffer=1` 在每一级查找**整个** IBM 流固界面并扩大 support band。因此 `max_level=3` 时，它会把长圆柱和 aft cap 也逐级推进 L3；问题不是用户 `user_tagging` 规则写得不够精细。

不要简单关闭 support buffer，也不能只修改其 interface seed。2026-08-03 的 tri/quad 严格门禁已经证明：对连续后圆柱表面使用 `x<=55 mm` 的一次性 seed 截断，即使再做标准 dilation，AMReX 矩形 boxing 仍会在截断外生成少量同级 GP；quad 实测 L1 有 `712` 个 target、`6408` 个非零 donor support 跨 coarse-fine 边界。把 `amr.n_error_buf` 从 `1` 提到 `2 2` 后缺失数完全不变，所以这不是普通 tag buffer 太小，而是局部截断没有形成“新增 fine GP -> 再补 fine support”的递归闭包。

因此统一的 body-frame `local_max_level(x,y,z)` 只能作为目标分布，**不能单独充当安全算法**。局部表面加密的 fine-grid 边界必须完全落在远离所有 GP/surface stencil 的纯流体或纯固体区、覆盖完整独立几何分量，或由新的迭代 support-closure 算法证明闭包；任何方案最终仍须保留：

```text
ib.amr_support_buffer = 1
ib.amr_support_audit  = 1
ib.amr_support_strict = 1
```

以当前 E2/I2、ratio=2 的典型设置，support tag 半径约为 7 个父级单元，即 L0/L1/L2 上约 `11.7/5.83/2.92 mm`；不能把它拍脑袋改成固定两格。level 切换放在直圆柱段并离喷口/肩部至少 10–20 mm，最终必须由 strict audit 给出零 missing support。

密度梯度 tag 只能在中心单元及其 `±1` 差分 stencil 全部是有效 fluid、且 stencil 不接触 IBM solid/ghost/jet imposed jump 时计算，并显式封顶到目标 L2。当前 [`CNS.cpp`](../../src/CNS.cpp) 已对 IBM 邻域中非网格独立的 gradient tag 导致 runaway regrid 给出警告；忽略这一条件会使后体重新被细化，甚至在粗细界面生成缺少 image-point support 的 child grid。

AMReX 是块结构 AMR。真实网格还会被 `n_error_buf`、proper nesting、Berger–Rigoutsos、`blocking_factor` 和 `max_grid_size` 扩成矩形 boxes，不会精确贴着设计色带。当前 tri/quad 的拓扑安全 bring-up baseline 已改为 `L0=(2.24,2.25,2.25) mm`、完整连接 IBM 界面 `L1=(1.12,1.125,1.125) mm`；初始化 regrid 后均为 `3,276,800 + 5,638,144 = 8,914,944` allocated cells，L0/L1 的 GP 与 surface missing 均为零并各完成一个 coarse step。后体外部体积仍为 L0，只有保持 donor 完整性所需的贴体薄带为 L1。这个 baseline 只有 10.85% 初始余量，且没有 0.5 mm nozzle level，只能用于集成/容量门禁，不能用于定量 SRP、载荷 PSD 或网格收敛结论。

### 19.4 局部 `dx` 还需要局部匹配的 STL

若把固定 STL 的生产建议 `hmax<=0.5 dx_local` 推广成经验证的局部门槛，上表对应建议硬上限约为 L3/L2/L1/L0 的 `0.104/0.208/0.417/0.833 mm`。当前 360k STL 不满足这一正式 L3 候选，只适合 bring-up/可视化。

按局部表面积作理想等边估算，这种前细后疏 STL 约为 5–6M faces，而不是全表面统一按本图 L3 `dx=0.2083 mm` 细化的约 26.5M；第 13.3 节的 28.75M 是 `dx=0.200 mm` benchmark。实际 Gmsh `lc` 不保证 hard max，最终面数会更高。约 5.7M 已接近 validation `MPI_Gatherv` 的 5.835M faces 硬上限，且 E2/I2 逐-rank IBM payload 粗估约 10–11 GB，所以在生成正式局部 STL 前仍应先完成第 17 节的 owned-only surface storage 和分布式输出。

可视化可重复生成：

```bash
.venv/bin/python geometry/test1853/visualize_test1853_mesh_distribution.py
```

### 19.5 job 862 的实际四层 hierarchy（不是设计图）

第 19.2 节的解析色带现已由 job 862 最终 plotfile 的真实 BoxArray 复核。新脚本
[`visualize_amr_hierarchy.py`](../../IBM/cases/test1853_run165_aws_gate/visualize_amr_hierarchy.py)
只读取 `Header` 和 `Level_0/Cell_H` 至 `Level_3/Cell_H`，不读取
`Cell_D_*` 流场二进制；它逐 box 核对 `Cell_H` 整数索引与 `Header` 物理边界，
并以“该位置实际覆盖的最细 level”生成切片：

- [`job862_plt00002_amr_hierarchy_slices.png`](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/amr_visualization/job862_plt00002_amr_hierarchy_slices.png)：实际 x-y、x-z 和 `x=10.25/25.25/100.25/260.25 mm` 的 y-z 截面；
- [`job862_plt00002_amr_hierarchy_forebody_zoom.png`](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/amr_visualization/job862_plt00002_amr_hierarchy_forebody_zoom.png)：同一 hierarchy 的喷管/前体/肩部放大；
- [`job862_plt00002_amr_hierarchy_stats.json`](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/amr_visualization/job862_plt00002_amr_hierarchy_stats.json)：metadata hash、每级 box/cell/composite-volume 和各截面面积统计。

| level | `dx`, mm | grids | allocated box cells | 去除细层覆盖后的 composite-active cells | composite 域体积占比 |
|---:|---:|---:|---:|---:|---:|
| L0 | 4.0 | 800 | 819,200 | 612,416 | 74.7578% |
| L1 | 2.0 | 878 | 1,654,272 | 902,976 | 13.7783% |
| L2 | 1.0 | 2,246 | 6,010,368 | 5,840,960 | 11.1407% |
| L3 | 0.5 | 805 | 1,355,264 | 1,355,264 | 0.3231% |

这里 `allocated box cells=9,839,104` 是计算/存储和 10M 门禁采用的数量；
`composite-active` 用于描述空间中最终由哪一级代表，不能拿后者替代容量预算。真实图还显示：
`x=10.25 mm` 的中心喷管/前体有 L3；到 `x=25.25 mm` 时外侧肩部仍是 L2，
只有中心喷管块为 L3；`x=100.25/260.25 mm` 的圆柱没有 L3，但完整 IBM
邻域仍保留大块 L2。因此当前门禁实现了“喷管 L3、全身至少 L2”，并没有实现
第 19.2 节所建议的肩部 L3 和后圆柱 L1/L0。仅余 `160,896 cells`（1.61%），
不能直接扩大肩部 L3；应先缩减后体 L2/固定流场包络并重新通过 strict support，
再用同一脚本和每次 regrid 的 10M guard 复核。

进一步按 `Cell_H` 的 FAB min/max 检查真实流固界面：L3 的 mixed-sld FAB 轴向并集只有 `x=0--16 mm`；延到约 `28 mm` 的其他 L3 boxes 是 `sld_min=sld_max=1` 的纯固体块，不能把它们算作肩部表面 L3。L2 mixed/ghost FAB 则连续到约 `x=280 mm`，所以人工 cap 也确实处于 L2。直接原因是当前 `ibm_support_refine_allowed()` 对 `parent_level<2` 无条件返回 true：全身同时 seed L0->L1 和 L1->L2，只限制了 L2->L3。

下面这个曾建议的 body-frame 上限现在只保留为**尚未实现的目标分布**，不得直接放入当前 one-pass interface seed：

```text
x <= 40 mm        local_max_level = 3
40 < x <= 70 mm   local_max_level = 2
70 < x <= 130 mm  local_max_level = 1
x > 130 mm        local_max_level = 0
```

`x=40 mm` 虽位于肩后直圆柱段，但连续表面并不会因为离肩部较远就允许硬切 level。原先按 job862 boxes 外推的 `8.1--8.8M` 只是 **[DERIVED] 容量估算**，已被上述 `6408` 个缺失 support 的实测门禁否决，不能作为可运行配置。以后若要恢复 0.5 mm nozzle/前体局部层，必须先实现或证明递归 support closure，并用 `max_step=1` 的初始化、真实 regrid 及后续动态 regrid 反复要求所有 GP/surface missing 为零、总 allocated cells 始终小于 10M。

## 20. Run 165 上游启动间断面与 AWS CPU 实测

### 20.1 已实现的流场启动方式

用户指出旧设置错误后，`IBM/cases/test1853_run165_aws_gate/prob.h` 已改为仓库 sphere/SRP 案例采用的 impulsive-start 模式。当前门禁只取轴向 `alpha=0`，不是 Run 165 的整段攻角扫描：

```text
x <= -0.052 m：Run 165 M=4.6，沿+x
x >  -0.052 m：同一p/T/rho，u=0
x-low：持续给定Run 165全状态超声速入口
x-high、y-low/high、z-low/high：一阶外推
```

`x=-0.052 m` 只依赖 `x`，覆盖整个 `y-z` 截面；它是当前 `xlo=-0.160 m, dx_L0=4 mm` 的精确 cell face，距 SI STL 的 `xmin=2.44126888 mm` 为 `54.441 mm`。这个位置是为启动距离和网格对齐选取的 **[ASSUMPTION] 数值参数**，不是 NASA 几何尺寸、风洞硬件位置或测得激波位置。

**[DERIVED]/[ASSUMPTION]** 唯一覆盖是实体中真实 N1 发散内通道 `2.441<=x<=14.415 mm` 内的 Run-165 quasi-1D 等熵启动状态。它避免从风洞静压直接跳到约 `2.18 MPa` 喉静压所导致的旧 HLLC 非有限状态；其作用域不包括 forebody 外部或“geom 周围”。如将来要求先建立无喷流外场、再开启喷流，必须给 throat wall 增加所有 MPI ranks 同步的物理时间和有限时间 ramp，不能从 checkpoint 瞬时开启。

等压、等密度、速度不同的两侧不是单一 Rankine--Hugoniot 激波。**[DERIVED]** 按一维理想气体 Riemann 问题估算，快/慢波约为 `497.5/244.5 m/s`：第一波约 `109 us` 到达 STL，慢波约 `223 us` 到达，约 `1.65 ms` 扫过 `x-high`。**[ASSUMPTION]** 因此先把 `2 ms` 作为第一轮全域 flushing 目标，再根据弓形激波/喷流激波位置、表面 `Cp`、Cartesian throat 质量流量和积分载荷判断是否进入可统计区间；`2 ms` 不是 NASA 测量的稳态时间。

### 20.2 AWS job 862 可复核结果

作业 `862` 使用 8 个 `hpc8a.96xlarge` CPU 节点、每节点 96 MPI ranks、合计 768 ranks；调度器为每 rank 保留 2 CPUs。计算节点用 GCC 11.5/OpenMPI 4.1.7 从实际 campaign source 重编译。保存的 shell 日志从 campaign 时间戳到 completion marker 为 `61 s`，覆盖编译、node-local STL 分发、初始化、两步、两次 plot 和后处理，但不覆盖节点 provisioning、首时间戳前和 marker 后工作。随后外部实时 `scontrol` 查询返回 `COMPLETED, ExitCode=0:0`；该 scheduler 快照不是原始 Slurm accounting 文件，因为集群关闭了 accounting storage。文件锁定值为：

| 对象 | SHA-256 |
|---|---|
| AWS executable | `aea85eb98be98660c6380b01e0d9bb1a4cc406e6bc9d1109d28f54eb7d70bbe1` |
| `prob.h` | `f3c041565590e7660be63a11404a0034a8d2a4aba4cdac3e0df3130fe086afe0` |
| `inputs_gate` | `bab1bd09cb62db91ed9c89a894027222a49fc05b4b3794a9154756c7d878cf0c` |
| Run165 SI STL | `f8d8095fa4d093fe421c1c6cc71712a1709839ddfeefd116cb8b91d91dd9df82` |

这些外部查询值、提交时脚本/analyzer hash、证据文件 hash 和事后复算命令另存于 [`job862_manifest.json`](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/job862_manifest.json)。AWS executable 没有下载进本地证据包，因此其 hash 是远端查询值；campaign source 也不是 Git checkout，日志中的 `source_sha256` 只标识 `src/CNS.cpp`，不能冒充整个 source snapshot hash。

实测层级是：

| level | `dx`, mm | BoxArray grids | box cells |
|---:|---:|---:|---:|
| L0 | 4.0 | 800 | 819,200 |
| L1 | 2.0 | 878 | 1,654,272 |
| L2 | 1.0 | 2,246 | 6,010,368 |
| L3 | 0.5 | 805 | 1,355,264 |
| **总计** | — | **4,729** | **9,839,104** |

以下两图直接解析 job862 最终两步启动 plotfile（`t=8.856607227e-8 s`）的 BoxArray；红线是同一 SHA 的 Run165 STL 截面。没有读取流场 `Cell_D_*`，也没有用 `user_tagging` 解析区域反推网格。启动前沿/激波继续发展后 tags 会变化，所以这不是稳态或 production hierarchy：

![job862 实际 AMReX hierarchy 切片](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/amr_visualization/job862_plt00002_amr_hierarchy_slices.png)

![job862 喷管、前体和肩部实际 AMR 放大](../../IBM/cases/test1853_run165_aws_gate/aws_results/job862/amr_visualization/job862_plt00002_amr_hierarchy_forebody_zoom.png)

图中 L3 确实集中在 N1 passage/近前体 guard；`x=10.25 mm` 有 L3，但 `x=25.25 mm` 外侧肩圆周仍为 L2，只有中心喷管 proper-nesting 小块延续 L3。`x=100.25 mm` 和 `260.25 mm` 的圆柱/近底面没有 L3，却仍有宽 L2 区域。L2 占全部 allocated cells 的约 `61.09%`，而 L3 只覆盖 composite 物理域体积的 `0.323%`。所以当前 hierarchy 是“低于10M的功能门禁”，还不是用户最终要求的“喷管/前体/肩最密、后圆柱/底面明显更粗”分配。不能在只剩 `160,896 cells` 时直接扩展肩部 L3；应先让 aft cylinder/cap 的 L1->L2 IBM baseline seed 服从安全的 axial local-max-level 并用 strict support 回归，再把释放的预算用于肩部/前体。

最终层级比 10M 上限只少 `160,896 cells`，余量 `1.61%`。当前脚本只从最终 plotfile 证明这个值；正式长算必须在每次 regrid 后于求解器内做 global BoxArray cell guard，否则中间层级可能越过上限而最终层级又降回来。

所有 L0--L3 GP 与 surface strict support 审计均为 `coarse-fine-missing-targets=0`、`coarse-fine-missing-supports=0`；日志没有 NaN、non-finite、positivity、signal 或 OOM 失败。两次 coarse step 分别为 `5.757810575 s` 与 `5.652851597 s`，平均 `5.705331086 s`；对应最终物理时间仅 `8.856607227e-8 s`。这证明新初值可以启动和推进，不证明流场已到模型，更不证明实验一致性。

IBM 计时要分成冷启动、regrid 重建和每步重复三类：

| 指标 | job 862 实测 | 判断 |
|---|---:|---|
| 完整 `read_geom` | 0.564004 s | 正常，不是当前瓶颈 |
| STL mesh I/O | 0.270194 s | 正常 |
| closed check | 0.154552 s | 正常 |
| BVH build | 0.157063 s | 正常 |
| in/out marker，68 calls critical-path sum | 5.415515 s | regrid 压力测试热点 |
| `initialiseGPs`，44 calls critical-path sum | 3.886411 s | regrid 压力测试热点 |
| image-point fused geometry，44 calls sum | 0.169604 s | 当前规模正常 |
| complete IBM rebuild，44 calls 的 MPI-rank maxima 之和 | 9.565132 s | **29.33% `Run Time total`，离线复算报警** |
| recurring `IBM::computeAllGPs` inclusive | 2.129 s，6.52% | 可接受但应继续跟踪 |

作业故意令 `amr.regrid_int=1`，所以 AMR subcycling 中每个 eligible level/substep 都可能 regrid；日志实测 `Amr::regrid=10 calls`、`rebuild_ibm_total=44 calls`。这适合暴露 GP/in-out/image-point 问题，却会高估正常 production schedule 的平均步成本。原 job862 使用旧分析器，原始 JSON/日志是 `IBM_TIMING_STATUS=OK` 和 `RUN165_CPU_GATE_PASS`；更新分析器在本地对同一日志离线复算才得到累计 `29.33%` WARN。未来重新运行同步后的脚本遇到同类警告时才会输出 `RUN165_CPU_GATE_FUNCTIONAL_PASS_PERFORMANCE_WARN`。

AMReX 只打印了 FArrayBox payload high-water：`7558--19307 KiB/rank ≈ 7.381--18.854 MiB/rank`；它不含 IBM geometry、进程 runtime 和通信内存。几何基础 payload 仍按第 17.4 节每 rank 全复制。此次没有 OOM，说明约 353k facets / 9.84M cells 在 8 节点 CPU 配置上可用于功能门禁；但该集群关闭 Slurm accounting storage，`MaxRSS/AveRSS` 不可得，所以不能声称已测得峰值 RSS。下次正式长算应由节点侧 sampler 记录每 rank 和每 node 的 RSS/PSS、cgroup memory.current/peak，并把缺失样本作为门禁失败。

### 20.3 30 分钟与正式启动策略

当前 `CFL=0.01` 的 coarse `dt` 约 `4.4283e-8 s`；两步门禁只走了 `0.0886 us`。可先以 `CFL=0.01--0.02` 跑 20--50 步检查喉口/IBM 正性，再逐级测试 `0.05`。当前运行路径尚无经验证的 checked-star/HLLE fallback；在补齐并验证之前不应冒进到 `CFL=0.1`。若只作乐观墙钟时间估算，按 `CFL=0.1` 和 job862 的压力测试步时，第一启动波需严格向上取整约 248 coarse steps、纯步进约 `23.6 min`；完整左侧来流约 503 步、约 `48 min`；`2 ms` 约 4,517 步、约 `7.2 h`。首段 CFL warm-up、构建、初始化、regrid 和输出都另有开销，因此不能承诺 30 分钟内波前到达；降低 production regrid 频率也必须重新实测，不能线性承诺。

启动前沿会移动。若 `CFL=0.1` 且每 20 coarse steps 才 regrid，前沿约移动 `4.4 mm`，大于当前 `n_error_buf=1` 对 L2 tag 的保护尺度；可先用 `CFL≈0.02, regrid_int≈20`，或在 `CFL=0.1` 时把前沿相关 regrid 缩到约 4 步。提高 `n_error_buf` 会进一步挤压仅 1.61% 的 cell 余量。x/y/z 外推边界满足“只有 x-low 注入”的启动要求；长时间计算若弓形激波接近侧边界，还必须做计算域/characteristic side boundary 敏感性测试。

job862 的原始日志、最终层级计数、IBM timing JSON 和最小 plot 元数据保存在 `IBM/cases/test1853_run165_aws_gate/aws_results/job862/`。此前全域直接初始化来流的 AWS 作业只能作为旧代码性能参考，不能用于这套修正后启动物理的结论。
