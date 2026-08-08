# NASA Langley UPWT Test 1853 几何重建

> IBM 交接入口：[`IBM_HANDOFF_REPORT_CN.md`](IBM_HANDOFF_REPORT_CN.md)。该报告汇总全部公开尺寸、重建站位、STL 面片/曲率/误差、SI 变换、throat 边界使用方法、正式网格门槛，以及当前 Cerisse `read_geom`/BVH 的兼容性、逐-rank 内存和通信改造优先级。

本目录生成三类实验硬件构型：

- **Run 165**：中心 Nozzle 1；外围三个位置为齐平 plug。
- **Runs 247/262/263**：中心齐平 plug；外围 Nozzles 2、3、4。它们使用同一 tri 实体硬件，flow/roll 条件按 run 区分。
- **Runs 307/315 quad**：中心 Nozzle 1 与外围 Nozzles 2、3、4 全部启用。输出仍以 `phi=0` 几何 master 为基准，准备器按所选 run 应用滚转。

此外还提供一个明确标注为非实验 Run 的参数化构型：

- **four-peripheral ring**：无中心喷口，四个相同名义外围喷口位于
  `63.5 mm` 节圆直径上，方位角为 `0°/90°/180°/270°`。它用于四重对称
  3-D 数值研究，不能标记为 Test 1853 Run 307/315。
- **single surface-normal**：一枚名义外围喷口位于 `theta=0°`、
  `R=31.75 mm`。virtual-exit 平面与局部 70° 前体相切，喷管轴沿局部表面
  法向，射流方向为 `(-sin70°, 0, +cos70°)`。这是参数化研究构型，不是
  NASA Test 1853 的实验硬件；“相切”指出口平面相切，不指喷管轴沿表面切线。

几何依据 NASA 官方 [NASA/TP-2014-218256](https://ntrs.nasa.gov/citations/20140006403)、报告中的 as-built Table A-2 以及 DWG 1168288–1168298。工程图原页的无损摘图与来源索引见 [`source_extract/SOURCE.md`](source_extract/SOURCE.md)。没有使用仓库中既有的简化 `trinozzle.stl` 或二手 CFD 图形反推尺寸。

## 交付内容

运行生成器后，`outputs/` 包含：

```text
cad/
  test1853_run165_single.{brep,step}
  test1853_run262_263_tri_n3at120.{brep,step}
  test1853_run262_263_tri_n3at240.{brep,step}
  test1853_quad_n3at120.{brep,step}
  test1853_quad_n3at240.{brep,step}
stl/
  对应的五个闭合、水密实体 STL
patches/<model>/
  body_wall_excluding_throat_caps.stl
  nozzle_N_throat_cap.stl
reference/<model>/
  nozzle_N_virtual_exit_REFERENCE_ONLY.stl
  physical_lip_reconstruction.csv
reference/run_roll_transforms.json
figures/
  主体轴向剖面、喷管前视布置、外围斜切唇线变化
  当前逐面 h_STL 热图、局部真三角线框、轴向密度曲线、建议 IBM dx 分区图
visualization/
  带逐面 hmax/hmean/area/x/zone 标量的可交互 VTP
reports/
  build_audit.json
  geometry_audit.md
  independent_validation.{json,md}
  mesh_quality_detailed.json
  mesh_distribution_visualization.json
  ibm_code_compatibility_audit.json
ibm_prepared/
  五个官方构型及参数化构型的 SI-metre master、patch、metadata 和 IBM audit
four_peripheral_ring90_build/
  四外围参数化构型的独立 CAD/STL/reference/reports 生成树
single_peripheral_surface_normal_build/
  局部法向单喷管的独立 CAD/STL/reference/reports/figures 生成树
```

三喷管和四喷管各输出两种编号映射，是因为 NASA 只明确 Nozzle 2 在 `θ=0°`；公开报告没有说明 Nozzle 3/4 与 `120°/240°` 的一一关系。两枚喷管的 as-built 尺寸略有差别，所以不能静默指定其中一种。四喷管不是四个外围喷管，而是中心 N1 加三个外围喷管 N2/N3/N4。

参数化 `test1853_four_peripheral_ring90` 与上述官方 quad 严格分开：四个
喷口均使用 `Dt=6.35 mm`、virtual `De=12.7 mm`、扩张半角 `15°` 的相同
名义外围轮廓，中心位置保持齐平连续 OML。

为保持已锁定的 tri 文件路径与 SHA，tri 文件名仍保留 `run262_263`，但其 `build_audit.json` 的 `run_ids` 已扩展为 `[247,262,263]`；三者使用相同实体硬件。quad 的 `run_ids` 为 `[307,315]`。

## 几何定义

所有 CAD 坐标使用 mm。NASA 图纸的 TSP（70° 理论尖锥点）为 `x=0`，`+x` 指向下游，射流向 `−x` 排出：

```text
r = sqrt(y² + z²)
y = r sin(theta)
z = r cos(theta)
theta=0 位于模型顶部 +z
NASA DWG 1168296 的下游视向前视图中，theta=90（+y）位于画面左侧
```

主体采用公开数据可支持的名义相切重建：

- 最大直径 `5.000 in = 127.000 mm`；
- 70° sphere-cone；名义鼻球半径 `1.000 in`；
- 图纸明确肩圆角 `R0.100 in = 2.540 mm`；
- 圆柱半径 `2.500 in = 63.500 mm`；
- aft-cover 起点 `x=1.0066 in`，长度 `9.550 in`，末端 `x=10.5566 in`；
- 后体 cover 壁厚 `0.075 in`，但外部 CFD 实体按实心封闭体输出；
- aft end 以平面封闭，供 STL/IBM 使用。

解析外形（英寸）为：

```text
sphere:   r = sqrt(1 - (x - 1.064177772)^2)
          0.064177772 <= x <= 0.124485152
cone:     r = x tan(70 deg)
          0.124485152 <= x <= 0.885977077
shoulder: r = 2.4 + sqrt(0.1^2 - (x - 0.979946339)^2)
          0.885977077 <= x <= 0.979946339
cylinder: r = 2.5
          0.979946339 <= x <= 10.5566
```

这里没有把图纸 `x=0.877 in` 错当成锥—肩切点；它是 B–B 剖切/结构站位。`Ø4.851 in` 是前后体配合直径，`x=1.0066 in` 是 aft-cover 前缘，也都不是 OML 切点。

## as-built 喷管

| Nozzle | 喉径 Dt, in | virtual De, in | At, in² | virtual Ae, in² | Ae/At | 扩张壁半角 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.2491 | 0.5014 | 0.048735 | 0.197451 | 4.051550 | 14.981° |
| 2 | 0.2498 | 0.4965 | 0.049009 | 0.193610 | 3.950514 | 14.977° |
| 3 | 0.2535 | 0.4989 | 0.050471 | 0.195487 | 3.873210 | 14.994° |
| 4 | 0.2503 | 0.4996 | 0.049205 | 0.196036 | 3.984035 | 14.985° |

外围三轴位于 `r=1.250 in`、相隔 120°，并全部平行于模型 x 轴。构造时从实体减去 as-built 扩张锥，真实外露唇线由喷管锥壁与 70° OML 的 Boolean 交线产生。因此外围物理开口不是 `virtual De` 的平面圆，也不是局部法向倾斜的喷管。

每条流道在 as-built 喉面结束，并输出独立平面 `throat_cap`，可直接作为 CFD 入口 patch。完整收敛段、plenum 和 compound throat 圆角不进入主模型，因为公开表只给直径/面积/出口角，无法唯一恢复其 as-built 轴向轮廓。

## 明确的重建假设

NASA 图纸多处写有 “THIS SURFACE DEFINED BY CAD FILE”，但原始 CAD 与最终装配 DWG 1168299 没有公开。因此以下内容必须区分：

1. **图纸/as-built 明确值**：主体直径、70°、R0.100、aftbody 尺寸、喷管 PCD、每枚喷管 Dt/De/At/Ae/角度等。
2. **唯一解析推导**：名义 1-in sphere 与 70° 锥相切、70° 锥/R0.100/5-in 圆柱相切位置。
3. **显式 CFD 替代**：中心 virtual-exit 圆与重建球面齐平；外围 virtual-exit 平面通过喷管轴与理想 70° 锥面的交点。
4. **保持未知**：四个真实 virtual-exit 全局 x、外围精确插入深度、Nozzle 3/4 方位映射、CAD-defined 插件表面、完整 manifold、sting 全局装配站位和外部线束/胶带。

所以这些文件是“公开资料所能支持的最高可审计外部 CFD 重建”，不能称作与未公开 NASA STEP 逐点相同。`geometry_parameters.yaml` 对未知量保留 `null`，不会用猜测数值填空。

## 参数和实验测点

[`geometry_parameters.yaml`](geometry_parameters.yaml) 是唯一几何参数源，采用 JSON 语法（同时是合法 YAML 1.2），不依赖 PyYAML。它还保留 plenum、interface plate、insert、sting 等未进入主外形的图纸尺寸及来源状态。

[`surface_ports_table_a1.csv`](surface_ports_table_a1.csv) 包含 Table A-1 全部 176 个测点的 `r, θ, x`，并给出由报告中三位小数 `r/θ` 重新计算的 `y/z`。这些点供 CFD probe/实验对比使用；默认不在 OML 上切 0.020-in 微孔，以免给 IBM 几何引入无意义小特征。P126、P127 标为装配损坏。

## 生成

环境中没有 SALOME/FreeCAD，但已有 Gmsh 4.15.2 的 OpenCASCADE 内核。批处理命令为：

```bash
cd /home/qiaoj/testcerisse/cerisse
.venv/bin/python geometry/test1853/build_test1853.py
MPLCONFIGDIR=/tmp/matplotlib-test1853 \
  .venv/bin/python geometry/test1853/render_test1853.py
.venv/bin/python geometry/test1853/validate_test1853.py
.venv/bin/python geometry/test1853/analyze_test1853_mesh_quality.py
.venv/bin/python geometry/test1853/visualize_test1853_mesh_distribution.py
```

开发时可加 `--quick` 生成较粗网格；正式交付不要使用。STEP/BREP 是 CAD 母版，STL 只是受控三角化结果。

四外围参数化构型及四构型 3-D 对比图的命令为：

```bash
.venv/bin/python geometry/test1853/build_test1853.py \
  --output geometry/test1853/outputs/four_peripheral_ring90_build \
  --config test1853_four_peripheral_ring90
.venv/bin/python geometry/test1853/validate_test1853.py \
  --output geometry/test1853/outputs/four_peripheral_ring90_build \
  --config test1853_four_peripheral_ring90
.venv/bin/python geometry/test1853/render_four_configurations.py
```

局部法向单喷管及其三正交视图加斜视图的命令为：

```bash
.venv/bin/python geometry/test1853/build_test1853.py \
  --output geometry/test1853/outputs/single_peripheral_surface_normal_build \
  --config test1853_single_peripheral_surface_normal
.venv/bin/python geometry/test1853/validate_test1853.py \
  --output geometry/test1853/outputs/single_peripheral_surface_normal_build \
  --config test1853_single_peripheral_surface_normal
.venv/bin/python geometry/test1853/render_surface_normal_configuration.py
```

正式生成固定 `reproducible_num_threads=1`，避免 Gmsh 并行表面网格在几何不变时产生 byte-different STL。当前锁定 hash 已通过连续两次全量生成复现。

## 验证门

生成和独立验证会检查：

- Boolean 后、BREP 重导入和 STEP 重导入均为一个正体积 solid；
- 包围盒、外径、总长和每枚喉面 CAD 面积；
- 每个配置只有正确的 active nozzle throat patches；
- STL 水密、manifold、统一外法向、单连通实体、正体积；
- 零面积面、重复三角面；
- 三角形边中点/重心到原 OCC 曲面的最大采样弦差，正式阈值 `0.02 mm`；
- GTS 的可定向 manifold 与自相交检查；
- 文件 SHA-256，便于计算任务锁定同一版几何。

Run 165 single 和 quad 的最上游实体点都是中心喷管圆唇，而不是被替换掉的球鼻尖；三喷管/中心 plug 构型则保留球鼻尖。这一差异已进入各自包围盒检查。

## IBM 的 SI 副本

当前 IBM 直接使用 STL 原始坐标，不从输入文件读取单位或 scale。米制算例应使用预变换副本，例如：

```bash
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration tri_n3at120 --run 262 \
  --tsp-m 0.275 0.275 0.275

# Run 307 quad；省略 --run 则输出 phi=0 几何 master
.venv/bin/python geometry/test1853/prepare_test1853_ibm.py \
  --configuration quad_n3at120 --run 307
```

若使用重新三角化的高分辨率 mm/TSP-frame STL，必须再传 `--input-stl /absolute/path/to/highres.stl` 并核对 metadata 的 source hash；否则工具默认读取 `outputs/stl/` 下当前约 36 万面的 master。匹配的开放辅助 patches 可用 `--patch-dir` 显式提供，未提供时不会混入旧 patches。

metadata 同时记录参数文件/准备器 SHA；指定 `--run` 时还会嵌入该 run 的完整核定条件。Run247 只接受 tri，Runs307/315 只接受 quad，错误的 run/configuration 组合会在写文件前失败关闭。

四外围构型使用 `--configuration ring4`，且禁止绑定官方 run ID。当前 SI
master 位于
`outputs/ibm_prepared/test1853_four_peripheral_ring90_master_SI/`。

局部法向单喷管使用 `--configuration single_surface_normal`，同样禁止绑定官方
run ID。当前 SI master 位于
`outputs/ibm_prepared/test1853_single_peripheral_surface_normal_master_SI/`。

本 Test 1853 构型只把生成目录中的一个完整闭合 STL 交给 `ib.filename`。`body_wall_excluding_throat_caps.stl` 和各 throat-cap STL 是开放的离线识别文件，不能作为独立 IBM 几何加载。正式亚毫米 IBM 计算应按仓库 fixed-STL 文档采用的推荐 production audit gate `max edge / finest dx <= 0.5` 重新三角化；这不是 solver 硬限制或充分精度证明。当前约 35–36 万面 STL 是接入/可视化版本。

tri 与 quad SI master 都已在各自独立 3-D case 中通过现有 CPU IBM reader、闭合校验、BVH2/BVH4、marker、严格同层 GP/surface 支撑审计和一个极短时间步。这个结果只认证当前约 36 万面 STL 的集成/容量门禁，不是 SRP 物理认证。对应 case 位于 `IBM/cases/test1853_tri_3d` 和 `IBM/cases/test1853_quad_3d`；报告第 17 节同时说明，18–29M 面的假设性超细 STL 在现有“每 rank 全复制 + 全面片 surface SoA + rank-0 gather”路径上仍不可直接生产运行，必须先完成所列 P0 改造。
