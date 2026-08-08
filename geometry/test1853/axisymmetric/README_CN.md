# Run 165 二维轴对称切面预览

这是 NASA Test 1853 Run 165 中心单喷管构型的二维轴对称交付。除原始中心
N1 切面外，现在还包含一套严格同坐标的几何 A/B：解析 N1 喷管，以及从名义
R25.4 鼻尖连续封闭的无喷管数值对照。两者已有独立 Cerisse case；当前仍只称为
几何和短时集成门禁，不宣称是 production 或统计稳定的 SRP 流场。

## 文件

- `outputs/run165_axisymmetric_meridional_section.png`：完整镜像剖面与正半径局部放大。
- `outputs/run165_axisymmetric_meridional_section.svg`：同一张矢量图。
- `outputs/run165_axisymmetric_section_preview.dat`：Cerisse 2-D IBM 预览多边形，每行严格为 `(r_m, x_TSP_m)`，最后一点由 reader 隐式连接回第一点。
- `outputs/run165_axisymmetric_section_profile.csv`：带 segment/source 标签的同一组点。
- `outputs/run165_axisymmetric_section_metadata.json`：参数、SHA-256、闭合/方向、OCC 体积和实际 STL 切面交叉校验。
- `build_run165_axisymmetric_section.py`：无 Gmsh/Trimesh 依赖的可复现生成器，仅需 NumPy 与 Matplotlib。

### 有/无喷管 A/B

- `build_run165_nozzle_geometry_ab.py`：同时生成、审计两套几何并绘制同尺度侧视图。
- `outputs/nozzle_geometry_ab/run165_rz_with_nozzle_profile.dat`：解析 N1 发散段；与审计 Case T 字节相同。
- `outputs/nozzle_geometry_ab/run165_rz_without_nozzle_profile.dat`：封闭球鼻，无开口、通道或喷流 patch。
- `outputs/nozzle_geometry_ab/run165_rz_nozzle_geometry_side_view_comparison.png`：完整侧视图与前体局部对比。
- `outputs/nozzle_geometry_ab/run165_rz_nozzle_geometry_ab_metadata.json`：语义、几何 QA、旋转体体积与 SHA-256。

对应算例目录为：

```text
IBM/cases/test1853_run165_rz_with_nozzle
IBM/cases/test1853_run165_rz_without_nozzle
```

这里的 `without_nozzle` 是明确标记的数值对照，不是 NASA “无喷管 Run165”硬件；
它也不是仍在 virtual exit 施加喷流的 Case E。

## 几何边界顺序

多边形按 CCW 顺序写出：

```text
(r=0, x=aft)
 -> axis closure
 -> N1 throat face（后续接入 sonic jet BC）
 -> N1 divergent wall
 -> nominal R25.4 sphere
 -> 70-deg cone
 -> NASA R2.54 shoulder
 -> R63.5 cylinder
 -> implicit artificial aft cap back to the first point
```

这与当前 3-D Run165 OCC/STL 实体的拓扑一致。解析旋转体积与已审计 OCC 体积约 `3,200,265.447 mm^3` 一致；生成器还独立读取 SHA 为 `f8d8095f...` 的 353,042 面 SI STL，检查 `y=0` 截线为单一闭环，并把实际截线叠加到图中。STL 只用于 QA，不作为二维曲线定义来源。

## 来源边界

- N1 `Dt=6.32714 mm`、virtual `De=12.73556 mm`、发散半角 `14.981 deg`：NASA Table A-2 as-built。
- 主体直径 `127 mm`、70°、肩部 `R2.54 mm` 和后体长度：NASA 报告/图纸。
- 鼻球 `R25.4 mm`：公开测点强支持的名义重建，不是公开图中明确给出的 as-built 点云。
- N1 virtual exit 与名义球面齐平及其全局 `x`：公开资料一致的重建假设。
- 平 aft cap 和 `r=0` 轴线：为了单个水密 IBM polygon 添加的数值闭合，不是完整 NASA base/sting。
- 没有擅自加入公开图无法唯一定位的收敛段、复合圆角、plenum、sting 或测压孔。

轴对称化只适用于 Run165 中心喷管、攻角 `alpha=0 deg`。Run165 的非零攻角点和外围三喷管构型不能用这个 RZ 模型。

## 后续接入约束

正式 Cerisse case 的编译配置必须采用：

```text
DIM = 2
USE_GPIBM = TRUE
```

运行输入中的关键项至少包括：

```text
geometry.coord_sys = 1
geometry.prob_lo[0] = 0
geometry.is_periodic[0] = 0
cns.lo_bc[0] = 3
ib.filename = run165_axisymmetric_section_preview.dat
```

在 RZ 中第一坐标/`UMX` 是径向，第二坐标/`UMY` 才是轴向；上游启动间断面必须写在第二坐标上，throat 的局部正法向速度对应全局 `-x` 喷流。

当前内置 2-D surface measure 是线长而不是环面积。总推力/阻力必须在后处理中采用
`dA=2*pi*r_mid*ds`；人工 aft cap 和 axis closure 必须排除。两个新 case 已固定
`prob.h` 合同、上游整面启动间断和几何文件哈希，但 RZ/AMR 的守恒共享面正性限制器
以及 mixed slip/blowing 的冻结 BI-CWLS closure 尚未完成，因此不能把短门禁改写成
production 几何或物理验证。
