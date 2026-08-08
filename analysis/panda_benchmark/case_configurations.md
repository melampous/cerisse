# U_jet (Panda & Seasholtz 1999) 案配置总档 — m119 / m142 / m180 (2D) + m142 (3D)

来源:各案目录的 `inputs` 与 `prob.h` 原文件逐项核对(2026-07-17),非记忆复述;评审修订 2026-07-17(统计分量数、探针计数、网格履历、模型定位)。
本档为 **as-built 配置档**(可复现算例、支撑平均激波场研究),不是结果档。
三个 2D 案的 `inputs` **逐字节完全相同**;`prob.h` 仅差喷压比一行(m180 另有一处 tagging 细节)。

**模型定位(正式名称)**:*Panda-condition, prescribed-exit, simplified flush-baffle model*。
不称 Panda-exact,原因(与实验的已知差异,均"未建模/非实验值"):
| 项 | 实验 | 本模型 |
|---|---|---|
| 喷管内部 | 收敛喷管内流动 | 未解析,出口面直接规定 p, T, u |
| 唇口 | 实体唇厚 | 无实体唇,tanh 速度剖面 δ/D=0.01(非实验测量值) |
| 流量系数 | 未报 | C_d,eff=0.8789 为剖面推论,非实验值 |
| 法兰 | 305 mm 金属喷口块(Panda 指出是稳定大幅 screech 的主要反射体) | 理想共面平板(slip wall) |
| 环境流 | ~8 m/s 横向装置流(密度实验) | 静止(首基线可忽略,须注明未建模) |

---

## 1. 共同物理设定(四案全同)

### 1.1 气体与输运
| 项 | 值 |
|---|---|
| 气体模型 | calorically perfect gas, γ = 1.4, M_w = 28.96 g/mol → R = 287.06 J/(kg·K) |
| 黏性 | Sutherland: μ(T) = 1.458e-6 · T^1.5/(T+110.4) Pa·s (`visc_suth_t`) |
| 导热 | Sutherland 型 (`cond_suth_t`);常数参考值 k=0.0262 W/(m·K), μ=1.85e-5 Pa·s (在 `methodparm_t` 中,供非 Sutherland 路径) |
| 方程 | 可压缩 Navier–Stokes(viscous_t,黏性项 2 阶中心),无源项 |

### 1.2 环境与喷口(Panda 实验日条件)
| 项 | 值 |
|---|---|
| 环境压力 p_amb | 99 780 Pa |
| 环境温度 T_amb = T0(不加热) | 297.15 K → ρ_amb = 1.1698 kg/m³ |
| 环境流动 | 静止(无 coflow) |
| 喷口 | 声速收敛喷嘴(choked, M_e = 1),D_e = 25.4 mm(r_jet = 12.7 mm) |
| 出口静温 T_e = T0/1.2 | 247.625 K → c_e = U_e = 315.5 m/s |
| 法兰(flange) | 半径 6D = 152.4 mm(Panda 305 mm 平板) |

### 1.3 入口速度剖面(Dirichlet,三区 x-low/z-low 边界的孔内区)
- tanh lip 剖面:f(r) = ½[1 − tanh((r − r0)/δ)],δ = 254 µm(δ/D = 0.01),**中心 r0 = R − 3δ**,使 u(R) ≈ 0(声学孔径 = 全 D)。
- **静压在整孔均匀 = p_e**(全压力孔径;2026-07-15 评审修订:此前把 ρ/T/p 一起 blend 会把有效压力孔径缩到 ~0.94D、削弱第一道 expansion fan)。
- 温度按绝热边界层关系 T = T0 − u²/(2c_p)(总温恒定),ρ = p_e/(RT)。
- 有效流量系数(banner 积分):**C_d,eff = 0.8789**(相对理想 top-hat)。
- 法兰区(R ≤ r ≤ 6D):adiabatic slip wall,exact mirror ghosts(法向动量奇、其余偶 → u_n = 0、零剪切、dT/dn = 0)。
- r > 6D:静止环境 Dirichlet。

### 1.4 各案喷压比与完全膨胀参数
p_jet(出口静压) = 99 780 × NPR/1.89293;Mj、T_j = T0/(1+0.2Mj²)、U_j = Mj·√(γRT_j)、ρ_j = p_amb/(RT_j);Re_j 按 Panda 惯例(完全膨胀态)。

| 案 | NPR (p0/p_amb) | p_e/p_amb | p_e [Pa] | ρ_e [kg/m³] | Mj | T_j [K] | U_j [m/s] | ρ_j [kg/m³] | Re_j |
|---|---|---|---|---|---|---|---|---|---|
| m119 | 2.393613 | 1.2645 | 126 173 | 1.775 | 1.19 | 231.6 | 363.0 | 1.5011 | 0.92×10⁶ |
| m142 | 3.273446 | 1.7293 | 172 550 | 2.428 | 1.42 | 211.75 | 414.2 | 1.6413 | 1.24×10⁶ |
| m180 | 5.745796 | 3.0354 | 302 871 | 4.261 | 1.80 | 180.3 | 484.5 | 1.9276 | 1.95×10⁶ |

(Panda 论文报 0.93×10⁶ / 1.26×10⁶,与 m119/m142 的**名义完全膨胀 Reynolds 数**一致;实际入口质量/动量通量受 C_d,eff = 0.8789 影响,不等同于实验出口边界层/实际流量。)

---

## 2. 共同数值方法(四案全同)

| 项 | 值 |
|---|---|
| 无黏通量 | WENO-Z5,characteristic-wise flux-vector splitting + face-local Lax–Friedrichs(`weno_t<ReconScheme::WenoZ5>`) |
| 黏性通量 | 2 阶中心(`viscous_t`, order=2, use_LES=false) |
| 时间推进 | SSP-RK3(`cns.order_rk=3, stages_rk=3`),CFL = 0.3 |
| AMR 框架 | AMReX,ref_ratio = 2(逐级),`amr.grid_eff = 0.75`,`cns.do_reflux = 1` |
| 在线统计 | dt 加权一/二阶矩;**2D:11 个分量**(u,v 的 MEAN+SQR、uv、p/T/ρ 的 MEAN+SQR);**3D:15 个分量**(u,v,w 的 MEAN+SQR、uv/uw/vw、p/T/ρ 的 MEAN+SQR);重启时自动清零,`cns.stats_start_time` 控制起算时刻 |
| GPU 编译 | gcc + CUDA 12.6, C++20, -O3, --use_fast_math, maxrregcount=255, AMREX_GPU_MAX_THREADS=256, RDC, AMREX_NO_PROBINIT=TRUE;现役可执行为 sm_80+sm_90 双架构 |

---

## 3. 2D 案(m119 / m142 / m180 通用配置)

### 3.1 几何与网格
| 项 | 值 |
|---|---|
| 坐标系 | 轴对称 RZ(`geometry.coord_sys = 1`),i→r,j→z;θ 动量恒等于 0(无旋流) |
| 域 | r ∈ [0, 0.3048 m] × z ∈ [0, 0.8128 m] = **12D × 32D** |
| 基础网格 | 192 × 512(dx0 = 1.5875 mm = **D/16**) |
| 层级 | max_level = 6 → 最细 **D/1024**(dx6 = 24.80 µm;每 D 1024 胞) |
| blocking_factor / max_grid_size | 32 / 256 |
| regrid_int / n_error_buf | 20 / 8(各级) |
| 边界码 | lo = (3, 1):r=0 对称轴,z-low 用户三区 BC;hi = (2, 2):r-high 与 z-high foextrap(零梯度外推) |

### 3.2 AMR tagging(prob.h `user_tagging`)
- 无量纲密度梯度 sensor:√[(∂ρ/∂r)² + (∂ρ/∂z)²]·(2Δ)/ρ;阈值 **0.05(lev<3)/ 0.03(lev≥3)**。
- 每级空间 caps(r 与 z 同时满足才允许 tag;tag 于 lev l 生成 lev l+1):
  r caps = 7D, 5D, 4D, 3D, 2D, 1.25D;z caps = 24D, 20D, 16D, 12D, 8D, 4D(lev0→5)。
- 强制区(与 sensor 取或,再与 caps 取与):
  - **A(→L4)**:z ≤ 8D 且 r ≤ 1.5D(射流核心);
  - **B(→L5)**:z ≤ 6D 且 r/D ∈ [0.35+0.04ζ, 0.65+0.12ζ](剪切层走廊,ζ=z/D);
  - **C(→L6)**:z ≤ 1D 且同走廊(唇口)。
- **m180 唯一差异**:sensor 加 `nt>0` 门控(fresh-start 初次建层仅用强制区,避免初始化伪影 seeding)。

### 3.3 探针(32 个,压力,level 0,每步采样;10+3+10+5+4=32)
- lp01–lp10:唇线 r ≈ 12.7 mm(=R)沿 z 10 站;
- ol1–ol3:外线 r ≈ 22.9 mm(0.9D)3 站;
- cl01–cl10:轴线 r ∈ [0, 3.2 mm] 沿 z 10 站(z 至 ~8D);
- f1a–f1e:r ≈ 38.1 mm(1.5D)5 站;f3a–f3d:r ≈ 76.2 mm(3D)4 站。
- 探针盒尺寸 3.2 mm(≈2 个 L0 胞)。输出 `probes.log`。

### 3.4 输出
plot_int = 50(Phase A 值;生产段以命令行覆盖,典型 250–500),check_int = 1000(生产段 250–500);plt/chk 落 `./plot/`。plt 含 20 分量 = 5 守恒 + 11 统计 + 4 派生(p/T/u/v)。

### 3.5 各案实际运行史(与配置区分)
| 案 | Phase A(发展段) | Phase B(统计段) | 统计窗 | 终帧 |
|---|---|---|---|---|
| m119 | 0→6.5 ms | 6.5→16.5 ms | 10.0 ms | plt32950 |
| m142 | 0→6.5 ms | 重启 chk13219(t=6.4998 ms)→15.75 ms | 9.25 ms | plt32080 |
| m180 | 0→6.5 ms(8×H100) | 6.5→10.0 ms | 3.5 ms | plt20700 |
- dt ≈ 0.49 µs(CFL 0.3, D/1024);~10M 胞级网格。8×H100 实测:无统计 ~0.95 s/step,统计+每 250 步 plt/chk I/O ~2.0–2.3 s/step。
- **网格履历**:各案发展段自较粗种子起步(如 m180 目录名 L3seed 所示的 L3 种子协议),此后抬升至 max_level=6;**最终 D/1024 分辨率并非从 t=0 启用**。
- 衍生分支 m142_teno5:`prob.h` 一行改 `ReconScheme::Teno5`,自 chk14000/chk14500 重启,统计窗 7.2→8.089 ms(box 到期截断)。

---

## 4. 3D 案(m142_3d,L5 生产 = l5_prod_cb_20260715)

### 4.1 几何与网格
| 项 | 值 |
|---|---|
| 坐标系 | Cartesian(coord_sys=0),**喷流轴 = x**;i 为内存连续方向 |
| 域 | x ∈ [0, 0.6096 m] × y,z ∈ [−0.2032, 0.2032 m] = **24D × (±8D)²** |
| 基础网格 | 192 × 128 × 128(dx0 = 3.175 mm = **D/8**) |
| 层级 | max_level = 5 → 最细 **D/256**(dx5 = 99.22 µm) |
| blocking_factor | 32 32 32 16 16 16(L0–L5) |
| max_grid_size | 64(L0)/ 512(L1+;静态大盒需大于默认上限) |
| n_error_buf | 8 8 8 6 4 |
| 边界码 | lo = (1,2,2):x-low 用户三区 BC,y/z-low foextrap;hi = (2,2,2) 全 foextrap |
| 注记 | 侧向 ±6D 恰等于法兰半径:对 mean-flow 分支平板等效无限大(已记录局限;声学分支需加宽)。**foextrap 不是声学无反射边界**:射流→侧边界→返回的传播时间 ≈16D/c∞ ≈ 1.18 ms,3.66 ms 统计窗内理论上可往返 ~3 次;平均激波场可接受,定量声学(screech 幅值/相位)需 NRBC/sponge 或 R(f,θ,m) 标定 |

### 4.2 静态层级(生产用 `amr.regrid_file`,冻结 BoxArray)
`fixed_grids_v5_L5.dat`(盒子按 level−1 粗指数写入,读取时 ×2 细化;`amr.regrid_int = 100000` 冻结,`amr.regrid_on_restart = 1` 建层;`DistributionMapping.strategy = SFC`):

| 级 | 盒数 | 胞数 | 分辨率 |
|---|---|---|---|
| L0 | (基础) | 3.15M | D/8 |
| L1 | 8 | 4.19M | D/16 |
| L2 | 8 | 4.19M | D/32 |
| L3 | 12 | 7.37M | D/64 |
| L4 | 44 | 19.40M | D/128 |
| L5 | 26 | 12.58M | D/256 |
| 合计 | 98(+基础) | **50.9M** | advance weight W = Σ2^ℓN_ℓ ≈ 0.80B |

层级布局(设计意图):L5 = 唇口/剪切层环形套筒(轴向至 x≈0.75D 附近),L4 = 核心+走廊(至 ~4D),L3 = 安全核心,L1/L2 = 声学/反馈区。**对称轴上无 L5 覆盖,shock-cell 列所在近轴区分辨率为 L4 = D/128,第 4 峰附近临 L4→L3 交界**(诊断口径注意)。

`prob.h` 内另有动态 `user_tagging`(固定 L1/L2/L3 区 + L6 lip collar 两种变体 + sensor 限 lev≤2),为 regrid_file 之前的发展阶段所用,生产段被静态层级取代。

### 4.3 入口宽带非轴对称 startup seed(一次性)
u′_x = ε·U_e·g(t)·h(r)·S(θ),其中:
- ε = 3.0e-4(**总方位角 RMS**/U_e,非单模态);
- S(θ) = Σ_{m=1..5} c_m cos(mθ+φ_m),c_m = {1, 0.7, 0.5, 0.35, 0.25}×√(2/1.925)(归一使 ⟨S²⟩=1),φ_m = {0.9273, 2.1416, 3.7851, 5.0119, 0.4636}(存档确定性种子,无手性、无载频);
- g(s) = 64s³(1−s)³ 的 C² 窗,**runtime 参数**:`prob.perturb_t_start = 5.0e-3`(生产实际),`prob.perturb_duration = 2.0e-4` s;
- h(r) = exp[−((r−R)/3δ)²] 唇口高斯包络。
一阶质量/动量净通量为零(方位角正交);离散孔残差 O(ε²) 由 banner 审计打印。重启不会自动重放窗口(需显式设 t_start)。

### 4.4 探针(定义 216,激活 210;压力,level 1,每 10 步)
- **12 个方位角环 × 16 探针 = 192**(s1–s4、o1–o4、s05、s6、n15、n30 站),供 azimuthal Fourier(m-mode)分解;
- **18 个单点**:q1–q5 = **近唇线**(y ≈ R = 12.7 mm)沿 x 的轴向阵列(非轴线);a1–a3 = 轴线组(y≈z≈0);r1a–d / r2a–d = 径向阵列;m1–m2;
- 6 支 legacy 探针(P_cl_0p5De/1De/2De/4De、P_sh_1De/2De)已定义但未列入 `cns.time_probes`,不采样。
输出 `probes3d_prod.log`。

### 4.5 输出与统计
plot_int = 500,check_int = 2500(生产链覆盖为 check_int=500/plot_int=500);`cns.stats_start_time = 6.5e-3`。plt 25 分量 = 5 守恒 + 15 统计 + 5 派生(p/T/u/v/w)。

### 4.6 生产运行史
| 段 | 区间 | 说明 |
|---|---|---|
| 种子 | v5_staged/chk04500 → chk04520(t = 5.0 ms) | 静态层级建立与校验(与 5 级网格表逐盒核对通过) |
| seed 窗 | 5.0–5.2 ms | 宽带非轴对称扰动(上表参数) |
| develop | 5.0 → 6.5 ms(chk05847),cfl = 0.30,stats off | 8×H100(CB 裸机) |
| statistics | 6.5 → 10.16 ms(chk09170/plt09170) | 统计窗 **3.66 ms**;末段中位速度 ≈ 4.67 s/step;dt ≈ 1.1 µs |
- **网格履历**:L5 静态层级自 t = 5.0 ms(chk04520)起建立;此前发展段来自较粗种子。**最终 D/256 层级实际运行 5.0→10.16 ms ≈ 5.16 ms**,并非从 t=0 启用。
- 衍生分支 teno5_run:`prob.h` 一行改 Teno5,单腿协议(chk05847 重启,stats_start_time=7.0e-3,stop 10.66e-3)备好未跑。

---

## 5. 位置索引
| 内容 | 路径(/shared = AWS 头节点导出) |
|---|---|
| 2D 案目录 | /shared/cerisse_reflux/exm/underexpanded_jet/2d/{m119_gpu_prod, m142_gpu_prod, m180_L3seed} |
| 2D TENO 分支 | .../2d/m142_teno5 |
| 3D 案目录 | /shared/cerisse_reflux/exm/underexpanded_jet/3d/m142_3d_gpu(生产链 l5_prod_cb_20260715;TENO 分支 teno5_run) |
| 静态网格文件 | .../3d/m142_3d_gpu/profile_20260715/fixed_grids_v5_L5.dat |
| 分析产物(本机) | ~/testcerisse/cerisse/analysis/panda_benchmark/ |
| 到手即算工具 | /shared/cerisse_reflux/tools/{go_h100.sh, new_p5_bootstrap.sh, gpubind_auto.sh} + /shared/cuda-12.6 |
