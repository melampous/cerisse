# Run165 RZ：HLLC–MUSCL 与 AFD–HLLC–WENO-Z5 逐点崩溃诊断

## 结论

本次失败分支是 `AFD–HLLC–WENO-Z5 + cns.afd_shock_llf=0`。当前默认的
`cns.afd_shock_llf=1` 没有在同一门禁中失败；它在问题面切换为高阶特征
LLF。HLLC–MUSCL O2 也稳定通过。

根因不是 IBM、RZ 轴线、物理边界、AFD 修正项或 HLLC 接触波分母奇异。
根因是 AFD 路径采用逐点变化的非线性“特征”变换，并对两个声学变量分别
使用不同 WENO 权重。在非光滑模板上，逆变换得到一个压力仍为正、但速度
严重越界的面状态。当前合法性检查只检查有限性以及 `rho>0, p>0`，因此这个
状态进入 HLLC。面波速远大于计算时间步时使用的单元中心波速，最终产生约
正常值一千倍的共享面能量通量，使 RK2 第二阶段的一个单元出现负内能。

## 完全一致的 A/B 条件

- 二维 RZ Euler，无黏性、无 LES、无 IBM、无 solid cell。
- uniform L1：`64 x 192 = 12,288` 个有效单元。
- `dr=dx=1.06129666667 mm`。
- 固定粗层时间步 `0.4 us`，L1 细层时间步 `0.2 us`。
- SSPRK2，Run165 虚拟喷口出口，jet ramp `100 us`。
- 两个格式从同一干净初始场开始；比较时间严格同步。
- 失败位置：L1 cell `(i,j)=(4,171)`，`r=4.775835 mm`，
  `x=-19.3153128611 mm`。它不在轴线、不在物理边界，也没有 IBM。

## 差异从哪里、什么时候开始

以 `AFD(LLF off) - MUSCL` 的内能密度差除以来流内能密度：

| 阈值 | 首次时间 | 单元数 | 首次区域 |
|---:|---:|---:|---|
| `1e-12` | `0.4 us` | 71 | `r=0.531–11.144 mm`, `x=-6.580–1.911 mm` |
| `1e-6` | `0.4 us` | 20 | `r=0.531–6.898 mm`, `x=-0.212–1.911 mm` |
| `1e-3` | `0.8 us` | 6 | `x=0.849 mm` 的喷口邻域 |
| `1e-2` | `5.2 us` | 6 | `x=1.911 mm` 的出口/唇口模板 |
| `1e-1` | `5.6 us` | 5 | `x=1.911 mm` 的出口/唇口模板 |

因此两种格式从第一步就在喷口出口/唇口 stencil 产生预期的离散差异。该差异
随逆向喷流及压缩结构向上游传播，并非两个初场不一致。

## 崩溃单元的局部前兆

cell `(4,171)` 的内能密度 `rho*e`（J/m3）：

| 时间 [us] | MUSCL | AFD, LLF off | AFD, LLF on |
|---:|---:|---:|---:|
| 60.0 | 89,866 | 88,577 | 91,805 |
| 68.0 | 73,427 | 81,086 | 72,167 |
| 70.0 | 48,758 | 58,612 | 47,462 |
| 72.0 | 29,276 | 29,828 | 28,826 |
| 72.4 | 27,741 | 24,199 | 26,998 |
| 72.8 | 26,711 | 19,240 | 25,801 |
| 73.2 | 26,055 | 14,359 | 25,287 |
| 73.6 | 25,788 | 10,846 | 24,785 |

明显的破坏性分离从约 `72.4 us` 开始。`73.6 us` 时 AFD-off 仍是合法状态，
但局部内能储备已经只有 MUSCL 的约 42%。

## 最后一个细时间步的逐阶段因果链

L1 step 369 的第一阶段从 `73.6 us` 推进到 `73.8 us`：

- RK2 输入 cell `(4,171)`：`rho*e=10,846.50 J/m3`。
- stage-1 Forward Euler 后：`rho*e=10,709.70 J/m3`，仍为正。
- stage-2 RHS 在共享轴向面 `face=(4,171), dir=1` 产生异常通量。
- RK2 最终 cell `(4,171)`：`rho=0.0739947 kg/m3`，
  `rhoE=-9.93217e6 J/m3`，动能密度 `2.64098e7 J/m3`，
  因而 `rho*e=-3.63420e7 J/m3`。

异常 RHS 在共享面的两个相邻单元中表现为近似等大反号：

| cell | density RHS | axial-momentum RHS | energy RHS |
|---|---:|---:|---:|
| `(4,170)` | `+2.15465e6` | `-2.24692e10` | `+1.00699e14` |
| `(4,171)` | `-2.37885e6` | `+2.26230e10` | `-1.00774e14` |

这证明它是一个守恒共享面通量的异常交换，不是非守恒 RZ source、BC source
或 IBM source。

## 问题面直接对照

第二阶段在 `face=(4,171), dir=1`：

| 量 | AFD–HLLC, LLF off | HLLC–MUSCL | AFD shock-LLF |
|---|---:|---:|---:|
| reconstructed right `rho` [kg/m3] | 0.307519 | 0.388782 | 不使用该 HLLC 状态 |
| reconstructed right `p` [Pa] | **56.1402** | 10,041.9 | — |
| reconstructed right `u_n` [m/s] | **-8,861.0** | -734.411 | — |
| maximum face wave speed [m/s] | 8,876.99 | 929.796 | cell-centred LLF alpha |
| face CFL `dt*a/dx` | **1.67286** | 0.175219 | 与 cell-centred CFL 同量级 |
| energy flux [W/m2] | **-1.07012e11** | -1.06249e8 | -1.34335e8 |

异常 HLLC 能量通量是 MUSCL 的 `1007.19` 倍，是 shock-LLF 的 `796.61` 倍。

HLLC 的接触波计算本身没有接近零的分母：

- `SL=-8876.99 m/s`, `SR=-645.508 m/s`, `S*= -4305.84 m/s`；
- `S*` 分母 `-5820.43`；
- `SL-S*=-4571.15`, `SR-S*=3660.33`。

因此不能把这次失败归因于 HLLC 分母除零。

## 为什么 WENO 会得到 -8861 m/s

AFD 当前把每个 stencil 点的 primitive state 分别变换为

```text
w1 = 0.5 * (p + sqrt(gamma*rho*p) * un)
w2 = 0.5 * (p - sqrt(gamma*rho*p) * un)
```

然后 `w1` 和 `w2` 各自计算 WENO-Z 权重。问题面右状态使用的五点模板中：

| 分量 | WENO 权重 `(omega0,omega1,omega2)` | 重构值 |
|---|---|---:|
| `w1` | `(0.19323, 0.58474, 0.22203)` | `-21753.52` |
| `w2` | `(0.80747, 0.11240, 0.08012)` | `+21809.66` |

两个声学变量选择了完全不同的 stencil 组合。逆变换时：

```text
p = w1 + w2 = 56.14 Pa
w1 - w2 = -43563.18
sqrt(gamma*rho*p) = 4.9163
un = (w1-w2)/sqrt(gamma*rho*p) = -8861.0 m/s
```

压力恰好仍为正，所以正密度/正压力检查不会拒绝它；但小正压力使逆变换分母
变小，速度被放大。标准 frozen characteristic WENO 应在整个面 stencil 上使用
同一个线性特征基及其逆矩阵。当前逐点变化的非线性变换没有这种一致性，
这是本次异常重构的核心。

此外，时间步由 cell-centred 最大特征速度计算，没有看到重构面的
`8877 m/s`，所以实际面 CFL 从目标约 0.2 跳到 1.67。

## AFD correction 是否负责

不是。追踪结果明确给出：

```text
smooth=0, shock_candidate=1
final_branch=base_hllc
afd_correction_applied=0 reason=nonsmooth_stencil
```

问题面已经被识别为非光滑/压缩面；当 `afd_shock_llf=0` 时，代码仍允许点值
WENO 状态进入 HLLC。当默认 `afd_shock_llf=1` 时，同一面选择高阶特征 LLF，
该步以及后续门禁通过。

## 修复优先级

1. 生产计算保持 `cns.afd_shock_llf=1`。这是本算例已经验证的最短止血措施。
2. 把 AFD primitive/characteristic reconstruction 改成面冻结的线性特征基，
   不要分别使用逐点变化的非线性声学尺度再做非线性逆变换。
3. 在所有 HLLC 路径加入 reconstructed-state envelope/face-CFL 守卫；若面状态的
   `|u|+c` 相对 stencil cell-centred 上界异常放大，守恒地回退到 HLLE/LLF。
4. 在重构后加入 Zhang–Shu 型 primitive/conservative positivity scaling；仅检查
   `rho,p>0` 不够，因为本次状态虽正但极端。
5. 为每个 SSP Forward-Euler bracket 加共享面守恒 positivity limiter，作为最终
   数学保险。不能用单元独立 clipping 修复，因为那会破坏守恒。

第 1 条是当前已通过 A/B 的配置结论；第 2–5 条需要独立回归后才能进入生产。
本次诊断没有修改 IBM 状态、GP、BI-CWLS、full-cell Cartesian RHS 或 IBM flux。

## 证据文件

- `hllc_muscl_vs_afd_collapse_mechanism.png`：局部前兆、异常面状态及通量图。
- `pointwise_internal_energy_difference_maps.png`：61 个同步时刻的逐单元差值传播。
- `difference_growth_timeline.png`：全域 L1/L2/Linf 差值增长。
- `precollapse_target_neighborhood.png`：崩溃前目标单元邻域。
- `pointwise_difference_fields.npz`：61 个同步时刻的完整逐单元差值数组。
- `field_difference_timeline.csv`：全域逐时刻范数、极值位置和阈值包围盒。
- `point_target_4_171_timeline.csv`：目标单元三分支时间序列。
- `afd_fallback_off/face_trace_chk184_hllc.log`：异常 HLLC 面状态和通量。
- `afd_fallback_off/face_trace_chk184_llf.log`：同一 checkpoint 的 LLF 回退 A/B。
- `muscl/face_trace_150_185.log`：MUSCL 对应面状态和通量。
