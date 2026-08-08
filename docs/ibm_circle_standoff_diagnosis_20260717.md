# Mach 4 二维圆柱 GPIBM 脱体弓形激波审计

日期：2026-07-17

## 强制停止规则

Mach 4 圆柱是所有固定壁面 GPIBM 修改的首要回归基准。任何 IBM 几何、ghost-state
extension、壁面边界条件、近壁通量、WENO/TENO、黏性通量、AMR 或载荷修改，在进入
其他验证前必须重新通过本报告定义的 `N=320` 稳态圆柱门槛。

出现以下任一情况时立即停止其他工作，只处理圆柱弓形激波：

1. 停滞线上不能识别独立的 detached bow shock；
2. 最后稳态窗口内的 `Delta/R` 未平台化；
3. 密度梯度与压力梯度给出的 shock location 相差超过 0.5 个网格；
4. `Delta/R` 相对冻结的 `N=320` 软件回归值偏差超过 2%；
5. 出现负密度、负压力、负内能或非有限状态；
6. 最终密度、压力和数值 schlieren 图不能显示连续、对称、与圆柱分离的弓形激波。

自动门槛由 `tools/check_ibm_circle_gate.py` 执行。第 5、6 项仍需由运行日志和最终场图
审计；只有中心线标量通过不足以证明 shock shape 正确。

## 当前结论

当前可以严格区分两个 GPIBM 边界算子：

- 原来的两阶段高阶路径
  `fluid support -> image-point WLS -> BI boundary condition -> GP extrapolation`
  在 `iorder=2, eorder=1, alpha=0.6` 时，即使运行到
  `t=0.00537`、`t*=t U_inf/D=37.3`，仍不能建立正确的圆柱脱体弓形激波。
- 新的 BI-centered constrained shared-GP 路径保持标准 ghost-cell GPIBM：每个真实
  solid-side ghost cell 只有一个 GP、一个 boundary intercept (BI)、一个壁面法向和一个
  primitive state；边界条件直接进入 BI-centered polynomial，然后在该真实 GP 中心直接
  求值。该路径能够形成稳定、对称、随网格单调移动的 detached bow shock。

因此，已经闭合的是“弓形激波完全消失”的灾难性故障；尚未闭合的是圆柱 shock
stand-off 的最终物理精度。当前方法通过软件回归门槛，但仍缺同方程、同通量、同外域的
body-fitted reference，不能据此宣称已经完成 SRP 生产认证。

## 固定数值协议

- 方程：二维可压缩 Euler 方程，量热完全理想气体，`gamma=1.4`；
- 来流：`M_inf=4`，uniform freestream initialization；
- 几何：固定圆柱，`R=0.1`，计算域 `[-0.5,0.5]^2`；
- 网格回归点：`N=320`，`D/h=64`，`max_level=0`，一个大 AMReX box；
- 空间离散：characteristic LLF flux splitting + WENO-Z5；
- IBM：shared-GP ghost-cell IBM，`ghost_layers=1`，stationary Euler-slip wall；
- 禁用：face-local state、direct crossing flux、cut-cell/EB、shifted shell；
- 时间积分：SSPRK(4,3)，`CFL=0.3`；
- 稳态判定：至少 `t*>=40`，最后四个 shock reports 的 `Delta/R` 相对极差不超过
  `2e-3`；
- shock extraction：停滞线密度与压力梯度峰值的亚网格拟合，并以 tanh transition fit
  交叉检查；
- 当前 BI-constrained production executable 的 source hash：`340a03a3ff5f8b81`。

所有后续圆柱回归必须输出上述完整配置。仅写“LLF-WENO-Z5”不是完整算法标识。

## 固定时刻与稳态结果

早先把 `t=0.002` 当作最终验收时刻是错误的。两种边界算子建立弓形激波的瞬态速度
不同，因此固定时刻会保留初始化记忆。

| 路径 | `Delta/R` at `t=0.002` | steady time | steady `Delta/R` |
|---|---:|---:|---:|
| legacy `i1/e1/alpha=0.6/pzero` | 0.5159359 | 0.008 | 0.5712879 |
| BI-constrained shared-GP | 约 0.560 | 0.008 | 0.5643135 |

稳态 BI-constrained 结果相对稳态 legacy 结果为 `-1.2208%`，密度/压力 shock-location
差为 `0.1427h`。两者都形成清晰弓形激波。旧 legacy 结果在 `t=0.002` 与 Billig 关系
接近，只是瞬态时刻的偶然相符，不能继续作为物理精度证明。

冻结的软件回归参考采用 BI-constrained `N=320` 稳态值：

\[
(\Delta/R)_{reg}=0.5643135034711624.
\]

它用于检测软件回归，不代表精确 Euler 解。

## BI-constrained 网格序列

| `N` | `D/h` | final time | steady `Delta/R` | 状态 |
|---:|---:|---:|---:|---|
| 192 | 38.4 | 0.008 | 0.5783158 | detached、对称、稳态 |
| 320 | 64.0 | 0.008 | 0.5643135 | detached、对称、稳态 |
| 512 | 102.4 | 0.006 | 0.5557441 | detached、对称、稳态 |

shock stand-off 随流体网格加密单调下降。三点自由阶拟合给出诊断性
`p=0.829`、`(Delta/R)_{h->0}=0.53776`。该拟合只有三个点和三个参数，没有剩余自由度，
不能作为严格的离散误差不确定度。固定 `p=1` 和 `p=2` 外推分别约为 0.54243 和
0.55348，说明当前序列尚不足以给出唯一连续极限。

圆柱折线几何随网格自动从约 600 个顶点增加到 800 个顶点。其 sagitta 远小于 `h`，
预计不是当前主误差，但严格的 fluid-grid convergence 仍应改为所有网格使用同一份高精度
圆柱几何。

## 经验关系的正确用途

Billig 圆柱经验关系为

\[
\Delta/R=0.386\exp(4.67/M_\infty^2),
\]

在 Mach 4 时为 0.516829。它来自实验数据相关，不是本计算域 Euler 方程的精确解。
[Billig 1967, DOI 10.2514/3.28969](https://doi.org/10.2514/3.28969)

Hornung 基于一组 perfect-gas Euler computations 给出的圆柱拟合为

\[
\Delta/R=2.14\epsilon(1+\epsilon/2),\qquad
\epsilon=\frac{\gamma-1+2/M_\infty^2}{\gamma+1}.
\]

本算例代入 `gamma=1.4` 得 0.519326，但原论文圆柱数据覆盖到约 `gamma=1.3`，这里包含
小幅参数外推。[Hornung 2021](https://authors.library.caltech.edu/records/chnt3-fpg55/latest)

因此，相关式只能用于异常检测。最终物理误差必须相对 matched body-fitted 或独立收敛
Euler reference 定义。

## 原 IP-WLS 路径为什么被否决

固定 LLF-WENO-Z5、shared-GP、`N=320` 后的单因素矩阵为：

| Case | `iorder` | `eorder` | `alpha` | pressure closure | fixed-time result |
|---|---:|---:|---:|---|---:|
| L | 1 | 1 | 0.6 | zero-gradient | `Delta/R=0.515936` |
| A | 2 | 1 | 0.6 | zero-gradient | 无 detached bow shock |
| B | 1 | 2 | 0.6 | zero-gradient | `Delta/R=0.560332` |
| C | 1 | 1 | 1.0 | zero-gradient | 无 detached bow shock |
| D | 1 | 1 | 0.6 | shock-aware | 与 L 相同 |

首步 GP 审计表明：

- `iorder=2` 的二次 3x3 image-point WLS 使用非凸权重，典型 `L1` 范数约
  `1.63--1.85`；
- 在弓形激波和停滞压缩层内，非凸插值会放大支撑点间的不连续变化；
- 即使把 IP placement 恢复到历史距离，二次 WLS 仍不能形成弓形激波；
- 线性 WLS 能恢复 shock blocking，但不构成一致的高阶 wall closure；
- 长时间复跑证明二次 IP-WLS 不是单纯建立激波较慢，而是没有产生正确的钝体阻塞。

这解释了 smooth MMS 可以给出良好形式阶，而 shock-wall interaction 仍失败：多项式
再现性不等于对非光滑数据的 monotonicity、bound preservation 或 nonlinear robustness。

## 当前 GPIBM 边界闭合

对每个真实 solid-side GP，使用该 GP 自己唯一的 BI、法向和可见流体支撑，求解

\[
\min_{\mathbf a}\|W(A\mathbf a-\mathbf q_f)\|_2^2,
\qquad C_{BI}\mathbf a=\mathbf g_B,
\]

并直接在真实 GP 中心求值

\[
q_G=P_{\mathbf a}(\mathbf x_G).
\]

对 stationary Euler-slip wall，在 BI 施加 `u_n=0`；压力、温度和切向速度采用相应的
齐次法向条件；primitive state 完成后由 EOS 重构密度和总能量。该方法仍是 shared-GP
ghost-cell IBM，不是 EB、cut-cell 或 face-local Riemann closure。

临时测试过仅对 `u_n` 使用一次 constrained polynomial 的低阶版本。它在 `N=320`
给出 `Delta/R=0.56591`，与当前二次 constrained 结果只差约 0.30%。因此高 `L1` 的
quadratic Dirichlet functional 会放大启动瞬态，但不是稳态 stand-off 偏差的主因；该临时
分支已经从源码删除。

## 自动回归命令

```bash
python3 tools/check_ibm_circle_gate.py \
  IBM/cases/ibm_tests/2d_bvh_gpu/runs/llf_wenoz5_sharedgp_20260717/history_t008_bicgp_g1_N320 \
  --output circle_gate.json
```

默认标准：

- `D/h=64`；
- `t*>=40`；
- 最后四个样本的 `Delta/R` 相对极差 `<=0.2%`；
- stand-off 至少 4 个网格；
- 密度/压力 shock location 差 `<=0.5h`；
- 相对冻结 `N=320` 参考差 `<=2%`。

脚本返回非零退出码即表示 **FAIL**。FAIL 后不允许继续其他 IBM 开发。

## 当前生产判定与下一步

- BI-constrained shared-GP + LLF-WENO-Z5：通过当前 Mach 4 圆柱软件回归门槛；
- 原二次 image-point WLS：强激波壁面问题明确不支持；
- LLF-TENO5、AFD、face-local、direct-crossing、shifted-shell：本报告不认证；
- 当前圆柱物理 stand-off 精度：尚未认证；
- 当前 SRP 生产状态：不能仅依据本报告放行。

下一步只允许围绕圆柱完成：

1. 固定同一高精度圆柱几何，补 `N=768/1024` 或达到明确 asymptotic trend；
2. 建立 matched body-fitted Euler reference；
3. 比较完整 shock shape、停滞线剖面、`Cp(theta)`、drag 和 BI 上的 `u_n`；
4. 对最终候选重复 grid-phase、CPU/GPU 和多 FAB/MPI 回归。

只有圆柱硬门槛持续通过并完成物理参考后，才恢复 compression corner、
Prandtl--Meyer expansion、double ramp、NS、AMR 和 SRP。

## 数据产物

主目录：

`IBM/cases/ibm_tests/2d_bvh_gpu/runs/llf_wenoz5_sharedgp_20260717/`

关键子目录和图：

- `history_t008_L_g1_N320/`：稳态 legacy 历史；
- `history_t008_bicgp_g1_N320/`：稳态 BI-constrained 历史；
- `history_t008_ipwls_i2e1a0p6_g1_N320/`：失败的二次 IP-WLS 长时间复跑；
- `convergence_bicgp/N192`, `N512`：网格序列；
- `standoff_history_t008_L_vs_bicgp_g1_N320.png`：瞬态与稳态对比；
- `convergence_bicgp/standoff_convergence.png`：网格趋势；
- `circle_gate_bicgp_vs_legacy_N320.json`：已通过的门槛记录。
