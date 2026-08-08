# RZ 轴邻径向面失稳修复方案:传感器门控的一阶度量一致 LLF 回退

**Sensor-gated first-order metric-consistent local Lax–Friedrichs (LLF) fallback at the innermost interior radial face**

日期:2026-08-07(UTC)。状态:已在 GPU A/B 对照中验证，并已合入当前
`src/rhs/Weno.h`。当前实现在保留轴邻专用门控的同时，将同一
局部 LLF 回退扩展到 Cartesian 和 R-Z，all-fluid 和 pure shared-GP，
以及 LLF-WENO-Z5、LLF-TENO5 和 LLF-TENO6。
本文档面向第三方实现者(含其他 AI 会话),给出问题诊断结论、修复设计、原型补丁的演变、参数语义、验证记录,以及向 AFD / HLLC / skew-symmetric 等其他 RZ 格式移植的原则。

---

## 1. 适用对象与问题概述

- **代码**:新求解器树 `src/rhs/Weno.h` 中的特征式通量矢量分裂(characteristic-wise flux-vector splitting)WENO-Z5 驱动器(`weno_t` / `compute_characteristic_llf_face_fluxes`),其 RZ(axisymmetric)分支采用**径向度量加权通量分裂**(metric-weighted flux split):径向面重构度量变量 `r·(U ± F/α)`,面数组存储 `h_(rF)/r_face`(轴面槽位存原始 h-通量),配对压力通道(paired-pressure)单独处理径向动量的压力部分。特征常量:`rz_radial_axis_face_flux_is_metric = true`。
- **故障现象**:2D RZ 欠膨胀射流(underexpanded jet)启动瞬态(starting transient)在 t ≈ 8×10⁻⁵ s 处触发负内能中止(negative internal energy abort),失稳单元恒为**轴邻第一格**(axis-adjacent cell,i = 0,r = Δr/2),位于接触面(contact front)越过轴线的轴向位置。CPU/GPU 崩溃点逐位一致;与 WENO ε 参数无关;粘性项无关(Euler 变体同样中止)。
- **根因(控制变量实验证实,非假设)**:
  1. 轴邻单元的更新由**唯一的内侧径向面**(face i = 1)主导,度量散度含 `1/r_c = 2/Δr` 放大因子;
  2. 间断越轴时,r-调制样本(r-modulated samples)的 WENO 非线性权重产生的面重构误差被该因子放大,并全部注入径向动量(实测虚假向轴速度 u_r ≈ −720 m/s,而对称性要求该处 u_r = O(Δr)),动能超过总能 → ρe < 0;
  3. 轴镜像模板臂(带符号半径 + 宇称 ghost)的度量变量延拓经代码审读**解析正确**——宇称 bug 假设已排除,脆弱性为结构性。
- **定位实验**:仅将 `reconstruct_radial_metric_flux` 与 `rz_radial_axis_face_flux_is_metric` 置 false(其余重构机制不变)即稳定 400 步 → 度量加权路径为失稳原因;仅在最内侧 K=1 径向面改用普通分裂即稳定(静态 collar 实验),但 RZ MMS 阶数从五阶退化到 ~2 阶 → **五阶精度的成败系于轴邻面,修复必须是动态门控而非静态降阶**。

## 2. 修复设计

在**最内侧内部径向面**(`face_index[0] == 1`,不含轴面本身)安装密度跳传感器(density-jump sensor):当面两侧相邻单元的相对密度跳超过阈值(0.3%)时,将该面的 r-调制 WENO 重构结果**整体替换**为一阶度量一致 LLF 通量:

```
F_face = 0.5 · (r_L·F_adv,L + r_R·F_adv,R) / r_face  −  0.5 · α · (U_R − U_L)
```

其中 `F_adv = F − p·e_UMX`(径向动量分量剥离压力,遵循 paired-pressure 约定),压力面值取算术平均 `p_face = 0.5(p_L + p_R)` 并写入 `radial_pressure_face_flux`;α 为该面已算好的 LLF 标量最大波速 `llf_max_wave_speed`。

**两条不可违背的性质**(均由失败实验反向确立):

1. **度量一致性(metric consistency)**:回退通量的对流部分必须保留 r-加权平均(`r_L·F_L + r_R·F_R` 除以 `r_face`),因为 `compute_rhs` 的散度装配按度量 h-通量约定读取该面。改用普通(非加权)分裂形式的回退在喷口区第 8 步即崩——通量约定不一致比不修更糟。
2. **耗散项用物理跳(physical jump),不加权**:LLF 耗散必须写 `α·(U_R − U_L)`,**不得**对 U 做 r-加权。若误写为 `α·(r_R U_R − r_L U_L)/r_face`,均匀流场下耗散项不为零,等效于持续抽取质量(实测均匀区 ρ 从 1 降至 0.12,collar 越宽崩得越早)。这是本修复中最隐蔽的错误来源。

**传感器阈值**:相对密度跳 `|ρ_L − ρ_R| / min(ρ_L, ρ_R) > 0.003`。阈值 0.01 曾止住首个崩溃点(step 348)但在 step 612 出现残余崩溃(ρe = −4256):阈值偏高使传感器间歇关闭,度量重构误差在"关闭步"内时间累积过冲;降至 0.003 使间断经过期间持续触发,残余崩溃消失。光滑流(MMS)下传感器零触发,五阶精度逐位保留。

## 3. 代码位置与逐处修改

> **历史记录，不要再次应用。** 本节保留最初只覆盖轴邻面的验证
> 补丁，便于追溯失稳诊断。当前生产实现已经重整到
> `characteristic_llf_flux_t::compute_characteristic_llf_face_fluxes` 的
> 共享面循环中，并在高阶重构之前完成传感器判断。当前代码直接读取
> 相邻状态并在触发后返回，因此不再使用下述四组临时数组和事后覆写
> 结构。实际实现与参数以当前 `src/rhs/Weno.h`、`docs/input.md` 和
> `docs/llf_local_fallback_validation_20260807.md` 为准。

**文件**:`src/rhs/Weno.h`
**函数**:`compute_characteristic_llf_face_fluxes(...)`(public 静态模板成员,含 `amrex::ParallelFor` 设备 lambda;四处插入**全部位于该 lambda 体内**)。
**基线**:补丁对当前本地工作树的 `Weno.h` 可无冲突应用(已核对:本地文件与验证版仅差这 63 行新增)。
**归档**:`patches/rz_axis_fix_20260807/rz_axis_fix.patch`(ed-diff 格式)与 `patches/rz_axis_fix_20260807/Weno_FIXED_verified.h`(验证版完整头文件,可直接 diff 定位)。

### 插入点 1:常量声明(lambda 体内,靠近入口)

锚点:lambda 内 `amrex::ignore_unused(ibm_markers, face_flux);` 之后、`Real radial_face = Real(1.0);`("In R-Z, reconstruct the metric flux…" 注释块)之前。

```cpp
// Near-axis metric-flux robustification constants (lambda-local
// to avoid a constexpr-if device capture): first-order metric LLF
// fallback on the innermost interior radial face. threshold < 0 =
// always on that face; a positive value gates by density jump.
constexpr int rz_axis_collar_faces = 1;
constexpr Real rz_axis_shock_sensor_threshold = Real(0.003);
```

⚠️ **NVCC 约束**:这两个常量必须声明在 lambda 体内(lambda-local)。若声明为函数级 `constexpr` 再由设备 lambda 经 constexpr-if 分支引用,nvcc 的扩展设备 lambda 捕获规则会报错。这与本树此前已修的"函数内 static ε 旋钮须以局部值捕获"是同一类限制。

### 插入点 2:回退存储数组

锚点:`sample_metric_scales` 初始化循环之后、主采样循环(`for (int sample_index = 0; …)` 遍历 `2 * required_ghost_cells` 个模板样本)之前。

```cpp
// First-order metric fallback storage (immediate L/R cells).
Real fo_cons[2][Closure::NCONS];
Real fo_flux[2][Closure::NCONS];
Real fo_density[2] = {Real(0.0), Real(0.0)};
Real fo_pressure[2] = {Real(0.0), Real(0.0)};
Real fo_radius[2] = {Real(0.0), Real(0.0)};
```

### 插入点 3:面两侧紧邻单元数据捕获(主采样循环体内)

锚点:采样循环内 `sample_metric_scales[sample_index] = radial_metric_scale;` 之后。此处 `conservative_state` / `physical_flux` 是刚为当前样本算出的**原始**(未分裂、未加权)单元中心守恒量与方向物理通量;`sample_index == required_ghost_cells − 1` / `required_ghost_cells` 恰为面左、右紧邻单元。

```cpp
if (reconstruct_radial_metric_flux &&
    (sample_index == required_ghost_cells - 1 ||
     sample_index == required_ghost_cells)) {
  const int fo_side = sample_index - (required_ghost_cells - 1);
  fo_density[fo_side] =
      primitive_states(stencil_cell_index, Closure::QRHO);
  fo_pressure[fo_side] =
      primitive_states(stencil_cell_index, Closure::QPRES);
  fo_radius[fo_side] = radial_origin +
      (Real(stencil_cell_index[0]) + Real(0.5)) * radial_spacing;
  for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
    fo_cons[fo_side][fo_c] = conservative_state[fo_c];
    fo_flux[fo_side][fo_c] = physical_flux[fo_c];
  }
}
```

### 插入点 4:传感器判定与通量覆写(paired-pressure 度量分支尾部)

锚点:`if (reconstruct_radial_metric_flux) { … }` 高阶度量分支(即计算 `reconstructed_advective_flux` / `reconstructed_pressure_flux` 并以 `return;` 结束的那个分支)内,轴槽位三元式写入
`face_flux(face_index, Closure::UMX) = radial_face > 0 ? metric_advective_flux_per_radius + pressure_face : 0;`
**之后、`return;` 之前**。高阶重构照常完整执行,回退仅在触发时整面覆写——光滑流下代码路径与原版完全一致(这是 MMS 逐位不变的原因)。

```cpp
if (radial_face > Real(0.0) && face_index[0] >= 1 &&
    face_index[0] <= rz_axis_collar_faces) {
  const Real fo_dmin = amrex::min(fo_density[0], fo_density[1]);
  const Real fo_djump = (fo_dmin > Real(0.0))
      ? std::abs(fo_density[0] - fo_density[1]) / fo_dmin
      : Real(1.0);
  if (fo_djump > rz_axis_shock_sensor_threshold) {
    const Real fo_inv_rface = Real(1.0) / radial_face;
    const Real fo_alpha = llf_max_wave_speed;
    Real fo_full[Closure::NCONS];
    for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
      Real fo_adv_l = fo_flux[0][fo_c];
      Real fo_adv_r = fo_flux[1][fo_c];
      if (fo_c == Closure::UMX) {
        fo_adv_l -= fo_pressure[0];
        fo_adv_r -= fo_pressure[1];
      }
      // Metric-weighted flux average with an UN-weighted
      // physical-jump LLF dissipation (vanishes for a uniform
      // state; a metric-weighted dissipation would not).
      fo_full[fo_c] = Real(0.5) * fo_inv_rface *
          (fo_radius[0] * fo_adv_l + fo_radius[1] * fo_adv_r) -
          Real(0.5) * fo_alpha *
          (fo_cons[1][fo_c] - fo_cons[0][fo_c]);
    }
    const Real fo_pressure_face =
        Real(0.5) * (fo_pressure[0] + fo_pressure[1]);
    for (int fo_c = 0; fo_c < Closure::NCONS; ++fo_c) {
      face_flux(face_index, fo_c) = fo_full[fo_c];
    }
    face_flux(face_index, Closure::UMX) =
        fo_full[Closure::UMX] + fo_pressure_face;
    radial_pressure_face_flux(face_index, 0) = fo_pressure_face;
  }
}
```

要点说明:

- `radial_face > 0` 排除轴面本身(face 0)——轴面的度量 h-通量处理保持原样,collar 只覆盖 face 1..K;
- 覆写包含三件事:全部守恒分量的 `face_flux`、径向动量分量补回 `p_face`、`radial_pressure_face_flux` 同步覆写——三者缺一即与 `compute_rhs` 的装配约定失配;
- `fo_alpha` 复用该面已算好的 LLF 标量波速,不另行估计。

## 4. 参数语义

下表记录原始轴邻原型中的编译期常量。当前实现已改为运行时参数，
并扩展到共享驱动的全部 Cartesian 和 R-Z 面。当前默认值及完整语义
见 `docs/input.md`。

| 常量 | 值 | 语义 |
|---|---|---|
| `rz_axis_collar_faces` | 1 | 回退覆盖 face 1..K;K=1 已充分(失稳完全局域于第一内侧径向面);K=2 亦稳定但无必要 |
| `rz_axis_shock_sensor_threshold` | 0.003 | 相对密度跳阈值;负值 = 该面无条件回退(静态 collar,应急开关);0.01 不足(间歇门控致时间累积残崩) |

## 5. 验证记录(2026-08-07,GPU A/B)

- **鲁棒性**:原失稳射流算例(baseline 2D RZ,启动瞬态)越过原崩溃点(step 114 / step 348)至 step 500+ 零中止;0→1.8 ms 长时程确认运行通过。旧代码(reflux 树)同算例作为健康对照。
- **精度**:RZ Euler MMS(`exm/mms/eulerrz`,轴正则制造解)修复版与全度量版 **N=64/128 逐位相同**(max|Δρ| = 0)——传感器在光滑流下零触发,含轴五阶完整保留(N=128 L2(ρ) = 3.96×10⁻¹¹;旧代码同 case 为二阶 4.19×10⁻⁷)。
- **已排除的替代方案**(均有数据否证):普通分裂形式 collar(喷口 step 8 崩,更早);α×5 附加耗散(过冲翻倍至 −177010);cfl 降至 0.1(仍崩 → 正性问题而非 CFL 问题);静态宽 collar + 加权耗散(抽质量致更早崩)。

## 6. 向其他 RZ 格式(AFD / HLLC / skew-symmetric)移植的原则

MMS 扫掠显示 RZ 下 AFD 与 skew-symmetric 存在同类轴线发散,任何在轴邻区重构/差分度量变量 `r·(·)` 的格式都可能需要同类处理。

**skew-symmetric 的直接实证 (2026-08-07, b2d_skew4_newsrc_L4)**: 新树 Skew.h 的 G=rF 移植版(4 阶 + JST)在同一 2D RZ 射流算例上**不崩溃**(JST 耗散维持正性, 完整跑完 0→6 ms),但产生**网格附着的轴线密度丝**(axis filament): 轴邻 1–2 个单元的时间平均密度高出离轴正常值约 +70–120%(2.2–2.6 ρ_j vs 离轴 1.0–1.4 ρ_j),丝宽随当地最细网格尺度缩放(L4 区 1 格 / L3 区 1 个 L3 单元);同一单元的平均径向速度 ≈ +57–64 m/s(轴正则性要求 u_r = O(r)),3–5 格内衰减到 ~0。**离轴 ≥3–8 格的射流本体与旧代码及实验一致**。即该缺陷在耗散型格式下以亚临界形式存在——不失稳但污染轴上解,进一步证明轴邻径向面处理必须修,而非仅 WENO 特有。证据图: analysis/panda_benchmark/skew4newsrc_axis_filament.png。

移植时不变的三条原则:

1. **局域化**:回退仅装在最内侧内部径向面(K=1),由间断传感器门控;全局或静态降阶会牺牲设计阶数(本例五阶→二阶),且没有必要——失稳可证明局域于该面。
2. **通量约定一致**:回退通量必须严格遵循宿主格式的面数组存储约定(本例:非轴面存 `h_(rF)/r_face`、轴面存原始 h、paired-pressure 通道分离压力)。覆写点要选在高阶通量写入之后、面循环退出之前,保证与散度装配读到的是同一套约定。
3. **耗散作用于物理跳**:任何回退/附加耗散项一律写 `α·(U_R − U_L)`(或宿主 Riemann 求解器的等价物理跳),不得对状态量做度量加权——否则均匀态不守恒。对自带耗散的近似 Riemann 求解器(HLLC),等价做法是把该面的重构输入降为一阶单元值,而对流部分的度量加权平均保持不变;对 skew-symmetric,注意其分裂形式的两半都要保持同一约定。

## 7. 相关资产位置

| 资产 | 位置 |
|---|---|
| 补丁(ed-diff)与验证版头文件 | `patches/rz_axis_fix_20260807/` |
| RZ Euler MMS 算例 | `exm/mms/eulerrz/`(制造解生成器脚本另存,阶梯 `fixed_dt ∝ h^{5/3}`) |
| 度量路径停用对照树 | AWS `/shared/cerisse_metrictest` |
| 静态 collar 实验树 | AWS `/shared/cerisse_collartest`(K 常量同名) |
| 新旧求解器 A/B 复现 | AWS `/shared/solver_ab/{run_old,run_new}` |
| 失稳单元诊断(argmin 输出) | `src/tim/advance.cpp` 状态检查(含 `ParallelDescriptor::ReduceRealMin` 全局归约——注意:秩局部分支内加集合调用会死锁) |

**当前状态提示**:归档补丁已合入当前 `Weno.h`，但没有直接
覆盖现文件。实现已重整到共享的特征 LLF 面通量驱动中，因此 WENO
和 TENO 共用同一逻辑。默认通用门控要求相邻密度跳变和密度二阶差分
同时超过阈值。第一内部径向面保留已验证的单密度跳变门控。轴面
`r=0` 不使用此回退。运行时参数和方法清单见 `docs/input.md`。
