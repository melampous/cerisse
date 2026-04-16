# Moving-Body IBM 修复工作总结 (v5-baseline)

> **定稿版 — 2026-04-16**
>
> 本文件记录 2D 菱形翼型 Mach 2 Euler 算例下, moving-body IBM 模块两类
> 数值失稳的定位与修复过程. 以 `v5-frozen-debug` (commit `dc70f8a40`)
> 为证据快照, 以 `baseline/ibm-v5-clean` 为后续 FSI 开发基线.

---

## 1. 问题陈述

2D 菱形翼型 Mach 2 Euler 算例:
- 720 μs settle 后开始恒定 pitch rate, ω = 2330 rad/s
- Pitch 轴位于 1/4 弦长 (−0.025, 0)
- IBM + rigid motion + AMR 三者耦合

**原始症状**:
- 静态 IBM (`ib.move=0`) 完全正常
- Pitch 一启动, moving-body 很快失稳
- `max_level=0` 在约 step 1933 崩 (~1° pitch)
- `max_level=1` 在约 step 2107 崩 (~6-9° pitch)

## 2. 修复时间线

### 2.1 初始假设: AMR / post_regrid / FillPatch 污染 — **证伪**

在 `CNS::post_regrid` 加跨级 fresh-cell 修复. 过早触发, 破坏 settle 阶段,
回滚.

**关键转折**: `max_level=0` 隔离测试后, level=0 仍在 step 1933 崩. 首发
机制**不**由 AMR 独有, 而是 moving-IBM 本身的 bug.

### 2.2 v3: 修复 exposed-cell 首发失稳

`fixExposedCells` Pass 1 由 "4 邻居保守量平均" 改为 **方向性 score 单 donor**:

    score = -(u · di + v · dj)

结果:
- Level=0 可稳定跑满 5000 steps, pitch 达 ~153°, schlieren 主结构合理
- Level=1 仍在 step 2107 失稳

结论: exposed-cell averaging 是真实 bug, 但不是 AMR 情况下的唯一根因.

### 2.3 v5: 修复首个 AMR donor 污染机制

进一步加入:
- `fixExposedCells` Pass 1 donor 合法性检查 (p > 0, finite)
- Pass 2/3 同样改为 single-donor 逻辑
- `Weno.h` 对 `near_ib = (gl < ng || gr < ng)` 的 face 启用 1st-order LLF fallback

结果:
- Level=1 首发坏状态**不再在 step 2107 出现**
- first BAD 推迟到约 step 2397 (比原始版延后 ~290 steps)

结论: fixExposedCells donor 污染导致的第一类失稳机制已显著削弱, 但后续
仍有第二类结构性失稳.

### 2.4 v6: 扩大 IB-band + positivity fallback — **回归, 放弃**

尝试进一步加入:
- widened near-IB mask (radius 2)
- face-level positivity fallback

结果: first BAD 反而提前到 step 2325 (比 v5 早 72 步).

结论: heuristic patching 开始互相打架, widened mask 过度耗散改变尾迹
演化. **不再沿此路线堆补丁**.

## 3. 最终冻结版本: v5 baseline

**保存位置**:
- GitHub: `melampous/cerisse`
- Branch: `baseline/ibm-v5-clean`
- Tag: `v5-frozen-debug`

> **用途限定**: `baseline/ibm-v5-clean` 的目的, 是**冻结 exposed-cell 修复
> 与 IBM-near-wall fallback 这两类已验证有效的改动**, 作为 moving-FSI
> 调试基线, **而不是最终生产级 AMR 求解器**.

### 3.1 保留的算法修改

- `src/ibm/ibm_solver.h`:
  `fixExposedCells` Pass 1/2/3 全改 single-donor
  (directional score + donor p > 0 / finite 检查)
- `src/rhs/Weno.h`:
  `eflux_ibm` 中 stencil 触 IBM solid/ghost 时
  fallback 1st-order LLF

> **Scope note**: `Weno.h` 中的 near_ib 1st-order fallback 作用于**所有**
> IBM case (不仅 airfoil pitch-up), 不只是 moving body. 此 branch 因此
> **有意定位为 moving-FSI 调试专用**, 不作为通用 IBM baseline.

### 3.2 保留但默认关闭的调试基础设施

- `src/ibm/ibm_tripwire.h` — T1/T2/T3/T4 四阶段 pipeline 扫描器
- `advance.cpp` / `CNS.cpp` — tripwire 调用点 (窗口 = (-1,-1) 关闭)
- `fixExposedCells` DIAG — per-cell donor/fresh-cell 记录
  (窗口 = (-1,-1) 关闭)

改常量即可重新激活.

## 4. 当前能力与限制

| 项目 | 状态 |
|---|---|
| 静态 IBM (`ib.move=0`) | ✓ 稳定 (未变动) |
| moving body, level=0 | ✓ 稳定到 5000 steps / ~153° pitch |
| moving body, level=1 AMR | △ **首发失稳由 step 2107 推迟至 step 2397; 在 ~24° AoA 之前局部结果可参考; 完整长时演化仍不可视为 fully physical** |
| moving body, level ≥ 2 | ✗ 未测试 |
| moving body + viscous (NS) | ✗ 未测试 |

## 5. 第二类根因: 已定位到模块, 未修复到 face

Tripwire 显示: **step 2373 附近 first BAD 出现在 `T1_step_entry`,
T1/T2/T3 数值完全相同**, 说明:

> 坏状态诞生于**前一步的 advance 过程 (即 RHS/reconstruction/update 路径)**,
> 不是当前步的 `rebuildIBM` 或 `fixExposedCells`.

> **证据上限**: 现有证据足以将第二类失稳归因到
> `compute_rhs / reconstruction / update` 这条路径, **但尚不足以证明是单一
> 公式、单一 face, 或单一 coarse-fine 交换操作的唯一责任**.

现象与以下机制一致 (相关, 非唯一):
- Airfoil 在较大 AoA 下, 吸力面出现强扩张区
- IBM/ghost 邻域 + coarse-fine 邻域中的高阶重构出现 overshoot
- 生成负压 / 病态高动量 cell
- 坏状态经 AMR 同步/传播扩散到其他 level

## 6. 下一阶段技术路线 (按优先级)

### 6.1 IBM 邻域改用 primitive-space MUSCL + limiter  *(最推荐下一步)*
- 远离壁面: 保持 WENO-Z5
- IBM/ghost 邻域: 不用 characteristic WENO, 改用 primitive-space MUSCL +
  minmod/MC
- 避免 IBM 邻域的病态 characteristic / Roe 变换, 用一致的稳健低阶方案
  替代 heuristic fallback

### 6.2 Face-local positivity-preserving fallback
- 不是 RK 后 clip cell, 而是**每个 face 上检查 reconstructed L/R states**
- ρ_L/ρ_R/p_L/p_R 不物理 → 该 face 直接回退到 LLF/Rusanov
- 比事后 clipping 更合理

### 6.3 更系统的 IBM/ghost 一致性处理 (中长期)
- primitive reconstruction for ghost-adjacent states
- moving-wall / GP consistency
- coarse-fine + IBM 邻域协同策略

### 6.4 其他补充选项 (仅作 backup)
- **Riemann-ghost-fluid** (Sambasivan-Udaykumar 2009): 物理最严谨,
  工程量最大
- **IBM geometry rebuild 频率控制**: 累积转角阈值驱动, 不是每步
- **Cut-cell IBM** (Seo-Mittal 2011): 等于重写 IBM 模块, 不现实

## 7. 工程经验教训

1. **先做 `level=0` 隔离测试** — 整个调试的转折点, 把 AMR 从嫌疑人名单移出
2. **Tripwire 比 log-inspection 值钱** — 2-3 个 stage 的定点诊断立刻锁定
   "compute_rhs 产生坏状态"
3. **Heuristic patch 边际收益递减很快** — v3→v5 大提升, v5→v6 回归.
   该停就停
4. **必须诚实保留每个版本的收益与回归** — v3 bought 290 steps, v5 bought 24
   more, v6 lost 72. 防止后人误以为 v6 是 v5 的改进
5. **外部 peer review 很值** — 多次纠正概念混淆 (ghost-cell vs fresh-cell
   trace) 和过度 patching

## 8. "目前解是否物理" 的最终判断

### Level 0
**物理上可接受, 工程上可用**.

适合:
- 可视化
- 趋势分析
- 方法调试 baseline

不宜直接宣称:
- 已严格收敛
- 已达到发表级验证完毕

### Level 1
**相比原始版大幅改善, 但长时 AMR 解仍不 fully reliable**.

适合:
- 证明 first-failure donor 污染机制已修复
- 指导下一阶段方法设计

**不适合**:
- 作为最终物理解
- 做定量结论
- 做正式生产结果

---

## Appendix A. Commit / Tag / Branch 对应

| 标记 | commit | 内容 |
|---|---|---|
| `v5-frozen-debug` (tag) | `dc70f8a40` | 完整 debug 快照 (含 tripwire + DIAG 启用) |
| `baseline/ibm-v5-clean` HEAD | `d0034f43f` | debug 关闭, 算法修改保留, FSI 专用 |
| `ibm_fsi` HEAD | `ec5351541` | 不含本次 v5 修改; 原版 Weno |

## Appendix B. 测试数据保留位置

- `temp/airfoil_pitchup/level0_step1patch_run.log` — level=0 v5 5000 步
- `temp/airfoil_pitchup/level1_*_run.log` — level=1 各版本 log
- `temp/airfoil_pitchup/step1_l0_plt*.png` — level=0 schlieren (clean)
- `temp/airfoil_pitchup/v3_l1_plt02200_schlieren.png` — level=1 pre-failure
- `temp/airfoil_pitchup/v3_l1_plt02400_schlieren.png` — level=1 post-failure
  (显示 silent cascade)
