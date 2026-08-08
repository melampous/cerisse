# Run 165 三维近壁低温故障修复审计

日期：2026-08-08

## 结论

当前 Run 165 任务在 Step 3367 因近壁单元跌破 10 K 而终止，已经找到并修复一个可确定复现的构建/AMR 分派缺陷：临时 H100 构建文件漏掉了
`CNS_IBM_PROBLEM_SUPPORT_CAP`，而旧 `CNS::errorEst()` 又只在该宏存在时才调用算例已经定义的
`ibm_support_refine_allowed()`。因此，算例要求的“完整固体界面只到 L1、L2 只加密纯流体喷流核心”被静默忽略，整个连接体表面被错误推到 L2。

修复后，生产源码通过 ADL/SFINAE 自动检测算例 hook，不再依赖复制到构建文件中的宏；没有 hook 的算例保持原来的整界面加密行为。使用最终源码、CUDA 12.6、H100，从同一 `chk03300` 严格重放并推进到 Step 3400：

- `cns.strict_positivity=1`；
- `cns.soft_positivity=0`；
- 正性通量 limiter 关闭；
- 退出码 0，正性失败数 0；
- Step 3400 全场最低温度为 21.9379 K，`T<20 K` 单元数为 0；
- 原故障邻域最近的 L1 流体单元温度为 104.068 K，而旧软下限延续在该邻域仍约为 13.383 K。

这证明**当前生产算例的 10 K 终止已经消除**，也证明该次 regrid 失败不是 coarse-to-fine 插值直接制造坏状态。但这不等价于证明“任意 L2 贴壁 GP 网格都稳定”：旧的、非设计层级上的 L2 贴壁解为何产生随时间增强的法向穿壁速度，仍应作为独立的分辨率一致性问题审计。

## 确定性因果链

1. 算例在 `IBM/cases/test1853_single_3d/prob.h` 中定义了 `ibm_support_refine_allowed()`：`parent_level==0` 时允许界面种子加密，其他层禁止。因此固体界面目标层级是 L1，L2 仅用于纯流体 jet core。
2. 正式 CPU/A100/H100 makefile 都显式添加了 `-DCNS_IBM_PROBLEM_SUPPORT_CAP`；故障任务使用的临时 makefile 漏掉该宏。
3. 旧源码把 hook 调用放在同一宏条件中。漏宏不会编译失败，也没有运行时告警，而是退回“整个界面继续加密”。
4. 故障层级的 L2 有 22,436,864 个有效单元、1,855 个 grids；L2 GP 目标为 739,216，L2 surface 目标为 706,084。原故障点也位于 L2、距 STL 约 0.349 mm。
5. 最终源码从同一 `chk03300` 重建层级后，L2 降为 372,736 个纯流体单元、56 个 grids；L2 GP 和 surface 目标都严格为 0，固体表面由 L1 承担。
6. regrid 诊断给出 L1 和 L2 的 `new_from_coarse=0`。新的 L2 core 单元全部与旧 L2 重叠，没有任何新有效单元由 coarse-to-fine prolongation 生成。因此本次修复不是靠插值重置流场，也排除了“该次 regrid 直接插出负内能”的解释。

## 源码修复

修改位于 `src/CNS.cpp`：

- 为 GP-IBM 构建增加 `HasIBMRefineSupportHook<ProbParm>` 检测；
- 用无资格名调用和 ADL 查找与 `PROB::ProbParm` 同命名空间的算例 hook；
- hook 存在时自动调用，不存在时返回 `true`，维持整界面默认行为；
- 启动时输出
  `[IBM-AMR-Support-Hook] detected=<0|1> dispatch=automatic default=whole_interface`，使构建结果可审计；
- 删除实际分派对 `CNS_IBM_PROBLEM_SUPPORT_CAP` 的依赖。

最终 `src/CNS.cpp` SHA-256：

```text
5c23101ec5ab55f175c827f67b883dd3214831ad7728987820c1785b24d1cd3d
```

H100 最终可执行文件 SHA-256：

```text
f433b965e52436654e973fc69eb5bb07290512743e9cb6ed61fc52d825dceec1
```

该 H100 构建**故意不定义**旧宏，运行时仍报告 `detected=1`，从而直接验证修复不是被遗留宏“假通过”。

## 同检查点对比

| 状态 | Step | 全场 Tmin (K) | 原故障邻域 T (K) | 邻域速度 (m/s) | 估计壁法向速度 (m/s) |
|---|---:|---:|---:|---:|---:|
| 旧错误层级 | 3310 | 约 21.9 | 65.707 | 427.46 | -147.44 |
| 修复层级 | 3310 | 21.923 | 94.788 | 362.72 | -63.44 |
| 旧错误层级 | 3360 | 10.669 | 15.535 | 509.89 | -219.49 |
| 修复层级 | 3360 | 21.845 | 97.913 | 386.74 | -80.25 |
| 修复层级 | 3370 | 21.862 | 95.423 | 392.19 | -85.28 |
| 旧软下限延续 | 3400 | 10.000 | 13.383 | 536.18 | -247.22 |
| 最终源码严格延续 | 3400 | 21.938 | 104.068 | 389.62 | -73.54 |

修复后全场最低温度位于自由喷流：`(x,y,z)=(-54.1875,25.6364,-26.7273) mm`，总温约 348.013 K、速度约 809.52 m/s，属于空气喷流膨胀产生的低静温，而不是近壁总能异常。

表中“原故障邻域”在修复层级上是距旧坐标约 0.907 mm 的最近 L1 叶层流体单元；由于层级改变，它不是同一个离散单元。法向速度使用旧故障点最近 STL 三角面的法向，仅用于趋势诊断，不能替代 BI 上的无穿透误差。

## 验证结果

### H100 实际算例

- Test B：从原始 `chk03300` 强制 regrid，严格推进到 Step 3370，退出码 0，失败数 0。
- Test C：使用最终源码和最终 H100 可执行文件，从修复后的 `chk03370` 推进到 Step 3400，退出码 0，失败数 0。
- 两次运行均报告 `family=pure_shared_gp_full_cartesian`、`cut_control=disabled`。
- Test C 报告 `detected=1`，L2 GP 目标为 0。

### 自动分派负路径

强制重编译无算例 hook 的 2-D smooth-wall MMS 可执行文件，并运行 AMR 支撑正/负门：

- 运行时报告 `detected=0`，证明无 hook 时能正常编译和回退；
- 正例：L1 GP 与 surface 的 `coarse-fine-missing-supports=0`；
- 负例：统计到 1,248 个缺失支撑并按设计 fail-closed；
- 方法清单仍为 `pure_shared_gp_full_cartesian`，`cut_control_compiled=0`。

该最终 MMS 可执行文件 SHA-256：

```text
4e30949bcaf6b4aa30725deb9f9c1e8294775f5681c78fe9b38ff46b78affe2c
```

另已完成 3-D Run 165 CPU MPI 编译，以及无 hook 的 R--Z pure-GP 两步 SSPRK(4,3) 正性 smoke；后者所有 Forward-Euler bracket 和同步后审计均通过。

### Step 3400 后续推进与能量检查

Test D 从修复后的 `chk03400` 继续推进到 Step 4000，恢复算例正常的
`amr.regrid_int=4 4`，并保持：

- `cns.do_reflux=1`；
- `cns.strict_positivity=1`；
- `cns.soft_positivity=0`；
- `cns.ibm_positivity_flux_limiter=0`；
- 不设置最低温度，不通过状态裁剪注入能量。

Step 3400 到首个新输出 Step 3425 的 AMR 叶层流体积分为：

| 量 | Step 3400 | Step 3425 | 相对变化 |
|---|---:|---:|---:|
| 流体体积 (m^3) | 0.08173226677685913 | 0.08173226677685906 | -8.49e-16 |
| 质量 (kg) | 0.004518623126305413 | 0.004520182223937732 | +3.4504e-4 |
| 总能 (J) | 1016.9187425981411 | 1017.4798195152259 | +5.5174e-4 |
| 动能 (J) | 540.9661075023664 | 541.2502786128329 | +5.2530e-4 |
| 内能 (J) | 475.95263509577455 | 476.22954090239494 | +5.8179e-4 |
| Tmin (K) | 21.937913 | 21.843972 | -- |

这段物理时间为 `5.60438e-7 s`。新增总能除以新增质量约为
`359.9 kJ/kg`，与来流和喷流约 `340--347 kJ/kg` 的总焓量级一致；动能和内能同步平滑变化，未发现 regrid 瞬时能量跳跃。由于 Run 165 是开放域，这些积分会随来流、喷流和出口净通量变化；在没有逐物理边界数值面通量账本前，该结果只能排除明显的异常能量注入，不能宣称机器精度的真域守恒残差。

## Pure GP 范围确认

本次修改只影响 **GP geometry / AMR interface support tagging**：

- 没有修改共享 ghost-point 状态；
- 没有修改 BI-CWLS、热力学延拓、WENO/LLF 通量、完整背景单元 RHS、正性 limiter 或表面恢复；
- 没有引入 cut cell、cut control volume、area/volume fraction、aggregate、redistribution 或 AMReX EB 流动离散；
- 每个活跃流体单元仍使用完整背景网格体积和完整笛卡尔面；
- 运行时方法清单保持 `pure_shared_gp_full_cartesian`；
- 既有 AMR 支撑严格审计仍 fail-closed。

由于这次行为变化只发生在 `max_level>0` 且算例定义 hook 的 `errorEst()` 路径，统一 L0 的 Mach-4 circle、PM expansion、compression corner、smooth MMS 和 inverse-MOC 数值算子没有变化。本轮没有重新跑完这些耗时的统一网格长时基准；这部分是回归覆盖缺口，不能写成已重新认证。

## 尚未关闭的问题

1. 旧 L2 贴壁层级上，首层流体单元法向速度随时间增至约 `-226.6 m/s`。当前修复恢复了算例设计层级并消除生产故障，但没有解释该非设计 L2 解的全部离散机制。
2. 修复后邻近壁面的 L1 单元仍有约 `-73.5 m/s` 的投影速度。它是有限距离处的单元中心，不要求严格为 0，但在定量表面流动结论前仍需做 BI 速度和法向剖面审计。
3. 下一项高价值测试应是受控的小域 L1/L2/L3 同物理时刻矩阵，冻结 AMR、记录共享 GP 状态、每个穿壁背景面通量以及每个 SSPRK Forward-Euler bracket 的内能变化。只有该测试才能判断 finer-wall 异常是瞬态/网格采样，还是 GP 延拓与背景通量的分辨率不一致。
4. 在上述审计完成前，不应把“整车表面强制到 L2”作为可靠生产配置；当前可继续使用已恢复的 L1 表面 + L2 纯流体 jet core 配置。

## 证据路径

- H100 Test B 日志：`h100_results/testB_run.log`
- H100 Test B 清单：`h100_results/run_manifest.txt`
- H100 Test C 日志：`h100_results/testC_final_source/run.log`
- H100 Test C 清单：`h100_results/testC_final_source/run_manifest.txt`
- 最终构建清单：`h100_results/testC_final_source/support_hook_build_manifest.txt`
- Step 3400 能量审计：`h100_results/final_source_energy_3400.json`
- Step 3400/3425 能量对比：`h100_results/energy_3400_3425.json`
- 无 hook AMR 门：`amr_support_gate/final_source/positive.log`、`negative.log`
- H100 远程运行目录：`/shared/cerisse_test1853_3d_20260807/runs/test1853_single_h100/run165_support_hook_repair_20260808/`
