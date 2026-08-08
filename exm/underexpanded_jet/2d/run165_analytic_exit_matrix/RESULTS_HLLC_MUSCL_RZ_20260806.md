# Run165 二阶 HLLC--MUSCL RZ 启动测试（2026-08-06）

## 结论

在与 Skew--JST 失败测试相同的二维 RZ、均匀 L1、Run165 外流和
100 microsecond 喷流 ramp 条件下，二阶 HLLC--MUSCL 在 `CFL=0.2`
下从干净初场推进到 500 microseconds，未出现负密度、负内能、NaN
或 Inf。这个结果越过了 Skew--JST 在约 13.37 microseconds 的重复失败点，
并覆盖了喷流完全开启后的约 400 microseconds。

这是一项稳定性和机制门禁，不是稳态 SRP 验证。该算例没有 IBM、
没有实体 forebody/nozzle，只在轴向高端面施加 Run165 虚拟喷口出口状态；
500 microseconds 时压缩结构仍在向上游移动。

## 离散与配置

- 方程：Euler。
- 坐标：二维轴对称 RZ。
- 对流通量：`riemann_t<false, ProbClosures>`，characteristic MUSCL 线性重构、
  有限斜率 limiter、HLLC Riemann 通量。
- 时间推进：沿用 `inputs` 的 `cns.order_rk=-2`、`cns.stages_rk=2`。
- AMR：L0 为 32 x 96；L1 覆盖全域，为 64 x 192，12,288 个 L1 cells；
  `dr=dx=1.06129666667 mm`。
- 喷口出口直径：L1 上 12 cells。
- 喷流：Run165 出口状态，100 microseconds 线性 ramp，质量流量目标
  `0.2822379563 kg/s`。
- IBM：关闭（`USE_GPIBM=FALSE`）；没有修改或调用 GP-IBM 生产路径。

## 运行结果

| 主机 | CFL | 终止时间 | coarse steps | advance wall time | 平均时间/step | 结果 |
|---|---:|---:|---:|---:|---:|---|
| 本地 CPU，1 OMP thread | 0.2 | 120 microseconds | 272 | 3.797 s | 0.01396 s | 通过 |
| Y9000X，1 OMP thread | 0.2 | 120 microseconds | 272 | 2.299 s | 0.00845 s | 通过 |
| Y9000X，1 OMP thread | 0.2 | 500 microseconds | 1178 | 10.631 s | 0.00902 s | 通过 |

Y9000X 的 500 microseconds profiler 中，Euler face flux 占总时间约 66.2%。

120 microseconds 场：`Mach=[0.0956,7.0403]`，
`p=[533.70,120570.05] Pa`，轴线上首次喷流 `M=1` 位置约 `x=-44.1 mm`。

500 microseconds 场：`Mach=[0.00255,8.9381]`，
`p=[204.49,120588.55] Pa`，相应轴向位置约 `x=-106.2 mm`。
所有报告的压力均为正。

## 与 Skew--JST 的直接比较

相同几何、网格、边界状态、100 microseconds ramp 和 `CFL=0.2` 下，
Skew--JST 在 `t=13.37030308 microseconds` 的 `(r,x)=(1.591945,-2.334566) mm`
产生负内能；HLLC--MUSCL 推进到 500 microseconds 仍保持 admissible state。
因此当前失败不是 RZ、边界状态或网格单独决定的，离散通量/重构的鲁棒性
是主导因素之一。

## 限制

HLLC--MUSCL 本身没有实现严格的数学正性证明或逐 stage 正性 limiter；本次
通过是经验稳定性证据。二阶 limiter 带来的耗散也会扩宽激波和剪切层，必须
用更细网格做激波位置、压力和质量流量收敛，才能判断这种稳定性是否以过度
耗散为代价。下一阶段应在完整 SRP/IBM 构型中复测，不能用本算例替代实验对比。

## 构建命令

```bash
make -j8 EULER_SCHEME=hllc-muscl JET_ENABLED=1 NS_ENABLED=0 LES_ENABLED=0 \
  USE_MPI=FALSE USE_OMP=TRUE \
  EBASE=run165_analytic_exit_hllc_muscl_o2_rz_omp \
  TMP_BUILD_DIR=tmp_build_run165_hllc_muscl_o2_rz_omp
```

