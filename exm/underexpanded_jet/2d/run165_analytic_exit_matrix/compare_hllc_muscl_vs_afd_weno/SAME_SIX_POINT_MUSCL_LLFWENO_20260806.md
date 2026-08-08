# 同一六点状态：MUSCL–HLLC、LLF–WENO 与 AFD–HLLC–WENO 对照

## 对照条件

- 同一 `RK2 stage-1 FE` 状态；
- 同一 `face=(4,171), dir=1`；
- 同一 `j168...j173` 六点数据；
- `dt=2e-7 s`, `dx=1.0612967e-3 m`；
- 不包含 IBM、BC 或不同格式此前的演化历史。

## MUSCL–HLLC

MUSCL 实际只使用：

```text
left slope/face state : j169,j170,j171
right slope/face state: j170,j171,j172
```

`j168,j173` 不参与该面的二阶 MUSCL 重构。MC limiter 得到：

```text
left : rho=0.381180628, p=7311.81283 Pa,
       un=-796.489209 m/s, ur=155.783982 m/s
right: rho=0.338058172, p=5073.72100 Pa,
       un=-925.018993 m/s, ur=158.614492 m/s
```

HLLC 波速：

```text
SL=-1066.71316 m/s, SR=-604.867481 m/s
```

`SR<0`，所以直接返回右物理通量：

```text
F_rho = -312.710230
F_mn  =  294336.623
F_mt  =  -49600.3744
F_E   = -1.54147017e8 W/m2
face CFL = 0.20102
```

## characteristic LLF–WENO-Z5

当前 shock fallback 的 LLF–WENO 使用这六点计算：

```text
alpha=max_j(|un_j|+c_j)=1048.05257 m/s
face CFL=dt*alpha/dx=0.197504
```

它以相邻单元的 Roe 状态冻结一个统一特征基，然后对每点形成：

```text
f_plus  = 0.5*(U + F/alpha)
f_minus = 0.5*(U - F/alpha)
```

正、负分支分别投影到同一冻结特征基，使用五阶 WENO-Z 重构，最后变回守恒通量。
它不构造左右 primitive face state，因此不存在可单独报告的
`p_face,u_face`，直接输出共享面数值通量。把指定六点直接送入当前算子得到：

```text
q=1:
F_rho = -316.847821
F_mn  =  273384.540
F_mt  =  -49264.3179
F_E   = -1.35883246e8 W/m2

q=2:
F_rho = -334.907356
F_mn  =  281019.275
F_mt  =  -52479.5503
F_E   = -1.37403833e8 W/m2
```

这里的 `q=2` 是用户指定的平方 smoothness-ratio 版本。原先记录的 C++
checkpoint face trace 在第一 RK 阶段已经切换到 LLF，因而它的第二阶段六点状态
不再严格等于本文指定六点；该 trace 不能冒充“固定六点”对照，已从下表移除。

## 三种通量直接比较

| 方法 | F_rho | F_mn | F_mt | F_E [W/m2] |
|---|---:|---:|---:|---:|
|MUSCL–HLLC|-312.710|294336.6|-49600.4|-1.54147e8|
|LLF–WENO, q=1|-316.848|273384.5|-49264.3|-1.35883e8|
|LLF–WENO, q=2|-334.907|281019.3|-52479.6|-1.37404e8|
|坏 AFD–HLLC–WENO|-2724.926|24145620.8|-428286.4|-1.07012e11|

MUSCL–HLLC 与两个 LLF–WENO 结果处于同一物理量级；坏 AFD 点状态 WENO 的
能量通量则约为 MUSCL 的 `694` 倍、LLF–WENO 的 `779--788` 倍。

## 仅替换这一面的更新反事实

保持原失败计算的其他面和 RZ RHS 不变，只把这一面换成对应通量：

|替换通量|第二 FE bracket rho|第二 FE bracket p [Pa]|最终 RK2 rho|最终 RK2 p [Pa]|
|---|---:|---:|---:|---:|
|MUSCL–HLLC|0.277526|4623.66|0.301284|4526.30|
|LLF–WENO, q=1|0.276746|4488.54|0.300894|4435.41|
|LLF–WENO, q=2|0.273343|4369.81|0.299193|4369.86|
|坏 AFD–HLLC–WENO|-0.177053|不可容许|0.073995|-1.45368e7|

因此只考虑该六点和该共享面，MUSCL–HLLC 与 LLF–WENO 都不会导致这次非物理
更新；失败专属于当前逐点非线性 acoustic-state WENO 后进入单侧 HLLC 的路径。
