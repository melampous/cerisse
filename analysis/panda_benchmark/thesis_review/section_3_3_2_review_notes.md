# §3.3.2 "Nozzle exit profiles" — 修订说明(2026-07-22)

替换文件:`section_3_3_2_revised.tex`。该版本已于 2026-07-22 合并到
`/home/qiaoj/Documents/thesis-overleaf/03_free_underexpanded_jet/sections/04_cases_and_diagnostics.tex`
的 §3.3.2。表格连续剖面数值均经独立复算；top-hat 的 \(y_{50}\) 从
`0` 改为 `--`，因为不连续剖面不存在严格的半速度点。

## 必须修(已落实)
1. **求积声明改为可核验的版本**(原 302–304 行):原句宣称"生产网格 99 µm 上质量与动量通量 0.01% 内"。现明确区分 \(D_e/512=49.6\) µm 求积检查、\(D_e/256=99.2\) µm 四剖面对比网格和 \(D_e/8\) 粗网格；只对孔面积与 \(C_d\) 给出离散误差，并注明 \(C_M\) 为解析值。
2. **删除内部注记残留**"not linearly blended with the ambient state"(原 203–204 行)，改为正面陈述温度和密度的径向变化由 \(f(r)\) 决定，而静压保持均匀；同时修正"only the axial velocity is prescribed"的不严谨(T、ρ 亦随 \(f\) 变化)。
3. **复合修饰语连字符**:finite-thickness / top-hat / shifted-tanh / wall-tanh(修饰语位)/ wall-attached / half-velocity / matched-pair / shock-cell / cell-centred / one-seventh-power。表格行名(非修饰语位)保持原样。⚠️ 全章尚有 ~28 处 "three dimensional" 等待统一,本文件只覆盖本节。
4. **wall-tanh 断言补论据**:原"更像贴壁出口边界层"为印象式；现引用表内 \(H=2.26\) 接近层流平板值 2.59（shifted 族为 6.02），并准确说明 shifted profile 的半速度位置位于唇口内 \(3\delta\)。
5. **prescribed vs achieved \(\delta_\omega\)**:配对句加 "prescribed"，表题注明解析剖面量，并明确演化后的剪切层不假定保持该名义涡量厚度。没有保留原草稿中“结果节已经报告实测厚度”的错误前向承诺。

## 压缩(净省约 18 行)
- top-hat 6 行分段公式 → 一句内联(−7 行);
- 上游边界三区描述与 §3.3.3 重复 → 一句交叉引用(−3 行);
- C_d,eff 引入句与免责声明合并(−4 行);
- θ₀≃δ/2 与 0.31a 推导句 → 括号并入定义段(−3 行);
- δ、速度分量两个 display 块 → 内联(−4 行);
- 正文与图题重复的 Pohlhausen–Blasius 说明已合并；Figure 3.3 现明确前两幅含六条剖面，梯度面板只含四条有限且非奇异的梯度曲线，并将额外两条剖面指向独立的 profile-shape 比较。
- 未引用的公式标签补引用:Eq.(wall-tanh) 与 Eq.(measures) 现于正文引用。

## 集成与编译检查
- 论文源中 §3.3.2 现为 142 行（原 153 行），删除了两个独立 display 和重复的边界描述。
- XeLaTeX 全文编译成功，共 87 页；本节位于第 37–40 页。
- Table 3.3 与 Figure 3.3 均位于页边距内，无 overfull box、未定义引用或未定义 citation。
- 实际合并版本与 `section_3_3_2_revised.tex` 已同步。

## 尚未落实、需要你决定的两项(不在本文件内)
- **敏感性问题的结果回收**(必须修):§3.3.2 提出的四剖面敏感性问题需在 10_results_profile.tex 回收。现有 \(S_1\) 与 \(S_1/\sqrt{C_d}\) 数据可用于讨论与积分通量亏损的一致相关性，但不能写成 \(C_d\) 的因果控制；候选 realised \(\delta_\omega\) 数值在正式写入前还需与结果节的统一提取方法核对（数据位于 analysis/panda_benchmark/bc_variants/）。
- **两套剖面战役的区分**(必须修):设置矩阵(top-hat/wall-tanh/δ100/δ254)与 10 节的另一套(δ=50/100/200+Blasius+1/7)需在 case 矩阵表中分开列行并命名区分。

## v2(2026-07-22,按用户审查意见收紧结论范围)
1. **定位改写**:开篇明确 "treated as a sensitivity study rather than as a strict separation of the underlying mechanisms";判别句改为 "would be consistent with / would instead support",并加 "Neither observation would provide a unique causal separation, because the complete profile shapes remain different"。全文不再出现 "indicate control by";不讨论 KH 增长率。
2. **分辨率限定**:新增 δω/Δx≈2、δ99≈5 格的明示,"local velocity gradients and their instability development are not considered quantitatively resolved",结论限定为"离散施加剖面下平均激波胞系统的敏感性"。
3. **δ*/θ0/H 加限定**:velocity-based、locally planar、无密度权与圆柱面积权;删除 H=2.26 vs Blasius 2.59 的论证句与 wall-tanh 边界层类比(按意见 #3)。
4. **C_M 更名**:convective momentum-flux coefficient,注明不含出口压力推力与黏性应力;四例同出口压力故适用于比较剖面动量亏损。表题同步注明。
5. **补正式定义并保持排版**:\(f(R_e-\delta_{99})=0.99\) 纳入量度公式，\(\Delta U=u_e-u_\infty=u_e\) 紧邻公式在正文中定义。整组量度采用居中的 (3.15a)–(3.15e) 子公式，避免将三个长定义挤在同一行。
6. **审计改到实际网格 D/256**:本地按与求解器一致的胞心求积重算(约定用日志 D/8=1.03450713、D/512=1.000022264 及各 Cd 逐位校验通过):D/256 离散孔面积 −0.0075%;Cd 与离散 C_M:top-hat −0.0075%(受面积主导),其余三剖面 ≤0.001%;D/512 ≤0.003% 作细网格参考。正文按此改写;D/8 的 +3.45% 从本节移除(建议放网格审计附录)。审计脚本数据见本目录 audit_d256.txt。
7. **整合时保留的两项一致性修正**:top-hat 的 \(y_{50}\) 继续记为双连字符，因为不连续剖面不存在严格半速度位置；未加入“结果节已报告 realised vorticity thickness”这一尚未兑现的前向声明。Figure 3.3 继续使用 H 浮动位置，保证 3.3.3 紧接图后。
