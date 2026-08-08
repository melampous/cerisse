#!/usr/bin/env python3
"""Plot the first hidden q=2 LLF-WENO RK-stage violation and its stencil."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np
import yt

from analyze_first_stage_violation import parse_trace
from analyze_hllc_muscl_afd_weno_pointwise import (
    load_fields,
    patch_yt_cylindrical_readonly_edge,
)


ROOT = Path(__file__).resolve().parent
BASE = ROOT / "compare_muscl_vs_llfweno_q2"
ANALYSIS = BASE / "analysis"
FLOW = BASE / "llfweno" / "first_event_flow"
AUDIT = BASE / "llfweno" / "stage_audit"
TARGET = (56, 133)
STENCIL_JS = np.arange(131, 137)
GAMMA = 1.4
PLAIN_TICK = FuncFormatter(lambda value, _: f"{value:g}")


def configure_plotting() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["AR PL SungtiL GB", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "font.size": 10.5,
        "axes.titlesize": 11.5,
        "figure.titlesize": 14,
    })


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def exact_stencil_rows() -> list[dict]:
    rows = read_csv(ANALYSIS / "first_violation_stencil_all_stages.csv")
    result = []
    for row in rows:
        if row["axis"] != "axial" or int(row["j"]) not in STENCIL_JS:
            continue
        result.append({
            "label": row["label"],
            "j": int(row["j"]),
            "rho": float(row["rho"]),
            "pressure": float(row["p_raw_Pa"]),
            "ur": float(row["ur_m_s"]),
            "uz": float(row["uz_m_s"]),
            "rhoe": float(row["rhoe_raw_J_m3"]),
            "energy": float(row["energy_J_m3"]),
        })
    return result


def draw_stencil(ax, r_edges, z_edges, r_centres, z_centres,
                 labels: bool = True) -> None:
    i, bad_j = TARGET
    for j in STENCIL_JS:
        bad = j == bad_j
        rect = Rectangle(
            (z_edges[j], r_edges[i]),
            z_edges[j + 1] - z_edges[j],
            r_edges[i + 1] - r_edges[i],
            fill=False,
            edgecolor="red" if bad else "#ffe45e",
            linewidth=2.5 if bad else 1.5,
            zorder=8,
        )
        ax.add_patch(rect)
        if labels:
            ax.text(
                z_centres[j], r_centres[i], str(j),
                color="white", ha="center", va="center",
                fontsize=8.5, weight="bold", zorder=10,
                path_effects=[patheffects.withStroke(linewidth=2, foreground="black")],
            )
    # axial-high face of cell j=133 is the edge between j=133 and j=134.
    ax.plot(
        [z_edges[bad_j + 1], z_edges[bad_j + 1]],
        [r_edges[i], r_edges[i + 1]],
        color="#00ffff", linewidth=4.0, zorder=11,
    )
    ax.scatter(
        [z_centres[bad_j]], [r_centres[i]], marker="x",
        s=90, linewidth=2.4, color="red", zorder=12,
    )


def plot_flow_context() -> None:
    snapshots = {}
    for step in (541, 542):
        snapshots[step] = load_fields(FLOW / f"plt{step:05d}")

    fields_542, time_542, r_edges, z_edges, rc, zc = snapshots[542]
    speed = np.hypot(fields_542["ur"], fields_542["uz"])
    sound = np.sqrt(GAMMA * np.maximum(fields_542["pressure"], 1.0e-12)
                    / fields_542["rho"])
    mach = speed / sound

    fig, axes = plt.subplots(2, 2, figsize=(15.2, 10.4), constrained_layout=True)
    ax = axes[0, 0]
    pcm = ax.pcolormesh(
        z_edges, r_edges, fields_542["pressure"], shading="flat",
        cmap="turbo", norm=LogNorm(vmin=10.0, vmax=1.3e5),
    )
    fig.colorbar(pcm, ax=ax, label="压力 p [Pa]", format=PLAIN_TICK)
    ax.scatter([zc[TARGET[1]]], [rc[TARGET[0]]], marker="x", s=100,
               color="red", linewidth=2.5, zorder=8)
    ax.add_patch(Rectangle((-66.0, 54.5), 13.5, 11.0, fill=False,
                           edgecolor="white", linewidth=1.8, linestyle="--"))
    ax.annotate(
        "最早坏点\n(i,j)=(56,133)\n(r,z)=(59.963,-59.645) mm",
        xy=(zc[TARGET[1]], rc[TARGET[0]]), xytext=(-145, 49),
        color="white", fontsize=10,
        arrowprops=dict(arrowstyle="->", color="white", lw=1.5),
        bbox=dict(facecolor="black", alpha=0.65, edgecolor="white"),
    )
    ax.set(title=f"完整压力场，Heun 接受态 t={time_542:.1f} μs",
           xlabel="轴向 z [mm]", ylabel="径向 r [mm]")

    ax = axes[0, 1]
    pcm = ax.pcolormesh(z_edges, r_edges, mach, shading="flat",
                        cmap="viridis", vmin=0.0, vmax=5.0)
    fig.colorbar(pcm, ax=ax, label="Mach", format=PLAIN_TICK)
    ax.scatter([zc[TARGET[1]]], [rc[TARGET[0]]], marker="x", s=100,
               color="red", linewidth=2.5, zorder=8)
    ax.add_patch(Rectangle((-66.0, 54.5), 13.5, 11.0, fill=False,
                           edgecolor="white", linewidth=1.8, linestyle="--"))
    ax.set(title=f"完整 Mach 场，Heun 接受态 t={time_542:.1f} μs",
           xlabel="轴向 z [mm]", ylabel="径向 r [mm]")

    for ax, step in zip(axes[1], (541, 542)):
        fields, time_us, re, ze, rcc, zcc = snapshots[step]
        pcm = ax.pcolormesh(
            ze, re, fields["pressure"], shading="flat",
            cmap="turbo", norm=LogNorm(vmin=10.0, vmax=2.5e4),
        )
        fig.colorbar(pcm, ax=ax, label="压力 p [Pa]", format=PLAIN_TICK)
        draw_stencil(ax, re, ze, rcc, zcc)
        ax.set_xlim(-66.0, -52.5)
        ax.set_ylim(54.5, 65.5)
        ax.set_aspect("equal")
        ax.grid(color="white", alpha=0.17, linewidth=0.5)
        ax.set(title=f"局部压力场，Heun 接受态 t={time_us:.1f} μs",
               xlabel="轴向 z [mm]", ylabel="径向 r [mm]")
    axes[1, 0].text(
        -65.6, 65.0,
        "黄框：轴向六点 stencil j=131…136\n红框/红叉：目标单元 j=133\n青线：危险 axial-high 面 j=133.5",
        color="white", va="top",
        bbox=dict(facecolor="black", alpha=0.65, edgecolor="white"),
    )
    axes[1, 1].text(
        -65.6, 65.0,
        "该 plot 是完整 Heun 后状态：p=27.06 Pa\n同一步 FE 中间候选：p=-8.38 Pa，ρe=-20.95 J/m^3",
        color="white", va="top",
        bbox=dict(facecolor="black", alpha=0.72, edgecolor="red"),
    )
    fig.suptitle("Run165 RZ：第一次隐性负内能的位置与轴向危险 stencil")
    fig.savefig(ANALYSIS / "first_violation_flow_context_annotated.png", dpi=210)
    plt.close(fig)


def plot_stencil_profiles() -> None:
    rows = exact_stencil_rows()
    stages = [
        ("heun_stage1_input", "216.6 μs：阶段输入", "o", "#1f77b4"),
        ("heun_fe_stage1_candidate", "216.8 μs：FE 候选", "x", "#d62728"),
        ("heun_final", "216.8 μs：Heun 最终态", "s", "#2ca02c"),
    ]
    variables = [
        ("rho", "密度 ρ [kg/m^3]", None),
        ("pressure", "压力 p [Pa]", "symlog"),
        ("ur", "径向速度 u_r [m/s]", None),
        ("uz", "轴向速度 u_z [m/s]", None),
        ("rhoe", "内能密度 ρe [J/m^3]", "symlog"),
        ("energy", "总能量密度 E [J/m^3]", None),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(14.5, 12.3), sharex=True,
                             constrained_layout=True)
    for ax, (key, ylabel, scale) in zip(axes.flat, variables):
        for label, legend, marker, color in stages:
            selected = sorted((r for r in rows if r["label"] == label),
                              key=lambda r: r["j"])
            ax.plot([r["j"] for r in selected], [r[key] for r in selected],
                    marker=marker, markersize=6.5, linewidth=1.8,
                    color=color, label=legend)
        ax.axvspan(132.72, 133.28, color="red", alpha=0.09)
        ax.axvline(133.5, color="#00a6a6", linestyle="--", linewidth=1.8)
        if key in ("pressure", "rhoe"):
            ax.axhline(0.0, color="black", linewidth=1.0)
        if scale == "symlog":
            ax.set_yscale("symlog", linthresh=100.0)
            ax.yaxis.set_major_formatter(PLAIN_TICK)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
        ax.set_xticks(STENCIL_JS)
    for ax in axes[-1]:
        ax.set_xlabel("轴向单元索引 j（左：更负 z；右：喷口方向）")
    axes[0, 0].legend(loc="best", fontsize=9)
    axes[0, 1].annotate(
        "首个负压\np=-8.38 Pa",
        xy=(133, -8.3792), xytext=(131.1, -500.0),
        arrowprops=dict(arrowstyle="->", color="red"), color="red",
    )
    axes[2, 0].annotate(
        "首个负内能\nρe=-20.95 J/m^3",
        xy=(133, -20.948), xytext=(131.0, -900.0),
        arrowprops=dict(arrowstyle="->", color="red"), color="red",
    )
    fig.suptitle(
        "轴向六点 stencil 的物理量：坏点前后空间点与 FE/Heun 状态\n"
        "红带为坏点 j=133，青虚线为危险轴向高侧面 j=133.5"
    )
    fig.savefig(ANALYSIS / "first_violation_six_point_physical_profiles.png", dpi=210)
    plt.close(fig)


def plot_bad_cell_history() -> None:
    rows = read_csv(ANALYSIS / "first_bad_cell_raw_history.csv")
    accepted = []
    for row in rows:
        t = float(row["time_us"])
        label = row["label"]
        if label == "heun_stage1_input" and abs(t - 216.0) < 1.0e-6:
            accepted.append(row)
        elif label == "heun_final":
            accepted.append(row)
    fe = [row for row in rows if row["label"] == "heun_fe_stage1_candidate"]

    def value(row, key):
        return float(row[key])

    variables = [
        ("rho", "密度 ρ [kg/m^3]"),
        ("p_raw_Pa", "压力 p [Pa]"),
        ("ur_m_s", "径向速度 u_r [m/s]"),
        ("uz_m_s", "轴向速度 u_z [m/s]"),
        ("rhoe_raw_J_m3", "内能密度 ρe [J/m^3]"),
        ("energy_J_m3", "总能量密度 E [J/m^3]"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(14.0, 11.8), sharex=True,
                             constrained_layout=True)
    for ax, (key, ylabel) in zip(axes.flat, variables):
        ax.plot([value(r, "time_us") for r in accepted],
                [value(r, key) for r in accepted], "o-", color="#1f77b4",
                linewidth=2.0, label="完整 Heun 接受态")
        ax.plot([value(r, "time_us") for r in fe],
                [value(r, key) for r in fe], "x--", color="#d62728",
                markersize=7, linewidth=1.5, label="第一阶段 FE 候选")
        ax.axvline(216.8, color="red", alpha=0.35, linewidth=1.0)
        if key in ("p_raw_Pa", "rhoe_raw_J_m3"):
            ax.axhline(0.0, color="black", linewidth=1.0)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("时间 [μs]")
    axes[0, 0].legend(loc="best", fontsize=9)
    axes[0, 1].annotate(
        "FE 首次越界\n-8.38 Pa",
        xy=(216.8, -8.3792), xytext=(216.15, 70),
        arrowprops=dict(arrowstyle="->", color="red"), color="red",
    )
    axes[2, 0].annotate(
        "FE 首次越界\n-20.95 J/m^3",
        xy=(216.8, -20.948), xytext=(216.05, 170),
        arrowprops=dict(arrowstyle="->", color="red"), color="red",
    )
    fig.suptitle(
        "最早坏点 (i,j)=(56,133) 的时间演化\n"
        "完整步仍可恢复，但 FE bracket 已经离开正性集合"
    )
    fig.savefig(ANALYSIS / "first_bad_cell_all_variables_time_history.png", dpi=210)
    plt.close(fig)


def plot_stencil_time_map() -> None:
    trace_specs = [
        (1080, AUDIT / "trace_1080_cell56_133.log"),
        (1081, AUDIT / "trace_1081_cell56_133.log"),
        (1082, AUDIT / "trace_1082_cell56_133.log"),
        (1083, AUDIT / "first_event_trace.log"),
    ]
    accepted_by_time = {}
    fe_bad = None
    for _, path in trace_specs:
        rows, _ = parse_trace(path)
        for row in rows:
            if row["axis"] != "axial" or row["j"] not in STENCIL_JS:
                continue
            t = round(row["time_us"], 7)
            if row["label"] == "heun_stage1_input" and abs(t - 216.0) < 1.0e-5:
                accepted_by_time.setdefault(t, {})[row["j"]] = row
            elif row["label"] == "heun_final":
                accepted_by_time.setdefault(t, {})[row["j"]] = row
            elif row["label"] == "heun_fe_stage1_candidate" and abs(t - 216.8) < 1.0e-5:
                fe_bad = fe_bad or {}
                fe_bad[row["j"]] = row
    columns = [(t, f"{t:.1f}\nHeun") for t in sorted(accepted_by_time)]
    columns.insert(-1, ("fe", "216.8\nFE"))
    data_sources = []
    for key, _ in columns:
        data_sources.append(fe_bad if key == "fe" else accepted_by_time[key])

    specs = [
        ("p_raw_Pa", "压力 p [Pa]", SymLogNorm(linthresh=100, vmin=-10, vmax=2.3e4), "coolwarm"),
        ("rhoe_raw_J_m3", "内能密度 ρe [J/m^3]", SymLogNorm(linthresh=250, vmin=-25, vmax=5.7e4), "coolwarm"),
        ("rho", "密度 ρ [kg/m^3]", Normalize(vmin=0.02, vmax=0.17), "viridis"),
        ("uz_m_s", "轴向速度 u_z [m/s]", Normalize(vmin=100, vmax=760), "plasma"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.8), constrained_layout=True)
    for ax, (key, title, norm, cmap) in zip(axes.flat, specs):
        array = np.array([
            [float(source[j][key]) for source in data_sources]
            for j in STENCIL_JS
        ])
        im = ax.imshow(array, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        fig.colorbar(im, ax=ax, label=title, format=PLAIN_TICK)
        ax.set_xticks(range(len(columns)), [label for _, label in columns])
        ax.set_yticks(range(len(STENCIL_JS)), STENCIL_JS)
        ax.set_xlabel("时间 [μs] / 状态")
        ax.set_ylabel("轴向 stencil 单元 j")
        ax.set_title(title)
        for yi, j in enumerate(STENCIL_JS):
            for xi in range(len(columns)):
                v = array[yi, xi]
                text = f"{v:.3g}"
                color = "white" if xi != len(columns) - 1 or j != 133 else "black"
                ax.text(xi, yi, text, ha="center", va="center",
                        fontsize=7.2, color=color)
        # Exact first negative candidate is column immediately before final Heun.
        fe_col = next(i for i, (key_col, _) in enumerate(columns) if key_col == "fe")
        bad_row = int(np.where(STENCIL_JS == TARGET[1])[0][0])
        ax.add_patch(Rectangle((fe_col - 0.49, bad_row - 0.49), 0.98, 0.98,
                               fill=False, edgecolor="#00ffff", linewidth=3.0))
    fig.suptitle(
        "轴向六点 stencil 的时空演化（数值为原始物理量）\n"
        "青框：216.8 μs FE 候选中首先变成非物理的 j=133"
    )
    fig.savefig(ANALYSIS / "first_violation_stencil_spacetime_values.png", dpi=210)
    plt.close(fig)


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    configure_plotting()
    patch_yt_cylindrical_readonly_edge()
    yt.funcs.mylog.setLevel(50)
    plot_flow_context()
    plot_stencil_profiles()
    plot_bad_cell_history()
    plot_stencil_time_map()
    for name in (
        "first_violation_flow_context_annotated.png",
        "first_violation_six_point_physical_profiles.png",
        "first_bad_cell_all_variables_time_history.png",
        "first_violation_stencil_spacetime_values.png",
    ):
        print(ANALYSIS / name)


if __name__ == "__main__":
    main()
