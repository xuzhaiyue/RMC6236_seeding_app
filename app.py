from __future__ import annotations

from datetime import datetime
from io import StringIO

import pandas as pd
import streamlit as st


TIMEPOINTS = ["0 h", "2 h", "6 h", "12 h", "24 h", "72 h", "Day 7", "Day 14"]

DEFAULT_DENSITIES = {
    "快生型": [400_000, 400_000, 400_000, 350_000, 300_000, 100_000, 30_000, 10_000],
    "中快/中速型": [500_000, 500_000, 500_000, 450_000, 400_000, 150_000, 50_000, 15_000],
    "慢生/分化型": [650_000, 650_000, 650_000, 550_000, 500_000, 220_000, 75_000, 25_000],
}

RESEEDING_GOALS = ["2 天左右长满", "3 天左右长满", "4–5 天左右长满"]

RESEEDING_DEFAULTS = {
    "快生型": {
        "2 天左右长满": {"T25": (400_000, 600_000), "T75": (1_200_000, 1_800_000)},
        "3 天左右长满": {"T25": (200_000, 300_000), "T75": (600_000, 900_000)},
        "4–5 天左右长满": {"T25": (100_000, 150_000), "T75": (300_000, 500_000)},
    },
    "中快/中速型": {
        "2 天左右长满": {"T25": (600_000, 800_000), "T75": (1_800_000, 2_400_000)},
        "3 天左右长满": {"T25": (300_000, 500_000), "T75": (900_000, 1_500_000)},
        "4–5 天左右长满": {"T25": (150_000, 250_000), "T75": (500_000, 800_000)},
    },
    "慢生/分化型": {
        "2 天左右长满": {"T25": (800_000, 1_000_000), "T75": (2_400_000, 3_000_000)},
        "3 天左右长满": {"T25": (500_000, 700_000), "T75": (1_500_000, 2_100_000)},
        "4–5 天左右长满": {"T25": (250_000, 400_000), "T75": (800_000, 1_200_000)},
    },
}

GOAL_DAYS = {"2 天左右长满": 2.0, "3 天左右长满": 3.0, "4–5 天左右长满": 4.5}
DEFAULT_DOUBLING_TIME_HOURS = {"快生型": 24.0, "中快/中速型": 30.0, "慢生/分化型": 40.0}
FULL_FLASK_CELLS = {"T25": 2_400_000, "T75": 7_200_000}


def build_default_table(preset: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "时间点": TIMEPOINTS,
            "cells_per_well": DEFAULT_DENSITIES[preset],
        }
    )


def calculate_plan(
    density_table: pd.DataFrame,
    cell_concentration_per_ml: float,
    final_volume_ml_per_well: float,
    wells_per_timepoint: int,
    extra_fraction: float,
    include_rmc6236: bool,
    rmc6236_final_nm: float,
    rmc6236_stock_mm: float,
) -> pd.DataFrame:
    result = density_table.copy()
    final_volume_ul = final_volume_ml_per_well * 1000
    multiplier = wells_per_timepoint * (1 + extra_fraction)

    result["细胞悬液_uL_per_well"] = (
        result["cells_per_well"] / cell_concentration_per_ml * 1000
    )

    if include_rmc6236 and rmc6236_stock_mm > 0:
        stock_nm = rmc6236_stock_mm * 1_000_000
        result["RMC6236_stock_uL_per_well"] = rmc6236_final_nm * final_volume_ul / stock_nm
    else:
        result["RMC6236_stock_uL_per_well"] = 0.0

    result["补培养基_uL_per_well"] = (
        final_volume_ul
        - result["细胞悬液_uL_per_well"]
        - result["RMC6236_stock_uL_per_well"]
    )
    result["wells"] = wells_per_timepoint
    result["extra_percent"] = extra_fraction * 100
    result["总细胞数_with_extra"] = result["cells_per_well"] * multiplier
    result["总终体积_mL_with_extra"] = final_volume_ml_per_well * multiplier
    result["总细胞悬液_mL_with_extra"] = result["细胞悬液_uL_per_well"] * multiplier / 1000
    result["总培养基_mL_with_extra"] = result["补培养基_uL_per_well"] * multiplier / 1000
    result["总RMC6236_stock_uL_with_extra"] = (
        result["RMC6236_stock_uL_per_well"] * multiplier
    )

    ordered_columns = [
        "时间点",
        "cells_per_well",
        "细胞悬液_uL_per_well",
        "补培养基_uL_per_well",
        "RMC6236_stock_uL_per_well",
        "wells",
        "extra_percent",
        "总细胞数_with_extra",
        "总终体积_mL_with_extra",
        "总细胞悬液_mL_with_extra",
        "总培养基_mL_with_extra",
        "总RMC6236_stock_uL_with_extra",
    ]
    return result[ordered_columns]


def recommend_reseed_container(remaining_cells: float) -> str:
    if remaining_cells < 200_000:
        return "建议 T25 低密度保种；不建议种 T75。"
    if remaining_cells < 800_000:
        return "建议 T25。"
    if remaining_cells < 2_000_000:
        return "建议 T75，或拆成 2 个 T25。"
    return "建议 T75；如果细胞很多，可以同时冻存一管。"


def calculate_reseeding_plan(
    growth_preset: str,
    target_goal: str,
    cell_concentration_per_ml: float,
    remaining_cells: float,
    flask_count: int,
    use_doubling_time: bool,
    doubling_time_hours: float,
) -> pd.DataFrame:
    rows = []
    culture_days = GOAL_DAYS[target_goal]

    for flask_type in ["T25", "T75"]:
        low, high = RESEEDING_DEFAULTS[growth_preset][target_goal][flask_type]
        experience_midpoint = (low + high) / 2
        if use_doubling_time and doubling_time_hours > 0:
            recommended_cells = FULL_FLASK_CELLS[flask_type] / (
                2 ** (culture_days * 24 / doubling_time_hours)
            )
        else:
            recommended_cells = experience_midpoint

        medium_ml = 15.0 if flask_type == "T75" and target_goal == "4–5 天左右长满" else 12.0
        if flask_type == "T25":
            medium_ml = 5.0

        cell_suspension_ul = recommended_cells / cell_concentration_per_ml * 1000
        medium_to_add_ml = medium_ml - cell_suspension_ul / 1000
        total_needed_cells = recommended_cells * flask_count
        max_flasks = int(remaining_cells // recommended_cells) if recommended_cells > 0 else 0

        rows.append(
            {
                "瓶型": flask_type,
                "目标长满时间": target_goal,
                "经验下限_cells_per_flask": low,
                "经验上限_cells_per_flask": high,
                "推荐_cells_per_flask": recommended_cells,
                "培养基终体积_mL_per_flask": medium_ml,
                "细胞悬液_uL_per_flask": cell_suspension_ul,
                "补培养基_mL_per_flask": medium_to_add_ml,
                "计划瓶数": flask_count,
                "需要总细胞数": total_needed_cells,
                "需要总细胞悬液_mL": cell_suspension_ul * flask_count / 1000,
                "剩余细胞数": remaining_cells,
                "剩余_minus_需要": remaining_cells - total_needed_cells,
                "按剩余细胞最多可种瓶数": max_flasks,
                "是否足够": "足够" if remaining_cells >= total_needed_cells else "不够",
            }
        )

    return pd.DataFrame(rows)


def to_csv_bytes(result: pd.DataFrame, settings: dict[str, object]) -> bytes:
    buffer = StringIO()
    buffer.write("# RMC-6236 cell seeding calculator\n")
    for key, value in settings.items():
        buffer.write(f"# {key}: {value}\n")
    buffer.write("\n")
    result.to_csv(buffer, index=False)
    return buffer.getvalue().encode("utf-8-sig")


st.set_page_config(
    page_title="RMC-6236 Seeding Calculator",
    page_icon="🧪",
    layout="wide",
)

st.title("RMC-6236 细胞铺板计算器")
st.caption("按计数后的 cells/mL 自动计算每孔细胞悬液、补培养基、总细胞量和可选 RMC-6236 加药体积。")

with st.sidebar:
    st.header("输入参数")
    cell_concentration = st.number_input(
        "计数后的细胞浓度 (cells/mL)",
        min_value=1.0,
        value=1_200_000.0,
        step=50_000.0,
        format="%.0f",
    )
    preset = st.selectbox("细胞生长类型", list(DEFAULT_DENSITIES.keys()), index=1)
    final_volume_ml = st.number_input(
        "终体积 (mL/well)",
        min_value=0.05,
        value=2.0,
        step=0.1,
        format="%.2f",
    )
    wells_per_timepoint = st.number_input(
        "每个时间点 wells 数",
        min_value=1,
        value=1,
        step=1,
    )
    extra_percent = st.slider("额外配液余量 (%)", min_value=0, max_value=50, value=15, step=1)

    st.divider()
    include_rmc6236 = st.checkbox("计算 RMC-6236 stock 加药体积", value=False)
    rmc6236_final_nm = st.number_input(
        "RMC-6236 终浓度 (nM)",
        min_value=0.0,
        value=100.0,
        step=10.0,
        disabled=not include_rmc6236,
    )
    rmc6236_stock_mm = st.number_input(
        "RMC-6236 stock 浓度 (mM)",
        min_value=0.001,
        value=10.0,
        step=1.0,
        disabled=not include_rmc6236,
    )

st.subheader("1. 自定义每个时间点 cells/well")
edited_density = st.data_editor(
    build_default_table(preset),
    use_container_width=True,
    hide_index=True,
    column_config={
        "时间点": st.column_config.TextColumn("时间点", disabled=True),
        "cells_per_well": st.column_config.NumberColumn(
            "cells/well",
            min_value=1,
            step=1_000,
            format="%d",
        ),
    },
)

plan = calculate_plan(
    density_table=edited_density,
    cell_concentration_per_ml=cell_concentration,
    final_volume_ml_per_well=final_volume_ml,
    wells_per_timepoint=int(wells_per_timepoint),
    extra_fraction=extra_percent / 100,
    include_rmc6236=include_rmc6236,
    rmc6236_final_nm=rmc6236_final_nm,
    rmc6236_stock_mm=rmc6236_stock_mm,
)

negative_medium = plan["补培养基_uL_per_well"] < 0
tiny_drug_volume = include_rmc6236 and (plan["RMC6236_stock_uL_per_well"] < 0.5).any()

if negative_medium.any():
    bad_timepoints = ", ".join(plan.loc[negative_medium, "时间点"].tolist())
    st.error(
        f"{bad_timepoints} 的细胞悬液/药物体积已经超过终体积。请降低 cells/well、提高细胞浓度或增加终体积。"
    )

if tiny_drug_volume:
    st.warning("单孔 RMC-6236 stock 体积低于 0.5 µL，建议先配中间工作液或 master mix 后再加药。")

st.subheader("2. 输出计算表")
display_plan = plan.copy()
rounding_columns = [
    "细胞悬液_uL_per_well",
    "补培养基_uL_per_well",
    "RMC6236_stock_uL_per_well",
    "总终体积_mL_with_extra",
    "总细胞悬液_mL_with_extra",
    "总培养基_mL_with_extra",
    "总RMC6236_stock_uL_with_extra",
]
display_plan[rounding_columns] = display_plan[rounding_columns].round(3)
display_plan["总细胞数_with_extra"] = display_plan["总细胞数_with_extra"].round(0).astype(int)

st.dataframe(display_plan, use_container_width=True, hide_index=True)

total_cells = int(plan["总细胞数_with_extra"].sum().round())
total_volume_ml = float(plan["总终体积_mL_with_extra"].sum())
total_cell_suspension_ml = float(plan["总细胞悬液_mL_with_extra"].sum())
total_medium_ml = float(plan["总培养基_mL_with_extra"].sum())

metric_cols = st.columns(4)
metric_cols[0].metric("总细胞数", f"{total_cells:,}")
metric_cols[1].metric("总终体积", f"{total_volume_ml:.2f} mL")
metric_cols[2].metric("总细胞悬液", f"{total_cell_suspension_ml:.2f} mL")
metric_cols[3].metric("总培养基", f"{total_medium_ml:.2f} mL")

settings = {
    "run_date": datetime.now().isoformat(timespec="seconds"),
    "cell_concentration_cells_per_ml": int(cell_concentration),
    "growth_preset": preset,
    "final_volume_ml_per_well": final_volume_ml,
    "wells_per_timepoint": int(wells_per_timepoint),
    "extra_percent": extra_percent,
    "include_rmc6236": include_rmc6236,
    "rmc6236_final_nm": rmc6236_final_nm if include_rmc6236 else "NA",
    "rmc6236_stock_mm": rmc6236_stock_mm if include_rmc6236 else "NA",
}

csv_name = f"{datetime.now().strftime('%Y-%m-%d')}_RMC6236_seeding_plan.csv"
st.download_button(
    "下载 CSV",
    data=to_csv_bytes(display_plan, settings),
    file_name=csv_name,
    mime="text/csv",
)

st.divider()
st.subheader("3. WB 后剩余细胞回种保种")

reseed_cols = st.columns([1.1, 1.0, 1.0, 1.0])
with reseed_cols[0]:
    remaining_cells = st.number_input(
        "WB 铺板后剩余细胞数",
        min_value=0.0,
        value=500_000.0,
        step=50_000.0,
        format="%.0f",
    )
with reseed_cols[1]:
    reseeding_goal = st.selectbox("希望几天左右长满", RESEEDING_GOALS, index=1)
with reseed_cols[2]:
    reseeding_flask_count = st.number_input(
        "计划瓶数",
        min_value=1,
        value=1,
        step=1,
    )
with reseed_cols[3]:
    use_doubling_time = st.checkbox("用 doubling time 估算", value=False)

if use_doubling_time:
    doubling_time_hours = st.number_input(
        "Doubling time (小时)",
        min_value=1.0,
        value=DEFAULT_DOUBLING_TIME_HOURS[preset],
        step=1.0,
        format="%.1f",
    )
else:
    doubling_time_hours = DEFAULT_DOUBLING_TIME_HOURS[preset]

reseeding_plan = calculate_reseeding_plan(
    growth_preset=preset,
    target_goal=reseeding_goal,
    cell_concentration_per_ml=cell_concentration,
    remaining_cells=remaining_cells,
    flask_count=int(reseeding_flask_count),
    use_doubling_time=use_doubling_time,
    doubling_time_hours=doubling_time_hours,
)

st.info(recommend_reseed_container(remaining_cells))

display_reseeding_plan = reseeding_plan.copy()
reseed_integer_columns = [
    "经验下限_cells_per_flask",
    "经验上限_cells_per_flask",
    "推荐_cells_per_flask",
    "需要总细胞数",
    "剩余细胞数",
    "剩余_minus_需要",
]
display_reseeding_plan[reseed_integer_columns] = (
    display_reseeding_plan[reseed_integer_columns].round(0).astype(int)
)
reseed_float_columns = [
    "培养基终体积_mL_per_flask",
    "细胞悬液_uL_per_flask",
    "补培养基_mL_per_flask",
    "需要总细胞悬液_mL",
]
display_reseeding_plan[reseed_float_columns] = display_reseeding_plan[
    reseed_float_columns
].round(3)

if (display_reseeding_plan["补培养基_mL_per_flask"] < 0).any():
    st.error("回种计算中细胞悬液体积超过瓶内终体积。请提高细胞浓度、降低回种细胞数或增加终体积。")

st.dataframe(display_reseeding_plan, use_container_width=True, hide_index=True)

reseeding_settings = {
    **settings,
    "reseeding_remaining_cells": int(remaining_cells),
    "reseeding_goal": reseeding_goal,
    "reseeding_flask_count": int(reseeding_flask_count),
    "reseeding_use_doubling_time": use_doubling_time,
    "reseeding_doubling_time_hours": doubling_time_hours if use_doubling_time else "NA",
}

reseeding_csv_name = f"{datetime.now().strftime('%Y-%m-%d')}_RMC6236_reseeding_plan.csv"
st.download_button(
    "下载回种 CSV",
    data=to_csv_bytes(display_reseeding_plan, reseeding_settings),
    file_name=reseeding_csv_name,
    mime="text/csv",
)

with st.expander("本次计算参数"):
    st.json(reseeding_settings)
