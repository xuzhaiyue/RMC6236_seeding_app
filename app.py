from __future__ import annotations

import math
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
WAN = 10_000


def build_default_table(preset: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "时间点": TIMEPOINTS,
            "cells_per_well": DEFAULT_DENSITIES[preset],
        }
    )


def build_display_density_table(preset: str) -> pd.DataFrame:
    table = build_default_table(preset)
    return table.assign(cells_per_well_万=table["cells_per_well"] / WAN)[
        ["时间点", "cells_per_well_万"]
    ]


def add_wan_columns(
    table: pd.DataFrame,
    cell_columns: list[str],
    drop_original: bool = True,
) -> pd.DataFrame:
    result = table.copy()
    for column in cell_columns:
        if column in result.columns:
            result[f"{column}_万"] = result[column] / WAN
    if drop_original:
        result = result.drop(columns=[column for column in cell_columns if column in result.columns])
    return result


def calculate_plan(
    density_table: pd.DataFrame,
    cell_concentration_per_ml: float,
    final_volume_ml_per_well: float,
    wells_per_timepoint: int,
    extra_fraction: float,
) -> pd.DataFrame:
    result = density_table.copy()
    final_volume_ul = final_volume_ml_per_well * 1000
    multiplier = wells_per_timepoint * (1 + extra_fraction)

    result["细胞悬液_uL_per_well"] = (
        result["cells_per_well"] / cell_concentration_per_ml * 1000
    )
    result["补培养基_uL_per_well"] = final_volume_ul - result["细胞悬液_uL_per_well"]
    result["wells"] = wells_per_timepoint
    result["extra_percent"] = extra_fraction * 100
    result["总细胞数_with_extra"] = result["cells_per_well"] * multiplier
    result["总终体积_mL_with_extra"] = final_volume_ml_per_well * multiplier
    result["总细胞悬液_mL_with_extra"] = result["细胞悬液_uL_per_well"] * multiplier / 1000
    result["总培养基_mL_with_extra"] = result["补培养基_uL_per_well"] * multiplier / 1000

    ordered_columns = [
        "时间点",
        "cells_per_well",
        "细胞悬液_uL_per_well",
        "补培养基_uL_per_well",
        "wells",
        "extra_percent",
        "总细胞数_with_extra",
        "总终体积_mL_with_extra",
        "总细胞悬液_mL_with_extra",
        "总培养基_mL_with_extra",
    ]
    return result[ordered_columns]


def round_up_to_increment(value: float, increment: float) -> float:
    return math.ceil(value / increment) * increment


def calculate_drug_medium_plan(
    stock_mm: float,
    stock_input_ul: float,
    intermediate_medium_ml: float,
    target_nm: float,
    drug_wells: int,
    volume_ml_per_well: float,
    loss_percent: float,
    round_to_ml: float,
) -> pd.DataFrame:
    intermediate_total_ml = intermediate_medium_ml + stock_input_ul / 1000
    intermediate_um = stock_mm * stock_input_ul / intermediate_total_ml
    raw_total_ml = drug_wells * volume_ml_per_well * (1 + loss_percent / 100)
    final_total_ml = round_up_to_increment(raw_total_ml, round_to_ml)
    intermediate_needed_ul = target_nm * final_total_ml * 1000 / (intermediate_um * 1000)
    medium_needed_ml = final_total_ml - intermediate_needed_ul / 1000

    return pd.DataFrame(
        [
            {
                "stock浓度_mM": stock_mm,
                "第一步_stock原液_uL": stock_input_ul,
                "第一步_培养基_mL": intermediate_medium_ml,
                "得到中间液_uM": intermediate_um,
                "目标终浓度_nM": target_nm,
                "换液孔数": drug_wells,
                "每孔换液_mL": volume_ml_per_well,
                "损耗_percent": loss_percent,
                "建议配制总量_mL": final_total_ml,
                "第二步_加中间液_uL": intermediate_needed_ul,
                "第二步_补培养基_mL": medium_needed_ml,
                "中间液是否足够": "足够" if intermediate_needed_ul <= intermediate_total_ml * 1000 else "不够",
            }
        ]
    )


def calculate_reseeding_plan(
    growth_preset: str,
    target_goal: str,
    cell_concentration_per_ml: float,
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
st.caption("按计数后的万 cells/mL 自动计算铺板、RMC-6236 含药培养基配置和 WB 后回种推荐。")

with st.sidebar:
    st.header("输入参数")
    cell_concentration_wan = st.number_input(
        "计数后的细胞浓度 (万 cells/mL)",
        min_value=0.0001,
        value=120.0,
        step=5.0,
        format="%.2f",
    )
    cell_concentration = cell_concentration_wan * WAN
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

st.subheader("1. 自定义每个时间点铺板量")
edited_density_display = st.data_editor(
    build_display_density_table(preset),
    use_container_width=True,
    hide_index=True,
    column_config={
        "时间点": st.column_config.TextColumn("时间点", disabled=True),
        "cells_per_well_万": st.column_config.NumberColumn(
            "万 cells/well",
            min_value=0.0001,
            step=1.0,
            format="%.2f",
        ),
    },
)
edited_density = edited_density_display.assign(
    cells_per_well=edited_density_display["cells_per_well_万"] * WAN
)[["时间点", "cells_per_well"]]

plan = calculate_plan(
    density_table=edited_density,
    cell_concentration_per_ml=cell_concentration,
    final_volume_ml_per_well=final_volume_ml,
    wells_per_timepoint=int(wells_per_timepoint),
    extra_fraction=extra_percent / 100,
)

negative_medium = plan["补培养基_uL_per_well"] < 0

if negative_medium.any():
    bad_timepoints = ", ".join(plan.loc[negative_medium, "时间点"].tolist())
    st.error(
        f"{bad_timepoints} 的细胞悬液体积已经超过终体积。请降低万 cells/well、提高细胞浓度或增加终体积。"
    )

st.subheader("2. 输出计算表")
display_plan = plan.copy()
display_plan = add_wan_columns(
    display_plan,
    ["cells_per_well", "总细胞数_with_extra"],
    drop_original=True,
)
display_plan = display_plan.rename(
    columns={
        "cells_per_well_万": "铺板量_万cells_per_well",
        "总细胞数_with_extra_万": "总细胞数_万_with_extra",
    }
)
rounding_columns = [
    "铺板量_万cells_per_well",
    "细胞悬液_uL_per_well",
    "补培养基_uL_per_well",
    "总细胞数_万_with_extra",
    "总终体积_mL_with_extra",
    "总细胞悬液_mL_with_extra",
    "总培养基_mL_with_extra",
]
display_plan[rounding_columns] = display_plan[rounding_columns].round(3)
display_plan = display_plan[
    [
        "时间点",
        "铺板量_万cells_per_well",
        "细胞悬液_uL_per_well",
        "补培养基_uL_per_well",
        "wells",
        "extra_percent",
        "总细胞数_万_with_extra",
        "总终体积_mL_with_extra",
        "总细胞悬液_mL_with_extra",
        "总培养基_mL_with_extra",
    ]
]

st.dataframe(display_plan, use_container_width=True, hide_index=True)

total_cells = int(plan["总细胞数_with_extra"].sum().round())
total_volume_ml = float(plan["总终体积_mL_with_extra"].sum())
total_cell_suspension_ml = float(plan["总细胞悬液_mL_with_extra"].sum())
total_medium_ml = float(plan["总培养基_mL_with_extra"].sum())

metric_cols = st.columns(4)
metric_cols[0].metric("总细胞数", f"{total_cells / WAN:.2f} 万")
metric_cols[1].metric("总终体积", f"{total_volume_ml:.2f} mL")
metric_cols[2].metric("总细胞悬液", f"{total_cell_suspension_ml:.2f} mL")
metric_cols[3].metric("总培养基", f"{total_medium_ml:.2f} mL")

settings = {
    "run_date": datetime.now().isoformat(timespec="seconds"),
    "cell_concentration_cells_per_ml": int(cell_concentration),
    "cell_concentration_万_cells_per_ml": cell_concentration_wan,
    "growth_preset": preset,
    "final_volume_ml_per_well": final_volume_ml,
    "wells_per_timepoint": int(wells_per_timepoint),
    "extra_percent": extra_percent,
}

csv_name = f"{datetime.now().strftime('%Y-%m-%d')}_RMC6236_seeding_plan.csv"
st.download_button(
    "下载 CSV",
    data=to_csv_bytes(display_plan, settings),
    file_name=csv_name,
    mime="text/csv",
)

st.divider()
st.subheader("3. RMC-6236 含药培养基配置")

drug_cols_top = st.columns(4)
with drug_cols_top[0]:
    drug_stock_mm = st.number_input(
        "RMC-6236 stock (mM)",
        min_value=0.001,
        value=10.0,
        step=1.0,
    )
with drug_cols_top[1]:
    drug_target_nm = st.number_input(
        "目标终浓度 (nM)",
        min_value=0.0,
        value=100.0,
        step=10.0,
    )
with drug_cols_top[2]:
    drug_wells = st.number_input(
        "要换液的孔数",
        min_value=1,
        value=6,
        step=1,
    )
with drug_cols_top[3]:
    drug_loss_percent = st.slider("加药配液损耗 (%)", min_value=0, max_value=50, value=15, step=5)

drug_cols_bottom = st.columns(4)
with drug_cols_bottom[0]:
    stock_input_ul = st.number_input(
        "第一步 stock 原液 (µL)",
        min_value=2.0,
        value=2.0,
        step=0.5,
    )
with drug_cols_bottom[1]:
    intermediate_medium_ml = st.number_input(
        "第一步培养基 (mL)",
        min_value=0.1,
        value=2.0,
        step=0.1,
    )
with drug_cols_bottom[2]:
    drug_volume_ml_per_well = st.number_input(
        "每孔换液体积 (mL)",
        min_value=0.05,
        value=final_volume_ml,
        step=0.1,
    )
with drug_cols_bottom[3]:
    round_to_ml = st.number_input(
        "总量向上取整到 (mL)",
        min_value=0.1,
        value=1.0,
        step=0.5,
    )

drug_plan = calculate_drug_medium_plan(
    stock_mm=drug_stock_mm,
    stock_input_ul=stock_input_ul,
    intermediate_medium_ml=intermediate_medium_ml,
    target_nm=drug_target_nm,
    drug_wells=int(drug_wells),
    volume_ml_per_well=drug_volume_ml_per_well,
    loss_percent=drug_loss_percent,
    round_to_ml=round_to_ml,
)

display_drug_plan = drug_plan.copy()
display_drug_plan[
    [
        "stock浓度_mM",
        "第一步_stock原液_uL",
        "第一步_培养基_mL",
        "得到中间液_uM",
        "目标终浓度_nM",
        "每孔换液_mL",
        "损耗_percent",
        "建议配制总量_mL",
        "第二步_加中间液_uL",
        "第二步_补培养基_mL",
    ]
] = display_drug_plan[
    [
        "stock浓度_mM",
        "第一步_stock原液_uL",
        "第一步_培养基_mL",
        "得到中间液_uM",
        "目标终浓度_nM",
        "每孔换液_mL",
        "损耗_percent",
        "建议配制总量_mL",
        "第二步_加中间液_uL",
        "第二步_补培养基_mL",
    ]
].round(3)

st.dataframe(display_drug_plan, use_container_width=True, hide_index=True)

drug_summary = display_drug_plan.iloc[0]
if drug_summary["中间液是否足够"] == "不够":
    st.error("第一步配出的 uM 中间液不够。请增加第一步 stock 原液和培养基体积。")
else:
    st.success(
        f"操作：先用 {drug_summary['第一步_stock原液_uL']:.1f} µL stock + "
        f"{drug_summary['第一步_培养基_mL']:.1f} mL 培养基配成 "
        f"{drug_summary['得到中间液_uM']:.3g} µM 中间液；再取 "
        f"{drug_summary['第二步_加中间液_uL']:.1f} µL 中间液，加培养基到 "
        f"{drug_summary['建议配制总量_mL']:.1f} mL。"
    )

drug_settings = {
    **settings,
    "drug_stock_mm": drug_stock_mm,
    "drug_target_nm": drug_target_nm,
    "drug_wells": int(drug_wells),
    "drug_loss_percent": drug_loss_percent,
    "drug_stock_input_ul": stock_input_ul,
    "drug_intermediate_medium_ml": intermediate_medium_ml,
    "drug_volume_ml_per_well": drug_volume_ml_per_well,
    "drug_round_to_ml": round_to_ml,
}

drug_csv_name = f"{datetime.now().strftime('%Y-%m-%d')}_RMC6236_drug_medium_plan.csv"
st.download_button(
    "下载加药配置 CSV",
    data=to_csv_bytes(display_drug_plan, drug_settings),
    file_name=drug_csv_name,
    mime="text/csv",
)

st.divider()
st.subheader("4. WB 后细胞回种推荐")

reseed_cols = st.columns(3)
with reseed_cols[0]:
    reseeding_goal = st.selectbox("希望几天左右长满", RESEEDING_GOALS, index=1)
with reseed_cols[1]:
    reseeding_flask_count = st.number_input(
        "计划瓶数",
        min_value=1,
        value=1,
        step=1,
    )
with reseed_cols[2]:
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
    flask_count=int(reseeding_flask_count),
    use_doubling_time=use_doubling_time,
    doubling_time_hours=doubling_time_hours,
)

display_reseeding_plan = reseeding_plan.copy()
display_reseeding_plan = add_wan_columns(
    display_reseeding_plan,
    [
        "经验下限_cells_per_flask",
        "经验上限_cells_per_flask",
        "推荐_cells_per_flask",
        "需要总细胞数",
    ],
    drop_original=True,
)
display_reseeding_plan = display_reseeding_plan.rename(
    columns={
        "经验下限_cells_per_flask_万": "经验下限_万cells_per_flask",
        "经验上限_cells_per_flask_万": "经验上限_万cells_per_flask",
        "推荐_cells_per_flask_万": "推荐_万cells_per_flask",
        "需要总细胞数_万": "需要总细胞数_万",
    }
)
reseed_wan_columns = [
    "经验下限_万cells_per_flask",
    "经验上限_万cells_per_flask",
    "推荐_万cells_per_flask",
    "需要总细胞数_万",
]
display_reseeding_plan[reseed_wan_columns] = display_reseeding_plan[reseed_wan_columns].round(3)
reseed_float_columns = [
    "培养基终体积_mL_per_flask",
    "细胞悬液_uL_per_flask",
    "补培养基_mL_per_flask",
    "需要总细胞悬液_mL",
]
display_reseeding_plan[reseed_float_columns] = display_reseeding_plan[
    reseed_float_columns
].round(3)
display_reseeding_plan = display_reseeding_plan[
    [
        "瓶型",
        "目标长满时间",
        "经验下限_万cells_per_flask",
        "经验上限_万cells_per_flask",
        "推荐_万cells_per_flask",
        "培养基终体积_mL_per_flask",
        "细胞悬液_uL_per_flask",
        "补培养基_mL_per_flask",
        "计划瓶数",
        "需要总细胞数_万",
        "需要总细胞悬液_mL",
    ]
]

if (display_reseeding_plan["补培养基_mL_per_flask"] < 0).any():
    st.error("回种计算中细胞悬液体积超过瓶内终体积。请提高细胞浓度、降低回种细胞数或增加终体积。")

st.dataframe(display_reseeding_plan, use_container_width=True, hide_index=True)

reseeding_settings = {
    **drug_settings,
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
