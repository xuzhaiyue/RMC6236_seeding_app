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

with st.expander("本次计算参数"):
    st.json(settings)
