"""Plot the participant-level association between UPDRS-III and no-response rate."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

DATA_SHEET = "main_qc_data"
FIXED_SHEET = "fixed_adjusted"
TARGET_ORDER = ["STN", "SN"]
COLORS = {"STN": "#1D4ED8", "SN": "#D62828"}
DEFAULT_EXCLUDED_SUBJECTS = ("sub36",)
DEFAULT_Y_MIN = 0.0
DEFAULT_Y_MAX = 30.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a participant-level scatter plot of mean UPDRS-III score "
            "against mean no-response rate across retained stimulation sessions."
        )
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Excel file containing the mixed-model output worksheets.",
    )
    parser.add_argument(
        "--data-sheet",
        default=DATA_SHEET,
        help=f"Worksheet containing session-level data (default: {DATA_SHEET}).",
    )
    parser.add_argument(
        "--fixed-sheet",
        default=FIXED_SHEET,
        help=f"Fixed-effects worksheet retained for input validation (default: {FIXED_SHEET}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory for output figures. By default, a folder named "
            "'updrs_no_response_association' is created beside the input file."
        ),
    )
    parser.add_argument(
        "--exclude-subject",
        action="append",
        default=None,
        help=(
            "Subject identifier to exclude. Repeat this option for multiple subjects. "
            f"Default: {', '.join(DEFAULT_EXCLUDED_SUBJECTS)}."
        ),
    )
    parser.add_argument(
        "--y-min",
        type=float,
        default=DEFAULT_Y_MIN,
        help=f"Lower y-axis limit in percent (default: {DEFAULT_Y_MIN:g}).",
    )
    parser.add_argument(
        "--y-max",
        type=float,
        default=DEFAULT_Y_MAX,
        help=f"Upper y-axis limit in percent (default: {DEFAULT_Y_MAX:g}).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure after saving it.",
    )
    return parser.parse_args()


def to_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def require_columns(data: pd.DataFrame, required: list[str], sheet_name: str) -> None:
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(
            f"Worksheet '{sheet_name}' is missing required columns: {missing}. "
            f"Available columns: {list(data.columns)}"
        )


def load_data(
    input_file: Path,
    data_sheet: str,
    fixed_sheet: str,
    excluded_subjects: list[str],
) -> pd.DataFrame:
    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    data = pd.read_excel(input_file, sheet_name=data_sheet)
    pd.read_excel(input_file, sheet_name=fixed_sheet)
    data.columns = data.columns.astype(str).str.strip()

    require_columns(data, ["subject", "Target", "condition", "UPDRS"], data_sheet)

    data["subject"] = data["subject"].astype(str).str.strip()
    data["Target"] = data["Target"].astype(str).str.strip()
    data["condition"] = data["condition"].astype(str).str.strip()

    excluded = {str(subject).strip().lower() for subject in excluded_subjects}
    data = data[~data["subject"].str.lower().isin(excluded)].copy()

    if "analysis_keep" in data.columns:
        data = data[to_bool(data["analysis_keep"])].copy()

    if "no_response_rate" not in data.columns:
        require_columns(data, ["responded_total", "no_response_total"], data_sheet)
        denominator = data["responded_total"] + data["no_response_total"]
        data["no_response_rate"] = np.where(
            denominator > 0,
            data["no_response_total"] / denominator,
            np.nan,
        )

    data["UPDRS"] = pd.to_numeric(data["UPDRS"], errors="coerce")
    data["no_response_rate"] = pd.to_numeric(
        data["no_response_rate"], errors="coerce"
    )
    return data


def make_subject_level_data(data: pd.DataFrame) -> pd.DataFrame:
    subject_data = (
        data.groupby(["subject", "Target"], observed=True)
        .agg(
            mean_UPDRS=("UPDRS", "mean"),
            mean_no_response_rate=("no_response_rate", "mean"),
            n_sessions=("condition", "nunique"),
        )
        .reset_index()
    )

    subject_data = subject_data.dropna(
        subset=["mean_UPDRS", "mean_no_response_rate"]
    ).copy()
    subject_data["mean_no_response_percent"] = (
        subject_data["mean_no_response_rate"] * 100
    )

    if len(subject_data) < 2:
        raise ValueError("At least two participants are required for the regression.")
    return subject_data


def fit_subject_level_regression(
    subject_data: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame, object]:
    x = subject_data["mean_UPDRS"].to_numpy(dtype=float)
    y = subject_data["mean_no_response_percent"].to_numpy(dtype=float)

    design = sm.add_constant(x)
    model = sm.OLS(y, design).fit()

    x_grid = np.linspace(x.min(), x.max(), 200)
    prediction_design = sm.add_constant(x_grid)
    predictions = model.get_prediction(prediction_design).summary_frame(alpha=0.05)
    return x_grid, predictions, model


def configure_plot_style() -> None:
    plt.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["axes.linewidth"] = 1.0
    plt.rcParams["xtick.major.width"] = 1.0
    plt.rcParams["ytick.major.width"] = 1.0
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["ytick.direction"] = "out"


def create_plot(
    subject_data: pd.DataFrame,
    x_grid: np.ndarray,
    predictions: pd.DataFrame,
    y_min: float,
    y_max: float,
) -> plt.Figure:
    configure_plot_style()
    figure, axis = plt.subplots(figsize=(7.2, 4.2))
    figure.subplots_adjust(left=0.18, right=0.98, bottom=0.18, top=0.98)

    for target in TARGET_ORDER:
        subset = subject_data[subject_data["Target"] == target]
        axis.scatter(
            subset["mean_UPDRS"],
            subset["mean_no_response_percent"],
            s=42,
            color=COLORS[target],
            edgecolor="white",
            linewidth=0.7,
            alpha=0.85,
            zorder=3,
        )

    axis.plot(
        x_grid,
        predictions["mean"],
        color="0.15",
        linewidth=2.0,
        zorder=2,
    )
    axis.fill_between(
        x_grid,
        predictions["mean_ci_lower"],
        predictions["mean_ci_upper"],
        color="0.70",
        alpha=0.28,
        linewidth=0,
        zorder=1,
    )

    axis.set_title("")
    axis.set_xlabel("")
    axis.set_ylabel("")
    axis.tick_params(axis="x", labelsize=9, pad=5)
    axis.tick_params(axis="y", labelsize=9, pad=5)
    axis.set_ylim(y_min, y_max)
    axis.set_yticks(np.arange(y_min, y_max + 1, 5))

    x_min = np.floor(subject_data["mean_UPDRS"].min() / 5) * 5 - 2
    x_max = np.ceil(subject_data["mean_UPDRS"].max() / 5) * 5 + 2
    axis.set_xlim(x_min, x_max)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    return figure


def main() -> None:
    args = parse_args()
    excluded_subjects = (
        args.exclude_subject
        if args.exclude_subject is not None
        else list(DEFAULT_EXCLUDED_SUBJECTS)
    )
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else args.input_file.parent / "updrs_no_response_association"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    data = load_data(
        input_file=args.input_file,
        data_sheet=args.data_sheet,
        fixed_sheet=args.fixed_sheet,
        excluded_subjects=excluded_subjects,
    )
    subject_data = make_subject_level_data(data)
    x_grid, predictions, _ = fit_subject_level_regression(subject_data)
    figure = create_plot(
        subject_data=subject_data,
        x_grid=x_grid,
        predictions=predictions,
        y_min=args.y_min,
        y_max=args.y_max,
    )

    output_png = output_dir / "updrs_no_response_association.png"
    output_pdf = output_dir / "updrs_no_response_association.pdf"
    output_svg = output_dir / "updrs_no_response_association.svg"

    figure.savefig(output_png, dpi=600, facecolor="white")
    figure.savefig(output_pdf, facecolor="white")
    figure.savefig(output_svg, facecolor="white")

    if args.show:
        plt.show()
    else:
        plt.close(figure)

    print(f"Saved to:\n{output_png}\n{output_pdf}\n{output_svg}")
    print(f"Display y-axis range: {args.y_min:g} to {args.y_max:g}")


if __name__ == "__main__":
    main()
