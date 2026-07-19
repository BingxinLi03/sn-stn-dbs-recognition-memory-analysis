"""Analyze UPDRS-III scores across DBS targets and stimulation conditions."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, mannwhitneyu, ttest_ind, ttest_rel, wilcoxon

SHEET_NAME = "Sheet1"
VALUE_COLUMN = "UPDRS"
VALUE_LABEL = "UPDRS-III score"
PLOT_TITLE = "UPDRS-III"

GROUP_ORDER = ["STN", "SN"]
CONDITION_ORDER = ["off", "10Hz", "130Hz"]

WITHIN_TEST = "wilcoxon"
BETWEEN_TEST = "mannwhitney"
P_ADJUST_METHOD = "fdr_bh"

SHOW_SIGNIFICANCE = True
SIGNIFICANCE_MODE = "auto"
SIGNIFICANCE_LABEL_STYLE = "stars"
P_VALUE_FOR_PLOT = "p_adj"
PLOT_ONLY_SIGNIFICANT = True
ALPHA = 0.05

PRINT_STATISTICS = True
EXPORT_STATISTICS = True

MANUAL_SIGNIFICANCE: list[dict[str, object]] = []

COLORS = {"STN": "#2F60D8", "SN": "#E3342F"}
X_OFFSET = {"STN": -0.035, "SN": 0.035}

FIGURE_SIZE = (10.5, 7.0)
DPI = 300
BACKGROUND_COLOR = "#eeeeee"
ZERO_LINE_COLOR = "#8a8a8a"

INDIVIDUAL_LINE_ALPHA = 0.14
INDIVIDUAL_POINT_ALPHA = 0.35
INDIVIDUAL_LINE_WIDTH = 1.4
INDIVIDUAL_POINT_SIZE = 46

MEAN_LINE_WIDTH = 4.0
MEAN_POINT_SIZE = 280
MEAN_EDGE_WIDTH = 1.6
ERRORBAR_LINE_WIDTH = 2.6
CAPSIZE = 7

TITLE_SIZE = 22
LABEL_SIZE = 21
TICK_SIZE = 17
SIGNIFICANCE_TEXT_SIZE = 13


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze UPDRS-III scores across OFF, 10 Hz, and 130 Hz conditions "
            "within STN-DBS and SN-DBS groups."
        )
    )
    parser.add_argument(
        "input_file",
        type=Path,
        help="Excel file containing subject, condition, Target, and UPDRS columns.",
    )
    parser.add_argument(
        "--sheet",
        default=SHEET_NAME,
        help=f"Worksheet name (default: {SHEET_NAME}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("updrs_analysis_results"),
        help="Directory for tables and figures (default: updrs_analysis_results).",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure after saving it.",
    )
    return parser.parse_args()


def normalize_condition(value: object) -> str:
    text = str(value).strip().replace(" ", "")
    normalized = text.lower()

    if normalized in {"off", "0", "dbsoff", "offdbs"}:
        return "off"
    if normalized in {"10hz", "10h", "10"}:
        return "10Hz"
    if normalized in {"130hz", "130h", "130"}:
        return "130Hz"
    return text


def sem(values: pd.Series) -> float:
    numeric = pd.Series(values).dropna().astype(float)
    if len(numeric) <= 1:
        return np.nan
    return float(numeric.std(ddof=1) / np.sqrt(len(numeric)))


def p_to_text(p_value: float, style: str = "stars") -> str:
    if pd.isna(p_value):
        return "n.s."

    if style == "p":
        return "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"

    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "n.s."


def adjust_pvalues(p_values: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full_like(p_values, np.nan, dtype=float)

    valid = ~np.isnan(p_values)
    valid_p = p_values[valid]
    count = len(valid_p)

    if count == 0:
        return adjusted

    if method == "none":
        adjusted[valid] = valid_p
        return adjusted

    if method == "bonferroni":
        adjusted[valid] = np.minimum(valid_p * count, 1.0)
        return adjusted

    if method == "fdr_bh":
        order = np.argsort(valid_p)
        ranked = valid_p[order] * count / np.arange(1, count + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        ranked = np.minimum(ranked, 1.0)

        restored = np.empty_like(ranked)
        restored[order] = ranked
        adjusted[valid] = restored
        return adjusted

    raise ValueError("P_ADJUST_METHOD must be 'fdr_bh', 'bonferroni', or 'none'.")


def add_adjusted_p_by_family(
    table: pd.DataFrame,
    p_column: str = "p_raw",
    family_column: str = "family",
    method: str = "fdr_bh",
) -> pd.DataFrame:
    output = table.copy()
    output["p_adj"] = np.nan

    for indices in output.groupby(family_column, observed=True).groups.values():
        output.loc[indices, "p_adj"] = adjust_pvalues(
            output.loc[indices, p_column].values,
            method=method,
        )
    return output


def round_for_output(table: pd.DataFrame, digits: int = 4) -> pd.DataFrame:
    output = table.copy()
    for column in output.columns:
        if pd.api.types.is_float_dtype(output[column]):
            output[column] = output[column].round(digits)
    return output


def add_significance_bracket(
    axis: plt.Axes,
    x1: float,
    x2: float,
    y: float,
    label: str,
    height: float,
    line_width: float = 1.6,
    color: str = "black",
    font_size: int = 12,
) -> None:
    axis.plot(
        [x1, x1, x2, x2],
        [y, y + height, y + height, y],
        lw=line_width,
        c=color,
        clip_on=False,
    )
    axis.text(
        (x1 + x2) / 2,
        y + height,
        label,
        ha="center",
        va="bottom",
        color=color,
        fontsize=font_size,
    )


def load_data(input_file: Path, sheet_name: str) -> pd.DataFrame:
    if not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    data = pd.read_excel(input_file, sheet_name=sheet_name)
    data.columns = data.columns.astype(str).str.strip()

    required_columns = ["subject", "condition", "Target", VALUE_COLUMN]
    missing_columns = [column for column in required_columns if column not in data.columns]
    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}. "
            f"Available columns: {list(data.columns)}"
        )

    data = data[required_columns].copy()
    data["subject"] = data["subject"].astype(str).str.strip()
    data["condition"] = data["condition"].apply(normalize_condition)
    data["Target"] = data["Target"].astype(str).str.strip().str.upper()
    data[VALUE_COLUMN] = pd.to_numeric(data[VALUE_COLUMN], errors="coerce")

    data["condition"] = pd.Categorical(
        data["condition"], categories=CONDITION_ORDER, ordered=True
    )
    data["Target"] = pd.Categorical(
        data["Target"], categories=GROUP_ORDER, ordered=True
    )
    data = data.dropna(
        subset=["subject", "condition", "Target", VALUE_COLUMN]
    ).copy()

    if data.empty:
        raise ValueError("No valid observations remain after data validation.")
    return data


def calculate_descriptive_statistics(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for group in GROUP_ORDER:
        for condition in CONDITION_ORDER:
            values = data.loc[
                (data["Target"] == group) & (data["condition"] == condition),
                VALUE_COLUMN,
            ].dropna().astype(float)

            rows.append(
                {
                    "Target": group,
                    "condition": condition,
                    "n": len(values),
                    "mean": values.mean() if len(values) else np.nan,
                    "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
                    "sem": sem(values),
                    "median": values.median() if len(values) else np.nan,
                    "q1": values.quantile(0.25) if len(values) else np.nan,
                    "q3": values.quantile(0.75) if len(values) else np.nan,
                    "min": values.min() if len(values) else np.nan,
                    "max": values.max() if len(values) else np.nan,
                }
            )

    return pd.DataFrame(rows)


def calculate_within_group_statistics(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for group in GROUP_ORDER:
        group_data = data[data["Target"] == group].copy()
        wide = group_data.pivot_table(
            index="subject",
            columns="condition",
            values=VALUE_COLUMN,
            aggfunc="mean",
            observed=True,
        )

        for condition1, condition2 in combinations(CONDITION_ORDER, 2):
            if condition1 in wide.columns and condition2 in wide.columns:
                paired = wide[[condition1, condition2]].dropna()
            else:
                paired = pd.DataFrame(columns=[condition1, condition2])

            n_pairs = len(paired)
            statistic = np.nan
            p_raw = np.nan
            test_df = np.nan
            mean_condition1 = np.nan
            mean_condition2 = np.nan
            mean_difference = np.nan
            median_difference = np.nan

            if WITHIN_TEST == "wilcoxon":
                test_name = "Wilcoxon signed-rank"
                statistic_name = "W"
            elif WITHIN_TEST == "paired_t":
                test_name = "paired t-test"
                statistic_name = "t"
            else:
                raise ValueError("WITHIN_TEST must be 'wilcoxon' or 'paired_t'.")

            if n_pairs >= 2:
                values1 = paired[condition1].astype(float).values
                values2 = paired[condition2].astype(float).values
                difference = values2 - values1

                mean_condition1 = np.mean(values1)
                mean_condition2 = np.mean(values2)
                mean_difference = np.mean(difference)
                median_difference = np.median(difference)

                if WITHIN_TEST == "wilcoxon":
                    if np.allclose(difference, 0):
                        statistic = 0.0
                        p_raw = 1.0
                    else:
                        result = wilcoxon(
                            values1,
                            values2,
                            alternative="two-sided",
                            zero_method="wilcox",
                        )
                        statistic = float(result.statistic)
                        p_raw = float(result.pvalue)
                else:
                    result = ttest_rel(values1, values2, nan_policy="omit")
                    statistic = float(result.statistic)
                    p_raw = float(result.pvalue)
                    test_df = n_pairs - 1

            rows.append(
                {
                    "family": f"within_{group}",
                    "comparison_type": "within_group",
                    "Target": group,
                    "cond1": condition1,
                    "cond2": condition2,
                    "n_pairs": n_pairs,
                    "test": test_name,
                    "statistic_name": statistic_name,
                    "statistic": statistic,
                    "df": test_df,
                    "p_raw": p_raw,
                    "mean_cond1": mean_condition1,
                    "mean_cond2": mean_condition2,
                    "mean_diff_cond2_minus_cond1": mean_difference,
                    "median_diff_cond2_minus_cond1": median_difference,
                }
            )

    output = add_adjusted_p_by_family(
        pd.DataFrame(rows),
        p_column="p_raw",
        family_column="family",
        method=P_ADJUST_METHOD,
    )
    output["sig_raw"] = output["p_raw"].apply(lambda p: p_to_text(p, "stars"))
    output["sig_adj"] = output["p_adj"].apply(lambda p: p_to_text(p, "stars"))
    return output


def calculate_between_group_statistics(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group1, group2 = GROUP_ORDER

    for condition in CONDITION_ORDER:
        values1 = data.loc[
            (data["Target"] == group1) & (data["condition"] == condition),
            VALUE_COLUMN,
        ].dropna().astype(float).values
        values2 = data.loc[
            (data["Target"] == group2) & (data["condition"] == condition),
            VALUE_COLUMN,
        ].dropna().astype(float).values

        statistic = np.nan
        p_raw = np.nan
        test_df = np.nan

        if BETWEEN_TEST == "mannwhitney":
            test_name = "Mann-Whitney U"
            statistic_name = "U"
        elif BETWEEN_TEST == "welch_t":
            test_name = "Welch t-test"
            statistic_name = "t"
        else:
            raise ValueError("BETWEEN_TEST must be 'mannwhitney' or 'welch_t'.")

        if len(values1) >= 2 and len(values2) >= 2:
            if BETWEEN_TEST == "mannwhitney":
                result = mannwhitneyu(values1, values2, alternative="two-sided")
                statistic = float(result.statistic)
                p_raw = float(result.pvalue)
            else:
                result = ttest_ind(
                    values1, values2, equal_var=False, nan_policy="omit"
                )
                statistic = float(result.statistic)
                p_raw = float(result.pvalue)

                variance1 = np.var(values1, ddof=1)
                variance2 = np.var(values2, ddof=1)
                test_df = (variance1 / len(values1) + variance2 / len(values2)) ** 2 / (
                    (variance1 / len(values1)) ** 2 / (len(values1) - 1)
                    + (variance2 / len(values2)) ** 2 / (len(values2) - 1)
                )

        rows.append(
            {
                "family": "between_groups",
                "comparison_type": "between_groups_same_freq",
                "condition": condition,
                "group1": group1,
                "group2": group2,
                "n_group1": len(values1),
                "n_group2": len(values2),
                "test": test_name,
                "statistic_name": statistic_name,
                "statistic": statistic,
                "df": test_df,
                "p_raw": p_raw,
                "mean_group1": np.mean(values1) if len(values1) else np.nan,
                "mean_group2": np.mean(values2) if len(values2) else np.nan,
                "mean_diff_group1_minus_group2": (
                    np.mean(values1) - np.mean(values2)
                    if len(values1) and len(values2)
                    else np.nan
                ),
                "median_group1": np.median(values1) if len(values1) else np.nan,
                "median_group2": np.median(values2) if len(values2) else np.nan,
                "median_diff_group1_minus_group2": (
                    np.median(values1) - np.median(values2)
                    if len(values1) and len(values2)
                    else np.nan
                ),
            }
        )

    output = add_adjusted_p_by_family(
        pd.DataFrame(rows),
        p_column="p_raw",
        family_column="family",
        method=P_ADJUST_METHOD,
    )
    output["sig_raw"] = output["p_raw"].apply(lambda p: p_to_text(p, "stars"))
    output["sig_adj"] = output["p_adj"].apply(lambda p: p_to_text(p, "stars"))
    return output


def calculate_omnibus_statistics(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for group in GROUP_ORDER:
        group_data = data[data["Target"] == group].copy()
        wide = group_data.pivot_table(
            index="subject",
            columns="condition",
            values=VALUE_COLUMN,
            aggfunc="mean",
            observed=True,
        )

        if all(condition in wide.columns for condition in CONDITION_ORDER):
            complete = wide[CONDITION_ORDER].dropna()
        else:
            complete = pd.DataFrame()

        statistic = np.nan
        p_raw = np.nan

        if len(complete) >= 2:
            result = friedmanchisquare(
                *[complete[condition].astype(float).values for condition in CONDITION_ORDER]
            )
            statistic = float(result.statistic)
            p_raw = float(result.pvalue)

        rows.append(
            {
                "Target": group,
                "test": "Friedman test",
                "statistic_name": "chi-square",
                "statistic": statistic,
                "df": len(CONDITION_ORDER) - 1,
                "n_complete": len(complete),
                "p_raw": p_raw,
                "sig_raw": p_to_text(p_raw, "stars"),
            }
        )

    return pd.DataFrame(rows)


def print_statistics(
    descriptive: pd.DataFrame,
    within_group: pd.DataFrame,
    between_groups: pd.DataFrame,
    omnibus: pd.DataFrame,
) -> None:
    pd.set_option("display.max_columns", 100)
    pd.set_option("display.width", 180)

    print("\n================ Descriptive statistics ================")
    print(round_for_output(descriptive).to_string(index=False))

    print("\n================ Within-group comparisons ================")
    print(round_for_output(within_group).to_string(index=False))

    print("\n================ Between-group comparisons ================")
    print(round_for_output(between_groups).to_string(index=False))

    print("\n================ Friedman omnibus tests ================")
    print(round_for_output(omnibus).to_string(index=False))


def export_statistics(
    output_dir: Path,
    descriptive: pd.DataFrame,
    within_group: pd.DataFrame,
    between_groups: pd.DataFrame,
    omnibus: pd.DataFrame,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    workbook_path = output_dir / "updrs_statistics.xlsx"
    csv_prefix = output_dir / "updrs_statistics"

    with pd.ExcelWriter(workbook_path) as writer:
        round_for_output(descriptive).to_excel(
            writer, sheet_name="descriptive", index=False
        )
        round_for_output(within_group).to_excel(
            writer, sheet_name="within_group", index=False
        )
        round_for_output(between_groups).to_excel(
            writer, sheet_name="between_groups", index=False
        )
        round_for_output(omnibus).to_excel(
            writer, sheet_name="friedman", index=False
        )

    round_for_output(descriptive).to_csv(
        f"{csv_prefix}_descriptive.csv", index=False, encoding="utf-8-sig"
    )
    round_for_output(within_group).to_csv(
        f"{csv_prefix}_within_group.csv", index=False, encoding="utf-8-sig"
    )
    round_for_output(between_groups).to_csv(
        f"{csv_prefix}_between_groups.csv", index=False, encoding="utf-8-sig"
    )
    round_for_output(omnibus).to_csv(
        f"{csv_prefix}_friedman.csv", index=False, encoding="utf-8-sig"
    )
    return workbook_path


def build_significance_items(
    within_group: pd.DataFrame,
    between_groups: pd.DataFrame,
) -> list[dict[str, object]]:
    if not SHOW_SIGNIFICANCE:
        return []
    if SIGNIFICANCE_MODE == "manual":
        return MANUAL_SIGNIFICANCE.copy()
    if SIGNIFICANCE_MODE != "auto":
        raise ValueError("SIGNIFICANCE_MODE must be 'auto' or 'manual'.")

    items: list[dict[str, object]] = []
    for _, row in within_group.iterrows():
        p_value = row[P_VALUE_FOR_PLOT]
        if PLOT_ONLY_SIGNIFICANT and (pd.isna(p_value) or p_value >= ALPHA):
            continue
        items.append(
            {
                "type": "within_group",
                "group": row["Target"],
                "cond1": row["cond1"],
                "cond2": row["cond2"],
                "label": p_to_text(p_value, SIGNIFICANCE_LABEL_STYLE),
            }
        )

    for _, row in between_groups.iterrows():
        p_value = row[P_VALUE_FOR_PLOT]
        if PLOT_ONLY_SIGNIFICANT and (pd.isna(p_value) or p_value >= ALPHA):
            continue
        items.append(
            {
                "type": "between_groups_same_freq",
                "cond": row["condition"],
                "group1": row["group1"],
                "group2": row["group2"],
                "label": p_to_text(p_value, SIGNIFICANCE_LABEL_STYLE),
            }
        )
    return items


def create_plot(
    data: pd.DataFrame,
    descriptive: pd.DataFrame,
    significance_items: list[dict[str, object]],
    output_dir: Path,
    show: bool,
) -> tuple[Path, Path]:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    figure, axis = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI)
    figure.patch.set_facecolor(BACKGROUND_COLOR)
    axis.set_facecolor(BACKGROUND_COLOR)

    x_base = np.arange(len(CONDITION_ORDER))
    x_map = {condition: index for index, condition in enumerate(CONDITION_ORDER)}

    for group in GROUP_ORDER:
        color = COLORS[group]
        group_data = data[data["Target"] == group].copy()

        for _, subject_data in group_data.groupby("subject"):
            subject_values: dict[str, float] = {}
            for condition in CONDITION_ORDER:
                condition_data = subject_data[subject_data["condition"] == condition]
                subject_values[condition] = (
                    condition_data[VALUE_COLUMN].iloc[0]
                    if len(condition_data) > 0
                    else np.nan
                )

            points_x: list[float] = []
            points_y: list[float] = []
            for condition in CONDITION_ORDER:
                value = subject_values[condition]
                if not pd.isna(value):
                    points_x.append(x_map[condition] + X_OFFSET[group])
                    points_y.append(value)

            axis.scatter(
                points_x,
                points_y,
                s=INDIVIDUAL_POINT_SIZE,
                color=color,
                alpha=INDIVIDUAL_POINT_ALPHA,
                edgecolors="none",
                zorder=2,
            )

            for condition1, condition2 in zip(
                CONDITION_ORDER[:-1], CONDITION_ORDER[1:]
            ):
                value1 = subject_values[condition1]
                value2 = subject_values[condition2]
                if not pd.isna(value1) and not pd.isna(value2):
                    axis.plot(
                        [
                            x_map[condition1] + X_OFFSET[group],
                            x_map[condition2] + X_OFFSET[group],
                        ],
                        [value1, value2],
                        color=color,
                        alpha=INDIVIDUAL_LINE_ALPHA,
                        linewidth=INDIVIDUAL_LINE_WIDTH,
                        zorder=1,
                    )

    for group in GROUP_ORDER:
        color = COLORS[group]
        group_summary = descriptive[descriptive["Target"] == group].copy()
        group_summary["x"] = (
            group_summary["condition"].map(x_map).astype(float) + X_OFFSET[group]
        )
        group_summary = group_summary.sort_values("condition")

        axis.errorbar(
            group_summary["x"],
            group_summary["mean"],
            yerr=group_summary["sem"],
            fmt="none",
            ecolor=color,
            elinewidth=ERRORBAR_LINE_WIDTH,
            capsize=CAPSIZE,
            capthick=ERRORBAR_LINE_WIDTH,
            zorder=4,
        )
        axis.plot(
            group_summary["x"],
            group_summary["mean"],
            color=color,
            linewidth=MEAN_LINE_WIDTH,
            zorder=5,
        )
        axis.scatter(
            group_summary["x"],
            group_summary["mean"],
            s=MEAN_POINT_SIZE,
            color=color,
            edgecolor="white",
            linewidth=MEAN_EDGE_WIDTH,
            zorder=6,
        )

    axis.axhline(
        0,
        color=ZERO_LINE_COLOR,
        linestyle=":",
        linewidth=2.0,
        zorder=0,
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(2.0)
    axis.spines["bottom"].set_linewidth(2.0)
    axis.tick_params(
        axis="x", width=2.0, length=9, pad=12, labelsize=TICK_SIZE
    )
    axis.tick_params(
        axis="y", width=2.0, length=9, pad=12, labelsize=TICK_SIZE
    )
    axis.set_xticks(x_base)
    axis.set_xticklabels(CONDITION_ORDER)
    axis.set_xlabel("Stimulation frequency", fontsize=LABEL_SIZE, labelpad=20)
    axis.set_ylabel(VALUE_LABEL, fontsize=LABEL_SIZE, labelpad=22)
    axis.set_title(PLOT_TITLE, fontsize=TITLE_SIZE, pad=35)

    y_min = float(np.nanmin(data[VALUE_COLUMN].values))
    y_max = float(np.nanmax(data[VALUE_COLUMN].values))
    y_range = y_max - y_min if y_max > y_min else 1.0
    lower_limit = min(0, y_min - 0.08 * y_range)
    upper_limit = y_max + 0.15 * y_range

    if SHOW_SIGNIFICANCE and significance_items:
        significance_y_start = y_max + 0.10 * y_range
        significance_y_step = 0.10 * y_range
        significance_bar_height = 0.03 * y_range

        for level, item in enumerate(significance_items):
            y_position = float(
                item.get("y", significance_y_start + level * significance_y_step)
            )

            if item["type"] == "within_group":
                group = str(item["group"])
                x1 = x_map[str(item["cond1"])] + X_OFFSET[group]
                x2 = x_map[str(item["cond2"])] + X_OFFSET[group]
            elif item["type"] == "between_groups_same_freq":
                condition = str(item["cond"])
                group1 = str(item["group1"])
                group2 = str(item["group2"])
                x1 = x_map[condition] + X_OFFSET[group1]
                x2 = x_map[condition] + X_OFFSET[group2]
            else:
                continue

            add_significance_bracket(
                axis,
                x1,
                x2,
                y_position,
                str(item["label"]),
                height=significance_bar_height,
                font_size=SIGNIFICANCE_TEXT_SIZE,
            )

        upper_limit = y_max + (len(significance_items) + 2) * significance_y_step

    axis.set_ylim(lower_limit, upper_limit)
    figure.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "updrs_condition_plot.png"
    pdf_path = output_dir / "updrs_condition_plot.pdf"
    figure.savefig(
        png_path,
        dpi=DPI,
        bbox_inches="tight",
        facecolor=figure.get_facecolor(),
    )
    figure.savefig(
        pdf_path,
        bbox_inches="tight",
        facecolor=figure.get_facecolor(),
    )

    if show:
        plt.show()
    plt.close(figure)
    return png_path, pdf_path


def main() -> None:
    args = parse_args()
    input_file = args.input_file.expanduser()
    output_dir = args.output_dir.expanduser()

    data = load_data(input_file, args.sheet)
    descriptive = calculate_descriptive_statistics(data)
    within_group = calculate_within_group_statistics(data)
    between_groups = calculate_between_group_statistics(data)
    omnibus = calculate_omnibus_statistics(data)

    if PRINT_STATISTICS:
        print_statistics(descriptive, within_group, between_groups, omnibus)

    workbook_path: Path | None = None
    if EXPORT_STATISTICS:
        workbook_path = export_statistics(
            output_dir,
            descriptive,
            within_group,
            between_groups,
            omnibus,
        )

    significance_items = build_significance_items(within_group, between_groups)
    png_path, pdf_path = create_plot(
        data,
        descriptive,
        significance_items,
        output_dir,
        show=args.show,
    )

    print(f"\nPNG figure: {png_path}")
    print(f"PDF figure: {pdf_path}")
    if workbook_path is not None:
        print(f"Statistics workbook: {workbook_path}")


if __name__ == "__main__":
    main()
