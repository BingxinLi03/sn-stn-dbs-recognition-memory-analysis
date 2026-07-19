"""Behavioral mixed-effects analysis with UPDRS-III adjustment.

The script applies the prespecified main and sensitivity quality-control
thresholds, fits unadjusted and UPDRS-III-adjusted mixed-effects models,
computes planned contrasts with Benjamini-Hochberg FDR correction, and writes
result tables, summaries, and diagnostic figures.
"""

import argparse
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

MAIN_QC = {"responded_total": 24, "old_sdt_total": 8, "new_sdt_total": 8}
SENS_QC = {"responded_total": 36, "old_sdt_total": 12, "new_sdt_total": 12}

PRIMARY_OUTCOME = "d_prime"
SECONDARY_OUTCOMES = [
    "pr",
    "hit_rate",
    "false_alarm_rate",
    "no_response_rate",
    "log_rt_hit_mean",
    "log_rt_correct_rejection_mean",
]

TARGET_ORDER = ["STN", "SN"]
CONDITION_ORDER = ["off", "10Hz", "130Hz"]

warnings.filterwarnings("ignore")


def prepare_behavior(path):
    df = pd.read_excel(path, sheet_name="summary")
    df = df.drop(columns=["file"], errors="ignore")
    df["subject"] = df["subject"].astype(str).str.strip()
    df["condition"] = df["condition"].astype(str).str.strip()
    df["Target"] = df["Target"].astype(str).str.strip()
    df["Target"] = pd.Categorical(df["Target"], categories=TARGET_ORDER, ordered=True)
    df["condition"] = pd.Categorical(df["condition"], categories=CONDITION_ORDER, ordered=True)

    denom = df["responded_total"] + df["no_response_total"]
    df["no_response_rate"] = np.where(denom > 0, df["no_response_total"] / denom, np.nan)

    df["log_rt_hit_mean"] = np.where(df["rt_hit_mean"] > 0, np.log(df["rt_hit_mean"]), np.nan)
    df["log_rt_correct_rejection_mean"] = np.where(
        df["rt_correct_rejection_mean"] > 0,
        np.log(df["rt_correct_rejection_mean"]),
        np.nan,
    )
    return df


def prepare_updrs(path):
    df = pd.read_excel(path)
    df["subject"] = df["subject"].astype(str).str.strip()
    df["condition"] = df["condition"].astype(str).str.strip()
    df["Target"] = df["Target"].astype(str).str.strip()
    keep_cols = ["subject", "condition", "Target", "UPDRS"]
    df = df[keep_cols].copy()
    df["UPDRS"] = pd.to_numeric(df["UPDRS"], errors="coerce")
    df["Target"] = pd.Categorical(df["Target"], categories=TARGET_ORDER, ordered=True)
    df["condition"] = pd.Categorical(df["condition"], categories=CONDITION_ORDER, ordered=True)
    return df


def merge_behavior_updrs(beh, updrs):
    outer = beh.merge(
        updrs,
        on=["subject", "condition"],
        how="outer",
        suffixes=("_beh", "_updrs"),
        indicator=True
    )

    outer["target_match"] = np.where(
        outer["_merge"] == "both",
        outer["Target_beh"].astype(str) == outer["Target_updrs"].astype(str),
        np.nan
    )

    merge_report = outer[[
        "subject", "condition", "Target_beh", "Target_updrs", "UPDRS", "_merge", "target_match"
    ]].copy()

    both = outer[outer["_merge"] == "both"].copy()
    both["Target"] = both["Target_beh"]
    target_mismatch = both[both["Target_beh"].astype(str) != both["Target_updrs"].astype(str)].copy()

    use_cols = [c for c in beh.columns if c not in ["Target"]]
    merged = both[use_cols + ["Target", "UPDRS"]].copy()
    merged["Target"] = pd.Categorical(merged["Target"], categories=TARGET_ORDER, ordered=True)
    merged["condition"] = pd.Categorical(merged["condition"], categories=CONDITION_ORDER, ordered=True)
    return merged, merge_report, target_mismatch


def apply_qc(df, qc):
    out = df.copy()
    out["qc_keep"] = (
        (out["responded_total"] >= qc["responded_total"]) &
        (out["old_sdt_total"] >= qc["old_sdt_total"]) &
        (out["new_sdt_total"] >= qc["new_sdt_total"])
    )
    session_counts = out.groupby("subject")["qc_keep"].sum()
    valid_subjects = session_counts[session_counts >= 2].index
    out["subject_keep"] = out["subject"].isin(valid_subjects)
    out["analysis_keep"] = out["qc_keep"] & out["subject_keep"]
    return out


def add_updrs_components(df):
    out = df.copy()
    subj_mean = out.groupby("subject")["UPDRS"].transform("mean")
    grand_mean = subj_mean.dropna().mean()
    out["updrs_subject_mean"] = subj_mean
    out["updrs_between_c"] = subj_mean - grand_mean
    out["updrs_within"] = out["UPDRS"] - subj_mean
    return out


def fit_mixedlm(data, formula, group_col="subject"):
    methods = ["lbfgs", "powell", "cg"]
    last_err = None
    for method in methods:
        try:
            model = smf.mixedlm(formula, data=data, groups=data[group_col])
            result = model.fit(reml=False, method=method, disp=False)
            return result, method
        except Exception as e:
            last_err = e
    raise RuntimeError(f"Model fitting failed: {last_err}")


def data_subject_count(result):
    try:
        return len(result.model.group_labels)
    except Exception:
        return np.nan


def fixed_effect_table(result, outcome, model_name, optimizer):
    ci = result.conf_int()
    fe_names = list(result.fe_params.index)
    rows = []
    for name in fe_names:
        rows.append({
            "outcome": outcome,
            "model": model_name,
            "optimizer": optimizer,
            "term": name,
            "coef": result.fe_params[name],
            "se": result.bse_fe[name],
            "z": result.tvalues[name],
            "p": result.pvalues[name],
            "ci_low": ci.loc[name, 0],
            "ci_high": ci.loc[name, 1],
            "aic": result.aic,
            "bic": result.bic,
            "logLik": result.llf,
            "n_obs": result.nobs,
            "n_subjects": data_subject_count(result)
        })
    return pd.DataFrame(rows)


def lr_compare(result_small, result_large, outcome):
    lr = 2 * (result_large.llf - result_small.llf)
    df_diff = len(result_large.fe_params) - len(result_small.fe_params)
    p = chi2.sf(lr, df_diff) if df_diff > 0 else np.nan
    return pd.DataFrame([{
        "outcome": outcome,
        "model_small": "total_effect",
        "model_large": "updrs_adjusted",
        "lr_chi2": lr,
        "df_diff": df_diff,
        "p": p,
        "ll_small": result_small.llf,
        "ll_large": result_large.llf
    }])


def build_contrast_vector(result, terms):
    names = list(result.fe_params.index)
    vec = np.zeros((1, len(names)))
    for term, weight in terms.items():
        if term in names:
            vec[0, names.index(term)] = weight
    return vec


def contrast_to_row(result, contrast_terms, label, outcome, model_name):
    vec = build_contrast_vector(result, contrast_terms)
    test = result.t_test(vec)
    coef = float(test.effect[0])
    se = float(test.sd[0][0]) if np.ndim(test.sd) == 2 else float(test.sd[0])
    z = float(test.tvalue[0][0]) if np.ndim(test.tvalue) == 2 else float(test.tvalue[0])
    p = float(test.pvalue)
    ci = test.conf_int()[0]
    return {
        "outcome": outcome,
        "model": model_name,
        "contrast": label,
        "estimate": coef,
        "se": se,
        "z": z,
        "p_unc": p,
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
    }


def add_fdr(df, p_col="p_unc", out_col="p_fdr"):
    out = df.copy()
    mask = out[p_col].notna()
    out[out_col] = np.nan
    if mask.sum() > 0:
        out.loc[mask, out_col] = multipletests(out.loc[mask, p_col], method="fdr_bh")[1]
    return out


def planned_contrasts(result, outcome, model_name):
    t_sn = 'C(Target, Treatment(reference="STN"))[T.SN]'
    c10 = 'C(condition, Treatment(reference="off"))[T.10Hz]'
    c130 = 'C(condition, Treatment(reference="off"))[T.130Hz]'
    i10 = 'C(Target, Treatment(reference="STN"))[T.SN]:C(condition, Treatment(reference="off"))[T.10Hz]'
    i130 = 'C(Target, Treatment(reference="STN"))[T.SN]:C(condition, Treatment(reference="off"))[T.130Hz]'

    rows = []
    rows.append(contrast_to_row(result, {c10: 1}, "STN: 10Hz - off", outcome, model_name))
    rows.append(contrast_to_row(result, {c130: 1}, "STN: 130Hz - off", outcome, model_name))
    rows.append(contrast_to_row(result, {c10: 1, i10: 1}, "SN: 10Hz - off", outcome, model_name))
    rows.append(contrast_to_row(result, {c130: 1, i130: 1}, "SN: 130Hz - off", outcome, model_name))
    rows.append(contrast_to_row(result, {t_sn: 1}, "SN - STN at off", outcome, model_name))
    rows.append(contrast_to_row(result, {t_sn: 1, i10: 1}, "SN - STN at 10Hz", outcome, model_name))
    rows.append(contrast_to_row(result, {t_sn: 1, i130: 1}, "SN - STN at 130Hz", outcome, model_name))
    rows.append(contrast_to_row(result, {i10: 1}, "Interaction difference at 10Hz", outcome, model_name))
    rows.append(contrast_to_row(result, {i130: 1}, "Interaction difference at 130Hz", outcome, model_name))

    out = pd.DataFrame(rows)
    out = add_fdr(out, p_col="p_unc", out_col="p_fdr")
    return out


def cell_predictions(result, outcome, model_name, adjusted=False):
    rows = []
    cells = [
        ("STN", "off"),
        ("STN", "10Hz"),
        ("STN", "130Hz"),
        ("SN", "off"),
        ("SN", "10Hz"),
        ("SN", "130Hz"),
    ]
    for target, condition in cells:
        new = pd.DataFrame({
            "Target": pd.Categorical([target], categories=TARGET_ORDER, ordered=True),
            "condition": pd.Categorical([condition], categories=CONDITION_ORDER, ordered=True),
        })
        if adjusted:
            new["updrs_between_c"] = [0.0]
            new["updrs_within"] = [0.0]
        pred = float(result.predict(new)[0])

        if target == "STN" and condition == "off":
            terms = {"Intercept": 1}
        elif target == "STN" and condition == "10Hz":
            terms = {"Intercept": 1, 'C(condition, Treatment(reference="off"))[T.10Hz]': 1}
        elif target == "STN" and condition == "130Hz":
            terms = {"Intercept": 1, 'C(condition, Treatment(reference="off"))[T.130Hz]': 1}
        elif target == "SN" and condition == "off":
            terms = {"Intercept": 1, 'C(Target, Treatment(reference="STN"))[T.SN]': 1}
        elif target == "SN" and condition == "10Hz":
            terms = {
                "Intercept": 1,
                'C(Target, Treatment(reference="STN"))[T.SN]': 1,
                'C(condition, Treatment(reference="off"))[T.10Hz]': 1,
                'C(Target, Treatment(reference="STN"))[T.SN]:C(condition, Treatment(reference="off"))[T.10Hz]': 1,
            }
        else:
            terms = {
                "Intercept": 1,
                'C(Target, Treatment(reference="STN"))[T.SN]': 1,
                'C(condition, Treatment(reference="off"))[T.130Hz]': 1,
                'C(Target, Treatment(reference="STN"))[T.SN]:C(condition, Treatment(reference="off"))[T.130Hz]': 1,
            }

        row = contrast_to_row(result, terms, f"{target}-{condition}", outcome, model_name)
        row["predicted_mean"] = pred
        rows.append(row)
    return pd.DataFrame(rows)


def run_outcome_models(df, outcome):
    use = df[df["analysis_keep"]].copy()
    use = use.dropna(subset=[outcome, "UPDRS"])
    use["Target"] = pd.Categorical(use["Target"], categories=TARGET_ORDER, ordered=True)
    use["condition"] = pd.Categorical(use["condition"], categories=CONDITION_ORDER, ordered=True)

    formula_base = (
        f"{outcome} ~ "
        f'C(Target, Treatment(reference="STN")) * C(condition, Treatment(reference="off"))'
    )
    formula_adj = (
        f"{outcome} ~ "
        f'C(Target, Treatment(reference="STN")) * C(condition, Treatment(reference="off"))'
        f" + updrs_between_c + updrs_within"
    )

    res_total, opt_total = fit_mixedlm(use, formula_base)
    res_adj, opt_adj = fit_mixedlm(use, formula_adj)

    fixed_total = fixed_effect_table(res_total, outcome, "total_effect", opt_total)
    fixed_adj = fixed_effect_table(res_adj, outcome, "updrs_adjusted", opt_adj)
    contrasts_total = planned_contrasts(res_total, outcome, "total_effect")
    contrasts_adj = planned_contrasts(res_adj, outcome, "updrs_adjusted")
    cells_total = cell_predictions(res_total, outcome, "total_effect", adjusted=False)
    cells_adj = cell_predictions(res_adj, outcome, "updrs_adjusted", adjusted=True)
    lr = lr_compare(res_total, res_adj, outcome)

    fit_info = pd.DataFrame([{
        "outcome": outcome,
        "n_obs": len(use),
        "n_subjects": use["subject"].nunique(),
        "optimizer_total": opt_total,
        "optimizer_adjusted": opt_adj,
        "aic_total": res_total.aic,
        "aic_adjusted": res_adj.aic,
        "bic_total": res_total.bic,
        "bic_adjusted": res_adj.bic,
        "ll_total": res_total.llf,
        "ll_adjusted": res_adj.llf
    }])

    return {
        "data_used": use,
        "fixed_total": fixed_total,
        "fixed_adjusted": fixed_adj,
        "contrasts_total": contrasts_total,
        "contrasts_adjusted": contrasts_adj,
        "cells_total": cells_total,
        "cells_adjusted": cells_adj,
        "fit_info": fit_info,
        "lr_test": lr,
        "res_total": res_total,
        "res_adjusted": res_adj,
    }


def plot_observed_primary(df, out_png):
    use = df[df["analysis_keep"]].copy()
    use = use.dropna(subset=[PRIMARY_OUTCOME])

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    x = np.arange(len(CONDITION_ORDER))
    styles = {
        "STN": {"color": "#4C72B0", "marker": "o", "linestyle": "--", "offset": -0.05},
        "SN": {"color": "#DD8452", "marker": "s", "linestyle": "-", "offset": 0.05},
    }

    for target in TARGET_ORDER:
        subdf = use[use["Target"] == target].copy()
        style = styles[target]

        for subject, g in subdf.groupby("subject"):
            g = g.sort_values("condition")
            if len(g) >= 2:
                xs = np.array([CONDITION_ORDER.index(c) for c in g["condition"]]) + style["offset"]
                ys = g[PRIMARY_OUTCOME].values
                ax.plot(xs, ys, color=style["color"], alpha=0.18, linewidth=0.9, zorder=1)
                ax.scatter(xs, ys, color=style["color"], alpha=0.20, s=16, zorder=2)

        mean_df = subdf.groupby("condition")[PRIMARY_OUTCOME].agg(["mean", "std", "count"]).reindex(CONDITION_ORDER)
        means = mean_df["mean"].values
        sems = mean_df["std"].values / np.sqrt(mean_df["count"].values)
        xs = x + style["offset"]
        ax.errorbar(xs, means, yerr=sems, color=style["color"], marker=style["marker"],
                    linestyle=style["linestyle"], linewidth=2.2, markersize=7, capsize=4,
                    label=f"{target} (n={subdf['subject'].nunique()})", zorder=5)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Observed d′")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_adjusted_primary(cells_df, out_png):
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    x = np.arange(len(CONDITION_ORDER))
    styles = {
        "STN": {"color": "#4C72B0", "marker": "o", "linestyle": "--", "offset": -0.05},
        "SN": {"color": "#DD8452", "marker": "s", "linestyle": "-", "offset": 0.05},
    }

    for target in TARGET_ORDER:
        sub = cells_df[cells_df["contrast"].str.startswith(target)].copy()
        sub["condition"] = sub["contrast"].str.split("-").str[1]
        sub = sub.set_index("condition").reindex(CONDITION_ORDER).reset_index()
        xs = x + styles[target]["offset"]
        means = sub["predicted_mean"].values
        lower = sub["ci_low"].values
        upper = sub["ci_high"].values
        yerr = np.vstack([means - lower, upper - means])
        ax.errorbar(xs, means, yerr=yerr, color=styles[target]["color"],
                    marker=styles[target]["marker"], linestyle=styles[target]["linestyle"],
                    linewidth=2.2, markersize=7, capsize=4, label=target)

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITION_ORDER)
    ax.set_xlabel("Condition")
    ax.set_ylabel("Adjusted d′ (UPDRS centered)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summarize_qc(df, label):
    return pd.DataFrame([{
        "qc_label": label,
        "n_rows_total": len(df),
        "n_rows_keep": int(df["analysis_keep"].sum()),
        "n_subjects_total": df["subject"].nunique(),
        "n_subjects_keep": df.loc[df["analysis_keep"], "subject"].nunique(),
        "mean_responded_total_keep": df.loc[df["analysis_keep"], "responded_total"].mean(),
        "mean_no_response_keep": df.loc[df["analysis_keep"], "no_response_total"].mean(),
    }])


def sensitivity_analysis(merged):
    rows = []
    for label, qc in [("main_qc", MAIN_QC), ("sensitivity_qc", SENS_QC)]:
        qcdf = add_updrs_components(apply_qc(merged, qc))
        use = qcdf[qcdf["analysis_keep"]].dropna(subset=[PRIMARY_OUTCOME, "UPDRS"]).copy()
        if len(use) == 0 or use["subject"].nunique() < 4:
            continue

        formula_base = (
            f"{PRIMARY_OUTCOME} ~ "
            f'C(Target, Treatment(reference="STN")) * C(condition, Treatment(reference="off"))'
        )
        formula_adj = formula_base + " + updrs_between_c + updrs_within"
        res_total, _ = fit_mixedlm(use, formula_base)
        res_adj, _ = fit_mixedlm(use, formula_adj)

        for model_name, res in [("total_effect", res_total), ("updrs_adjusted", res_adj)]:
            for term in res.fe_params.index:
                rows.append({
                    "qc_label": label,
                    "model": model_name,
                    "term": term,
                    "coef": res.fe_params[term],
                    "se": res.bse_fe[term],
                    "z": res.tvalues[term],
                    "p": res.pvalues[term],
                    "n_obs": len(use),
                    "n_subjects": use["subject"].nunique(),
                })
    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run behavioral mixed-effects analyses with UPDRS-III adjustment "
            "and predefined quality-control thresholds."
        )
    )
    parser.add_argument(
        "behavior_file",
        type=Path,
        help="Excel workbook containing the behavioral summary sheet.",
    )
    parser.add_argument(
        "updrs_file",
        type=Path,
        help="Excel workbook containing subject, condition, Target, and UPDRS columns.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Directory for result tables and figures. By default, results are written "
            "to 'updrs_mixed_model_results_fdr' beside the behavioral workbook."
        ),
    )
    return parser.parse_args()


def main():
    args = parse_args()
    beh_file = args.behavior_file.expanduser().resolve()
    updrs_file = args.updrs_file.expanduser().resolve()

    if not beh_file.is_file():
        raise FileNotFoundError(f"Behavioral workbook not found: {beh_file.name}")
    if not updrs_file.is_file():
        raise FileNotFoundError(f"UPDRS workbook not found: {updrs_file.name}")

    out_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else beh_file.parent / "updrs_mixed_model_results_fdr"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    beh = prepare_behavior(beh_file)
    updrs = prepare_updrs(updrs_file)
    merged, merge_report, target_mismatch = merge_behavior_updrs(beh, updrs)

    main_df = add_updrs_components(apply_qc(merged, MAIN_QC))
    sens_df = add_updrs_components(apply_qc(merged, SENS_QC))

    all_outcomes = [PRIMARY_OUTCOME] + SECONDARY_OUTCOMES

    fixed_total_all = []
    fixed_adj_all = []
    contrasts_total_all = []
    contrasts_adj_all = []
    cells_total_all = []
    cells_adj_all = []
    fit_info_all = []
    lr_all = []

    model_objects = {}

    for outcome in all_outcomes:
        results = run_outcome_models(main_df, outcome)
        model_objects[outcome] = results
        fixed_total_all.append(results["fixed_total"])
        fixed_adj_all.append(results["fixed_adjusted"])
        contrasts_total_all.append(results["contrasts_total"])
        contrasts_adj_all.append(results["contrasts_adjusted"])
        cells_total_all.append(results["cells_total"])
        cells_adj_all.append(results["cells_adjusted"])
        fit_info_all.append(results["fit_info"])
        lr_all.append(results["lr_test"])

    fixed_total_df = pd.concat(fixed_total_all, ignore_index=True)
    fixed_adj_df = pd.concat(fixed_adj_all, ignore_index=True)
    contrasts_total_df = pd.concat(contrasts_total_all, ignore_index=True)
    contrasts_adj_df = pd.concat(contrasts_adj_all, ignore_index=True)
    cells_total_df = pd.concat(cells_total_all, ignore_index=True)
    cells_adj_df = pd.concat(cells_adj_all, ignore_index=True)
    fit_info_df = pd.concat(fit_info_all, ignore_index=True)
    lr_df = pd.concat(lr_all, ignore_index=True)

    qc_summary = pd.concat([
        summarize_qc(main_df, "main_qc_24_8_8"),
        summarize_qc(sens_df, "sensitivity_qc_36_12_12"),
    ], ignore_index=True)

    sensitivity_df = sensitivity_analysis(merged)

    observed_plot = out_dir / "dprime_observed_main_qc.png"
    adjusted_plot = out_dir / "dprime_adjusted_main_qc.png"
    plot_observed_primary(main_df, observed_plot)
    plot_adjusted_primary(cells_adj_df[cells_adj_df["outcome"] == PRIMARY_OUTCOME], adjusted_plot)

    results_xlsx = out_dir / "updrs_mixed_model_results_fdr.xlsx"
    with pd.ExcelWriter(results_xlsx, engine="openpyxl") as writer:
        merge_report.to_excel(writer, sheet_name="merge_report", index=False)
        target_mismatch.to_excel(writer, sheet_name="target_mismatch", index=False)
        main_df.to_excel(writer, sheet_name="main_qc_data", index=False)
        sens_df.to_excel(writer, sheet_name="sensitivity_qc_data", index=False)
        qc_summary.to_excel(writer, sheet_name="qc_summary", index=False)
        fixed_total_df.to_excel(writer, sheet_name="fixed_total", index=False)
        fixed_adj_df.to_excel(writer, sheet_name="fixed_adjusted", index=False)
        contrasts_total_df.to_excel(writer, sheet_name="contrasts_total", index=False)
        contrasts_adj_df.to_excel(writer, sheet_name="contrasts_adjusted", index=False)
        cells_total_df.to_excel(writer, sheet_name="cells_total", index=False)
        cells_adj_df.to_excel(writer, sheet_name="cells_adjusted", index=False)
        fit_info_df.to_excel(writer, sheet_name="fit_info", index=False)
        lr_df.to_excel(writer, sheet_name="lr_compare", index=False)
        sensitivity_df.to_excel(writer, sheet_name="sensitivity_models", index=False)

        notes = pd.DataFrame({
            "item": [
                "main_qc",
                "sensitivity_qc",
                "primary_outcome",
                "secondary_outcomes",
                "updrs_between_c",
                "updrs_within",
                "total_effect_model",
                "updrs_adjusted_model",
                "multiple_comparison_correction",
                "important_note",
            ],
            "value": [
                "responded_total >= 24, old_sdt_total >= 8, new_sdt_total >= 8; subject requires >=2 usable sessions",
                "responded_total >= 36, old_sdt_total >= 12, new_sdt_total >= 12; subject requires >=2 usable sessions",
                "d_prime",
                "pr, hit_rate, false_alarm_rate, no_response_rate, log_rt_hit_mean, log_rt_correct_rejection_mean",
                "subject mean UPDRS minus grand mean, representing between-person motor severity",
                "session UPDRS minus subject mean UPDRS, representing within-person motor fluctuation",
                "outcome ~ Target * condition + (1|subject)",
                "outcome ~ Target * condition + updrs_between_c + updrs_within + (1|subject)",
                "Planned contrasts are corrected with Benjamini-Hochberg FDR (p_fdr). Fixed effects in the main mixed model are not FDR-corrected.",
                "For the adjusted model, cell predictions are evaluated at average UPDRS (updrs_between_c=0, updrs_within=0).",
            ]
        })
        notes.to_excel(writer, sheet_name="notes", index=False)

    summary_txt = out_dir / "analysis_summary_with_updrs_fdr.txt"
    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("UPDRS-adjusted mixed model analysis finished (FDR-corrected planned contrasts).\n\n")
        f.write(f"Behavior file: {beh_file.name}\n")
        f.write(f"UPDRS file: {updrs_file.name}\n")
        f.write(f"Output directory: {out_dir.name}\n\n")
        f.write("QC thresholds (main): responded_total>=24, old_sdt_total>=8, new_sdt_total>=8\n")
        f.write("QC thresholds (sensitivity): responded_total>=36, old_sdt_total>=12, new_sdt_total>=12\n")
        f.write("Planned contrasts corrected with Benjamini-Hochberg FDR (p_fdr).\n\n")
        f.write("Merge report:\n")
        f.write(f"- rows in behavior summary: {len(beh)}\n")
        f.write(f"- rows in UPDRS file: {len(updrs)}\n")
        f.write(f"- merged rows for analysis candidate: {len(merged)}\n")
        f.write(f"- target mismatches after merge: {len(target_mismatch)}\n\n")

        for outcome in [PRIMARY_OUTCOME] + SECONDARY_OUTCOMES:
            f.write("=" * 72 + "\n")
            f.write(f"Outcome: {outcome}\n\n")
            fit_sub = fit_info_df[fit_info_df["outcome"] == outcome]
            f.write("[Fit info]\n")
            f.write(fit_sub.to_string(index=False))
            f.write("\n\n[Total effect fixed effects]\n")
            f.write(fixed_total_df[fixed_total_df["outcome"] == outcome].to_string(index=False))
            f.write("\n\n[UPDRS-adjusted fixed effects]\n")
            f.write(fixed_adj_df[fixed_adj_df["outcome"] == outcome].to_string(index=False))
            f.write("\n\n[Adjusted planned contrasts with FDR]\n")
            f.write(contrasts_adj_df[contrasts_adj_df["outcome"] == outcome].to_string(index=False))
            f.write("\n\n")

    print(f"Results workbook: {results_xlsx.name}")
    print(f"Analysis summary: {summary_txt.name}")
    print(f"Figures: {observed_plot.name}, {adjusted_plot.name}")


if __name__ == "__main__":
    main()
