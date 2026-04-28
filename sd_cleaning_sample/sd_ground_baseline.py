# -*- coding: utf-8 -*-
"""Flat-segment baseline correction: ``detect_flat_offsets``.

This repository is a **minimal sample** of the full snow-depth processing chain: it
isolates **one** step—piecewise baseline (ground-level) adjustment in low-variance
segments. Other pipeline stages (spike removal, physical limits, multi-station I/O,
etc.) live in the main project; they are **not** included here.

The module is meant to be imported from the full pipeline or used with
``baseline_viz.ipynb`` and the bundled example CSV ``SD_ValleOlivares_05706003.csv``.
"""

from __future__ import annotations

import argparse
import os
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if not os.environ.get("DISPLAY", "").strip():
    matplotlib.use("Agg")


def detect_flat_offsets(
    data: pd.DataFrame,
    cols: List[str],
    *,
    window: int = 7,
    var_thr: float = 1,
    tol: float = 0.5,
    drift_thr: float = 0.5,
    snow_thr: float = 5.0,
    backfill_first: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, pd.Series]]:
    """Detect and correct step-like baseline offsets in low-variance segments.

    Uses a rolling min–max range; stable runs whose median (unless snow is likely)
    become subtracted piecewise-constant offsets.
    """
    corr = data.copy()
    offs: Dict[str, pd.Series] = {}

    for col in cols:
        s = corr[col].astype(float)
        roll_rng = s.rolling(window, min_periods=window).apply(
            lambda x: x.max() - x.min(), raw=True
        )
        stable = roll_rng.le(var_thr).fillna(False)
        seg_id = stable.ne(stable.shift()).cumsum()

        off = pd.Series(0.0, index=s.index, name=f"flat_{col}")
        current = 0.0
        first_applied = False

        for seg in seg_id[stable].unique():
            mask = (seg_id == seg) & stable
            if mask.sum() < window:
                continue

            est = s[mask].median()
            if abs(est) > snow_thr:
                continue
            if abs(est) <= tol:
                est = 0.0
            elif abs(est - current) < drift_thr:
                continue

            current = est
            first_idx = mask.idxmax()
            off.loc[first_idx:] = current
            if backfill_first and not first_applied:
                off.loc[:first_idx] = current
                first_applied = True

        corr[col] = s - off
        offs[col] = off

    return corr, offs


def _apply_ylim(
    ax: matplotlib.axes.Axes,
    ylim: Optional[Tuple[Optional[float], Optional[float]]],
) -> None:
    """Set y-limits; use ``None`` for one end to keep autoscaling there."""
    if ylim is None:
        return
    lo, hi = ylim
    if lo is not None and hi is not None:
        ax.set_ylim(float(lo), float(hi), auto=False)
    elif lo is not None:
        ax.set_ylim(bottom=float(lo), auto=False)
    elif hi is not None:
        ax.set_ylim(top=float(hi), auto=False)


def trim_dataframe_to_first_valid_value(
    d: pd.DataFrame, col: str
) -> pd.DataFrame:
    """Drop leading rows where ``col`` is null."""
    if col not in d.columns:
        return d
    first = d[col].first_valid_index()
    if first is None:
        return d
    return d.loc[first:].copy().reset_index(drop=True)


def filter_by_date_range(
    d: pd.DataFrame,
    date_col: str,
    start: Optional[Union[str, pd.Timestamp]] = None,
    end: Optional[Union[str, pd.Timestamp]] = None,
) -> pd.DataFrame:
    """Keep rows with ``date_col`` in the inclusive [start, end] window (date-level)."""
    if start is None and end is None:
        return d
    ts = d[date_col]
    m = ts.notna()
    if start is not None:
        m = m & (ts >= pd.Timestamp(start))
    if end is not None:
        m = m & (ts < pd.Timestamp(end) + pd.Timedelta(days=1))
    out = d.loc[m].copy().reset_index(drop=True)
    if out.empty:
        raise ValueError("No rows left after date range filter; widen start/end.")
    return out


def load_and_prepare_raw_table(
    csv_path: str,
    date_col: str = "date",
) -> pd.DataFrame:
    """Read CSV and coerce numeric value columns (same rules as visualization)."""
    raw = pd.read_csv(csv_path, encoding="utf-8-sig")
    if date_col not in raw.columns:
        raise ValueError(f"Missing date column: {date_col!r}")
    raw[date_col] = pd.to_datetime(raw[date_col], errors="coerce")
    num_cols = [
        c
        for c in raw.select_dtypes(include=[np.number]).columns
        if c != date_col
    ]
    for c in raw.columns:
        if c == date_col:
            continue
        if c not in num_cols and raw[c].dtype == object:
            raw[c] = pd.to_numeric(raw[c], errors="coerce")
    return raw


def _resolve_explicit_column(
    raw: pd.DataFrame,
    date_col: str,
    col: str,
) -> str:
    """Map user ``col`` to a real index label (string vs int station ids, e.g. ``'06020001'`` vs ``6020001``)."""
    if col in raw.columns:
        if col == date_col:
            raise ValueError("Value column cannot be the date column.")
        return col
    # Match digit-only names by integer (pandas may use int 6020001; CSV may say 06020001)
    try:
        want = int(str(col).strip())
    except (TypeError, ValueError):
        want = None
    if want is not None:
        for c in raw.columns:
            if c == date_col:
                continue
            try:
                if int(str(c).strip()) == want:
                    return c
            except (TypeError, ValueError):
                continue
    avail = [c for c in raw.columns if c != date_col]
    raise ValueError(
        f"Column {col!r} not in CSV. Available value columns: {avail}. "
        "Use col=None to pick the only numeric series, or pass the exact column name from the file."
    )


def infer_value_column(
    raw: pd.DataFrame,
    date_col: str = "date",
    col: Optional[str] = None,
) -> str:
    """Return the series column name: ``col`` if set, else the only numeric column besides ``date_col``."""
    num_cols = [
        c
        for c in raw.columns
        if c != date_col and pd.api.types.is_numeric_dtype(raw[c])
    ]
    if not num_cols:
        raise ValueError("No numeric value column besides date.")
    if col is not None:
        return _resolve_explicit_column(raw, date_col, col)
    if len(num_cols) > 1:
        raise ValueError(
            "Multiple numeric columns "
            f"{num_cols}; pass col= or use a single-station CSV."
        )
    return num_cols[0]


def visualize_baseline_from_csv(
    csv_path: str,
    *,
    col: Optional[str] = None,
    date_col: str = "date",
    window: int = 7,
    var_thr: float = 1.0,
    tol: float = 0.5,
    drift_thr: float = 0.5,
    snow_thr: float = 5.0,
    backfill_first: bool = True,
    title: Optional[str] = None,
    savefig: Optional[str] = None,
    y_lim_snow: Optional[Tuple[Optional[float], Optional[float]]] = None,
    y_lim_offset: Optional[Tuple[Optional[float], Optional[float]]] = None,
    trim_to_first_valid: bool = True,
    date_start: Optional[Union[str, pd.Timestamp]] = None,
    date_end: Optional[Union[str, pd.Timestamp]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, pd.Series]]:
    """Load a CSV, run ``detect_flat_offsets`` only, plot raw / corrected and offset.

    - ``y_lim_snow`` / ``y_lim_offset``: exact axis limits in cm; ``(None, None)`` to autoscale.
    - ``trim_to_first_valid``: drop leading nulls in the value column.
    - ``date_start`` / ``date_end``: optional inclusive window to subset rows before
      correction (better zoom for inspection; algorithm runs on the subset only).
    - If ``col`` is omitted and the file has exactly one numeric column (besides ``date_col``),
      that column is used automatically.
    """
    raw = load_and_prepare_raw_table(csv_path, date_col=date_col)
    col = infer_value_column(raw, date_col=date_col, col=col)

    d = raw[[date_col, col]].copy().sort_values(date_col).reset_index(drop=True)
    if trim_to_first_valid:
        d = trim_dataframe_to_first_valid_value(d, col)
    d = filter_by_date_range(d, date_col, date_start, date_end)

    t = d[date_col]
    work = d[[col]].copy()

    kwargs: Dict[str, Any] = dict(
        window=window,
        var_thr=var_thr,
        tol=tol,
        drift_thr=drift_thr,
        snow_thr=snow_thr,
        backfill_first=backfill_first,
    )
    corrected, offsets = detect_flat_offsets(work, [col], **kwargs)
    off_series = offsets[col]
    s_raw = d[col].astype(float)
    s_corr = corrected[col].astype(float)

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True, height_ratios=[2, 1])
    ax0, ax1 = axes

    ax0.plot(
        t,
        s_raw,
        ".-",
        alpha=0.45,
        color="gray",
        label="Input (raw)",
        markersize=2,
    )
    ax0.plot(
        t,
        s_corr,
        "-",
        color="steelblue",
        linewidth=1.2,
        label="After detect_flat_offsets",
    )
    ax0.set_ylabel("Snow depth (cm)")
    range_note = ""
    if date_start is not None or date_end is not None:
        range_note = f" (subset: {date_start or '…'} to {date_end or '…'})"
    ax0.set_title(
        title
        or f"Baseline correction — {os.path.basename(csv_path)} — {col}{range_note}"
    )
    ax0.legend(loc="upper right")
    ax0.grid(True, linestyle="--", alpha=0.5)
    _y_snow_any = y_lim_snow is not None and (
        y_lim_snow[0] is not None or y_lim_snow[1] is not None
    )
    ax0.margins(
        x=0.01,
        y=0.0 if _y_snow_any else 0.05,
    )
    _apply_ylim(ax0, y_lim_snow)

    ax1.fill_between(
        t,
        0.0,
        off_series,
        color="coral",
        alpha=0.4,
        step=None,
    )
    ax1.plot(
        t,
        off_series,
        color="darkred",
        linewidth=0.8,
        label="Cumulative offset removed",
    )
    ax1.set_ylabel("Offset (cm)")
    ax1.set_xlabel("Date")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    _y_off_any = y_lim_offset is not None and (
        y_lim_offset[0] is not None or y_lim_offset[1] is not None
    )
    ax1.margins(y=0.0 if _y_off_any else 0.05)
    _apply_ylim(ax1, y_lim_offset)
    fig.autofmt_xdate(rotation=30)
    fig.tight_layout()

    if savefig:
        fig.savefig(savefig, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    default_csv = os.path.join(here, "SD_ValleOlivares_05706003.csv")

    p = argparse.ArgumentParser(
        description="Sample CLI: baseline step (detect_flat_offsets) from a demo CSV.",
    )
    p.add_argument(
        "csv",
        nargs="?",
        default=default_csv,
        help=f"Path to CSV (default: {default_csv})",
    )
    p.add_argument(
        "--col",
        default=None,
        help="Value column (e.g. 05703014). Inferred if there is a single numeric column.",
    )
    p.add_argument("--window", type=int, default=7)
    p.add_argument("--var-thr", type=float, default=1.0)
    p.add_argument("--tol", type=float, default=0.5)
    p.add_argument("--drift-thr", type=float, default=0.5)
    p.add_argument("--snow-thr", type=float, default=5.0)
    p.add_argument(
        "--no-backfill-first",
        action="store_true",
        help="Use backfill_first=False in detect_flat_offsets",
    )
    p.add_argument(
        "-o",
        "--output",
        default=None,
        help="Save figure to this file (e.g. headless)",
    )
    p.add_argument(
        "--y-snow",
        metavar=("MIN", "MAX"),
        nargs=2,
        type=float,
        default=None,
        help="Y axis limits for depth panel (cm), e.g. 0 25",
    )
    p.add_argument(
        "--y-offset",
        metavar=("MIN", "MAX"),
        nargs=2,
        type=float,
        default=None,
        help="Y axis limits for offset panel (cm)",
    )
    p.add_argument(
        "--no-trim-start",
        action="store_true",
        help="Keep leading rows with no value in the series column",
    )
    p.add_argument(
        "--date-from",
        default=None,
        help="Start date (inclusive), e.g. 2018-04-01; subset data before processing",
    )
    p.add_argument(
        "--date-to",
        default=None,
        help="End date (inclusive), e.g. 2019-03-31",
    )
    args = p.parse_args()

    if not os.path.isfile(args.csv):
        raise SystemExit(f"File not found: {args.csv}")

    y_s = tuple(args.y_snow) if args.y_snow is not None else None
    y_o = tuple(args.y_offset) if args.y_offset is not None else None

    visualize_baseline_from_csv(
        args.csv,
        col=args.col,
        window=args.window,
        var_thr=args.var_thr,
        tol=args.tol,
        drift_thr=args.drift_thr,
        snow_thr=args.snow_thr,
        backfill_first=not args.no_backfill_first,
        savefig=args.output,
        y_lim_snow=y_s,
        y_lim_offset=y_o,
        trim_to_first_valid=not args.no_trim_start,
        date_start=args.date_from,
        date_end=args.date_to,
    )


if __name__ == "__main__":
    main()
