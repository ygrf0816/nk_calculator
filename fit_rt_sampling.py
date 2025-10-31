import argparse
from typing import Tuple, List, Optional, Callable

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt


def _find_columns(df: pd.DataFrame) -> Tuple[Optional[str], str, Optional[str], str]:
    """
    返回 (wl_R_col, R_col, wl_T_col, T_col)
    - 若存在单一波长列（wavelength/lambda/wl），则 wl_R_col=wl_T_col=该列
    - 若存在分离列（wl_R, wl_T），则分别返回
    """
    # 先找强匹配：分离列
    wl_r_col = "wl_R" if "wl_R" in df.columns else None
    wl_t_col = "wl_T" if "wl_T" in df.columns else None

    r_candidates: List[str] = ["R", "r", "Reflectance", "reflectance"]
    t_candidates: List[str] = ["T", "t", "Transmittance", "transmittance"]
    R_col = next((c for c in r_candidates if c in df.columns), None)
    T_col = next((c for c in t_candidates if c in df.columns), None)

    # 单一波长列候选
    wl_candidates: List[str] = [
        "wavelength", "lambda", "wl", "Wavelength", "Lambda"
    ]
    wl_single = next((c for c in wl_candidates if c in df.columns), None)

    if R_col is None or T_col is None:
        raise ValueError(f"CSV需包含 R 与 T 列。实际列: {list(df.columns)}")

    if wl_r_col is None and wl_t_col is None:
        if wl_single is None:
            raise ValueError(
                f"未找到波长列。需包含 wl_R/wl_T 或通用 wavelength/lambda/wl。实际列: {list(df.columns)}"
            )
        return wl_single, R_col, wl_single, T_col

    # 若至少一个分离列存在，缺失的一侧可回退为单一波长列
    if wl_r_col is None:
        if wl_single is None:
            raise ValueError("缺少 wl_R 且无通用波长列可用")
        wl_r_col = wl_single
    if wl_t_col is None:
        if wl_single is None:
            raise ValueError("缺少 wl_T 且无通用波长列可用")
        wl_t_col = wl_single

    return wl_r_col, R_col, wl_t_col, T_col


def _fit_spline(x: np.ndarray, y: np.ndarray, s_factor: float = 0.002) -> Callable[[np.ndarray], np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    # 排序
    sorter = np.argsort(x)
    x = x[sorter]
    y = y[sorter]
    # 去重：对重复x取均值，确保严格递增
    unique_x, idx_starts = np.unique(x, return_index=True)
    if len(unique_x) != len(x):
        # 对每个唯一x求平均y
        means = []
        for ux in unique_x:
            mask = (x == ux)
            means.append(float(np.mean(y[mask])))
        x = unique_x
        y = np.asarray(means, dtype=float)
    # 至少需要两个点
    if x.size < 2:
        raise ValueError("数据点不足以拟合样条（至少需要2个不同波长点）")
    # 使用单调PCHIP，避免振荡，支持外推
    interpolator = PchipInterpolator(x, y, extrapolate=True)
    return interpolator


def process_file(csv_path: str, out_csv: str, out_png: str,
                 wl_min_nm: float = 220.0, wl_max_nm: float = 2200.0, step_nm: float = 5.0,
                 s_factor: float = 0.002) -> None:
    df = pd.read_csv(csv_path)
    wl_r_col, r_col, wl_t_col, t_col = _find_columns(df)

    # R通道
    r_df = df[[wl_r_col, r_col]].dropna()
    wl_R = r_df[wl_r_col].to_numpy(dtype=float)
    R = r_df[r_col].to_numpy(dtype=float)

    # T通道
    t_df = df[[wl_t_col, t_col]].dropna()
    wl_T = t_df[wl_t_col].to_numpy(dtype=float)
    T = t_df[t_col].to_numpy(dtype=float)

    # 限制R/T在[0,1]
    R = np.clip(R, 0.0, 1.0)
    T = np.clip(T, 0.0, 1.0)

    # 拟合样条
    spline_R = _fit_spline(wl_R, R, s_factor=s_factor)
    spline_T = _fit_spline(wl_T, T, s_factor=s_factor)

    # 均匀采样
    wl_grid = np.arange(wl_min_nm, wl_max_nm + 0.5 * step_nm, step_nm, dtype=float)
    R_fit = np.clip(spline_R(wl_grid), 0.0, 1.0)
    T_fit = np.clip(spline_T(wl_grid), 0.0, 1.0)

    # 保存采样CSV
    out_df = pd.DataFrame({
        "wavelength_nm": wl_grid,
        "R_fit": R_fit,
        "T_fit": T_fit,
    })
    out_df.to_csv(out_csv, index=False)

    # 作图：单一子图，原始散点 + 拟合曲线（R 与 T 同图）
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.scatter(wl_R, R, s=14, alpha=0.6, label="R data", color="tab:blue")
    ax.plot(wl_grid, R_fit, color="tab:blue", label="R fit")
    ax.scatter(wl_T, T, s=14, alpha=0.6, label="T data", color="tab:orange")
    ax.plot(wl_grid, T_fit, color="tab:orange", label="T fit")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("R / T")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2)
    fig.suptitle(f"Fit and Sampling: {csv_path}")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def cli():
    parser = argparse.ArgumentParser(description="拟合R/T曲线并在220–2200 nm按5 nm采样，保存CSV与图像")
    parser.add_argument("--csv", nargs=2, required=True,
                        help="两个输入CSV路径，列需包含 wavelength,R,T 或其常见别名")
    parser.add_argument("--out_prefix", type=str, default="fit_",
                        help="输出文件前缀（将自动附加文件基名）")
    parser.add_argument("--wl_min", type=float, default=220.0)
    parser.add_argument("--wl_max", type=float, default=2200.0)
    parser.add_argument("--step", type=float, default=5.0)
    parser.add_argument("--s_factor", type=float, default=0.002,
                        help="UnivariateSpline 平滑系数比例，越大越平滑")
    args = parser.parse_args()

    for csv_path in args.csv:
        base = csv_path.rsplit("/", 1)[-1].rsplit(".csv", 1)[0]
        out_csv = f"{args.out_prefix}{base}_sampled.csv"
        out_png = f"{args.out_prefix}{base}_fit.png"
        process_file(
            csv_path,
            out_csv,
            out_png,
            wl_min_nm=args.wl_min,
            wl_max_nm=args.wl_max,
            step_nm=args.step,
            s_factor=args.s_factor,
        )


if __name__ == "__main__":
    cli()


