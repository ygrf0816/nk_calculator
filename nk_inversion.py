import argparse
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import logging
import os
from datetime import datetime

def tm_inversion(wl: np.ndarray, R: np.ndarray, T: np.ndarray, t: float, debug_logger=None, debug_csv_path=None):
    # 辅助变量
    a = T**2-(1-R)**2
    b = np.sqrt(a**2+4*T**2)

    #Equation (25) Cal attenuation factor
    eta = (a+b)/2*T

    #Equation (26) Cal absorption coeff
    # alpha = (-1/t)*np.log(eta)

    #Equation (27) Cal imag part
    k = (-wl/(4*math.pi*t))*np.log(eta)

    #Equation (28) Cal Ras
    R_as = 2*R/(2+a+b)

    #Equation (29) Cal real part
    sqrt_arg = (4*R_as/(1-R_as)**2)-k**2
    n = (1+R_as)/(1-R_as) + np.sqrt(sqrt_arg)

    # 调试输出：保存逐点数据
    if debug_logger and debug_csv_path:
        debug_df = pd.DataFrame({
            'wavelength_nm': wl * 1e9,
            'R': R,
            'T': T,
            'a': a,
            'b': b,
            'eta': eta,
            'R_as': R_as,
            'sqrt_arg': sqrt_arg,
            'k': k,
            'n': n,
            'k_isnan': np.isnan(k),
            'k_isinf': np.isinf(k),
            'n_isnan': np.isnan(n),
            'n_isinf': np.isinf(n)
        })
        debug_df.to_csv(debug_csv_path, index=False, float_format='%.8e')
        debug_logger.info(f"Saved detailed debug data to: {debug_csv_path}")
        debug_logger.debug("=== TM Inversion Debug Info ===")
        debug_logger.debug(f"Input: wl (nm) range: [{wl.min()*1e9:.2f}, {wl.max()*1e9:.2f}], R range: [{R.min():.6f}, {R.max():.6f}], T range: [{T.min():.6f}, {T.max():.6f}], t={t*1e9:.2f} nm")
        debug_logger.debug(f"Intermediate variables - a: min={a.min():.6e}, max={a.max():.6e}, mean={a.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - b: min={b.min():.6e}, max={b.max():.6e}, mean={b.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - eta: min={eta.min():.6e}, max={eta.max():.6e}, mean={eta.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - R_as: min={R_as.min():.6e}, max={R_as.max():.6e}, mean={R_as.mean():.6e}")
        debug_logger.debug(f"Results - n: min={n.min():.6e}, max={n.max():.6e}, mean={n.mean():.6e}, nan_count={np.isnan(n).sum()}, inf_count={np.isinf(n).sum()}")
        debug_logger.debug(f"Results - k: min={k.min():.6e}, max={k.max():.6e}, mean={k.mean():.6e}, nan_count={np.isnan(k).sum()}, inf_count={np.isinf(k).sum()}")

    return n,k

def rt_inversion(wl: np.ndarray, R: np.ndarray, T: np.ndarray, t: float, debug_logger=None, debug_csv_path=None):
    # 辅助变量
    a = 2+T**2-(1-R)**2
    b = np.sqrt(a**2-4*(2-R)*R)
    c = a-b

    #Equation (20) Cal Ras
    R_as = c / (2*(2-R))
    
    # 检查分母
    denominator_k = (2*(2-R)*R+b-a)
    
    #Equation (21) Cal imag part
    log_arg = (T*c)/denominator_k
    k = (wl/(4*math.pi*t))*np.log(log_arg)

    #Equation (23) Cal real part
    sqrt_arg = (4*R_as/(1-R_as)**2)-k**2
    n = (1+R_as)/(1-R_as) + np.sqrt(sqrt_arg)

    # 调试输出：保存逐点数据
    if debug_logger and debug_csv_path:
        debug_df = pd.DataFrame({
            'wavelength_nm': wl * 1e9,
            'R': R,
            'T': T,
            'a': a,
            'b': b,
            'c': c,
            'R_as': R_as,
            'denominator_k': denominator_k,
            'log_arg': log_arg,
            'sqrt_arg': sqrt_arg,
            'k': k,
            'n': n,
            'denominator_k_iszero': denominator_k == 0,
            'denominator_k_isneg': denominator_k < 0,
            'log_arg_isneg_or_zero': log_arg <= 0,
            'sqrt_arg_isneg': sqrt_arg < 0,
            'k_isnan': np.isnan(k),
            'k_isinf': np.isinf(k),
            'n_isnan': np.isnan(n),
            'n_isinf': np.isinf(n)
        })
        debug_df.to_csv(debug_csv_path, index=False, float_format='%.8e')
        debug_logger.info(f"Saved detailed debug data to: {debug_csv_path}")
        debug_logger.debug("=== RT Inversion Debug Info ===")
        debug_logger.debug(f"Input: wl (nm) range: [{wl.min()*1e9:.2f}, {wl.max()*1e9:.2f}], R range: [{R.min():.6f}, {R.max():.6f}], T range: [{T.min():.6f}, {T.max():.6f}], t={t*1e9:.2f} nm")
        debug_logger.debug(f"Intermediate variables - a: min={a.min():.6e}, max={a.max():.6e}, mean={a.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - b: min={b.min():.6e}, max={b.max():.6e}, mean={b.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - c=(a-b): min={c.min():.6e}, max={c.max():.6e}, mean={c.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - R_as: min={R_as.min():.6e}, max={R_as.max():.6e}, mean={R_as.mean():.6e}")
        debug_logger.debug(f"Intermediate variables - denominator_k: min={denominator_k.min():.6e}, max={denominator_k.max():.6e}, mean={denominator_k.mean():.6e}, zero_count={(denominator_k==0).sum()}, neg_count={(denominator_k<0).sum()}")
        debug_logger.debug(f"Intermediate variables - log argument (T*w/denominator_k): min={log_arg.min():.6e}, max={log_arg.max():.6e}, neg_count={(log_arg<=0).sum()}")
        debug_logger.debug(f"Intermediate variables - sqrt argument: min={sqrt_arg.min():.6e}, max={sqrt_arg.max():.6e}, neg_count={(sqrt_arg<0).sum()}")
        debug_logger.debug(f"Results - n: min={n.min():.6e}, max={n.max():.6e}, mean={n.mean():.6e}, nan_count={np.isnan(n).sum()}, inf_count={np.isinf(n).sum()}")
        debug_logger.debug(f"Results - k: min={k.min():.6e}, max={k.max():.6e}, mean={k.mean():.6e}, nan_count={np.isnan(k).sum()}, inf_count={np.isinf(k).sum()}")

    return n,k

def roughness_correction(R_m: np.ndarray, T_m: np.ndarray, wl: np.ndarray, sigma: float):
    # Equation (37) Cal R under roughness
    R = R_m*np.exp(-16*math.pi**2*sigma**2/(wl**2))

    # Equation (38) Cal R under roughness
    T = T_m*np.exp(-16*math.pi**2*sigma**2/(wl**2))

    return R,T

def _maybe_convert_wavelengths(wavelengths: np.ndarray, unit: str) -> np.ndarray:
    unit = unit.lower()
    if unit in ["m", "meter", "meters"]:
        return wavelengths
    if unit in ["nm", "nanometer", "nanometers"]:
        return wavelengths * 1e-9
    if unit in ["um", "µm", "micrometer", "micron", "micrometers", "microns"]:
        return wavelengths * 1e-6
    if unit in ["ang", "angstrom", "angströms", "a"]:
        return wavelengths * 1e-10
    raise ValueError(f"Unsupported wavelength unit: {unit}")


def _maybe_convert_length(length_value: float, unit: str) -> float:
    unit = unit.lower()
    if unit in ["m", "meter", "meters"]:
        return length_value
    if unit in ["nm", "nanometer", "nanometers"]:
        return length_value * 1e-9
    if unit in ["um", "µm", "micrometer", "micron", "micrometers", "microns"]:
        return length_value * 1e-6
    if unit in ["mm", "millimeter", "millimeters"]:
        return length_value * 1e-3
    if unit in ["ang", "angstrom", "angströms", "a"]:
        return length_value * 1e-10
    raise ValueError(f"Unsupported length unit: {unit}")

def plot_fig(n,k,wl,figure_path):
    # 作图：单一子图，原始散点 + 拟合曲线（R 与 T 同图）
    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    ax[0].plot(wl, n, color="tab:blue", label="n")
    ax[1].plot(wl, k, color="tab:orange", label="k")
    ax[0].set_xlabel("Wavelength (nm)")
    ax[0].set_ylabel("n")
    ax[1].set_xlabel("Wavelength (nm)")
    ax[1].set_ylabel("k")
    fig.tight_layout()
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

def cli():
    parser = argparse.ArgumentParser(
        description=(
            "Invert n,k from thickness and R/T spectra using single-layer TMM at normal incidence."
        )
    )
    parser.add_argument("--csv", type=str, required=True,
                        help="Input CSV with columns: wavelength,R,T (wavelength unit set by --wl_unit)")
    parser.add_argument("--out", type=str, required=True,
                        help="Output CSV path with columns: wavelength(nm),n,k,status")
    parser.add_argument("--thickness", type=float, required=True, help="Sample thickness value")
    parser.add_argument("--th_unit", type=str, default="nm", help="Thickness unit (nm|um|m|mm|ang)")
    parser.add_argument("--wl_unit", type=str, default="nm", help="Wavelength unit in input CSV (nm|um|m|ang)")
    parser.add_argument("--sigma", type=float, required=False, help="root mean square roughness")
    parser.add_argument("--fig_path", type=str, required=False, help="figure output path")
    parser.add_argument("--method", type=str, required=True, help="calculate method")
    parser.add_argument("--debug_log", type=str, required=False, help="Enable debug logging and save to specified log file path")
    args = parser.parse_args()

    # 读入数据
    df = pd.read_csv(args.csv)
    if not {"wavelength_nm", "R", "T"}.issubset(df.columns):
        raise ValueError("CSV must contain columns: wavelength,R,T")

    wavelengths = df["wavelength_nm"].to_numpy(dtype=float)
    R = df["R"].to_numpy(dtype=float)
    T = df["T"].to_numpy(dtype=float)

    # 验证 R+T <= 1 (能量守恒)
    invalid_mask = (R + T) > 1
    if invalid_mask.any():
        invalid_indices = np.where(invalid_mask)[0]
        invalid_wavelengths = wavelengths[invalid_indices]
        invalid_R = R[invalid_indices]
        invalid_T = T[invalid_indices]
        error_msg = f"错误：发现 {invalid_mask.sum()} 个数据点的 R+T > 1（违反能量守恒定律）：\n"
        error_msg += f"前10个错误数据点：\n"
        for i in range(min(10, len(invalid_indices))):
            error_msg += f"  波长={invalid_wavelengths[i]:.2f} nm, R={invalid_R[i]:.6f}, T={invalid_T[i]:.6f}, R+T={invalid_R[i]+invalid_T[i]:.6f}\n"
        raise ValueError(error_msg)

    # 处理0值：如果T或R中有0，则替换为1e-8以避免计算时出现无穷
    # R = np.where(R == 0, 1e-8, R)
    # T = np.where(T == 0, 1e-8, T)
    R = np.where(R == 0, 1e-3, R)
    T = np.where(T == 0, 1e-3, T)
    # 设置调试日志
    debug_logger = None
    if args.debug_log:
        debug_logger = logging.getLogger('nk_inversion_debug')
        debug_logger.setLevel(logging.DEBUG)
        file_handler = logging.FileHandler(args.debug_log, mode='w', encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        debug_logger.addHandler(file_handler)
        debug_logger.info("="*60)
        debug_logger.info("NK Inversion Debug Log")
        debug_logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        debug_logger.info("="*60)
        debug_logger.info(f"Input file: {args.csv}")
        debug_logger.info(f"Output file: {args.out}")
        debug_logger.info(f"Method: {args.method}")
        debug_logger.info(f"Thickness: {args.thickness} {args.th_unit}")
        debug_logger.info(f"Wavelength unit: {args.wl_unit}")
        if args.sigma:
            debug_logger.info(f"Roughness correction sigma: {args.sigma}")
        debug_logger.info(f"Total data points: {len(wavelengths)}")
        debug_logger.info(f"Wavelength range: [{wavelengths.min():.2f}, {wavelengths.max():.2f}] nm")
        debug_logger.info(f"R range: [{R.min():.6f}, {R.max():.6f}]")
        debug_logger.info(f"T range: [{T.min():.6f}, {T.max():.6f}]")
        debug_logger.info(f"R+T range: [{(R+T).min():.6f}, {(R+T).max():.6f}]")

    # 根据输入参数转换单位
    wavelengths_m = _maybe_convert_wavelengths(wavelengths, args.wl_unit)
    thickness_m = _maybe_convert_length(args.thickness, args.th_unit)
    
    if debug_logger:
        debug_logger.info(f"Converted thickness: {thickness_m*1e9:.6f} nm")
        debug_logger.info(f"Converted wavelength range: [{wavelengths_m.min()*1e9:.6f}, {wavelengths_m.max()*1e9:.6f}] nm")
    
    # 粗糙度修正
    if args.sigma != None:
        R_original = R.copy()
        T_original = T.copy()
        R,T = roughness_correction(R_original,T_original,wavelengths,args.sigma)
        if debug_logger:
            debug_logger.info(f"Applied roughness correction (sigma={args.sigma})")
            debug_logger.info(f"After correction - R range: [{R.min():.6f}, {R.max():.6f}]")
            debug_logger.info(f"After correction - T range: [{T.min():.6f}, {T.max():.6f}]")

    # 准备调试CSV文件路径
    debug_csv_path = None
    if args.debug_log:
        # 生成调试CSV文件路径：将.log扩展名替换为_debug.csv，如果没有.log则添加_debug.csv
        base_name, ext = os.path.splitext(args.debug_log)
        debug_csv_path = base_name + '_debug.csv'

    # TM法计算n,k
    if args.method == "tm":
        n, k = tm_inversion(wavelengths_m, R, T, thickness_m, debug_logger, debug_csv_path)
    elif args.method == "rt":
        n, k = rt_inversion(wavelengths_m, R, T, thickness_m, debug_logger, debug_csv_path)
    
    if debug_logger:
        debug_logger.info("="*60)
        debug_logger.info("Final Results Summary")
        debug_logger.info(f"n - min: {n.min():.6e}, max: {n.max():.6e}, mean: {n.mean():.6e}")
        debug_logger.info(f"n - NaN count: {np.isnan(n).sum()}, Inf count: {np.isinf(n).sum()}")
        debug_logger.info(f"k - min: {k.min():.6e}, max: {k.max():.6e}, mean: {k.mean():.6e}")
        debug_logger.info(f"k - NaN count: {np.isnan(k).sum()}, Inf count: {np.isinf(k).sum()}")
        debug_logger.info("="*60)
        debug_logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 输出结果
    out = pd.DataFrame()
    out["wavelength"] = wavelengths_m * 1e9 #默认输出nm单位
    out["n"] = n
    out["k"] = k
    out.to_csv(args.out, index=False)

    if args.fig_path != None:
        plot_fig(n,k,wavelengths_m* 1e9, args.fig_path)



if __name__ == "__main__":
    cli()





