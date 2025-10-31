# Material Complex Refractive Index Inversion Tool

This tool is used to invert and calculate the complex refractive index n(λ) + ik(λ) of semi-transparent or transparent materials from measured reflectance (R) and transmittance (T) spectra. Note that it does not support negative refractive index inversion for opaque (T=0 or very low T) cases. It supports two inversion methods: **TM method** (Transmission Method) and **RT method** (Ray Tracing Method), suitable for optical property analysis of single-layer thin film materials under normal incidence conditions.

## Features

- ✅ Supports both TM and RT inversion methods
- ✅ Automatic data validation (energy conservation check)
- ✅ Roughness correction function (optional)
- ✅ Debug mode: Detailed recording of intermediate variables for each data point
- ✅ Automatic visualization generation
- ✅ Supports multiple units (wavelength: nm, μm, m, Å; thickness: nm, μm, mm, m, Å)

## Installation

```bash
pip install -r requirements.txt
```

Dependencies:
- numpy >= 1.24
- scipy >= 1.11
- pandas >= 2.0
- matplotlib >= 3.7

## Code Structure

- nk_inversion.py (main program for nk value conversion)
- fit_rt_sampling.py (code for data fitting and sampling to corresponding wavelength data points)
- requirements.txt
- README.md

Other files are reference papers and R, T curve data cases extracted from the papers.
R, T curve data are obtained by taking screenshots from papers, extracting data using chart extraction software, then fitting curves and sampling to unify to the same wavelength, so there may be significant errors (personal testing found that compared to the values in the paper, the differences are quite large, but the shapes are correct).
When using FDTD method to model and simulate R, T curves to inversely solve material nk values, both methods show acceptable accuracy.

## Input Data Format

The CSV file must contain the following three columns:

| Column Name | Description | Range |
|-------------|-------------|-------|
| `wavelength_nm` | Wavelength value (unit specified by `--wl_unit`) | > 0 |
| `R` | Reflectance | 0 ≤ R ≤ 1 |
| `T` | Transmittance | 0 ≤ T ≤ 1 |

**Notes:**
- Data must satisfy energy conservation: R + T ≤ 1
- If there are data points with R+T > 1, the program will report an error and list the erroneous data points
- R or T values of 0 will be automatically replaced with 1e-3 to avoid calculation errors

## Usage

### Basic Command Format

```bash
python nk_inversion.py \
  --csv <input CSV file path> \
  --out <output CSV file path> \
  --thickness <thickness value> \
  --th_unit <thickness unit> \
  --wl_unit <wavelength unit> \
  --method <inversion method> \
  [optional parameters]
```

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--csv` | Input CSV file path | `input/data.csv` |
| `--out` | Output CSV file path | `output/nk_results.csv` |
| `--thickness` | Sample thickness value | `5.19` |
| `--th_unit` | Thickness unit (nm/um/mm/m/ang) | `mm` |
| `--wl_unit` | Wavelength unit in input CSV (nm/um/m/ang) | `nm` |
| `--method` | Inversion method: `tm` or `rt` | `rt` |

### Optional Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--sigma` | Roughness correction parameter (RMS roughness) | `0.5` |
| `--fig_path` | Output image file path | `output/result.png` |
| `--debug_log` | Debug log file path | `output/debug.log` |

### Usage Examples

#### Example 1: RT Method Inversion (with debugging)

```bash
python nk_inversion.py \
  --csv input/fit_test_borate_RT_sampled.csv \
  --out output/borate_nk_rt.csv \
  --thickness 5.19 \
  --th_unit mm \
  --wl_unit nm \
  --method rt \
  --fig_path output/borate_rt.png \
  --debug_log output/debug_borate.log
```

#### Example 2: TM Method Inversion (with roughness correction)

```bash
python nk_inversion.py \
  --csv input/fit_test_corning_RT_sampled.csv \
  --out output/corning_nk_tm.csv \
  --thickness 0.92 \
  --th_unit mm \
  --wl_unit nm \
  --method tm \
  --sigma 0.5 \
  --fig_path output/corning_tm.png \
  --debug_log output/debug_corning_tm.log
```

## Output Description

### Main Output File

After running, a CSV file (specified by `--out` parameter) will be generated, containing the following columns:

| Column Name | Description | Unit |
|-------------|-------------|------|
| `wavelength` | Wavelength | nm |
| `n` | Real part of refractive index | - |
| `k` | Imaginary part of refractive index (extinction coefficient) | - |

### Debug Output (when `--debug_log` is enabled)

When debug mode is enabled, two files will be generated:

1. **Log file** (`*.log`): Contains summary information of the calculation process
   - Input parameter statistics
   - Intermediate variable statistics (minimum, maximum, average values)
   - Anomaly statistics (NaN, Inf counts)

2. **Detailed data file** (`*_debug.csv`): Contains detailed intermediate variables for each data point
   - TM method includes: wavelength, R, T, a, b, eta, R_as, sqrt_arg, n, k, and anomaly markers
   - RT method includes: wavelength, R, T, a, b, c, R_as, denominator_k, log_arg, sqrt_arg, n, k, and detailed problem markers

### Image Output (when `--fig_path` is specified)

Generates an image containing two subplots:
- Top plot: Real part of refractive index n vs wavelength
- Bottom plot: Extinction coefficient k vs wavelength

## Theory

### TM Method (Transmission Method)

The TM method is based on the Transfer Matrix Method (TMM) and is suitable for samples with high transmittance. The calculation formulas are as follows:

**Intermediate variables:**
- a = T² - (1-R)²
- b = √(a² + 4T²)
- η = (a+b)/(2T)  [Attenuation factor, Equation (25)]

**Calculation results:**
- k = -[λ/(4πt)] · ln(η)  [Equation (27)]
- R_as = 2R/(2+a+b)  [Equation (28)]
- n = (1+R_as)/(1-R_as) + √[(4R_as/(1-R_as)²) - k²]  [Equation (29)]

### RT Method (Ray Tracing Method)

The RT method is based on ray tracing and is suitable for samples with strong reflection and transmission signals. The calculation formulas are as follows:

**Intermediate variables:**
- a = 2 + T² - (1-R)²
- b = √(a² - 4(2-R)R)
- c = a - b

**Calculation results:**
- R_as = c/[2(2-R)]  [Equation (20)]
- k = [λ/(4πt)] · ln[(Tc)/denominator_k]  [Equation (21)]
  where denominator_k = 2(2-R)R + b - a
- n = (1+R_as)/(1-R_as) + √[(4R_as/(1-R_as)²) - k²]  [Equation (23)]

### Roughness Correction (Optional)

When the sample surface has roughness, the `--sigma` parameter can be used for correction:

- R_corrected = R_measured · exp(-16π²σ²/λ²)  [Equation (37)]
- T_corrected = T_measured · exp(-16π²σ²/λ²)  [Equation (38)]

where σ is the root mean square roughness (RMS roughness).

## Citation

This tool implements methods based on the following paper:

El-Zaiat, S.Y., Determination of the complex refractive index of a thick slab material from its spectral reflectance and transmittance at normal incidence. Optik, 2013. 124(2): p. 157-161. "https://doi.org/https://doi.org/10.1016/j.ijleo.2011.11.039"

The paper provides analytical methods for inverting complex refractive index from reflectance and transmittance spectra, including:
- Theoretical derivation and formulas of the Transfer Matrix Method (TM method)
- Theoretical derivation and formulas of the Ray Tracing Method (RT method)
- Effects of surface roughness on measurement results and correction methods

If you use this tool, please cite the relevant original paper.

## Notes

1. **Energy Conservation Validation**: The program automatically checks R+T ≤ 1. If there are data points violating energy conservation, the program will report an error and stop.

2. **Zero Value Handling**: When R or T is 0, it will be automatically replaced with 1e-3 to avoid division-by-zero errors in calculations.

3. **Anomaly Handling**:
   - In the RT method, if denominator_k approaches zero or is negative, it may cause abnormal k values
   - When sqrt_arg is negative, it will result in NaN values for n
   - The CSV file in debug mode will mark these anomaly conditions

4. **Method Selection Recommendations**:
   - **TM method**: Suitable for high transmittance samples (large T)
   - **RT method**: Suitable for samples with strong reflection and transmission signals

5. **Unit Consistency**: Ensure that the units of input data match the `--wl_unit` and `--th_unit` parameters.

## Troubleshooting

If you encounter calculation anomalies (NaN or Inf values), it is recommended to:

1. Enable debug mode: Use the `--debug_log` parameter
2. Check the detailed debug CSV file: Examine which data points have abnormal intermediate variables
3. Check input data: Ensure R+T ≤ 1 and data quality is good
4. Try roughness correction: If the sample surface is rough, use the `--sigma` parameter
5. Try the other method: If one method fails, try the other method

## License

Please refer to the project license file.

---

# 材料复折射率反演工具

本工具用于从测量的反射率（R）和透射率（T）光谱，反演计算半透明或透明材料的复折射率 n(λ) + ik(λ)。注意不支持不透明（T=0或T非常低）情况下的负折射率反演。支持两种反演方法：**TM方法**（Transmission Method）和**RT方法**（Ray Tracing Method），适用于单层薄膜材料在法向入射条件下的光学特性分析。

## 功能特点

- ✅ 支持 TM 和 RT 两种反演方法
- ✅ 自动数据验证（能量守恒检查）
- ✅ 粗糙度修正功能（可选）
- ✅ 调试模式：详细记录每个数据点的中间变量
- ✅ 自动生成可视化图表
- ✅ 支持多种单位（波长：nm, μm, m, Å；厚度：nm, μm, mm, m, Å）

## 安装依赖

```bash
pip install -r requirements.txt
```

依赖包：
- numpy >= 1.24
- scipy >= 1.11
- pandas >= 2.0
- matplotlib >= 3.7

## 代码库结构
- nk_inversion.py （主程序，用于nk值转换）
- fit_rt_sampling.py (用于数据拟合采样成对应波长数据点的代码)
- requirements.txt
- README.md

除此以外的文件为对应的参考论文，和我从论文中获取的R，T曲线数据案例。
R,T曲线数据是通过在论文中截图，使用图表提取软件获取数据后，经过曲线拟合后采样，统一成同一个波长的，所以可能存在较大的误差（我个人测试发现对比论文中的数值发现差的蛮大的，不过形状是对的）。
用FDTD方法建模仿真的R，T曲线来反向求解材料nk值，发现两种方法的精度都还可以。

## 输入数据格式

CSV 文件必须包含以下三列：

| 列名 | 说明 | 范围 |
|------|------|------|
| `wavelength_nm` | 波长值（单位由 `--wl_unit` 指定） | > 0 |
| `R` | 反射率 | 0 ≤ R ≤ 1 |
| `T` | 透射率 | 0 ≤ T ≤ 1 |

**注意事项：**
- 数据必须满足能量守恒：R + T ≤ 1
- 如果数据中存在 R+T > 1 的点，程序会报错并列出错误数据点
- R 或 T 为 0 的值会被自动替换为 1e-3 以避免计算异常

## 使用方法

### 基本命令格式

```bash
python nk_inversion.py \
  --csv <输入CSV文件路径> \
  --out <输出CSV文件路径> \
  --thickness <厚度值> \
  --th_unit <厚度单位> \
  --wl_unit <波长单位> \
  --method <反演方法> \
  [可选参数]
```

### 必需参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--csv` | 输入CSV文件路径 | `input/data.csv` |
| `--out` | 输出CSV文件路径 | `output/nk_results.csv` |
| `--thickness` | 样品厚度数值 | `5.19` |
| `--th_unit` | 厚度单位（nm/um/mm/m/ang） | `mm` |
| `--wl_unit` | 输入CSV中波长单位（nm/um/m/ang） | `nm` |
| `--method` | 反演方法：`tm` 或 `rt` | `rt` |

### 可选参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `--sigma` | 粗糙度修正参数（RMS roughness） | `0.5` |
| `--fig_path` | 输出图像文件路径 | `output/result.png` |
| `--debug_log` | 调试日志文件路径 | `output/debug.log` |

### 使用示例

#### 示例1：RT方法反演（带调试）

```bash
python nk_inversion.py \
  --csv input/fit_test_borate_RT_sampled.csv \
  --out output/borate_nk_rt.csv \
  --thickness 5.19 \
  --th_unit mm \
  --wl_unit nm \
  --method rt \
  --fig_path output/borate_rt.png \
  --debug_log output/debug_borate.log
```

#### 示例2：TM方法反演（带粗糙度修正）

```bash
python nk_inversion.py \
  --csv input/fit_test_corning_RT_sampled.csv \
  --out output/corning_nk_tm.csv \
  --thickness 0.92 \
  --th_unit mm \
  --wl_unit nm \
  --method tm \
  --sigma 0.5 \
  --fig_path output/corning_tm.png \
  --debug_log output/debug_corning_tm.log
```

## 输出说明

### 主输出文件

运行后会生成一个CSV文件（`--out` 参数指定），包含以下列：

| 列名 | 说明 | 单位 |
|------|------|------|
| `wavelength` | 波长 | nm |
| `n` | 折射率实部 | - |
| `k` | 折射率虚部（消光系数） | - |

### 调试输出（启用 `--debug_log` 时）

当启用调试模式时，会生成两个文件：

1. **日志文件**（`*.log`）：包含计算过程的摘要信息
   - 输入参数统计
   - 中间变量的统计信息（最小值、最大值、平均值）
   - 异常值统计（NaN、Inf数量）

2. **详细数据文件**（`*_debug.csv`）：包含每个数据点的详细中间变量
   - TM方法包含：波长、R、T、a、b、eta、R_as、sqrt_arg、n、k 及异常标记
   - RT方法包含：波长、R、T、a、b、c、R_as、denominator_k、log_arg、sqrt_arg、n、k 及详细问题标记

### 图像输出（指定 `--fig_path` 时）

生成包含两个子图的图像：
- 上图：折射率实部 n 随波长的变化
- 下图：消光系数 k 随波长的变化

## 原理说明

### TM方法（Transmission Method）

TM方法基于传输矩阵法（TMM），适用于透射率较大的样品。计算公式如下：

**中间变量：**
- a = T² - (1-R)²
- b = √(a² + 4T²)
- η = (a+b)/(2T)  [衰减因子，Equation (25)]

**计算结果：**
- k = -[λ/(4πt)] · ln(η)  [Equation (27)]
- R_as = 2R/(2+a+b)  [Equation (28)]
- n = (1+R_as)/(1-R_as) + √[(4R_as/(1-R_as)²) - k²]  [Equation (29)]

### RT方法（Ray Tracing Method）

RT方法基于光线追迹法，适用于反射和透射信号较强的样品。计算公式如下：

**中间变量：**
- a = 2 + T² - (1-R)²
- b = √(a² - 4(2-R)R)
- c = a - b

**计算结果：**
- R_as = c/[2(2-R)]  [Equation (20)]
- k = [λ/(4πt)] · ln[(Tc)/denominator_k]  [Equation (21)]
  其中 denominator_k = 2(2-R)R + b - a
- n = (1+R_as)/(1-R_as) + √[(4R_as/(1-R_as)²) - k²]  [Equation (23)]

### 粗糙度修正（可选）

当样品表面存在粗糙度时，可以使用 `--sigma` 参数进行修正：

- R_corrected = R_measured · exp(-16π²σ²/λ²)  [Equation (37)]
- T_corrected = T_measured · exp(-16π²σ²/λ²)  [Equation (38)]

其中 σ 为均方根粗糙度（RMS roughness）。

## 引用信息

本工具实现的方法基于如下论文：

El-Zaiat, S.Y., Determination of the complex refractive index of a thick slab material from its spectral reflectance and transmittance at normal incidence. Optik, 2013. 124(2): p. 157-161."https://doi.org/https://doi.org/10.1016/j.ijleo.2011.11.039"

论文提供了从反射率和透射率光谱反演复折射率的解析方法，包括：
- 传输矩阵法（TM方法）的理论推导与公式
- 光线追迹法（RT方法）的理论推导与公式
- 表面粗糙度对测量结果的影响及修正方法

如使用本工具，请引用相关原始论文。

## 注意事项

1. **能量守恒验证**：程序会自动检查 R+T ≤ 1。如果存在违反能量守恒的数据点，程序会报错并停止。

2. **零值处理**：当 R 或 T 为 0 时，会被自动替换为 1e-3 以避免计算中的除零错误。

3. **异常值处理**：
   - RT方法中，如果 denominator_k 接近零或为负，可能导致 k 值异常
   - sqrt_arg 为负时，会导致 n 值为 NaN
   - 调试模式下的 CSV 文件会标记这些异常情况

4. **方法选择建议**：
   - **TM方法**：适用于高透射率样品（T较大）
   - **RT方法**：适用于反射和透射信号都较强的样品

5. **单位一致性**：确保输入数据的单位与 `--wl_unit` 和 `--th_unit` 参数匹配。

## 错误排查

如果遇到计算异常（NaN或Inf值），建议：

1. 启用调试模式：使用 `--debug_log` 参数
2. 查看详细调试CSV文件：检查哪些数据点的中间变量异常
3. 检查输入数据：确保 R+T ≤ 1，且数据质量良好
4. 尝试粗糙度修正：如果样品表面粗糙，使用 `--sigma` 参数
5. 尝试另一种方法：如果一种方法失败，可以尝试另一种方法

## 许可证

请参考项目许可证文件。
