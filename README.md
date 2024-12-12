# HPLC Data Analysis

This project contains a Python script for analyzing High-Performance Liquid Chromatography (HPLC) data from CSV files. The script provides baseline correction, peak detection, area integration, and various visualization options for HPLC chromatograms. It also supports batch analysis and normalization of multiple files.

## Features

- **Baseline Correction**: Implements a "rolling ball" method to correct the baseline of the data.
- **Peak Detection**: Identifies peaks, calculates their prominence, width, and area.
- **Normalization**: Allows normalization of spectra based on a specified time or maximum intensity.
- **Batch Processing**: Analyzes multiple CSV files and generates results for each.
- **3D Waterfall Plot**: Visualizes normalized spectra from multiple files in a 3D plot.
- **Output Reports**: Generates PNG plots, CSV files, and TXT reports for peak data and parameters.

## Requirements

This script requires the following Python libraries:

- `numpy`
- `matplotlib`
- `scipy`
- `pandas`
- `argparse`

Install the required libraries using pip:

```bash
pip install numpy matplotlib scipy pandas
```

## Usage

The script supports both single-file analysis and batch processing of all CSV files in a directory.

### Command-Line Arguments

| Argument                     | Description                                                                                       | Default        |
|------------------------------|---------------------------------------------------------------------------------------------------|----------------|
| `--all`                      | Analyze all `.csv` files in the current directory.                                               | Off            |
| `--nt <float>`               | Normalization time for spectra.                                                                 | 0.0            |
| `--x_start <float>`          | Start of x-range for analysis.                                                                  | Entire dataset |
| `--x_end <float>`            | End of x-range for analysis.                                                                    | Entire dataset |
| `--baseline_radius <int>`    | Radius for rolling ball baseline correction.                                                    | 50             |
| `--peak_distance <int>`      | Minimum distance between peaks.                                                                 | 20             |
| `--rel_height <float>`       | Relative height for measuring peak width (e.g., 0.5 for full-width at half maximum).            | 0.5            |
| `--threshold <float>`        | Minimum height required for a peak to be considered.                                            | 0.1            |
| `--no_baseline_correction`   | Disable baseline correction.                                                                    | Enabled        |
| `--no_plot`                  | Disable plotting of results.                                                                    | Enabled        |
| `--raw_offset <float>`       | Offset raw data by a specified amount to aid visualization.                                     | 0.0            |
| `--cmap <str>`               | Colormap for the 3D waterfall plot (e.g., `viridis`, `plasma`).                                 | `viridis`      |
| `--height_max <float>`       | Maximum height for clipping spectra.                                                            | None           |
| `--filename <str>`           | Path to a single HPLC data CSV file for analysis.                                               | None           |

### Examples

#### Analysis on a ***Single*** File

This does the analysis on a ***single*** .csv file; the output is only saved after closing the matplotlib window.
```bash
python analyze_hplc_data.py --filename example.csv --x_start 0 --x_end 10 --baseline_radius 30 --peak_distance 15
```

#### Batch Analysis

Analyses ***all*** `.csv` files in the current directory:

```bash
python analyze_hplc_data.py --all --nt 2.0 --x_start 0 --x_end 10 --baseline_radius 30 --peak_distance 15
```

#### Generate a 3D Waterfall Plot

```bash
python analyze_hplc_data.py --all --nt 2.0 --cmap plasma
```

## Output

For each analyzed file, the script generates the following:

1. **PNG Plot**: A plot of the chromatogram with identified peaks and baseline correction.
   - Saved as `<filename>_report.png`.
2. **CSV Report**: A file containing peak positions, heights, areas, and percent areas.
   - Saved as `<filename>_report.csv`.
3. **TXT Report**: A summary of analysis parameters and settings.
   - Saved as `<filename>_report.txt`.

## Script Functions

### `load_hplc_data`
Loads HPLC data from a CSV file, handling null bytes and encoding issues.

### `rolling_ball_baseline`
Applies a "rolling ball" method for baseline correction.

### `analyze_hplc_data`
Main function for analyzing a single HPLC data file. Includes baseline correction, peak detection, and plotting.

### `analyze_all_hplc_data`
Analyzes multiple HPLC data files and generates a 3D waterfall plot.

### `normalize_spectrum`
Normalizes HPLC data to a specific time or maximum intensity.

## Example Output:
A processed HPLC chromatogram:
![LUTIDINE-COMU_report](https://github.com/user-attachments/assets/e9b2a086-ef93-43cc-980e-baf2cc6a058c)
A .csv file contains the peak positions, areas etc.
A .txt file contains the parameters used to generate the report
