import numpy as np  # module for math operations
import matplotlib.pyplot as plt  # module for plotting
import scipy.signal as signal  # module for finding peaks
import csv  # module for reading csv files
from scipy.ndimage import minimum_filter, maximum_filter  # for baseline correction
from scipy.integrate import simpson  # for calculating the area under peaks
from scipy.signal import find_peaks, peak_widths, peak_prominences  # for finding peaks and their properties
import pandas as pd  # module for creating and saving tables
import argparse  # module for handling command-line arguments
from datetime import datetime  # module for working with dates and times
from mpl_toolkits.mplot3d import Axes3D # for waterfall plots
import os  # module for handling file paths

# Wanhe: this bit loads data from the .csv file. The ones you sent are a little odd, otherwise this would be simpler...
def load_hplc_data(file_path, delimiter='\t', encoding='ISO-8859-1'):
    """
    load hplc data from a csv file, handle encoding issues, remove null bytes, and return a numpy array.
    - file_path: path to the csv file.
    - delimiter: delimiter used in the csv file (default is tab '\t').
    - encoding: encoding of the file (default is 'iso-8859-1').
    """
    
    data = []  # list to store the data
    try:
        with open(file_path, 'rb') as file:  # open file in binary mode
            raw_data = file.read()  # read the file content
            
            # remove null values
            cleaned_data = raw_data.replace(b'\x00', b'').decode(encoding)

            decoded_lines = cleaned_data.splitlines()

            csv_reader = csv.reader(decoded_lines, delimiter=delimiter)
            for row in csv_reader:
                # we want each row to have 2 columns, and we convert them to floats
                if len(row) == 2:
                    try:
                        data.append([float(row[0]), float(row[1])])
                    except ValueError:
                        # skip rows that can't be converted to floats
                        continue
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    
    data_array = np.array(data)  # convert list to numpy array

    return data_array


# function for "rolling ball" baseline correction; makes a nice baseline.
def rolling_ball_baseline(y, radius=50):
    """
    baseline correction using a "rolling ball" method.
    - y: hplc detector response data as a numpy array.
    - radius: radius of the ball, which determines the baseline "smoothness".
    """
    # apply minimum and maximum filters to smooth the baseline
    min_filtered = minimum_filter(y, size=radius)
    max_filtered = maximum_filter(min_filtered, size=radius)
    return max_filtered

# this is the part that actually analyses the data
def analyze_hplc_data(filename, 
                      x_start=None, x_end=None, nt = 0.0,
                      baseline_radius=50, 
                      peak_distance=20, 
                      rel_height=0.5, 
                      threshold=0.1,  
                      baseline_correction=True, 
                      plot=True, 
                      raw_offset = 0.0):
    """
    analyze hplc data by detecting peaks and calculating peak areas.
    - filename: path to the csv file containing the hplc data.
    - x_start, x_end: specify the range of x-values to analyze. default is the entire dataset.
    - baseline_radius: radius for rolling ball baseline correction. default is 50.
    - peak_distance: minimum distance between peaks for detection. default is 20.
    - rel_height: relative height for measuring peak width. default is 0.5 (full width at half maximum).
    - threshold: minimum height for a peak to be considered. default is 0.1.
    - baseline_correction: apply baseline correction? default is true.
    - plot: whether to plot the results. default is true.
    """

    # load the data
    data = load_hplc_data(filename)
    if data is None:
        print(f"Failed to load data from {filename}. Please check the file path and format.")
        return
        
    original_x_data = data[:, 0]  # time data (original, uncropped)
    original_y_data = data[:, 1]  # detector response data (original, uncropped)

    # clip data to the specified x-window if provided
    if x_start is not None or x_end is not None:
        mask = (original_x_data >= (x_start if x_start is not None else original_x_data.min())) & \
               (original_x_data <= (x_end if x_end is not None else original_x_data.max()))
        x_data = original_x_data[mask]
        y_data = original_y_data[mask]
    else:
        x_data = original_x_data  # No cropping, use original data
        y_data = original_y_data

    # apply Rolling Ball Baseline Correction if enabled
    if baseline_correction:
        baseline = rolling_ball_baseline(y_data, radius=baseline_radius)
        corrected_data = y_data - baseline  # subtract baseline from data
    else:
        corrected_data = y_data
        
    if nt != 0.0 and nt != 3.14:
        print(f'Normalizing data to time {nt} min')
        corrected_data = normalize_spectrum(x_data, corrected_data, nt)
        y_data = normalize_spectrum(x_data, y_data, nt)
        
    if nt == 3.14:
        print(f'Normalizing data to y-max')
        corrected_data = corrected_data / np.max(corrected_data)
        y_data = y_data / np.max(y_data)
        
    # create a plot for the data
    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(x_data, y_data + raw_offset, label='Raw data')  # Plot the cropped data
        if baseline_correction:
            plt.plot(x_data, corrected_data, label='Corrected data')
        plt.legend()
        plt.xlabel('Time / min')
        plt.ylabel('Normalized Detector Response / a.u.')
        
    # get peaks in the corrected data
    peaks, _ = find_peaks(corrected_data, distance=peak_distance, height=threshold)

    # filter out peaks based on the provided x_start and x_end using the cropped time data
    peak_mask = (x_data[peaks] >= (x_start if x_start is not None else x_data.min())) & \
                (x_data[peaks] <= (x_end if x_end is not None else x_data.max()))
    peaks = peaks[peak_mask]  # Only keep peaks within the specified range

    prominences = peak_prominences(corrected_data, peaks)[0]  # get peak prominences
    widths = peak_widths(corrected_data, peaks, rel_height=rel_height)  # get peak widths

    # record peak data
    peak_areas = []
    peak_data = []
    colors = ['orange', 'green']  # alternating colors for peaks

    for i, peak in enumerate(peaks):
        # find the width of the peak
        width = widths[0][i]

        # define the region for integration based on the peak width
        left = max(0, int(peaks[i] - width / 2))
        right = min(len(corrected_data), int(peaks[i] + width / 2))

        # adjust the integration region based on the baseline
        while left > 0 and corrected_data[left] > 0:
            left -= 1

        while right < len(corrected_data) and corrected_data[right] > 0:
            right += 1

        # get the area under the peak using Simpsons rule
        area = simpson(corrected_data[left:right])
        peak_areas.append(area)

        # store the peak information (use x_data for the cropped time)
        peak_data.append({
            "Peak #": i + 1,
            "Position (x)": x_data[peak],  # Use cropped time data
            "Height": prominences[i],
            "Area": area
        })

        # Plot the peak if requested
        if plot:
            plt.annotate(f'{i+1}', 
                         (x_data[peak], corrected_data[peak]), 
                         textcoords="offset points", 
                         xytext=(0, 10), 
                         ha='center', fontsize=8, color='red')
            plt.fill_between(x_data[left:right], corrected_data[left:right], color=colors[i % len(colors)], alpha=0.3)

    # get the total area to find the % area for each peak
    total_area = sum(peak_areas)
    for peak in peak_data:
        peak["% Area"] = (peak["Area"] / total_area) * 100

    peak_df = pd.DataFrame(peak_data)

    # adjust plot limits and display the plot if requested
    base_filename = os.path.splitext(filename)[0]  # get the base filename without extension
    if plot:
        plt.xlim([x_data.min(), x_data.max()])
        plt.ylim([corrected_data.min() - 5, corrected_data.max() + 5])
        
        # save the plot as a PNG image
        output_png = f"{base_filename}_report.png"  # create the new filename
        plt.title(f'HPLC data for: {base_filename}')
        plt.savefig(output_png)
        print(f"Plot saved as {output_png}")
        
        # show the plot!!!
        plt.show()

    # save peak data to a CSV file
    output_csv = f"{base_filename}_report.csv"
    peak_df.to_csv(output_csv, index=False)
    print(f"Peak data saved to {output_csv}")

    # store the function parameters and save them to a TXT file
    output_txt = f"{base_filename}_report.txt"
    with open(output_txt, "w") as log_file:
        log_file.write("HPLC Data Analysis Report\n")
        log_file.write(f"Filename: {filename}\n")
        log_file.write(f"Date/Time: {datetime.now()}\n\n")
        log_file.write(f"x_start: {x_start}\n")
        log_file.write(f"x_end: {x_end}\n")
        log_file.write(f"baseline_radius: {baseline_radius}\n")
        log_file.write(f"peak_distance: {peak_distance}\n")
        log_file.write(f"rel_height: {rel_height}\n")
        log_file.write(f"threshold: {threshold}\n")  # Log the threshold parameter
        log_file.write(f"baseline_correction: {baseline_correction}\n")
        log_file.write(f"plot: {plot}\n")
    print(f"Parameters saved to {output_txt}")

    return peak_df
    
# Wanhe: this new function to normalises data based peak height at time `nt`
def normalize_spectrum(x_data, y_data, nt):
    # Find the index of the closest time point to the specified normalization time `nt`
    idx = (np.abs(x_data - nt)).argmin()
    normalization_factor = y_data[idx]
    if normalization_factor != 0:
        return y_data / normalization_factor  # Normalize the spectrum
    else:
        return y_data  # In case of zero at normalization point, return original data
        
# Wanhe: this function lets you analyse and normalize HPLC data from multiple files

from mpl_toolkits.mplot3d import Axes3D

def analyze_all_hplc_data(files, nt, x_start=None, x_end=None, baseline_radius=50, peak_distance=20, rel_height=0.5, threshold=0.1, baseline_correction=True, plot=True, raw_offset=0.0, cmap='viridis', height_max=None):
    spectra = []
    labels = []

    for filename in files:
        data = load_hplc_data(filename)
        if data is None:
            print(f"Failed to load data from {filename}. Skipping file.")
            continue
        
        x_data = data[:, 0]
        y_data = data[:, 1]

        baseline = rolling_ball_baseline(y_data, radius=baseline_radius)
        corrected_data = y_data - baseline

        normalized_data = normalize_spectrum(x_data, corrected_data, nt)

        # Clip data to height_max if specified
        if height_max is not None:
            normalized_data = np.clip(normalized_data, None, height_max)

        # Append both x_data and normalized_data to spectra as a tuple
        spectra.append((x_data, normalized_data))
        labels.append(os.path.basename(filename))
        
        # Also run individual analysis on the file
        analyze_hplc_data(filename, x_start=x_start, x_end=x_end, nt=nt, baseline_radius=baseline_radius, 
                          peak_distance=peak_distance, rel_height=rel_height, threshold=threshold, 
                          baseline_correction=baseline_correction, plot=plot, raw_offset=raw_offset)

    # Create 3D waterfall plot
    fig = plt.figure(figsize=(16, 8), constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d')

    # Use a color gradient for better visual differentiation
    color_map = plt.get_cmap(cmap)
    num_files = len(spectra)

    # Plot each spectrum in 3D space
    for i, (x_data, spectrum) in enumerate(spectra):
        # clip data to the specified x-window if provided
        if x_start is not None or x_end is not None:
            mask = (x_data >= (x_start if x_start is not None else x_data.min())) & \
                   (x_data <= (x_end if x_end is not None else x_data.max()))
            x_data = x_data[mask]
            y_data = spectrum[mask]
        else:
            x_data = x_data  # no cropping; use original data
            y_data = spectrum
        
        z_offset = i * 0.5  # Offset each plot on the z-axis for the waterfall effect
        ax.plot(x_data, [z_offset] * len(x_data), y_data, label=labels[i], color=color_map(i / num_files))

        # Add labels directly to the plot
        ax.text(x_data[-1] + 1, z_offset, spectrum[-1], labels[i], fontsize=10, ha='left', va='bottom')

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('')  # Remove the Y-axis label as requested
    ax.set_zlabel('Normalized Intensity')
    ax.set_title('3D Waterfall Plot of HPLC Data')
    ax.grid(True)  # Add grid for better visibility

    plt.show()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Analyze HPLC data by detecting peaks and calculating peak areas.')
    
    parser.add_argument('--all', action='store_true', help='Process all .csv files in the current directory.')
    parser.add_argument('--nt', type=float, default=0.0, help='Normalization time for spectra.')
    parser.add_argument('--x_start', type=float, default=None, help='Start of x-range for analysis (default: entire dataset).')
    parser.add_argument('--x_end', type=float, default=None, help='End of x-range for analysis (default: entire dataset).')
    parser.add_argument('--baseline_radius', type=int, default=50, help='Radius for rolling ball baseline correction (default: 50).')
    parser.add_argument('--peak_distance', type=int, default=20, help='Minimum distance between peaks for detection (default: 20).')
    parser.add_argument('--rel_height', type=float, default=0.5, help='Relative height for measuring peak width (default: 0.5, i.e., FWHM).')
    parser.add_argument('--threshold', type=float, default=0.1, help='Minimum height required for a peak to be considered (default: 0.1).')
    parser.add_argument('--no_baseline_correction', action='store_true', help='Disable baseline correction (default: enabled).')
    parser.add_argument('--no_plot', action='store_true', help='Disable plotting of results (default: enabled).')
    parser.add_argument('--raw_offset', type=float, default=0.0, help='Offset raw data by this many units in y- to aid visualization (defaults to 0.0)')
    parser.add_argument('--cmap', type=str, default='viridis', help='Colormap for the waterfall plot (default: viridis).')
    parser.add_argument('--height_max', type=float, default=None, help='Maximum height for clipping the spectra (default: None).')
    parser.add_argument('--filename', type=str, default=None, help='Filename of raw data to load, in .csv format')
        
    args = parser.parse_args()    
    
    # Wanhe - this makes a list of all .csv files, provided the --all flag is set
    if args.all:
        files = [f for f in os.listdir('.') if f.lower().endswith('.csv') and not f.lower().endswith('_report.csv')]
        if args.nt == 0.0:
            print("Please provide a normalization time (--nt) when using the --all flag; defaulting to t = 0.0.")
        analyze_all_hplc_data(files, nt=args.nt, x_start=args.x_start, x_end=args.x_end, 
                              baseline_radius=args.baseline_radius, peak_distance=args.peak_distance, 
                              rel_height=args.rel_height, threshold=args.threshold, 
                              baseline_correction=not args.no_baseline_correction,  # Fixing logic here
                              plot=not args.no_plot,  # Fixing logic here
                              raw_offset=args.raw_offset, height_max = args.height_max)
    else:
        if args.filename:
            analyze_hplc_data(
                filename=args.filename,
                x_start=args.x_start,
                x_end=args.x_end,
                nt=args.nt,
                baseline_radius=args.baseline_radius,
                peak_distance=args.peak_distance,
                rel_height=args.rel_height,
                threshold=args.threshold,
                baseline_correction=not args.no_baseline_correction,  # Fixing logic here
                plot=not args.no_plot,  # Fixing logic here
                raw_offset=args.raw_offset
            )

if __name__ == '__main__':
    main()