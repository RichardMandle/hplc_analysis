# hplc_analysis.py
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy.ndimage import minimum_filter, maximum_filter
from scipy.integrate import simpson
from scipy.signal import find_peaks, peak_widths, peak_prominences
import pandas as pd

def load_hplc_data(file_path, delimiter='\t', encoding='ISO-8859-1'):
    """
    load HPLC data from a CSV file. WE need to removes null bytes, handles encoding,
    and returns a numpy array with two columns (time and detector response). The UoL
    HPLC .csv files are a little odd, Pandas seems to throw errors when trying to read them
    """
    data = []
    try:
        with open(file_path, 'rb') as file:
            raw_data = file.read()
            cleaned_data = raw_data.replace(b'\x00', b'').decode(encoding)
            decoded_lines = cleaned_data.splitlines()
            csv_reader = csv.reader(decoded_lines, delimiter=delimiter)
            for row in csv_reader:
                if len(row) == 2:
                    try:
                        data.append([float(row[0]), float(row[1])])
                    except ValueError:
                        continue
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    data_array = np.array(data)
    return data_array

def rolling_ball_baseline(y, radius=50):
    """
    Use rolling ball (minimum then maximum filter) baseline correction from scipy.ndimage.
    """
    min_filtered = minimum_filter(y, size=radius)
    max_filtered = maximum_filter(min_filtered, size=radius)
    return max_filtered

def normalize_spectrum(x_data, y_data, nt):
    """
    Normalize y_data so that its value at the time closest to nt is 1.
    """
    idx = (np.abs(x_data - nt)).argmin()
    normalization_factor = y_data[idx]
    if normalization_factor != 0:
        return y_data / normalization_factor
    else:
        return y_data

def analyze_hplc_file(filename, 
                      x_start=None, x_end=None, nt=0.0,
                      baseline_radius=50, 
                      peak_distance=20, 
                      rel_height=0.5, 
                      threshold=0.1,  
                      baseline_correction=True, 
                      raw_offset=0.0):
    """
    this is our major analysis function. Here, we are going to:
      - load and (optionally) cropping the data,
      - apply baseline correction (and normalization if requested),
      - find peaks,
      - make a plot of both the raw data (with an optional offset) and the corrected data,
      - annotate peaks with numbers and fill them with.
    
    Returns: a tuple containing loads of things:
      fig           : matplotlib Figure with the plot.
      peak_df       : pandas DataFrame containing peak information.
      x_data        : time (x) data (cropped if x_start/x_end provided).
      y_raw         : raw detector response data (cropped and, if normalized, normalized).
      corrected_data: baseline-corrected (and normalized) data.
      peaks         : array of indices (into x_data) where peaks were detected.
      prominences   : detected peak prominences.
      widths        : detected peak widths (as returned by peak_widths).
    """
    data = load_hplc_data(filename)
    if data is None:
        raise Exception(f"Failed to load data from {filename}")
    
    original_x_data = data[:, 0]
    original_y_data = data[:, 1]
    
    # here we'll "crop" the data if needed.
    if x_start is not None or x_end is not None:
        mask = (original_x_data >= (x_start if x_start is not None else original_x_data.min())) & \
               (original_x_data <= (x_end if x_end is not None else original_x_data.max()))
        x_data = original_x_data[mask]
        y_raw = original_y_data[mask]
    else:
        x_data = original_x_data
        y_raw = original_y_data

    #  baseline correction.
    if baseline_correction:
        baseline = rolling_ball_baseline(y_raw, radius=baseline_radius)
        corrected_data = y_raw - baseline
    else:
        corrected_data = y_raw

    # normalize - either by a specified time or by maximum.
    if nt != 0.0 and nt != 3.14:
        corrected_data = normalize_spectrum(x_data, corrected_data, nt)
        y_raw = normalize_spectrum(x_data, y_raw, nt)
    elif nt == 3.14:
        corrected_data = corrected_data / np.max(corrected_data)
        y_raw = y_raw / np.max(y_raw)
    
    # make the plot.
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_data, y_raw + raw_offset, label='Raw data')
    if baseline_correction:
        ax.plot(x_data, corrected_data, label='Corrected data')
    ax.set_xlabel('Time / min')
    ax.set_ylabel('Normalized Detector Response / a.u.')
    ax.legend()
    
    # look for peaks in the corrected data
    peaks, _ = find_peaks(corrected_data, distance=peak_distance, height=threshold)
    peak_mask = (x_data[peaks] >= (x_start if x_start is not None else x_data.min())) & \
                (x_data[peaks] <= (x_end if x_end is not None else x_data.max()))
    peaks = peaks[peak_mask]
    prominences = peak_prominences(corrected_data, peaks)[0]
    widths = peak_widths(corrected_data, peaks, rel_height=rel_height)
    
    peak_areas = []
    peak_data = []
    colors = ['orange', 'green'] # these are just hard-coded defaults; the GUI might override these {if I can figure out how}
    for i, peak in enumerate(peaks):
        width = widths[0][i]
        left = max(0, int(peak - width / 2))
        right = min(len(corrected_data), int(peak + width / 2))
        while left > 0 and corrected_data[left] > 0:
            left -= 1
        while right < len(corrected_data) and corrected_data[right] > 0:
            right += 1
        area = simpson(corrected_data[left:right], x=x_data[left:right])
        peak_areas.append(area)
        peak_data.append({
            "Peak #": i + 1,
            "Position (x)": x_data[peak],
            "Height": prominences[i],
            "Area": area
        })
        ax.annotate(f'{i+1}', 
                    (x_data[peak], corrected_data[peak]), 
                    textcoords="offset points", 
                    xytext=(0, 10), 
                    ha='center', fontsize=8, color='red')
        ax.fill_between(x_data[left:right], corrected_data[left:right], 
                        color=colors[i % len(colors)], alpha=0.3)
    
    total_area = sum(peak_areas)
    for d in peak_data:
        d["% Area"] = (d["Area"] / total_area) * 100 if total_area != 0 else 0

    peak_df = pd.DataFrame(peak_data)
    
    ax.set_xlim([x_data.min(), x_data.max()])
    ax.set_ylim([corrected_data.min() - 5, corrected_data.max() + 5])
    ax.set_title(f'HPLC Data for: {os.path.basename(filename)}')
    
    return fig, peak_df, x_data, y_raw, corrected_data, peaks, prominences, widths
