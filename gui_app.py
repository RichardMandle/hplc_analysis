# gui_app.py
import sys
import os
import numpy as np
from datetime import datetime

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QPushButton, QListWidget,
    QVBoxLayout, QHBoxLayout, QFileDialog, QGroupBox, QFormLayout,
    QLineEdit, QSpinBox, QDoubleSpinBox, QCheckBox, QTableWidget,
    QTableWidgetItem, QSplitter, QMessageBox, QColorDialog
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)

# get key functions from the analysis module.
from hplc_analysis import analyze_hplc_file, peak_widths, peak_prominences

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("HPLC Analysis Tool")
        self.resize(1400, 900)
        self.files = []          # full paths for loaded CSV files
        self.current_file = None
        self.last_peak_df = None  # pd.df for peak table
        self.last_parameters = {} # the fitting options for saving
        
        # store analysis data (x_data, y_raw, corrected_data, peaks, etc.)
        self.last_analysis_data = {}
        
        # we need some default plot options:
        self.raw_line_color = "blue"
        self.corrected_line_color = "green"
        self.peak_fill_color1 = "orange"
        self.peak_fill_color2 = "green"
        self.dashed_line_color = "gray"
        self.peak_fill_alpha = 0.3
        self.dashed_lines_enabled = True
        
        # we want a mode for deleting peaks by clicking
        self.interactive_mode = False  # now: "Click to Delete Peak(s)"
        self.cid_click = None  # Connection ID for click events
        
        self._create_widgets()
        self._create_layouts()
        self._create_connections()
    
    def _create_widgets(self):
        # this is the left panel: load button at top, file list in the middle, and options below.
        self.loadButton = QPushButton("Load CSV Files")
        self.fileListWidget = QListWidget()
        
        # set a list of analysis/fitting options
        self.optionsGroup = QGroupBox("Fitting Options")
        self.xStartEdit = QLineEdit()
        self.xEndEdit = QLineEdit()
        self.ntEdit = QLineEdit()  # Normalization time
        self.baselineRadiusEdit = QSpinBox()
        self.baselineRadiusEdit.setRange(1, 1000)
        self.baselineRadiusEdit.setValue(50)
        self.peakDistanceEdit = QSpinBox()
        self.peakDistanceEdit.setRange(1, 1000)
        self.peakDistanceEdit.setValue(20)
        self.relHeightEdit = QDoubleSpinBox()
        self.relHeightEdit.setRange(0.0, 1.0)
        self.relHeightEdit.setSingleStep(0.1)
        self.relHeightEdit.setValue(0.5)
        self.thresholdEdit = QDoubleSpinBox()
        self.thresholdEdit.setRange(0.0, 1000.0)
        self.thresholdEdit.setSingleStep(0.1)
        self.thresholdEdit.setValue(0.1)
        self.rawOffsetEdit = QDoubleSpinBox()
        self.rawOffsetEdit.setRange(-1000.0, 1000.0)
        self.rawOffsetEdit.setValue(0.0)
        self.baselineCorrectionCheck = QCheckBox("Baseline Correction")
        self.baselineCorrectionCheck.setChecked(True)
        self.toggleDashedLinesCheck = QCheckBox("Show Dashed Lines")
        self.toggleDashedLinesCheck.setChecked(True)
        
        # Incase we need it, lets have some options to change the look of the plot.
        # These should work OK; pops out an RGB selection window
        self.rawLineColorButton = QPushButton("Select Raw Line Color")
        self.correctedLineColorButton = QPushButton("Select Corrected Line Color")
        self.peakFillColor1Button = QPushButton("Select Peak Fill Color #1")
        self.peakFillColor2Button = QPushButton("Select Peak Fill Color #2")
        self.dashedLineColorButton = QPushButton("Select Dashed Line Color")
        
        # initial background colours.
        self._update_color_button(self.rawLineColorButton, self.raw_line_color)
        self._update_color_button(self.correctedLineColorButton, self.corrected_line_color)
        self._update_color_button(self.peakFillColor1Button, self.peak_fill_color1)
        self._update_color_button(self.peakFillColor2Button, self.peak_fill_color2)
        self._update_color_button(self.dashedLineColorButton, self.dashed_line_color)
        
        # peak fill alpha
        self.peakFillAlphaEdit = QDoubleSpinBox()
        self.peakFillAlphaEdit.setRange(0.0, 1.0)
        self.peakFillAlphaEdit.setSingleStep(0.05)
        self.peakFillAlphaEdit.setValue(self.peak_fill_alpha)
        
        # this click box makes peaks deletable 
        self.interactiveModeCheck = QCheckBox("Click to Delete Peak")
        self.interactiveModeCheck.setChecked(False)
        
        self.runButton = QPushButton("Run Analysis")
        self.saveButton = QPushButton("Save Analysis")
        
        # main panel: plot the data (as canvas) on top, with peak table below.
        self.figure = plt.figure(figsize=(10, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        self.peakTable = QTableWidget()
        self.peakTable.setColumnCount(5)
        self.peakTable.setHorizontalHeaderLabels(["Peak #", "Position (x)", "Height", "Area", "% Area"])
        self.peakTable.horizontalHeader().setStretchLastSection(True)
    
    def _create_layouts(self):
        # layout for options.
        optionsLayout = QFormLayout()
        optionsLayout.addRow("x_start:", self.xStartEdit)
        optionsLayout.addRow("x_end:", self.xEndEdit)
        optionsLayout.addRow("Normalization Time (nt):", self.ntEdit)
        optionsLayout.addRow("Baseline Radius:", self.baselineRadiusEdit)
        optionsLayout.addRow("Peak Distance:", self.peakDistanceEdit)
        optionsLayout.addRow("Relative Height:", self.relHeightEdit)
        optionsLayout.addRow("Threshold:", self.thresholdEdit)
        optionsLayout.addRow("Raw Offset:", self.rawOffsetEdit)
        optionsLayout.addRow(self.baselineCorrectionCheck)
        optionsLayout.addRow(self.toggleDashedLinesCheck)
        
        # colour choice bits (using colour selection buttons)
        optionsLayout.addRow("Raw Line Color:", self.rawLineColorButton)
        optionsLayout.addRow("Corrected Line Color:", self.correctedLineColorButton)
        optionsLayout.addRow("Peak Fill Color #1:", self.peakFillColor1Button)
        optionsLayout.addRow("Peak Fill Color #2:", self.peakFillColor2Button)
        optionsLayout.addRow("Dashed Line Color:", self.dashedLineColorButton)
        optionsLayout.addRow("Peak Fill Alpha:", self.peakFillAlphaEdit)
        optionsLayout.addRow(self.interactiveModeCheck)
        
        optionsLayout.addRow(self.runButton)
        optionsLayout.addRow(self.saveButton)
        self.optionsGroup.setLayout(optionsLayout)
        
        # set the layout of the left panel 
        leftLayout = QVBoxLayout()
        leftLayout.addWidget(self.loadButton)
        leftLayout.addWidget(self.fileListWidget)
        leftLayout.addWidget(self.optionsGroup)
        leftWidget = QWidget()
        leftWidget.setLayout(leftLayout)
        
        #set the layout of the right panel 
        rightLayout = QVBoxLayout()
        rightLayout.addWidget(self.toolbar)
        rightLayout.addWidget(self.canvas, stretch=3)
        rightLayout.addWidget(self.peakTable, stretch=1)
        rightWidget = QWidget()
        rightWidget.setLayout(rightLayout)
        
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(leftWidget)
        splitter.addWidget(rightWidget)
        splitter.setSizes([400, 1000])
        
        self.setCentralWidget(splitter)
    
    def _create_connections(self):
        # connections for updating on click and so on
        self.loadButton.clicked.connect(self.load_files)
        self.fileListWidget.itemSelectionChanged.connect(self.file_selected)
        self.runButton.clicked.connect(self.run_analysis)
        self.saveButton.clicked.connect(self.save_analysis)
        self.toggleDashedLinesCheck.stateChanged.connect(self.update_plot)
        self.interactiveModeCheck.stateChanged.connect(self.toggle_interactive_mode)
        
        # Connect colour selection buttons.
        self.rawLineColorButton.clicked.connect(lambda: self.choose_color("raw"))
        self.correctedLineColorButton.clicked.connect(lambda: self.choose_color("corrected"))
        self.peakFillColor1Button.clicked.connect(lambda: self.choose_color("peak1"))
        self.peakFillColor2Button.clicked.connect(lambda: self.choose_color("peak2"))
        self.dashedLineColorButton.clicked.connect(lambda: self.choose_color("dashed"))
    
    def _update_color_button(self, button, color):
        # just sets the colour of a colour selection button so it matches the choice (clever)
        button.setStyleSheet("background-color: %s" % color)
    
    def choose_color(self, target):
        """make a QColorDialog to choose a colour for the specified target."""
        current = ""
        if target == "raw":
            current = self.raw_line_color
        elif target == "corrected":
            current = self.corrected_line_color
        elif target == "peak1":
            current = self.peak_fill_color1
        elif target == "peak2":
            current = self.peak_fill_color2
        elif target == "dashed":
            current = self.dashed_line_color
        col = QColorDialog.getColor(QColor(current), self, "Select Color")
        if col.isValid():
            hex_color = col.name()
            if target == "raw":
                self.raw_line_color = hex_color
                self._update_color_button(self.rawLineColorButton, hex_color)
            elif target == "corrected":
                self.corrected_line_color = hex_color
                self._update_color_button(self.correctedLineColorButton, hex_color)
            elif target == "peak1":
                self.peak_fill_color1 = hex_color
                self._update_color_button(self.peakFillColor1Button, hex_color)
            elif target == "peak2":
                self.peak_fill_color2 = hex_color
                self._update_color_button(self.peakFillColor2Button, hex_color)
            elif target == "dashed":
                self.dashed_line_color = hex_color
                self._update_color_button(self.dashedLineColorButton, hex_color)
            self.update_plot()
    
    def load_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select CSV Files", "", "CSV Files (*.csv)")
        if files:
            self.files.extend(files)
            for f in files:
                self.fileListWidget.addItem(os.path.basename(f))
    
    def file_selected(self):
        selectedItems = self.fileListWidget.selectedItems()
        if selectedItems:
            selectedName = selectedItems[0].text()
            for f in self.files:
                if os.path.basename(f) == selectedName:
                    self.current_file = f
                    break
            self.run_analysis()  # Auto-run analysis on file selection.
    
    def run_analysis(self):
        if not self.current_file:
            QMessageBox.warning(self, "No File Selected", "Please select a CSV file.")
            return
        
        try:
            x_start = float(self.xStartEdit.text()) if self.xStartEdit.text() else None
        except ValueError:
            x_start = None
        try:
            x_end = float(self.xEndEdit.text()) if self.xEndEdit.text() else None
        except ValueError:
            x_end = None
        try:
            nt = float(self.ntEdit.text()) if self.ntEdit.text() else 0.0
        except ValueError:
            nt = 0.0
        
        baseline_radius = self.baselineRadiusEdit.value()
        peak_distance = self.peakDistanceEdit.value()
        rel_height = self.relHeightEdit.value()
        threshold = self.thresholdEdit.value()
        raw_offset = self.rawOffsetEdit.value()
        baseline_correction = self.baselineCorrectionCheck.isChecked()
        self.dashed_lines_enabled = self.toggleDashedLinesCheck.isChecked()
        
        # write the current analysis parameters so we can reproduce.
        self.last_parameters = {
            "x_start": x_start,
            "x_end": x_end,
            "nt": nt,
            "baseline_radius": baseline_radius,
            "peak_distance": peak_distance,
            "rel_height": rel_height,
            "threshold": threshold,
            "raw_offset": raw_offset,
            "baseline_correction": baseline_correction
        }
        
        try:
            fig, peak_df, x_data, y_raw, corrected_data, peaks, prominences, widths = analyze_hplc_file(
                self.current_file,
                x_start=x_start,
                x_end=x_end,
                nt=nt,
                baseline_radius=baseline_radius,
                peak_distance=peak_distance,
                rel_height=rel_height,
                threshold=threshold,
                baseline_correction=baseline_correction,
                raw_offset=raw_offset
            )
        except Exception as e:
            QMessageBox.critical(self, "Analysis Error", f"Error: {e}")
            return
        
        # store the analysis data for later replotting.
        self.last_analysis_data = {
            "x_data": x_data,
            "y_raw": y_raw,
            "corrected_data": corrected_data,
            "peaks": list(peaks)  # store as list for editing
        }
        
        self.last_peak_df = peak_df
        self.update_plot()
        self.update_peak_table()
    def update_plot(self):
        """replot the data and peaks using current customization options and non-overlapping integration regions."""
        if not self.last_analysis_data:
            return
        
        self.figure.clf()
        ax = self.figure.add_subplot(111)
        
        x_data = self.last_analysis_data["x_data"]
        y_raw = self.last_analysis_data["y_raw"]
        corrected_data = self.last_analysis_data["corrected_data"]
        peaks = self.last_analysis_data["peaks"]
        raw_offset = self.last_parameters.get("raw_offset", 0.0)
        
        # plot both the raw and corrected data using the selected colours.
        ax.plot(x_data, y_raw + raw_offset, label='Raw data', color=self.raw_line_color)
        if self.last_parameters.get("baseline_correction", True):
            ax.plot(x_data, corrected_data, label='Corrected data', color=self.corrected_line_color)
        ax.set_xlabel('Time / min')
        ax.set_ylabel('Normalized Detector Response / a.u.')
        
        # if requested, draw those dashed lines.
        if self.dashed_lines_enabled:
            thr = self.last_parameters.get("threshold", 0.1)
            ax.axhline(y=thr, linestyle='--', color=self.dashed_line_color, label='Threshold')
            x_start = self.last_parameters.get("x_start", None)
            x_end = self.last_parameters.get("x_end", None)
            if x_start is not None:
                ax.axvline(x=x_start, linestyle='--', color=self.dashed_line_color, label='x_start')
            if x_end is not None:
                ax.axvline(x=x_end, linestyle='--', color=self.dashed_line_color, label='x_end')
        
        # calculate the integration boundaries
        n_peaks = len(peaks)
        boundaries = []
        if n_peaks == 0 or n_peaks == 1:
            boundaries = [x_data[0], x_data[-1]]
        else:
            boundaries.append(x_data[0])
            for i in range(n_peaks - 1):
                # For the boundary between peak i and peak i+1, find the index where corrected_data is minimum.
                start_idx = peaks[i]
                end_idx = peaks[i + 1]
                # Ensure there is a valid region.
                if end_idx > start_idx:
                    region = corrected_data[start_idx:end_idx + 1]
                    valley_index = np.argmin(region) + start_idx
                    boundaries.append(x_data[valley_index])
                else:
                    boundaries.append(x_data[start_idx])
            boundaries.append(x_data[-1])
        
        # Draw peaks and fill their integration regions.
        if n_peaks > 0:
            # For plotting purposes, also store each peak's x-position.
            for i, peak in enumerate(peaks):
                ax.plot(x_data[peak], corrected_data[peak], marker='o', color='red')
                ax.annotate(f'{i+1}', (x_data[peak], corrected_data[peak]),
                            textcoords="offset points", xytext=(0, 10),
                            ha='center', fontsize=8, color='red')
                # Determine the left and right boundaries for the integration region.
                left_bound = boundaries[i]
                right_bound = boundaries[i+1]
                left_index = np.searchsorted(x_data, left_bound)
                right_index = np.searchsorted(x_data, right_bound)
                # Alternate between the two peak fill colours.
                fill_color = self.peak_fill_color1 if i % 2 == 0 else self.peak_fill_color2
                ax.fill_between(x_data[left_index:right_index],
                                corrected_data[left_index:right_index],
                                color=fill_color, alpha=self.peak_fill_alpha)
        
        ax.legend()
        self.canvas.draw()
     
    def update_peak_table(self):
        """Recalculate and update the peak table using non-overlapping integration regions based on valley minima."""
        if not self.last_analysis_data:
            return
        
        x_data = self.last_analysis_data["x_data"]
        corrected_data = self.last_analysis_data["corrected_data"]
        peaks = self.last_analysis_data["peaks"]
        peak_data = []
        peak_areas = []
        n_peaks = len(peaks)
        
        # Compute integration boundaries.
        boundaries = []
        if n_peaks == 0 or n_peaks == 1:
            boundaries = [x_data[0], x_data[-1]]
        else:
            boundaries.append(x_data[0])
            for i in range(n_peaks - 1):
                start_idx = peaks[i]
                end_idx = peaks[i + 1]
                if end_idx > start_idx:
                    region = corrected_data[start_idx:end_idx + 1]
                    valley_index = np.argmin(region) + start_idx
                    boundaries.append(x_data[valley_index])
                else:
                    boundaries.append(x_data[start_idx])
            boundaries.append(x_data[-1])
        
        # Compute peak properties within their unique integration regions.
        for i, peak in enumerate(peaks):
            left_bound = boundaries[i]
            right_bound = boundaries[i+1]
            left_index = np.searchsorted(x_data, left_bound)
            right_index = np.searchsorted(x_data, right_bound)
            # Compute area using the trapezoidal rule.
            area = np.trapz(corrected_data[left_index:right_index], x_data[left_index:right_index])
            height = corrected_data[peak]
            peak_data.append({
                "Peak #": i + 1,
                "Position (x)": x_data[peak],
                "Height": height,
                "Area": area
            })
            peak_areas.append(area)
        total_area = sum(peak_areas)
        for d in peak_data:
            d["% Area"] = (d["Area"] / total_area * 100) if total_area != 0 else 0
        
        import pandas as pd
        self.last_peak_df = pd.DataFrame(peak_data)
        
        self.peakTable.setRowCount(0)
        for row in self.last_peak_df.itertuples(index=False):
            rowPosition = self.peakTable.rowCount()
            self.peakTable.insertRow(rowPosition)
            self.peakTable.setItem(rowPosition, 0, QTableWidgetItem(str(row[0])))
            self.peakTable.setItem(rowPosition, 1, QTableWidgetItem(f"{row[1]:.3f}"))
            self.peakTable.setItem(rowPosition, 2, QTableWidgetItem(f"{row[2]:.3f}"))
            self.peakTable.setItem(rowPosition, 3, QTableWidgetItem(f"{row[3]:.3f}"))
            self.peakTable.setItem(rowPosition, 4, QTableWidgetItem(f"{row[4]:.3f}"))

    
    def toggle_interactive_mode(self):
        # just toggles interactive mode (peak delete)
        # originally this was meant to let you add peaks, but I couldn't make it work nicely
        self.interactive_mode = self.interactiveModeCheck.isChecked()
        if self.interactive_mode:
            self.cid_click = self.canvas.mpl_connect("button_press_event", self.on_click)
        else:
            if self.cid_click is not None:
                self.canvas.mpl_disconnect(self.cid_click)
                self.cid_click = None
    
    def on_click(self, event):
        # handles click events incase a peak is **almost** clicked
        if event.inaxes is None or not self.last_analysis_data:
            return
        tol_x = 0.02 * (self.last_analysis_data["x_data"].max() - self.last_analysis_data["x_data"].min())
        tol_y = 0.02 * (self.last_analysis_data["corrected_data"].max() - self.last_analysis_data["corrected_data"].min())
        x_click = event.xdata
        y_click = event.ydata
        remove_index = None
        for i, peak in enumerate(self.last_analysis_data["peaks"]):
            x_peak = self.last_analysis_data["x_data"][peak]
            y_peak = self.last_analysis_data["corrected_data"][peak]
            if abs(x_click - x_peak) < tol_x and abs(y_click - y_peak) < tol_y:
                remove_index = i
                break
        if remove_index is not None:
            del self.last_analysis_data["peaks"][remove_index]
            self.update_peak_table()
            self.update_plot()
    
    def save_analysis(self):
        # saves the current status to txt
        if not self.current_file:
            # throw an error incase someone clicks save before loading
            QMessageBox.warning(self, "No File Selected", "Please select a CSV file first.")
            return
        if self.last_peak_df is None:
            # throw an error rather than crash
            QMessageBox.warning(self, "No Analysis Data", "Run an analysis before saving.")
            return
        
        base_filename = os.path.splitext(os.path.basename(self.current_file))[0]
        output_png = f"{base_filename}_report.png"
        self.figure.savefig(output_png)
        
        output_txt = f"{base_filename}_report.txt"
        try:
            with open(output_txt, "w") as log_file:
                log_file.write("HPLC Data Analysis Report\n")
                log_file.write(f"Filename: {self.current_file}\n")
                log_file.write(f"Date/Time: {datetime.now()}\n\n")
                log_file.write("Fitting Options:\n")
                for key, value in self.last_parameters.items():
                    log_file.write(f"  {key}: {value}\n")
                log_file.write("\nPeak Table:\n")
                log_file.write(self.last_peak_df.to_string(index=False))
            QMessageBox.information(self, "Analysis Saved", f"Saved:\n{output_png}\n{output_txt}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save analysis:\n{e}")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
