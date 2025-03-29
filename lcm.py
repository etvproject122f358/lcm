import sys
import os
import re
import webbrowser
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.dates import DateFormatter
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.widgets import RectangleSelector
import mplcursors

from astropy.time import Time
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, get_body
import astropy.units as u
from astropy.visualization import quantity_support
from astroplan import moon_illumination
from scipy.special import erfcinv  # for Chauvenet's criterion

from PyQt5.QtCore import Qt, QSize, pyqtSlot, QLocale
from PyQt5.QtWidgets import (
    QApplication, QWidget, QPushButton, QLabel, QLineEdit, QTextEdit,
    QVBoxLayout, QHBoxLayout, QGroupBox, QComboBox, QSpinBox, QDoubleSpinBox,
    QCheckBox, QFileDialog, QDialog, QDialogButtonBox, QMessageBox, QTabWidget,
    QTableWidget, QTableWidgetItem, QSpacerItem, QSizePolicy
)
from PyQt5.QtGui import QIcon

from tess_util import (
    search_tess_sector, tess_analysis, get_tess_sector_data,
    phase_fold_and_normalize, phase_fold_manually
)
from simbad_util import create_simbad_link
from astro_app import AstroApp

# If you want to write file path manually, comment out below line
from dropbox_util import find_dropbox_folder

try:
    import dropbox_util

    dropbox_available = True
except ImportError:
    dropbox_available = False


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.data_file_path = None
        self.minima_file_path = None
        self.minima_file_loaded = False
        self.fig0 = None
        self.txt_window = None
        self.star_name = None
        self.sector_numbers = []
        self.tess_window = None
        self.simbad_window = None
        self.data_window = None  # Data window için referans
        self.tess_window = None
        self.initui()

    def initui(self):
        self.setWindowTitle("lcm")
        self.setWindowIcon(QIcon("images/star.ico"))

        self.data_select_btn = QPushButton()
        self.data_select_btn.setIcon(QIcon("images/lc1.png"))
        self.data_select_btn.setIconSize(QSize(50, 40))
        self.data_select_btn.clicked.connect(self.select_data_file)

        self.minima_select_btn = QPushButton()
        self.minima_select_btn.setIcon(QIcon("images/minima1.png"))
        self.minima_select_btn.setIconSize(QSize(50, 40))
        self.minima_select_btn.clicked.connect(self.select_minima_file)

        self.tess_btn = QPushButton()
        self.tess_btn.setIcon(QIcon("images/tess.png"))
        self.tess_btn.setIconSize(QSize(50, 40))
        self.tess_btn.clicked.connect(self.open_tess_window)

        self.simbad_btn = QPushButton()
        self.simbad_btn.setIcon(QIcon("images/simbad.png"))
        self.simbad_btn.setIconSize(QSize(50, 40))
        self.simbad_btn.clicked.connect(self.open_simbad_window)

        self.astro_app_btn = QPushButton()
        self.astro_app_btn.setIcon(QIcon("images/icon.ico"))
        self.astro_app_btn.setIconSize(QSize(50, 40))
        self.astro_app_btn.clicked.connect(self.open_astro_app)

        self.warning_label = QLabel()

        layout = QHBoxLayout()
        layout.addWidget(self.data_select_btn)
        layout.addWidget(self.minima_select_btn)
        layout.addWidget(self.tess_btn)
        layout.addWidget(self.simbad_btn)
        layout.addWidget(self.astro_app_btn)

        main_layout = QVBoxLayout()
        main_layout.addLayout(layout)
        main_layout.addWidget(self.warning_label)

        self.setLayout(main_layout)

        self.screen_resolution = QApplication.desktop().screenGeometry()
        self.screen_width = self.screen_resolution.width()
        self.screen_height = self.screen_resolution.height()

        self.obs = QWidget()
        self.txt_window = QWidget()

        if self.screen_width == 2560 and self.screen_height == 1440:
            self.setGeometry(300, 600, 200, 10)
        elif self.screen_width == 1920 and self.screen_height == 1080:
            self.setGeometry(110, 500, 200, 10)

    def open_astro_app(self):
        self.astro_app = AstroApp()
        self.astro_app.show()

        if self.star_name:
            self.astro_app.star_input.setText(self.star_name)

    def close_txt_window(self):
        if self.txt_window and self.txt_window.isVisible():
            self.txt_window.close()

    @pyqtSlot()
    def select_data_file(self):
        if dropbox_available:
            dropbox_path = dropbox_util.find_dropbox_folder()
            initial_dir = os.path.join(dropbox_path, "PCEB_project/observations/reduction_outputs")
        else:
            initial_dir = "C:/Users/furka/Documents/data"

        self.data_file_path, _ = QFileDialog.getOpenFileName(self, "Select Data File", initial_dir,
                                                             "Data Files (*.dat)")

        if self.data_file_path:
            self.warning_label.setText("Press Enter to save the plot")
            self.plot_data(self.data_file_path)
            self.info(self.data_file_path)

            self.star_name, _ = get_star_name_and_telescope(self.data_file_path)

            if self.tess_window and hasattr(self.tess_window, 'star_name_input'):
                self.tess_window.star_name_input.setText(self.star_name)

            if self.simbad_window and hasattr(self, 'simbad_star_name_input'):
                self.simbad_star_name_input.setText(self.star_name)

    def plot_data(self, file_path):
        if file_path:
            self.minima_file_loaded = True
            try:
                if self.fig0:
                    self.fig0.clf()
                    plt.close(self.fig0)
                    plt.close()
                    self.close_txt_window()

                data = pd.read_csv(file_path, sep='\s+', skiprows=1, names=["time", "flux", "error"])

                time = data["time"].values
                flux = data["flux"].values
                error = data["error"].values

                date = get_observation_date(file_path)
                star_name, telescope = get_star_name_and_telescope(file_path)
                filter_color, filter_name = get_filter(file_path)

                self.fig0 = plt.figure("Reduced Photometric Light Curve Data")
                manager = plt.get_current_fig_manager()
                manager.window.setWindowIcon(QIcon("images/star.ico"))

                if self.screen_width == 2560 and self.screen_height == 1440:
                    self.fig0.canvas.manager.window.setGeometry(600, 150, 1100, 700)
                elif self.screen_width == 1920 and self.screen_height == 1080:
                    self.fig0.canvas.manager.window.setGeometry(430, 30, 970, 540)

                self.ax = self.fig0.add_subplot(111)

                self.ax.errorbar(time, flux, yerr=error, marker=".", markersize=5, capsize=2, elinewidth=1.5, ls="None",
                                 color=filter_color)
                self.ax.set_xlabel('Time [BJD-TDB]')
                self.ax.set_ylabel('Relative Flux')
                self.ax.set_title(f'Date: {date}, Object: {star_name}, Telescope: {telescope}, Filter: {filter_name}')

                self.fig0.subplots_adjust(left=0.086, right=0.975, top=0.925, bottom=0.116)

                self.fig0.canvas.mpl_connect('key_press_event', self.on_key_press)
                self.ax.minorticks_on()
                self.ax.grid(True, which='major', linestyle='-', linewidth=0.1)
                self.fig0.canvas.draw()

            except Exception as e:
                print("Error reading data and plotting:", e)
                self.warning_label.setText("An error occurred while reading data and plotting.")
        else:
            self.warning_label.setText("You should select a .dat file first.")

    def on_key_press(self, event):
        if event.key == 'enter':
            self.save_plot()

    def save_plot(self):
        if not self.minima_file_loaded:
            QMessageBox.warning(self, "Warning", "Minima Data File is not selected.")
            return

        if not self.txt_window:
            QMessageBox.warning(self, "Warning", "Text window is not open. Can't save the plot.")
            return

        if self.fig0:
            if self.minima_file_path:
                output_filename = self.minima_file_path.replace('.txt', '-profile.png')
            else:
                QMessageBox.warning(self, "Warning", "Minima Data File is not selected.")
                return
        else:
            QMessageBox.warning(self, "Warning", "Data File is not selected.")
            return

        self.fig0.savefig(output_filename, dpi=300, bbox_inches='tight')
        QMessageBox.information(self, "Information", f"Plot saved as {output_filename}.")
        self.warning_label.clear()

    def info(self, file_path):
        global m_file_path
        if file_path:
            data_folder = os.path.dirname(file_path)
            m_file_path = None

            for file_name in os.listdir(data_folder):
                if file_name.endswith(".xls"):
                    m_file_path = os.path.join(data_folder, file_name)
                    break

            if not m_file_path:
                for file_name in os.listdir(data_folder):
                    if file_name.endswith(".tbl"):
                        m_file_path = os.path.join(data_folder, file_name)
                        break

            self.obs.close()

            self.obs = QWidget()

            if self.screen_width == 2560 and self.screen_height == 1440:
                self.obs.setGeometry(600, 900, 190, 400)
            elif self.screen_width == 1920 and self.screen_height == 1080:
                self.obs.setGeometry(430, 650, 200, 380)

            self.obs.setWindowTitle("Info")
            self.obs.setWindowIcon(QIcon("images/star.ico"))

            if m_file_path is not None:
                try:
                    layout = QVBoxLayout()
                    df = pd.read_csv(m_file_path, delimiter="\t")

                    exp_time = int(df["EXPTIME"].mean())
                    start_time = Time(df['JD_UTC'].iloc[0], format='jd', scale='tai').to_datetime()
                    end_time = Time(df['JD_UTC'].iloc[-1], format='jd', scale='tai').to_datetime()

                    start_time += pd.Timedelta(hours=3)
                    end_time += pd.Timedelta(hours=3)

                    rel_flux_SNR_T1 = [round(value) for value in df['rel_flux_SNR_T1'].iloc[:3]]
                    peak_T1 = round(df['Peak_T1'].iloc[:3])
                    sky = [round(value) for value in df['Sky/Pixel_T1'].iloc[:3]]

                    cursor = mplcursors.cursor(self.ax, hover=True, bindings={"toggle_visible": "h"})
                    cursor.connect("add", lambda sel: (sel.annotation.set_text(
                        "{:^12s}: {:^13g}\n{:^10s}: {:^10g}\n{:^10s}: {:^10g}".format(
                            "SNR",
                            round(df['rel_flux_SNR_T1'][sel.index]),
                            "Star Count",
                            round(df['Peak_T1'][sel.index]),
                            "BG Count",
                            round(df['Sky/Pixel_T1'][sel.index])
                        )
                    ),
                                                       sel.annotation.set_bbox({'facecolor': 'white', 'alpha': 0.7})
                    ))

                    self.warning_label.setText("<font color='blue'>Press 'h' for the hover on and off.</font>")

                    lab1 = QLabel("Exposure Time")
                    lab1.setAlignment(Qt.AlignCenter)
                    txt1 = QTextEdit(str(exp_time))
                    txt1.setAlignment(Qt.AlignCenter)
                    txt1.setReadOnly(True)
                    txt1.setFixedHeight(25)

                    lab2 = QLabel("Start Time")
                    lab2.setAlignment(Qt.AlignCenter)
                    txt2 = QTextEdit(start_time.strftime('%H:%M'))
                    txt2.setAlignment(Qt.AlignCenter)
                    txt2.setReadOnly(True)
                    txt2.setFixedHeight(25)

                    lab3 = QLabel("End Time")
                    lab3.setAlignment(Qt.AlignCenter)
                    txt3 = QTextEdit(end_time.strftime('%H:%M'))
                    txt3.setAlignment(Qt.AlignCenter)
                    txt3.setReadOnly(True)
                    txt3.setFixedHeight(25)

                    lab4 = QLabel("F3 SNR Values")
                    lab4.setAlignment(Qt.AlignCenter)
                    txt4 = QTextEdit(", ".join(map(str, rel_flux_SNR_T1)))
                    txt4.setAlignment(Qt.AlignCenter)
                    txt4.setReadOnly(True)
                    txt4.setFixedHeight(25)

                    lab5 = QLabel("F3 Count Values")
                    lab5.setAlignment(Qt.AlignCenter)
                    txt5 = QTextEdit(", ".join(map(str, peak_T1)))
                    txt5.setAlignment(Qt.AlignCenter)
                    txt5.setReadOnly(True)
                    txt5.setFixedHeight(25)

                    lab6 = QLabel("F3 BG Count Values")
                    lab6.setAlignment(Qt.AlignCenter)
                    txt6 = QTextEdit(", ".join(map(str, sky)))
                    txt6.setAlignment(Qt.AlignCenter)
                    txt6.setReadOnly(True)
                    txt6.setFixedHeight(25)

                    lay1 = QVBoxLayout()
                    lay1.addWidget(lab1)
                    lay1.addWidget(txt1)

                    lay2 = QVBoxLayout()
                    lay2.addWidget(lab2)
                    lay2.addWidget(txt2)

                    lay3 = QVBoxLayout()
                    lay3.addWidget(lab3)
                    lay3.addWidget(txt3)

                    lay4 = QVBoxLayout()
                    lay4.addWidget(lab4)
                    lay4.addWidget(txt4)

                    lay5 = QVBoxLayout()
                    lay5.addWidget(lab5)
                    lay5.addWidget(txt5)

                    lay6 = QVBoxLayout()
                    lay6.addWidget(lab6)
                    lay6.addWidget(txt6)

                    hlayout = QVBoxLayout()
                    hlayout.addLayout(lay1)
                    hlayout.addLayout(lay2)
                    hlayout.addLayout(lay3)
                    hlayout.addLayout(lay4)
                    hlayout.addLayout(lay5)
                    hlayout.addLayout(lay6)

                    self.ra_deg = df["RAOBJ2K"][0] * u.hour
                    self.dec_deg = df["DECOBJ2K"][0] * u.deg

                    obj = SkyCoord(ra=self.ra_deg, dec=self.dec_deg)

                    star_name, telescope = get_star_name_and_telescope(file_path)

                    if telescope == "T80" or telescope == "T35":
                        location = EarthLocation(lat=39.84 * u.deg, lon=32.782 * u.deg, height=1256.69 * u.m)
                    elif telescope == "T100":
                        location = EarthLocation(lat=36.825 * u.deg, lon=30.335 * u.deg, height=1256.69 * u.m)

                    utcoffset = 3 * u.hour
                    time = Time(df['JD_UTC'], format='jd', scale='tai') + utcoffset

                    new_time = time[0] - 12 * u.hour
                    midnight_new_time_str = new_time.datetime.strftime('%Y-%m-%d')
                    midnight = Time(midnight_new_time_str + ' 21:00:00')

                    delta_midnight = np.arange(-12 * 60, 13 * 60, 1) * u.minute
                    times = midnight + delta_midnight

                    local_times = times + utcoffset

                    frame = AltAz(obstime=times, location=location)
                    sunaltaz = get_sun(times).transform_to(frame)

                    moon = get_body("moon", times)
                    moonaltaz = moon.transform_to(frame)

                    objaltaz = obj.transform_to(frame)

                    illumination = moon_illumination(Time(df['JD_UTC'].iloc[0], format='jd', scale='tai'))

                    nearest_start_time_index = np.argmin(
                        np.abs((local_times.datetime - start_time).astype('timedelta64[s]')))
                    nearest_end_time_index = np.argmin(
                        np.abs((local_times.datetime - end_time).astype('timedelta64[s]')))

                    start_time_moon_alt = moonaltaz.alt.deg[nearest_start_time_index]
                    end_time_moon_alt = moonaltaz.alt.deg[nearest_end_time_index]

                    start_time_obj_alt = objaltaz.alt.deg[nearest_start_time_index]
                    end_time_obj_alt = objaltaz.alt.deg[nearest_end_time_index]

                    start_alt_difference = np.abs(start_time_obj_alt - start_time_moon_alt)
                    end_alt_difference = np.abs(end_time_obj_alt - end_time_moon_alt)

                    lab7 = QLabel("Init/End Moon-Object Separation")
                    lab7.setAlignment(Qt.AlignCenter)
                    txt7 = QTextEdit(
                        str(int(np.round(start_alt_difference))) + " - " + str(int(np.round(end_alt_difference))))

                    txt7.setAlignment(Qt.AlignCenter)
                    txt7.setReadOnly(True)
                    txt7.setFixedHeight(25)

                    lay7 = QVBoxLayout()
                    lay7.addWidget(lab7)
                    lay7.addWidget(txt7)

                    hlayout.addLayout(lay7)

                    layout.addLayout(hlayout)

                    quantity_support()

                    self.fig2 = plt.figure("Observation Visualization")
                    manager = plt.get_current_fig_manager()
                    manager.window.setWindowIcon(QIcon("images/star.ico"))
                    self.ax = self.fig2.add_subplot(111)

                    if self.screen_width == 2560 and self.screen_height == 1440:
                        self.fig2.canvas.manager.window.setGeometry(800, 900, 900, 400)
                    elif self.screen_width == 1920 and self.screen_height == 1080:
                        self.fig2.canvas.manager.window.setGeometry(660, 650, 740, 380)

                    self.fig2.subplots_adjust(left=0.086, right=0.981, top=0.898, bottom=0.159, hspace=0)
                    self.ax.plot(local_times.datetime, sunaltaz.alt, color="r", label="Sun")
                    self.ax.plot(local_times.datetime, moonaltaz.alt, color=[0.75] * 3, ls="--",
                                 label="Moon ({:.2f}%)".format(illumination * 100))

                    self.ax.scatter(
                        local_times.datetime,
                        objaltaz.alt,
                        c=objaltaz.az.value,
                        label=star_name,
                        lw=0,
                        s=8,
                        cmap="viridis",
                    )

                    self.ax.fill_between(
                        local_times.datetime,
                        0 * u.deg,
                        90 * u.deg,
                        sunaltaz.alt < -0 * u.deg,
                        color="k",
                        alpha=0.05,
                        zorder=0,
                    )
                    self.ax.fill_between(
                        local_times.datetime,
                        0 * u.deg,
                        90 * u.deg,
                        sunaltaz.alt < -6 * u.deg,
                        color="k",
                        alpha=0.15,
                        zorder=0,
                    )

                    self.ax.fill_between(
                        local_times.datetime,
                        0 * u.deg,
                        90 * u.deg,
                        sunaltaz.alt < -12 * u.deg,
                        color="k",
                        alpha=0.25,
                        zorder=0,
                    )

                    self.ax.fill_between(
                        local_times.datetime,
                        0 * u.deg,
                        90 * u.deg,
                        sunaltaz.alt < -18 * u.deg,
                        color="k",
                        alpha=0.35,
                        zorder=0,
                    )
                    self.ax.axvspan(start_time, end_time, color='yellow', alpha=0.3)
                    plt.colorbar(self.ax.collections[0], ax=self.ax, label="Azimuth [deg]")
                    plt.legend(loc="upper left")
                    plt.ylim(0 * u.deg, 90 * u.deg)
                    plt.xlabel("EEST")
                    plt.ylabel("Altitude [deg]")
                    plt.title("Visualization of Observation Duration Highlighted in Yellow")
                    self.ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))

                    self.obs.setLayout(layout)

                    self.obs.show()
                    self.fig2.show()
                    self.fig0.show()

                except Exception as e:
                    print("Error displaying observation info:", e)
                    self.warning_label.setText("An error occurred while displaying observation info.")
        else:
            print("No .xls file found in the folder of the .dat file.")
            QMessageBox.warning(self, "Warning", "No .xls file found in the folder of the .dat file.")

    def select_minima_file(self):
        minima_file_path, _ = QFileDialog.getOpenFileName(self, "Select Minima Data File", "", "Data Files (*.txt)")

        if minima_file_path:
            self.minima_file_loaded = True
            self.minima_file_path = minima_file_path

            self.txt_window = QWidget()
            self.txt_window.setWindowTitle("Minima Data")
            self.txt_window.setWindowIcon(QIcon("images/star.ico"))

            with open(minima_file_path, 'r') as file:
                lines = file.readlines()

                if len(lines) == 0 or any(len(line.strip().split()) != 3 for line in lines):
                    self.warning_label.setText("<font color='red'>The minima file is not in the expected format</font>")
                else:
                    if self.screen_width == 2560 and self.screen_height == 1440:
                        self.txt_window.setGeometry(1800, 550, 600, 200)
                    elif self.screen_width == 1920 and self.screen_height == 1080:
                        self.txt_window.setGeometry(1450, 440, 250, 200)

                    layout = QVBoxLayout()

                    lay1 = QVBoxLayout()
                    label1 = QLabel("Minima Time and Error [BJD-TDB]")
                    label1.setAlignment(Qt.AlignCenter)
                    txt_edit1 = QTextEdit()
                    txt_edit1.setReadOnly(True)
                    txt_edit1.setAlignment(Qt.AlignCenter)
                    txt_edit1.setFixedHeight(40)
                    copy_button1 = QPushButton("Copy")
                    copy_button1.clicked.connect(lambda: self.copy_text(txt_edit1.toPlainText()))
                    lay1.addWidget(label1)
                    lay1.addWidget(txt_edit1)
                    lay1.addWidget(copy_button1)

                    lay2 = QVBoxLayout()
                    label2 = QLabel("Minima Time in EEST")
                    label2.setAlignment(Qt.AlignCenter)
                    txt_edit2 = QTextEdit()
                    txt_edit2.setReadOnly(True)
                    txt_edit2.setAlignment(Qt.AlignCenter)
                    txt_edit2.setFixedHeight(40)
                    copy_button2 = QPushButton("Copy")
                    copy_button2.clicked.connect(lambda: self.copy_text(txt_edit2.toPlainText()))
                    lay2.addWidget(label2)
                    lay2.addWidget(txt_edit2)
                    lay2.addWidget(copy_button2)

                    lay3 = QVBoxLayout()
                    label3 = QLabel("Error in Seconds")
                    label3.setAlignment(Qt.AlignCenter)
                    txt_edit3 = QTextEdit()
                    txt_edit3.setReadOnly(True)
                    txt_edit3.setAlignment(Qt.AlignCenter)
                    txt_edit3.setFixedHeight(40)
                    copy_button3 = QPushButton("Copy")
                    copy_button3.clicked.connect(lambda: self.copy_text(txt_edit3.toPlainText()))
                    lay3.addWidget(label3)
                    lay3.addWidget(txt_edit3)
                    lay3.addWidget(copy_button3)

                    layout.addLayout(lay1)
                    layout.addLayout(lay2)
                    layout.addLayout(lay3)

                    on = get_star_name_and_telescope(minima_file_path)[0]

                    star_coord = SkyCoord.from_name(on)

                    for line in lines:
                        parts = line.strip().split()
                        if len(parts) == 3:
                            minima_time = float(parts[0])
                            minima_error = float(parts[2])

                            loc = EarthLocation(0, 0, 0, unit='m')
                            time = Time(minima_time, format='jd', scale='tdb', location=loc)
                            ltt_bary = time.light_travel_time(star_coord)
                            jd_tdb = time.tdb - ltt_bary
                            jd_utc = Time(jd_tdb, format='jd', scale='utc')
                            jd_utc += pd.Timedelta(hours=3)
                            dt = jd_utc.to_datetime()

                            txt_edit1.append(f"{minima_time:.6f} +/- {minima_error:.6f}")
                            txt_edit2.append(f"{dt.strftime('%H:%M')}")
                            txt_edit3.append(f"{minima_error * 86400:.2f}")

                    copy_button1 = QPushButton("Copy Minima Time and Error")
                    copy_button1.clicked.connect(lambda: self.copy_text(txt_edit1.toPlainText()))

                    copy_button2 = QPushButton("Copy Minima Time in EEST")
                    copy_button2.clicked.connect(lambda: self.copy_text(txt_edit2.toPlainText()))

                    copy_button3 = QPushButton("Copy Error in Seconds")
                    copy_button3.clicked.connect(lambda: self.copy_text(txt_edit3.toPlainText()))

                    self.txt_window.setLayout(layout)
                    self.txt_window.show()

            if not self.fig0:
                self.warning_label.setText("Data File is not selected.")
            else:
                if self.data_file_path:
                    minima_data_list = []

                    for line in lines:
                        minima_data = line.strip().split()
                        if len(minima_data) != 3:
                            self.warning_label.setText(
                                "<font color='red'>The minima file is not in the expected format</font>")
                        else:
                            minima_time = float(minima_data[0])
                            minima_type = int(minima_data[1])
                            minima_error = float(minima_data[2])
                            upper_limit = minima_time + minima_error
                            lower_limit = minima_time - minima_error
                            minima_data_list.append((minima_time, minima_type, minima_error, lower_limit, upper_limit))

                    for minima_time, minima_type, _, lower_limit, upper_limit in minima_data_list:
                        if self.fig0:
                            ax = self.fig0.gca()
                            if minima_type == 1:
                                ax.axvline(x=minima_time, color='black', linestyle='-', label='Min I')
                            elif minima_type == 2:
                                ax.axvline(x=minima_time, color='red', linestyle='-', label='Min II')

                        ax.axvline(x=lower_limit, color='grey', linestyle='--', lw=1)
                        ax.axvline(x=upper_limit, color='grey', linestyle='--', lw=1)

                        ax.set_xlabel('Time (BJD(TDB))')
                        ax.set_ylabel('Relative Flux')
                        ax.legend(loc="best")
                        self.warning_label.setText("<font color='green'>Press Enter to save the plot.</font>")
                        self.fig0.canvas.draw()
                        self.fig0.show()

                if not self.txt_window:
                    self.minima_file_loaded = False
        else:
            self.warning_label.setText("You should select a .txt file first.")

    def bjd_tdb_to_time_string(self, x, pos):
        jd = x
        star_coord = SkyCoord(ra=self.ra_deg, dec=self.dec_deg)
        loc = EarthLocation(0, 0, 0, unit='m')
        time = Time(jd, format='jd', scale='tdb', location=loc)
        ltt_bary = time.light_travel_time(star_coord)
        bjd_tdb = time.tdb + ltt_bary
        bjd_tdb += pd.Timedelta(hours=3)
        dt = bjd_tdb.to_datetime()
        return dt.strftime('%H:%M')

    def copy_text(self, text):
        QApplication.clipboard().setText(text)

    def open_tess_window(self):
        try:
            if self.tess_window is None or not self.tess_window.isVisible():
                star_name = self.star_name if self.star_name else ""
                self.tess_window = TESSWindow(self.screen_width, self.screen_height, star_name=star_name,
                                              sector_numbers=self.sector_numbers)
                self.tess_window.show()
            else:
                self.tess_window.raise_()
                self.tess_window.activateWindow()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while opening TESS window: {e}")

    def open_simbad_window(self):
        try:
            self.simbad_window = QWidget()
            self.simbad_window.setWindowTitle("SIMBAD Query")
            self.simbad_window.setWindowIcon(QIcon("images/star.ico"))

            self.simbad_star_name_input = QLineEdit()
            self.simbad_star_name_input.setPlaceholderText("Enter star name")

            self.simbad_search_btn = QPushButton("Open SIMBAD Link")
            self.simbad_search_btn.clicked.connect(
                lambda: self.open_simbad_link(self.simbad_star_name_input.text().strip()))

            layout = QVBoxLayout()
            layout.addWidget(self.simbad_star_name_input)
            layout.addWidget(self.simbad_search_btn)

            self.simbad_window.setLayout(layout)
            self.simbad_window.show()

            if self.star_name:
                self.simbad_star_name_input.setText(self.star_name)

            if self.screen_width == 2560 and self.screen_height == 1440:
                self.simbad_window.setGeometry(1800, 550, 250, 70)
            elif self.screen_width == 1920 and self.screen_height == 1080:
                self.simbad_window.setGeometry(610, 500, 250, 70)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")

    def open_simbad_link(self, star_name):
        if star_name:
            simbad_link = create_simbad_link(star_name)
            webbrowser.open(simbad_link)
        else:
            QMessageBox.warning(self, "Warning", "Please enter a star name.")


def get_observation_date(filename):
    match = re.search(r'(\d{8})', filename)
    if match:
        date_str = match.group(1)
        date = datetime.strptime(date_str, '%Y%m%d')
        date_n = date.strftime('%d %B %Y')
        return date_n
    return None


def get_star_name_and_telescope(filename):
    global star_name
    star_name = None

    if os.path.isfile(filename):
        filename = os.path.basename(filename)

    match = re.search(r'_(.*?)_(T35|T80|T100|ATA50|IST60)', filename)
    if match:
        star = match.group(1)
        telescope = match.group(2)
    else:
        star = filename
        telescope = None

    if '_' in star:
        star_name = star.replace('_', ' ')
    else:
        if len(star) > 3 and star[-3:] in constellation_abbr:
            star_name = star[:-3] + " " + star[-3:]
        elif star.startswith("GAIA"):
            parts = re.split(r'(\d+)', star, 1)
            star_name = f"{parts[0]} {parts[1][:3]} {parts[1][3:]}{parts[2]}"
        elif star.startswith("FBS"):
            parts = re.split(r'(\d+)', star, 1)
            if parts[2] == "":
                star_name = f"{parts[0]} {parts[1][:-3]}+{parts[1][-3:]}{parts[2]}"
            else:
                star_name = f"{parts[0]} {parts[1][:-3]}{parts[1][-3:]}{parts[2]}"
        elif star.startswith("KPD"):
            parts = re.split(r'(\d+)', star, 1)
            if parts[2] == "":
                star_name = f"{parts[0]} {parts[1][:-4]}+{parts[1][-4:]}{parts[2]}"
            else:
                star_name = f"{parts[0]} {parts[1][:-4]}{parts[1][-4:]}{parts[2]}"
        elif star.startswith("HS"):
            parts = re.split(r'(\d+)', star, 1)
            star_name = f"{parts[0]} {parts[1]}{parts[2]}"
        elif star.startswith("ATO"):
            parts = re.split(r'(\d+)', star, 1)
            star_name = f"{parts[0][0:3]} {parts[0][-1]}{parts[1]}{parts[2]}"
        elif star.startswith("GALEX"):
            # GALEX'ten sonra bir boşluk ekle
            star_name = f"GALEX {star[5:]}"
        elif star.startswith("ZTF"):
            # ZTF'ten sonra boşluk ekle ve J ile koordinatları koru
            star_name = f"ZTF {star[3:]}"
        elif star.startswith(("TIC", "KIC", "CHSS", "LAMOST")):
            parts = re.split(r'(\d+)', star, 1)
            star_name = f"{parts[0]} {parts[1]}{parts[2]}"
        else:
            star_name = star


    return star_name, telescope


constellation_abbr = [
    "And", "Ant", "Aps", "Aql", "Aqr", "Ara", "Ari", "Aur", "Boo", "Cae", "Cam", "Cap", "Car", "Cas", "Cen", "Cep",
    "Cet", "Cha", "Cir", "CMa", "CMi", "Cnc", "Col", "Com", "CrA", "CrB", "Crt", "Cru", "Crv", "CVn", "Cyg", "Del",
    "Dor", "Dra", "Equ", "Eri", "For", "Gem", "Gru", "Her", "Hor", "Hya", "Hyi", "Ind", "Lac", "Leo", "Lep", "Lib",
    "LMi", "Lup", "Lyn", "Lyr", "Men", "Mic", "Mon", "Mus", "Nor", "Oct", "Oph", "Ori", "Pav", "Peg", "Per", "Phe",
    "Pic", "PsA", "Psc", "Pup", "Pyx", "Ret", "Scl", "Sco", "Sct", "Ser", "Sex", "Sge", "Sgr", "Tau", "Tel", "TrA",
    "Tri", "Tuc", "UMa", "UMi", "Vel", "Vir", "Vol", "Vul"
]


def get_filter(filename):
    if os.path.isfile(filename):
        filename = os.path.basename(filename)

    match = re.search(r'_(.*?)_(LC|lc)', filename)
    prefix = match.group(1)
    filter_color = None
    filter_name = None

    if prefix.endswith('V'):
        filter_color = "green"
        filter_name = "V"
    elif prefix.endswith('g'):
        filter_color = "green"
        filter_name = "g"
    elif prefix.endswith('R'):
        filter_color = "red"
        filter_name = "R"
    elif prefix.endswith('r'):
        filter_color = "red"
        filter_name = "r"
    elif prefix.endswith('B'):
        filter_color = "blue"
        filter_name = "B"
    elif prefix.endswith('I'):
        filter_color = "saddlebrown"
        filter_name = "I"
    elif prefix.endswith('i'):
        filter_color = "brown"
        filter_name = "i"
    elif prefix.endswith('u'):
        filter_color = "purple"
        filter_name = "u"
    elif prefix.endswith('z'):
        filter_color = "gray"
        filter_name = "z"
    elif prefix.endswith('C'):
        filter_color = "black"
        filter_name = "Clear"

    return filter_color, filter_name


class TESSWindow(QWidget):
    def __init__(self, screen_width, screen_height, star_name=None, sector_numbers=None):
        super().__init__()
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.star_name = star_name
        self.sector_numbers = sector_numbers
        self.canvas = None
        self.canvas_outliers = None
        self.toolbar = None
        self.time = None
        self.flux = None
        self.flux_err = None
        self.selected_start_time = None
        self.selected_end_time = None
        self.canvas_detrending = None
        # Cleaned data attributes initialized to None
        self.cleaned_time = None
        self.cleaned_flux = None
        self.cleaned_flux_err = None

        # Initialize original data attributes
        self.original_time = None
        self.original_flux = None
        self.original_flux_err = None

        # Initialize current data attributes
        self.current_time = None
        self.current_flux = None
        self.current_flux_err = None
        self.automated_segments = []
        self.manual_segments = []

        self.segment_mode = "automate"

        self.init_ui()

    def init_ui(self):
        try:
            self.setWindowTitle("TESS Query")
            self.setWindowIcon(QIcon("images/star.ico"))

            tabs = QTabWidget()

            # ---- Query Tab ----
            query_tab = QWidget()
            self.query_layout = QVBoxLayout(query_tab)

            self.star_name_input = QLineEdit()
            self.star_name_input.setPlaceholderText("Enter star name")

            self.tess_search_btn = QPushButton("Search Sectors")
            self.tess_search_btn.clicked.connect(self.search_tess_sectors)

            self.display_data_btn = QPushButton("Query More")
            self.display_data_btn.clicked.connect(self.display_tess_sector_data)

            button_layout = QHBoxLayout()
            button_layout.addWidget(self.tess_search_btn)
            button_layout.addWidget(self.display_data_btn)

            self.query_layout.addWidget(self.star_name_input)
            self.query_layout.addLayout(button_layout)

            dropdown_layout = QHBoxLayout()
            self.sector_combo = QComboBox()
            self.sector_combo.setEditable(True)
            self.sector_combo.lineEdit().setPlaceholderText("Enter sector number")

            self.exptime_combo = QComboBox()
            self.exptime_combo.addItems(["120", "20"])

            self.flux_type_combo = QComboBox()
            self.flux_type_combo.addItems(["SAP FLUX", "PDCSAP FLUX"])

            dropdown_layout.addWidget(self.sector_combo)
            dropdown_layout.addWidget(self.exptime_combo)
            dropdown_layout.addWidget(self.flux_type_combo)

            self.query_layout.addLayout(dropdown_layout)

            download_layout = QHBoxLayout()
            self.quality_nan_checkbox = QCheckBox("Filter Data")  # Shorter label
            self.quality_nan_checkbox.setToolTip("Filter out data where QUALITY == 0 and drop NaN values.")
            self.quality_nan_checkbox.setChecked(True)

            spacer = QSpacerItem(0, 0, QSizePolicy.Fixed, QSizePolicy.Minimum)

            self.download_data_btn = QPushButton("Download and Plot Data")
            self.download_data_btn.clicked.connect(self.download_and_plot_data)

            download_layout.addWidget(self.quality_nan_checkbox)
            download_layout.addItem(spacer)
            download_layout.addWidget(self.download_data_btn)

            self.query_layout.addLayout(download_layout)

            light_curve_group = QGroupBox("Light Curve")
            light_curve_layout = QVBoxLayout()

            self.figure_query = Figure()
            self.canvas = FigureCanvas(self.figure_query)
            self.toolbar = NavigationToolbar(self.canvas, self)

            self.query_error_checkbox = QCheckBox("Display with Errors")
            self.query_error_checkbox.setChecked(True)
            self.query_error_checkbox.stateChanged.connect(self.update_query_plot)

            toolbar_layout = QHBoxLayout()
            toolbar_layout.addWidget(self.toolbar)
            toolbar_layout.addWidget(self.query_error_checkbox)
            toolbar_layout.setAlignment(Qt.AlignLeft)  # Align to the left side

            light_curve_layout.addWidget(self.canvas)
            light_curve_layout.addLayout(toolbar_layout)

            light_curve_group.setLayout(light_curve_layout)
            self.query_layout.addWidget(light_curve_group)

            query_button_layout = QHBoxLayout()

            self.detach_query_plot_button = QPushButton("Detached")
            self.detach_query_plot_button.clicked.connect(
                self.detach_query_plot)
            query_button_layout.addWidget(self.detach_query_plot_button)

            self.save_query_data_button = QPushButton("Save Data")
            self.save_query_data_button.clicked.connect(
                self.save_query_data)
            query_button_layout.addWidget(self.save_query_data_button)

            self.query_layout.addLayout(query_button_layout)

            query_tab.setLayout(self.query_layout)

            # ---- Outlier Tab ----
            outlier_tab = QWidget()
            outlier_layout = QVBoxLayout()

            # Manual Detection Group
            manual_detection_group = QGroupBox("Manual Detection")
            manual_layout = QVBoxLayout()

            self.select_region_button = QPushButton("Select Region")
            self.select_region_button.clicked.connect(self.enable_region_selection)
            manual_layout.addWidget(self.select_region_button)

            self.remove_selected_button = QPushButton("Remove Region")
            self.remove_selected_button.clicked.connect(self.remove_selected_interval)
            manual_layout.addWidget(self.remove_selected_button)

            manual_layout.addStretch()
            manual_detection_group.setLayout(manual_layout)
            manual_detection_group.setFixedHeight(100)

            # Automated Detection Group
            automated_detection_group = QGroupBox("Automated Detection")
            automated_layout = QVBoxLayout()

            detection_method_layout = QHBoxLayout()
            self.outlier_method_combo = QComboBox()
            self.outlier_method_combo.addItems(["Sigma Clipping", "Boxplot", "Chauvenet"])
            self.outlier_method_combo.currentIndexChanged.connect(self.toggle_sigma_input)
            self.outlier_method_combo.setToolTip("Select the method for detecting outliers")
            detection_method_layout.addWidget(self.outlier_method_combo)

            self.sigma_value_spinbox = QSpinBox()
            self.sigma_value_spinbox.setRange(1, 10)
            self.sigma_value_spinbox.setValue(3)
            self.sigma_value_spinbox.setFixedWidth(50)
            self.sigma_value_spinbox.setToolTip("Select the sigma value for clipping outliers")
            detection_method_layout.addWidget(self.sigma_value_spinbox)

            automated_layout.addLayout(detection_method_layout)

            apply_remove_layout = QHBoxLayout()
            self.apply_detection_button = QPushButton("Apply Detection")
            self.apply_detection_button.clicked.connect(self.apply_outlier_detection)
            apply_remove_layout.addWidget(self.apply_detection_button)

            self.remove_outliers_button = QPushButton("Remove Outliers")
            self.remove_outliers_button.clicked.connect(self.remove_outliers)
            apply_remove_layout.addWidget(self.remove_outliers_button)

            automated_layout.addLayout(apply_remove_layout)
            automated_layout.addStretch()
            automated_detection_group.setLayout(automated_layout)
            automated_detection_group.setFixedHeight(100)

            # Normalization Group
            normalization_group = QGroupBox("Normalization")
            normalization_layout = QVBoxLayout()

            self.apply_normalization_button = QPushButton("Apply Normalization")
            self.apply_normalization_button.clicked.connect(self.apply_normalization)
            normalization_layout.addWidget(self.apply_normalization_button)

            normalization_group.setLayout(normalization_layout)
            normalization_group.setFixedHeight(100)

            top_layout = QHBoxLayout()
            top_layout.addWidget(manual_detection_group)
            top_layout.addWidget(automated_detection_group)
            top_layout.addWidget(normalization_group)
            outlier_layout.addLayout(top_layout)

            restore_original_data_button = QPushButton("Restore Original Data")
            restore_original_data_button.clicked.connect(self.plot_original_data)
            outlier_layout.addWidget(restore_original_data_button)

            outlier_light_curve_group = QGroupBox("Light Curve")
            outlier_light_curve_layout = QVBoxLayout()

            self.figure_outliers = Figure()
            self.canvas_outliers = FigureCanvas(self.figure_outliers)
            self.toolbar_outliers = NavigationToolbar(self.canvas_outliers, self)
            self.outlier_error_checkbox = QCheckBox("Display with Errors")
            self.outlier_error_checkbox.setChecked(True)
            self.outlier_error_checkbox.stateChanged.connect(self.update_outlier_plot)

            toolbar_layout_outliers = QHBoxLayout()
            toolbar_layout_outliers.addWidget(self.toolbar_outliers)
            toolbar_layout_outliers.addWidget(self.outlier_error_checkbox)
            toolbar_layout_outliers.setAlignment(Qt.AlignLeft)

            outlier_light_curve_layout.addWidget(self.canvas_outliers)
            outlier_light_curve_layout.addLayout(toolbar_layout_outliers)
            outlier_light_curve_group.setLayout(outlier_light_curve_layout)

            outlier_layout.addWidget(outlier_light_curve_group)

            bottom_button_layout = QHBoxLayout()
            self.detach_plot_button = QPushButton("Detached")
            self.detach_plot_button.clicked.connect(self.detach_plot)
            self.detach_plot_button.setToolTip("Click to detach the plot for viewing separately.")
            bottom_button_layout.addWidget(self.detach_plot_button)

            self.save_cleaned_data_button = QPushButton("Save Cleaned Data")
            self.save_cleaned_data_button.clicked.connect(self.save_cleaned_data)
            self.save_cleaned_data_button.setToolTip("Click to save the cleaned data to a file.")
            bottom_button_layout.addWidget(self.save_cleaned_data_button)

            outlier_layout.addLayout(bottom_button_layout)
            outlier_tab.setLayout(outlier_layout)

            # ---- Detrending Tab ----
            detrending_tab = QWidget()
            detrending_layout = QVBoxLayout()

            # Manual Detection Group
            manual_detrending_group = QGroupBox("Manual Detection")
            manual_detrending_layout = QVBoxLayout()

            self.select_detrending_region_button = QPushButton("Select Region for Manual Detrending")
            self.select_detrending_region_button.clicked.connect(self.enable_detrending_region_selection)

            manual_buttons_layout = QHBoxLayout()
            self.show_manual_trend_button = QPushButton("Show Trend")
            self.show_manual_trend_button.clicked.connect(self.show_manual_trend)
            self.remove_manual_trend_button = QPushButton("Remove Trend")
            self.remove_manual_trend_button.clicked.connect(self.remove_manual_trend)

            manual_buttons_layout.addWidget(self.show_manual_trend_button)
            manual_buttons_layout.addWidget(self.remove_manual_trend_button)

            manual_detrending_layout.addWidget(self.select_detrending_region_button)
            manual_detrending_layout.addLayout(manual_buttons_layout)

            manual_detrending_group.setLayout(manual_detrending_layout)
            manual_detrending_group.setFixedHeight(100)

            # Polynomial Fit Group
            polynomial_fit_group = QGroupBox("Polynomial Fit")
            polynomial_fit_layout = QVBoxLayout()

            self.degree_combo = QComboBox()
            self.degree_combo.addItems(["1st Degree Polynomial", "2nd Degree Polynomial", "3rd Degree Polynomial", "4th Degree Polynomial"])
            self.degree_combo.setCurrentIndex(1)

            polynomial_fit_layout.addWidget(self.degree_combo)
            polynomial_fit_group.setLayout(polynomial_fit_layout)
            polynomial_fit_group.setFixedHeight(100)

            # Automated Detection Group
            automated_detection_group = QGroupBox("Automated Detection")
            automated_detection_layout = QVBoxLayout()
            threshold_detection_layout = QHBoxLayout()

            self.detect_gaps_button = QPushButton("Detect Gaps for Automated Detrending")
            self.detect_gaps_button.clicked.connect(self.detect_gaps)

            self.threshold_input = QDoubleSpinBox()
            self.threshold_input.setRange(0, 10)
            self.threshold_input.setValue(1.0)
            self.threshold_input.setSingleStep(0.1)
            self.threshold_input.setDecimals(2)
            self.threshold_input.setLocale(QLocale(QLocale.English))

            threshold_detection_layout.addWidget(self.detect_gaps_button)
            threshold_detection_layout.addWidget(self.threshold_input)

            automated_buttons_layout = QHBoxLayout()
            self.show_automated_trend_button = QPushButton("Show Trend")
            self.show_automated_trend_button.clicked.connect(self.show_automated_trend)
            self.remove_automated_trend_button = QPushButton("Remove Trend")
            self.remove_automated_trend_button.clicked.connect(self.remove_automated_trend)

            automated_buttons_layout.addWidget(self.show_automated_trend_button)
            automated_buttons_layout.addWidget(self.remove_automated_trend_button)

            automated_detection_layout.addLayout(threshold_detection_layout)
            automated_detection_layout.addLayout(automated_buttons_layout)

            automated_detection_group.setLayout(automated_detection_layout)
            automated_detection_group.setFixedHeight(100)

            detrending_controls_layout = QHBoxLayout()
            detrending_controls_layout.addWidget(manual_detrending_group)
            detrending_controls_layout.addWidget(polynomial_fit_group)
            detrending_controls_layout.addWidget(automated_detection_group)
            detrending_layout.addLayout(detrending_controls_layout)

            self.restore_button = QPushButton("Restore Data")
            self.restore_button.clicked.connect(self.restore_data)  # Connect to restore function
            detrending_layout.addWidget(self.restore_button)

            detrending_light_curve_group = QGroupBox("Light Curve")
            detrending_light_curve_layout = QVBoxLayout()

            self.figure_detrending = Figure()
            self.canvas_detrending = FigureCanvas(self.figure_detrending)
            self.toolbar_detrending = NavigationToolbar(self.canvas_detrending, self)
            self.detrending_error_checkbox = QCheckBox("Display with Errors")
            self.detrending_error_checkbox.setChecked(True)
            self.detrending_error_checkbox.stateChanged.connect(self.update_detrending_plot)

            toolbar_layout_detrending = QHBoxLayout()
            toolbar_layout_detrending.addWidget(self.toolbar_detrending)
            toolbar_layout_detrending.addWidget(self.detrending_error_checkbox)
            toolbar_layout_detrending.setAlignment(Qt.AlignLeft)

            detrending_light_curve_layout.addWidget(self.canvas_detrending)
            detrending_light_curve_layout.addLayout(toolbar_layout_detrending)
            detrending_light_curve_group.setLayout(detrending_light_curve_layout)
            detrending_layout.addWidget(detrending_light_curve_group)

            detrending_button_layout = QHBoxLayout()
            self.detach_detrending_plot_button = QPushButton("Detached")
            self.detach_detrending_plot_button.clicked.connect(self.detach_detrending_plot)
            detrending_button_layout.addWidget(self.detach_detrending_plot_button)

            self.save_detrending_data_button = QPushButton("Save")
            self.save_detrending_data_button.clicked.connect(self.save_detrending_data)
            detrending_button_layout.addWidget(self.save_detrending_data_button)

            detrending_layout.addLayout(detrending_button_layout)

            detrending_tab.setLayout(detrending_layout)

            # ---- Phase Fold Tab ----
            phase_tab = QWidget()
            phase_layout = QVBoxLayout()

            # Manual Phase Folding Group
            manual_group = QGroupBox("Manual Phase Folding")
            manual_layout = QVBoxLayout()

            manual_inputs_layout = QHBoxLayout()
            self.t0_input = QLineEdit()
            self.t0_input.setPlaceholderText("e.g., 2451350.024259")
            self.period_input = QLineEdit()
            self.period_input.setPlaceholderText("e.g., 0.2375")
            manual_inputs_layout.addWidget(QLabel("Minima Time:"))
            manual_inputs_layout.addWidget(self.t0_input)
            manual_inputs_layout.addWidget(QLabel("Period Value:"))
            manual_inputs_layout.addWidget(self.period_input)

            manual_button_layout = QHBoxLayout()
            manual_phase_button = QPushButton("Phase Fold Manually")
            manual_phase_button.clicked.connect(self.manual_phase_folding)
            manual_button_layout.addWidget(manual_phase_button)

            manual_layout.addLayout(manual_inputs_layout)
            manual_layout.addLayout(manual_button_layout)
            manual_group.setLayout(manual_layout)

            # Automate Phase Folding Group
            automate_group = QGroupBox("Automate Phase Folding")
            automate_layout = QVBoxLayout()

            automate_inputs_layout = QHBoxLayout()
            self.best_t0_display = QLineEdit("<Not calculated>")
            self.best_t0_display.setReadOnly(True)
            self.best_t0_display.setStyleSheet("background-color: #f0f0f0;")  # Grey background
            self.best_period_display = QLineEdit("<Not calculated>")
            self.best_period_display.setReadOnly(True)
            self.best_period_display.setStyleSheet("background-color: #f0f0f0;")  # Grey background
            automate_inputs_layout.addWidget(QLabel("Adjusted Minima Time:"))
            automate_inputs_layout.addWidget(self.best_t0_display)
            automate_inputs_layout.addWidget(QLabel("Best Period:"))
            automate_inputs_layout.addWidget(self.best_period_display)

            automate_button_layout = QHBoxLayout()
            self.period_dropdown = QComboBox()
            self.period_dropdown.setEditable(True)
            self.period_dropdown.addItem("Lomb-Scargle Period")
            self.period_dropdown.addItem("Period Estimate")
            automate_button = QPushButton("Automate Phase Folding")
            automate_button.clicked.connect(self.automate_phase_folding)
            automate_button_layout.addWidget(self.period_dropdown)
            automate_button_layout.addWidget(automate_button)

            automate_layout.addLayout(automate_inputs_layout)
            automate_layout.addLayout(automate_button_layout)
            automate_group.setLayout(automate_layout)

            phase_controls_layout = QHBoxLayout()
            manual_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            automate_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            phase_controls_layout.addWidget(manual_group)
            phase_controls_layout.addWidget(automate_group)

            phase_light_curve_group = QGroupBox("Light Curve")
            phase_light_curve_layout = QVBoxLayout()

            self.figure_phase = Figure()
            self.canvas_phase = FigureCanvas(self.figure_phase)
            self.toolbar_phase = NavigationToolbar(self.canvas_phase, self)

            self.phase_error_checkbox = QCheckBox("Display with Errors")
            self.phase_error_checkbox.setChecked(True)
            self.phase_error_checkbox.stateChanged.connect(self.update_phase_plot)

            toolbar_layout_phase = QHBoxLayout()
            toolbar_layout_phase.addWidget(self.toolbar_phase)
            toolbar_layout_phase.addWidget(self.phase_error_checkbox)
            toolbar_layout_phase.setAlignment(Qt.AlignLeft)

            phase_light_curve_layout.addWidget(self.canvas_phase)
            phase_light_curve_layout.addLayout(toolbar_layout_phase)
            phase_light_curve_group.setLayout(phase_light_curve_layout)

            phase_layout.addLayout(phase_controls_layout)
            phase_layout.addWidget(phase_light_curve_group)

            # Detached and Save Buttons
            phase_button_layout = QHBoxLayout()
            self.detach_phase_plot_button = QPushButton("Detached")
            self.detach_phase_plot_button.clicked.connect(self.detach_phase_plot)
            phase_button_layout.addWidget(self.detach_phase_plot_button)

            self.save_phase_data_button = QPushButton("Save Data")
            self.save_phase_data_button.clicked.connect(self.save_phase_data)
            phase_button_layout.addWidget(self.save_phase_data_button)
            phase_layout.addLayout(phase_button_layout)

            phase_tab.setLayout(phase_layout)

            tabs.addTab(query_tab, "Query")
            tabs.addTab(outlier_tab, "Outlier")
            tabs.addTab(detrending_tab, "Detrending")
            tabs.addTab(phase_tab, "Phase Fold")

            main_layout = QVBoxLayout()
            main_layout.addWidget(tabs)
            self.setLayout(main_layout)

            if self.star_name:
                self.star_name_input.setText(self.star_name)
            if self.sector_numbers:
                self.sector_combo.addItems(map(str, self.sector_numbers))
                self.sector_combo.setCurrentIndex(len(self.sector_numbers) - 1)

            self.setGeometry(300, 200, 1000, 800)  # Adjusted size

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")

    def update_query_plot(self):
        """Update the query plot based on the checkbox state."""
        self.plot_data_check(self.figure_query, self.canvas, self.query_error_checkbox.isChecked())

    def update_outlier_plot(self):
        """Update the outlier plot based on the checkbox state."""
        self.plot_data_check(self.figure_outliers, self.canvas_outliers, self.outlier_error_checkbox.isChecked())

    def update_detrending_plot(self):
        """Update the detrending plot based on the checkbox state."""
        self.plot_data_check(self.figure_detrending, self.canvas_detrending, self.detrending_error_checkbox.isChecked())


    def plot_data_check(self, figure, canvas, display_errors):
        """Plot data on the given figure and canvas with optional error bars."""
        ax = figure.gca()
        previous_title = ax.get_title()
        previous_xlabel = ax.get_xlabel()
        previous_ylabel = ax.get_ylabel()

        figure.clear()
        ax = figure.add_subplot(111)

        # Check if time and flux data are available
        if self.time is not None and self.flux is not None:
            if display_errors and self.flux_err is not None:
                ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue", markersize=5,
                            capsize=1, elinewidth=1.5, ls="none")
            else:
                ax.plot(self.time, self.flux, '.', color='dodgerblue', markersize=5)

            ax.set_title(previous_title)
            ax.set_xlabel(previous_xlabel)
            ax.set_ylabel(previous_ylabel)

        canvas.draw()


    def save_query_data(self):
        """Save the query data to a file in the chosen file format."""
        try:
            if self.time is None or self.flux is None or self.flux_err is None:
                QMessageBox.warning(self, "Warning", "No query data available to save.")
                return

            dialog = QDialog(self)
            dialog.setWindowTitle("Save Query Data Options")

            layout = QVBoxLayout(dialog)
            layout.addWidget(QLabel("Select file format to save:"))

            file_format_combo = QComboBox()
            file_format_combo.addItems(["XLSX", "TXT", "CSV"])
            layout.addWidget(file_format_combo)

            add_2457000_checkbox = QCheckBox("Add 2457000 to Time")
            add_2457000_checkbox.setChecked(True)  # Checked by default
            layout.addWidget(add_2457000_checkbox)

            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)

            if dialog.exec_() != QDialog.Accepted:
                return

            file_format = file_format_combo.currentText().lower()  # Convert to lowercase for consistency

            time_data = self.time + 2457000 if add_2457000_checkbox.isChecked() else self.time
            data = pd.DataFrame({
                "Time": time_data,
                "Flux": self.flux,
                "Error": self.flux_err
            })

            file_name, _ = QFileDialog.getSaveFileName(
                self, f"Save Query Data as {file_format.upper()}", "", f"{file_format.upper()} Files (*.{file_format})"
            )
            if not file_name:
                return

            if not file_name.endswith(f".{file_format}"):
                file_name += f".{file_format}"

            try:
                if file_format == "csv":
                    data.to_csv(file_name, index=False, float_format="%.10f")
                elif file_format == "xlsx":
                    with pd.ExcelWriter(file_name, engine="openpyxl", mode="w") as writer:
                        data.to_excel(writer, sheet_name="Query Data", index=False, float_format="%.10f")
                elif file_format == "txt":
                    np.savetxt(
                        file_name,
                        data.values,
                        fmt="%.10f",
                        header="\t".join(data.columns),
                        delimiter="\t",
                        comments=""
                    )

                QMessageBox.information(self, "Data Saved", f"Query data saved successfully to {file_name}.")

            except Exception as e:
                QMessageBox.critical(self, "Error",
                                     f"An error occurred while saving query data as {file_format.upper()}: {e}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}")


    def detach_query_plot(self):
        """Detach the query plot and display it in a separate window."""
        if self.time is None or self.flux is None or self.flux_err is None:
            QMessageBox.warning(self, "Warning", "No data available for plotting.")
            return

        fig, ax = plt.subplots()

        ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color='dodgerblue', markersize=5, capsize=1,
                    elinewidth=1.5, ls='none')
        ax.set_ylabel("Flux ($electrons^{-1}$)")
        ax.set_xlabel("Time - 2457000 ($day$)")
        ax.set_title("Detached Light Curve - Query Data")

        fig.tight_layout()
        fig.show()

    def detach_detrending_plot(self):
        """Detach the detrended plot and display it in a separate window."""
        try:
            if self.cleaned_time is not None and len(self.cleaned_time) > 0:
                time = self.cleaned_time
                flux = self.cleaned_flux
                flux_err = self.cleaned_flux_err
            elif self.time is not None and len(self.time) > 0:
                time = self.time
                flux = self.flux
                flux_err = self.flux_err
            else:
                QMessageBox.warning(self, "Warning", "No data available for plotting.")
                return

            fig, ax = plt.subplots()

            ax.errorbar(time, flux, yerr=flux_err, fmt=".", color="dodgerblue", markersize=5, capsize=1,
                        elinewidth=1.5, ls="none")
            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrended Light Curve - Detached View")
            fig.tight_layout()
            fig.show()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while detaching the detrending plot: {e}")

    def toggle_sigma_input(self):
        method = self.outlier_method_combo.currentText()
        if method == "Sigma Clipping":
            self.sigma_value_spinbox.setVisible(True)
        else:
            self.sigma_value_spinbox.setVisible(False)

    def apply_normalization(self):
        """Apply normalization to the current data based on the median value and update the plots."""
        try:
            flux_median = np.median(self.flux)

            normalized_flux = self.flux / flux_median
            normalized_flux_err = self.flux_err / flux_median  # Normalize errors as well

            self.time = self.time
            self.flux = normalized_flux
            self.flux_err = normalized_flux_err

            self.update_outlier_data(self.time, self.flux, self.flux_err)

            self.plot_data_in_detrending()

            QMessageBox.information(self, "Normalization Applied",
                                    "Data has been normalized based on median successfully.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred during normalization: {e}")

    def download_and_plot_data(self):
        try:
            self.download_data_btn.setText("Downloading and plotting...")
            self.download_data_btn.setStyleSheet("color: red;")
            QApplication.processEvents()  # Refresh UI immediately

            star_name = self.star_name_input.text().strip() or self.star_name
            sector = self.sector_combo.currentText()
            exptime = int(self.exptime_combo.currentText())
            flux_type = self.flux_type_combo.currentText().replace(" ", "_")

            if not star_name:
                QMessageBox.warning(self, "Warning", "Please enter a star name.")
                self.download_data_btn.setText("Download and Plot Data")
                self.download_data_btn.setStyleSheet("color: black;")
                return

            if not sector.isdigit():
                QMessageBox.warning(self, "Warning", "Sector number must be an integer.")
                self.download_data_btn.setText("Download and Plot Data")
                self.download_data_btn.setStyleSheet("color: black;")
                return

            sector = int(sector)

            apply_filter = self.quality_nan_checkbox.isChecked()

            time, flux, flux_err, quality = tess_analysis(
                star_name=star_name,
                sector=sector,
                exptime=exptime,
                flux_type=flux_type,
                parent=self,
                apply_filter=apply_filter
            )

            if time is None:
                self.download_data_btn.setText("Download and Plot Data")
                self.download_data_btn.setStyleSheet("color: black;")
                return

            self.time = time
            self.flux = flux
            self.flux_err = flux_err

            self.original_time = self.time.copy()
            self.original_flux = self.flux.copy()
            self.original_flux_err = self.flux_err.copy()

            self.download_data_btn.setText("Download and Plot Data")
            self.download_data_btn.setStyleSheet("color: black;")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while downloading and plotting data: {e}")
            self.download_data_btn.setText("Download and Plot Data")
            self.download_data_btn.setStyleSheet("color: black;")

    def on_select(self, eclick, erelease):
        """Callback when a region is selected based on both time and flux values."""
        self.selected_start_time = eclick.xdata
        self.selected_end_time = erelease.xdata
        self.selected_start_flux = eclick.ydata
        self.selected_end_flux = erelease.ydata

        ax = self.figure_outliers.gca()
        xlim = ax.get_xlim()  # Get current x-axis limits
        ylim = ax.get_ylim()  # Get current y-axis limits

        ax.axvspan(self.selected_start_time, self.selected_end_time, color='yellow', alpha=0.3)
        ax.axhspan(self.selected_start_flux, self.selected_end_flux, color='yellow', alpha=0.3)
        self.canvas_outliers.draw()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.canvas_outliers.draw()

    def remove_selected_data(self, time, flux, flux_err, selected_range):
        """Remove data within the selected time range for time, flux, and flux_err."""
        start_time, end_time = selected_range

        if start_time > end_time:
            start_time, end_time = end_time, start_time

        mask = (time < start_time) | (time > end_time)

        time_filtered = time[mask]
        flux_filtered = flux[mask]
        flux_err_filtered = flux_err[mask]  # Flux errorları da filtreliyoruz

        return time_filtered, flux_filtered, flux_err_filtered

    def remove_selected_time_flux_data(self, time, flux, flux_err, selected_range_time, selected_range_flux):
        """Remove data within the selected time and flux range and return the filtered data."""
        start_time, end_time = selected_range_time
        start_flux, end_flux = selected_range_flux

        if start_time > end_time:
            start_time, end_time = end_time, start_time
        if start_flux > end_flux:
            start_flux, end_flux = end_flux, start_flux

        mask = ((time < start_time) | (time > end_time)) | ((flux < start_flux) | (flux > end_flux))

        time_filtered = time[mask]
        flux_filtered = flux[mask]
        flux_err_filtered = flux_err[mask]

        return time_filtered, flux_filtered, flux_err_filtered

    def update_plot(self, time, flux):
        if self.canvas is not None:
            self.plot_layout.removeWidget(self.canvas)
            self.canvas.deleteLater()

        # Yeni canvas oluştur
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.plot_layout.addWidget(self.canvas)

        ax = self.figure.add_subplot(111)
        ax.plot(time, flux, '.', color='blue', markersize=5)
        ax.set_ylabel("Flux")
        ax.set_xlabel("Time")
        self.canvas.draw()

    def save_cleaned_data(self):
        try:
            if self.time is None or self.flux is None or self.flux_err is None:
                QMessageBox.warning(self, "Warning", "No cleaned data available to save.")
                return

            dialog = QDialog(self)
            dialog.setWindowTitle("Save Options")

            layout = QVBoxLayout(dialog)
            layout.addWidget(QLabel("Select file format to save:"))

            self.file_format_combo = QComboBox()
            self.file_format_combo.addItems(["XLSX", "TXT", "CSV"])
            layout.addWidget(self.file_format_combo)

            self.add_2457000_checkbox = QCheckBox("Add 2457000 to Time")
            self.add_2457000_checkbox.setChecked(True)
            layout.addWidget(self.add_2457000_checkbox)

            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)

            if dialog.exec_() != QDialog.Accepted:
                return

            file_format = self.file_format_combo.currentText().lower()

            time_data = self.time + 2457000 if self.add_2457000_checkbox.isChecked() else self.time
            data = pd.DataFrame({
                "Time": time_data,
                "Flux": self.flux,
                "FluxErr": self.flux_err
            })

            file_name, _ = QFileDialog.getSaveFileName(
                self, f"Save Data as {file_format.upper()}", "", f"{file_format.upper()} Files (*.{file_format})"
            )
            if not file_name:
                return

            if not file_name.endswith(f".{file_format}"):
                file_name += f".{file_format}"

            try:
                if file_format == "csv":
                    data.to_csv(file_name, index=False, float_format="%.10f")
                elif file_format == "xlsx":
                    with pd.ExcelWriter(file_name, engine="openpyxl", mode="w") as writer:
                        data.to_excel(writer, sheet_name="Cleaned Data", index=False, float_format="%.10f")
                elif file_format == "txt":
                    np.savetxt(
                        file_name,
                        data.values,
                        fmt="%.10f",
                        header="\t".join(data.columns),
                        delimiter="\t",
                        comments=""
                    )

                QMessageBox.information(self, "Data Saved", f"Cleaned data saved successfully to {file_name}.")

            except Exception as e:
                QMessageBox.critical(self, "Error",
                                     f"An error occurred while saving data as {file_format.upper()}: {e}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}")

    def plot_data_in_outliers(self):
        """Plot data in the Outliers tab based on the current state of self.time, self.flux, and self.flux_err."""
        ax = self.figure_outliers.gca()
        ax.clear()

        if self.outlier_error_checkbox.isChecked() and self.flux_err is not None:
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")
        else:
            ax.plot(self.time, self.flux, ".", color="dodgerblue", markersize=5)

        ax.set_ylabel("Flux ($electrons^{-1}$)")
        ax.set_xlabel("Time - 2457000 ($day$)")
        ax.set_title("Light Curve without Outliers")

        self.canvas_outliers.draw()

    def restore_data_detrending(self):
        """Restore data in the Detrending tab to match the latest cleaned data from the Outlier tab."""
        if self.cleaned_time is not None and self.cleaned_flux is not None and self.cleaned_flux_err is not None:
            self.time = self.cleaned_time
            self.flux = self.cleaned_flux
            self.flux_err = self.cleaned_flux_err
            data_type = "Cleaned"
        else:
            self.time = self.original_time
            self.flux = self.original_flux
            self.flux_err = self.original_flux_err
            data_type = "Original"

        self.plot_data_in_detrending()
        QMessageBox.information(self, "Data Restored", f"Detrending data has been restored to {data_type} data.")

    def remove_selected_interval(self):
        """Remove the selected time and flux range from the data and update the plots."""
        try:
            if hasattr(self, 'selected_start_time') and hasattr(self, 'selected_end_time') and \
                    hasattr(self, 'selected_start_flux') and hasattr(self, 'selected_end_flux'):

                selected_range_time = (self.selected_start_time, self.selected_end_time)
                selected_range_flux = (self.selected_start_flux, self.selected_end_flux)

                if selected_range_time[0] is None or selected_range_time[1] is None or \
                        selected_range_flux[0] is None or selected_range_flux[1] is None:
                    QMessageBox.warning(self, "Warning", "No valid region selected.")
                    return

                time_filtered, flux_filtered, flux_err_filtered = self.remove_selected_time_flux_data(
                    self.time, self.flux, self.flux_err, selected_range_time, selected_range_flux)

                self.time = time_filtered
                self.flux = flux_filtered
                self.flux_err = flux_err_filtered

                ax = self.figure_outliers.gca()
                for patch in ax.patches:
                    patch.remove()
                for collection in ax.collections:
                    collection.remove()
                self.canvas_outliers.draw()

                self.plot_data_in_outliers()

                self.selected_start_time = None
                self.selected_end_time = None
                self.selected_start_flux = None
                self.selected_end_flux = None

                if hasattr(self, 'selector'):
                    self.selector.set_active(False)
                    del self.selector  # Delete the selector to remove the highlight

                QMessageBox.information(self, "Region Removed",
                                        f"Data between time {selected_range_time[0]:.2f} and {selected_range_time[1]:.2f} "
                                        f"and flux {selected_range_flux[0]:.2f} and {selected_range_flux[1]:.2f} has been removed.")

            else:
                QMessageBox.warning(self, "Warning", "No valid region selected.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")

    def get_selected_range_from_plot(self):
        """Use the RectangleSelector to select a region and return the start and end times of selection."""

        def onselect(eclick, erelease):
            """Callback to get the selected range when the selection is made."""
            self.selected_start_time = eclick.xdata
            self.selected_end_time = erelease.xdata
            print(f"Selected region: {self.selected_start_time:.2f} to {self.selected_end_time:.2f}")

            ax = self.figure_outliers.gca()
            ax.axvspan(self.selected_start_time, self.selected_end_time, color='yellow', alpha=0.3)
            self.canvas_outliers.draw()

        ax = self.figure_outliers.gca()

        self.selector = RectangleSelector(ax, onselect, useblit=True, interactive=True,
                                          button=[1], minspanx=0.1, minspany=0.1, spancoords='data')

        self.canvas_outliers.draw()

        return getattr(self, 'selected_start_time', None), getattr(self, 'selected_end_time', None)

    def apply_sigma_clipping(self, sigma_value):
        """Apply sigma clipping using the provided sigma value."""
        try:
            if not hasattr(self, 'current_flux') or self.current_flux is None:
                self.current_flux = self.original_flux.copy()
                self.current_time = self.original_time.copy()
                self.current_flux_err = self.original_flux_err.copy()

            overall_mean = np.mean(self.current_flux)
            overall_std = np.std(self.current_flux)

            lower_bound = overall_mean - sigma_value * overall_std
            upper_bound = overall_mean + sigma_value * overall_std

            if hasattr(self, 'lower_bound_line') and hasattr(self, 'upper_bound_line'):
                self.lower_bound_line.remove()
                self.upper_bound_line.remove()

            ax = self.figure_outliers.gca()
            self.lower_bound_line = ax.axhline(lower_bound, color='red', linestyle='--',
                                               label=f'{sigma_value} sigma lower bound')
            self.upper_bound_line = ax.axhline(upper_bound, color='green', linestyle='--',
                                               label=f'{sigma_value} sigma upper bound')

            self.canvas_outliers.draw()

            QMessageBox.information(self, "Sigma Clipping",
                                    f"Applied {sigma_value}-sigma clipping with mean={overall_mean:.2f} and std={overall_std:.2f}.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred during sigma clipping: {e}")

    def remove_sigma_clipped_data(self):
        """Sigma clipping ile belirlenen sınırlar dışındaki verileri sil ve grafikleri güncelle."""
        try:
            if hasattr(self, 'lower_bound_line') and hasattr(self, 'upper_bound_line'):
                lower_bound = self.lower_bound_line.get_ydata()[0]
                upper_bound = self.upper_bound_line.get_ydata()[0]

                mask = (self.current_flux >= lower_bound) & (self.current_flux <= upper_bound)

                if np.any(~mask):
                    self.current_time = self.current_time[mask]
                    self.current_flux = self.current_flux[mask]
                    self.current_flux_err = self.current_flux_err[mask]

                    self.time = self.current_time
                    self.flux = self.current_flux
                    self.flux_err = self.current_flux_err

                    self.plot_data_in_outliers()

                    self.plot_data_in_detrending()

                    QMessageBox.information(self, "Outliers Removed",
                                            f"Data outside of the sigma clipping bounds has been removed and graph updated.")
                else:
                    QMessageBox.information(self, "Sigma Clipping", "No data points were removed.")
            else:
                QMessageBox.warning(self, "Warning", "No sigma clipping has been applied yet.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while removing sigma clipped data: {e}")

    def apply_outlier_detection(self):
        method = self.outlier_method_combo.currentText()
        if method == "Sigma Clipping":
            sigma_val = self.sigma_value_spinbox.value()
            self.apply_sigma_clipping(sigma_val)
        elif method == "Boxplot":
            self.apply_boxplot_outliers()
        elif method == "Chauvenet":
            self.apply_chauvenet_outliers()

    def remove_outliers(self):
        method = self.outlier_method_combo.currentText()
        if method == "Sigma Clipping":
            self.remove_sigma_clipped_data()
        elif method == "Boxplot":
            self.remove_boxplot_outliers()
        elif method == "Chauvenet":
            self.remove_chauvenet_outliers()

    def detach_plot(self):
        try:
            if self.time is None or self.flux is None or self.flux_err is None:
                QMessageBox.warning(self, "Warning", "No data available for plotting.")
                return

            fig1, ax = plt.subplots()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue", markersize=5, capsize=1, elinewidth=1.5, ls="none")
            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detached Light Curve")
            fig1.tight_layout()
            fig1.show()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while detaching plot: {e}")



    def apply_boxplot_outliers(self):
        try:
            if not hasattr(self, 'current_flux') or self.current_flux is None:
                self.current_flux = self.original_flux.copy()
                self.current_time = self.original_time.copy()
                self.current_flux_err = self.original_flux_err.copy()

            q1 = np.percentile(self.current_flux, 25)
            q3 = np.percentile(self.current_flux, 75)
            iqr = q3 - q1

            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr

            if hasattr(self, 'boxplot_lower_line') and hasattr(self, 'boxplot_upper_line'):
                self.boxplot_lower_line.remove()
                self.boxplot_upper_line.remove()

            ax = self.figure_outliers.gca()
            self.boxplot_lower_line = ax.axhline(lower_bound, color='red', linestyle='--',
                                                 label="Boxplot lower bound")
            self.boxplot_upper_line = ax.axhline(upper_bound, color='green', linestyle='--',
                                                 label="Boxplot upper bound")
            self.canvas_outliers.draw()

            QMessageBox.information(self, "Boxplot Outliers",
                                    f"Boxplot thresholds applied.\nLower: {lower_bound:.3f}, Upper: {upper_bound:.3f}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"An error occurred while applying Boxplot outliers: {e}")

    def remove_boxplot_outliers(self):
        try:
            if hasattr(self, 'boxplot_lower_line') and hasattr(self, 'boxplot_upper_line'):
                lower_bound = self.boxplot_lower_line.get_ydata()[0]
                upper_bound = self.boxplot_upper_line.get_ydata()[0]

                mask = (self.current_flux >= lower_bound) & (self.current_flux <= upper_bound)

                if np.any(~mask):
                    self.current_time = self.current_time[mask]
                    self.current_flux = self.current_flux[mask]
                    self.current_flux_err = self.current_flux_err[mask]

                    self.time = self.current_time
                    self.flux = self.current_flux
                    self.flux_err = self.current_flux_err

                    self.plot_data_in_outliers()

                    self.plot_data_in_detrending()

                    QMessageBox.information(self, "Outliers Removed",
                                            "Data outside of Boxplot bounds has been removed.")
                else:
                    QMessageBox.information(self, "Boxplot Outliers",
                                            "No data points were removed. (All within Boxplot bounds.)")
            else:
                QMessageBox.warning(self, "Warning",
                                    "Boxplot thresholds have not been applied yet.")
        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"An error occurred while removing Boxplot outliers: {e}")


    def apply_chauvenet_outliers(self):
        try:
            if not hasattr(self, 'current_flux') or self.current_flux is None:
                self.current_flux = self.original_flux.copy()
                self.current_time = self.original_time.copy()
                self.current_flux_err = self.original_flux_err.copy()

            mean_flux = np.mean(self.current_flux)
            std_flux = np.std(self.current_flux)
            n = len(self.current_flux)

            probability = 1.0 / (2 * n)
            critical_value = np.sqrt(2) * erfcinv(2 * probability)


            lower_bound = mean_flux - critical_value * std_flux
            upper_bound = mean_flux + critical_value * std_flux

            if hasattr(self, 'chauvenet_lower_line') and hasattr(self, 'chauvenet_upper_line'):
                self.chauvenet_lower_line.remove()
                self.chauvenet_upper_line.remove()

            ax = self.figure_outliers.gca()
            self.chauvenet_lower_line = ax.axhline(lower_bound, color='red', linestyle='--',
                                                   label="Chauvenet lower bound")
            self.chauvenet_upper_line = ax.axhline(upper_bound, color='green', linestyle='--',
                                                   label="Chauvenet upper bound")
            self.canvas_outliers.draw()

            QMessageBox.information(self, "Chauvenet Outliers",
                                    f"Chauvenet limit applied.\nLower: {lower_bound:.3f}, Upper: {upper_bound:.3f}")

        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"An error occurred while applying Chauvenet outliers: {e}")

    def remove_chauvenet_outliers(self):
        try:
            if hasattr(self, 'chauvenet_lower_line') and hasattr(self, 'chauvenet_upper_line'):
                lower_bound = self.chauvenet_lower_line.get_ydata()[0]
                upper_bound = self.chauvenet_upper_line.get_ydata()[0]

                mask = (self.current_flux >= lower_bound) & (self.current_flux <= upper_bound)

                if np.any(~mask):
                    self.current_time = self.current_time[mask]
                    self.current_flux = self.current_flux[mask]
                    self.current_flux_err = self.current_flux_err[mask]

                    self.time = self.current_time
                    self.flux = self.current_flux
                    self.flux_err = self.current_flux_err

                    self.plot_data_in_outliers()
                    self.plot_data_in_detrending()

                    QMessageBox.information(self, "Outliers Removed",
                                            "Data outside of Chauvenet bounds has been removed.")
                else:
                    QMessageBox.information(self, "Chauvenet Outliers",
                                            "No data points were removed. (All within Chauvenet bounds.)")
            else:
                QMessageBox.warning(self, "Warning",
                                    "Chauvenet thresholds have not been applied yet.")
        except Exception as e:
            QMessageBox.critical(self, "Error",
                                 f"An error occurred while removing Chauvenet outliers: {e}")


    def plot_original_data(self):
        """Redraw the original data and update the plots."""
        try:
            # Update the data to original
            self.time = self.original_time.copy()
            self.flux = self.original_flux.copy()
            self.flux_err = self.original_flux_err.copy()

            # Redraw the Outlier tab plot
            ax = self.figure_outliers.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")
            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Original Light Curve from Query Tab")
            self.canvas_outliers.draw()

            self.plot_data_in_detrending()

            QMessageBox.information(self, "Data Restored",
                                    "The original data has been restored in both Outlier and Detrending sections.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while restoring the original data: {e}")

    def search_tess_sectors(self):
        star_name = self.star_name_input.text().strip() or self.star_name
        if not star_name:
            QMessageBox.warning(self, "Warning", "Please enter a star name.")
            return

        self.tess_search_btn.setText("Searching...")
        self.tess_search_btn.setStyleSheet("color: red;")
        QApplication.processEvents()  # Update UI immediately

        message, sector_numbers = search_tess_sector(star_name)

        self.tess_search_btn.setText("Search Sectors")
        self.tess_search_btn.setStyleSheet("color: black;")

        QMessageBox.information(self, "Information", message)

        if sector_numbers:
            self.sector_numbers = sector_numbers
            self.sector_combo.clear()  # Clear the current list
            self.sector_combo.addItems(map(str, self.sector_numbers))
            self.sector_combo.setCurrentIndex(len(self.sector_numbers) - 1)  # Focus on the last sector

    def display_tess_sector_data(self):
        star_name = self.star_name_input.text().strip() or self.star_name
        if not star_name:
            QMessageBox.warning(self, "Warning", "Please enter a star name.")
            return

        self.display_data_btn.setText("Querying...")
        self.display_data_btn.setStyleSheet("color: red;")
        QApplication.processEvents()  # Update UI immediately

        lc_data, ffi_data = get_tess_sector_data(star_name)

        self.display_data_btn.setText("Query More")
        self.display_data_btn.setStyleSheet("color: black;")

        if not lc_data and not ffi_data:
            QMessageBox.information(self, "Information", "No TESS data found for this star.")
            return

        self.data_window = QWidget()
        self.data_window.setWindowTitle("TESS Sector Data")
        self.data_window.setWindowIcon(QIcon("images/star.ico"))

        star_name_label = QLabel(f"<b>Star Name: <font color='blue'>{star_name}</font></b>")
        star_name_label.setAlignment(Qt.AlignCenter)

        unique_sectors = set(item[0] for item in lc_data)
        sector_count_label = QLabel(f"<b>Total number of sectors: <font color='red'>{len(unique_sectors)}</font></b>")
        sector_count_label.setAlignment(Qt.AlignCenter)

        timeseries_label = QLabel("<b>Timeseries Data</b>")
        timeseries_label.setAlignment(Qt.AlignCenter)

        lc_table = QTableWidget()
        lc_table.setRowCount(len(lc_data))
        lc_table.setColumnCount(2)
        lc_table.setHorizontalHeaderLabels(["Sector Number", "Exposure Time"])

        for row, (sector, exptime) in enumerate(lc_data):
            sector_item = QTableWidgetItem(str(sector))
            sector_item.setTextAlignment(Qt.AlignCenter)
            exptime_item = QTableWidgetItem(str(exptime))
            exptime_item.setTextAlignment(Qt.AlignCenter)
            lc_table.setItem(row, 0, sector_item)
            lc_table.setItem(row, 1, exptime_item)

        ffi_label = QLabel("<b>Full Frame Image Data</b>")
        ffi_label.setAlignment(Qt.AlignCenter)

        ffi_table = QTableWidget()
        ffi_table.setRowCount(len(ffi_data))
        ffi_table.setColumnCount(2)
        ffi_table.setHorizontalHeaderLabels(["Sector Number", "Exposure Time"])

        for row, (sector, data_type, exptime) in enumerate(ffi_data):
            sector_item = QTableWidgetItem(str(sector))
            sector_item.setTextAlignment(Qt.AlignCenter)
            exptime_item = QTableWidgetItem(str(exptime))
            exptime_item.setTextAlignment(Qt.AlignCenter)
            ffi_table.setItem(row, 0, sector_item)
            ffi_table.setItem(row, 1, exptime_item)

        # Layouts
        timeseries_layout = QVBoxLayout()
        timeseries_layout.addSpacing(10)
        timeseries_layout.addWidget(timeseries_label)
        timeseries_layout.addWidget(lc_table)

        ffi_layout = QVBoxLayout()
        ffi_layout.addSpacing(10)
        ffi_layout.addWidget(ffi_label)
        ffi_layout.addWidget(ffi_table)

        table_layout = QHBoxLayout()
        table_layout.addLayout(timeseries_layout)
        table_layout.addLayout(ffi_layout)

        layout = QVBoxLayout()
        layout.addWidget(star_name_label)
        layout.addWidget(sector_count_label)
        layout.addLayout(table_layout)

        self.data_window.setLayout(layout)
        self.data_window.setGeometry(100, 100, 480, 300)
        self.data_window.show()


## Manual Detrending

    def enable_detrending_region_selection(self):
        """Enable region selection on the graph for detrending."""
        try:
            if hasattr(self, 'detrending_selector'):
                self.detrending_selector.set_active(False)
                del self.detrending_selector

            ax = self.figure_detrending.gca()
            self.original_xlim = ax.get_xlim()
            self.original_ylim = ax.get_ylim()

            # Create a RectangleSelector for selecting detrending region
            self.detrending_selector = RectangleSelector(
                ax, self.on_select_detrending, useblit=True, interactive=True,
                button=[1], minspanx=0.1, minspany=0.1, spancoords='data'
            )

            QMessageBox.information(self, "Detrending Selection",
                                    "Drag to select the region where you want to apply detrending.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")
            raise

    def on_select_detrending_region(self, eclick, erelease):
        """Callback when a region is selected for detrending based on both time and flux values."""
        start_time = eclick.xdata
        end_time = erelease.xdata

        if start_time is not None and end_time is not None:
            self.manual_segments.append((min(start_time, end_time), max(start_time, end_time)))
            QMessageBox.information(self, "Manual Segment Selected",
                                    f"Manual segment added from {start_time:.2f} to {end_time:.2f}.")
        else:
            QMessageBox.warning(self, "Warning", "No valid region selected.")

    def on_select_detrending(self, eclick, erelease):
        """Callback when a region is selected for detrending based on both time and flux values."""
        self.detrending_start_time = eclick.xdata
        self.detrending_end_time = erelease.xdata

        if self.detrending_start_time and self.detrending_end_time:
            self.manual_segments.append((self.detrending_start_time, self.detrending_end_time))
            print(f"Manual segment added: {self.detrending_start_time} to {self.detrending_end_time}")  # Debugging için
            QMessageBox.information(self, "Manual Segment Selected",
                                    f"Manual segment added from {self.detrending_start_time:.2f} to {self.detrending_end_time:.2f}.")
        else:
            QMessageBox.warning(self, "Warning", "No valid region selected.")

    def enable_region_selection(self):
        """Enable region selection on the graph and set mode to manual."""
        try:
            self.segment_mode = "manual"
            if hasattr(self, 'selector'):
                self.selector.set_active(False)
                del self.selector

            ax = self.figure_outliers.gca()
            self.selector = RectangleSelector(ax, self.on_select, useblit=True, interactive=True,
                                              button=[1], minspanx=0.1, minspany=0.1, spancoords='data')
            QMessageBox.information(self, "Region Selection", "Drag to select the region you want to remove.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred: {e}")

    def show_manual_trend(self):
        try:
            degree = int(self.degree_combo.currentText().split()[0][0])

            if not self.manual_segments:
                QMessageBox.warning(self, "Warning", "No manual regions selected for detrending.")
                return

            ax = self.figure_detrending.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color="dodgerblue", markersize=5, capsize=1,
                        elinewidth=1.5)

            for start_time, end_time in self.manual_segments:
                mask = (self.time >= start_time) & (self.time <= end_time)
                time_segment = self.time[mask]
                flux_segment = self.flux[mask]
                flux_err_segment = self.flux_err[mask]

                if len(time_segment) == 0:
                    continue

                pdc_coef = np.polyfit(time_segment, flux_segment, w=1 / flux_err_segment, deg=degree)
                pdc_parabola = np.poly1d(pdc_coef)
                pdc_par_val = pdc_parabola(time_segment)

                ax.errorbar(time_segment, pdc_par_val, linestyle="--", color="red")


            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrending of Light Curve")
            self.canvas_detrending.draw()

            QMessageBox.information(self, "Manual Detrending Applied",
                                    f"Polynomial fit of degree {degree} applied to manual regions and displayed as a dashed line.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while showing manual trend: {e}")

    def remove_manual_trend(self):
        """Remove polynomial fit (detrend) for manually selected segments."""
        try:
            if not hasattr(self, 'pre_detrend_time'):
                self.pre_detrend_time = self.time.copy()
                self.pre_detrend_flux = self.flux.copy()
                self.pre_detrend_flux_err = self.flux_err.copy()

            if not self.manual_segments:
                QMessageBox.warning(self, "Warning", "No regions selected for manual detrending.")
                return

            degree = int(self.degree_combo.currentText().split()[0][0])

            updated_flux = self.flux.copy()
            updated_flux_err = self.flux_err.copy()

            for start_time, end_time in self.manual_segments:
                # Seçilen segment için mask oluştur
                mask = (self.time >= start_time) & (self.time <= end_time)
                time_segment = self.time[mask]
                flux_segment = self.flux[mask]
                flux_err_segment = self.flux_err[mask]

                if len(time_segment) == 0:
                    continue

                pdc_coef = np.polyfit(time_segment, flux_segment, w=1 / flux_err_segment, deg=degree)
                pdc_parabola = np.poly1d(pdc_coef)
                pdc_par_val = pdc_parabola(time_segment)

                detrended_flux = flux_segment - pdc_par_val + 1

                updated_flux[mask] = detrended_flux
                updated_flux_err[mask] = flux_err_segment

            self.flux = updated_flux
            self.flux_err = updated_flux_err

            median_flux = np.median(self.flux)
            self.flux /= median_flux
            self.flux_err /= median_flux


            self.manual_segments.clear()

            ax = self.figure_detrending.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5)

            ax.set_ylabel("Normalized Flux")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrended and Normalized Light Curve - Manual")
            self.canvas_detrending.draw()

            QMessageBox.information(self, "Manual Detrending Applied",
                                    "Polynomial trend removed from manual segments, and data is detrended and normalized.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while removing manual detrending: {e}")

    ## Automate Detrending
    def detect_gaps(self):
        """Detect gaps and set segment mode to automate."""
        try:
            self.segment_mode = "automate"
            self.segment_data()

            if not self.automated_segments:
                QMessageBox.warning(self, "Warning", "No segments detected.")
                return

            ax = self.figure_detrending.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")

            ax.set_ylabel("Normalized Flux")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Light Curve with Segments")

            for i, (start_time, end_time) in enumerate(self.automated_segments):
                ax.axvspan(start_time, end_time, color="gray", alpha=0.3, label="Detected Segment" if i == 0 else "")

            ax.legend(loc="best")
            self.canvas_detrending.draw()

            QMessageBox.information(self, "Gaps Detected",
                                    f"{len(self.automated_segments)} segment(s) detected based on the threshold of {self.threshold_input.value()}.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while detecting gaps: {e}")

    def segment_data(self):
        """Segment the data based on the selected gap threshold."""
        threshold = self.threshold_input.value()
        time_diffs = np.diff(self.time)
        gap_indices = np.where(time_diffs > threshold)[0] + 1

        self.time_segments = np.split(self.time, gap_indices)
        self.flux_segments = np.split(self.flux, gap_indices)
        self.flux_err_segments = np.split(self.flux_err, gap_indices)

        self.automated_segments = [(segment[0], segment[-1]) for segment in self.time_segments if len(segment) > 0]

    def show_automated_trend(self):
        try:
            degree = int(self.degree_combo.currentText().split()[0][0])

            if not self.automated_segments:
                QMessageBox.warning(self, "Warning", "No automated regions selected for detrending.")
                return

            ax = self.figure_detrending.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color="dodgerblue", markersize=5, capsize=1,
                        elinewidth=1.5)

            for start_time, end_time in self.automated_segments:
                mask = (self.time >= start_time) & (self.time <= end_time)
                time_segment = self.time[mask]
                flux_segment = self.flux[mask]
                flux_err_segment = self.flux_err[mask]

                if len(time_segment) == 0:
                    continue

                pdc_coef = np.polyfit(time_segment, flux_segment, w=1 / flux_err_segment, deg=degree)
                pdc_parabola = np.poly1d(pdc_coef)

                pdc_par_val = pdc_parabola(time_segment)


                ax.errorbar(time_segment, pdc_par_val, linestyle="--", color="red")


            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrending of Light Curve")
            self.canvas_detrending.draw()

            QMessageBox.information(self, "Automated Detrending Applied",
                                    f"Polynomial fit of degree {degree} applied to automated regions and displayed as a dashed line.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while showing automated trend: {e}")

    def remove_automated_trend(self):
        """Remove polynomial fit (detrend) for automatically detected segments."""
        try:
            if not hasattr(self, 'pre_detrend_time'):
                self.pre_detrend_time = self.time.copy()
                self.pre_detrend_flux = self.flux.copy()
                self.pre_detrend_flux_err = self.flux_err.copy()

            if not self.automated_segments:
                QMessageBox.warning(self, "Warning", "No segments detected for automated detrending.")
                return

            degree = int(self.degree_combo.currentText().split()[0][0])

            updated_flux = self.flux.copy()
            updated_flux_err = self.flux_err.copy()

            for idx, (start_time, end_time) in enumerate(self.automated_segments):
                mask = (self.time >= start_time) & (self.time <= end_time)
                time_segment = self.time[mask]
                flux_segment = self.flux[mask]
                flux_err_segment = self.flux_err[mask]

                if len(time_segment) == 0:
                    continue

                pdc_coef = np.polyfit(time_segment, flux_segment, w=1 / flux_err_segment, deg=degree)

                pdc_parabola = np.poly1d(pdc_coef)
                pdc_par_val = pdc_parabola(time_segment)
                detrended_flux = flux_segment - pdc_par_val + 1

                updated_flux[mask] = detrended_flux
                updated_flux_err[mask] = flux_err_segment

            self.flux = updated_flux
            self.flux_err = updated_flux_err

            median_flux = np.median(self.flux)
            self.flux /= median_flux
            self.flux_err /= median_flux

            self.automated_segments.clear()

            ax = self.figure_detrending.gca()
            ax.clear()
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5)

            ax.set_ylabel("Normalized Flux")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrended and Normalized Light Curve")
            self.canvas_detrending.draw()

            QMessageBox.information(self, "Automated Detrending Applied",
                                    "Polynomial trend removed from automated segments, and data is detrended and normalized.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while removing automated detrending: {e}")

    def initialize_cleaned_data(self):
        """Initialize cleaned data attributes."""
        self.cleaned_time = None
        self.cleaned_flux = None
        self.cleaned_flux_err = None

    def restore_data(self):
        """Restore data in the Detrending tab to match the current data in the Outlier tab."""
        try:
            # If pre-detrend data exists, restore it
            if hasattr(self, 'pre_detrend_flux'):
                self.time = self.pre_detrend_time.copy()
                self.flux = self.pre_detrend_flux.copy()
                self.flux_err = self.pre_detrend_flux_err.copy()
                # Remove the pre-detrend data attributes
                del self.pre_detrend_time
                del self.pre_detrend_flux
                del self.pre_detrend_flux_err
            else:
                # If no detrending was applied, use the current data from the Outlier tab
                pass  # self.time, self.flux, self.flux_err are already the current data

            # Plot the data in the Detrending tab
            self.plot_data_in_detrending()
            QMessageBox.information(self, "Data Restored",
                                    "Detrending data has been restored to match the current data in the Outlier tab.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while restoring data: {e}")

    def update_outlier_data(self, time, flux, flux_err):
        """Update the outlier data and save the latest cleaned data."""
        self.cleaned_time = time
        self.cleaned_flux = flux
        self.cleaned_flux_err = flux_err


        ax = self.figure_outliers.gca()
        ax.clear()
        ax.errorbar(self.cleaned_time, self.cleaned_flux, yerr=self.cleaned_flux_err, fmt=".", color="dodgerblue",
                    markersize=5, capsize=1, elinewidth=1.5, ls="none")
        ax.set_ylabel("Normalized Flux ($electrons^{-1}$)")
        ax.set_xlabel("Time - 2457000 ($day$)")
        ax.set_title("Normalized Light Curve")
        self.canvas_outliers.draw()

        self.plot_data_in_detrending()

    def restore_original_data(self):
        try:
            if self.cleaned_time is None or self.cleaned_flux is None or self.cleaned_flux_err is None:
                QMessageBox.warning(self, "Warning", "No cleaned data available to restore.")
                return

            self.outlier_plot_layout.removeWidget(self.canvas_outliers)
            self.canvas_outliers.deleteLater()
            self.toolbar_outliers.deleteLater()
            self.outlier_error_checkbox.setParent(None)  # Checkbox'ı da yerinden kaldır

            self.figure_outliers = Figure()
            self.canvas_outliers = FigureCanvas(self.figure_outliers)
            self.toolbar_outliers = NavigationToolbar(self.canvas_outliers, self)

            toolbar_layout_outliers = QHBoxLayout()
            toolbar_layout_outliers.addWidget(self.toolbar_outliers)
            toolbar_layout_outliers.addWidget(self.outlier_error_checkbox)
            toolbar_layout_outliers.setAlignment(Qt.AlignLeft)

            self.outlier_plot_layout.addLayout(toolbar_layout_outliers)
            self.outlier_plot_layout.addWidget(self.canvas_outliers)

            ax = self.figure_outliers.add_subplot(111)
            ax.errorbar(self.original_time, self.original_flux, yerr=self.original_flux_err, fmt=".",
                        color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")
            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Original Light Curve from Query Tab")
            self.canvas_outliers.draw()

            self.time = self.original_time
            self.flux = self.original_flux
            self.flux_err = self.original_flux_err

            self.plot_data_in_detrending()

            QMessageBox.information(self, "Data Restored",
                                    "The cleaned data from the Outlier section has been restored in the detrending plot.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while restoring the cleaned data: {e}")


    def plot_data_in_detrending(self):
        """Plot data in the Detrending tab based on the current state of self.time, self.flux, and self.flux_err."""
        ax = self.figure_detrending.gca()
        ax.clear()  # Clear any previous plots

        if self.detrending_error_checkbox.isChecked() and self.flux_err is not None:
            ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt=".", color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")
        else:
            ax.plot(self.time, self.flux, ".", color="dodgerblue", markersize=5)

        ax.set_ylabel("Flux ($electrons^{-1}$)")
        ax.set_xlabel("Time - 2457000 ($day$)")
        ax.set_title("Detrending of Light Curve")

        self.canvas_detrending.draw()

    def plot_data_in_phase_folding(self, detrended_data=None):
        """Plot data in the Phase Folding tab, including errors if provided."""
        try:
            ax = self.figure_phase.gca()
            ax.clear()

            if detrended_data:
                time_segments, flux_segments, flux_err_segments = detrended_data
                for time_segment, flux_segment, flux_err_segment in zip(time_segments, flux_segments,
                                                                        flux_err_segments):
                    ax.errorbar(time_segment, flux_segment, yerr=flux_err_segment, fmt='.', color="dodgerblue",
                                markersize=5, capsize=1, elinewidth=1.5, label="Detrended and Normalized Data")
            else:
                ax.errorbar(self.time, self.flux, yerr=self.flux_err, fmt='.', color="dodgerblue", markersize=5,
                            capsize=1, elinewidth=1.5, label="Original Data")

            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Phase")
            ax.set_title("Phase Folded Light Curve")
            ax.legend(loc="best")
            self.canvas_phase.draw()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while plotting in Phase Folding tab: {e}")

    def detach_detrending_plot(self):
        try:
            if self.cleaned_time is None or self.cleaned_flux is None or self.cleaned_flux_err is None:
                QMessageBox.warning(self, "Warning", "No data available for plotting.")
                return

            fig, ax = plt.subplots()
            ax.errorbar(self.cleaned_time, self.cleaned_flux, yerr=self.cleaned_flux_err, fmt=".", color="dodgerblue",
                        markersize=5, capsize=1, elinewidth=1.5, ls="none")
            ax.set_ylabel("Flux ($electrons^{-1}$)")
            ax.set_xlabel("Time - 2457000 ($day$)")
            ax.set_title("Detrended Light Curve - Detached View")
            fig.tight_layout()
            fig.show()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while detaching the detrending plot: {e}")

    def save_detrending_data(self):
        """Save only the detrended data in the chosen file format."""
        try:
            dialog = QDialog(self)
            dialog.setWindowTitle("Save Detrended Data Options")

            layout = QVBoxLayout(dialog)
            layout.addWidget(QLabel("Select file format to save:"))

            file_format_combo = QComboBox()
            file_format_combo.addItems(["XLSX", "CSV", "TXT"])  # Default order with XLSX first
            layout.addWidget(file_format_combo)

            add_2457000_checkbox = QCheckBox("Add 2457000 to Time")
            add_2457000_checkbox.setChecked(True)  # Checked by default
            layout.addWidget(add_2457000_checkbox)

            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)

            if dialog.exec_() != QDialog.Accepted:
                return

            file_format = file_format_combo.currentText().lower()  # Convert to lowercase for consistency

            file_name, _ = QFileDialog.getSaveFileName(
                self, f"Save Detrended Data as {file_format.upper()}", "",
                f"{file_format.upper()} Files (*.{file_format})"
            )
            if not file_name:
                return

            if not file_name.endswith(f".{file_format}"):
                file_name += f".{file_format}"

            if self.time is not None:
                time_data = self.time + 2457000 if add_2457000_checkbox.isChecked() else self.time
                detrended_data = pd.DataFrame({
                    'Time': time_data,
                    'Flux': self.flux,
                    'Error': self.flux_err
                })
            else:
                QMessageBox.warning(self, "No Data Available", "Detrended data is not available to save.")
                return

            try:
                if file_format == "xlsx":
                    with pd.ExcelWriter(file_name, engine="openpyxl", mode="w") as writer:
                        detrended_data.to_excel(writer, sheet_name="Detrended Data", index=False, float_format="%.10f")
                elif file_format == "csv":
                    detrended_data.to_csv(file_name, index=False, float_format="%.10f")
                elif file_format == "txt":
                    np.savetxt(
                        file_name,
                        detrended_data.values,
                        fmt="%.10f",
                        header="\t".join(detrended_data.columns),
                        delimiter="\t",
                        comments=""
                    )

                QMessageBox.information(self, "Data Saved",
                                        f"Detrended data has been saved successfully as {file_format.upper()}.")

            except Exception as e:
                QMessageBox.critical(self, "Error",
                                     f"An error occurred while saving detrended data as {file_format.upper()}: {e}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}")

## Phase Tab

    def manual_phase_folding(self):
        try:
            t0 = float(self.t0_input.text().strip())
            period = float(self.period_input.text().strip())

            if not t0 or not period:
                raise ValueError("Both t0 and period must be provided.")

            if self.time.size == 0 or self.flux.size == 0 or self.flux_err.size == 0:
                raise ValueError("Loaded data is empty. Please ensure the TESS data is loaded properly.")

            phase, flux, flux_err = phase_fold_manually(
                time=self.time,
                flux=self.flux,
                flux_err=self.flux_err,
                t0=t0,
                period=period,
                canvas=self.canvas_phase
            )

            self.phase_time = phase
            self.phase_flux = flux
            self.phase_flux_err = flux_err

            QMessageBox.information(self, "Manual Phase Folding", "Manual phase folding completed successfully.")

        except ValueError as ve:
            QMessageBox.warning(self, "Input Error", f"Invalid input for t0 or period: {ve}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while manually phase folding: {e}")

    def automate_phase_folding(self):
        try:
            t0_final, period_estimate, lomb_scargle_period, phase, flux, flux_err = phase_fold_and_normalize(
                time=self.time,
                flux=self.flux,
                flux_err=self.flux_err,
                canvas=self.canvas_phase,  # Phase fold tabındaki canvas'ı gönderiyoruz
                display_fields={
                    "best_t0_display": self.best_t0_display,
                    "best_period_display": self.best_period_display,
                    "period_dropdown": self.period_dropdown
                }
            )

            self.phase_time = phase
            self.phase_flux = flux
            self.phase_flux_err = flux_err

            QMessageBox.information(self, "Automate Phase Folding", "Automated phase folding completed successfully.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred during automatic phase folding: {e}")

    def update_phase_plot(self):
        """Update the phase plot based on the checkbox state."""
        self.plot_data_check(self.figure_phase, self.canvas_phase, self.phase_error_checkbox.isChecked())

    def detach_phase_plot(self):
        """Detach the phase-folded plot and display it in a separate window."""
        try:
            if not hasattr(self, "phase_time") or self.phase_time is None:
                QMessageBox.warning(self, "Warning", "No data available for plotting.")
                return

            fig, ax = plt.subplots()

            ax.errorbar(
                self.phase_time, self.phase_flux, yerr=self.phase_flux_err, fmt='.', color='dodgerblue',
                markersize=5, capsize=1, elinewidth=1.5, ls='none'
            )
            ax.set_ylabel("Normalized Flux")
            ax.set_xlabel("Phase")
            ax.set_title("Detached Light Curve - Phase-Folded Data")

            ax.set_xlim(0.2, 1.2)  # Only show phase values between 0.2 and 1.2

            fig.tight_layout()
            fig.show()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An error occurred while detaching the phase plot: {e}")

    def save_phase_data(self):
        """Save the phase-folded data in the chosen file format."""
        try:
            dialog = QDialog(self)
            dialog.setWindowTitle("Save Phase-Folded Data Options")

            layout = QVBoxLayout(dialog)
            layout.addWidget(QLabel("Select file format to save:"))

            file_format_combo = QComboBox()
            file_format_combo.addItems(["XLSX", "CSV", "TXT"])  # Default order with XLSX first
            layout.addWidget(file_format_combo)

            button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
            button_box.accepted.connect(dialog.accept)
            button_box.rejected.connect(dialog.reject)
            layout.addWidget(button_box)

            if dialog.exec_() != QDialog.Accepted:
                return

            file_format = file_format_combo.currentText().lower()  # Convert to lowercase for consistency

            file_name, _ = QFileDialog.getSaveFileName(
                self, f"Save Phase-Folded Data as {file_format.upper()}", "",
                f"{file_format.upper()} Files (*.{file_format})"
            )
            if not file_name:
                return

            if not file_name.endswith(f".{file_format}"):
                file_name += f".{file_format}"

            if hasattr(self, "phase_time") and self.phase_time is not None:
                phase_data = pd.DataFrame({
                    'Phase': self.phase_time,
                    'Flux': self.phase_flux,
                    'Error': self.phase_flux_err
                })
            else:
                QMessageBox.warning(self, "No Data Available", "Phase-folded data is not available to save.")
                return

            try:
                if file_format == "xlsx":
                    with pd.ExcelWriter(file_name, engine="openpyxl", mode="w") as writer:
                        phase_data.to_excel(writer, sheet_name="Phase-Folded Data", index=False, float_format="%.10f")
                elif file_format == "csv":
                    phase_data.to_csv(file_name, index=False, float_format="%.10f")
                elif file_format == "txt":
                    np.savetxt(
                        file_name,
                        phase_data.values,
                        fmt="%.10f",
                        header="\t".join(phase_data.columns),
                        delimiter="\t",
                        comments=""
                    )

                QMessageBox.information(self, "Data Saved",
                                        f"Phase-folded data has been saved successfully as {file_format.upper()}.")

            except Exception as e:
                QMessageBox.critical(self, "Error",
                                     f"An error occurred while saving phase-folded data as {file_format.upper()}: {e}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}")


def main():
    app = QApplication(sys.argv)
    app.setAttribute(Qt.AA_DisableWindowContextHelpButton, True)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()