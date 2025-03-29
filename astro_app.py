import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from astroquery.simbad import Simbad
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body
from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astroplan import Observer
import astropy.units as u
from datetime import timedelta
from zoneinfo import ZoneInfo
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QDateEdit, QComboBox, \
    QStyledItemDelegate, QTabWidget, QSpacerItem, QSizePolicy, QHBoxLayout, QTableWidget, QTableWidgetItem, QDialog, \
    QFormLayout, QMenu, QAbstractItemView, QShortcut, QFileDialog, QCheckBox, QMessageBox, QHeaderView, QAbstractScrollArea
from PyQt5.QtCore import QDate, Qt
from PyQt5.QtGui import QKeySequence, QIcon
import pandas as pd
import shutil
import os


class CenterAlignDelegate(QStyledItemDelegate):
    def paint(self, painter, option, index):
        option.displayAlignment = Qt.AlignCenter
        super().paint(painter, option, index)


class SaveEphemerisDialog(QDialog):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Save Ephemeris Data')
        layout = QFormLayout()

        self.object_input = QLineEdit(self)
        self.ref_minima_input = QLineEdit(self)
        self.period_input = QLineEdit(self)

        layout.addRow('Object:', self.object_input)
        layout.addRow('Reference Minima Time:', self.ref_minima_input)
        layout.addRow('Period:', self.period_input)

        # Pre-fill inputs with existing data if available
        if self.parent.star_input.text():
            self.object_input.setText(self.parent.star_input.text())
        if self.parent.t0_input.text():
            self.ref_minima_input.setText(self.parent.t0_input.text())
        if self.parent.P_input.text():
            self.period_input.setText(self.parent.P_input.text())

        save_button = QPushButton('Save', self)
        save_button.clicked.connect(self.save_data)
        layout.addWidget(save_button)

        self.setLayout(layout)


    def save_data(self):
        obj = self.object_input.text()
        ref_minima = self.ref_minima_input.text()
        period = self.period_input.text()

        with open('last_ephemeris_data.txt', 'a') as file:
            file.write(f"{obj},{ref_minima},{period}\n")

        # Update the ephemeris data in the parent
        self.parent.ephemeris_data.append([obj, ref_minima, period])

        # Save the updated data to 'last_ephemeris_data.csv'
        df_to_save = pd.DataFrame(self.parent.ephemeris_data, columns=['Object Name', 'Reference Minima Time', 'Period'])
        df_to_save.to_csv('last_ephemeris_data.csv', index=False)

        self.accept()



class AddLocationDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Add New Location')
        layout = QFormLayout()

        self.name_input = QLineEdit(self)
        self.lat_input = QLineEdit(self)
        self.lon_input = QLineEdit(self)
        self.alt_input = QLineEdit(self)
        self.timezone_input = QLineEdit(self)
        self.timezone_input.setPlaceholderText("e.g., Europe/Istanbul")

        layout.addRow('Location Name:', self.name_input)
        layout.addRow('Latitude (deg):', self.lat_input)
        layout.addRow('Longitude (deg):', self.lon_input)
        layout.addRow('Altitude (m):', self.alt_input)
        layout.addRow('Time Zone:', self.timezone_input)

        save_button = QPushButton('Add Location', self)
        save_button.clicked.connect(self.add_location)
        layout.addWidget(save_button)

        self.setLayout(layout)

    def add_location(self):
        name = self.name_input.text()
        lat = self.lat_input.text()
        lon = self.lon_input.text()
        alt = self.alt_input.text()
        timezone = self.timezone_input.text()
        if name and lat and lon and alt and timezone:
            try:
                lat = float(lat)
                lon = float(lon)
                alt = float(alt)
                location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
                self.parent.add_new_location(name, location, timezone)

                # Yeni gözlemevini 'observatories.txt' dosyasına alt satıra ekleyelim
                with open('observatories.txt', 'a') as file:
                    file.write(f"{name}\t{lat}\t{lon}\t{alt}\t{timezone}\n")

                self.accept()
            except ValueError:
                print("Please enter valid numerical values for latitude, longitude, and altitude.")
        else:
            print("All fields must be filled.")



class LoadEphemerisDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.loaded_columns = None  # To store selected columns
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Saved Ephemeris Data')
        layout = QVBoxLayout()
        # self.setGeometry(800, 400, 500, 800)  # Bu satırı kaldırıyoruz
        self.resize(800, 600)  # Başlangıç boyutu
        self.setMinimumSize(600, 400)  # Minimum boyut

        # Search Box
        search_layout = QHBoxLayout()
        self.search_input = QLineEdit(self)
        self.search_input.setPlaceholderText("Search by star name...")
        self.search_input.textChanged.connect(self.search_table)
        search_layout.addWidget(self.search_input)

        layout.addLayout(search_layout)

        self.table = QTableWidget(self)

        # Tablo boyutlandırma ayarları
        self.table.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        # Sütunları otomatik genişleyecek şekilde ayarla
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        # Initialize with three columns
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(['Object Name', 'Reference Minima Time', 'Period'])

        self.load_data()

        self.table.setContextMenuPolicy(Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self.open_context_menu)
        self.table.cellDoubleClicked.connect(self.load_selected_data)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        layout.addWidget(self.table)

        save_button = QPushButton('Save Changes', self)
        save_button.clicked.connect(self.save_changes)
        layout.addWidget(save_button)

        # Load from file button
        load_file_button = QPushButton('Load from File', self)
        load_file_button.clicked.connect(self.load_from_file)
        layout.addWidget(load_file_button)

        self.table.resizeColumnsToContents()
        self.setLayout(layout)

        # Ctrl+F shortcut
        self.shortcut = QShortcut(QKeySequence("Ctrl+F"), self)
        self.shortcut.activated.connect(self.focus_search)


    def focus_search(self):
        self.search_input.setFocus()

    def load_data(self):
        """Clears the table and updates it with new data."""
        self.table.setRowCount(0)  # Clear the table
        if self.parent.ephemeris_data:
            # Set column headers
            column_headers = ['Object Name', 'Reference Minima Time', 'Period']
            self.table.setColumnCount(3)
            self.table.setHorizontalHeaderLabels(column_headers)

            # Load ephemeris data into the table
            for row_data in self.parent.ephemeris_data:
                row_position = self.table.rowCount()
                self.table.insertRow(row_position)
                for column, data in enumerate(row_data):
                    item = QTableWidgetItem(str(data))
                    item.setTextAlignment(Qt.AlignCenter)  # Metni ortala
                    self.table.setItem(row_position, column, item)
        else:
            print("No ephemeris data loaded.")

    def search_table(self, text):
        """Tabloyu arama kutusuna girilen metne göre filtreler"""
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)  # Object sütunu
            if text.lower() in item.text().lower():
                self.table.setRowHidden(row, False)  # Eğer aranan metin varsa satırı göster
            else:
                self.table.setRowHidden(row, True)  # Yoksa satırı gizle

    def load_selected_data(self, row, column):
        # Seçilen satırdaki ilgili hücrelerden verileri alıyoruz
        obj = self.table.item(row, 0).text()  # Object sütunu
        ref_minima = self.table.item(row, 1).text()  # Reference Minima Time sütunu
        period = self.table.item(row, 2).text()  # Period sütunu

        # Ana ekrandaki alanlara bu değerleri yazıyoruz
        self.parent.star_input.setText(obj)
        self.parent.t0_input.setText(ref_minima)
        self.parent.P_input.setText(period)

        # Pencereyi kapatıyoruz
        self.accept()

    def open_context_menu(self, position):
        context_menu = QMenu(self)
        edit_action = context_menu.addAction("Edit Cell")
        delete_action = context_menu.addAction("Delete Cell")
        action = context_menu.exec_(self.table.viewport().mapToGlobal(position))

        if action == edit_action:
            item = self.table.itemAt(position)
            if item:
                self.table.editItem(item)
        elif action == delete_action:
            row = self.table.indexAt(position).row()
            if row != -1:
                self.table.removeRow(row)

    def save_changes(self):
        # Tablo verilerini alalım
        updated_data = []
        for row in range(self.table.rowCount()):
            row_data = []
            for column in range(self.table.columnCount()):
                item = self.table.item(row, column)
                if item is not None:
                    row_data.append(item.text())
                else:
                    row_data.append('')
            updated_data.append(row_data)

        # self.parent.ephemeris_data'yı güncelleyelim
        self.parent.ephemeris_data = updated_data

        # ephemeris_data.txt dosyasına kaydedelim
        try:
            with open('ephemeris_data.txt', 'w') as file:
                for row_data in updated_data:
                    file.write(','.join(row_data) + '\n')
        except Exception as e:
            QMessageBox.warning(self, "File Save Error", f"An error occurred while saving the file:\n{str(e)}")
            return

        # last_ephemeris_data.csv dosyasına kaydedelim
        try:
            df_to_save = pd.DataFrame(updated_data, columns=['Object Name', 'Reference Minima Time', 'Period'])
            df_to_save.to_csv('last_ephemeris_data.csv', index=False)
        except Exception as e:
            QMessageBox.warning(self, "File Save Error", f"An error occurred while saving the CSV file:\n{str(e)}")
            return

        # Kullanıcıya başarı mesajı gösterelim
        QMessageBox.information(self, "Success", "Changes have been saved successfully.")
    def preview_file(self, file_name, file_extension):
        preview_dialog = QDialog(self)
        preview_dialog.setWindowTitle('Preview and Configure File')
        preview_layout = QVBoxLayout()
        table = QTableWidget(preview_dialog)
        self.preview_table = table  # Tabloyu diğer fonksiyonlarda kullanmak için saklıyoruz

        # Kullanıcıya sırayla sütun seçmesini hatırlatan bir uyarı ekle
        info_label = QLabel('Please select columns in the order: Object Name, Reference Minima Time, Period',
                            preview_dialog)
        info_label.setStyleSheet("color: red; font-weight: bold;")
        preview_layout.addWidget(info_label)

        # İlk satırı atlayıp atlamama seçeneği
        self.first_row_checkbox = QCheckBox('Ignore first row (header)', preview_dialog)
        self.first_row_checkbox.setChecked(False)  # Varsayılan olarak işaretli değil
        preview_layout.addWidget(self.first_row_checkbox)

        # Reference Minima Time'a 2400000 eklemek için bir checkbox
        self.add_2400000_checkbox = QCheckBox('Add 2400000 to Reference Minima Time', preview_dialog)
        preview_layout.addWidget(self.add_2400000_checkbox)

        # Tabloya Ctrl ile birden fazla sütun seçebilme
        table.setSelectionBehavior(QAbstractItemView.SelectColumns)
        table.setSelectionMode(QAbstractItemView.MultiSelection)
        preview_layout.addWidget(table)

        # Tabloyu güncelleyen fonksiyonu tanımlayalım
        def update_preview_table():
            try:
                header_option = 0 if not self.first_row_checkbox.isChecked() else None
                if file_extension in ['.txt', '.csv']:
                    df = pd.read_csv(file_name, header=header_option)
                elif file_extension == '.xlsx':
                    df = pd.read_excel(file_name, header=header_option, engine='openpyxl')
                else:
                    print(f"Unsupported file format: {file_extension}")
                    return

                data = df.values.tolist()
                headers = df.columns.tolist()
                table.setRowCount(len(data))
                table.setColumnCount(len(headers))
                table.setHorizontalHeaderLabels(headers)

                for row_idx, row_data in enumerate(data):
                    for col_idx, cell_data in enumerate(row_data):
                        table.setItem(row_idx, col_idx, QTableWidgetItem(str(cell_data)))

                # Pencereyi dinamik boyutlandırmak için
                table.resizeColumnsToContents()
                table.resizeRowsToContents()

            except Exception as e:
                QMessageBox.warning(self, "File Error", f"An error occurred while reading the file:\n{str(e)}")
                return

        # Checkbox'ın durum değişikliğinde tabloyu güncelle
        self.first_row_checkbox.stateChanged.connect(update_preview_table)

        # İlk kez tabloyu yükle
        update_preview_table()

        load_button = QPushButton('Load Data', preview_dialog)
        load_button.clicked.connect(
            lambda: self.load_data_with_settings(file_name, file_extension, table, preview_dialog))
        preview_layout.addWidget(load_button)

        preview_dialog.setLayout(preview_layout)
        preview_dialog.exec_()

    def load_data_with_settings(self, file_name, file_extension, table, preview_dialog):
        selected_columns = table.selectionModel().selectedColumns()
        selected_columns = sorted([index.column() for index in selected_columns])

        if len(selected_columns) != 3:
            QMessageBox.warning(self, "Selection Error",
                                "Please select exactly 3 columns: Object Name, Reference Minima Time, and Period.")
            return

        # Seçilen sütunları alalım
        object_column = selected_columns[0]
        ref_minima_column = selected_columns[1]
        period_column = selected_columns[2]

        # Seçilen sütunları kaydedelim
        self.loaded_columns = (object_column, ref_minima_column, period_column)

        # Dosyayı pandas ile okuyalım
        try:
            header_option = 0 if not self.first_row_checkbox.isChecked() else None
            if file_extension in ['.txt', '.csv']:
                df = pd.read_csv(file_name, header=header_option)
            elif file_extension == '.xlsx':
                df = pd.read_excel(file_name, header=header_option, engine='openpyxl')
            else:
                QMessageBox.warning(self, "File Error", "Unsupported file format.")
                return
        except Exception as e:
            QMessageBox.warning(self, "File Error", f"An error occurred while reading the file:\n{str(e)}")
            return

        # Sütun isimlerini alalım
        if header_option == 0:
            column_headers = [df.columns[object_column],
                              df.columns[ref_minima_column],
                              df.columns[period_column]]
        else:
            column_headers = ['Object Name', 'Reference Minima Time', 'Period']

        # Tablo başlıklarını güncelle
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(column_headers)

        invalid_rows = []

        # Ephemeris verilerini sıfırlayalım
        self.parent.ephemeris_data = []

        self.table.setRowCount(0)

        for row_idx in range(len(df)):
            obj_name = ''  # obj_name değişkenini burada tanımlıyoruz
            try:
                obj_name = str(df.iloc[row_idx, object_column])

                ref_minima_value = df.iloc[row_idx, ref_minima_column]
                if ref_minima_value == '-' or pd.isnull(ref_minima_value):
                    ref_minima_value = None
                else:
                    ref_minima_value = float(ref_minima_value)
                    if self.add_2400000_checkbox.isChecked():
                        ref_minima_value += 2400000  # Değeri burada güncelliyoruz

                period_value = df.iloc[row_idx, period_column]
                if period_value == '-' or pd.isnull(period_value):
                    period_value = None
                else:
                    period_value = float(period_value)

                # Verileri bellekte saklayalım
                self.parent.ephemeris_data.append([
                    obj_name,
                    str(ref_minima_value) if ref_minima_value is not None else '-',
                    str(period_value) if period_value is not None else '-'
                ])

                # Verileri tabloya ekleyelim
                row_position = self.table.rowCount()
                self.table.insertRow(row_position)

                # Object Name
                item_obj = QTableWidgetItem(obj_name)
                item_obj.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(row_position, 0, item_obj)

                # Reference Minima Time
                item_ref = QTableWidgetItem(str(ref_minima_value) if ref_minima_value is not None else '-')
                item_ref.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(row_position, 1, item_ref)

                # Period
                item_period = QTableWidgetItem(str(period_value) if period_value is not None else '-')
                item_period.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(row_position, 2, item_period)
            except ValueError:
                invalid_rows.append((row_idx + 1, obj_name))

        # Geçersiz satırlar varsa kullanıcıya uyarı göster
        if invalid_rows:
            invalid_rows_text = "\n".join([f"Row {idx}: {obj}" for idx, obj in invalid_rows])
            QMessageBox.warning(self, "Invalid Rows Removed",
                                f"Some rows could not be loaded due to non-numeric values in the selected columns:\n\n{invalid_rows_text}")

        # Ephemeris verilerini dosyaya kaydedelim
        try:
            df_to_save = pd.DataFrame(self.parent.ephemeris_data, columns=column_headers)
            df_to_save.to_csv('last_ephemeris_data.csv', index=False)
        except Exception as e:
            print(f"Failed to save ephemeris data to file: {e}")

        preview_dialog.accept()

        # Tabloyu yeniden boyutlandır
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()

    def load_from_file(self):
        """Kullanıcının bir dosya seçmesini sağlar ve dosyayı önizleme için açar."""
        file_name, _ = QFileDialog.getOpenFileName(self, 'Open File', '',
                                                   'All Files (*.txt *.csv *.xlsx);;Text Files (*.txt);;CSV Files (*.csv);;Excel Files (*.xlsx)')
        if file_name:
            file_extension = '.' + file_name.split('.')[-1].lower()

            # Dosya türünü kontrol et ve ön izleme ekranını aç
            if file_extension in ['.txt', '.csv', '.xlsx']:
                self.preview_file(file_name, file_extension)  # Dosyayı önizleme için aç
            else:
                print(f"Unsupported file format: {file_extension}")


class AstroApp(QWidget):
    def __init__(self):
        super().__init__()
        self.constellation_abbr = [
            "And", "Ant", "Aps", "Aql", "Aqr", "Ara", "Ari", "Aur", "Boo", "Cae", "Cam", "Cap", "Car", "Cas", "Cen",
            "Cep", "Cet", "Cha", "Cir", "CMa", "CMi", "Cnc", "Col", "Com", "CrA", "CrB", "Crt", "Cru", "Crv", "CVn",
            "Cyg", "Del", "Dor", "Dra", "Equ", "Eri", "For", "Gem", "Gru", "Her", "Hor", "Hya", "Hyi", "Ind", "Lac",
            "Leo", "Lep", "Lib", "LMi", "Lup", "Lyn", "Lyr", "Men", "Mic", "Mon", "Mus", "Nor", "Oct", "Oph", "Ori",
            "Pav", "Peg", "Per", "Phe", "Pic", "PsA", "Psc", "Pup", "Pyx", "Ret", "Scl", "Sco", "Sct", "Ser", "Sex",
            "Sge", "Sgr", "Tau", "Tel", "TrA", "Tri", "Tuc", "UMa", "UMi", "Vel", "Vir", "Vol", "Vul"
        ]
        self.custom_locations = {}
        self.basic_plot_created = False
        self.fig = None
        self.ax = None
        self.ephemeris_data = [] # Yüklü dosyayı saklayacak değişken

        # Uygulama başlarken last_ephemeris_data.csv dosyasından verileri yükleyelim
        if os.path.exists('last_ephemeris_data.csv'):
            try:
                df = pd.read_csv('last_ephemeris_data.csv')
                self.ephemeris_data = df.values.tolist()
            except Exception as e:
                print(f"Failed to load ephemeris data from file: {e}")


        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        self.setGeometry(800, 400, 250, 250)
        self.setWindowIcon(QIcon("images/icon.ico"))

        tabs = QTabWidget()
        layout.addWidget(tabs)

        main_tab = QWidget()
        tabs.addTab(main_tab, "Altitude Calculation")

        phase_tab = QWidget()
        tabs.addTab(phase_tab, "Phase Calculation")

        main_layout = QVBoxLayout()
        main_tab.setLayout(main_layout)

        phase_layout = QVBoxLayout()
        phase_tab.setLayout(phase_layout)

        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.star_label = QLabel('Object', self)
        self.star_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.star_label)
        self.star_label.setAlignment(Qt.AlignCenter)

        self.star_input = QLineEdit(self)
        self.star_input.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.star_input)

        self.date_label = QLabel('Date', self)
        self.date_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.date_label)
        self.date_label.setAlignment(Qt.AlignCenter)

        self.date_input = QDateEdit(self)
        self.date_input.setAlignment(Qt.AlignCenter)
        self.date_input.setCalendarPopup(True)
        self.date_input.setDate(QDate.currentDate())
        main_layout.addWidget(self.date_input)

        self.location_label = QLabel('Location', self)
        self.location_label.setStyleSheet("color: black; font-weight: bold;")
        main_layout.addWidget(self.location_label)
        self.location_label.setAlignment(Qt.AlignCenter)

        self.location_combo = QComboBox(self)
        self.location_combo.setEditable(True)
        self.load_locations()
        self.location_combo.lineEdit().setAlignment(Qt.AlignCenter)
        self.location_combo.setItemDelegate(CenterAlignDelegate(self.location_combo))
        self.location_combo.lineEdit().setReadOnly(True)
        main_layout.addWidget(self.location_combo)

        add_location_btn = QPushButton('Add New Location', self)
        add_location_btn.clicked.connect(self.open_add_location_dialog)
        main_layout.addWidget(add_location_btn)

        main_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.plot_basic_btn = QPushButton('Plot Graph', self)
        self.plot_basic_btn.clicked.connect(self.plot_basic_star)
        main_layout.addWidget(self.plot_basic_btn)

        phase_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.t0_label = QLabel('Reference Minima Time', self)
        self.t0_label.setStyleSheet("color: black; font-weight: bold;")
        phase_layout.addWidget(self.t0_label)
        self.t0_label.setAlignment(Qt.AlignCenter)

        self.t0_input = QLineEdit(self)
        self.t0_input.setAlignment(Qt.AlignCenter)
        phase_layout.addWidget(self.t0_input)

        self.P_label = QLabel('Period', self)
        self.P_label.setStyleSheet("color: black; font-weight: bold;")
        phase_layout.addWidget(self.P_label)
        self.P_label.setAlignment(Qt.AlignCenter)

        self.P_input = QLineEdit(self)
        self.P_input.setAlignment(Qt.AlignCenter)
        phase_layout.addWidget(self.P_input)

        button_layout = QHBoxLayout()

        self.save_ephemeris_btn = QPushButton('Save Ephemeris', self)
        self.save_ephemeris_btn.clicked.connect(self.open_save_ephemeris_dialog)
        button_layout.addWidget(self.save_ephemeris_btn)

        self.load_ephemeris_btn = QPushButton('Load Ephemeris', self)
        self.load_ephemeris_btn.clicked.connect(self.open_load_ephemeris_dialog)
        button_layout.addWidget(self.load_ephemeris_btn)

        phase_layout.addLayout(button_layout)

        phase_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.submit_btn = QPushButton('Show Phase Info', self)
        self.submit_btn.clicked.connect(self.plot_star_with_phases)
        phase_layout.addWidget(self.submit_btn)

        self.setLayout(layout)
        self.setWindowTitle('AstPlan')
        self.show()

    def load_locations(self):
        try:
            with open('observatories.txt', 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        name = parts[0]
                        self.location_combo.addItem(name)
        except FileNotFoundError:
            print(f"Observatories file not found.")

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.submit_btn.click()

    def format_star_name(self, star_name):
        star_name = star_name.strip()
        parts = star_name.split()
        if len(parts) > 1 and parts[1].upper() in [abbr.upper() for abbr in self.constellation_abbr]:
            return f"{parts[0].upper()} {self.constellation_abbr[[abbr.upper() for abbr in self.constellation_abbr].index(parts[1].upper())]}"
        return star_name.upper()

    def plot_basic_star(self):
        star_name = self.star_input.text()
        formatted_star_name = self.format_star_name(star_name)
        date = self.date_input.date().toString("yyyy-MM-dd")

        location_name = self.location_combo.currentText()
        location_data = self.get_selected_location(location_name)
        if location_data is None:
            print("Location not found. Please select a valid location.")
            return
        location, timezone_str = location_data

        if location is None or timezone_str is None:
            print("Location or time zone not found. Please select a valid location.")
            return

        try:
            timezone = ZoneInfo(timezone_str)
        except Exception as e:
            print(f"Invalid time zone: {timezone_str}")
            return

        observer = Observer(location=location, timezone=timezone)

        simbad = Simbad()
        simbad.add_votable_fields('ra', 'dec')
        result_table = simbad.query_object(star_name)

        if result_table is None:
            print("Star not found. Please make sure the name is correct.")
            return

        ra = result_table['ra'][0]
        dec = result_table['dec'][0]
        obj = SkyCoord(f"{ra} {dec}", unit=(u.deg, u.deg))

        start_time = observer.sun_set_time(Time(date), which='next')
        end_time = observer.sun_rise_time(Time(date) + 1 * u.day, which='next')

        observing_time = start_time + (end_time - start_time) * np.linspace(0, 1, 1000)
        altaz = AltAz(obstime=observing_time, location=location)
        objaltaz = obj.transform_to(altaz)
        sunaltaz = get_sun(observing_time).transform_to(altaz)

        sun_below_horizon = sunaltaz.alt < 0
        observing_time_filtered = observing_time[sun_below_horizon]
        objaltaz_filtered = objaltaz[sun_below_horizon]

        moonaltaz = get_body("moon", observing_time_filtered).transform_to(altaz[sun_below_horizon])

        obj_altitude = objaltaz_filtered.alt.to_value(u.deg)
        obj_azimuth = objaltaz_filtered.az.to_value(u.deg)
        sun_altitude = sunaltaz.alt[sun_below_horizon].to_value(u.deg)
        moon_altitude = moonaltaz.alt.to_value(u.deg)

        # Step 1: Convert to local time with timezone info
        observing_time_local = observing_time_filtered.to_datetime(timezone)

        # Step 2: Remove timezone info to make naive datetime objects
        observing_time_local_naive = np.array([dt.replace(tzinfo=None) for dt in observing_time_local])

        self.fig, self.ax = plt.subplots(figsize=(11, 7))
        sc = self.ax.scatter(observing_time_local_naive, obj_altitude, c=obj_azimuth, cmap='viridis', lw=0, s=8,
                             label="Object")
        self.ax.plot(observing_time_local_naive, moon_altitude, color='lightgray', label='Moon')
        self.ax.fill_between(observing_time_local_naive, 0, 90, sun_altitude < -0, color="k", alpha=0.05, zorder=0)
        self.ax.fill_between(observing_time_local_naive, 0, 90, sun_altitude < -6, color="k", alpha=0.15, zorder=0)
        self.ax.fill_between(observing_time_local_naive, 0, 90, sun_altitude < -12, color="k", alpha=0.25, zorder=0)
        self.ax.fill_between(observing_time_local_naive, 0, 90, sun_altitude < -18, color="k", alpha=0.35, zorder=0)
        self.ax.axhline(20, color="red", ls="-", alpha=0.8)

        self.ax.set_xlim(observing_time_local_naive[0], observing_time_local_naive[-1])
        self.ax.set_xlabel(f'Time ({timezone_str})')
        self.ax.set_ylabel('Altitude (degree)')
        self.ax.set_ylim(0, 90)

        self.ax.set_title(f"Altitude-Time Graph of {formatted_star_name} on {date}", fontstyle='italic', fontsize=13)

        self.ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))

        plt.colorbar(sc, label='Azimuth (degree)', ax=self.ax)
        self.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.075), fancybox=True, ncol=5)
        self.ax.grid()

        self.basic_plot_created = True
        self.fig.show()

    def plot_star_with_phases(self):
        # Mevcut figürü kapatıp bayrağı sıfırla
        if self.fig is not None and plt.fignum_exists(self.fig.number):
            plt.close(self.fig)
            self.basic_plot_created = False

        # Yeni temel grafiği oluştur (bu, self.fig ve self.ax değerlerini ayarlamalı)
        self.plot_basic_star()

        # Girdi verilerini al
        star_name = self.star_input.text().strip()
        formatted_star_name = self.format_star_name(star_name)
        date = self.date_input.date().toString("yyyy-MM-dd")
        t0 = self.t0_input.text().strip()
        P = self.P_input.text().strip()

        # Konum bilgilerini al
        location_name = self.location_combo.currentText()
        location_data = self.get_selected_location(location_name)
        if location_data is None:
            print("Location not found. Please select a valid location.")
            return
        location, timezone_str = location_data
        if location is None or timezone_str is None:
            print("Location or time zone not found. Please select a valid location.")
            return

        try:
            timezone = ZoneInfo(timezone_str)
        except Exception as e:
            print(f"Invalid time zone: {timezone_str}")
            return

        observer = Observer(location=location, timezone=timezone)

        # Simbad sorgusunu optimize etmek için caching (aynı yıldız tekrar sorgulanırsa yeniden sorgulanmaz)
        simbad = Simbad()
        simbad.add_votable_fields('ra', 'dec')
        if hasattr(self, 'last_star_name') and self.last_star_name == star_name:
            result_table = self.last_simbad_result
        else:
            result_table = simbad.query_object(star_name)
            self.last_star_name = star_name
            self.last_simbad_result = result_table

        if result_table is None:
            print("Star not found. Please make sure the name is correct.")
            return

        ra = result_table['ra'][0]
        dec = result_table['dec'][0]
        obj = SkyCoord(f"{ra} {dec}", unit=(u.deg, u.deg))

        # Gözlem zaman aralığı: gün batımından gün doğumuna kadar
        start_time = observer.sun_set_time(Time(date), which='next')
        end_time = observer.sun_rise_time(Time(date) + 1 * u.day, which='next')
        # 1000 nokta üzerinden hesaplama
        observing_time = start_time + (end_time - start_time) * np.linspace(0, 1, 1000)

        altaz = AltAz(obstime=observing_time, location=location)
        objaltaz = obj.transform_to(altaz)
        sunaltaz = get_sun(observing_time).transform_to(altaz)

        # Güneşin ufkun altında olduğu zamanları filtrele
        sun_below_horizon = sunaltaz.alt < 0
        observing_time_filtered = observing_time[sun_below_horizon]
        objaltaz_filtered = objaltaz[sun_below_horizon]

        # Gözlem zamanlarını yerel zamana çevir ve tz bilgisini kaldır
        observing_time_local = observing_time_filtered.to_datetime(timezone)
        observing_time_local_naive = [dt.replace(tzinfo=None) for dt in observing_time_local]

        obj_altitude = objaltaz_filtered.alt.to_value(u.deg)
        obj_azimuth = objaltaz_filtered.az.to_value(u.deg)
        sun_altitude = sunaltaz.alt[sun_below_horizon].to_value(u.deg)

        # Faz hesaplaması (t0 ve P girilmişse)
        phases = None
        if t0 and P:
            try:
                t0_float = float(t0)
                P_float = float(P)
            except ValueError:
                print("Reference Minima Time ve Period sayısal değerler olmalı.")
                return
            P_time = P_float * u.day
            phases = (observing_time_filtered.jd - t0_float) / P_time.value
            phases = phases - np.floor(phases)

        # Faz çizgilerini ekle (while döngüsü yerine vektörleştirilmiş hesaplama)
        if phases is not None:
            first_phase_min_0 = None
            first_phase_min_0_5 = None
            for i, phase in enumerate(phases):
                if first_phase_min_0 is None and np.isclose(phase, 0.0, atol=0.01):
                    first_phase_min_0 = observing_time_filtered[i]
                if first_phase_min_0_5 is None and np.isclose(phase, 0.5, atol=0.01):
                    first_phase_min_0_5 = observing_time_filtered[i]
                if first_phase_min_0 and first_phase_min_0_5:
                    break

            if first_phase_min_0:
                n_steps_0 = int(np.floor((observing_time_filtered[-1].jd - first_phase_min_0.jd) / P_time.value))
                min_0_times = first_phase_min_0 + TimeDelta(np.arange(n_steps_0 + 1) * P_time.value, format='jd')
                min_0_times_local_naive = [t.to_datetime(timezone).replace(tzinfo=None) for t in min_0_times]
                self.ax.axvline(min_0_times_local_naive[0], color='blue', linestyle='--', label="Primary")
                self.ax.text(min_0_times_local_naive[0] - timedelta(minutes=8), 82,
                             min_0_times_local_naive[0].strftime('%H:%M'),
                             color='blue', fontsize=10, ha='right', rotation=90,
                             bbox=dict(facecolor='yellow', alpha=0.8))
                for time_local in min_0_times_local_naive[1:]:
                    self.ax.axvline(time_local, color='blue', linestyle='--')
                    self.ax.text(time_local - timedelta(minutes=8), 82,
                                 time_local.strftime('%H:%M'),
                                 color='blue', fontsize=10, ha='right', rotation=90,
                                 bbox=dict(facecolor='yellow', alpha=0.8))

            if first_phase_min_0_5:
                n_steps_0_5 = int(np.floor((observing_time_filtered[-1].jd - first_phase_min_0_5.jd) / P_time.value))
                min_0_5_times = first_phase_min_0_5 + TimeDelta(np.arange(n_steps_0_5 + 1) * P_time.value, format='jd')
                min_0_5_times_local_naive = [t.to_datetime(timezone).replace(tzinfo=None) for t in min_0_5_times]
                self.ax.axvline(min_0_5_times_local_naive[0], color='purple', linestyle='--', label="Secondary")
                self.ax.text(min_0_5_times_local_naive[0] - timedelta(minutes=8), 82,
                             min_0_5_times_local_naive[0].strftime('%H:%M'),
                             color='purple', fontsize=10, ha='right', rotation=90,
                             bbox=dict(facecolor='yellow', alpha=0.8))
                for time_local in min_0_5_times_local_naive[1:]:
                    self.ax.axvline(time_local, color='purple', linestyle='--')
                    self.ax.text(time_local - timedelta(minutes=8), 82,
                                 time_local.strftime('%H:%M'),
                                 color='purple', fontsize=10, ha='right', rotation=90,
                                 bbox=dict(facecolor='yellow', alpha=0.8))

        # Grafik başlığını ve legend'i güncelle
        self.ax.set_title(f"Phase Information of {formatted_star_name} on {date}", fontstyle='italic', fontsize=13)
        self.ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.075), fancybox=True, ncol=5)

        # Figürü güncelle veya göster
        if plt.fignum_exists(self.fig.number):
            self.fig.canvas.draw()
        else:
            self.fig.show()

    def open_save_ephemeris_dialog(self):
        dialog = SaveEphemerisDialog(self)
        dialog.exec_()

    def open_load_ephemeris_dialog(self):
        dialog = LoadEphemerisDialog(self)
        dialog.exec_()

    def open_add_location_dialog(self):
        dialog = AddLocationDialog(self)
        dialog.exec_()

    def add_new_location(self, name, location, timezone):
        self.custom_locations[name] = (location, timezone)
        self.location_combo.addItem(name)

    def get_selected_location(self, name):
        try:
            with open('observatories.txt', 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5 and parts[0] == name:
                        lat, lon, alt, timezone = float(parts[1]), float(parts[2]), float(parts[3]), parts[4]
                        location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=alt * u.m)
                        return location, timezone
        except FileNotFoundError:
            print("Observatories file not found.")
        return self.custom_locations.get(name, (None, None))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = AstroApp()
    sys.exit(app.exec_())