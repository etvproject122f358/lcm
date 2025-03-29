from astroquery.mast import Observations, Tesscut
from PyQt5.QtWidgets import QMessageBox
from astropy.io import fits
import re
from astropy.timeseries import LombScargle
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.mixture import GaussianMixture
from scipy.interpolate import UnivariateSpline

def get_tess_sector_data(star_name):
    query_result = Observations.query_criteria(obs_collection="TESS", objectname=star_name, dataproduct_type=["timeseries", "image"], radius=0)

    lc_data = []
    ffi_data = []

    for i in query_result:
        if i["dataproduct_type"] == "timeseries":
            if i['dataURL'][-9:-5] == "s_lc" or i['dataURL'][-12:-5] == "fast-lc":
                lc_data.append([i["sequence_number"], i["t_exptime"]])
        elif i["dataproduct_type"] == "image":
            if i["target_name"] == "TESS FFI":
                ffi_data.append([i["sequence_number"], "FFI", round(i["t_exptime"], 1)])

    lc_data.sort(key=lambda x: (x[0], x[1]))
    ffi_data.sort(key=lambda x: (x[0], x[2]))

    return lc_data, ffi_data

def search_tess_sector(star_name):
    star_name = star_name.strip()
    if not star_name:
        return "Please enter a star name.", []

    try:
        query_result = Tesscut.get_sectors(objectname=star_name)
    except Exception as e:
        return f"Error querying TESS sectors: {e}", []

    if len(query_result) == 0:
        return f"No TESS observations found for star: {star_name}", []

    sector_numbers = []
    for sector in query_result['sectorName']:
        parts = sector.split('-')
        if len(parts) > 1:
            try:
                sector_number = int(re.findall(r'\d+', parts[1])[0])
                sector_numbers.append(sector_number)
            except (ValueError, IndexError):
                continue

    sector_numbers = sorted(sector_numbers)
    sectors_text = "<br>".join(map(str, sector_numbers))

    message = (
        f"<div style='text-align: center;'>"
        f"TESS observations for {star_name} found in sectors:<br><br>{sectors_text}</div>"
    )

    return message, sector_numbers

def tess_analysis(star_name, sector, exptime, flux_type, parent=None, apply_filter=True):
    try:
        # Query TESS data
        query = Observations.query_criteria(obs_collection="TESS",
                                            objectname=star_name,
                                            dataproduct_type="timeseries",
                                            sequence_number=sector,
                                            radius=0)
        if len(query) == 0:
            QMessageBox.warning(parent, "Warning", f"No data found for {star_name} in sector {sector}.")
            return None, None, None, None

        product_type = "LC" if exptime == 120 else "FAST-LC"
        products = Observations.get_product_list(query)
        filtered = Observations.filter_products(products, productSubGroupDescription=product_type)
        if len(filtered) == 0:
            QMessageBox.warning(parent, "Warning", f"No {product_type} data found for {star_name} in sector {sector}.")
            return None, None, None, None

        download_table = Observations.download_products(filtered, mrp_only=True)
        path = download_table["Local Path"][0]
        hdul = fits.open(path)

        time = hdul[1].data["TIME"]
        flux = hdul[1].data[flux_type]
        flux_err = hdul[1].data[f"{flux_type}_ERR"]
        quality = hdul[1].data["QUALITY"]

        if apply_filter:
            mask = (quality == 0) & ~np.isnan(time) & ~np.isnan(flux) & ~np.isnan(flux_err)

        else:
            mask = ~np.isnan(time) & ~np.isnan(flux) & ~np.isnan(flux_err)


        time, flux, flux_err, quality = time[mask], flux[mask], flux_err[mask], quality[mask]

        parent.time, parent.flux, parent.flux_err = time, flux, flux_err
        parent.canvas.figure.clear()
        ax_query = parent.canvas.figure.add_subplot(111)
        ax_query.errorbar(time, flux, yerr=flux_err, fmt=".", color="dodgerblue", markersize=5, capsize=1,
                          elinewidth=1.5, ls="none")
        ax_query.set_ylabel("Flux ($electrons^{-1}$)")
        ax_query.set_xlabel("Time - 2457000 ($day$)")
        ax_query.set_title(f"{flux_type.replace('_', ' ')} Light Curve", fontstyle="italic")
        parent.canvas.draw()

        parent.canvas_outliers.figure.clear()
        ax_outliers = parent.canvas_outliers.figure.add_subplot(111)
        ax_outliers.errorbar(time, flux, yerr=flux_err, fmt=".", color="dodgerblue",
                             markersize=5, capsize=1, elinewidth=1.5, ls="none")
        ax_outliers.set_ylabel("Flux ($electrons^{-1}$)")
        ax_outliers.set_xlabel("Time - 2457000 ($day$)")
        ax_outliers.set_title(f"{flux_type.replace('_', ' ')} Light Curve", fontstyle="italic")
        parent.canvas_outliers.draw()

        parent.canvas_detrending.figure.clear()
        ax_detrending = parent.canvas_detrending.figure.add_subplot(111)
        ax_detrending.errorbar(time, flux, yerr=flux_err, fmt=".", color="dodgerblue",
                               markersize=5, capsize=1, elinewidth=1.5, ls="none")
        ax_detrending.set_ylabel("Flux ($electrons^{-1}$)")
        ax_detrending.set_xlabel("Time - 2457000 ($day$)")
        ax_detrending.set_title(f"{flux_type.replace('_', ' ')} Light Curve", fontstyle="italic")
        parent.canvas_detrending.draw()

        return time, flux, flux_err, quality

    except Exception as e:
        QMessageBox.critical(parent, "Error", f"An error occurred while downloading and plotting data: {e}")
        return None, None, None, None

def phase_fold_and_normalize(time, flux, flux_err=None, t0=None, period=None, auto=False, canvas=None, display_fields=None):
    # 1. Clean Data
    mask = np.isfinite(time) & np.isfinite(flux)
    time = time[mask]
    flux = flux[mask]
    if flux_err is not None:
        flux_err = flux_err[mask]

    # 2. Normalize Flux
    median_flux = np.median(flux)
    flux = flux / median_flux
    if flux_err is not None:
        flux_err = flux_err / median_flux

    # 3. Invert Flux and Detect Minima
    inverted_flux = -flux
    prominence_value = 0.05
    height_value = -0.98
    peaks, properties = find_peaks(inverted_flux, prominence=prominence_value, height=height_value)
    minima_times = time[peaks]
    minima_fluxes = flux[peaks]

    # 5. Separate Minima Using Gaussian Mixture Model
    gmm = GaussianMixture(n_components=2)
    gmm.fit(minima_fluxes.reshape(-1, 1))
    means = gmm.means_.flatten()
    sorted_means = np.sort(means)

    # 6. Classify Minima by Flux Value
    flux_threshold = (sorted_means[0] + sorted_means[1]) / 2
    primary_indices = np.where(minima_fluxes <= flux_threshold)[0]
    primary_minima_times = minima_times[primary_indices]

    # 7. Sort Primary Minima by Time
    sorted_indices = np.argsort(primary_minima_times)
    primary_minima_times_sorted = primary_minima_times[sorted_indices]

    # 8. Set Initial t0
    t0_initial = primary_minima_times_sorted[0]

    # 9. Calculate Time Differences Between Primary Minima
    time_diffs = np.diff(primary_minima_times_sorted)
    period_estimate = np.median(time_diffs)

    # 10. Refine t0 Using Cubic Spline
    delta_t = 0.05
    mask_fit = (time >= t0_initial - delta_t) & (time <= t0_initial + delta_t)
    time_fit = time[mask_fit]
    flux_fit = flux[mask_fit]
    flux_err_fit = flux_err[mask_fit]

    spline = UnivariateSpline(time_fit, flux_fit, w=1.0 / flux_err_fit, k=3, s=0)
    time_fine = np.linspace(time_fit.min(), time_fit.max(), 10000)
    flux_fine = spline(time_fine)
    min_index = np.argmin(flux_fine)
    t0_refined = time_fine[min_index]

    # 11. Lomb-Scargle Period Analysis
    frequency_guess = 1 / period_estimate
    frequency_min = frequency_guess * 0.5
    frequency_max = frequency_guess * 1.5
    frequency = np.linspace(frequency_min, frequency_max, 10000)
    ls = LombScargle(time, flux)
    power = ls.power(frequency)
    best_frequency = frequency[np.argmax(power)]
    best_period = 1 / best_frequency

    # 13. Fourier Series Fit for Improved t0 Calculation
    def fourier_series(x, a0, a1, b1, a2, b2, a3, b3):
        return (a0 +
                a1 * np.cos(2 * np.pi * 1 * x) + b1 * np.sin(2 * np.pi * 1 * x) +
                a2 * np.cos(2 * np.pi * 2 * x) + b2 * np.sin(2 * np.pi * 2 * x) +
                a3 * np.cos(2 * np.pi * 3 * x) + b3 * np.sin(2 * np.pi * 3 * x))

    phase = ((time - t0_refined) / best_period) % 1
    extended_phase = np.concatenate([phase, phase + 1])
    extended_flux = np.concatenate([flux, flux])

    mask_phase = (extended_phase >= 0.9) & (extended_phase <= 1.1)
    phase_fit = extended_phase[mask_phase]
    flux_fit = extended_flux[mask_phase]

    initial_params = [np.mean(flux_fit)] + [0] * 6
    params, _ = curve_fit(fourier_series, phase_fit, flux_fit, p0=initial_params)

    phase_fine = np.linspace(phase_fit.min(), phase_fit.max(), 1000)
    flux_fine = fourier_series(phase_fine, *params)
    min_index = np.argmin(flux_fine)
    phase_minimum = phase_fine[min_index]
    time_minimum = t0_refined + (phase_minimum * best_period)
    t0_final = time_minimum

    # 14. Final Phase Folding and Plotting on Provided Canvas
    phase_final = ((time - t0_final) / best_period) % 1
    extended_phase_final = np.concatenate([phase_final, phase_final + 1])
    extended_flux_final = np.concatenate([flux, flux])
    extended_flux_err_final = np.concatenate([flux_err, flux_err])

    if canvas:
        canvas.figure.clear()
        ax = canvas.figure.add_subplot(111)
        ax.plot(extended_phase_final, extended_flux_final, 'k.', markersize=1)
        ax.errorbar(extended_phase_final, extended_flux_final, yerr=extended_flux_err_final, fmt='.',
                    color="dodgerblue", markersize=5, capsize=1, elinewidth=1.5)
        ax.axvline(1, color='g', linestyle='--', label="Primary Minima")
        ax.set_xlabel('Phase')
        ax.set_ylabel("Normalized Flux ($electrons^{-1}$)")
        ax.set_title('Phase-Folded Light Curve')
        ax.set_xlim(0.2, 1.2)
        ax.legend(loc=3)
        canvas.draw()

    if display_fields:
        display_fields["best_t0_display"].setText(f"{t0_refined:.6f}")
        display_fields["best_period_display"].setText(f"{best_period:.6f}")
        display_fields["period_dropdown"].clear()
        display_fields["period_dropdown"].addItem(f"Lomb-Scargle: {best_period:.6f}")
        display_fields["period_dropdown"].addItem(f"Estimate: {period_estimate:.6f}")

    return t0_refined, period_estimate, best_period, extended_phase_final, extended_flux_final, extended_flux_err_final

def phase_fold_manually(time, flux, flux_err, t0, period, canvas):

    median_flux = np.median(flux)
    flux = flux / median_flux
    if flux_err is not None:
        flux_err = flux_err / median_flux

    phase = ((time - t0) / period) % 1

    extended_phase = np.concatenate([phase, phase + 1])
    extended_flux = np.concatenate([flux, flux])
    extended_flux_err = np.concatenate([flux_err, flux_err])


    canvas.figure.clear()
    ax = canvas.figure.add_subplot(111)
    ax.errorbar(extended_phase, extended_flux, yerr=extended_flux_err, fmt='.',
                color="dodgerblue", markersize=5, capsize=1, elinewidth=1.5)
    ax.axvline(1, color='g', linestyle='--', label='Phase 1')
    ax.set_xlabel('Phase')
    ax.set_ylabel('Normalized Flux')
    ax.set_title('Phase-Folded Light Curve')
    ax.set_xlim(0.2, 1.2)
    ax.legend()
    canvas.draw()

    return extended_phase, extended_flux, extended_flux_err