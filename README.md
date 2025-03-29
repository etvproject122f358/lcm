<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
</head>
<body>

<h1>lcm: Photometric Observation and Data Analysis Tool</h1>

<p><strong>lcm</strong> is a scientific software developed for the visualization and interpretation of photometric light curves, integrating multiple observational analysis tools into a user-friendly structure. The program includes modules for TESS space telescope data analysis, ground-based observation assessment using AstroImageJ (AIJ) outputs, and an astronomical observation planning system based on ephemeris data. It is designed to assist both in analyzing existing data and planning future observations.</p>

<h2>Main Features</h2>
<ul>
  <li>Visualization of photometric light curves from TESS and ground-based telescopes</li>
  <li>Automatic identification and marking of primary and secondary minima on light curves</li>
  <li>Calculation of signal-to-noise ratios (SNR), source and background counts, and exposure diagnostics</li>
  <li>Outlier detection and removal using statistical methods (Sigma Clipping, Boxplot, Chauvenet's criterion)</li>
  <li>Normalization and detrending of light curves using polynomial fits</li>
  <li>Manual and automated phase folding using Gaussian Mixture Models, Fourier analysis, cubic spline interpolation, and Lomb-Scargle periodogram</li>
  <li>Observation planning with altitude-time plots and Moon angular separation visualizations</li>
</ul>

<h2>Modules</h2>

<h3>1. TESS Analysis Module</h3>
<p>This module allows the user to query whether a target star has been observed by the TESS telescope and in which sectors. Data products (short/long cadence or Full Frame Images - FFIs) are listed and categorized by cadence type (20s, 120s) and exposure duration (158.4s, 475.2s, 1425.6s). The user can filter by quality flags, visualize light curves, and optionally remove outliers. Polynomial detrending (up to 4th degree) can be applied. Phase folding is available via both manual inputs and automatic period estimation techniques. The reference time and period are displayed on the light curve, with minima positions indicated by markers.</p>

<h3>2. Photometric Observation Analysis Module</h3>
<p>This module uses measurement files reduced with AstroImageJ (AIJ). Users must upload the measurement file (.csv/.xls) and may optionally provide a file containing predicted minima. The tool displays the light curve with overlaid SNR, source and background counts, and exposure time diagnostics. The minimum file should contain lines in the following format:</p>
<pre>2460758.266568452585489    1    0.000008184855781</pre>
<p>Where the columns indicate: Minimum time [BJD], Minimum type (1=Primary, 2=Secondary), and Uncertainty [days].</p>
<p>If provided, the minima are displayed on the plot and compared to the observational data. The angular separation between the target star and the Moon is also calculated and presented at the beginning and end of the observation session.</p>

<h3>3. Astronomical Observation Planning Module</h3>
<p>This module calculates the visibility of stars on a given date and location. Users can enter ephemeris data manually or load from a file. Uploaded ephemeris files are parsed dynamically and allow searching by star name. For each target, the tool calculates the times of primary and secondary minima based on reference minima and orbital period, and plots altitude-time diagrams. Twilight times are also considered. The target's azimuth values are color-mapped to better visualize spatial orientation, enabling clearer understanding of angular separation from the Moon. The system allows saving observatory locations for later use.</p>

<h2>Time and Units</h2>
<p>The software currently uses EEST (Eastern European Summer Time) as the default time reference. All time-based outputs (including minima, plots, and observation windows) are adjusted accordingly. Support for alternate time zones or UTC-based outputs will be considered in future versions.</p>

<h2>Input Requirements</h2>
<ul>
  <li><strong>AIJ Measurement File:</strong> A .csv or .xls file exported from AstroImageJ, containing columns for time and flux.</li>
  <li><strong>Minima File (optional):</strong> A .txt file with each line containing:
    <pre>Minimum Time [BJD]    Min Type (1/2)    Uncertainty [days]</pre>
  </li>
  <li><strong>Ephemeris File:</strong> A list of target names with associated reference minima and orbital period.</li>
  <li><strong>Observatory Location:</strong> Coordinates can be entered manually or selected from pre-defined templates.</li>
</ul>

<h2>License</h2>
<p>This software is released under the <strong>Apache License 2.0</strong>. Please refer to the LICENSE file for detailed terms of use.</p>

</body>
</html>
