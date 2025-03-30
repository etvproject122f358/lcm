<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
</head>
<body>

<h1>lcm: A Multi-Module System for Photometric Data Analysis and Planning</h1>

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
  <li>Quick access to SIMBAD for star name resolution and basic astrophysical data</li>
</ul>

<h2>Modules</h2>

<h3>1. Photometric Observation Analysis Module</h3>
<p>This module processes measurement files reduced with AstroImageJ (AIJ). Users must upload the measurement file (.csv/.xls) and may optionally provide a file containing predicted minima. The light curve is visualized with overlays showing the signal-to-noise ratio (SNR), source/background counts, and exposure time for each observation. A dedicated diagnostics section helps identify ideal exposure and count values.</p>
The software currently visualizes minimum times based on EEST (Eastern European Summer Time), with conversion applied to match local time. Future updates will include configurable time zones.</p>
<p>Additionally, the module calculates the angular separation between the target and the Moon at the beginning and end of the observation, providing further insights for planning and evaluation.</p>


<h3>2. TESS Analysis Module</h3>
<p>This module queries whether a target has been observed by the TESS space telescope and retrieves data across available sectors. The user is presented with a categorized list of observations based on cadence: 20s and 120s short/long cadence data, as well as Full Frame Images (FFIs) with exposures of 158.4s, 475.2s, and 1425.6s. Quality flags can be applied to filter the data. Users may visualize light curves, identify and remove outliers (using manual and automatic techniques), and apply polynomial detrending up to the 4th degree.</p>
<p>The phase folding section enables both manual entry of reference time and period, as well as automatic detection using techniques such as Lomb-Scargle periodograms, Gaussian Mixture Models, Fourier analysis, and cubic spline interpolation.</p>

<h3>3. Astronomical Observation Planning Module</h3>
<p>This module calculates the observability of target stars for a given date and observatory location. Observers can manually input ephemeris data or load from a formatted file. Ephemeris files are dynamically parsed and searchable by star name. For each target, the tool uses the reference minimum and period to compute the times of primary and secondary minima on the selected date.</p>
<p>The module dynamically generates altitude-time plots, taking twilight times into account. The azimuth position of the target is visualized using a color bar, aiding in the interpretation of sky position and angular separation from the Moon. Observatory locations can be manually saved and reused across sessions. The graph also includes the exact timing of predicted minima, shown in hour format for observation planning.</p>

<h2>Input Requirements</h2>
<ul>
  <li><strong>AIJ Measurement File:</strong> lcm expects measurement file exported from AstroImageJ.</li>
  <li><strong>Minima File (optional):</strong> A .txt file with each line containing:
    <pre>Minimum Time [BJD]    Min Type (1/2)    Uncertainty [days]</pre>
  </li>
  <li><strong>Ephemeris File:</strong> A list of target names with associated reference minima and orbital periods. The ephemeris information can either be entered manually or loaded from a file.</li>
  <li><strong>Observatory Location:</strong> Coordinates can be entered manually or selected from pre-defined templates.</li>
<li><strong>.dat File:</strong> Each <code>.dat</code> file must contain BJD, relative flux, and relative flux error in that exact order, followed by a corresponding measurement table (if applicable) that details the observation parameters (e.g., exposure time, star counts, background counts). Additionally, its filename must follow the format <code>YYYYMMDD_StarName_Telescope_Filter_LC</code>. For example, <code>20240101_DD-CrB_T100_u_LC.dat</code>. If a minima file is provided for the same dataset, its name should be identical but end with <code>_minima</code>, e.g., <code>20240101_DD-CrB_T100_u_LC_minima.dat</code>.</li>

</ul>

<h2>License</h2>
<p>This software is released under the <strong>Apache License 2.0</strong>. Please refer to the LICENSE file for detailed terms of use.</p>

</body>
</html>
