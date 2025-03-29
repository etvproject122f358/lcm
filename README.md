<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
</head>
<body>

<h1>lcm: Photometric Observation and Data Analysis Tool</h1>

<p>
  <strong>lcm</strong> is a scientific software developed to interpret photometric light curves,
  integrate TESS data analysis, and facilitate efficient observation planning. Built with Python,
  it provides a user-friendly interface for both ground-based and space-based photometric data.
</p>

<hr>

<h2>1. Photometric Observation Analysis</h2>
<p>
  <strong>lcm</strong> processes measurement files reduced using <strong>AstroImageJ (AIJ)</strong> and visualizes them
  as light curves. For each observation, it displays:
</p>
<ul>
  <li><strong>SNR (signal-to-noise ratio)</strong>, source counts, and background counts</li>
  <li>Exposure time and basic observational details (e.g., date, filter, telescope)</li>
  <li>The angular distance to the Moon at the start/end of the observation</li>
</ul>
<p>
  To help observers understand ideal observational conditions, the software generates a 
  horizon altitude-time graph that accounts for twilight intervals (astronomical, nautical, 
  and civil). Because expressing the target-Moon distance only by altitude is essentially 
  one-dimensional, <strong>lcm</strong> visualizes the target’s azimuth with a color scale, providing
  a clearer view of the star’s position relative to the Moon.
</p>
<p>
  If a minima file is provided, previously calculated minimum times are overlaid on the 
  light curve, enabling a visual check of when these minima occur. This also helps validate 
  ephemeris accuracy. Minimum times are currently plotted in <strong>EEST</strong> (Eastern European 
  Summer Time); future versions may allow configurable time zones. The minimum time format:
</p>
<pre>2460758.266568452585489   1   0.000008184855781</pre>
<p>
  Here, columns represent BJD (minimum time), type (1=primary, 2=secondary), and uncertainty 
  (days). The software further computes the error in seconds to guide the observer.
</p>

<hr>

<h2>2. SIMBAD Integration</h2>
<p>
  <strong>lcm</strong> provides a direct link to the <strong>SIMBAD</strong> database, allowing users
  to quickly access detailed information about the target star. This ensures seamless
  integration of archival data (e.g., spectral type, magnitudes, coordinates) without
  leaving the <strong>lcm</strong> interface.
</p>

<hr>

<h2>3. TESS Analysis Module</h2>
<p>
  This module focuses on analyzing <strong>TESS space telescope</strong> data. It checks whether
  the target star has been observed by TESS, listing the relevant sectors in which it appears.
  Data products are categorized by cadence (<em>20s, 120s</em>) and Full Frame Images (FFIs) 
  with exposure times of <em>158.4s, 475.2s, and 1425.6s</em>. Users can:
</p>
<ul>
  <li><strong>Query</strong> TESS observations for a target</li>
  <li>Filter and visualize light curves based on <em>quality flags</em> and <em>NaN</em> values</li>
  <li>Remove outliers <em>manually</em> or using <em>Sigma Clipping</em>, <em>Boxplot analysis</em>, or <em>Chauvenet’s criterion</em></li>
  <li>Apply <strong>detrending</strong> up to 4th-degree polynomial fits</li>
  <li><strong>Phase Fold</strong> the light curves manually or using automatic methods 
      (Gaussian Mixture Models, cubic spline, Fourier analysis, Lomb-Scargle periodograms)</li>
</ul>
<p>
  At each step, data can be saved in <code>csv</code>, <code>txt</code>, or <code>xlsx</code> format for
  further analysis. 
</p>

<hr>

<h2>4. Observation Planning Module (AstPlan)</h2>
<p>
  This module helps users determine target visibility (altitude) and primary/secondary
  eclipse events for a specified date and observatory location. Observers can manually
  enter or load ephemeris data (with star names, reference minima, and orbital periods).
  Once a location is selected (from a saved list or via manual coordinates), a dynamic 
  altitude-time graph is generated, including twilight intervals. The target’s azimuth 
  is shown using a color scale, clarifying its position relative to the Moon. 
</p>
<p>
  If the user provides reference minima (T<sub>0</sub>) and a period (P), the software calculates
  the times of any minima crossing on that date, marking them on the graph in hour format 
  (currently EEST). These can help observers pinpoint the best observing window for 
  precise eclipse timing.
</p>

<hr>

<h2>Time Reference</h2>
<p>
  All observation times, plots, and minima calculations are displayed in <strong>EEST</strong>. 
  Future releases may offer localization to different time zones or direct use of UTC.
</p>

<hr>

<h2>Input Requirements</h2>
<ul>
  <li><strong>AIJ Measurement File</strong>: <code>.csv</code> or <code>.xls</code> exported from AstroImageJ, containing
    time and flux columns.</li>
  <li><strong>Minima File (optional)</strong>: A text file listing minima times in the format:
    <pre>2459223.1234   1   0.000010</pre>
    The columns represent <em>minimum time [BJD]</em>, <em>minima type (1 or 2)</em>, 
    and <em>error in days</em>.</li>
  <li><strong>Ephemeris File</strong>: A file detailing star names, reference minimum times, 
    and periods. This can be entered manually or loaded from disk.</li>
  <li><strong>Observatory Location</strong>: Coordinates can be entered manually or selected 
    from predefined templates. Additional locations can be saved for future reuse.</li>
</ul>

<hr>

<h2>License</h2>
<p>
  <strong>lcm</strong> is released under the <strong>Apache License 2.0</strong>. Refer to the <code>LICENSE</code> 
  file for details on usage and redistribution.
</p>

</body>
</html>
