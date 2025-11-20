# StressStrain
Simple stress-strain and hardening curve fitting in Julia. Mainly intended for metals/ elasto-plastic materials.

# Installation

Make a new empty folder and open Julia into it (or navigate via cd("path/to/folder") in Julia).
Then open the package promp with ']' and type:
```
    activate .
    add https://github.com/barabule/StressStrain.git

```
After installation finishes, exit package mode (with Backspace) and type:
```
    using StressStrain
```

# Usage

Type 
```
    main()
```
 to start the main GUI:
![Main GUI](/assets/main_GUI.png)

It's best to maximize the window.


Alternatively you can already add the 'data' argument to main:
```
    data = (;strain = ..., stress = ...) #data needs the strain and stress fields
    main(data)

```

You can drag and drop a file onto the main window and am import subwindow will appear:
![Import GUI](/assets/import_GUI.png)

After changing the textbox defaults, click on Import! and a plot should appear.
Red point indicate abnormal points, which you can exclude by clicking on Clean!.
After you're satisfied with the data click on Done! and the main window with open the new data, 
closing the import window.

The controls are grouped into several categories:
## Bottom
On the bottom is a slider with 2 values: toein and cut off.
*Toein is meant to eliminate the soft portion at the start of the stress strain curve where the slope
is too small due to insufficient clamping etc at the start of the measurement.
*Cutoff cut off the last portion of the true stress curve to not  use it for fitting the hardening 
portion. In general a hardening curve for metals should not have non-monotonic rising portions, 
with this option you can eliminate segments from the hardening portion that would make a worse fit to
the data.

Both values can be reset non-destructively.
## Overview
Change the material name and check if your data is true stress (checked) or 
engineering stress (unchecked).

## True Stress Curve
Here you can resample the true stress function by several methods and also change the number of points.
It's intended mainly to decimate large datasets (1000s of points) to something easier to handle (~100)
without loss of precision.

Usually on sparse  data you don't need to touch any of these. Don't touch this if you're ok with the
quality of the data (smoothness etc.).
### Manual Fit
A new window will launch if you click on *Manual:
![Manual Edit](/assets/manual_edit.png)
Here you can manually edit a piecewise cubic Bezier curve
(similar to Inkscape) over your data. Close the window when satisfied with your fit, the true stress curve
will change according to the new Bezier curve. 
After this you don't need to fit a function anymore to the hardening, because the curve should be already smooth enough.

The reset button will reset the true stress curve back to the original data, the settings will stay the same.

## E modulus
Here you can change the E modulus (automatically computed) to some desired value and freeze it once satisfied.

## Hardening Curve
Change the fitting type for the extracted hardening portion of the true stress curve.
You can choose from linear, cubic spline and several analytic functions (Swift, Voce etc).
Whenever you change the fitting, the status label (on the bottom) will change reflecting the fitted parameters.

All fitted hardening curves are extrapolatable (linear) to some desired plastic strain and you can change how many
points to export (default 100).

## Export
True stress and hardening curves are exported to comma seprated csv.
The plot is exported to PNG.
