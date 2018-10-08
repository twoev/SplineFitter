# SplineFitter
Fit B-Splines to a set of data, and optionally estimate a signal

Author: James Monk

This uses basis-splines to fit an input background data.  It can optionally include a signal template in the fit to estimate signal/background.
Splines are defined by a user-given knot vector, and the control points of the splines, as well as fit uncertainties, can be estimated.
The idea is to be a model-agnostic way of estimating backgrounds (or signals) on a smooth curve.

To build you need to do

./configure --prefix=/path/to/install <other options>
make
make install


./configure --help tells you the full set of options.  

You will need GSL and BOOST to build the library, and ROOT to build the example
