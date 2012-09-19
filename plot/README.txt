LOADING EXPERIMENTS

Multiple experiments can be divided into two groups for uni- or bi-variate
analysis.  Use the "Add" button to invoke a dialog box to browse directories
for experiments.  Remove selected experiments with the "Delete" button.
Transfer experiments from one group to the other with the "Move Below" and
"Move Above" buttons.  Designate loaded experiments to analyze by selecting
(a subset of) them with <shift/ctrl> left button.

The features and behaviors analyzed are only those common to all of the
current experiments.

Experiment lists and other settings are maintained from one session to
the next by automatically saving the configuration to disk when exiting
("most_recent_config.mat").  Multiple different configurations can be
maintained using the "Load" and "Save" buttons.  The "Reset" button does
the obvious.


BEHAVIOR STATISTICS

Create a table of the fraction of time each behavior is performed using the
"Behavior Stats" button.  Each row is a behavior;  the columns break it
down by sex and individual.  If both experiment lists are populated then
the data are alternately striped and colored.  A contextual menu on the
"Behavior Stats" button controls whether averages are done on a per-fly
or per-frame basis.  Mean, median, and variance can be set by a contextual
menu on the "Prefs" button.  Data can be exported to a tab-delimited file
using the "Export Table" button.

Logical combinations of behaviors can be analyzed using the pull-down menu
in the lower-left "Behavior" box.  So for example, if "And" "Behavior X"
is chosen, then the rows in the table show the fraction of time when each
behavior and Behavior X are performed simultaneously.

The fraction of time each behavior is performed as a function of time
during the experiment can be plotted by selecting a cell in the table for
the particular behavior and sex/individual of interest.  Multiple cells
selected simultaneously will be overlayed.  Zoom and pan using the buttons
provided.  Use the "Stats" button to overlay mean, median, and variance data.


BOUT STATISTICS

Create a table of bout lengths (BL) and inter-bout lengths (IBL) using the
"Bout Stats" button.  As for behavior stats, each row is a behavior;  the
columns break it down by sex and individual;  the experiment groups are
striped;  logical operators can be applied; and the Prefs contextual menu
controls statistics.

Selecting a cell plots a histogram, normalized to unit area, of bout and
inter-bout lengths.  The "LogY" and "Stats" buttons can be used to scale
the axis and overlay summary statistics, respectively.  Multiple selected
cells are overlayed.


SOCIAL STATISTICS

Create a table of the closest individual using the "Social Stats" button.
Choose which behavior to analysis using the pull-down menu in the Behavior
panel, and which metric to use for closest using the pull-down menu in the
Feature panel.  In the resulting table, each row is an individual, and the
columns show the mode and percentage of frames of the closest individual
for each bout, as well as the overall mode of the per-bout modes.

Selecting a cell in the table plots a histogram of the closest indidividual.


FEATURE STATISTICS

Histograms of feature values during behaviors of interest can be created by
selecting the desired behavior, feature, and individuals using the pull-down
menus in the lower left panels and clicking on the "Feature Histogram" button.
The thick line corresponds to feature values during the behavior;  the
thin line to values not during the behavior;  the two experiment lists
are color coded.  A contextual menu in this button controls whether the
calculations are done on a per-frame or per-bout basis.  Logical operations
on behaviors can be performed as described earlier.  The "LogY" and "Stats"
buttons work in this context as well.

A time series of feature values near the beginning and end of bouts can be
similarly plotted with the pull-down menus and "Feature Time Series" button.
A contextual menu on this button controls various parameters.  The thick
colored line corresponds to the central tendency;  the thin colored lines
to the variance;  and the black lines to the raw data.  


INTERESTING DIFFERENCES

Create a table of behavior-feature pairs whose histograms are most different
between the two experiment lists using the "Interesting Histograms" button.
The metric used is the difference in the means normalized by the combined
standard deviations (i.e. d-prime of signal detection theory).  Selecting
a cell in the table plots its histogram as with bout stats.

Still in alpha:  a similar table of behavior-feature pairs whose time-series
are most different can be created using the "Interesting Time Series" button.
For univariate data, the metric is the during vs. not during the behavior
difference between the root-mean-square of the feature values.  For bivariate
data, the mean of the feature values not during behavior is first subtracted
from each experiment list and then the difference of the root-mean-square
of the feature values during the behavior between the two lists is used.
Selecting a cell plots the time series below.

Both of these calculations take awhile and so the results are cached to disk.


KNOWN ISSUES



TODO

normality and significance tests
implement export graph button
add logX button for histograms
