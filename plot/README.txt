LOADING EXPERIMENTS

Multiple experiments can be divided into an arbitrary number of groups for
multi-variate analysis.  Use the "New" button to create a new group and
give it a name, then the "Add" button to browse directories for experiments.
Remove selected experiments with the "Delete" button.  Deleting all of the
experiments in a group will also delete the group.  Transfer experiments from
one group to another with the "Move" button.  Designate loaded experiments
to analyze by selecting (a subset of) them with <shift/ctrl> left button.

The features and behaviors analyzed are only those common to all of the
current experiments.

Experiment lists and other settings are maintained from one session to
the next by automatically saving the configuration to disk when exiting
("most_recent_config.mat").  Multiple different configurations can be
maintained using the "Load" and "Save" buttons.  The "Reset" button deletes
all groups and experiments and returns all settings to their defaults.


BEHAVIOR STATISTICS

Create bar charts of the fraction of time each behavior is performed using
the "Behavior Bar Chart" button.  Each panel is a behavior.  A contextual
menu on the "Behavior Bar Chart" button controls whether averages are done
on a per-fly or per-frame basis.

Logical combinations of behaviors can be analyzed using the pull-down menu
in the "Behavior" box to the left.  So for example, if "And" "Behavior <X>"
is chosen, then the bars in the panels show the fraction of time when each
behavior and Behavior <X> are performed simultaneously.

Plots for just a single individual or a specific sex can be chosen using
the pull-down mneu in the "Individual" box.  If "All" is specified here,
then all individuals from the experiments selected in the "Experiments"
box are analyzed.

The table is filled with the p-values from K-S normality tests and two-sided
Wilcoxen rank sum tests.


The fraction of time a selected behavior is performed as a function of time
during the experiment can be plotted using the "Behavior Time Series" button.
A contextual menu controls whether error bars are plotted or per-experiment
data overlayed.  Flies are first pooled within an experiment and then combined
across experiments according to the central tendency and dispersion measures
specified in a shared contextual menu on the "Prefs" button.  The experiment
groups are color coded to match that in text of the group pull-down menu.


FEATURE STATISTICS

Histograms of feature values during behaviors of interest can be created by
selecting the desired behavior, feature, and individuals using the pull-down
menus on the left and clicking on the "Feature Histogram" button.  As with
behavior statistics, the experiment groups are color coded, a contextual
menu controls the the specific calculations, and logical operations on
behaviors can be performed.

A time series of feature values can be plotted with the "Feature Time
Series" button.  Colors and contextual buttons are as before, but not logical
operators.  The time axis can either be absolute time during the experiment,
or relative to the beginning or end of a bout of the specified behavior.


INTERESTING DIFFERENCES

Create a table of behavior-feature pairs whose feature histograms are most
different using the "Interesting Feature Histograms" button.  The metric used
is the difference in the means normalized by the combined standard deviations
(i.e. d-prime of signal detection theory).  Differences between both during
versus not-during the behavior, as well as during the behavior for each
possible pair of groups are considered.  Selecting a cell in the table changes
the "Behavior" and "Feature" boxes and plots the corresponding histograms.
A contextual menu provides options to include or exclude Nan and Inf d-primes.

Still in alpha:  a similar table of behavior-feature pairs whose feature
time-series are most different can be created using the "Interesting
Feature Time Series" button.  For univariate data, the metric is the during
vs. not during the behavior difference between the root-mean-square of the
feature values.  For bivariate data, the mean of the feature values not
during behavior is first subtracted from each experiment list and then the
difference of the root-mean-square of the feature values during the behavior
between the two lists is used.  Selecting a cell plots the time series below.

Both of these calculations take a (really long) while and so the results
are cached to disk.


BOUT STATISTICS -- not converted to multivariate yet

Create a table of bout lengths (BL) and inter-bout lengths (IBL) using the
"Bout Stats" button.  As for behavior stats, each row is a behavior;  the
columns break it down by sex and individual;  the experiment groups are
striped;  logical operators can be applied; and the Prefs contextual menu
controls statistics.

Selecting a cell plots a histogram, normalized to unit area, of bout and
inter-bout lengths.  The "LogY" and "Stats" buttons can be used to scale
the axis and overlay summary statistics, respectively.  Multiple selected
cells are overlayed.


SOCIAL STATISTICS -- not converted to multivariate yet

Create a table of the closest individual using the "Social Stats" button.
Choose which behavior to analysis using the pull-down menu in the Behavior
panel, and which metric to use for closest using the pull-down menu in the
Feature panel.  In the resulting table, each row is an individual, and the
columns show the mode and percentage of frames of the closest individual
for each bout, as well as the overall mode of the per-bout modes.

Selecting a cell in the table plots a histogram of the closest indidividual.


KNOWN ISSUES



TODO
