** Starting up

1) Change matlab directory to JAABA/miceCode.
2) Run "setuppaths".


** Creating JAABA directories

To use JAABA, you need to create a JAABA directory for each movie. 
makeMovies(dirname) command creates the directories. The function runs
recursively, so it'll create irrespective of how deep inside the dirname
the movie exists. It'll create a directory with the same name as the movie
(which has to be an .avi) and create a copy of the movie with the name 
movie.avi in the directory.


** Preparing a set of mice videos with the same setup

Run command annotateAllDirs(dirname), where dirname is the directory that 
has the jaab folders for the movies. Make sure to generate the features (below) 
after annotating the videos

** Generating features for mice videos.

After annotating a set of videos, generate features for them by running
genAllFeatures(dirname). This is the most time consuming step (around 5 min 
for each video), so it is better done overnight.

** Creating, annotating and generating features in one step

All the above 3 steps can also be done in one step using setUpDir(dirname).


** To Start JAABA

1) Change matlab directory to JAABA/miceCode.
2) Run StartJAABA

** Creating a new project in JAABA 

1) Select "New" when you start JAABA.
2) Select "Adam_mice" from the animal type drop down menu.
3) Give the behavior name and say Done.
4) Once done with labeling, save the project as a jab file.

** Opening an existing project

1) After starting JAABA, select "Open an existing project" 
and browse to previously saved jab file

** Labeling in JAABA

1) Details are in JAABA documentation. (http://jaaba.sourceforge.net/)

** Training/Predicting

1) Use the Train/Predict buttons.

** To plot ethograms

1) ethogram_plot(expdirs,jabfiles,nframesplot)
 You can double click on the ethogram to open JAABA with the corresponding
experiment and behavior detector.

** Classifying (Predicting) a set of movies

1) classifyAllDirs(rootdir,jabfile)

This function classifies all the folders inside the rootdir.

** Closing All Windows when something fails

 Run closeEverything.m

** Clearing a directory to start fresh

Run clearDirectory(rootdir)

This will remove all the subfolders that have movie.avi within the root directory. Be careful while using it.
