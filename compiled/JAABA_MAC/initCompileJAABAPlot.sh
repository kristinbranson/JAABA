#!/bin/bash
rm -rf JAABAPlot.app
mv JAABAPlot/src/JAABAPlot.app .
cp preprelaunch_JAABAPlot JAABAPlot.app/Contents/MacOS/preprelaunch
sed -i "" 's/prelaunch/preprelaunch/g' JAABAPlot.app/Contents/Info.plist 
