#!/bin/bash
cp preprelaunch JAABA.app/Contents/MacOS
sed -i "" 's/prelaunch/preprelaunch/g' JAABA.app/Contents/Info.plist 
