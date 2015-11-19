#!/bin/bash
#rm -rf PrepareJAABAData.app
#mv PrepareJAABAData/src/PrepareJAABAData.app .
cp preprelaunch_PrepareJAABAData PrepareJAABAData.app/Contents/MacOS/preprelaunch
sed -i "" 's/prelaunch/preprelaunch/g' PrepareJAABAData.app/Contents/Info.plist 
