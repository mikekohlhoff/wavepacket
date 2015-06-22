#!/bin/sh

# run wpp code for a range of velocities
# range of velocities

declare -a velRange=('v0 = -3d-4' 'v0 = -3d-5' 'v0 = -3d-6') 
echo velRange[1] 

echo 'old velocity:'
sed '31q;d' config.info
sed -i '' "31s/.*/testvel/" config.info
echo 'new velocity:'
sed '31q;d' config.info

echo ${velRange[@]}
echo ${#velRange[@]}
for i in ${#velRange[@]}
do echo line
done

