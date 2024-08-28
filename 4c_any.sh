#!/bin/bash

# written by MarkEwert03 on Aug 28, 2024

# List of parameter values to iterate over
c1values=(1 0 0 -1 0 0) #(1.1x, 1.1y, 1.1y 2.1x, 2.1y, 2.1y)
c2values=(2 0 0 -2 0 0) #(1.2x, 1.2y, 1.2z 2.2x, 2.2y, 2.2z)
c3values=(3 0 0 -3 0 0) #(1.3x, 1.3y, 1.3z 2.3x, 2.3y, 2.3z)
c4values=(4 0 0 -4 0 0) #(1.4x, 1.4y, 1.4z 2.4x, 2.4y, 2.4z)
alphavalues=(1 1 1 1) #(1.1a 1.2a 1.3a 1.4a 2.1a 2.2a 2.3a 2.4a)

# SHOULD NOT NEED TO CHANGE ANYTHING BELOW THIS LINE
# -------------------------

# n = length of list. # centers = n/3
n=${#c1values[@]}

# Log file where outputs will be saved
logfile="test/output_CHANGENAME.log"

# Clear the logfile if it exists (optional)
#> "$logfile"
echo "" >> "$logfile"

# Loop through values and execute the command

for ((i=0; i<=n-3; i=i+3)) do
    # j is for indexing alphas
    j=$(((i/3) * 4))

    # center 1
    c1x=${c1values[i]}
    c1y=${c1values[i+1]}
    c1z=${c1values[i+2]}
    a1=${alphavalues[j]}

    # center 2
    c2x=${c2values[i]}
    c2y=${c2values[i+1]}
    c2z=${c2values[i+2]}
    a2=${alphavalues[j+1]}

    # center 3
    c3x=${c3values[i]}
    c3y=${c3values[i+1]}
    c3z=${c3values[i+2]}
    a3=${alphavalues[j+2]}
    
    # center 4
    c4x=${c4values[i]}
    c4y=${c4values[i+1]}
    c4z=${c4values[i+2]}
    a4=${alphavalues[j+3]}

    echo "Running: ./simple4c_any -a1 $a1 -c1 $c1x $c1y $c1z -a2 $a2 -c2 $c2x $c2y $c2z -a3 $a3 -c3 $c3x $c3y $c3z -a4 $a4 -c4 $c4x $c4y $c4z"
    ./simple4c_any -a1 "$a1" -c1 "$c1x" "$c1y" "$c1z" -a2 "$a2" -c2 "$c2x" "$c2y" "$c2z" -a3 "$a3" -c3 "$c3x" "$c3y" "$c3z" -a4 "$a4" -c4 "$c4x" "$c4y" "$c4z" >> "$logfile"
done