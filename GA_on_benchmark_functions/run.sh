#!/bin/bash

for (( i=1; i<=3; i++ ))
do
    for (( j=1; j<=12; j++ ))
    do
    ./main.exe &
    done
    wait
done
