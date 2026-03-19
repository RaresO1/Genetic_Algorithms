#!/bin/bash

for (( k=1; k<=3; k++ ))
do

    for (( i=1; i<=3; i++ ))
    do
        for (( j=1; j<=12; j++ ))
        do
            ./main.exe $k &
            sleep 1;
        done
        wait
    done

done
