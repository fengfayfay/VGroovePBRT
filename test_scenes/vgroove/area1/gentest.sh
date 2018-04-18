#!/bin/sh
b=64
for a in 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9 
do
    for bounce in 1 3
    do
        #val=`expr $b \* $bounce`
        val=32
        sed -e '1,$s/bounceCount/'$bounce'/g' -e '1,$s/alpha/'$a'/g' -e '1,$s/sampleCount/'$val'/g' area1_beck.pbrt > test\_b$bounce\_a$a.pbrt
        ~/work/bin/VGroove/pbrt test\_b$bounce\_a$a.pbrt &
    done
done
