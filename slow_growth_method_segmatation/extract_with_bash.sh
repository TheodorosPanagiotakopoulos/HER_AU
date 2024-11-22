#!/bin/bash

rm grad.dat

i=1

while [ $i -le 1000 ]
do
  if test -f REPORT.$i
  then
    grep cc REPORT.$i |awk '{print $3}' >xxx
    grep b_m REPORT.$i |awk '{print $2}' >fff
    paste xxx fff >> grad.dat
  fi
  let i=i+1
done


sort -n grad.dat >grad_.dat
mv grad_.dat grad.dat
rm xxx
rm fff
