#!/bin/bash
$path="Input/"
for i in $(ls Input)
do
   file=$path$i
   time Rscript runIt.R $file params.csv
done
