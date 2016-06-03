#!/bin/bash
$path="Input/"
for i in $(ls Intput)
do
   file=$path$i
   Rscript runIt.R $file params.csv
done
