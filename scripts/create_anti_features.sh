#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "[ERROR] unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <feature BED> <chromosome lengths> <output BED>"
    exit
fi


bedtools complement -i $1 -g $2 > $3
