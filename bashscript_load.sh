#!/bin/bash

source $(dirname $(readlink -f $0))/build/config.sh 
module load gnu/gcc/6.5
module load fairroot/18.00
STRING="environment loaded $(readlink -f $0)"
echo $STRING
