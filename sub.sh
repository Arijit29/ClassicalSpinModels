#!/bin/bash

# max processors to be used :
# c12 - 240
# c11 - 240
# c09 - 192
# c10 - 70

for L in 24
do
	qsub -v l=$L -N ann.L$L -o ./log/ parm.sh 
done
