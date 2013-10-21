#!/bin/bash
# Author: Aram Avila-Herrera
# Email: Aram.Avila-Herrera@ucsf.edu

#export __SRC_PATH="/path/to/infstats"
export __SRC_PATH="/home/aram/dev-src/infstats_calc"

python ${__SRC_PATH}/infCalc.py ${@}


