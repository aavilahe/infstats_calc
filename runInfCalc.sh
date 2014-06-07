#!/bin/bash
# Author: Aram Avila-Herrera
# Email: Aram.Avila-Herrera@ucsf.edu

#export __SRC_PATH="/path/to/infCalc"
export __SRC_PATH="${HOME}/dev-src/infCalc"

python ${__SRC_PATH}/infCalc.py ${@}


