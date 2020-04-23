#!/bin/bash

NUMSAMPLES="4096"
DESTINATION="dataset/"


XS=("0.5" "0.4924" "0.4698" "0.433" "0.383" "0.3214" "0.25" "0.171" "0.0868" "0.0")
YS=("0.866" "0.866" "0.866" "0.866" "0.866" "0.866" "0.866" "0.866" "0.866" "0.866")
ZS=("0.0" "0.0868" "0.171" "0.25" "0.3214" "0.383" "0.433" "0.4698" "0.4924" "0.5")
ENVMAPS=("0" "10" "20" "30" "40" "50" "60" "70" "80" "90")
MESHMODELS=("cube" "bun" "bunny" "lucy" "star_smooth" "armadillo" "dragon" "soap" "buddha" "cap" "bust")
SIGMATS=("30" "36" "44" "54" "67" "82" "100" "122" "150" "184" "225" "276")
ALBEDOES=("0.39" "0.59" "0.74" "0.87" "0.95")
GS=("0.0" "0.2" "0.4" "0.6" "0.8")

for MESHMODEL in "${MESHMODELS[@]}"; do
	COUNTER=0
        for ENVMAP in "${ENVMAPS[@]}"; do
                X="${XS[COUNTER]}"
                Y="${YS[COUNTER]}"
                Z="${ZS[COUNTER]}"
		let COUNTER=COUNTER+1
                for SIGMAT in "${SIGMATS[@]}"; do
                        for ALBEDO in "${ALBEDOES[@]}"; do
                                for G in "${GS[@]}"; do
                                        FILE="${MESHMODEL}_e${ENVMAP}_d${SIGMAT}_a${ALBEDO}_g${G}_q${NUMSAMPLES}"
                        		qsub -V -o "${FILE}.o${JOB_ID}" -e "${FILE}.e${JOB_ID}" -b y -cwd -v X="${X}",Y="${Y}",Z="${Z}",DESTINATION="${DESTINATION}",NUMSAMPLES="${NUMSAMPLES}",MESHMODEL="${MESHMODEL}",ENVMAP="${ENVMAP}",SIGMAT="${SIGMAT}",ALBEDO="${ALBEDO}",G="${G}" "./run_${MESHMODEL}_job.sh"
                                done
                        done
                done
                
        done
done

