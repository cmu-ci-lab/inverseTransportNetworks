#!/bin/bash
#DESTINATION="${1}"
#NUMSAMPLES="${2}"
#HEIGHTMAP="${3}"
#HEIGHTSCALE="${4}"
#DENSITY="${5}"
FILENAME="${DESTINATION}${MESHMODEL}_e${ENVMAP}_d${SIGMAT}_a${ALBEDO}_g${G}_q${NUMSAMPLES}.exr"
echo -e "${FILENAME}"

mitsuba ../scenes/lucy_sunsky.xml -Dx="${X}" -Dy="${Y}" -Dz="${Z}" -Dmeshmodel="${MESHMODEL}" -DsigmaT="${SIGMAT}" -Dalbedo="${ALBEDO}" -Dg="${G}" -DnumSamples="${NUMSAMPLES}" -o "${FILENAME}"
