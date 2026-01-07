#!/bin/bash
# ======================================
#
# Start_Casper_jobs.sh
#
# This scripts submits individual ASR MCS simulations
# to perform feature tracking on the cloud and precipitation
# field
#
#
# Dependencies:
# none
#
# ======================================
# ======================================
# USER Input

DataDir='/glade/derecho/scratch/bblanc/ERA5_hail_model/ERA5-hailpredictors/'

ST_DT='1990-01-01'
EN_DT='2022-12-01'   
endt=$(date '+%s' -d "$EN_DT")
i="$ST_DT"


while [[ $(date +%s -d $i) -le $endt ]]; do
   # echo "${i%-*}"
   MM=$(date -d "$i" '+%m')
   YYYY=$(date -d "$i" '+%Y')

   # check if file exists?
   FILE=$DataDir$YYYY$MM'_ERA-5_HailPredictors_newSRH03.nc'
   if [ ! -f $FILE ]; then
       echo $YYYY$MM
       sed "s/MM/$MM/g" CasperSubmit_newSRH03.sh > CasperSubmit_fin_newSRH03.sh
       sed -i "s/YYYY/$YYYY/g" CasperSubmit_fin_newSRH03.sh
       qsub CasperSubmit_fin_newSRH03.sh
       # sleep 1s
   fi
   i=$(date '+%Y-%m-%d' -d "$i +1 month")
done



exit



