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

declare -a area=("4")
declare -a Plow=("4")
declare -a Phig=("26")



for yy in {2021..2022};
do
    for aa in "${area[@]}"
    do
        for pl in "${Plow[@]}"
        do
            for ph in "${Phig[@]}"
            do
                echo "$aa"'-'"$pl"'-'"$ph"

                sed "s/area/$aa/g" CasperSubmit.sh > CasperSubmit_fin.sh
                sed -i "s/qT/$pl/g" CasperSubmit_fin.sh
                sed -i "s/qD/$ph/g" CasperSubmit_fin.sh
                sed -i "s/YYYY/$yy/g" CasperSubmit_fin.sh
                qsub CasperSubmit_fin.sh

            done
        done
    done
done

exit



