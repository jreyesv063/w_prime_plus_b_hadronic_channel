#!/bin/bash

# Before running this script, be sure to grant execution permissions with the following command:
# chmod +x run.sh
# Script is used with the command ./run.sh

# Variables:
processor="ttbar_h"
executor="futures"
year="2017"
nfiles="-1"
redirector="xcache"
tag="tt_CR_topXfinder"
output_type="hist"
syst="nominal"

# Samples:
samples=(
  "SingleTau"
  "TTToSemiLeptonic"
  "TTToHadronic"
  "TTTo2L2Nu"
  "ST_s-channel_4f_leptonDecays"
  "ST_t-channel_antitop_5f_InclusiveDecays"
  "ST_t-channel_top_5f_InclusiveDecays"
  "ST_tW_antitop_5f_inclusiveDecays"
  "ST_tW_top_5f_inclusiveDecays"
  "WJetsToLNu_HT-100To200"
  "WJetsToLNu_HT-200To400"
  "WJetsToLNu_HT-600To800"
  "WJetsToLNu_HT-800To1200"
  "WJetsToLNu_HT-1200To2500"
  "WJetsToLNu_HT-2500ToInf"
  "WW"
  "WZ"
  "ZZ"
  "DYJetsToLL_M-50"
)

# For run over all the samples
for sample in "${samples[@]}"; do
  if [ "$sample" == "TTTo2L2Nu" ] || [ "$sample" == "TTToHadronic" ] || [ "$sample" == "TTToSemiLeptonic" ]; then
    nsplit="10"
  else
    nsplit="5"
  fi

  python submit.py --processor "$processor" --executor "$executor" --year "$year" --nfiles "$nfiles" --redirector "$redirector" --tag "$tag" --output_type "$output_type" --syst "$syst" --sample "$sample" --nsplit "$nsplit"
done

