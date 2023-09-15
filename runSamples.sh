#!/bin/bash

# Before running this script, be sure to grant execution permissions with the following command:
# chmod +x runSamples.sh
# Script is used with the command ./runSamples.sh

# Define las variables con las opciones
processor="ttbar_h"
executor="futures"
year="2017"
nfiles="-1"
redirector="xcache"
tag="tt_CR_1b1tau1top"
output_type="hist"
syst="nominal"

# Array con los nombres de los samples
samples=(

  TTToSemiLeptonic
  TTTo2L2Nu
  TTToHadronic
  
  SingleTau
  
  ST_s-channel_4f_leptonDecays
  ST_t-channel_antitop_5f_InclusiveDecays
  ST_t-channel_top_5f_InclusiveDecays
  ST_tW_antitop_5f_inclusiveDecays
  ST_tW_top_5f_inclusiveDecays
  
  WW
  WZ
  ZZ
  
  DYJetsToLL_M-50_HT-100to200
  DYJetsToLL_M-50_HT-200to400
  DYJetsToLL_M-50_HT-400to600
  DYJetsToLL_M-50_HT-600to800
  DYJetsToLL_M-50_HT-800to1200
  DYJetsToLL_M-50_HT-1200to2500
  DYJetsToLL_M-50_HT-2500ToInf


  WJetsToLNu_HT-100To200
  WJetsToLNu_HT-200To400
  WJetsToLNu_HT-400To600
  WJetsToLNu_HT-600To800
  WJetsToLNu_HT-800To1200
  WJetsToLNu_HT-1200To2500
  WJetsToLNu_HT-2500ToInf
  
)

# Loop sobre los samples y ejecutar el comando Python para cada uno de ellos
for sample in "${samples[@]}"; do
  echo "Ejecutando para el sample: $sample"
  # Establece el valor predeterminado de nsplit
  nsplit=10

  python submit.py --processor "$processor" --executor "$executor" --year "$year" --nfiles "$nfiles" --redirector "$redirector" --tag "$tag" --output_type "$output_type" --syst "$syst" --sample "$sample" --nsplit "$nsplit"
done
