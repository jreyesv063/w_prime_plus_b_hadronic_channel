#!/bin/bash

# Before running this script, be sure to grant execution permissions with the following command:
# chmod +x runSamples.sh
# Script is used with the command ./runSamples.sh

# Define las variables con las opciones
processor="ttbar_h"
executor="dask"
year="2017"
nfiles="-1"
redirector="xcache"
tag="tt_CR_1b1tau1top"
output_type="hist"
syst="nominal"

# Array con los nombres de los samples
samples=(
  TTToSemiLeptonic
  DYJetsToLL_M-50_HT-70to100
  DYJetsToLL_M-50_HT-100to200
  DYJetsToLL_M-50_HT-200to400
  DYJetsToLL_M-50_HT-400to600
  DYJetsToLL_M-50_HT-600to800
  DYJetsToLL_M-50_HT-800to1200
  DYJetsToLL_M-50_HT-1200to2500
  DYJetsToLL_M-50_HT-2500ToInf
  ST_s-channel_4f_leptonDecays
  ST_t-channel_antitop_5f_InclusiveDecays
  ST_t-channel_top_5f_InclusiveDecays
  ST_tW_antitop_5f_inclusiveDecays
  ST_tW_top_5f_inclusiveDecays
  TTTo2L2Nu
  TTToHadronic
  WJetsToLNu_HT-100To200
  WJetsToLNu_HT-200To400
  WJetsToLNu_HT-400To600
  WJetsToLNu_HT-600To800
  WJetsToLNu_HT-800To1200
  WJetsToLNu_HT-1200To2500
  WJetsToLNu_HT-2500ToInf
  WW
  WZ
  ZZ
  SingleTau
)

# Loop sobre los samples y ejecutar el comando Python para cada uno de ellos
for sample in "${samples[@]}"; do
  echo "Ejecutando para el sample: $sample"
  # Establece el valor predeterminado de nsplit
  nsplit=5

  # Ajusta nsplit según el nombre del sample
  case "$sample" in
    DYJets*|WJets*|WW|WZ|ZZ)
      ;;
    TTTo*|SingleTau)
      ;;
    ST*)
      ;;
  esac

  python submit.py --processor "$processor" --executor "$executor" --year "$year" --nfiles "$nfiles" --redirector "$redirector" --tag "$tag" --output_type "$output_type" --syst "$syst" --sample "$sample"
done
