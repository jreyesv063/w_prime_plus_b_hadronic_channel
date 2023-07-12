#!/bin/bash

# Definición de variables comunes
processor="ttbar_cr1"
first_executor="iterative"
channel="mu"
nfiles=-1
nsplit=5
tag="tt_CR1"

date="2023-07-10"
year="2017"

output="outfiles/$tag/$date/$processor/$2017/$channel"

second_executor="futures"


echo -e '\n'

echo ':::::::::::::::::::::::::::::::::::::::::::::'
echo '  Running samples   '
echo ':::::::::::::::::::::::::::::::::::::::::::::'


# Lista de muestras Single Top
samples_ST=(
    "ST_s-channel_4f_leptonDecays"
    "ST_t-channel_antitop_5f_InclusiveDecays"
    "ST_t-channel_top_5f_InclusiveDecays"
    "ST_tW_antitop_5f_inclusiveDecays"
    "ST_tW_top_5f_inclusiveDecays"
)


# Lista de muestras Diboson
samples_VV=(
    "WW"
    "WZ"
    "ZZ"
)


# Lista de muestras w+jets
samples_WJ=(
    "WJetsToLNu_HT-100To200"
    "WJetsToLNu_HT-200To400"
    "WJetsToLNu_HT-400To600"
    "WJetsToLNu_HT-600To800"
    "WJetsToLNu_HT-800To1200"
    "WJetsToLNu_HT-1200To2500"
    "WJetsToLNu_HT-2500ToInf"       
)

# Lista de muestras drell-yan
samples_DY=(
    "DYJetsToLL_M-50_HT-70to100"
    "DYJetsToLL_M-50_HT-100to200"
    "DYJetsToLL_M-50_HT-200to400"
    "DYJetsToLL_M-50_HT-400to600"
    "DYJetsToLL_M-50_HT-600to800"
    "DYJetsToLL_M-50_HT-800to1200"
    "DYJetsToLL_M-50_HT-1200to2500"
    "DYJetsToLL_M-50_HT-2500toInf"
)


# Lista de muestras tt
samples_TT=(
    "TTTo2L2Nu"
    "TTToHadronic"
    "TTToSemiLeptonic" 
)

# Lista de muestras de data
samples_DATA=(
    "SingleMuon"
)


# Función para ejecutar el comando de Python y verificar la existencia del archivo .pkl
run_and_check() {
    local sample="$1"
    local executor="$2"
    local file="$sample.pkl"

    run_python_command "$sample" "$executor"

    if [ -f "$output/$file" ]; then
        echo -e "\t The $file file was successfully created"
    else
        run_python_command "$sample" "$second_executor"
    fi
}

# Función para ejecutar el comando de Python
run_python_command() {
    local sample="$1"
    local executor="$2"
    python3 submit.py --processor "$processor" --executor "$executor" --channel "$channel" --sample "$sample" --nsplit "$nsplit" --nfiles "$nfiles" --tag "$tag"
}


# Bucle for para ejecutar y corroborar los comandos para cada muestra Single Top

echo -e '\n'
echo '*** Running Single Top samples ***'

for sample in "${samples_ST[@]}"
do
    run_and_check "$sample" "$first_executor"
done


####


echo -e '\n'
echo '*** Running Diboson samples ***'

# Bucle for para ejecutar y corroborar los comandos para cada muestra Diboson

for sample in "${samples_VV[@]}"
do
    run_and_check "$sample" "$first_executor"
done


####


echo -e '\n'
echo '*** Running W+Jets samples ***'

# Bucle for para ejecutar y corroborar los comandos para cada muestra Diboson

for sample in "${samples_WJ[@]}"
do
    run_and_check "$sample" "$first_executor"
done


####


echo -e '\n'
echo '*** Running Drell-Yan+Jets samples ***'

# Bucle for para ejecutar y corroborar los comandos para cada muestra Diboson

for sample in "${samples_DY[@]}"
do
    run_and_check "$sample" "$first_executor"
done


####


echo -e '\n'
echo '*** Running ttbar samples ***'

# Bucle for para ejecutar y corroborar los comandos para cada muestra Diboson

for sample in "${samples_TT[@]}"
do
    run_and_check "$sample" "$first_executor"
done



####


echo -e '\n'
echo '*** Running ttbar samples ***'

# Bucle for para ejecutar y corroborar los comandos para cada muestra Diboson

for sample in "${samples_DATA[@]}"
do
    run_and_check "$sample" "$first_executor"
done




echo ' All the samples have been completed'
echo ':::::::::: Thank you! :::::::::::::::::::'