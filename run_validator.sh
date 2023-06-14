#!/bin/bash

# you will need biopython and gffutils in your conda environment
conda_env='your_conda_env'

source activate $conda_env

assembly='input_files/ST7B_GCA_00151665.1.fna_shortened_headers.fna'
prediction='input_files/1_ST7C_to_ST7B_transfer_additional_orfs.gff3'
output='output.gff3'
mode='all'

python false_prediction_validator.py -f $assembly -p $prediction -o $output -m $mode

conda deactivate
