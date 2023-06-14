#!/bin/bash

source activate blastocyctisValidator

python false_prediction_validator.py -f input_files/ST7B_GCA_00151665.1.fna_shortened_headers.fna -p input_files/1_ST7C_to_ST7B_transfer_additional_orfs.gff3 -o test.gff3

conda deactivate
