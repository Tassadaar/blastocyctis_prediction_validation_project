#!/bin/bash

source activate blastocyctisValidator

python false_prediction_validator.py -f input_files/ST7B_GCA_00151665.1.fna_shortened_headers.fna -p input_files/test.gff3

conda deactivate
