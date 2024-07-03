#!/bin/bash

# pharokka.py
# Miranda Tomas

# Este script utiliza el programa Pharokka para anotar genomas de fagos.
# Para ello, hay que pasar como argumentos: (1) ficheros de entrada en formato
# FASTA con los genomas de los fagos que se quieren anotar, (2) directorio en el
# que se quieren guardar los resultados.

# Numero de argumentos que requiere el script
n=2

# Verificar si se proporcionan un numero de argumentos correct
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt> <directorio_salida>"
        exit 1
else

	# Ficheros de entrada
	input_files=$1

	# Directorio de salida 
	outdir1=$2

fi

# Base de datos de pharokka
database="/home/principal/anaconda3/envs/pharokka_env/databases"

# Bucle que recorre todos los ficheros de entrada
for file in $input_files; 
do
	name=$(basename $file .fasta)
    	outdir2=$outdir1"/"$name
    
    	echo "Analizando" $file
    	
	# Análisis con pharokka 
    	pharokka.py -i $file -o $outdir2 -d $database -t 4

done

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then

        echo "El programa ha finalido. Los ficheros resultantes de anotación se encuentran en el directorio: $outdir2"

else
        echo "El programa ha finalizado con errores."
fi

