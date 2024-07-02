#!/bin/bash
#SBATCH -n 1            	# Number of cores requested
#SBATCH -N 1            	# Number of nodes requested
#SBATCH -t 2:00:00    		# Runtime in minutes
#SBATCH --qos medium   		# The QoS to submit the job
#SBATCH --mem=10G      		# Memory per cpu in G
#SBATCH -o virsorter2.out  	# Standard output goes to this file
#SBATCH -e virsorter2.err  	# Standard error goes to this file

# virsorter2.sh
# Miranda Tomas 

# Este script utiliza el programa VirSorter2 para clasificar los contigs virales
# de fagos de DNA de doble cadena a partir de ficheros de entrada en formato FASTA.
# Para ello, hay que pasar como argumentos: (1) ficheros de entrada que se quieren
# analizar y (2) nombre del directorio en el que guardar los resultados.

# Numero de argumentos que requiere el script
n=2

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
	echo "Uso: $0 <archivo_input.txt>, <outdir>"
        exit 1
else
        # Ficheros de entrada
	input_files=$1
	
	# Directorio de salida
	outdir=$2
fi

# Bucle que recorre los ficheros y utiliza VirSorter para analizarlos
for file in $input_files;
do
	name=$(basename $file _filtered.fasta)

	echo "Procesando la muestra " $name

	virsorter run -i $file -w $name".out" --include-groups "dsDNAphage" -j 4 all

done	
