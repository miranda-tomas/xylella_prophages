#!/bin/bash
#SBATCH -n 1            # Number of cores requested
#SBATCH -N 1            # Number of nodes requested
#SBATCH -t 2:00:00    	# Runtime in minutes
#SBATCH --qos medium   	# The QoS to submit the job
#SBATCH --mem=10G      	# Memory per cpu in G
#SBATCH -o checkv.out  	# Standard output goes to this file
#SBATCH -e checkv.err  	# Standard error goes to this file

# checkv.sh
# Miranda Tomas

# Este script utiliza el programa CheckV para evaluar la calidad de genomas virales y 
# la estimación de la integridad de los fragmentos del genoma. 

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


# Bucle que recorre los ficheros y usa la pipeline de CheckV
for file in $input_files; do

        outdi2r=$outdir/$(basename $file .fasta)

        checkv contamination $file $outdir2 -t 16
        checkv completeness $file $outdir2 -t 16
        checkv complete_genomes $file $outdir2
        checkv quality_summary $file $outdir2

done

