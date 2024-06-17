#!/bin/bash
#SBATCH -n 10           # Number of cores requested
#SBATCH -N 1            # Number of nodes requested
#SBATCH -t 24:00:00     # Runtime in minutes
#SBATCH --qos medium    # The QoS to submit the job
#SBATCH --mem=100G     	# Memory per cpu in G 
#SBATCH -o blastn.out  	# Standard output goes to this file
#SBATCH -e blastn.err  	# Standard error goes to this file

# blastn.sh
# Miranda Tomas

# Este script compara mediante Blastn las secuencias input frente a la base de datos de virus o procaritas, según elija el usuario.
# Para ello, hay que pasar como argumentos: (1) fichero txt con los nombres de los ficheros input, (2) base de datos que se quiere
# utilizar (procariotas o virus) y (3) nombre del directorio en el que guardar los ficheros resultantes.

# Numero de argumentos que requiere el script
n=3
     
# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt> <base_de_datos_(virus_o_prok)> <directorio_salida>"
        exit 1
else
        # Fichero con los nombres de las muestras
        file_inputs=$1

	# Base de datos que se quiere utilizar (viruses o prok)
	db=$2
	database=/home/principal/blast_db_$db/nt_$db

        # Nombre del directorio en el que se quieren guardar los fastq filtrados
	outdir=$3
	mkdir -p $outdir
        
fi

# Verificar si el fichero con los nombres de los inputs existe
if [ ! -f "$file_inputs" ]; then
        echo "Error: El fichero $file_inputs no existe."
        exit 1
fi

# Variable para almacenar el estado de error
error_occurred=false

# Funcion para verificar errores durante el script
check_error() {
        if [ $? -ne 0 ]; then
                error_occured=true
                echo "Error: $1"
        fi
}

# Valor de corte del e-value que determina la significancia de las coincidencias
evalue=0.005
	
# Formato de los ficheros de salida
formato="6 qseqid saccver sscinames stitle pident length evalue mismatch qstart qend sstart send"

# Numero maximo de coincidencias para cada query
max_target_seqs=5

# Numero maximo de alineamientos para cada pareja query/target
max_hsps=1

# Bucle que recorre los ficheros, leyendo el fihcero txt con los nombres de las muestras
while IFS= read -r name; 
do
        echo "Procesando la muestra $name"
        
        in_file="all_assembly/"$name"_filtered.fasta"
        out_file=$outdir"/"$name"_blast.csv"
        
	blastn -query $in_file -db $database -evalue $evalue -outfmt "$formato" -max_target_seqs $max_target_seqs -max_hsps $max_hsps > $out_file
	echo >> $out_file

done < $file_inputs
