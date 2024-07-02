#!/bin/bash
#SBATCH -n 10             # Number of cores requested
#SBATCH -N 1              # Number of nodes requested
#SBATCH -t 48:00:00       # Runtime in minutes
#SBATCH --qos medium      # The QoS to submit the job
#SBATCH --mem=30G         # Memory per cpu in G 
#SBATCH -o unicycler.out  # Standard output goes to this file
#SBATCH -e unicycler.err  # Standard error goes to this file

# unicycler.sh
# Miranda Tomas

# Este script utiliza el programa UNICYCLER para ensamblar reads de ficheros fastq resultantes 
# de secuenciacion con Illumina y con Nanopore que han sido previamente procesados. 
# Para ello, hay que proporcionar como argumentos: (1) el fichero en el que están los nombres 
# de las muestras que se quieren ensamblar, el directorio en el que se encuentran los
# ficheros filtrados de Illumina (2) y de ONT (3)

# Ruta al programa 
path="/storage/enbivir/software/Unicycler"
spades_path="/storage/enbivir/software/SPAdes-3.15.4-Linux/bin/"

# Numero de argumentos que requiere el script
n=3

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt>, <directorio_input_illumina>, <directorio_input_ONT>"
        exit 1
else
        # Fichero con los nombres de las muestras
        file_inputs=$1
        
        # Directorio donde se encuentran los ficheros fastq de Illumina (los filtrados con bbmap)
	indir_illumina=$2

	# Directorio donde se encuentra el fichero fastq de Nanopore (el filtrado con nanofilt)
	indir_ONT=$3
        
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

# Bucle que recorre los ficheros, leyendo el fihcero txt con los nombres de las muestras
while IFS= read -r name; 
do
        echo "Procesando la muestra $name"

        # Nombre de los ficheros
        read1=$indir_illumina"/"$name"_R1_filtered.fastq"
        read2=$indir_illumina"/"$name"_R2_filtered.fastq"
        nanopore=$indir_ONT"/"$name"_ONT_filtered.fastq"
        
        # Nombre del directorio en el que se guardan los resultados de ensamblaje
	outdir="res_unicycler_"$name
        mkdir -p $outdir

	$path/unicycler-runner.py --threads 10 --spades_path $spades_path/spades.py -1 $read1 -2 $read2 -l $nanopore -o $outdir

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then

        echo "El programa ha finalido. Los ficheros resultantes se encuentran en el directorio: $outdir"

else
        echo "El programa ha finalizado con errores."
fi

