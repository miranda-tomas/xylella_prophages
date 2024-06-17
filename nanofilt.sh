#!/bin/bash
#SBATCH -n 1             # Number of cores requested
#SBATCH -N 1             # Number of nodes requested
#SBATCH -t 2:00:00       # Runtime in minutes
#SBATCH --qos short      # The QoS to submit the job
#SBATCH --mem=10G        # Memory per cpu in G 
#SBATCH -o nanofilt.out  # Standard output goes to this file
#SBATCH -e nanofilt.err  # Standard error goes to this file

# nanofilt.sh
# Miranda Tomas

# Este script utiliza el programa NanoFilt para filtrar los ficheros fastq obtenidos por secuenciación ONT.
# Para ello, hay que pasar como argumentos: (1) los ficheros de entrada, (2) la calidad y (3) lo longitud
# mínima a la que se quiere filtrar y (4) nombre del directorio en el que guardar los fastq filtrados.

# Verificar el número correcto de argumentos
if [ "$#" -ne 2 ]; then
    echo "Error: Número incorrecto de argumentos."
    echo "Uso: $0 <archivo_input.txt> <directorio_salida>"
    exit 1

else
	# Fichero con los nombres de los ficheros fastq input
	file_inputs=$1
	
	# Calidad mínima a la que se quiere filtrar
	quality=$2
	
	# Longitud mínima a la que se quiere filtrar
	length=$3

	# Directorio en el que se quieren guardar los resultados
	outdir=$4
	mkdir -p $outdir
fi

# Verificar si el fichero con los nombres de los inputs existe
if [ ! -f "$1" ]; then
        echo "El fichero $1 no existe"
        exit 1
else
        file_inputs=$1
fi

# Variable para almacenar el estado de error
error_occurred=false

# Funcion para verificar errores durante el script
check_error() {
        if [$? -ne 0 ]; then
                error_occured=true
        fi
}

# Bucle que recorre los ficheros, leyendo el fihcero txt con los nombres de las muestras
while IFS= read -r name;
do
        input_file=$name"_ONT.fastq"
        output_file=$name"_ONT_filtered.fastq"

        NanoFilt -q $quality -l $length $input_file > $outdir/$output_file
        check_error "Fallo en ejecutar NanoFilt para $input_file"

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then
	echo "El programa ha terminado. Los ficheros resultantes se encuentran en el directorio $outdir"
else
	echo "El programa ha finalizado con errores."
fi



