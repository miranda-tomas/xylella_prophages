#!/bin/bash
#SBATCH -n 1		# Number of cores requested
#SBATCH -N 1		# Number of nodes requested
#SBATCH -t 2:00:00	# Runtime in minutes
#SBATCH --qos short	# The QoS to submit the job
#SBATCH --mem=10G	# Memory per cpu in G 
#SBATCH -o fastqc.out	# Standard output goes to this file
#SBATCH -e fastqc.err	# Standard error goes to this file

# fastqc.sh
# Miranda Tomas 

# Este script utiliza el programa fastqc para obtener el control de calidad de ficheros fastq paired-ends 
# obtenidos por secuenciación Illumina. 
# Para ello, hay que pasar como argumentos: (1) fichero txt con los nombres de los ficheros fastq que se 
# quieren procesar y (2) el nombre del directorio en el que guardar los resultados.

# Ruta para el programa
path=/storage/enbivir/software/FastQC

# Verificar el número correcto de argumentos
if [ "$#" -ne 2 ]; then
    echo "Error: Número incorrecto de argumentos."
    echo "Uso: $0 <archivo_input.txt> <directorio_salida>"
    exit 1

else
	# Fichero con los nombres de los ficheros fastq input
	file_inputs=$1

	# Directorio en el que se quieren guardar los resultados
	outdir=$2
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
	echo "Procesado muestra $name"

        # Nombre de los ficheros 
	file1=$name"_R1.fastq"
	file2=$name"_R2.fastq"

        # Fastqc
        $path/fastqc -o $outdir $file1
	check_error "Fallo en ejecutar FastQC para $file1"
        $path/fastqc -o $outdir $file2
	check_error "Fallo en ejecutar FastQC para $file2"

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then
	echo "El programa ha terminado. Los ficheros resultantes se encuentran en el directorio $outdir"
fi
else
	echo "El programa ha finalizado con errores."
fi
