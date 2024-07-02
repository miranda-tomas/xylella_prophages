#!/bin/bash

# spades.sh
# Miranda Tomás

# Este script utiliza el programa SPADES para ensamblar reads de ficheros fastq resultantes
# de secuenciacion Illumina.
# Para ello, hay que proporcionar como argumentos: (1) el fichero en el que están los nombres 
# de las muestras que se quieren ensamblar, (2) el directorio en el que se encuentran los ficheros 
# fastq y (3) el directorio en el que se quieren guardar los resultados

# Numero de argumentos que requiere el script
n=3

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt>, <directorio_input>, <directorio_output>"
        exit 1
else
        # Fichero con los nombres de las muestras
        file_inputs=$1
        
        # Directorio donde se encuentran los ficheros fastq de Illumina 
	indir_fastq=$2

	# Directorio donde se quieren guardar los resultados
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

# Bucle que recorre los ficheros, leyendo el fihcero txt con los nombres de las muestras
while IFS= read -r name; 
do
	file1=$indir_fastq/$name"_R1_mapped.fastq"
	file2=$indir_fastq/$name"_R2_mapped.fastq"

	spades.py --only-assembler -1 $file1 -2 $file2 -o $outdir/$name

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then

        echo "El programa ha finalido. Los ficheros resultantes se encuentran en el directorio: $outdir"

else
        echo "El programa ha finalizado con errores."
fi
