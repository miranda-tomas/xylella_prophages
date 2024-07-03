#!/bin/bash

# nanoplot.sh
# Miranda Tomás

# Este script utiliza el programa NanoPlot para obtener el control de calida de ficheros fastq
# procedentes de secuenciación ONT.
# Para ello, hay que pasar como argumentos: (1) los ficheros fastq que se quieren analizar y 
# (2) el nombre del directorio en el que se quieren guardar los resultados.

# Numero de argumentos que requiere el script
n=2

# Verificar el número correcto de argumentos
if [ "$#" -ne "$n" ]; then
    echo "Error: Número incorrecto de argumentos."
    echo "Uso: $0 <archivo_input.txt> <directorio_salida>"
    exit 1

else
	# Fichero con los nombres de los ficheros fastq input
	input_files=$1

	# Directorio en el que se quieren guardar los resultados
	outdir=$2
	mkdir -p $outdir
fi

# Variable para almacenar el estado de error
error_occurred=false

# Funcion para verificar errores durante el script
check_error() {
        if [$? -ne 0 ]; then
                error_occured=true
        fi
}

# Bucle que recorre todos los ficheros de entrada y los analiza con NanoPlot
for file in $input_files; 
do
	name=$(basename $file .fastq)
	
	NanoPlot --fastq $file --loglength -o $outdir --title $name
done

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then
	echo "El programa ha terminado. Los ficheros resultantes se encuentran en el directorio $outdir"
fi
else
	echo "El programa ha finalizado con errores."
fi
