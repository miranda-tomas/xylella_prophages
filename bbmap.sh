#!/bin/bash
#SBATCH -n 1            # Number of cores requested
#SBATCH -N 1            # Number of nodes requested
#SBATCH -t 12:00:00     # Runtime in minutes
#SBATCH --qos medium    # The QoS to submit the job
#SBATCH --mem=10G    	# Memory per cpu in G 
#SBATCH -o bbmap.out  	# Standard output goes to this file
#SBATCH -e bbmap.err  	# Standard error goes to this file

# bbmap.sh
# Miranda Tomas

# Este script utiliza el programa bbmap para limpiar los ficheros fastq paired-ends obtenidos por secuenciación Illumina. En primer lugar, se eliminan las reads correspondientes al fago phiX y después se elimina el adaptador y se filtran las reads que tengan una calidad mayor a la indicada.
# Para ello, hay que pasar como argumentos: (1) fichero txt con los nombres de los ficheros fastq que se quieren procesar, (2) la posicion del adaptador, (3) la longitud y (4) la calidad minima que se quieren en las reads resultantes, (5) lado de las secuencias por el que se quiere filtrar y (6) nombre del directorio en el que guardar los fastq filtrados.

# Ruta a la carpeta del programa 
path=/storage/enbivir/software/BBMap_38.95/bbmap

# Numero de argumentos que requiere el script
n=6

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
	echo "Uso: $0 <archivo_input.txt> <posicion_adaptador> <longitud_mínima> <calidad_mínima> <lado_filtrado> <directorio_salida>"
        exit 1
else
        # Fichero con los nombres de las muestras
	file_inputs=$1
	
	# Posicion del adaptador: izq (l), dcha (r) o no hay (f)
	adapter_place=$2

	# Longitud minima de las reads resultantes
        min_len=$3

	# Calidad minima a la que se quiere filtrar
        quality=$4
	
	# Indicar si se quiere filtrar por la inzquierda (l), derecha (r), ambos (rl) o por ninguno (f)	
	side=$5
        
	# Nombre del directorio en el que se quieren guardar los fastq filtrados
	outdir=$6
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
	echo "Procesando la muestra $name"

        # Nombre de los ficheros
        file1=$name"_R1.fastq"
        file2=$name"_R2.fastq"

        # Limpiar si hay reads de fago phiX
        $path/bbduk.sh in1=$file1 in2=$file2 out1=$outdir/$name"_R1_nophi.fastq" out2=$outdir/$name"_R2_nophi.fastq" ref=$path/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=$outdir/phiX_stats.txt
	check_error "Fallo en limpiar reads de phiX"

        # Eliminar el adaptador y filtrar las reads
        $path/bbduk.sh in1=$outdir/$name"_R1_nophi.fastq" in2=$outdir/$name"_R2_nophi.fastq" out1=$outdir/$name"R1_filtered.fastq" out2=$outdir/$name"R2_filtered.fastq" ref=$path/resources/adapters.fa ktrim=$adapter_place hdist=1 tpe tbo minlen=$min_len trimq=$quality qtrim=$side
	check_error "Fallo en eliminar el adaptar y filtrar las reads"

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then	
	
	# Eliminar ficheros temporales
	rm $outdir/*nophi.fastq
	
	echo "El programa ha finalido. Los ficheros resultantes (filtered) se encuentran en el directorio: $outdir"

else
	echo "El programa ha finalizado con errores."
fi
