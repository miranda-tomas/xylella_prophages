#!/bin/bash

# mapeo_virales.sh
# Miranda Tomás

# Este script utiliza el programa bwa mem para alinear ficheros fastq frente a un fichero fasta que 
# contiene los contigs clasificados como virales, samtools para obtener solo las reads virales 
# y bedtools para guardarlas en nuevos ficheros fastq.

# Para ello, hay que pasar como argumentos: (1) fichero txt con los nombres de los ficheros fastq que se
# quieren procesar, (2) nombre del directorio con los ficheros fastq que se quieren alinear, (3) nombre 
# del directorio con los ficheros fasta que se utilizan como referencia para alinear y (4) nombre del 
# directorio en el que guardar los ficheros fastq con las lecturas que mapean

# Numero de argumentos que requiere el script
n=4

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt>, <directorio_salida>"
        exit 1
else
	# Fichero con los nombres de las muestras
	file_inputs=$1

	# Nombre del directorio con los ficheros fastq de Illumina
	indir_illumina=$2

	# NOmbre del directorio con el fichero fasta con las secuencias virales
	indir_fasta=$3

	# Nombre del directorio en el que se quieren guardar los resultados
	outdir=$4
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

# Bucle que recorre los ficheros, leyendo el fichero txt con los nombres de las muestras
while IFS= read -r name;
do
	echo "Analizando muestra " $name
	
	# Ficheros fastq con todas las reads 
	file1=$indir_illumina/$name"_R1_filtered.fastq"
	file2=$indir_illumina/$name"_R2_filtered.fastq"

	# Fichero de referencia
	reference=$indir_fasta/$name"_viral.fasta"
	
	dir_index="index"
	mkdir -p $dir_index
		
	# Indexar genoma referencia
	if [ ! -f $dir_index/$reference".sa" ]; then 
		bwa index $reference 
		mv *fasta.* $dir_index
	fi
		
	# Mapear fastq contra el fago que se usa de referencia
	bwa mem $dir_index/$reference $file1 $file2 > $outdir/$name".sam"
		
	# Convertir sam a bam 
	samtools view -S -b $outdir/$name".sam" > $outdir/$name".bam"
	rm $outdir/$name".sam"
		
	# Obtener las reads que han mapeado
	samtools view -b -F 4 $outdir/$name".bam" > $outdir/$name"_mapped.bam"
	rm $outdir/$name".bam"

	# Pasar las reads mapeadas a nuevos fastq
	bedtools bamtofastq -i $outdir/$name"_mapped.bam" \
				-fq $outdir/$name"_R1_mapped.fastq" \
				-fq2 $outdir/$name"_R2_mapped.fastq"
done < $file_inputs
