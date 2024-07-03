#!/bin/bash

# pilon.sh
# Miranda Tomas

# Este script corrige genomas en formato FASTA. Primero utiliza Samtools para alinear las lecturas,
# y ordenar e indexar el fichero BAM resultante. Este lo utiliza Pilon para detectar errores de
# ensamblaje y los soluciona, creando un nuevo fasta con el genoma corregido.

# Para ello, hay que pasar como argumentos: (1) fichero txt con el nombre de los ficheros fasta
# que se quieren corregir, (2) directorio que contiene los ficheros fasta, (3) directorio que 
# contiene los ficheros fastq correspondientes y (4) nombre del directorio en el que guardar los 
# ficheros fasta corregidos.

# Numero de argumentos que requiere el script
n=4

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
	echo "Uso: $0 <archivo_input.txt> <directorio_entrada_fasta> <directorio_entrada_fastq> <directorio_salida>"
        exit 1
else
      
	# Fichero con el nombre de las muestras
	file_inputs=$1

	# Directorio con los ficheros fasta que se quieren corregir
	$indir_fasta=$2

	# Direcotiro con los ficheros fastq necesarios
	$indir_fastq=$3

	# Directorio donde se guardan ficheros de salida
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
	echo "Procesando la muestra $name"

	# Fichero con el contig a corregir
	file=$indir_fasta"/"$name"_prophage.fasta"

	# Ficheros fastq con las reads con las que se ha construido ese contig
	fastq1=$indir_fastq"/"$name"_R1_filtered.fastq"
	fastq2=$indir_fastq"/"$name"_R2_filtered.fastq"

	# Nombre del fichero con el contig tras la correcion 
	output_file=$name"_prophage_pilon"

	#Extraer coverage por contig con BBMap
	~/bbmap/bbmap.sh ref=$file in1=$fastq1 in2=$fastq2 covstats=$outdir/$name"_covstats.txt" out=$outdir/$name"_mapped.sam"

		#covstats contiene info sobre el coverage de los contigs ("avg_fold")
		#contig_mapped.sam contiene los reads alineados con el contig del input file

	# Crear el fichero BAM con las reads alineadas
	samtools view -bS -F4 $outdir/$name"_mapped.sam" | samtools sort - -o $outdir/$name"_mapped_sorted.bam"
	samtools index $outdir/$name"_mapped_sorted.bam"
	
	# Correcion del ensamblaje
	pilon --genome $file --frags $outdir/$name"_mapped_sorted.bam" --output $outdir/$output_file --verbose --changes

done < $file_inputs

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then	
	
	# Eliminar ficheros temporales
	rm $outdir/*nophi.fastq
	
	echo "El programa ha finalido. Los ficheros resultantes se encuentran en el directorio: $outdir"
	echo "Los detalles de los cambios se muestra en .changes"
	echo "Volver a repetir hasta que no haya cambios (fichero .changes vacio)"

else
	echo "El programa ha finalizado con errores."
fi

