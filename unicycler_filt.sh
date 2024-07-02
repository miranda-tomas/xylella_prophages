#!/bin/bash

# filtrar_contigs.sh
# Miranda Tomas

# Este script filtra los contigs de un fasta resultante de un ensamblaje por Unicycler y 
# selecciona aquellos que tengan una longitud y un coverage mínimo, indicado por el usuario. 
# Para ello, hay que proporcionar como argumentos: (1) fichero txt con los nombres de las muestras,
# (2) longitud y (3) cobertura mínimos a las que se quiere filtrar. 

# Numero de argumentos que requiere el script
n=3

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
        echo "Uso: $0 <archivo_input.txt> <longitud_minima> <cobertura_minima>"
        exit 1
else
        # Fichero con los nombres de las muestras
        input_files=$1
        
        # Parámetros por los que se quiere filtrar (longitud y coverage mínimos)
	len_min=$2
	cov_min=$3
        
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


out_file_txt="contigs_filtrados.txt"

for file in $input_files;
do
	name=$(basename $file .fasta)
	file2=$name"_lineas.fasta"
	out_file=$name"_filtered.fasta"

	# Eliminar el salto de linea entre una misma secuencia
        awk '$0~/>/ {print "\n" $0} $0!~ />/ {printf $0}' $file > $file2
        sed -i "1d" $file2

        siguiente_linea=""

        while IFS= read -r linea
        do
                if [[ -n "$siguiente_linea" ]]; then
                        echo "$siguiente_linea" >> $out_file
                        siguiente_linea=""
                fi

                if [[ $linea == ">"* ]]; then

                        len=$(echo "$linea" | awk -F 'length=' '{print $2}' | awk '{print $1}')
                        cov=$(echo "$linea" | awk -F 'depth=|x' '{print $2}' | awk '{print $1}')
			
                        if (( $(echo "$len >= $len_min") )) && (( $(echo "$cov >= $cov_min" | bc -l) )); then

                                echo $linea >> $out_file
                                siguiente_linea=$(IFS= read -r; echo "$REPLY")

                        fi
                fi

        done < "$file2"

	if [[ -n "$siguiente_linea" ]]; then
		echo "$siguiente_linea" >> $out_file
	fi
	
	rm $file2
	
	echo $out_file >> $out_file_txt
	
	num_contigs_ini=$(grep ">" $file | wc -l)
	num_contigs_filt=$(grep ">" $out_file | wc -l)
	echo "Contigs iniciales: " $num_contigs_ini >> $out_file_txt
	echo "Contigs filtrados: " $num_contigs_filt >> $out_file_txt

	grep ">" $out_file | while IFS= read -r line; do
        	contig=$(echo $line | awk '{split($2, a, "="); print substr($1, 2)}')
                length=$(echo "$line" | grep -oP '(?<=length=)[0-9]+')
                cov=$(echo "$line" | grep -oP '(?<=depth=)[0-9]+\.[0-9]+')

                echo -e "${contig}\t${length}\t${cov}"
        
	done >> $out_file_txt
	echo >> $out_file_txt
done

# Verificar si ha ocurrido algún error
if [ "$error_occurred" = false ]; then
        echo "El programa ha finalido. Los ficheros resultantes se encuentran en el directorio: $outdir"
        echo "Se ha creado el fichero $out_file_txt con un resumen de los contigs filtrados.

else
        echo "El programa ha finalizado con errores."
fi

	
