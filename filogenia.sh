#SBATCH -n 1        		# Number of cores requested
#SBATCH -N 1           		# Number of nodes requested
#SBATCH -t 24:00:00  		# Runtime in minutes
#SBATCH --qos medium     	# The QoS to submit the job
#SBATCH --mem=10G       	# Memory per cpu in G 
#SBATCH -o filogenia.out	# Standard output goes to this file
#SBATCH -e filogenia.err	# Standard error goes to this file

# filogenia.sh
# Miranda Tomás

# Este script realiza un alinemiento de secuencias fasta incluidas en un fichero FASTA utilizando
# MAFFT y construye una filogenia utilizando iqtree, calculando el mejor modelo nucleotídico para ello.

# Numero de argumentos que requiere el script
n=4

# Verificar si se proporciona un numero de argumentos correcto
if [ "$#" -ne "$n" ]; then
        echo "Error: Numero incorrecto de argumentos."
	echo "Uso: $0 <fichero_input> <fichero_alineamiento> <fichero_arbol> <numero_bootstrap>"
	exit 1
else
        # Fichero de entrada con las secuencias problema
	input_file=$1
	
	# Fichero de salida donde guardar el alineamiento
	alignment_file=$2

	# Fichero de salida donde guardar la filogenia
        tree_file=$3

	# Numero de bootstrap utilizado
        bootstrap=$4	

fi

# Alinear secuencias con MAFFT
mafft --auto $input_file > $alignment_file

# Determinar el mejor modelo de sustitución con IQ-TREE
iqtree2 -s $alignment_file -m TEST -nt AUTO

# Ejecutar IQ-TREE con el mejor modelo encontrado
best_model=$(grep "Best-fit model according to BIC:" $alignment_file.iqtree | awk '{print $NF})'
iqtree2 -s $alignment_file -m $best_model -bb $bootstrap -nt AUTO -pre $tree_file
