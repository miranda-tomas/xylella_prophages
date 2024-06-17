# -*- coding: utf-8 -*-
#*****************************************************************************
# 
# genome_reordering.py
#
# Este programa reordena un conjunto de ficheros fasta tomando como punto de
# inicio una secuencia de referencia, obtenienod las coordenadas a través de 
# un alineamiento Blastn de ambas
# 
# Para ello se pasa como referencia los ficheros fasta input, el fichero fasta
# con la secuencia de referencia y el nombre del fichero con el output del 
# alineamiento Blast 
#
# 21/01/2024
# Miranda Tomas
#
#*****************************************************************************


from Bio import SeqIO
import sys, os, subprocess

def Blast(query_file:str, subject_file:str, blast_file:str):
    """
    Esta funcion utiliza la linea de comandos de bash para realizar un alineamiento con blastn
    entre dos ficheros fasta, guardando el resultado en un fichero csv, todos pasados como
    parámetros

    Parameters
    ----------
    query_file : str
        Nombre del fichero que se utiliza como query en el alineamiento por Blast

    subject_file : str
        Nombre del fichero que se utiliza como subject en el alineamiento por Blast

    blast_file : str
        Nombre del fichero en el que se guardan los resultados de los alineamientos por Blast
    """

    comando_blast = f"blastn -query {query_file} -subject {subject_file} -max_target_seqs 1 -max_hsps 1 -outfmt 6 >> {blast_file} "

    process = subprocess.Popen(comando_blast, shell=True)
    process.wait()


def Reverse(in_file:str, out_file:str):
    """
    Esta funcion hace el reverso de una secuencia fasta y la guarda en un nuevo fichero fasta 
    
    Parameters
    ----------
    in_file : str
        Nombre del fichero con 
    
    out_file : str
        Nombre del fichero con 

    """
    
    for record in SeqIO.parse(open(in_file), 'fasta'):
        reverse_seq = str(record.seq.reverse_complement())
    with open(out_file, "w") as out:
        out.write('>' + record.id + '_reversed' + '\n' + reverse_seq)


def Reorder(in_file:str, out_file:str, pos:int):
    """
    Esta funcion reordena una secuencia fasta en funcion a unas coordenadas dadas y la guarda 
    en un nuevo fichero fasta 
    
    Parameters
    ----------
    in_file : str
        Nombre del fichero con 
    
    out_file : str
        Nombre del fichero con 

    pos : int
        Entero con 

    """
    
    for contig_record in SeqIO.parse(open(in_file), 'fasta'):
        contig = str(contig_record.seq)
        str1 = contig[pos - 1:]
        str2 = contig[:pos - 1]
        str_final = '>' + contig_record.id + '\n' + str1 + str2
    with open(out_file, "w") as out:
        out.write(str_final)

            
def main():

    # Leer inputs
    in_files = sys.argv[1]
    referencia = sys.argv[2]
    blast_file = sys.argv[3]

    # Bucle que recorre los ficheros de entrada 
    for in_file in in_files.split():

        basename = os.path.splitext(in_file)[0]

        # Blast de la secuencia de referencia dada frente al genoma del fago
        Blast(referencia, in_file, blast_file)

    
    # Leer el fichero con los resultados del Blast
    with open(blast_file, "r") as in_file:

        for line in in_file:
        
            # Obtener nombre del fichero, posicion de inicio y fin del alineamiento 
            file = line.split()[1]
            ini = int(line.split()[8])
            fin = int(line.split()[9])

            # Si la posicion de inicio es mayor que la de fin, se hace el reverso
            if ini > fin:
            
                in_file_reverse = file+".fasta"
                out_file_reverse = file+"_reversed.fasta"

                Reverse(in_file_reverse, out_file_reverse)

                # Y se vuelven a sacar las coordenadas del alineamiento con Blast
                Blast(referencia, out_file_reverse, blast_file)
                

            # Si la posicion de inicio es menor que la de fin, se reordena el genoma  
            if ini < fin:
                
                in_file_reorder = file+".fasta"
                out_file_reorder = file+"_reordered.fasta"

                Reorder(in_file_reorder, out_file_reorder, ini)


if __name__ == '__main__':
    main()

    
