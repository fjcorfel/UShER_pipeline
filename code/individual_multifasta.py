'''Este script toma un archivo .snp y genera un alineamiento de todo el genoma de nuestra referencia contra el genoma de dicha muestra
teniendo en cuenta los SNPs, de manera que podamos sacar un VCF para sacar el archivo diff'''

import sys

seq_file = sys.argv[1] # Referencia .fasta
snp_file = sys.argv[2] # Archivo annoF
sample_name = snp_file.split(".")[0]
ref_seq = ""

# Importamos la secuencia y la convertimos en lista para poder modificarla
with open(seq_file, "r+") as fichero:
    lines = fichero.readlines()
    lines = lines[1:]
    for line in lines:
        ref_seq += line.rstrip()
    new_seq = list(ref_seq) # Secuencia que sera modificada

# Importamos el archivo de SNP y, para cada posicion, comprobamos que la referencia sea la correcta y en tal caso modificamos con la variante
with open(snp_file, "r+") as fichero:
    lines = fichero.readlines()
    lines = lines[1:]

    for line in lines:
        tokens = line.split("\t")
        position,ref,alt = tokens[1:4]

        if new_seq[int(position)-1] == ref:
            new_seq[int(position)-1] = alt
        else:
            print("La posicion {0} no coincide".format(position))

# Guardamos el ouput
with open("{0}_and_ref_multifasta.fas".format(sample_name), "w+") as output:
    output.write(">MTB_ancestor\n")
    output.write(ref_seq+"\n")
    output.write(">{0}\n".format(sample_name))
    output.write(("").join(new_seq))
