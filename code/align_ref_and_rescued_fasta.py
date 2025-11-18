import sys

ref_seq = sys.argv[1]
snp_seq = sys.argv[2] # .fasta individual sin header obtenido del multifasta a partir del consensus (se puede dividir el multifasta con el script zsplit_multifasta.sh)
snp_table = sys.argv[3]
sample_name = snp_seq.split(".")[0]
print(sample_name)
positions_snptable = []

# Importamos la SNP_table y nos quedamos con las posiciones
with open(snp_table, "r+") as fichero:
    lines = fichero.readlines()
    lines = lines[1:]

    for line in lines:
        positions_snptable.append(int(line.split("\t")[0].replace(" ", "")))

with open(ref_seq, "r+") as fichero:
    lineas = fichero.readlines()
    lineas = lineas[1:]
    lineas_strip = [line.strip() for line in lineas]
    ref_seq = "".join(lineas_strip)

    newseq_rescatada = list(ref_seq) # hacemos lista la secuencia de referencia para reemplazar los nucleótidos

with open(snp_seq, "r+") as fichero:
    lineas = fichero.readlines()
    linea = lineas[0].strip()
    snp_seq = list(linea) # hacemos lista la secuencia de SNPs

count = 0

for position in positions_snptable:
    #print("Position: ",position)
    pos_in_ref = position-1
    #print(newseq_rescatada[pos_in_ref])
    #print(snp_seq[count])
    newseq_rescatada[pos_in_ref] = snp_seq[count]
   # print(count, snp_seq[count])
    count += 1


with open("{}_aligned.fasta".format(sample_name), "w+") as output_file:
    output_file.write(">MTB_anc\n{0}\n>{1}\n{2}\n".format(ref_seq,sample_name,"".join(newseq_rescatada)))
