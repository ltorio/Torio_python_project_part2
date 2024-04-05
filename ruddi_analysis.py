from Bio import SeqIO
import os

# ruddi_analysis.py
# L Torio
# downloads specific genbank files and creates a csv of their
# accession id, family, genus, and species, the number of features in the genbank file,
# and the source of the record

# sets working directory to the the same location as this file
# from https://stackoverflow.com/a/1432949
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

f = open("ruddi.csv", "w")


# getting the sequence from the fna file
for rec in SeqIO.parse("GCA_000287275.1_ASM28727v1_genomic.fna", "fasta"):
    sequence = rec.seq  # the sequence
    rev = rec.seq.reverse_complement()  # the reverse compliment sequence

f.write("Length_of_genome,{}\n".format(len(sequence)))
print("Length_of_genome,{}\n".format(len(sequence)))

gc_content = (sequence.count("C") + sequence.count("G")) / len(sequence)
f.write("GC_content,{:.4f}\n".format(gc_content))
print("GC_content,{:.4f}\n".format(gc_content))

f.write("ATG_forward,{}\n".format(sequence.count("ATG")))
print("ATG_forward,{}\n".format(sequence.count("ATG")))

f.write("ATG_reverse,{}\n".format(rev.count("ATG")))
print("ATG_reverse,{}\n".format(rev.count("ATG")))

f.close()
# downloading genbank files for:
# NZ_CALPCP010000001.1 NZ_CALPCY010000130.1 NZ_BHVZ01000001.1 NZ_SRYA01000017.1 NZ_CAJTFZ010000019.1
