from Bio import SeqIO
from Bio import Entrez

"""
fasta_parse.py
L Torio

downloads specific protein fasta files and creates a csv of their 
fasta id, first 10 amino acids in the protein, the length of the protein, and 
the number of cysteines in the protein
"""

# %% Biopython Exercise 2: Fasta fun with proteins!

# downloading genbank files for:
# AGI40145.1, AGJ87295.1, WVV45440.1, WVS05366.1
with Entrez.efetch(
    db="protein",
    rettype="fasta",
    retmode="text",
    id="AGI40145.1, AGJ87295.1, WVV45440.1, WVS05366.1",
) as fastaFiles:
    # writing fasta id, first 10 amino acids in the protein,
    # the length of the protein, and the number of cysteines in the protein
    # in a csv for each fasta file

    f = open("protein_info.csv", "w")
    f.write("ID,First_10_AA,Length,Number_Cs\n")  # column titles for csv

    # loops thru each fasta file to write info to protein_info.csv
    for seq_record in SeqIO.parse(fastaFiles, "fasta"):
        print(
            "writing %s to csv" % (seq_record.id)
        )  # printing which file is being worked on
        f.write(
            "%(ID)s,%(First_10_AA)s,%(Length)s,%(Number_Cs)s, \n"
            % {
                "ID": seq_record.id,
                "First_10_AA": seq_record.seq[0:10],
                "Length": len(seq_record.seq),
                "Number_Cs": seq_record.seq.count("C"),
            }
        )
    f.close()
# %% Test code
# import os
# os.chdir("/Users/lyndo/OneDrive/Desktop/2024Spring/BIOL668/Python/pythonProjectPart2/")
# subsub = SeqIO.read("AGI40145.fasta","fasta")
# print(subsub.id)
# print(subsub.seq[0:10])
# print(len(subsub.seq))
# print(subsub.seq.count("C"))
