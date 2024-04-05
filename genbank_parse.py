from Bio import SeqIO
from Bio import Entrez

"""
genbank_parse.py
L Torio

downloads specific genbank files and creates a csv of their 
accession id, family, genus, and species, the number of features in the genbank file, 
and the source of the record
"""
# %% Biopython Exercise 1: Parse Genbank files!
Entrez.email = "ltorio7987@sdsu.edu"

# downloading genbank files for:
# NZ_CALPCP010000001.1 NZ_CALPCY010000130.1 NZ_BHVZ01000001.1 NZ_SRYA01000017.1 NZ_CAJTFZ010000019.1
with Entrez.efetch(
    db="nucleotide",
    rettype="gb",
    retmode="text",
    id="NZ_CALPCP010000001.1,NZ_CALPCY010000130.1,NZ_BHVZ01000001.1,NZ_SRYA01000017.1,NZ_CAJTFZ010000019.1",
) as gbFiles:
    # writing family, genus, and species, the number of features in the genbank , and the source of the record
    # in a csv for each genbank file
    f = open("genbank_parse.csv", "w")
    f.write(
        "Accession,Family,Genus,Species,Num_Features,Source\n"
    )  # column titles for csv

    # loops thru each genbank file to write to genbank_parse.csv
    for seq_record in SeqIO.parse(gbFiles, "gb"):
        print(
            "writing %s to csv" % (seq_record.id)
        )  # printing which file is being worked on
        f.write(
            "%(accession)s,%(family)s,%(genus)s,%(species)s,%(num_Features)s,%(source)s, \n"
            % {
                "accession": seq_record.id,
                "family": seq_record.annotations["taxonomy"][-2],
                "genus": seq_record.annotations["taxonomy"][-1],
                "species": seq_record.annotations["organism"],
                "num_Features": len(seq_record.features),
                "source": seq_record.annotations["source"],
            }
        )
    f.close()
