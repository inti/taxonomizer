import pysam
import numpy as np
import pandas as pd 
import subprocess
import h5py
import os

dt_str_vlen = h5py.special_dtype(vlen=str)

######################################################################################


def create_genome_group_and_Q_tables(h5_group,specie,R,G):
        if specie not in h5_group:
                h5_specie_folder = h5_group.create_group(specie)
                h5_specie_folder.create_dataset("Q",(R,G),dtype='float16')
                h5_specie_folder.create_dataset("read_names",(R,2),dtype=dt_str_vlen)
                h5_specie_folder.create_dataset("reference_names",(R,2),dtype=dt_str_vlen )
		
        else:
                print specie, "already exists on the DB"



#######################################################################################


aln_folder = "/media/TeraData/ipedroso/METAGENOMICS/ANALYSES/ALN/simLC/"

proc = subprocess.Popen(["".join(("ls ",aln_folder,"/Z*bam"))], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()

bam_files = out.split()
all_genomes=pd.concat([ pd.DataFrame([pysam.AlignmentFile(bam, "rb").count(),
			pysam.AlignmentFile(bam, "rb").nreferences, 
			os.path.basename(bam).replace(".sorted.bam","")],
		      index=["R","G","species_name"]).T	
		 for bam in bam_files])

h5_file = h5py.File("mytestfile.hdf5", "w")

genomes = h5_file.create_group("genomes")

#all_genomes.apply(lambda (nR,nG,species_name): create_genome_group_and_Q_tables(genomes,species_name,nR,nG),axis=1)

for bam in bam_files:
	samfile = pysam.AlignmentFile(bam, "rb")
	species_name = os.path.basename(bam).replace(".sorted.bam","")
	nR = samfile.count()
	nG = samfile.nreferences
	
	create_genome_group_and_Q_tables(genomes,species_name,nR,nG)

	read_counter = 0
	for read in samfile.fetch():
		genomes[species_name]["Q"][read_counter,read.rname] = read.mapq
		genomes[species_name]["read_names"][read_counter,:] = read.qname, read_counter
		read_counter += 1




