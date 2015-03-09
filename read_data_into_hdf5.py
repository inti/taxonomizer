from time import ctime
import pysam
import numpy as np
import pandas as pd 
import subprocess
import h5py
import os
from progressbar import ProgressBar, Percentage, Bar, ETA, FileTransferSpeed, RotatingMarker

widgets = ['Proggres: ', Percentage(), ' ', Bar(), ' ', ETA(), ' ', FileTransferSpeed()] 

dt_str_vlen = h5py.special_dtype(vlen=str)

dtype_stringVlen_float16 =[('name1', dt_str_vlen), ('name2', 'float16')]


######################################################################################
def openBamfile(bam):
	samfile = pysam.AlignmentFile(bam, "rb")
        species_name = os.path.basename(bam).replace(".sorted.bam","")
        nR = samfile.count()
        nG = samfile.nreferences
	return samfile, species_name, nR, nG

def rescale_samMAPQ(mapq):
	return 1 - np.power(10,mapq/-10)

def create_genome_group_and_Q_tables(h5_group,specie,R,G,ref_names):
        if specie not in h5_group:
                h5_specie_folder = h5_group.create_group(specie)
                h5_specie_folder.create_dataset("Q",(R,G),dtype='float')
                h5_specie_folder.create_dataset("read_names",(R,1),dtype=dt_str_vlen) 
                h5_specie_folder.create_dataset("reference_names",(G,),dtype=dt_str_vlen)
		h5_specie_folder["reference_names"][:] = [ rname for i,rname in enumerate(ref_names)]	
        else:
                print specie, "already exists on the DB"



#######################################################################################


aln_folder = "/media/TeraData/ipedroso/METAGENOMICS/ANALYSES/ALN/simLC/"

proc = subprocess.Popen(["".join(("ls ",aln_folder,"/Esch*bam"))], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
bam_files = out.split()

h5_file = h5py.File("mytestfile.hdf5", "w")

genomes = h5_file.create_group("genomes")

#all_genomes.apply(lambda (nR,nG,species_name): create_genome_group_and_Q_tables(genomes,species_name,nR,nG),axis=1)

print ctime(), "Gathering read mapped to all genomes"
master_index = dict()
genome_index = dict()
j = 0
pbar = ProgressBar(widgets=widgets, maxval=len(bam_files)*10).start()
for bam_counter,bam in enumerate(bam_files):
	#print ctime(), bam
	try:
		samfile, species_name, nR, nG = openBamfile(bam)
		
		create_genome_group_and_Q_tables(genomes,species_name,nR,nG,samfile.references)
		genome_index[j] = species_name	
		#read_counter = 0
		for read_counter,read in enumerate(samfile.fetch()):
			genomes[species_name]["Q"][read_counter,read.rname] = rescale_samMAPQ( read.mapq )
			genomes[species_name]["read_names"][read_counter,:] = read.qname
			#read_counter += 1
			if read.qname not in master_index.keys():
                        	master_index[ read.qname ] = dict()
                	if species_name in master_index[read.qname]:
                        	master_index[ read.qname ][ j ] = np.append( master_index[ read.qname ][ species_name ] , read_counter)
                	else:
                        	master_index[ read.qname ][ j ] = read_counter

		j += 1
	except:
		print "not working"
	pbar.update(10*bam_counter+1)
pbar.finish()


nR, nG = len(master_index), len(genome_index)

Q_all      = h5_file.create_dataset("Q_all"     ,(nR,nG),dtype='float')
Q_all_rows = h5_file.create_dataset("Q_all_rows",(nR,),dtype=dt_str_vlen)
Q_all_cols = h5_file.create_dataset("Q_all_cols",(nG,),dtype=dt_str_vlen)


print ctime(), "Combining mapping results across genomes"
pbar = ProgressBar(widgets=widgets, maxval=len(master_index.keys())*10).start()
for i,read in enumerate(master_index):
	for j in master_index[read]: 
		Q_all[i,j] = genomes[ genome_index[j] ]["Q"][  master_index[ read ][ j ]  , 0 ] 
	pbar.update(10*i+1)
pbar.finish()

print ctime(), "Done"

