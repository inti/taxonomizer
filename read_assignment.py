import pysam
import numpy as np
#import tables as hdf5
import pandas as pd 
import subprocess


aln_folder = "/media/TeraData/ipedroso/METAGENOMICS/ANALYSES/ALN/simLC/"

proc = subprocess.Popen(["".join(("ls ",aln_folder,"/Z*bam"))], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()

bam_files = out.split()
all_genomes = np.hstack([pysam.AlignmentFile(bam, "rb").references for bam in bam_files])


Q = pd.DataFrame(columns=all_genomes,dtype=np.float16)

flush = 100

for bam in bam_files:
	samfile = pysam.AlignmentFile(bam, "rb")
	nR = samfile.count()
	local_Q = pd.DataFrame(np.zeros((flush , samfile.nreferences)),columns=samfile.references,dtype=np.float16)
	reads_id = {}
	read_counter = 0
	flush_counter = 0
	for read in samfile.fetch():
    		if read.qname not in reads_id:
			reads_id[read.qname] = read_counter
			read_counter += 1

		local_Q.ix[flush_counter,read.rname] = read.mapq
		flush_counter += 1
		
		if flush_counter > 99:
			# here we need to add local_Q to Q with the align methods of 
			flush_counter = 0
			print local_Q.head()
			local_Q = pd.DataFrame(np.zeros((flush , samfile.nreferences)),columns=samfile.references,dtype=np.float16)
			continue



