import pysam
import numpy as np


samfile = pysam.AlignmentFile("/media/TeraData/ipedroso/METAGENOMICS/ANALYSES/ALN/simLC/Zymomonas_mobilis_CP4___NRRL_B_14023_uid229874.sorted.bam", "rb")

# get number of genomes and of reads on file
nG = len(samfile.header["SQ"])
info_seq_G = samfile.header["SQ"]
nR = samfile.count()

Q = np.zeros((nR,nG),dtype=np.float16)


reads_id = {}
read_counter = 0
for read in samfile.fetch():
    if read.qname not in reads_id:
	reads_id[read.qname] = read_counter
	read_counter += 1
    Q[reads_id[read.qname],read.rname] = read.mapq

