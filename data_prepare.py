import re
import os
import gzip

def split_read_by_barcode(filename, max_workers):
	"""
	Classify the barcode's read name

	Args:
		filename(the fastq file): input 
		max_workers(specified number of separations)(int): input

	Return :
		the specified numbers file of read's name
	"""
	of = gzip.open(filename,'rb')
	barcode_dict = {}
	barcode_list = []
	barcode_count = 0
	index = 1
	for line in of.readlines():
		if(index%4==1):
			line = line.decode()
			read = line.split('\t')[0].strip()
			read = read.lstrip('@')
			read_list = re.split('[#,/]',read)
			barcode = read_list[1]
			if barcode not in barcode_dict.keys():
				barcode_count += 1
				barcode_dict[barcode]=[]
			barcode_dict[barcode].append(read)
		index += 1
	if barcode_count < max_workers:
		max_workers = barcode_count	
	for i in range(max_workers):
		barcode_list.append([])
	i=0
	for bar,read_list in barcode_dict.items():
		barcode_list[i].extend(read_list)
		i += 1
		if i%max_workers == 0:
			i=0
	for i in range(max_workers):
		of = open("read_list/"+str(i)+'_read_list.txt','w')
		for read in barcode_list[i]:
			of.write(read+'\n')
		of.close
	
	return max_workers


def create_work_floder():
	"""
	Create the work floders

	Return:
		the work floders
	"""
	os.system("mkdir read_list")
	os.system("mkdir read_1_floder")
	os.system("mkdir read_2_floder")
	os.system("mkdir r1_bam")
	os.system("mkdir r2_bam")
	os.system("mkdir merge_bam")
	os.system("mkdir r1_pos_info")
	os.system("mkdir r2_pos_info")
	os.system("mkdir Solved_Bar_info")
	os.system("mkdir Update_bam_floder")


def extract_read_by_name(tool_path, read_1, read_2, max_workers, threads):
	"""
	Extract read form the fastq file to new fastq file

	Args:
		tool_path (the tools path): input 
		read_1,read_2 (the fastq file of sequencing data): input
 		max_workers (the specified number of copies): input
		thread (the specified number of thread): input

	Return :
		The specified number fastq file of paired end
	"""
	os.system("bash "+tool_path+"/shell/extract_read.sh"+" -r "+read_1+" -n read_list"+" -o read_1_floder"+" -m "+str(max_workers)+" -t "+str(threads)+" -i 1")
	os.system("bash "+tool_path+"/shell/change_read_name.sh"+" -n read_list"+" -m "+str(max_workers)+" -t "+str(max_workers))
	os.system("bash "+tool_path+"/shell/extract_read.sh"+" -r "+read_2+" -n read_list"+" -o read_2_floder"+" -m "+str(max_workers)+" -t "+str(threads)+" -i 2")


def bwa_align_and_merge(tool_path, samtools, bwa, ref_file, read1, read2, index):
    """
    Using bwa to Alignment the pair end to reference sequence

    Args:
        tool_path (the tools path): input
        samtools (the samtools path) : input
        bwa (the bwa path) : input
        ref_file (the reference sequence): input
        read1,read2(the fastq file of paired end): input
        index (the number of the fastq)(int): input
	
    Return:
        the bam file of the alignment result
    """
    os.system("bash "+tool_path+"/shell/single_bwa.sh"+" -s "+samtools+" -b "+bwa+" -r "+ref_file+" -f "+read1+" -o r1_bam"+" -m "+str(index)+" -i 1")
    os.system("bash "+tool_path+"/shell/single_bwa.sh"+" -s "+samtools+" -b "+bwa+" -r "+ref_file+" -f "+read2+" -o r2_bam"+" -m "+str(index)+" -i 2")
    os.system("bash "+tool_path+"/shell/merge_bam.sh"+" -s "+samtools+" -1 r1_bam"+" -2 r2_bam"+" -o merge_bam"+" -m "+str(index))


def merge_bam_after_sloved(tool_path, samtools ,threads):
	"""
	Merge the bam file after solve multiple alignment process

	Args:
		tool_path (the tools path): input
		threads (the specified number of thread): input

	Return:
		A bam file of the final result of alignment
	"""
	os.system("bash "+tool_path+"/shell/merge_bam_after_update.sh"+" -s "+samtools+" -f Update_bam_floder -o stLFR_SMA -t "+str(threads))


def detele_intermediate_file():
	"""
	Detele the intermediate file
	"""
	os.system("rm -rf read_list")
	os.system("rm -rf read_1_floder")
	os.system("rm -rf read_2_floder")
	os.system("rm -rf r1_bam")
	os.system("rm -rf r2_bam")
	os.system("rm -rf Update_bam_floder")
	os.system("rm -rf merge_bam")
