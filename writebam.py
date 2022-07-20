########################
##writebam.py		  ##
########################
#!/user/bin/env/python3
import pysam
import gzip

def add_align_info(align_1, align_2, read_seq_quality_dict):
	"""
	Add other main information for single end alignment transfrom to paired end

	Args:
		align_1(the read1's alignment): input
		align_2(the read2's alignment): input
		read_seq_quality_dict(a dict of the read's sequnce and mapquality): input

	Return:
		The paired end alignment ,which add main information
	"""
	comefile = align_1.get_tag("RG")
	index = comefile.split(".")[1]
	if index == 'r1':
		align_1.cigar = read_seq_quality_dict['r1'][0]
		align_1.query_sequence = read_seq_quality_dict['r1'][1]
		align_1.query_qualities = read_seq_quality_dict['r1'][2]
		align_2.cigar = read_seq_quality_dict['r2'][0]
		align_2.query_sequence = read_seq_quality_dict['r2'][1]
		align_2.query_qualities = read_seq_quality_dict['r2'][2]
		align_1.is_read1=True
		align_2.is_read2=True
	else:
		align_1.cigar = read_seq_quality_dict['r2'][0]
		align_1.query_sequence = read_seq_quality_dict['r2'][1]
		align_1.query_qualities = read_seq_quality_dict['r2'][2]
		align_2.cigar = read_seq_quality_dict['r1'][0]
		align_2.query_sequence = read_seq_quality_dict['r1'][1]
		align_2.query_qualities = read_seq_quality_dict['r1'][2]
		align_1.is_read2=True
		align_2.is_read1=True

	align_1.next_reference_id = align_2.reference_id
	align_1.next_reference_start = align_2.reference_start
	align_2.next_reference_id = align_1.reference_id
	align_2.next_reference_start = align_1.reference_start

	read_2_len = align_2.query_length

	if  align_2.reference_start == align_1.reference_start:
		insert = 0
	else:
		insert = align_2.reference_start + read_2_len - align_1.reference_start

	align_1.template_length	= insert
	align_2.template_length = -insert
	align_1.mapping_quality = 60
	align_2.mapping_quality = 60
	align_1.is_paired = True
	align_2.is_paired = True
    


def get_info_form_barlist(bar_readpos_list):
	"""
	Collect each barcode marked reads to dictionaries

	Args :
		bar_redpos_list : input

	Return :
		three dictionaries.
		count_read_numdict is count the read apper numbers(key is readname,value is the frequnecy),
		read_pos_dict is collect one read's alignment positions(key is readname,value is the positions list),
		bar_dict is marked the barccode lacation in bar_readpos_list.
	"""
	count_read_numdict = {}
	bar_dict = {}
	read_pos_dict = {}

	for i in range(len(bar_readpos_list)):
		bar = bar_readpos_list[i]
		bar_name = bar.bar_name
		bar_dict[bar_name] = i
		for readpos in bar.readpos_list :
			read_name = readpos.read_name
			pos_num = readpos.pos_num
			count_read_numdict[read_name] = int(pos_num)
			read_pos_dict[read_name] = []
			for onepos in readpos.pos_list:
				pos_info= onepos.get_pos_str()
				read_pos_dict[read_name].append(pos_info)

	return count_read_numdict, bar_dict, read_pos_dict


def update_pos_info(bar_readpos_list, bam_file, prefix):
	"""
	Update the read alignment information to new_bam according to bar_readpos_list
		
	Args :
		bar_readpos_list (new bar_readpos_list after solve multip align ): input
		bamfile (the megered align bam file ): input
		prefix (str) : the prefix of output bamfile

	Return :
		the update bam file after run solve_multi_barcode_v2.py
	"""
	count_read_numdict, bar_dict, read_pos_dict = get_info_form_barlist(bar_readpos_list)

	alignments = pysam.AlignmentFile(bam_file, 'rb')
	outfile = pysam.AlignmentFile(prefix+".bam", "wb", template = alignments)

	new_read_numdict = {}
	read_seq_quality = {}
	read_index = {}
	align_1 = pysam.AlignedSegment()
	align_2 = pysam.AlignedSegment()

	for align in alignments.fetch(until_eof=True):
		query_name = align.query_name
		rname = query_name.split('#')[0]
		bname = query_name.split('#')[1]

		if rname in read_index.keys():
			continue

		supple_index = align.is_supplementary
		dup_index = align.is_duplicate
		unmap_index = align.is_unmapped

		if supple_index == True or dup_index == True or unmap_index == True:
			outfile.write(align)
			continue

		if bname not in bar_dict.keys():  #read's barcode not in bar_dict output directly
			outfile.write(align)
			continue

		if rname not in count_read_numdict.keys():  # read not in count_read_numdict output directly
			outfile.write(align)
			continue

		if count_read_numdict[rname] == 0 :  # read alignment information not in bar_readpos_list after data processing
			outfile.write(align)
			continue
		
		come_from = align.get_tag("RG")
		readfile = come_from.split('.')[1]

		if count_read_numdict[rname] == 1:
			if (rname not in read_seq_quality.keys()) and (align.query_sequence != None ):
				read_seq_quality[rname]={}
				read_seq_quality[rname][readfile]=[]
				read_seq_quality[rname][readfile].append(align.cigar)
				read_seq_quality[rname][readfile].append(align.query_sequence)
				read_seq_quality[rname][readfile].append(align.query_qualities)
			elif (rname in read_seq_quality.keys()) and (readfile not in read_seq_quality[rname].keys()) and (align.query_sequence != None):
				read_seq_quality[rname][readfile]=[]
				read_seq_quality[rname][readfile].append(align.cigar)
				read_seq_quality[rname][readfile].append(align.query_sequence)
				read_seq_quality[rname][readfile].append(align.query_qualities)
	
			if is_correct_pe(align, read_pos_dict[rname]) == True :
				if rname not in new_read_numdict.keys():
					new_read_numdict[rname] = 1
					align_1 = align
				else :
					new_read_numdict[rname] += 1
					align_2 = align

			if rname in read_seq_quality.keys() and len(read_seq_quality[rname].keys()) == 2:
				if rname in new_read_numdict and new_read_numdict[rname] == 2 :
					add_align_info(align_1, align_2, read_seq_quality[rname])
					outfile.write(align_1)
					outfile.write(align_2)
					read_index[rname] = 1
					read_seq_quality.pop(rname)
					new_read_numdict.pop(rname)
		else :					
			outfile.write(align)		# the mutiple align reads after processing output directly
	outfile.close()
	alignments.close()	
	return


def is_correct_pe(align, read_pos_list):
	"""
	Determine if the read matches correctly

	Args:
		align (the read alignment information): input
		read_pos_list (list stored one read's align) : input

	Return:
		True or False
	"""
	index = 0

	ref_name = align.reference_name
	pos_1 = align.pos

	if align.is_reverse == True :
		direct_1 = '-'
	else:
		direct_1 = '+'

	onepos = read_pos_list[0]
	paired_end = onepos.split('|')

	read_info_1 = paired_end[0]
	ref_1 = read_info_1.split(':')[0]
	pos_info_1 = int(read_info_1.split(':')[1])
	pos_direct_1 = read_info_1.split(':')[2]
	read_info_2 = paired_end[1]
	ref_2 = read_info_2.split(':')[0]
	pos_info_2 = int(read_info_2.split(':')[1])
	pos_direct_2 = read_info_2.split(':')[2]

	if (ref_name == ref_1 and pos_1 == pos_info_1 and direct_1 == pos_direct_1 ) or (ref_name == ref_2 and pos_1 == pos_info_2 and direct_1 == pos_direct_2 ) :
		index = 1

	if index == 1 :
		return True
	else :
		return False
