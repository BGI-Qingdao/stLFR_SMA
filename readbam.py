########################
##readbam.py		  ##
########################
import pysam

def get_all_read_pos_set_v2(file_name, mhit):
    """
    Getting all read alignment positions(unique and multiple) of a barcode for all barcodes from a bam file

    Args:
        file_name (the bam file) :input
        mhit (int ,limit number of positions) : input

    Return:
        A dictionary ,which stored all barcode's read positions information(key is barcode name ,value is the reads position dictionary)
        this dictionary nested reads information dictionary
    """
    allreadpos_barcode = {}
    alignments = pysam.AlignmentFile(file_name, 'rb')
    for align in alignments.fetch(until_eof=True):
        anal_index = 1
        query_name = align.query_name
        rname = query_name.split('#')[0]
        bname = query_name.split('#')[1]
        file_index = align.get_tag("RG")
        index = file_index.split(".")[1]
        supple_index = align.is_supplementary
        dup_index = align.is_duplicate
        unmap_index = align.is_unmapped
        if supple_index == True or dup_index == True or unmap_index == True :
            anal_index = 0

        if anal_index == 1 :
            if bname not in allreadpos_barcode.keys() :
                allreadpos_barcode[bname]={}

            if rname not in allreadpos_barcode[bname].keys():
                allreadpos_barcode[bname][rname]={}
                allreadpos_barcode[bname][rname][index]={'hitnum':0,'pos': [],'deal':True}
            elif index not in allreadpos_barcode[bname][rname].keys():
                allreadpos_barcode[bname][rname][index]={'hitnum':0,'pos': [],'deal':True}

            if allreadpos_barcode[bname][rname][index]['hitnum'] == mhit+1:
               allreadpos_barcode[bname][rname][index]['hitnum'] = mhit
               allreadpos_barcode[bname][rname][index]['pos'].pop()
               allreadpos_barcode[bname][rname][index]['deal'] = False
            if allreadpos_barcode[bname][rname][index]['deal'] == False:
               continue

            pos_info = get_main_info(align)
            allreadpos_barcode[bname][rname][index]['hitnum'] += 1
            allreadpos_barcode[bname][rname][index]['pos'].append(pos_info)

    return allreadpos_barcode


def get_main_info(align):
    """
    get the main information of a alignment

    Args:
        align(a read's alignment): input

    Return:
        A list of the alignment's mian information ,which include the chromosome name ,the pos and the direct
    """
    pos_info = []

    ref_name = align.reference_name

    pos_info.append(ref_name)

    pos = align.reference_start

    pos_info.append(pos)

    reverse_index = align.is_reverse

    if reverse_index == True :
        direct = '-'
    else:
        direct = '+'

    pos_info.append(direct)

    return pos_info

