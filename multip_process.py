########################
##multip_process.py	  ##
########################
from dealpairedlib import *
from data_prepare import bwa_align_and_merge
from bar_lfr_read import *
from bar_read_pos_v2 import *
from readbam import *
from writebam import *
import multiprocessing
import os

def multip_run(tools_path, samtools, bwa, ref_file, num_file, max_workers, mhit, length, lfr_thrshd):
    """
    Parallel processing of solve multiple alignment process

    Args:
        tools_path (the tools path) :input
        samtools (the samtools path): input
        bwa (the bwa path) : input
        ref_file (the reference sequence): input
        num_file (the specified number of files)(int): input
        max_workers (the specified number of threads)(int): input
        mhit (the max hit numbers of a read pos)(int): input
        length (the max length of PE insert distant can be accepted ):input
        lfr_thrshd (the distance threshold is used to determine whether two neighboring reads are in the same lfr): input

    Return:
        the all bam file of solve multiple alignment process
    """
    pool = multiprocessing.Pool(processes = max_workers)
    for i in range(num_file):
        pool.apply_async(multip_process_all_barcode,(tools_path, samtools, bwa, ref_file, "read_1_floder/"+str(i)+"_1.fq", "read_2_floder/"+str(i)+"_2.fq", str(i), mhit, length, lfr_thrshd, ))
    pool.close()
    pool.join()


def multip_process_all_barcode(tools_path, samtools, bwa, ref_file, read1, read2, number, mhit, length, lfr_thrshd):
    """
    One process of solve multiple alignment process

    Args:
        tools_path (the tools path): input
        samtools (the samtools path) : input
        bwa (the bwa path) : input
        ref_file (the reference sequence): input
        read_1,read_2(the pair end fastq file) : input
        number (the processing number): input
        mhit (as above): input
        length (as above)：input
        lfr_thrshd (as above) :input

    Return:
        One bam file of solve multiple alignment process   
    """
    bwa_align_and_merge(tools_path, samtools, bwa, ref_file, read1, read2, number)
    new_readpos_barcode = {}
    allreadpos_barcode = get_all_read_pos_set_v2("merge_bam/"+number+".sorted.bam", mhit)
    output_allreadpos_v3('r1_pos_info/'+number+'.r1',allreadpos_barcode, 'r1')
    output_allreadpos_v3('r2_pos_info/'+number+'.r2',allreadpos_barcode, 'r2')

    index = list(allreadpos_barcode.keys())
    bar_readpos_list = []
    for barcode in index:
         new_readpos_barcode[barcode] = deal_pe_for_one_barcode(allreadpos_barcode[barcode], length)
         allreadpos_barcode.pop(barcode)
         bar_readpos = transfrom_dict_to_class(new_readpos_barcode[barcode],barcode)
         new_readpos_barcode.pop(barcode)
         new_bar_readpos = solve_multi_1bar(bar_readpos, lfr_thrshd)
         bar_readpos_list.append(new_bar_readpos)

    output_update_info("Solved_Bar_info/"+number, bar_readpos_list)
    update_pos_info(bar_readpos_list, "merge_bam/"+number+".sorted.bam", "Update_bam_floder/"+number)

    del index
    del new_readpos_barcode
    del bar_readpos_list

    return 


def transfrom_dict_to_class(one_bar_dict , bar_name):
    """
    Transfrom the barcode's pos information dict to the BarRead class

    Args:
        one_bar_dict(the dict of a barcode's read pos information): input 
        bar_name(the barcode name)(str): input

    Return:
        A BarRead class
    """
    bar_readpos = BarRead(bar_name,[])
    for rname ,pos_info in one_bar_dict.items():
        pos_num = len(pos_info)
        pos_list = []
        for i in pos_info:
            i_new = [str(x) for x in i]
            pos_i = ':'.join(i_new)
            pos_i = pos_i.split(':')
            pos_1 = pos_i[0]+':'+pos_i[1]+':'+pos_i[2]
            pos_2 = pos_i[3]+':'+pos_i[4]+':'+pos_i[5]
            pos_i = pos_1+'|'+pos_2
            one_pos = OnePos(pos_i)
            pos_list.append(one_pos)

        readpos = ReadPos(rname,pos_num,pos_list)
        bar_readpos.readpos_list.append(readpos)

    return bar_readpos
