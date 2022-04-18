import os,sys,time
import gzip
import logging
from numpy import *
from bar_lfr_read import *
import multiprocessing

"""
This version takes a read with multiple alignments as an unit to solve the multiple aligment read.
"""


def get_cddt_lfr_1seq(uniq_readset_1seq, lfr_thrshd):
    """
    clustering the alignments on the same sequence to get the candidate long fragments.

	Args:
        uniq_readset_1seq (a list of readpos on the same sequence and these reads are with the same barcode) : input
        lfr_thrshd (The distance threshold is used to determine whether two neighboring reads are in the same lfr) : input

    Return:
        A list of long fragments.
    """

    def get_pos(readpos):
        tmp_pos = readpos.pos_list[0].get_pos_value()
        return tmp_pos

    uniq_readset_1seq.sort(key=get_pos)

    lfr_read_list = []
    lfr_read = LfrRead([], lfr_thrshd)
    lfr_count = 0
    lfr_read_num = 0
    for readpos in uniq_readset_1seq:
        if lfr_read_num ==0:
            lfr_read.readpos_list.append(readpos)
            lfr_read_num = lfr_read_num + 1
            lfr_count = lfr_count + 1
        else:
            bf_readpos = lfr_read.readpos_list[lfr_read_num-1]
            pos_bf = get_pos(bf_readpos)
            pos_self = get_pos(readpos)

            pos_diff = int(pos_self) - int(pos_bf)

            if pos_diff <= lfr_thrshd:
                lfr_read.readpos_list.append(readpos)
                lfr_read_num = lfr_read_num + 1
            else:
                lfr_read_list.append(lfr_read)
                lfr_read = LfrRead([], lfr_thrshd)
                lfr_read_num = 0

                lfr_read.readpos_list.append(readpos)
                lfr_read_num = lfr_read_num + 1
                lfr_count = lfr_count + 1

    lfr_read_list.append(lfr_read)

    return lfr_read_list 


def get_readset_byname(uniq_readset):
    """
    Divide a readset with uniqe alignments into subsets by the sequence name in OnePos format.

    Args:
        uniq_readset(a list of the barcode's unique alignment read position and name): input
	
    Return:
        A dict about barcode's unique read alignment chromosome and position information
    """
    name_dict ={}
    for readpos in uniq_readset:
        tmp_name = readpos.pos_list[0].get_seq_name()
        if tmp_name not in name_dict.keys():
            name_dict[tmp_name] = []
            name_dict[tmp_name].append(readpos)
        else:
            name_dict[tmp_name].append(readpos)

    return name_dict
            


def get_cddt_lfr_1bar(uniq_readset, lfr_thrshd):
    """
    Divide the reads of same barcode into different clusters which refer to candidate LFR. Firstly, dividing by sequence name; secondly, dividing by postion.
	
    Args:
        uniq_readset : (a list of the barcode's all unique alignment read position and name) : input
        lfr_thrshd: the same to that in get_cddt_lfr_1seq()
    
    Return:
        A list of the candidate long fragments
    """
    lfr_read_list_all = []
    name_dict = get_readset_byname(uniq_readset)
    for seq_name in name_dict.keys():
        tmp_uniq_readset_1seq = name_dict[seq_name]
        tmp_lfr_read_list = get_cddt_lfr_1seq(tmp_uniq_readset_1seq, lfr_thrshd)
        lfr_read_list_all.extend(tmp_lfr_read_list)

    return lfr_read_list_all


def get_readset(lfr_read_list):
    """
    Get the a list of read position from a list of LfrRead.

    Args:
        lfr_read_list(the candidate long fragments list): input 

    Return:
        A list of all long fragments's read position
    """
    readset = []
    for lfr_read in lfr_read_list:
        readset.extend(lfr_read.readpos_list)

    return readset


def solve_multi_1bar(bar_readpos, lfr_thrshd):
    """
    A complete process to determine the solid aligment of the multiple reads in one barcode which contains both uniqe and multiple reads.
    
    Args:
        bar_readpos(a BarRead class contains all the position information of all reads with a same barcode) : input
        lfr_thrshd (the same to above) : input
   
    Return:
        A new BarRaed class after solved multiple alignment
    """
    valid_index = bar_readpos.get_valid_index()
    if valid_index == False:
        return bar_readpos

    uniq_readset = bar_readpos.get_uniq_readset()
    multi_readset = bar_readpos.get_multi_readset()

    new_readset = []
    new_multi_readset = multi_readset
    lfr_read_list_all = get_cddt_lfr_1bar(uniq_readset, lfr_thrshd)
    new_lfr_read_list_all = lfr_read_list_all

    p_index = 1

    while p_index == 1:
        p_index, update_multi_readset, update_lfr_read_list_all = solve_multi_1loop(new_multi_readset, new_lfr_read_list_all, lfr_thrshd)

        new_multi_readset = update_multi_readset
        new_lfr_read_list_all = update_lfr_read_list_all

    tmp_uniq_readset = get_readset(new_lfr_read_list_all)
    new_readset = tmp_uniq_readset
    new_readset.extend(new_multi_readset)


    new_bar_name = bar_readpos.bar_name
    new_readpos_list = new_readset
    new_bar_readpos = BarRead(new_bar_name, new_readpos_list)
    
    return new_bar_readpos


def solve_multi_1loop(multi_readset, lfr_read_list ,lfr_thrshd):
    """
    A loop of solve multiple alignment process
    
    Args:
        multi_readset (the barcode's all multiple alignment read list): input
        lfr_read_list (the candidate long fragments list): input
        lfr_thrshd (as above): input

    Return:
        index(0 or 1): Judge whether run next loop
        update_multi_readset: the mulitple alignment reads' list after a loop process
        new_lfr_read_list: the unqiue reads' long fragments after a loop process
    """
    p_index = 1
    sol_num = 0
    update_multi_readset = []

    new_lfr_read_list = lfr_read_list

    for readpos in multi_readset:
        s_index,update_lfr_read_list = solve_multi_1read(readpos, new_lfr_read_list,lfr_thrshd)

        if s_index == 1:
            sol_num = sol_num +1
            new_lfr_read_list = update_lfr_read_list
        else:
            readpos.m2u_index = s_index
            update_multi_readset.append(readpos)

    if sol_num > 0:
        p_index = 1
    else:
        p_index = 0

    return p_index, update_multi_readset, new_lfr_read_list


def get_pos_by_same_chr(readpos):
    """
    Classify the read pos by chromosome 

    Args:
        readpos(class of readpos): input
    
    Return:
        A dict of chromosomes name and read's pos list
    """
    chr_pos_dict = {}
    for pos in readpos.pos_list:
        chr_name = pos.get_seq_name()
        pos_value = int(pos.get_pos_value())
        if chr_name not in chr_pos_dict.keys():
           chr_pos_dict[chr_name]=[]
        chr_pos_dict[chr_name].append(pos_value)

    return chr_pos_dict


def is_inter_read(pos_list, lfr_thrshd):
    """
    Judge a read's all pos is inter in lfr

	Args:
        pos_list(a read's pos list): input
        lfr_thrshd(the distance threshold is used to determine whether two neighboring reads are in the same lfr): input

    Return:
        True or False
    """
    min_distant = -1
    num = 0
    temp_value = -1
    for pos in pos_list:
        if num > 0 :
            close_distant = abs(pos - temp_value)
            if close_distant < lfr_thrshd/2 :
               return True
            else :
               old_value = temp_value
        if num > 1 :
           if abs(temp_value - old_value) < abs(pos - temp_value):
               min_distant = abs(temp_value - old_value)
           else :
               min_distant = abs(pos - temp_value)

           if min_distant < lfr_thrshd/2 :
               return True

        temp_value = pos
        num += 1

    return False


def solve_multi_1read(readpos, lfr_read_list, lfr_thrshd):
    """
    Deal one read for all multiple alignment reads
    
    Args:
        readpos (A multiple alignment read's all positions): input
        lfr_read_list (all long fragments’ list of a barcode): input
        lfr_thrshd (as above): input

    Return:
        s_index(the index after a process)
        update_lfr_read_list(all long fragments' list of a barcode after update)
    """
    update_lfr_read_list = []
    update_readpos = readpos
    pos_lfr_read_dict = {}

    s_index = None

    chr_pos_dict = get_pos_by_same_chr(readpos)
    for chr_name,pos_list in chr_pos_dict.items():
        if is_inter_read(pos_list, lfr_thrshd) == True:
            s_index = 3
            update_lfr_read_list = lfr_read_list
            return s_index, update_lfr_read_list

    if len(readpos.pos_list) != int(readpos.pos_num):
        s_index = 6
        update_lfr_read_list = lfr_read_list

        return s_index, update_lfr_read_list

    else:
        for pos in readpos.pos_list:
            for lfr_read in lfr_read_list:
                lfr_index = is_pos_in_lfr(pos, lfr_read)

                if lfr_index == 1:
                    if lfr_read not in pos_lfr_read_dict.keys():
                        pos_lfr_read_dict[lfr_read]=[]
                    pos_lfr_read_dict[lfr_read].append(pos)
                    break

        
        s_index = get_s_index(pos_lfr_read_dict)

        if s_index == 1:
            tmp_update_lfr_read_list = list(pos_lfr_read_dict.keys())
            update_lfr_read = tmp_update_lfr_read_list[0]

            update_readpos.pos_num  = 1
            update_readpos.pos_list = pos_lfr_read_dict[update_lfr_read]
            update_readpos.m2u_index = 1

            for lfr_read in lfr_read_list:
                if lfr_read == update_lfr_read:
                    update_lfr_read.readpos_list.append(update_readpos)
                    update_lfr_read_list.append(update_lfr_read)

                else:
                    update_lfr_read_list.append(lfr_read)
        else:
            update_lfr_read_list = lfr_read_list

        return s_index, update_lfr_read_list


def is_pos_in_lfr(pos, lfr_read):
    """
    Judge the read's position is in the long fragments' internal

    Args:
        pos (read all alignment positiona) : input
        lfr_read(as lfr_read_list) : input
        
    Return:
        0 or 1 
    """
    pos_seq_name = pos.get_seq_name()
    pos_pos_value = pos.get_pos_value()

    lfr_seq_name = lfr_read.get_seq_name()
    lfr_start =lfr_read.get_start()
    lfr_end = lfr_read.get_end()
    
    lfr_index = None
    if pos_seq_name == lfr_seq_name:
        if float(pos_pos_value) >= float(lfr_start) and float(pos_pos_value) <= float(lfr_end):
            lfr_index = 1
        else:
            lfr_index = 0
    else:
        lfr_index = 0

    return lfr_index


def get_s_index(pos_lfr_read_dict):
    """
    To analyze the distribution of multiple alignments on one group LFR. The map between postions and LFR are recorded in a dictionary.

    Args:
        pos_lfr_read_dict : input

    Return:
        A index

    Notes:
        s_index: 6 refers to unsolved case for the reads with incomplete postions; 1 refers to solved case; 2 refers to unsolved case for intra-lfr repeats; 3 refers to unsolved case for inter-lfr repeats; 4 refers to unsolved case for mixture repeats; 5 refers to unsolved case for the read whose positions are fall out from all the LFRs; 0 refers to the reads whose s_index has no change in the process.  
    """
    max_pos_num = 0
    lfr_num = len(list(pos_lfr_read_dict.keys()))
    for lfr_read in pos_lfr_read_dict.keys():
        tmp_pos_num = len(pos_lfr_read_dict[lfr_read])
        if tmp_pos_num > max_pos_num:
            max_pos_num = tmp_pos_num


    if lfr_num == 1 and max_pos_num ==1:
        s_index = 1
    elif lfr_num ==1 and max_pos_num > 1:
        s_index = 2
    elif lfr_num >1 and max_pos_num ==1:
        s_index = 3
    elif lfr_num > 1 and max_pos_num > 1 :
        s_index = 4
    else:
        s_index = 5

    return s_index


def output_update_info(prefix, new_bar_readpos_list):
    """
    Output all the position information of a barcode 

    Args：
        prefix(the output file prefix)(str) : input
        new_bar_readpos_list(a list of the barcode all read's position infromation): input

    Return:
        A gz file of the solved multiple alignment position information for a barcode
    """
    pout = gzip.open(prefix+".update.bar.info.gz", "wb")
    for bar in new_bar_readpos_list:
        pout_bar(bar, 1, pout)
    pout.close()


def pout_bar(bar_readpos, index, pout):
    """
    Output all the position information of a barcode in two modes controlled by index.

    Args:
        bar_readpos (a BarRead class record all the position information): input
        index (1 indicates outputing m2u_index, 0 is not): input
        pout (output stream class of gz file) : input

    Return:
        A gz file of the solved multiple alignment position information for a barcode
    """
    b_content = "#"+bar_readpos.bar_name+"\n"
    pout.write(b_content.encode())
    for readpos in bar_readpos.readpos_list:
        pout_readpos(readpos, index, pout)

    return


def pout_readpos(readpos, index, pout): 
    """
    Output all the alignment information of a read in two modes controlled by index.
  
    Args:
        readpos(a ReadPos class record all the aligment information) : input
        index (1 indicates outputing m2u_index, 0 is not) : input
        pout (output stream class of gz file) :input

    Return:
        A gz file of the solved multiple alignment position information for a read
    """
    pos_str_list = []
    for pos_info in readpos.pos_list:
        pos_str_i = pos_info.get_pos_str()
        pos_str_list.append(pos_str_i)
            
    allpos = "\t".join(pos_str_list)
    if index == 1:
        tmp_content = readpos.read_name + "\t" + str(readpos.pos_num) + "\t"+ allpos + "\t" + str(readpos.m2u_index)+ "\n"
        pout.write(tmp_content.encode())
    elif index ==0:
        tmp_content = readpos.read_name + "\t" + str(readpos.pos_num) + "\t"+ allpos + "\n"
        pout.write(tmp_content.encode())
    else: 
        print("Wrong index is set.")


def pout_lfr_read(lfr_read, prefix, pout):
    """
    Output all the position information of a lfr without m2u_index.

    Args: 
        lfr_read (a LfrRead class record all the position information): input
        pout (output stream class of gz file) : input

    """
    seq_name = lfr_read.get_seq_name()
    start = lfr_read.get_start()
    end = lfr_read.get_end()
    thrshd = lfr_read.thrshd_len

    p_content = ">"+seq_name + str(start)+ str(end)+ str(thrshd)+"\n"
    pout.write(p_content.endcode())
    for readpos in lfr_read.readpos_list:
        pout_readpos(readpos, 0)

    return

