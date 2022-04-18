import os,sys,time
import gzip
from numpy import *

class OnePos:
    def __init__(self, pos_str):
        """
        Description: this class record the position information of an aligment in the form "seq_name:pos_value:direct"
        """
        self.pos_info = pos_str
        
        if self.check_valid() == False :
            print("there exist error in format of recording position information")
            exit(-1)

    def check_valid(self):
        pos_info = self.pos_info.split("|")
        pos_1 = pos_info[0]
        pos_2 = pos_info[1]
        pos_info = pos_1+":"+pos_2
        pos_info = pos_info.split(":")
        if len(pos_info) != 6 :
            return False
        else:
            return True

    def get_pos_str(self):
        pos_info = self.pos_info
        return pos_info

    def get_seq_name(self):
        pos_info = self.pos_info.split(":")
        return pos_info[0]

    def get_pos_value(self):
        pos_info = self.pos_info.split(":")
        return pos_info[1]

    def get_direct(self):
        pos_info = self.pos_info.split("|")
        pos_info = pos_info[1]
        pos_info = pos_info.split(":")
        return pos_info[2]


class ReadPos:
    def __init__(self, read_name, pos_num, pos_list):
        """
        Description: This class record all the alignment postion information of a read.
        read_name:
        pos_num: number of alignment position.
        pos_list: it must be a list of OnePos class. 
        """
        self.read_name = read_name
        self.pos_num = pos_num
        self.pos_list = pos_list
        self.m2u_index = 0


class BarRead:
    def __init__(self, bar_name, readpos_list):
        """
        Description: This class record all the aligment information of reads with the same barcode name.
        readpos_list: it must be a list of ReadPos class.
        """
        self.bar_name = bar_name
        self.readpos_list = readpos_list

    def get_uniq_read_num(self):
        uniq_count = 0
        for readpos in self.readpos_list:
            if int(readpos.pos_num) == 1:
                uniq_count = uniq_count +1

        return uniq_count

    def get_read_num(self):
        tot_num = len(self.readpos_list)
        return tot_num

    def get_multi_read_num(self):
        multi_count = 0
        for readpos in self.readpos_list:
            if int(readpos.pos_num) > 1:
                multi_count = multi_count + 1

        return multi_count


    def get_valid_index(self):
        uniq_num = self.get_uniq_read_num()
        multi_num = self.get_multi_read_num()

        if uniq_num >0 and multi_num >0 :
            index = True
        else:
            index = False

        return index


    def get_uniq_readset(self):
        uniq_readset = []
        for readpos in self.readpos_list:
            if int(readpos.pos_num) == 1 :
                uniq_readset.append(readpos)

        return uniq_readset

    
    def get_multi_readset(self):
        multi_readset  = []
        for readpos in self.readpos_list:
            if int(readpos.pos_num) > 1 :
                multi_readset.append(readpos)

        return multi_readset


class LfrRead:
    def __init__(self, readpos_list, thrshd_len):
        """
        Description: This class can be used to get LFR postion information (seq_name, start, end) by the reads with unique alignment position.
        readpos_list: it must be a list of ReadPos class.
        thrshd_len: it can be set by the mean length of LFR in experiments. 
        """

        self.readpos_list = readpos_list
        self.thrshd_len  = thrshd_len

    def get_seq_name(self):
        name_list = []
        for readpos in self.readpos_list:
            tmp_name = readpos.pos_list[0].get_seq_name()
            name_list.append(tmp_name)

        if len(set(name_list)) == 1:
            seq_name = name_list[0]
        else:
            print("the read set of LFR come from two and more different sequnce")
            exit(-1)

        return seq_name


    def get_start(self):
        pos_list = []
        for readpos in self.readpos_list:
            tmp_pos = readpos.pos_list[0].get_pos_value()
            pos_list.append(float(tmp_pos))
        
        average_pos = round(mean(pos_list), 3)
        start = average_pos - self.thrshd_len
        if start <= 0:
            start = 0

        return start


    def get_end(self):
        pos_list = []
        for readpos in self.readpos_list:
            tmp_pos = readpos.pos_list[0].get_pos_value()
            pos_list.append(float(tmp_pos))
        
        average_pos = round(mean(pos_list), 3)
        end = average_pos + self.thrshd_len

        return end

    def get_direct(self):
        direct = "+"

        return direct
