########################
##dealpairedlib.py	  ##
########################
#!/use/bin/env/python3
import os
import pysam
import gzip

def deal_pe_for_one_barcode(barcode_read_dict, length):
    """
    determine the Highly credible paired end by the barcode's read pos information and the 
    appropriate length

    Args :
        barcode_read_dict (a dictionary) : input
        length (int) the length of the max distant of paired end (bp): input

    Return :
        A dict of the barcode's paired end numbers and pos information
    """

    temp = {}
    for rname,pos_info in barcode_read_dict.items():
        if len(pos_info.keys()) < 2 :
            temp[rname] = []
            continue

        if (pos_info['r1']['deal'] == False or pos_info['r2']['deal'] == False):
            temp[rname] = []
            continue

        pos_info_1 = pos_info['r1']['pos']
        pos_info_2 = pos_info['r2']['pos']
        if len(pos_info_2) == 0 or len(pos_info_1) == 0:
            output_info = []
        else:
            output_info = anchored_paired(pos_info_1,pos_info_2,length)
        
        temp[rname] = output_info
        
    return temp


def anchored_paired(pos_info_1,pos_info_2,length):
    """
    determine the Highly credible paired end by diffrent read pos information list and the appropriate
    length
    Args:
        pos_info_1,pos_info_2(list) : input
        length 	(int) as above : input

    Return:
        A list of the read's credible paired end 
    """
    unfindchr = []
    output_info = []
    chr_pos_dict_1 = get_chr_read_all_pos(pos_info_1)
    chr_pos_dict_2 = get_chr_read_all_pos(pos_info_2)
    
    for chr_key,chr_pos_info_1 in chr_pos_dict_1.items():
        if chr_key not in chr_pos_dict_2.keys():
            chrstr = str(chr_key)
            unfindchr.append(chrstr)
            continue
        
        chr_pos_info_2 = chr_pos_dict_2[chr_key]
        
        if len(chr_pos_info_1) <= len(chr_pos_info_2):
            paired_end_list = determine_combination(chr_pos_info_1,chr_pos_info_2,length,1)
        else :
            paired_end_list = determine_combination(chr_pos_info_2,chr_pos_info_1,length,2)
        
        for paired in paired_end_list :
            pos_info = outpaired(chr_key, paired)
            output_info.append(pos_info)

    return output_info


def outpaired(chr_key, paired):
    """
    output the chromosome's one paired end pos infromation
    
    Args:
        chr_key (the chromosome's name)(str) : input
        paired(the paired's pos infromation)(str): input

    Retrun :
       A list of the paired end pos infromation
    """
    pos_info= []
    paired_1 = chr_key+':'+paired.split('#')[0]
    str1 = paired_1.split(':')
    paired_2 = chr_key+':'+paired.split('#')[1]
    str2 = paired_2.split(':')
    pos_info.append(str1[0])
    pos_info.append(int(str1[1]))
    pos_info.append(str1[2])

    pos_info.append(str2[0])
    pos_info.append(int(str2[1]))
    pos_info.append(str2[2])

    return pos_info


def determine_combination(chr_pos_info_1, chr_pos_info_2, length, model):
    """
    determine the chromosome's all Highly credible paired end

    Args:
        chr_pos_info_1 (the read1's all pos information)(list): input
        chr_pos_info_2 (the_read2's all pos information)(list): input
        length (the appropriate length of a paired end) (int)(bp): input
        model (the index of deal way)(int): input

    Return:
        A list chromosome's all Highly credible paired end
    """
    count = 0
    flag = len(chr_pos_info_1)
    paired_end_list = []
    for pos_1 in chr_pos_info_1:
        pos_set_1 = pos_1.split(':')[0]
        direct_1 = pos_1.split(':')[1]
        for pos_2 in chr_pos_info_2 :
            pos_set_2 = pos_2.split(':')[0]
            direct_2 = pos_2.split(':')[1]
            if is_right_direct(direct_1,direct_2) == False :
                continue
            if is_proper_distant(pos_set_1,pos_set_2,length) == False :
                continue
            if count < flag :

                if model == 1 :
                    paired_end = pos_1+'#'+pos_2
                if model == 2 :
                    paired_end = pos_2+'#'+pos_1

                paired_end_list.append(paired_end)
                count += 1

            else : break

    return paired_end_list


def is_right_direct(direct_1, direct_2):
    """
    Judge whether the two direct match
    
    Args:
        direct_1,direct_2(the single end's direct)(str): input
    
    Return:
        True or False
    """
    if (direct_1 == '+' and direct_2 == '-' ) or (direct_1 == '-' and direct_2 == '+'):
        return True
    else :
        return False


def is_proper_distant(pos_set_1, pos_set_2, length):
    """
    Juge whether the two pos is appropriate

    Args:
        pos_set_1,pos_set_2(the pos for the read's)(str): input

    Return:
        True or False
    """
    pos_set_1 = int(pos_set_1)
    pos_set_2 = int(pos_set_2)
    
    if abs(pos_set_1-pos_set_2) <= length:
        return True
    else :
        return False


def get_chr_read_all_pos(pos_info_all):
    """
    Classify pos information according to chromosome names

    Args:
        pos_info_all(the read's all blast pos): input

    Return:
        A dict of the read's all balst chromosome and the pos list
    """

    chr_pos_dict ={}
    for pos in pos_info_all :
        chr_key = pos[0]
        pos_info = str(pos[1])+':'+pos[2]
        if chr_key not in chr_pos_dict.keys() :
            chr_pos_dict[chr_key] = []
        chr_pos_dict[chr_key].append(pos_info)

    return chr_pos_dict


def output_allreadpos_v3(prefix, allreadpos_barcode, index):
    """
    Output the barcode's reads position information according to allreadpos_barcode

    Args :
        allreadpos_barcode (a dictionary) : input
        prefix (str) :the prefix of output file

    Return :
        A file collect the read position information

    """
    f_out = gzip.open(prefix+".all.pos.bar.gz", "wb")
    for barcode, readset in allreadpos_barcode.items():
        b_content = '#'+barcode+'\n'
        f_out.write(b_content.encode())
        for rname, pos_info in readset.items():
            if index in pos_info.keys():
                pos_info_list = list(pos_info[index]['pos'])
                pos_info_all = []
                hit_num = str(pos_info[index]['hitnum'])
                pos_info_all.append(hit_num)                
                for i in pos_info_list :
                    i_new = [str(x) for x in i]
                    pos_i = ':'.join(i_new)
                    pos_info_all.append(pos_i)
 
                pos_all = '\t'.join(pos_info_all)
                r_content = str(rname)+'\t'+pos_all+'\n'
                f_out.write(r_content.encode())

    f_out.close()


def output_allreadpos_paired_end_v2(prefix, allreadpos_barcode):
    """
    This model  as above
    """
    f_out = gzip.open(prefix+".all.pos.bar.gz", "wb")
    for barcode, readset in allreadpos_barcode.items():
        b_content = '#'+barcode+'\n'
        f_out.write(b_content.encode())
        for rname, pos_info in readset.items():
            pos_info_all = []
            hit_num = str(len(pos_info))
            pos_info_all.append(hit_num)
            for i in pos_info :
                i_new = [str(x) for x in i]
                pos_i = ':'.join(i_new)
                pos_info = pos_i.split(':')
                pos_1 = pos_info[0]+':'+pos_info[1]+':'+pos_info[2]
                pos_2 = pos_info[3]+':'+pos_info[4]+':'+pos_info[5]
                pos_i = pos_1+'|'+pos_2
                pos_info_all.append(pos_i)
            pos_all = '\t'.join(pos_info_all)
            r_content = rname+'\t'+pos_all+'\n'
            f_out.write(r_content.encode())

    f_out.close()

