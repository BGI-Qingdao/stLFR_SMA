import sys
import re
from multip_process import *
import os,sys,time,argparse
import gzip
from us_log import *
from data_prepare import *

if __name__=="__main__":
	tool_path = sys.path[0]
	max_workers = 100
	parser = argparse.ArgumentParser(description='Sloved the mutilp alignment reads')
	parser.add_argument('-ref', required=True, dest='reference', type=str, help='the reference file')
	parser.add_argument('-1', required=True, dest='read_1', type=str, help='the paired_end of read1 must be .gz')
	parser.add_argument('-2', required=True, dest='read_2', type=str, help='the paired_end of read2 must be .gz')
	parser.add_argument('-mhit', required=True, dest='max_hit', type=int, help='The max hit numbers of a read pos')
	parser.add_argument('-length', required=True, dest='insert_distant', type=int, help='The max length of PE insert distant can be accepted')
	parser.add_argument('-lfr_thrshd', required=True, dest='lfr_thrshd', type=int, help='The distance threshold is used to determine whether two neighboring reads are in the same lfr')
	parser.add_argument('-thread', required=True, dest='thread', type=int, help='the number of threads')
	args = parser.parse_args()
	ref = args.reference
	read_1 = args.read_1
	read_2 = args.read_2
	length = args.insert_distant
	max_hit = args.max_hit
	lfr_thrshd = args.lfr_thrshd
	threads = args.thread

	log_setting("run_stLFR_SMA")
	cmd_info = " ".join(sys.argv)
	logger.info(cmd_info)

	create_work_floder()
	logger.info("Start to split read with the same barcode")
	max_workers = split_read_by_barcode(read_1,max_workers)
	logger.info("***End split")

	logger.info("Start to extract read by the name")
	extract_read_by_name(tool_path, read_1, read_2, max_workers, threads)
	logger.info("***End extract")

	logger.info("Start to slove multip alignment read")
	multip_run(tool_path, ref, max_workers, threads, max_hit, length, lfr_thrshd)
	logger.info("***End sloved")

	logger.info("merge the result")
	merge_bam_after_sloved(tool_path, threads)
	logger.info("***END merge")

	logger.info("delete the intermediate file")
	detele_intermediate_file()
	logger.info("The process successfully!")
