stLFR_SMA (solve multiple alignment)
====================================

Install stLRF_SMA
------------------------------------
```
$git clone https://github.com/sumoii/stLFR_SMA.git
```
Dependencies
------------------------------------
1.python3.6+
2.gcc 4.4+

Usage
-----------------------------------
```
usage: solve_multip_align.py [-h] -ref REFERENCE -1 READ_1 -2 READ_2 -mhit
                             MAX_HIT -length INSERT_DISTANT -lfr_thrshd
                             LFR_THRSHD -thread THREAD

Sloved the mutilp alignment reads

optional arguments:
  -h, --help            show this help message and exit
  -ref REFERENCE        the reference file
  -1 READ_1             the paired_end of read1 must be .gz
  -2 READ_2             the paired_end of read2 must be .gz
  -mhit MAX_HIT         The max hit numbers of a read pos
  -length INSERT_DISTANT
                        The max length of PE insert distant can be accepted
  -lfr_thrshd LFR_THRSHD
                        The distance threshold is used to determine whether two
                        neighboring reads are in the same lfr
  -thread THREAD        the number of threads
```

Examples
----------------------------------
an example of simulation reference data and the paired end stLFR sequencing data (fastq)
```
python3 /dellfsqd2/ST_OCEAN/USER/xiaogaohong/software/stLFR_SMA/solve_multip_align.py -ref /inputs/chr8/chr8_diploid.fa -1 /inputs/chr8/chr8_400k.r1.fq.gz -2 /inputs/chr8/chr8_400k.r2.fq.gz -mhit 10 -length 2000 -lfr_thrshd 100000 -thread 10
```

Output
----------------------------------
```
.
├── name_list.txt
├── r1_pos_info
├── r2_pos_info
├── run_stLFR_SMA_log.txt
├── Sloved_Bar_info
└── stLFR_SMA.bam
```
