#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/17 3:49 PM
    @Usage:
"""
import sys


def add_newline(seq, maxchar_in_one_line=60, endnewline=True):
    newseq = ''
    seqlen = len(seq)
    for i in range(0, seqlen, maxchar_in_one_line):
        newseq += seq[i:min(i + maxchar_in_one_line, seqlen)] + '\n'
    if endnewline:
        return newseq
    else:
        return newseq.rstrip('\n')

def filter_gapfill(update_file, ctg_file, method='max', min_gap_len=1, max_gap_len=1000000):
    # updated_scaff_infos
    gap_seq = dict()
    chr_head = dict()
    chr_first_ctg_id = dict()
    with open(update_file) as f:
        for line in f:
            if line.startswith('>'):
                continue
            temp = line.rstrip().split('\t')
            if len(temp) < 7:
                continue
            ctg_id, strand, fill_len, ctg_len, chr_start, ctg_chr_id, chr_id = temp[:7]
            if len(temp) <= 7:
                if ctg_id not in gap_seq:
                    gap_seq[ctg_id] = 'N' * int(fill_len)
                if ctg_chr_id == '1':
                    chr_head[chr_id] = ''
                    chr_first_ctg_id[ctg_id] = chr_id
            else:
                if ctg_chr_id == '1':
                    chr_head[chr_id] = ''
                    chr_first_ctg_id[ctg_id] = chr_id
                    if 'PREV_N=' in line:
                        gap_seq[ctg_id] = 'N' * int(fill_len)
                        head_n = line.rstrip().split('PREV_N=')[1].split('\t')[0]
                        chr_head[chr_id] = 'N' * int(head_n)

                if '_FILL=' in line:
                    seq = line.rstrip().split('_FILL=')[1].split('\t')[0]
                    if ctg_id not in gap_seq:
                        gap_seq[ctg_id] = seq
                    else:
                        if gap_seq[ctg_id].startswith('N') and len(gap_seq[ctg_id]) >= min_gap_len:
                            gap_seq[ctg_id] = seq
                        elif method == 'max' and len(seq) > len(gap_seq[ctg_id]) and len(seq) < max_gap_len:
                            gap_seq[ctg_id] = seq
                        elif method == 'min' and (len(seq) < len(gap_seq[ctg_id]) or len(gap_seq[ctg_id]) == 0):
                            gap_seq[ctg_id] = seq
    #for ctg_id in gap_seq:
    #    print(ctg_id, str(len(gap_seq[ctg_id])))
    with open(ctg_file) as fa:
        chr_n = 0
        current_seq = ''
        tail_seq = ''
        for line in fa:
            if line.startswith('>'):
                ctg_id = line.rstrip()[1:].split()[0]
                current_seq = current_seq + tail_seq

                tail_seq = gap_seq[ctg_id]
                if ctg_id in chr_first_ctg_id:
                    if current_seq != '':
                        print('>Chr' + str(chr_n))
                        print(add_newline(current_seq, maxchar_in_one_line=100, endnewline=False))
                    #print(ctg_id)
                    chr_n += 1
                    current_seq = chr_head[str(chr_n)]
            else:
                current_seq += line.rstrip()
        # final
        current_seq = current_seq + tail_seq
        print('>Chr' + str(chr_n))
        print(add_newline(current_seq, maxchar_in_one_line=100, endnewline=False))





filter_gapfill(sys.argv[1],sys.argv[2],method=sys.argv[3])