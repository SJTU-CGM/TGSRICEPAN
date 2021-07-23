#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/3/10 2:50 PM
    @Usage: python elongate_seq.py cut.fa elongated_length rawctg.fa elongated.fa > seq_long.table
"""
import sys


def union_interval(interval_list):
    union_i = []
    interval_list = sorted(interval_list)
    list_length = len(interval_list)
    if list_length == 0:
        return None
    elif list_length == 1:
        return [(interval_list[0][0], interval_list[0][1])]
    else:
        # >=2
        last_interval = interval_list[0]
        last_start = last_interval[0]
        last_end = last_interval[1]
        for i in range(1, list_length):
            current_interval = interval_list[i]
            current_start = current_interval[0]
            current_end = current_interval[1]
            if last_end < current_start:
                # xxx---
                # ----xx
                union_i += [(last_start, last_end)]
                last_start = max(last_start, current_start)
            else:
                last_start = min(last_start, current_start)
            last_end = max(last_end, current_end)
        union_i += [(last_start, last_end)]
        return union_i


def sub2parent_interval(pattern_interval, interval_list):
    for region in interval_list:
        if region[0] <= pattern_interval[0] and region[1] >= pattern_interval[1]:
            return region
    return tuple()


def subseq(seq, region):
    return seq[region[0] - 1:region[1]]


def batch_subseq(ctg, seq, regions):
    subseqs = dict()
    # regions may have same interval
    for region in set(regions):
        subseq_name = ctg + ':' + str(region[0]) + '-' + str(region[1])
        subseqs[subseq_name] = subseq(seq, region)
    return subseqs


# input
cutfa = sys.argv[1]
elongate_length = sys.argv[1]   #5000
rawfa = sys.argv[3]
elongatedfa = sys.argv[4]

# read raw ctg length
rawctg_lenth_dict = dict()
with open(rawfa) as rf:
    ctg = ''
    ctglen = 0
    for line in rf:
        if line.startswith('>'):
            # check
            if ctglen != 0:
                rawctg_lenth_dict[ctg] = ctglen
            ctg = line.rstrip()[1:]
            ctglen = 0
        else:
            ctglen += len(line.rstrip())
    # check last one
    if ctglen != 0:
        rawctg_lenth_dict[ctg] = ctglen
print("# Load raw seqs ok!")

# read cut ctg region
fulllength_seq_n = 0
partlength_seq_n = 0
cutctg_region_dict = dict()  # {'ctg1':[(1,2),(444,3333)]}
with open(cutfa) as cf:
    for line in cf:
        if line.startswith('>'):
            # >xxxxxxdddd:222-222222 or xxxxxxdddd
            keyword = "RagTag"
            temp = line.rstrip()[1:].split(keyword)
            region_str = temp[1].replace(':', '')
            sample_ctg = temp[0] + keyword
            if region_str == '':
                cutctg_region_dict[sample_ctg] = []
                fulllength_seq_n += 1
                continue
            try:
                region_start, region_end = int(region_str.split('-')[0]), int(region_str.split('-')[1])
            except Exception as e:
                print("# Error", temp)
                exit()
            if sample_ctg not in cutctg_region_dict:
                cutctg_region_dict[sample_ctg] = []
            cutctg_region_dict[sample_ctg] += [(region_start, region_end)]
            partlength_seq_n += 1
        else:
            continue
print("# Load {cutseq_n} cut seqs from {cutctg_n} ctgs, load {fulllength_ctg_n} full length seqs(ctgs)!".format(
    cutseq_n=str(partlength_seq_n),
    cutctg_n=str(len(cutctg_region_dict) - fulllength_seq_n),
    fulllength_ctg_n=str(fulllength_seq_n)))

# elongation
elongationctg_region_dict = dict()
elongationctg_length_dict = dict()
# count (only left/right)
left_elong_n = 0
right_elong_n = 0
both_elong_n = 0
for cutctg in cutctg_region_dict:
    if not cutctg_region_dict[cutctg]:
        continue
    for region in cutctg_region_dict[cutctg]:
        # init
        left_elongate_length = elongate_length
        right_elongate_length = elongate_length

        try:
            # left
            if region[0] - elongate_length <= 0:
                left_elongate_length = region[0] - 1
            # right
            if region[1] + elongate_length >= rawctg_lenth_dict[cutctg]:
                # 100 -1999
                right_elongate_length = rawctg_lenth_dict[cutctg] - region[1]
        except Exception as e:
            print("# Fail elongating at " + cutctg + str(region) + ' !')
            left_elongate_length = 0
            right_elongate_length = 0

        if left_elongate_length > 0 and right_elongate_length == 0:
            left_elong_n += 1
        elif left_elongate_length == 0 and right_elongate_length > 0:
            right_elong_n += 1
        elif left_elongate_length > 0 and right_elongate_length > 0:
            both_elong_n += 1

        if cutctg not in elongationctg_length_dict:
            elongationctg_length_dict[cutctg] = []
        if cutctg not in elongationctg_region_dict:
            elongationctg_region_dict[cutctg] = []
        elongationctg_length_dict[cutctg] += [(max(0, left_elongate_length), max(0, right_elongate_length))]
        elongationctg_region_dict[cutctg] += [(region[0] - left_elongate_length, region[1] + right_elongate_length)]

print(
    '''# Elongate {left_elong_n} cut seqs only left, elongate {right_elong_n} cut seqs only right, elongate {both_elong_n} cut seqs both sides!'''.format(
        left_elong_n=str(left_elong_n),
        right_elong_n=str(right_elong_n),
        both_elong_n=str(both_elong_n)))

# merge overlap
mergedctg_region_dict = dict()
for cutctg in cutctg_region_dict:
    if cutctg in elongationctg_region_dict:
        mergedctg_region_dict[cutctg] = union_interval(elongationctg_region_dict[cutctg])

# print all
print('\t'.join(['# RawSeq', 'SourceRegion', 'ElongatedLength(left,right)', 'ElongatedRegion', 'MergedRegion']))
for cutctg in cutctg_region_dict:
    if cutctg not in elongationctg_region_dict:
        # full unaligned
        if cutctg_region_dict[cutctg] == []:
            print('\t'.join([cutctg, '-', '-', '-', '-']))
        else:
            # part but not found parent ctg
            for i in range(len(cutctg_region_dict[cutctg])):
                fromregion = cutctg_region_dict[cutctg][i]
                print('\t'.join([str(x) for x in [cutctg, fromregion, (0, 0), fromregion, fromregion]]))
    else:
        elong_region_n = len(elongationctg_length_dict[cutctg])
        for i in range(elong_region_n):
            fromregion = cutctg_region_dict[cutctg][i]
            elen = elongationctg_length_dict[cutctg][i]
            eregion = elongationctg_region_dict[cutctg][i]
            mregion = sub2parent_interval(eregion, mergedctg_region_dict[cutctg])

            # ctg, source_region, elongated_length(left,right), elongated_region, merged_region
            print('\t'.join([str(x) for x in [cutctg, fromregion, elen, eregion, mregion]]))

# write fasta
success_write_part_n = 0
success_write_full_n = 0
ctg = ''
seq = ''
with open(elongatedfa, 'w') as ef:
    with open(rawfa) as rf:
        for line in rf:
            if line.startswith('>'):
                # check
                if seq != '':
                    if ctg in cutctg_region_dict:
                        # need write
                        # need part
                        if ctg in mergedctg_region_dict:
                            seqs_dict = batch_subseq(ctg, seq, mergedctg_region_dict[ctg])
                            for subregion_name in seqs_dict:
                                subregion_seq = seqs_dict[subregion_name]
                                if subregion_seq != '':
                                    ef.write('>' + subregion_name + '\n' + subregion_seq + '\n')
                                    success_write_part_n += 1
                        # need full
                        else:
                            ef.write('>' + ctg + '\n' + seq + '\n')
                            success_write_full_n += 1

                ctg = line.rstrip()[1:]
                seq = ''
            else:
                seq += line.rstrip()
        # last one
        if seq != '':
            if ctg in cutctg_region_dict:
                # need write
                # need part
                if ctg in mergedctg_region_dict:
                    seqs_dict = batch_subseq(ctg, seq, mergedctg_region_dict[ctg])
                    for subregion_name in seqs_dict:
                        subregion_seq = seqs_dict[subregion_name]
                        if subregion_seq != '':
                            ef.write('>' + subregion_name + '\n' + subregion_seq + '\n')
                            success_write_part_n += 1
                # need full
                else:
                    ef.write('>' + ctg + '\n' + seq + '\n')
                    success_write_full_n += 1
print("# Write {part} part seqs, {full} full seqs!".format(part=success_write_part_n,
                                                           full=success_write_full_n))
