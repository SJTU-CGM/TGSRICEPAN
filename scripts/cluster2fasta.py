#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/3/10 2:50 PM
    @Usage: python3 cluster2fasta.py raw.fa glust.out output.fa
"""
# get the clusters' representative sequences from fasta and glust_cluster_output(X,cd-hit cover head of seq)

import sys
# input:
'''
>Cluster 0
0       5888230nt, >seq1... *
>Cluster 1
0       4800869nt, >seq2... *
>Cluster 2
0       3906592nt, >seq3... *
1       20nt, >seq4... at -/100.00%
'''

fasta_file = sys.argv[1]
# glust.out
clust_file = sys.argv[2]

#
outfa = sys.argv[3]
if fasta_file == outfa:
    exit()

representative_ctgs = dict()

i = 0
with open(clust_file) as f:
    for line in f:
        if line.startswith('>'):
            i += 1
        else:
            temp = line.rstrip().split()
            if temp[-1] == '*':
                representative = 1
            else:
                representative = 0
            if representative == 1:
                ctgname = temp[2].rstrip('.')
                representative_ctgs[ctgname] = ''

print("Representative number: " + str(i))

outFlag = 0
with open(outfa, 'w') as fout:
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if line.rstrip() in representative_ctgs:
                    outFlag = 1
                else:
                    outFlag = 0
            if outFlag == 1:
                fout.write(line)
