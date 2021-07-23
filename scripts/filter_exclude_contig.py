# python3 filter_exclude_contig.py blast.out taxfile cutoff_pident cutoff_coverage cutoff_length > exclude_contig.txt
# exclude_contig.txt: the contigs not in white list
# white list: Viridiplantae


import sys


ctg_dict = {}

if __name__ == "__main__":
    try:
        blastout = sys.argv[1]
        taxfile = sys.argv[2]
        cutoff_pident = float(sys.argv[3]) # 0
        cutoff_coverage = float(sys.argv[4]) # 0
        cutoff_length = float(sys.argv[5]) # 0
        # 10 viridiplantae, 5 oryza
        used_col = 10
        white_list = 'Viridiplantae'
        #used_col = 5
        #white_list = 'Oryza'
    except:
        pass
        exit()

    acc_tax = dict()
    with open(taxfile) as f:
        for line in f:
            temp = line.rstrip().split('|')
            acc_tax[temp[0]] = temp[1:]

    evalue_idx = 2
    pident_idx = 3
    length_idx = 4

    queryname_idx = 0
    querystart_idx = 5
    queryend_idx = 6
    refname_idx = 1
    refstart_idx = 7
    refend_idx = 8


    exclude_ctg = []

with open(blastout) as f:
    for line in f:
        temp = line.rstrip().split('\t')
        queryname = temp[queryname_idx]
        pident = float(temp[pident_idx])
        refname = temp[refname_idx]
        length = float(temp[length_idx])
        if pident < cutoff_pident or length < cutoff_length :
            continue
        try:
            if acc_tax[refname][used_col - 2] != white_list:
                exclude_ctg += [queryname]
        except:
            pass


for x in set(exclude_ctg):
    print(x)


