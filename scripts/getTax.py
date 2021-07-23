# python3 getTax.py <accession_file> <accession_col> <acc2taxid> <ranklineage> > tax.txt
# OUTPUT: tax.txt
# accession     taxid   linkage_information

import sys


if __name__ == "__main__":
    try:
        input_acc_file = sys.argv[1]
        acc_col = int(sys.argv[2])
        acc2taxid_file = sys.argv[3]
        ranklineage_file = sys.argv[4]

    except:
        print('''
        python3 getTax.py <accession_file> <accession_col> <acc2taxid> <ranklineage>
    
        <accession_file> - output such as blast
        <accession_col>  - the colnum of accession (e.g. blast: 2)
        <acc2taxid> - accesiontotax file
        <ranklineage> - rank lingeage file
        
        
        ''')
        exit()

    # used acc add
    used_accession = dict()
    with open(input_acc_file) as f:
        for line in f:
            temp = line.rstrip().split()
            #print(temp)
            acc = temp[acc_col - 1].split('.')[0]
            if not acc in used_accession:
                # ['taxid','']
                used_accession[acc] = ['','']

    # acc to tax
    used_tax = dict()
    with open(acc2taxid_file) as f:
        taxidx = 2
        if acc2taxid_file.endswith('FULL'):
            taxidx = 1
        for line in f:
            temp = line.rstrip().split()
            acc = temp[0].split('.')[0]
            if acc in used_accession:
                used_accession[acc][0] = temp[taxidx]
                if not temp[taxidx] in used_tax:
                    used_tax[temp[taxidx]] = [acc]
                else:
                    used_tax[temp[taxidx]] += [acc]


    # tax to linkage
    with open(ranklineage_file) as f:
        for line in f:
            temp = [x.strip() for x in  line.split('|')]
            if temp[0] in used_tax:
                for acc in used_tax[temp[0]]:
                    used_accession[acc][1] = '|'.join(temp[1:])

    for key in used_accession:
        print(key+"|"+used_accession[key][0]+"|"+used_accession[key][1])


