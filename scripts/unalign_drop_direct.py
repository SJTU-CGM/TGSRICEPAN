import sys

def trans_contig_name(stringA):
    return stringA.replace(':','_').replace('+','_').replace('-','_').replace(')','_').replace('(','_')


def string2interval(stringB,cutoff,sep1=',',sep2='-'):
    coords = list([coord for coord in stringB.split(sep1)])
    coords = list([tuple([int(y) for y in x.split(sep2)]) for x in coords])
    for i in range(len(coords)):
        if coords[i][0] > coords[i][1]:
            coords[i] = tuple((coords[i][1], coords[i][0]))
        if len(cutoff) == 2:
            if (coords[i][1] - coords[i][0] < cutoff[0]) or (coords[i][1] - coords[i][0] > cutoff[1]):
                # not pass
                coords[i] = tuple((-1, -1))
        else:
            if coords[i][1] - coords[i][0] < cutoff[0]:
                # not pass
                coords[i] = tuple((-1, -1))
    return coords


try:
    fasta_file = sys.argv[1]
    unalign_table = sys.argv[2]
    output_fasta = sys.argv[3]

except:
    print('''
    python3 unalign_drop_direct.py <assembly_fasta> <unalign_info> <output_fasta>  [length_cufoff] [add_sample_title] 
    
    Filter unalign scaffolds/contigs, including fully and partially 

Need:
    <assembly_fasta> - assembly of contig or scaffold
    <unalign_info> - contig_report_xxx.unaligned.info generated from quast in contigs_reports
    <output_fasta> - write the unaligned region sequencea
    
Optional:    
    [length_cutoff] int - the min length to output unaligned regions from assembly (default=500)
                     int:int - an interval of length also is ok
    [add_sample_title] string - add sample name before each contig/scaffold
    
''')
    exit()

try:
    if len(sys.argv[4].split(':')) == 2:
        length_cutoff = list([int(x) for x in sys.argv[4].split(':')])
    else:
        length_cutoff = [int(sys.argv[4])]
except:
    length_cutoff = 500
print('# unalign length at least ' + str(length_cutoff))

try:
    add_sample_title = sys.argv[5]
    print('# add sample title = ' + add_sample_title)
except:
    add_sample_title = ''



unalign_dict= {}
with open(unalign_table) as f:
    # contig: total_length unaligned_length type unaligned_parts
    for line in f:
        temp = line.rstrip().split('\t')
        if temp[3] == 'none':
            unalign_dict[trans_contig_name(temp[0].rstrip())] = tuple((temp[1], temp[2], temp[3], '-'))
        else:
            unalign_dict[trans_contig_name(temp[0].rstrip())] = tuple((temp[1],temp[2],temp[3],temp[4]))

# output
# OUTPUT = 0
# init
input_record_num = 0
output_record_num = 0
partial_queryname = ''
partial_queryseq = ''
output_interval = []

with open(output_fasta,'w') as fout:
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                input_record_num += 1
                # check partial record
                if partial_queryname != '':
                    for interval in output_interval:
                        if interval[0] >= 0:
                            output_record_num += 1
                            fout.write('>'+add_sample_title+'_'+partial_queryname+':'+str(interval[0])+'-'+str(interval[1])+'\n')
                            fout.write(partial_queryseq[interval[0]-1:interval[1]]+'\n')

                # init
                tempquery = line.rstrip()[1:]
                partial_queryname = ''
                partial_queryseq = ''
                # full output
                FULL_OUTPUT_FLAG = 0


                if trans_contig_name(tempquery) in unalign_dict:
                    #if unalign_dict[trans_contig_name(tempquery)][2] == 'full':
                    #    # full
                    #    OUTPUT = 1
                    #else:
                        # partial
                    if unalign_dict[trans_contig_name(tempquery)][2] == 'none':
                        continue
                    partial_queryname = tempquery
                    output_interval = string2interval(unalign_dict[trans_contig_name(tempquery)][3],length_cutoff)
                # no mapped region in blast
                else:
                    FULL_OUTPUT_FLAG = 1
                    output_record_num += 1
            else:
                if partial_queryname != '':
                    partial_queryseq += line.rstrip()


            if FULL_OUTPUT_FLAG == 1:
                #fout.write(line)
                if line.startswith('>'):
                    fout.write('>'+add_sample_title+'_'+line[1:])
                else:
                    fout.write(line)

        # last record
        if partial_queryname != '':
            for interval in output_interval:
                if interval[0] >= 0:
                    output_record_num += 1
                    fout.write('>' +add_sample_title+'_' + partial_queryname + ':' + str(interval[0]) + '-' + str(interval[1]) + '\n')
                    # [1,maxlength]
                    fout.write(partial_queryseq[interval[0]-1:interval[1]]+'\n')

print('# Input: {in_r} records, Output: {out_r} records'.format(in_r=input_record_num,out_r=output_record_num))