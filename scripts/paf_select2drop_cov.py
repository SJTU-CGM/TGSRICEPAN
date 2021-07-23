import sys


def search_leaf(b_dict, name):
    visited = set()
    if not name in b_dict:
        return list(visited)
    else:
        sons = b_dict[name]
        unvisit = sons
        while unvisit != []:
            for son in unvisit:
                visited.add(son)
                if son in b_dict:
                    grandsons = b_dict[son]
                    for grandson in grandsons:
                        if not grandson in visited:
                            visited.add(grandson)
                            unvisit += [grandson]
                unvisit.remove(son)
        return visited


if __name__ == "__main__":
    try:
        func = sys.argv[1]
        paf = sys.argv[2]
        cutoff_cov_percent = float(sys.argv[3])
    except:
        print('''
        python paf_select.py command <paf_file>  <identity>

        command:
            drop - 
                paf file - output of tools like minimap2
                coverage - match coverage (0.8)      

        ''')

        exit()
    if func != 'drop':
        exit()

    target_dict = {}
    remove_set = set()

    with open(paf) as f:
        for line in f:
            # paf format
            # 1 string, query name
            # 2 int, query length
            # 3 int, query start (0 base,closed)
            # 4 int, query end (0 base,open)
            # 5 char, strand (+/-)
            # 6 string, target name
            # 7 int, target length
            # 8 int, target start
            # 9 int, target end
            # 10 int, number of residue matched
            # 11 int, alignment block length
            # 12 int, mapping quality (0~255,255:missing)
            temp = line.rstrip().split('\t')
            query_name = temp[0]
            target_name = temp[5]
            # self
            if query_name == target_name:
                continue

            query_length = float(temp[1])
            matchedlength = float(temp[9])
            # not self
            # less than cutoff
            if matchedlength / query_length < cutoff_cov_percent:
                continue

            # query_start = min(int(temp[2]),int(temp[3]))
            # query_end = max(int(temp[2]),int(temp[3]))

            target_length = float(temp[6])
            # target_start = min(int(temp[7]),int(temp[8]))
            # target_end = max(int(temp[7]), int(temp[8]))
            # identity = float(temp[9]) / float(temp[10])

            if query_length <= target_length:
                # no circle
                leaves = search_leaf(target_dict, query_name)
                if not target_name in leaves:
                    if not target_name in target_dict:
                        target_dict[target_name] = [query_name]
                    else:
                        target_dict[target_name] += [query_name]
                    remove_set.add(query_name)

        for query in remove_set:
            print(query)

