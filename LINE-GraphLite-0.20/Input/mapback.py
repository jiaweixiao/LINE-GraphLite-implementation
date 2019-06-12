import sys
import os
fname = "../Output_order1_thread40_e5/"
fi = open(fname + "out_1", 'r')
fo = open("../experiment/vec_1nd_wo_norm.txt", 'wb')
f = open("mapdict",'r')
dic = eval(f.read())
rdic = dict(zip([str(x) for x in dic.values()], dic.keys()))
for line in fi:
    items = line.strip().split(' ', 1 )
    if items[0] in rdic:
        fo.write('{0} {1}\n'.format(rdic[items[0]], items[1]))
    else:
        print items[0]
        print(" dict error!\n")
        exit(1)
fi.close()
fo.close()
