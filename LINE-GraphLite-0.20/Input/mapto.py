import sys
import os

fi = open("../experiment/net_youtube_dense.txt", 'r')
fo = open("net_youtube_dense_t", 'w')
fo.write("336934\n")   # number of vertices
fo.write("46598619\n") # number of edges
dic = dict()
count = 0
for line in fi:
    items = line.strip().split()
    if items[0] in dic:
        items0 = dic.get(items[0])
    else:
        items0 = count
        dic[items[0]] = count
        count += 1
    if items[1] in dic:
        items1 = dic.get(items[1])
    else:
        items1 = count
        dic[items[1]]= count
        count += 1
    #fo.write('{0} {1} {2}\n'.format(items0, items1, items[2]))
fi.close()
fo.close()

with open("mapdict", 'w') as f:
    f.write(str(dic))
