import os,sys

for i in range(len(sys.argv)):
    "sys.argv[%d] = %s" % (i, sys.argv[i])

fin = open(sys.argv[1], "r")
fout = open(sys.argv[2], "w")

count = 0
for line in fin:
    line = line.rstrip(os.linesep)
    s = line.split("\t")
    if(s[10]=="HERVL-int"):
        fout.write(s[5]+"\t"+s[6]+"\t"+s[7]+"\t"+s[9]+"\n")





fin.close()
fout.close()

