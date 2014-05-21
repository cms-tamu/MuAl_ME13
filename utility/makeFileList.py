import os
import commands
import sys

if(len(sys.argv) < 2):
    print "usage: python makeFileList.py [path to eos folder]"
    sys.exit()

lines = []

for dataset in sys.argv[1:]:
    status,output = commands.getstatusoutput("cmsLs %s | awk \'{print $5}\'" % (dataset))

    lines += [line for line in output.split("\n") if len(line) > 2]


print "fileNames = ["
for i,line in enumerate(lines):
    line = line.strip()

    if(i == len(lines)-1):
        print "'"+line+"'"
    else:
        print "'"+line+"',"

print "]"

