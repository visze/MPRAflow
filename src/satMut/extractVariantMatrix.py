#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *18.07.2016
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import gzip

def selColumnsToVector(selColumns,lenColumns):
  selColumns.sort()
  last = 0
  res = ""
  for elem in selColumns:
    res += (elem-last)*"0\t"+"1\t"
    last=elem+1
  res += (lenColumns-last)*"0\t"
  return res.rstrip()
    
parser = OptionParser("%prog [options] filename")
parser.add_option("-a","--assignment", dest="assignment", help="Assignment of variants and tags (def '../assignment/LDLR.variants.txt.gz')",default="../assignment/LDLR.variants.txt.gz")
parser.add_option("-i","--InDels", dest="InDels", help="Keep barcodes with InDels (def Off)",default=False,action="store_true")
parser.add_option("-e","--IgnoreInDels", dest="IgnoreInDels", help="Ignore InDel annotations for barcodes (def Off)",default=False,action="store_true")
parser.add_option("-v","--verbose", dest="verbose", help="Turn on verbose output (def Off)",default=False,action="store_true")
(options, args) = parser.parse_args()

if not os.path.exists(options.assignment):
  sys.stderr.write("Error: Assignment file does not exist\n")
  sys.exit()

if len(args) != 1:
  sys.stderr.write("Error: Only one count file allowed as input.\n")
  sys.exit()
else:
  filename = args[0]
  
if not os.path.exists(filename):
  sys.stderr.write("Error: Input file does not exist.\n")
  sys.exit()
  
variants = defaultdict(set)
allVariants = set()
poolCounts = {}
isWT = set()

lPrefix = None
infile = gzip.open(options.assignment) if options.assignment.endswith(".gz") else open(options.assignment)
for line in infile:
  if line.startswith('#'): continue
  fields = line.rstrip().split()
  haveIndel = False
  if not options.InDels and not options.IgnoreInDels:
    for elem in fields[1:]:
      var=elem.split(":")[-1]
      if (len(var) != 3) or ("." in var):
        haveIndel = True
        break
  if haveIndel: continue
  if len(fields) == 1:
    isWT.add(fields[0])
    #fields.append("WT:0:X>X")
  for elem in fields[1:]:
    if options.IgnoreInDels:
      var=elem.split(":")[-1]
      if (len(var) != 3) or ("." in var): continue
    fkey = elem.split(":")
    if lPrefix == None: lPrefix = len(fkey[0])
    elif len(fkey[0]) < lPrefix: lPrefix = len(fkey[0])
    allVariants.add((int(fkey[1]),elem))
    variants[fields[0]].add(elem)
  tags = fields[0].split(',')
  if len(tags) > 1:
    for tag in tags:
      poolCounts[tag]=fields[0]
infile.close()

if options.verbose: sys.stderr.write("Read %d variant to barcode assignments\n"%(len(variants)))

poolCountVal1 = defaultdict(int)
poolCountVal2 = defaultdict(int)

infile = gzip.open(filename) if filename.endswith(".gz") else open(filename)
for line in infile:
  if line.startswith('#'): continue
  fields = line.rstrip().split()
  if len(fields) > 2:
    if fields[0] in poolCounts:
      val1 = int(fields[1])
      val2 = int(fields[2])
      poolCountVal1[poolCounts[fields[0]]]+=val1
      poolCountVal2[poolCounts[fields[0]]]+=val2
infile.close()

if options.verbose: sys.stderr.write("Read %d variants\n"%(len(allVariants)))

lPrefix+=1
columnAssignment = {}
header = "#Barcode\tDNA\tRNA\t"
for ind,(pos,variant) in enumerate(sorted(allVariants)):
  columnAssignment[variant]=ind
  header += variant[lPrefix:].replace(":","_").replace(">",".")+"\t"
header = header.rstrip()
sys.stdout.write(header+"\n")

totalVarColumns = len(allVariants)
del allVariants

infile = gzip.open(filename) if filename.endswith(".gz") else open(filename)
for line in infile:
  if line.startswith('#'): continue
  fields = line.rstrip().split()
  if len(fields) > 2 and ((fields[0] in variants) or (fields[0] in isWT)):
    if fields[0] in poolCounts:
      pname = poolCounts[fields[0]]
      selColumns = []
      for varName in variants[pname]:
        if varName in columnAssignment: selColumns.append(columnAssignment[varName])
      sys.stdout.write("%s\t%s\t%s\t%s\n"%(pname,poolCountVal1[pname],poolCountVal1[pname],selColumnsToVector(selColumns,totalVarColumns)))
    else:
      selColumns = []
      for varName in variants[fields[0]]:
        if varName in columnAssignment: selColumns.append(columnAssignment[varName])
      sys.stdout.write("%s\t%s\t%s\t%s\n"%(fields[0],fields[1],fields[2],selColumnsToVector(selColumns,totalVarColumns)))
infile.close()
