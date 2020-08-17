#!/usr/bin/env python3

# TO DO
## check if output (esp rev) is in frame

import pandas as pd
import subprocess
from Bio import SeqIO, Seq
import argparse
import os

parser = argparse.ArgumentParser(description='Extract genes from genes based on ABRicate output.')

parser.add_argument("-a", "--abricatefile", dest="abricatefile", help="ABRicate file to parse genes", metavar="ABRICATE FILE", type=str)
parser.add_argument("-g", "--genomedir", dest="genomedir", help="directory containing genomes", metavar="GENOMES DIR", type=str)
parser.add_argument("-o", "--output", dest="outdir", help="directory for output", metavar="OUTPUT DIR", type=str)
parser.add_argument("-s", "--suffix", dest="suffix", default=".fasta", help="Genome assembly file suffix (default: .fasta)")

args = parser.parse_args()

if not os.path.exists(args.outdir):
  os.makedirs(args.outdir)

df = pd.read_csv(args.abricatefile, sep = '\t')

for index, row in df.iterrows():
  strain = str(os.path.basename(row['#FILE'])).replace(args.suffix, '')
  genome = args.genomedir  + '/' + strain + args.suffix
  combination = strain + '_' + row['GENE']
  output = args.outdir + '/' + combination + '.out'
  if row['STRAND'] != '-':
    writestring = str(row['SEQUENCE']) + " " + str(row['START'] - 1) + " " + str(row['END'])
    f = open('tmp.txt', 'w')
    f.write(writestring)
    f.close()

    out = open(output, "w")
    subprocess.call(["seqtk", "subseq", genome, "tmp.txt"], stdout=out)
    out.close()
    record = SeqIO.read(output, "fasta")
    record.id = combination
    record.description = ''
    SeqIO.write(record, output, "fasta")

  else:
    writestring = str(row['SEQUENCE']) + " " + str(row['START'] - 1) + " " + str(row['END'])
    f = open('tmp.txt', 'w')
    f.write(writestring)
    f.close()

    rev = open("rev.fasta", "w")
    subprocess.call(["seqtk", "subseq", genome, "tmp.txt"], stdout = rev)
    rev.close()
    record = SeqIO.read("rev.fasta", "fasta").reverse_complement()
    record.id = combination
    record.description = ''
    SeqIO.write(record, output, "fasta")

