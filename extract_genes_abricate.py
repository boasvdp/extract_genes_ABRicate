#!/usr/bin/env python3

import pandas as pd
import subprocess
from Bio import SeqIO, Seq
import argparse
import os
import re
import sys

def check_outdir(outdir):
  if not os.path.exists(outdir):
    os.makedirs(outdir)

def read_abricatefile(abricatefile, csv):
  if csv == True:
    df = pd.read_csv(abricatefile, sep = ',')
  else:
    df = pd.read_csv(abricatefile, sep = '\t')
  return df

def check_combination(combination, combinations_passed):
  found_duplicate = False
  if combination not in combinations_passed:
    combinations_passed[combination] = 1
    checked_combination = combination
  else:
    combinations_passed[combination] += 1
    combination_count = combinations_passed[combination]
    checked_combination = '_'.join([combination, str(combination_count)])
    found_duplicate = True
  return checked_combination, combinations_passed, found_duplicate

def parse_row(row, combinations_passed, suffix, genomedir):
  strain = str(os.path.basename(row['#FILE'])).replace(suffix, '')
  genome = genomedir  + '/' + strain + suffix
  combination = strain + '_' + row['GENE']
  combination = re.sub('[^\w\-_\. ]', '_', combination)
  checked_combination, combinations_passed, found_duplicate = check_combination(combination, combinations_passed)
  if found_duplicate:
    printstring = ''.join(['INFO: ', combination, ' is found more than once in ', strain, '. Writing output sequence to ', checked_combination])
    print(printstring, file=sys.stderr)
  output = args.outdir + '/' + checked_combination + '.out'
  return genome, checked_combination, output

def process_sense(row, genome, output):
  writestring = str(row['SEQUENCE']) + " " + str(row['START'] - 1) + " " + str(row['END'])
  f = open('tmp.txt', 'w')
  f.write(writestring)
  f.close()

  out = open(output, "w")
  subprocess.call(["seqtk", "subseq", genome, "tmp.txt"], stdout=out)
  out.close()
  record = SeqIO.read(output, "fasta")

  return record

def process_antisense(row, genome, output):
  writestring = str(row['SEQUENCE']) + " " + str(row['START'] - 1) + " " + str(row['END'])
  f = open('tmp.txt', 'w')
  f.write(writestring)
  f.close()

  rev = open("rev.fasta", "w")
  subprocess.call(["seqtk", "subseq", genome, "tmp.txt"], stdout = rev)
  rev.close()
  record = SeqIO.read("rev.fasta", "fasta").reverse_complement()

  return record

def update_record(record, combination):
  record.id = combination
  record.description = ''
  return record

def cleanup():
  for file in ['rev.fasta', 'tmp.txt']:
    if os.path.isfile(file):
      os.remove(file)

def main(args):
  check_outdir(args.outdir)
  df = read_abricatefile(args.abricatefile, args.csv)
  cleanup()
  combinations_passed = {}
  for index, row in df.iterrows():
    genome, checked_combination, output = parse_row(row, combinations_passed, args.suffix, args.genomedir)
    if row['STRAND'] != '-':
      record = process_sense(row, genome, output)
    else:
      record = process_antisense(row, genome, output)
    updated_record = update_record(record, checked_combination)
    SeqIO.write(updated_record, output, "fasta")
  cleanup()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract genes from genes based on ABRicate output.')
  
  parser.add_argument("-a", "--abricatefile", dest="abricatefile", help="ABRicate file to parse genes", metavar="ABRICATE FILE", required=True, type=str)
  parser.add_argument("-g", "--genomedir", dest="genomedir", help="directory containing genomes", metavar="GENOMES DIR", required=True, type=str)
  parser.add_argument("-o", "--output", dest="outdir", help="directory for output", metavar="OUTPUT DIR", required=True, type=str)
  parser.add_argument("-s", "--suffix", dest="suffix", default=".fasta", help="Genome assembly file suffix (default: .fasta)")
  parser.add_argument("--csv", dest="csv", action="store_true", help="Use this option if your ABRicate output file is comma-separated (default: parse as tab-separated file).")
  
  args = parser.parse_args()
  
  main(args)
