#!/usr/bin/env python3

import pandas as pd
import subprocess
from Bio import SeqIO, Seq
import argparse
import os
import re
import sys
import logging
import tempfile

def check_outdir(outdir):
  logging.debug(f"Calling function check_outdir")
  logging.debug(f"Checking whether {outdir} exists")
  if not os.path.exists(outdir):
    logging.debug(f"{outdir} does not exist, making this now")
    os.makedirs(outdir)

def read_abricatefile(abricatefile, csv):
  logging.debug(f"Calling function read_abricatefile")
  logging.debug(f"Checking whether ABRicate file should be read as a csv file: {csv}")
  if csv == True:
    logging.debug(f"Reading {abricatefile} as a csv file")
    df = pd.read_csv(abricatefile, sep = ',')
  else:
    logging.debug(f"Reading {abricatefile} as a tsv file")
    df = pd.read_csv(abricatefile, sep = '\t')
  return df

def check_combination(combination, combinations_passed):
  logging.debug(f"Calling function check_combination")
  found_duplicate = False
  logging.debug(f"Checking whether {combination} is in {combinations_passed}")
  if combination not in combinations_passed:
    logging.debug(f"{combination} was not found. Is added to checked combinations")
    combinations_passed[combination] = 1
    checked_combination = combination
  else:
    logging.debug(f"{combination} was found in checked combinations. Marked as duplicate")
    combinations_passed[combination] += 1
    combination_count = combinations_passed[combination]
    logging.debug(f"Joining {combination} with {combination_count}")
    checked_combination = '_'.join([combination, str(combination_count)])
    found_duplicate = True
  return checked_combination, combinations_passed, found_duplicate

def parse_row(row, combinations_passed, suffix, genomedir):
  logging.debug(f"Calling function parse_row")
  logging.debug(f"Constructing strain name")
  strain = str(os.path.basename(row['#FILE'])).replace(suffix, '')
  logging.debug(f"Genome name is constructed based on {genomedir} and {strain} and {suffix}")
  genome = genomedir  + '/' + strain + suffix
  combination = strain + '_' + row['GENE']
  logging.debug(f"Filtering {combination} using regex")
  combination = re.sub('[^\w\-_\. ]', '_', combination)
  logging.debug(f"Name is now {combination}")
  logging.debug(f"{combination} is checked against previously constructed combinations")
  checked_combination, combinations_passed, found_duplicate = check_combination(combination, combinations_passed)
  if found_duplicate:
    logging.warning(f"{combination} is found more than once in {strain}. Writing output sequence to {checked_combination}")
  logging.debug(f"Constructing output file name for {checked_combination}")
  output = args.outdir + '/' + checked_combination + '.out'
  return genome, checked_combination, output

def process_sense(row, genome, output):
  logging.debug(f"Calling function process_sense")
  START = row['START'] - 1
  END = row['END']
  writestring = ' '.join([str(row['SEQUENCE']), str(START), str(END)])
  logging.debug(f"{writestring} is written to NamedTemporaryFile")
  with tempfile.NamedTemporaryFile(mode='w+') as tf:
    tf.write(writestring)
    tf.seek(0)
    logging.debug(f"seqtk is called using subprocess for {genome}")
    with open(output, 'w+') as out:
      subprocess.call(["seqtk", "subseq", genome, tf.name], stdout=out)

  logging.debug(f"{output} is read using SeqIO")
  record = SeqIO.read(output, "fasta")

  return record

def process_antisense(row, genome, output):
  logging.debug(f"Calling function process_antisense")
  START = row['START'] - 1
  END = row['END']
  writestring = ' '.join([str(row['SEQUENCE']), str(START), str(END)])
  logging.debug(f"{writestring} is written to NamedTemporaryFile")
  with tempfile.NamedTemporaryFile(mode='w+') as tf:
    tf.write(writestring)
    tf.seek(0)
    logging.debug(f"seqtk is called using subprocess for {genome}")
    with tempfile.NamedTemporaryFile(mode='w+') as rev:
      subprocess.call(["seqtk", "subseq", genome, tf.name], stdout = rev)
      rev.seek(0)
      logging.debug(f"NamedTemporaryFile is read using SeqIO and reverse complemented")
      record = SeqIO.read(rev.name, "fasta").reverse_complement()

  return record

def parse_multiple_rows(df, suffix, genomedir):
  logging.debug(f"Calling function parse_multiple_rows")
  logging.debug(f"Assert whether only one FILE is found")
  original_file_name = df['#FILE'].unique()
  assert len(original_file_name) == 1
  logging.debug(f"Constructing strain name based on {original_file_name}")
  strain = str(os.path.basename(original_file_name[0]).replace(suffix, ''))
  logging.debug(f"Constructing genome name based on {genomedir} and {strain} and {suffix}")
  genome = genomedir  + '/' + strain + suffix
  logging.debug(f"Asser whether only one DATABASE is found")
  original_db = df['DATABASE'].unique()
  assert len(original_db) == 1
  combination = strain + '_' + original_db[0]
  logging.debug(f"Filtering {combination} using regex")
  combination = re.sub('[^\w\-_\. ]', '_', combination)
  logging.debug(f"Name is now {combination}")
  logging.debug(f"Constructing output file name for {combination}")
  output = args.outdir + '/' + combination + '.out'
  return genome, combination, output

def find_gene_boundary_extremes(df):
  logging.debug(f"Calling function find_gene_boundary_extremes")
  logging.debug(f"Finding gene boundaries based on whole file")
  list_extreme_gene_boundaries = [df['START'].min(), df['START'].max(), df['END'].min(), df['END'].max()]
  return min(list_extreme_gene_boundaries), max(list_extreme_gene_boundaries)

def update_record(record, combination):
  logging.debug(f"Calling function update_record")
  logging.debug(f"Updating record ID using {combination}")
  record.id = combination
  logging.debug(f"Removing description from record")
  record.description = ''
  return record

def main_genes(df, args):
  logging.debug(f"Calling function main_genes")
  combinations_passed = {}
  logging.debug(f"Looping through df")
  for index, row in df.iterrows():
    logging.debug(f"Parsing row {index}")
    genome, checked_combination, output = parse_row(row, combinations_passed, args.suffix, args.genomedir)
    logging.debug(f"Assert that STRAND column is found in abricatefile")
    assert 'STRAND' in row, "STRAND column is not found in ABRicate file, was ABRicate version 0.9.8 or later?"
    if row['STRAND'] != '-':
      logging.debug(f"Processing row with gene in sense")
      record = process_sense(row, genome, output)
    else:
      logging.debug(f"Processing row with gene in antisense")
      record = process_antisense(row, genome, output)
    logging.debug(f"Updating record with {checked_combination}")
    updated_record = update_record(record, checked_combination)
    logging.debug(f"Writing updated record to {output} as fasta file")
    SeqIO.write(updated_record, output, "fasta")

def main_genecluster(df, args):
  logging.debug(f"Calling function main_genecluster")
  logging.debug(f"Parse multiple rows at once for gene cluster processing")
  genome, combination, output = parse_multiple_rows(df, args.suffix, args.genomedir)
  logging.debug(f"Assert that STRAND column is found in abricatefile")
  assert 'STRAND' in df, "STRAND column is not found in ABRicate file, was ABRicate version 0.9.8 or later?"
  logging.debug(f"Checking the number of contigs on which hits are found")
  value_count_contigs = df['SEQUENCE'].value_counts()
  total_nr_contigs = len(value_count_contigs)
  total_nr_hits = value_count_contigs.sum()
  most_common_contig = value_count_contigs.index[0]
  logging.debug(f"Total number of contigs in file are checked")
  if total_nr_contigs > 1:
    logging.warning(f"Hits were found on {total_nr_contigs}. Contig {most_common_contig} has most hits and will be selected")
    if (value_count_contigs[0] / total_nr_hits) < 0.5:
      logging.warning(f"Contig {most_common_contig} has most hits, but contains less than half of all hits")
    elif value_count_contigs[0] == value_count_contigs[1]:
      logging.warning(f"The same number of hits were found on (at least) the top two contigs. Top contig is selected based on value_counts() sorting")
  df_most_common_contig = df[df['SEQUENCE'] == most_common_contig]
  logging.debug(f"Find gene boundary extremes")
  START, END = find_gene_boundary_extremes(df_most_common_contig)

  logging.debug(f"Construct pandas Series mimicking an object from df.iterrows")
  combined_row = pd.Series({'SEQUENCE': most_common_contig, 'START': START, 'END': END})

  logging.debug(f"Check how many hits are sense or antisense")
  strand_series = df_most_common_contig['STRAND'].value_counts()
  if ('+' in strand_series) and ('-' in strand_series):
    logging.debug(f"Sense and antisense hits in abricate output")
    if strand_series.loc['-'] > strand_series.loc['+']:
      logging.debug(f"More antisense hits are found")
      decision_strand = 'antisense'
    else:
      logging.debug(f"More sense hits are found")
      decision_strand = 'sense'
  elif '-' not in strand_series:
    logging.debug(f"No antisense hits in file: processing as sense")
    decision_strand = 'sense'
  elif '+' not in strand_series:
    logging.debug(f"No sense hits in file: processing as antisense")
    decision_strand = 'antisense'

  logging.debug(f"decision_strand has value {decision_strand}")
  if decision_strand == 'sense':
    logging.debug(f"Processing combined_row with gene cluster in sense")
    record = process_sense(combined_row, genome, output)
  elif decision_strand == 'antisense':
    logging.debug(f"Processing combined_row with gene cluster in antisense")
    record = process_antisense(combined_row, genome, output)
  else:
    logging.critical(f"decision_strand is {decision_strand}. This is an invalid value")
  logging.debug(f"Updating record with {combination}")
  updated_record = update_record(record, combination)
  logging.debug(f"Writing updated record to {output} as fasta file")
  SeqIO.write(updated_record, output, "fasta")

def main(args):
  logging.debug(f"Calling function main")
  check_outdir(args.outdir)
  df = read_abricatefile(args.abricatefile, args.csv)
  if args.genecluster == False:
    logging.debug(f"Executing main_genes for extraction of single genes")
    main_genes(df, args)
  else:
    logging.debug(f"Executing main_genecluster for extraction of gene clusters")
    main_genecluster(df, args)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Extract genes from genes based on ABRicate output.')
  
  parser.add_argument("-a", "--abricatefile", dest="abricatefile", help="ABRicate file to parse genes", metavar="ABRICATE FILE", required=True, type=str)
  parser.add_argument("-g", "--genomedir", dest="genomedir", help="directory containing genomes", metavar="GENOMES DIR", required=True, type=str)
  parser.add_argument("-o", "--output", dest="outdir", help="directory for output", metavar="OUTPUT DIR", required=True, type=str)
  parser.add_argument("-s", "--suffix", dest="suffix", default=".fasta", help="Genome assembly file suffix (default: .fasta)")
  parser.add_argument("--genecluster", dest="genecluster", action="store_true", help="Extract all genes to a single fasta if located on a single contig (default: false)")
  parser.add_argument("--csv", dest="csv", action="store_true", help="Use this option if your ABRicate output file is comma-separated (default: parse as tab-separated file).")
  parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="Increase verbosity")
  
  args = parser.parse_args()
  
  if not args.verbose:
    logging.basicConfig(level=logging.WARNING)
  elif args.verbose == 1:
    ## This can be set to INFO level if more levels are needed
    logging.basicConfig(level=logging.DEBUG)
  else:
    logging.basicConfig(level=logging.DEBUG)

  logging.debug(f"Executing main function")
  main(args)
