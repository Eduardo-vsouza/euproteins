import sys
import os

import pandas as pd
from Bio import SeqIO


class SieveFiltering(object):
    def __init__(self, filetype, folder):
        """ Filters the database excluding protein sequences that do not contain a proteotypic peptide after
        running Peptide Sieve. """
        self.filetype = filetype
        self.folder = folder
        self.dependencies = f'{sys.path[0]}/dependencies/peptide_sieve'
        self.sievePath = f'{self.dependencies}/PeptideSieve.linux.i386'
        self.propertiesFile = f'{self.dependencies}/properties.txt'

    def run_peptide_sieve(self, probability_cutoff=0.8, miss_cleavages=2):
        cmd = f'{self.sievePath} -O {self.folder}/PeptideSieve -p {probability_cutoff} -P {self.propertiesFile} -c {miss_cleavages} ' \
              f'-s {self.filetype}_database.fasta'
        os.system(cmd)

    def __add_header(self):
        with open(f'{self.folder}/PeptideSieve/{self.filetype}_database.ptps.out', 'r') as handler, \
                open(f'{self.folder}/PeptideSieve/{self.filetype}_database_ptps_fixed.txt', 'w') as out:
            to_write = ['Experiment\tProtein\tPeptide\tProbability_score\n']
            lines = handler.readlines()
            to_write = to_write + lines
            # for line in lines:
            #     to_write.append(line)
            out.writelines(to_write)

    def filter_database(self):
        self.__add_header()
        df = pd.read_csv(f'{self.folder}/PeptideSieve/{self.filetype}_database_ptps_fixed.txt', sep='\t')
        proteins = df["Protein"].tolist()
        records = SeqIO.parse(f'{self.filetype}_database.fasta', 'fasta')
        fasta = []
        for record in records:
            if str(record.description) in proteins:
                fasta.append(f'>{str(record.description)}\n{str(record.seq)}\n')
        with open(f'{self.filetype}_database_sieve_filtered.fasta', 'w') as out:
            out.writelines(fasta)