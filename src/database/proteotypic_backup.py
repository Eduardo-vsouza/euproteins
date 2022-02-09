import sys
import os
import multiprocessing as mp
import pandas as pd
from Bio import SeqIO
import ahocorasick


class DeepMSPeptide(object):
    def __init__(self, filetype, folder):
        self.executable = f'{sys.path[0]}/dependencies/DeepMSPeptide-master/DeepMSPeptide/DeepMSPeptide.py'
        self.filetype = filetype
        self.folder = folder
        self.digestionFolder = f'{self.folder}/digestions'

    def predict(self):
        peptides = f'{self.digestionFolder}/db_peptides.txt'
        cmd = f'python3 {self.executable} {peptides}'
        os.system(cmd)

    def filter_proteotypic(self, peps_and_orfs):
        filtered_database = []
        # df = pd.read_csv(f'{self.digestionFolder}/db_peptides_Predictions.txt', sep='\t')
        df = pd.read_csv(f'{self.digestionFolder}/db_peptides_Predictions.txt', sep='\t')

        df = df[df["Detectability"] == 1]
        ptps = list(set(df["Peptide"].tolist()))
        i = 0
        entries = []
        for ptp in ptps:
            i += 1
            print('ptps', i/len(ptps)*100, end='\r')
            if ptp in peps_and_orfs:
                orfs = peps_and_orfs[ptp]
                for orf in orfs:
                    entries.append(orf)
        with open(f'orfs_with_proteotypic.txt', 'w') as proteot:
            proteot.writelines([f'{entry}\n' for entry in entries])

    def remove_dupli(self):
        with open('orfs_with_proteotypic.txt', 'r') as prtps, open('orfs_no_dupli_proteot.txt', 'w') as out:
            lines = list(set(prtps.readlines()))
            out.writelines(lines)

    def filter_orf_db(self):
        entries = []
            # for line in lines:
            #     if line not in entries:
            #         entries.append(line.rstrip())
        #
        entries = []
        with open('orfs_no_dupli_proteot.txt', 'r') as handler:
            lines = handler.readlines()
            for line in lines:
                line = line.rstrip()
                entries.append(line)
        fa = []
        rcs = SeqIO.parse(f'{self.filetype}_database.fasta', 'fasta')
        recs = ['a' for rec in rcs]
        records = SeqIO.parse(f'{self.filetype}_database.fasta', 'fasta')
        i = 0
        total_recs = {}
        for record in records:
            # i += 1
            # print(i/len(recs)*100, end='\r')
            # if str(record.description) in entries:
            #     fa.append(record)
            total_recs[str(record.description)] = str(record.seq)
        def checker(entry, seq):
            if entry in entries:
                return f'>{entry}\n{seq}\n'
            else:
                return ''
        list_res = list(map(checker, list(total_recs.keys()), list(total_recs.values())))
        with open('filtered_transcriptome_database.fasta', 'w') as out:
            out.writelines(list_res)
        # SeqIO.write(fa, 'filtered_transcriptome_database.fasta', 'fasta')
