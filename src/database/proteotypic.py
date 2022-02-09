import sys
import os
import multiprocessing as mp
import pandas as pd
from Bio import SeqIO





manager = mp.Manager()
entries = []
list_res = manager.list()


def checker(protein):
    if protein[0] in entries:
        list_res.append(f'>{protein[0]}\n{protein[1]}\n')

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

    def filter_orf_db(self, threads):



            # for line in lines:
            #     if line not in entries:
            #         entries.append(line.rstrip())
        #
        # entries = []
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

        total_recs = []
        for record in records:
            # i += 1
            # print(i/len(recs)*100, end='\r')
            # if str(record.description) in entries:
            #     fa.append(record)
            total_recs.append((str(record.description), str(record.seq)))

            # else:
            #     return ''
        jobs = []
        p = mp.Pool(threads)
        p.map(checker, total_recs)
        # p.start()
        # p.join()
        # print(p)
        # jobs.append(p)
        # p.start()
        # for proc in jobs:
        #     proc.join()
        # p.join()
        # return_dict = manager

        # print(return_dict.values())
        # list_res = list(map(checker, list(total_recs.keys()), list(total_recs.values())))
        # list_res = p.starmap(checker, total_recs)

        with open('filtered_transcriptome_database.fasta', 'w') as out:

            out.writelines(list_res)
        # SeqIO.write(fa, 'filtered_transcriptome_database.fasta', 'fasta')

    def readd_annotation(self, proteome):
        records = SeqIO.parse('filtered_transcriptome_database.fasta', 'fasta')
        seqs = [str(record.seq) for record in records]
        proteome = SeqIO.parse(proteome, 'fasta')
        recs = SeqIO.parse('filtered_transcriptome_database.fasta', 'fasta')

        for record in proteome:
            if str(record.seq) not in seqs:
                recs.append(record)
        SeqIO.write(format='fasta', handle='transcriptome_database_filtered_anno_readded.fasta', sequences=recs)
