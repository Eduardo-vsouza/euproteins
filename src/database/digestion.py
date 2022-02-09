import os.path

import pandas as pd


class Digestor(object):
    def __init__(self, filetype, folder, threads):
        self.filetype = filetype
        self.folder = folder
        self.digestionFolder = f'{self.folder}/digestions'
        self.__check_folder()
        self.threads = threads

    def __check_folder(self):
        if not os.path.exists(self.digestionFolder):
            os.system(f'mkdir {self.digestionFolder}')

    def digest(self, enzyme=42):
        cmd = f'rpg -e {enzyme} --fmt tsv --inputdata {self.filetype}_database.fasta --outputfile ' \
              f'{self.digestionFolder}/db_digested --quiet --processes {self.threads}'
        os.system(cmd)

    def filter_by_size(self, min_size=6, max_size=40):
        df = pd.read_csv(f'{self.digestionFolder}/db_digested.tsv', sep='\t')
        df = df[df["Peptide_size"] >= min_size]
        df = df[df["Peptide_size"] <= max_size]
        peptides = df["Sequence"].tolist()
        orfs = df["Original_header"].tolist()
        peps_and_orfs = {}
        i = 0
        for pep, orf in zip(peptides, orfs):
            i += 1
            print('pep', i/len(peptides)*100, end='\r')
            if 'X' not in pep and '*' not in pep:
                if pep not in peps_and_orfs:
                    peps_and_orfs[pep] = [orf]
                else:
                    peps_and_orfs[pep].append(orf)

        to_write = [f'{pep}\n' for pep in peptides if 'X' not in pep and '*' not in pep]
        with open(f'{self.digestionFolder}/db_peptides.txt', 'w') as peps:
            peps.writelines(to_write)
        return peps_and_orfs