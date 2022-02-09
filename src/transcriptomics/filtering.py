import pandas as pd
import os

from ..sequtils import StringTieGFF


class ExpressionFilter(object):
    def __init__(self, gtf_intermediate_folder, assembled_gtf):
        self.intermediatesFolder = gtf_intermediate_folder
        self.intermediateFiles = os.listdir(gtf_intermediate_folder)
        self.filteredGenes = []
        self.assembled = assembled_gtf

    def __behead(self, gtf, file):
        with open(gtf, 'r') as handler, open(f'{self.intermediatesFolder}/{file[:-4]}_filtered.gtf', 'r') as filtered, \
                open(f'{self.intermediatesFolder}/{file[:-4]}_reheaded.gtf', 'w') as out:

            head = handler.readlines()[0]
            to_write = []
            to_write.append(head)
            to_write.append('# StringTie version 2.2.0\n')

            # to_write.append(head)
            filtered_lines = filtered.readlines()[1:]
            for line in filtered_lines:
                to_write.append(line)
            # total = head+filtered.readlines()[1:]
            # with open(f'{self.intermediatesFolder}/{file[:-4]}_filtered.gtf', 'w') as out:
            out.writelines(to_write)

    def filter_intermeds_before_assembly(self, cutoff):
        i = 0
        for file in self.intermediateFiles:
            if 'filtered' not in file and 'rehead' not in file:
                i += 1
                gtf = StringTieGFF(f'{self.intermediatesFolder}/{file}')
                df_init = gtf.get_gtf()
                df = df_init[df_init["feature"] == "transcript"]
                attrs = df["attributes"].tolist()
                tpms = []
                genes = []
                for a in attrs:
                    splat = a.split(";")
                    gene = splat[0].replace("\"", "").replace("gene_id ", "").replace("\"", "")
                    index = [idx for idx, s in enumerate(splat) if 'TPM ' in s and idx >= 4][0]
                    tpm = splat[index].replace("TPM", "").replace("\"", "").replace(" ", "").replace("\"", "")
                    tpms.append(float(tpm))
                    genes.append(gene)
                df.insert(5, "TPM", tpms)
                df.insert(6, "gene", genes)
                df = df[df["TPM"] > cutoff]
                df = df_init[df_init["attributes"].str.contains('|'.join(df["gene"].tolist())).groupby(level=0).any()]
                df.to_csv(f'{self.intermediatesFolder}/{file[:-4]}_filtered.gtf', sep='\t', index=False)
                self.__behead(f'{self.intermediatesFolder}/{file}', file)



    def move_outside(self):
        """ Moves the old intermediante GTF files to the other folder, so the filtered ones can be merged.
        """
        path = f"{self.intermediatesFolder}/unfiltered_gtf"
        if not os.path.exists(path):
            os.system(f'mkdir {path}')
        new = os.listdir(self.intermediatesFolder)
        for file in new:
            if 'rehead' not in file:
                if f'{self.intermediatesFolder}/{file}' not in os.listdir(path) and f'{self.intermediatesFolder}/{file}' != path:
                    cmd = f'mv {self.intermediatesFolder}/{file} {path}/.'
                    os.system(cmd)

    def filter_intermeds(self, cutoff):
        for file in self.intermediateFiles:
            gtf = StringTieGFF(f'{self.intermediatesFolder}/{file}')
            df = gtf.get_gtf()
            df = df[df["feature"] == "transcript"]
            attrs = df["attributes"].tolist()
            tpms = []
            genes = []
            for a in attrs:
                splat = a.split(";")
                gene = splat[0].replace("\"", "").replace("gene_id ", "")
                index = [idx for idx, s in enumerate(splat) if 'TPM ' in s and idx >= 4][0]
                tpm = splat[index].replace("TPM", "").replace("\"", "").replace(" ", "")

                # if 'ref_gene_name' not in a:
                #     # print(splat)
                #     # print(splat[3])
                #     # print(splat[4])
                #     # print(splat[5])
                #     tpm = splat[4].replace("TPM", "").replace("\"", "").replace(" ", "")
                # else:
                #     tpm = splat[7].replace("TPM", "").replace("\"", "").replace(" ", "")
                tpms.append(float(tpm))
                genes.append(gene)
            df.insert(5, "TPM", tpms)
            df.insert(6, "gene", genes)
            df = df[df["TPM"] > cutoff]
            gene_list = df["gene"].tolist()
            # print(gene_list)
            self.__add_filtered_genes(gene_list)

    def __add_filtered_genes(self, gene_list):
        for gene in gene_list:
            if gene not in self.filteredGenes:
                self.filteredGenes.append(gene)

    def __get_assembly(self):
        self.assembly = StringTieGFF(self.assembled).get_gtf()


    def filter_assembly(self, output):
        """ Filters the transcriptome assembly based on genes whose transcripts have a TPM value higher than the
        specified cutoff in at least one replicate. """
        self.__get_assembly()
        assembly = self.assembly
        attrs = assembly["attributes"].tolist()
        genes = []
        for a in attrs:
            splat = a.split(";")
            gene = splat[0].replace("\"", "").replace("gene_id ", "")
            print(gene)
            genes.append(gene)
        # print(genes)
        assembly.insert(4, "gene", genes)
        assembly = assembly[assembly["gene"].isin(self.filteredGenes)]
        assembly.to_csv(output, sep='\t', index=False)