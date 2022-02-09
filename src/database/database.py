from ..sequtils import TranscriptomeTranslator, GenomeTranslator, ORFCollection, DatabaseGenerator
from Bio import SeqIO


class Database(object):
    def __init__(self, args):
        self.args = args
        self.genome = args.genome

    def translate(self):
        rna = TranscriptomeTranslator(sequence="HISAT/transcripts.fasta", form='fasta',
                                             minsize=int(self.args.minsize), maxsize=int(self.args.maxsize))
        rna_orfs = rna.parse_frames(starts=self.args.starts.split(","), stops=self.args.stops.split(","), seqtype='aa', entry='full')

        db = DatabaseGenerator(name="transcriptome", db_type="sql")
        db.add_orfs(rna_orfs)
        db.to_fasta(filename="transcriptome_ORFs.fasta", identifier='t')
