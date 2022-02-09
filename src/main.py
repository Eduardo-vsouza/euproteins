import os
import sys
from src.database import database_generator as dg
from . import peptide_search as ps
from . import results_new_approach as pms
#from . import the_visualizer as vis
# from . import orthologs as phylo
from . import f_translation as trans
#from .working_runs import OrganizePlot
#from .Digestion_sets import Digestion, Digested, PlotData
# import prowser.gene_organizer as gorg
# import prowser.browser_gui_v16 as prsr
# from .testing import uproteins_testing as test
from .percolator import Decoy
from .assembly import TranscriptAssembly, CompareTranscripts, ReadMapper
from .master import Archives
from .database import Database, SieveFiltering, DeepMSPeptide, Digestor
from .postprocess import PostPercolator, ExtendedInformation
from .sequtils.orflib import AltCodons
from .upstream import SDInspection
from .testing import PipelineTesting
from .postms import TSVConverter
from .pipelines import PostMSPipeline, ValidatePipeline
from .forest import ProteinFixer, PreFiltering, FeatureFishing, SpectralForest
from .transcriptomics import ExpressionFilter


pypath = sys.path[0]


def run_workflow(args):
    run = True
    # argfix = ArgumentsFixer(parser)
    # args = argfix.get_abspath()
    mode = args.mode
    args.Transcriptome = "YES"
    if run:
        if not os.path.exists(args.outdir):
            cmd_dir = "mkdir %s" % args.outdir
            os.system(cmd_dir)
        os.chdir(args.outdir)
        # if args.smorfs is not None:
        #     if args.smorfs != "YES" or args.smorfs != "NO":
        #         print("\nInvalid argument. smorfs argument must be YES or NO.\n")
        #         return os.EX_USAGE

        if mode == "assembly":
            arxivs = Archives()
            if args.step <= 1:
                mapping = ReadMapper(args)
                mapping.create_indexes()
                mapping.align_reads()
            assembly = TranscriptAssembly(args)
            if args.step <= 2:
                assembly.assemble()
            if args.step <= 3:
                assembly.create_gtf_list()
                gtf = assembly.merge_transcripts(tpm_cutoff=10, fpkm_cutoff=10)
                arxivs.assembledGTF = gtf
                comparisons = CompareTranscripts(args, arxivs.assembledGTF)
                compare_dir = comparisons.run_gffcompare()
                arxivs.comparisonsDirectory = compare_dir
                sequences = comparisons.extract_sequences()
                arxivs.RNASequences = sequences
            # if args.step <= 4:
                # assembly_filtering = ExpressionFilter(gtf_intermediate_folder="StringTie/gtf_intermediates",
                #                                       assembled_gtf='assembled.gtf')
                # assembly_filtering.filter_intermeds(cutoff=0.5)
                # assembly_filtering.filter_assembly(output='StringTie/filtered_assembly.gtf')

        elif mode == "testing":
            testing = PipelineTesting(args)
            testing.run()

        elif mode == "database":
            filetype = 'transcriptome'
            folder = 'Transcriptome'
            if args.step <= 1:
                db = Database(args)
                db.translate()
                print("Generating the transcriptome database.")
            if args.step <= 2:
                transcriptome_db = dg.Database("transcriptome_ORFs.fasta", args.proteome, "transcriptome")
                transcriptome_db.mark_annotated()
                transcriptome_db.unify()
                print("Transcriptome database generated. Now filtering out sequences that do not have a proteotypic "
                      "peptide. \n")
            # if args.step <= 3:
            #     print("Performing an in silico digestion of the database.\n")
            #     digests = Digestor(filetype=filetype, folder=folder, threads=args.threads)
            #     digests.digest()
            #     peps_and_orfs = digests.filter_by_size()
            # elif args.step <= 4:
            #     print("Predicting proteotypic peptides.\n")
            #     proteotypic = DeepMSPeptide(filetype='transcriptome', folder='Transcriptome')
            #     proteotypic.predict()
            #     proteotypic.filter_proteotypic(peps_and_orfs=peps_and_orfs)
            #     proteotypic.remove_dupli()
            #     proteotypic.filter_orf_db(args.threads)
            #     proteotypic.readd_annotation(args.proteome)
            print("Database generation step complete. Look for databases in %s" % args.outdir)

        elif mode == "ms":
            # transcriptome = ps.PeptideSearch("Transcriptome", args.Mass_spec, "transcriptome_database.fasta", args)
            # transcriptome.peptide_identification()
            transcriptome_decoy = Decoy(db="transcriptome_database.fasta", db_type="Transcriptome")
            transcriptome_decoy.reverse_sequences().to_fasta()
            # transcriptome_decoy_search = ps.PeptideSearch("Transcriptome", args.Mass_spec, "Transcriptome/Percolator/Transcriptome_decoy.fasta", args, decoy=True)
            # transcriptome_decoy_search.peptide_identification()

        elif mode == "postms":
            """ newest method """
            # genome = PostMSPipeline(args=args, filetype='genome', folder='Genome')
            # genome.run()
            # if args.Transcriptome == 'YES':
            transcriptome = PostMSPipeline(args=args, filetype='transcriptome', folder='Transcriptome')
            transcriptome.run()

            print("DONE. Results are inside Genome or Transcriptome folders.")

        elif mode == "validate":
            validation = ValidatePipeline(args=args)
            validation.validate_genome()
            validation.validate_transcriptome()

        else:
            print("Invalid mode. Please choose a valid mode.")
