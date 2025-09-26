### kmerize.py

"""
https://stackoverflow.com/questions/44587669/python-multiprocessing-how-to-close-the-multiprocessing-pool-on-exception
https://stackoverflow.com/questions/47903791/how-to-terminate-a-multiprocess-in-python-when-a-given-condition-is-met
"""

import os
import sys
import shutil
import subprocess
import multiprocessing
import copy
import random
import string

from color import *

MAX_CHARS = os.pathconf('/', 'PC_NAME_MAX') - 40   # with --fasta-file arg: avoid too large file name

class SpecificKmers:
    """ Class doc """

    def __init__(self, args, report, items, transcriptome_dict=None, geneinfo_dict=None):
        """
        1. Get sequence
        2. write on file
        3. From sequence, build abundance for each kmer, using jellyfish againt the genome (output: dict)
        4. From sequence, build abundance for each kmer, uding jellyfish againt the transcriptome (output: dict)
        5. filter the specific kmers, according to the given arguments
        6. write specific kmers and contigs in files
        """

        ### create a shared dict among multiple processes with Manager()
        manager = multiprocessing.Manager()
        self.args = manager.dict(args.__dict__)
        if args.selection:
            self.transcriptome_dict = manager.dict(transcriptome_dict)
            self.geneinfo_dict = manager.dict(geneinfo_dict)

        ### launch workers
        self.transcriptome_file = os.path.join(args.datadir, f"{args.specie}.{args.assembly}.{args.release}.k{args.kmer_length}.transcriptome.jf")
        try:
            with multiprocessing.Pool(processes=args.thread) as pool:
                if args.selection:
                    messages = pool.map(self.worker_selection, items)
                elif args.fasta_file:
                    messages = pool.map(self.worker_fasta_file, items)
            for type,mesg in messages:
                report[type].append(mesg)
        except KeyError as e:
            pool.close()
            shutil.rmtree(args.tmpdir)
            sys.exit(f"{RED}Error: {e.args[0]}")


    def worker_selection(self, item):
        """ Function doc """
        ### unpack item
        globals().update(item)
        global level
        level = 'transcript' if type == 'transcript' else 'gene'

        ### 1. Get sequence of the transcript (from the transcriptome)
        try:
            seq = self.transcriptome_dict[ENST]
        except KeyError:
            mesg = f"{given}: transcript not found in transcriptome ({ENST})."
            return 'failed', mesg
        if len(seq) < self.args['kmer_length']:
            mesg = f"{given}: sequence to short ({len(seq)} < {self.args['kmer_length']})."
            return 'failed', mesg

        ### 2. Write sequences in temporary file (needed by jellyfish)
        f_id = f"{item['given']}.{item['ENST']}"
        seq_file, fasta_name = self.dump_seq(f_id, seq)

        ### 3. From sequence, compute jellyfish againt the genome/transcriptome and convert results as dict ()
        kmercounts_genome_dict = self.jellyfish(seq_file, self.args['genome'])
        kmercounts_transcriptome_dict = self.jellyfish(seq_file, self.transcriptome_file)
        ### 4. filter the specific kmers, according to the arguments
        mesg = self.get_specific_kmers(item, kmercounts_transcriptome_dict, kmercounts_genome_dict, fasta_name)
        return mesg


    def worker_fasta_file(self, item):
        ### unpack item
        globals().update(item)
        ### 1. write separate file for each item
        seq_file, fasta_name = self.dump_seq(item['f_id'], item['seq'])

        ### 2. From sequence, compute jellyfish against the genome/transcriptome and convert results as dict ()
        kmercounts_transcriptome_dict = self.jellyfish(seq_file, self.transcriptome_file)
        kmercounts_genome_dict = self.jellyfish(seq_file, self.args['genome'])

        ### 3. find to specific kmers
        mesg = self.get_specific_kmers(item, kmercounts_transcriptome_dict, kmercounts_genome_dict, fasta_name)

        return mesg


    def dump_seq(self, f_id, seq):
        '''
        Write sequences in temporary file (needed by jellyfish)
        '''
        seq_dir = os.path.join(self.args['tmpdir'], 'sequences')
        os.makedirs(seq_dir, exist_ok=True)
        ### replace forbiden characters
        fasta_name = f_id.replace('/', '--')
        ### handle max file lenght
        if len(fasta_name) > MAX_CHARS:
            fasta_name = f"{f_id[:MAX_CHARS-15]}...{''.join(random.choices(string.ascii_letters, k=11))}"
        ### write sequence as fasta
        seq_file = os.path.join(seq_dir, f"{fasta_name}.fa")
        with open(seq_file, 'w') as fh:
            fh.write(f">{fasta_name}\n{seq}\n")
        return seq_file, fasta_name


    def jellyfish(self, seq_file, jf_file):
        """
        From sequence, compute jellyfish againt the genome/transcriptome and convert results
        as dict ()
        """
        if self.args['debug']: print(f"{YELLOW}start jellyfish on {os.path.basename(seq_file)} "
                                     f"against {os.path.basename(jf_file)}{ENDCOL}")
        cmd = f'jellyfish query -s "{seq_file}" {jf_file}'
        try:
            result = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT) \
                                             .rstrip('\n').split('\n')
        except subprocess.CalledProcessError as err:
            sys.exit(f"{ERROR}Error: executing jellyfish:\n"
                     f"  {ERROR}command: {ENDCOL}{cmd}\n"
                     f"  {ERROR}returned: {ENDCOL}{err.output}")
        if not result[0]:
            return {}
        result_dict = dict([(b[0], int(b[1])) for b in [a.split() for a in result]])
        return result_dict


    def get_specific_kmers(self, item, kmercounts_transcriptome_dict, kmercounts_genome_dict, fasta_name):
        '''
        Keep only specific kmers, according to the arguments
        Launched for each gene, transcript, or unanotated sequence
        '''
        args = copy.deepcopy(self.args)
        ### Define some variables: gene_name, transcript_name, variants_dic and output file names
        '''
        d = {
            'specific_kmers': [],              # specific kmers list
            'specific_contigs': [],            # specific contigs list
            masked_kmers = []                  # masked kmers list
            'contig': "",                      # initialize contig sequence
            'knb': 0,                          # kmer number (selected kmer)
            'c_nb': 1,                         # contig number
            'kmer_pos_prev': 0,                # position of retained kmer
            'contig_pos': 0,                   # position of retained contig
        }
        '''
        level = 'gene' if item['type'] != 'transcript' else item['type']
        specific_kmers = []              # specific kmers list
        masked_kmers = []                # masked kmers list
        specific_contigs = []            # specific contigs list
        contig = ""                      # initialize contig sequence
        knb = 0                          # kmer number (selected kmer)
        c_nb = 1                         # contig number
        kmer_pos_prev = 0                # position of retained kmer
        contig_pos = 0                   # position of retained contig

        ### initialization of count variables
        total_kmers = len(kmercounts_transcriptome_dict)
        if args['selection']:
            given_up = given.upper()
            isoforms = self.geneinfo_dict['gene'][ENSG]['transcript']
            isoforms_nb = len(isoforms)
            isoforms_dict = { isoform:self.transcriptome_dict[isoform] for isoform in isoforms }

            ### Define conditionnal variables
            ## When '--selection' option is set
            kmer_outfile = f"{given_up}-{ENST}-{level}-specific_kmers.fa"
            contig_outfile = f"{given_up}-{ENST}-{level}-specific_contigs.fa"
            masked_outfile = f"{given_up}-{ENST}-{level}-masked_kmers.fa"
        # ~ elif args['chimera']:
            # ~ ## When '--chimera' option is set
            # ~ kmer_outfile = f"{fasta_name}-chimera-specific_kmers.fa"
            # ~ contig_outfile = f"{fasta_name}-chimera-specific_contigs.fa"
        elif args['fasta_file']:
            ## When '--fasta-file' option is set
            kmer_outfile = f"{fasta_name}-transcript-specific_kmers.fa"
            contig_outfile = f"{fasta_name}-transcript-specific_contigs.fa"
            masked_outfile = f"{fasta_name}-masked_kmers.fa"


        i = 1
        for mer, abund_in_tr in kmercounts_transcriptome_dict.items():

            kmer_pos = i
            try:
                canonical = self.__canonical(mer)
                abund_in_ge = int(kmercounts_genome_dict[canonical])    # abundance in genome for this kmer
            except KeyError as err:
                if len(mer) != len(next(iter(kmercounts_genome_dict))):
                    raise KeyError(f"ErrorIndexLength: length of kmer expected: {self.args['kmer_length']}\n"
                                   f"  Genome kmer index length: {len(next(iter(kmercounts_genome_dict)))}\n"
                                   f"  Transcriptome kmer index length (dataset): {len(mer)}")
                raise KeyError(f"Error: kmer not found in genome: {err}")

            if level == 'gene':

                isoforms_with_mer = [enst for enst, seq in isoforms_dict.items() if mer in seq]
                isoforms_with_mer_nb = len(isoforms_with_mer)

                ### if the kmer is present/unique or does not exist (splice/chimera?) on the genome
                if abund_in_ge <= 1:

                    # knb, specific_kmers, given_up, transcript, kmer_pos, isoforms_with_mer_nb,
                    # isoforms_nb, mer, contig, kmer_pos_prev, kmer_pos, contig_pos

                    ### all kmers found (abund_in_tr) must be in isoforms
                    if not args['stringent'] and abund_in_tr and abund_in_tr == isoforms_with_mer_nb:
                        ## contigs case
                        if knb == 0:                                    # first kmer in contig
                            contig = mer
                            kmer_pos_prev = kmer_pos
                            contig_pos = kmer_pos
                        elif knb>0 and kmer_pos == kmer_pos_prev+1:     # add last base to existing contig
                            contig += mer[-1]
                            kmer_pos_prev = kmer_pos
                        else:                                           # store contig and create new
                            specific_contigs.append(f">{given_up}:{ENST}.contig_{c_nb} (at position {contig_pos})\n{contig}")
                            c_nb += 1
                            contig = mer
                            kmer_pos_prev = kmer_pos
                            contig_pos = kmer_pos
                        ## kmers case
                        knb += 1
                        specific_kmers.append(f">{given_up}:{ENST}.kmer{kmer_pos} ct:{c_nb} tr:{isoforms_with_mer_nb}/{isoforms_nb}\n{mer}")

                    ### for gene level only, the stringent argument implies retaining kmers present in ALL isoforms of the gene
                    elif args['stringent'] and abund_in_tr == isoforms_nb == isoforms_with_mer_nb:
                        ## contigs case
                        if knb == 0:                                    # first kmer in contig
                            contig = mer
                            kmer_pos_prev = kmer_pos
                            contig_pos = kmer_pos
                        elif knb>0 and kmer_pos == kmer_pos_prev+1:     # add last base to existing contig
                            contig += mer[-1]
                            kmer_pos_prev = kmer_pos
                        else:                                           # store contig and create new
                            specific_contigs.append(f">{given_up}:{transcript}.contig_{c_nb} (at position {kmer_pos})\n{contig}")
                            c_nb += 1
                            contig = mer
                            kmer_pos_prev = kmer_pos
                            contig_pos = kmer_pos
                        ## kmers case
                        knb += 1
                        specific_kmers.append(f">{given_up}:{ENST}.kmer{kmer_pos} ct:{c_nb} tr:{isoforms_with_mer_nb}/{isoforms_nb}\n{mer}")
                    else:
                        masked_kmers.append(f">{given_up}:{ENST}.kmer{kmer_pos} "
                                        f"tr:{isoforms_with_mer_nb}/{isoforms_nb} "
                                        f"transcriptome:{abund_in_tr}\n{mer}"
                                        )

                ### the kmer count exceeded 1 in the genome
                else:
                    masked_kmers.append(f">{given_up}:{ENST}.kmer{kmer_pos} "
                                    f"tr:{isoforms_with_mer_nb}/{isoforms_nb} "
                                    f"genome:{abund_in_ge}\n{mer}"
                                    )

            ### Cases of transcripts 1) annotated, 2) unannotated.
            elif level == 'transcript':
                ### Case of annotated transcripts
                if args['selection'] and abund_in_tr == 1 and abund_in_ge <= 1:
                    ## contigs case
                    if knb == 0:                                    # first kmer in contig
                        contig = mer
                        kmer_pos_prev = kmer_pos
                        contig_pos = kmer_pos
                    elif knb>0 and kmer_pos == kmer_pos_prev+1:     # add last bp to existing contig
                        contig += mer[-1]
                        kmer_pos_prev = kmer_pos
                    else:                                           # store contig and create new
                        specific_contigs.append(f">{ENST}.contig_{c_nb} (at position {contig_pos})\n{contig}")
                        c_nb += 1
                        contig = mer
                        kmer_pos_prev = kmer_pos
                        contig_pos = kmer_pos
                    ## kmers case
                    knb += 1
                    specific_kmers.append(f">{given_up}:{ENST}.kmer{kmer_pos} ct:{c_nb}\n{mer}")

                ### Case of unannotated transcripts
                elif args['fasta_file'] and abund_in_tr <= self.args['max_on_transcriptome'] and abund_in_ge <= self.args['max_on_genome']:     # max_on_transcriptome = 0 by default
                    ### contigs case
                    if knb == 0:                                    # first kmer in contig
                        contig = mer
                        kmer_pos_prev = kmer_pos
                        contig_pos = kmer_pos
                    elif knb>0 and kmer_pos == kmer_pos_prev+1:     # add last bp to existing contig
                        contig += mer[-1]
                        kmer_pos_prev = kmer_pos
                    else:
                        specific_contigs.append(f">{f_id}.contig_{c_nb} (at position {contig_pos})\n{contig}")
                        c_nb += 1
                        contig = mer
                        kmer_pos_prev = kmer_pos
                        contig_pos = kmer_pos
                    ### kmers case
                    knb += 1
                    specific_kmers.append(f">{f_id}.kmer{kmer_pos} ct:{c_nb}\n{mer}")

                ### the kmer count exceeded 1 in the genome
                else:
                    fasta_id = f"{given_up}:{ENST}.kmer{kmer_pos}" if args['selection'] else f"{f_id}.contig_{c_nb}"
                    masked_kmers.append(f">{fasta_id} "
                                    f"genome:{abund_in_ge} "
                                    f"transcriptome:{abund_in_tr}\n{mer}"
                                    )

            # ~ ### Case of chimera
            # ~ elif level == 'chimera':
                # ~ if abund_in_tr == abund_in_ge == 0:
                    # ~ ### contig case
                    # ~ if knb == 0:
                        # ~ contig = mer
                        # ~ kmer_pos_prev = kmer_pos
                        # ~ contig_pos = kmer_pos
                    # ~ elif knb>0 and kmer_pos == kmer_pos_prev+1:
                        # ~ contig += mer[-1]
                        # ~ kmer_pos_prev = kmer_pos
                    # ~ else:
                        # ~ specific_contigs.append(f">{f_id}.contig_{c_nb} (at position {contig_pos})\n{contig}")
                        # ~ c_nb += 1
                        # ~ contig = mer
                        # ~ kmer_pos_prev = kmer_pos
                        # ~ contig_pos = kmer_pos
                    # ~ ### kmers case
                    # ~ knb += 1
                    # ~ specific_kmers.append(f">{f_id}.kmer{kmer_pos} ct:{c_nb}\n{mer}")
            ### not a gene or transcript
            else:
                raise KeyError(f"{RED}Error: level {level!r} unknown.{ENDCOL}")
            i += 1


        ### append last contig in list
        if args['fasta_file'] and contig:
            specific_contigs.append(f">{f_id}.contig_{c_nb} (at position {contig_pos})\n{contig}")
        elif level == "gene" and contig:
            specific_contigs.append(f">{given_up}:{ENST}.contig_{c_nb} (at position {contig_pos})\n{contig}")
        elif level == "transcript" and contig:
            specific_contigs.append(f">{ENST}.contig_{c_nb} (at position {contig_pos})\n{contig}")
        if args['debug']:
            if args['selection']:
                print(f"{YELLOW} {ENST} kmers/contig: {len(specific_kmers)}/{len(specific_contigs)} ({given}){ENDCOL}")
            else:
                print(f"{YELLOW} {f_id} kmers/contig: {len(specific_kmers)}/{len(specific_contigs)}{ENDCOL}")

        ### masked kmers
        if masked_kmers:
            self.write(args, masked_outfile, masked_kmers, 'masked')

        ### write kmer/contig files
        if specific_kmers:
            self.write(args, kmer_outfile, specific_kmers, 'kmers')
            self.write(args, contig_outfile, specific_contigs, 'contigs')
        else:
            if args['selection']:
                mesg = f"{given}: no specific kmers found for {ENST} (level: {level})"
            else:
                mesg = f"{f_id}: no specific kmers found"
            return 'failed', mesg

        ### report
        if args['selection']:
            mesg = f"{given}: {item['symbol']}:{ENST} - kmers/contigs: {len(specific_kmers)}/{len(specific_contigs)} (level: {level})"
        else:
            mesg = f"{f_id} - kmers/contigs: {len(specific_kmers)}/{len(specific_contigs)} (level: {level})"

        return 'done', mesg


    def write(self, args, outfile, specific_seq, type):
        """ Write results as fasta file"""
        outdir = os.path.join(args['tmpdir'], type)
        os.makedirs(outdir, exist_ok=True)
        with open(os.path.join(outdir, outfile.replace('/', '--')), 'w') as fh:
            fh.write("\n".join(specific_seq) + '\n')


    def __canonical(self, mer):
        revcomp = mer[::-1].translate(str.maketrans('ATCG','TAGC'))
        return sorted((mer, revcomp))[0]


'''
def handle_kmers(**kwargs):
    print("LOCALS():", locals())
    ## kmers case
    knb += 1
    specific_kmers.append(f">{given_up}:{given_up}.kmer{kmer_pos} ({isoforms_with_mer_nb}/{isoforms_nb})\n{mer}")
    ## contigs case
    if knb == 1:                                    # first kmer in contig
        contig = mer
        kmer_pos_prev = kmer_pos
        print(f"{given_up}:{transcript}.contig{c_nb}\n{contig}")
    elif knb>1 and kmer_pos == kmer_pos_prev+1:     # add last pb to existing contig
        contig = contig + mer[-1]
        kmer_pos_prev = kmer_pos
    else:                                           # store contig and create new
        specific_contigs.append(f"{given_up}:{transcript}.contig{c_nb}\n{contig}")
        c_nb += 1
        contig = mer
        kmer_pos_prev = kmer_pos
'''
