
import os
import sys
import requests
from bs4 import BeautifulSoup
from urllib.request import urlretrieve
import gzip
import pickle
import shutil
import tempfile
import threading

from color import *


class TranscriptomeBuilder:
    """
    - download 'cdna' and 'ncrna' from Ensembl ftp site
    - remove transcripts not located in regular chromosomes (avoid duplicates)
    - concatene a temporary file with results
    - build a jellyfish index
    - delete temporary files
    """

    def __init__(self, args, base_url, transcriptome_jf, report, chr):
        """ Class initialiser """
        self.args = args
        self.base_url = base_url
        self.report = report
        self.chr = chr
        self.transcript_count = 0

        if self.args.debug: print(f"{DEBUG}Build transcriptome, please wait...{ENDCOL}")
        ### select fasta file to download
        release = self.args.release

        self.tmpdir = tempfile.mkdtemp(prefix="kmerator_")
        self.downloaded_fasta = []
        self.temp_fasta_files = []

        if self.args.debug: print(f"{DEBUG}Download transcriptome files (release {release}), please wait...{ENDCOL}")
        ### download fasta gzipped files and clean them
        for item in ('cdna', 'ncrna'):
            self.build_fasta(item)

        ### concatene filtered cDNA and ncRNA (and remove temp fasta files)
        transcriptome_fa = self.mk_transcriptome()

        ### make index of transcriptome (kmc or jellyfish)
        self.mk_index(transcriptome_fa, transcriptome_jf)

        ### create a pickle file of transcriptome
        self.mk_pickle(transcriptome_fa)

        '''
        ### make index of transcriptome (kmc or jellyfish)
        th1_1 = threading.Thread(target=self.mk_index,
                               args=(transcriptome_fa, transcriptome_jf))
        ### create a pickle file of transcriptome
        th1_2 = threading.Thread(target=self.mk_pickle,
                               args=(transcriptome_jf,))
        th1_1.start()
        th1_2.start()
        th1_1.join()
        th1_2.join()
        '''

        ### Add to report
        self._to_report()

        ### remove intermediate files
        shutil.rmtree(self.tmpdir, ignore_errors=True)
        if not args.keep:
            os.remove(transcriptome_fa)


    def build_fasta(self, item):
        """
        for cdna and ncrna:
          - define url
          - download file
          - clean file
        """
        ### get url links for cDNA and ncRNA fasta files
        url = os.path.join(self.base_url, f"release-{self.args.release}", "fasta", self.args.specie, item)
        file_name = self.get_link(url, item)
        link = os.path.join(url, file_name)

        ### Download fasta files
        fasta_path = self.wget_fasta(link, file_name)
        self.downloaded_fasta.append(fasta_path)

        ### check fasta file
        self.check_fasta(fasta_path)

        ### create a temp file, with alternate chromosome removed
        temp_fasta_path = self.filtered_fasta(fasta_path)
        self.temp_fasta_files.append(temp_fasta_path)


    def get_link(self, url, item):
        if item == 'cdna':
            pattern = 'cdna.all.fa.gz'
        elif item == 'ncrna':
            pattern = 'ncrna.fa.gz'
        try:
            response = requests.get(url, timeout=10)
        except requests.exceptions.ConnectionError:
            sys.exit("Error: unable to connect to Ensembl API.")
        if response.ok:
            response_text = response.text
        else:
            return response.raise_for_status()
        soup = BeautifulSoup(response_text, 'html.parser')
        href = [link.get('href') for link in soup.find_all('a') if pattern in link.get('href')][0]
        return href


    def wget_fasta(self, link, file_name):
        gz_file = file_name.split('.')
        gz_file = f"{'.'.join(gz_file[0:2])}.{self.args.release}.{'.'.join(gz_file[2:])}"
        fasta_path = os.path.join(self.tmpdir, gz_file)
        urlretrieve(link, fasta_path)
        return fasta_path


    def check_fasta(self, fasta_file):
        try :
            with gzip.open(fasta_file) as fh:
                first_line = fh.readline().rstrip().decode()
                if not first_line.startswith(">"):
                    sys.exit(f"{RED}Error: Are you sure {fasta} is a fasta file ?{ENDCOL}")
                return first_line
        except FileNotFoundError:
            sys.exit(f"{RED}Error: File '{fasta}' not found.{ENDCOL}")


    def filtered_fasta(self, fasta_file_path):
        ## Handle sequences
        sequences = []
        seq = []
        current_header = ""
        end_head =  '\n'            # '\t' if args.tsv else '\n'
        sep_seq  =  '\n'            # '' if args.uniq or args.tsv else '\n'
        with gzip.open(fasta_file_path, 'rt') as fh:
            previous_header = fh.readline()
            for line in fh:
                if line.rstrip():
                    if line.startswith('>'):
                        current_header = line
                        # ~ if not 'CHR_' in previous_header.split()[2]:
                        ### keep only transcripts in regular chromosomes
                        if previous_header.split()[2].split(':')[2] in self.chr:
                            ## add only the transcript name in the header (without version number)
                            sequences.append(f"{previous_header.split('.')[0]}\n{''.join(seq)}\n")
                        seq = []
                        previous_header = current_header
                    else:
                        seq.append(line.rstrip())
        ### last fasta sequence is not printed in the loop
        if previous_header.split()[2].split(':')[2] in self.chr:
        # ~ if not 'CHR_' in previous_header.split()[2]:
            sequences.append(f"{previous_header.split('.')[0]}\n{''.join(seq)}\n")
        ### write temp file
        fasta_file_basename = ".".join(os.path.basename(fasta_file_path).split('.')[:-2])
        temp_fasta_path = os.path.join(self.tmpdir, f'{fasta_file_basename}-altCHR.fa')
        with open(temp_fasta_path, 'w') as fh:
            for sequence in sequences:
                fh.write(sequence)
        ### check if temporary file done
        if not os.path.isfile(temp_fasta_path):
            sys.exit(f"Error: temporary {temp_fasta_path!r} not found.")
        return temp_fasta_path


    def mk_transcriptome(self):
        ### define output for filtered transcriptome
        specie, version, release = self.temp_fasta_files[0].split('.')[:3]
        specie = os.path.basename(specie).lower()
        transcriptome_fa = f"{specie}.{version}.{release}.transcriptome.fa"
        transcriptome_fa = os.path.join(self.args.datadir, transcriptome_fa)
        ### concatene cDNA and ncRNA
        if self.args.debug: print(f"{DEBUG}creating transcriptome {os.path.basename(transcriptome_fa)!r}.{ENDCOL}")
        with open(transcriptome_fa, 'w') as outfile:
            for fasta in self.temp_fasta_files:
                with open(fasta, 'r') as infile:
                    # ~ outfile.write(infile.read())        # small files
                    for line in infile:                     # big files
                        outfile.write(line)
        ''' Other method
        with open(transcriptome_fa, 'w') as fout, fileinput.input(temp_fasta_files, 'rb') as fin:
            for line in fin:
                fout.write(line)
        '''
        return transcriptome_fa


    def mk_index(self, transcriptome_fa, transcriptome_jf):
        """ Function doc """
        ### Select best tool to indexing transcriptome
        tools = ['_kmc', 'jellyfish']
        num = None
        for i,bin in enumerate(tools):
            if shutil.which(bin):
                num = i
                break
        if not num:
            print(f"Warning: no tool found to index the transcriptome")
            return None

        ### Define command line according to the tool used
        ''' in the future (python 3.10)
        match tools[num]:
            case 'kmc':
                print(f"prefered tool: {os.path.basename(tool)}")
            case 'jellyfish':
                cmd = f"{tool} count -m self.args.kmer_length -s 100000 {transcriptome_fa} -o {transcriptome_idx}"
        '''
        tool = tools[num]
        ## kmc case
        if tool == 'kmc': print(f"prefered tool: {os.path.basename(tool)}")
        ## jellyfish case
        elif tool == 'jellyfish':
            filename = os.path.basename(f"{os.path.splitext(transcriptome_fa)[0]}.jf")
            cmd = (f"{tool} count -t {self.args.thread} -m {self.args.kmer_length} -s 100000 "
                   f"{transcriptome_fa} -o {transcriptome_jf}")

        ### Build index
        if self.args.debug: print(f"{DEBUG}Build index with {tool!r}, please wait...{ENDCOL}")
        ret = os.system(cmd)


    def mk_pickle(self, transcriptome_fa):
        """ Build pickle file instead the fasta file (faster to load) """
        fa_dict = {}
        with open(transcriptome_fa) as fh:
            header = fh.readline()[1:].split()[0].rstrip()
            for line in fh:
                if line.startswith('>'):
                    fa_dict[header] = seq
                    header = line[1:].rstrip().split()[0].split('.')[0]
                else:
                    seq = line.rstrip()
        ### Don't forget the last record
        fa_dict[header] = seq
        ### how many records ?
        self.transcript_count = len(fa_dict)
        ### write as pickle file
        transcriptome_pkl = f"{os.path.splitext(transcriptome_fa)[0]}.pkl"
        with open(transcriptome_pkl, 'wb') as fh:
            pickle.dump(fa_dict, fh)


    def _to_report(self):
        """ Function doc """
        self.report.append("\n## Transcriptome Info\n")
        self.report.append(f"- nb transcripts: {self.transcript_count}\n")
