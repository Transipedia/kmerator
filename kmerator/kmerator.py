#!/usr/bin/env python3

"""
Decomposition of transcript or gene sequences and extraction of their specific k-mers
Once:
  - Download some file of Ensembl database, regarding of requested release of transcriptome
    (once per release) to build files for kmerator
Per request:
  1. Load kmerator files of the release resquested
  3. Build sequences for each transcript
  4. Get specific kmers/contigs
  5. Concatene specific kmers/contigs files
  6. write and show report
  7. exit removing temporary file unless the --keep option is set
"""

import sys
import os
import signal
from functools import partial
import shutil
from datetime import datetime
import getpass
import pickle

### local packages
import info
from config import Config
from dataset import Dataset
from kmerize import SpecificKmers
import geneinfo
from options import usage
from color import *
import exit



def main():
    ### Handle arguments
    conf = Config(info.APPNAME)
    args = usage(conf)
    if args.debug: print(f"{DEBUG}Args:\n {args}\n{ENDCOL}")

    ### Stops running indexes on exit (Ctrl C, kill -15)
    global sigint_handler
    sigint_handler = partial(sigint_handler, args=args)
    signal.signal(signal.SIGINT, sigint_handler)

    ### check files of the dataset (download if needed)
    dataset = Dataset(args)

    ### Handle dataset options: '--list-dataset', '--rm-dataset'
    if args.list_dataset: dataset.list()
    if args.rm_dataset: dataset.remove()
    if args.mk_dataset: dataset.make()
    if args.last_avail: dataset.last_available()
    if args.update_dataset: dataset.update_last()
    if args.info: geneinfo.info(args)


    ### load kmerator geneinfo file as dict
    print(f" 🧬 Load dataset {args.release!r}.")

    ### Important objects
    items = []  # {given: xxx, ENST: xxx, ENSG: xxx, symbol: xxx, type: xxx, f_id: xxx, seq: xxx}

    report = {'failed': [], 'done': [], 'multiple': [], 'warning': []}
    transcriptome_dict = None       # not used with --fasta-file option
    geneinfo_dict = None            # not used with '--fasta-file' option

    ### load transcriptome if --selection option is set
    if args.selection:
        ### create shared objects among multiple processes with Manager()
        transcriptome_dict = dataset.load_transcriptome()

    ### load geneinfo
    geneinfo_dict = dataset.load_geneinfo()
    ### Find transcripts according the selection
    items = find_items(args, report, geneinfo_dict)

    ### get specific kmers (using multithreading)
    print(f" 🧬 Extract specific kmers, please wait..")
    SpecificKmers(args, report, items, transcriptome_dict, geneinfo_dict)

    ### Concatene results
    print(f" 🧬 Build finals results and report")
    merged_results(args)

    ### show some final info in the prompt
    show_report(args, report)

    ### set markdown report
    markdown_report(args, report)

    ### Bye
    exit.gracefully(args)


def find_items(args, report, geneinfo_dict=None):
    """
    for each name given (symbol/alias/ENSG/ENST), get information from geneinfo_dict
    """
    items = []
    ### With '--selection' option
    if args.selection:
        for given in args.selection:
            found_count = 0
            found_transcripts = []
            item = given.upper()

            ### sometimes, canonical transcript is not define (old releases), get longest transcript
            # ~ canonical = lambda args, ENSG: (geneinfo_dict['gene'][ENSG]['canonical']
                                 # ~ if 'canonical' in geneinfo_dict['gene'][ENSG]
                                 # ~ else longest_transcript(args, ENSG))

            if item in geneinfo_dict['symbol']:
                ENSGs = geneinfo_dict['symbol'][item]
                for ENSG in ENSGs:
                    found_count += 1
                    found_transcripts.append(geneinfo_dict['gene'][ENSG]['canonical'])
                    items.append({'given': given, 'ENSG': ENSG, 'type': 'symbol',
                                       'ENST':   geneinfo_dict['gene'][ENSG]['canonical'],
                                       'symbol': geneinfo_dict['gene'][ENSG]['symbol'],
                                       })
            elif item in geneinfo_dict['alias']:
                ENSGs = geneinfo_dict['alias'][item]
                for ENSG in ENSGs:
                    found_count += 1
                    found_transcripts.append(geneinfo_dict['gene'][ENSG]['canonical'])
                    items.append({'given':  given, 'ENSG': ENSG, 'type': 'alias',
                                        'ENST':   geneinfo_dict['gene'][ENSG]['canonical'],
                                        'symbol': geneinfo_dict['gene'][ENSG]['symbol'],
                                        })
            elif item in geneinfo_dict['gene']:
                found_count += 1
                ENSG = item
                items.append({'given':  given, 'ENSG': ENSG, 'type': 'gene',
                                    'ENST':   geneinfo_dict['gene'][ENSG]['canonical'],
                                    'symbol': geneinfo_dict['gene'][ENSG].get('symbol', 'N/A'),
                                    })
            elif item in geneinfo_dict['transcript']:
                found_count += 1
                ENSG = geneinfo_dict['transcript'][item]
                items.append({'given':  given, 'ENST': item, 'type': 'transcript',
                                    'ENSG':   ENSG,
                                    'symbol': geneinfo_dict['gene'][ENSG].get('symbol', 'N/A'),
                                    })
            else:
                report['failed'].append(f"{given}: not found in transcriptome")

            if found_count > 1:
                report['multiple'].append(f"{given}: {found_count} ({', '.join([i for i in found_transcripts])})")


    ### With '--fasta-file' option
    if args.fasta_file:
        type = 'transcript'
        with open(args.fasta_file) as fh:
            seq = ""
            f_id = fh.readline()[1:].split(' ')[0].rstrip()
            for raw in fh:
                if raw.startswith('>'):
                    if len(seq) >= args.kmer_length:
                        items.append({'f_id': f_id, 'seq': seq, 'type': type})
                    else:
                        report['failed'].append(f"{f_id}: sequence to short ({len(seq)} < {args.kmer_length})")
                    f_id = raw[1:].split(' ')[0].rstrip()
                    seq = ""
                else:
                    seq += raw.rstrip()
            ### last f_id/sequence
            if len(seq) >= args.kmer_length:
                items.append({'f_id': f_id, 'seq': seq, 'type': type})
            else:
                report['failed'].append(f"{f_id}: sequence to short ({len(seq)} < {args.kmer_length})")

    return items


def longest_transcript(args, ENSG):
    print(args)
    print(ENSG)
    if args.debug: print(f"{YELLOW}Finding the longest transcript for the gene {ENSG}.{ENDCOL}")
    sys.exit(f"EOF (Dev)")
    transcripts_dict = { k:len(v) for (k,v) in transcriptome_dict.items() if k.startswith(f"{gene_name}:")}
    # ~ print(*[k for k in variants_dict], sep='\n')
    nb_transcripts = len(transcript_dict)
    if args.verbose: print(f"{Color.YELLOW}Number of variants: {nb_transcripts}")
    longest_transcript = None
    length = 0
    for k,v in transcript_dict.items():
        if v > length:
            length = v
            longest_transcript = ':'.join(k.split(':')[1:2])
    return longest_transcript


def merged_results(args):
    if not os.path.isdir(os.path.join(args.tmpdir, 'kmers')):
        return None
    for item in ['kmers', 'contigs', 'masked']:
        try:
            files = os.listdir(os.path.join(args.tmpdir, item))
        except FileNotFoundError:
            return
        if files:
            merged_file = os.path.join(args.output, f"{item}.fa")
            os.makedirs(args.output, exist_ok=True)
            with open(merged_file,'wb') as mergefd:
                for file in files:
                    with open(os.path.join(args.tmpdir, item, file),'rb') as fd:
                        shutil.copyfileobj(fd, mergefd)


def show_report(args, report):
    ### show some final info in the prompt
    if report['done']:
        print(f"{CYAN}\n Done ({len(report['done'])}):")
        for i,mesg in enumerate(report['done']):
            if i == 15:
                print("  - ... (more responses)")
                break
            print(f"  - {mesg}")

    if report['multiple']:
        print(f"{BLUE}\n Multiple responses ({len(report['multiple'])}):")
        for i,mesg in enumerate(report['multiple']):
            if i == 15:
                print("  - ... (more responses)")
                break
            print(f"  - {mesg}")
            # ~ for k,v in mesg.items():
                # ~ print(f"  - {k}: {' '.join(v)}")

    if report['failed']:
        print(f"{PURPLE}\n Failed ({len(report['failed'])}):")
        for i,mesg in enumerate(report['failed']):
            if i == 15:
                print("  - ... (more responses)")
                break
            print(f"  - {mesg}")

    if report['warning']:
        print(f"{RED}\n Warning ({len(report['warning'])}):")
        for i,mesg in enumerate(report['warning']):
            if i == 15:
                print("  - ... (more responses)")
                break
            print(f"  - {mesg}")

    print(f"{ENDCOL}")


def markdown_report(args, report):
    sel_or_fa = 'selection' if args.selection else 'fasta_file'
    to_report = [
        sel_or_fa, 'datadir', 'genome', 'specie', 'kmer_length', 'release', 'stringent',
        'max_on_transcriptome', 'max_on_genome', 'output', 'thread', 'keep', 'assembly',
    ]
    os.makedirs(args.output, exist_ok=True)
    with open(os.path.join(args.output, 'report.md'), 'w') as fh:
        fh.write('# kmerator report\n')
        fh.write(f"*date: {datetime.now().strftime('%Y-%m-%d %H:%M')}*  \n")
        fh.write(f'*login: {getpass.getuser()}*\n\n')
        fh.write(f"**kmerator version:** {info.VERSION}\n\n")
        ### report command line and args, included defaults args
        cmd_args = ''
        for k,v in vars(args).items():
            if k in to_report:
                k = k.replace('_', '-')
                if isinstance(v, list):
                    v = ' '.join(v)
                cmd_args += f" \\\n  --{k} {v}"
        # ~ cmd_args = ' \\\n  '.join([f"--{k} {v}" for k,v in vars(args).items() if v])
        # ~ command = f"command: \n{__appname__}{cmd_args}"
        command = f"{info.APPNAME}{cmd_args}"
        # ~ command = ' '.join(sys.argv).replace(' -', ' \\\n  -')
        fh.write(f"**Command:**\n\n```\n{command}\n```\n\n")
        fh.write(f"**Working directory:** `{os.getcwd()}`\n\n")
        fh.write(f"**Specie:** `{args.specie}`\n\n")
        fh.write(f"**Assembly:** `{args.assembly}`\n\n")
        fh.write(f"**Transcriptome release:** `{args.release}`\n\n")
        if report['done']:
            fh.write(f"**Genes/transcripts succesfully done ({len(report['done'])})**\n\n")
            for mesg in report['done']:
                fh.write(f"- {mesg}\n")
        if report['multiple']:
            fh.write(f"\n**Multiple Genes returned for one given by Ensembl API ({len(report['multiple'])})**\n\n")
            for mesg in report['multiple']:
                fh.write(f"- {mesg}\n")
                # ~ for k,v in mesg.items():
                    # ~ fh.write(f"- {k}: {' '.join(v)}\n")
        if report['failed']:
            fh.write(f"\n\n**Genes/transcripts missing ({len(report['failed'])})**\n\n")
            for mesg in report['failed']:
                fh.write(f"- {mesg}\n")
        if report['warning']:
            fh.write(f"\n\n**Warnings ({len(report['warning'])})**\n\n")
            for mesg in report['warning']:
                fh.write(f"- {mesg}\n")


def html_report(args, report):
    '''TODO (plus option to change report format)'''
    print("Work in progress...")


def yaml_report(args, report):
    '''TODO'''
    print("Work in progress...")


def json_report(args, report):
    '''TODO'''
    print("Work in progress...")


def sigint_handler(signal, frame, args):
    exit.gracefully(args)


if __name__ == '__main__':
    main()
