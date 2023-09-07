#!/usr/bin/env python3

"""
Give some information about gene
examples:
  kmerator --info npm1 ensg00000268555 ENSG00000258154
  kmerator --info file.txt (list of items separated by space/tab/newline. comment: #)
"""

import sys
import os
import pickle

### Local modules
from dataset import Dataset
import color as col
import exit


def info(args):
    """ Function doc """

    ### load geneinfo pickle file as dict
    dataset = Dataset(args)
    geneinfo_dict = dataset.load_geneinfo()
    transcriptome_dict = dataset.load_transcriptome() if args.all else None
    ### get info for each item
    info, not_found = get_info(args, geneinfo_dict)
    ### print info
    output(args, geneinfo_dict, info, not_found, transcriptome_dict)
    ### quit removing temporary files
    exit.gracefully(args)


def get_info(args, geneinfo_dict):
    """
    according to the type of info (symbol, Ensembl gene, alias, Ensembl transcript),
    get information on the gene (or parent gene) from geneinfo_dict
    """
    info = {}
    not_found = []

    for given in args.info:
        ### if item given is a symbol
        if given.upper() in geneinfo_dict['symbol']:
            info[given] = {'_type': 'gene symbol'}
            ENSGs = geneinfo_dict['symbol'][given.upper()]
            for ensg in ENSGs:
                info[given][ensg] = geneinfo_dict['gene'][ensg]
        ### if item given is an alias
        elif given.upper() in geneinfo_dict['alias']:
            info[given] = {'_type': 'alias'}
            ENSGs = geneinfo_dict['alias'][given.upper()]
            for ensg in ENSGs:
                info[given][ensg] = geneinfo_dict['gene'][ensg]
        ### if item given is a ensembl gene name (ex: ENSG00000117)
        elif given.upper() in geneinfo_dict['gene']:
            info[given] = {'_type': 'ensembl gene name'}
            info[given][given] = geneinfo_dict['gene'][given.upper()]
        ### if item given is a ensembl transcript name (ex:ENST00000633522 )
        elif given.upper() in geneinfo_dict['transcript']:
            info[given] = {'_type': 'transcript'}
            ensg = geneinfo_dict['transcript'][given.upper()]
            info[given][given.upper()] = geneinfo_dict['gene'][ensg]
        else:
            not_found.append(given)

    return info, not_found


    for given in args.info:
        found_count = 0
        found_transcripts = []
        item = given.upper()

        if item in geneinfo_dict['transcript']:
            found_count += 1
            ENSG = geneinfo_dict['transcript'][item]
            items.append({'given':  given, 'ENST': item, 'type': 'transcript',
                                'ENSG':   ENSG,
                                'symbol': geneinfo_dict['gene'][ENSG]['symbol'],
                                })
        else:
            report['failed'].append(f"{given}: not found in transcriptome")

        if found_count > 1:
            report['multiple'].append(f"{given}: {found_count} ({', '.join([i for i in found_transcripts])})")


def output(args, geneinfo_dict, info, not_found, transcriptome_dict):
    """ Function doc """

    ### output as readable format
    for given, ENSGs in info.items():
        number_found = "" if ENSGs['_type'] == "ensembl gene name" else f"({len(ENSGs)-1} found) "

        print(f"\n=== {given} ({ENSGs['_type']}) {number_found} ===")
        given_type = ENSGs.pop('_type')     # gene symbol, alias, ENSG, transcript

        for ensg, val in ENSGs.items():
            blank = "\n" + 23 * " "
            print(f"  Ensembl ID           {ensg}")
            print(f"  Gene Name            {val.get('symbol', '')}")
            if args.all: print( "  Aliases              {}".format(blank.join(val.get('aliases', ''))))
            print(f"  Specie               {args.specie}")
            print(f"  Assembly             {geneinfo_dict['assembly']}")
            print(f"  Coordinates          {val['chr']}:{val['start']}-{val['end']}")
            print(f"  Strand               {val['strand']}")
            print(f"  Canonical transcript {val['canonical']}")
            print(f"  Biotype              {val.get('biotype', 'unknown')}")
            print(_text_format(23, '  Description', val.get('desc', 'unknown')))
            ### Print transcrpts (or not)
            if args.all:
                if given_type == "transcript":
                    transcript = next(iter(info))
                    seq = transcriptome_dict.get(transcript, f"{col.RED}✘✘✘{col.ENDCOL} (sequence not found)")
                    length = f" ({len(seq)})" if seq[0] in ['A', 'T', 'C', 'G'] else ""
                    print(_text_format(23, f"{transcript}{length}", seq, pos_key='top'))
                else:
                    print(f"  Transcripts of {val.get('symbol', ensg)}")
                    transcripts = val['transcript']
                    transcripts.remove(val['canonical'])    # I want canonical at first
                    transcripts.insert(0, val['canonical']) # I want canonical at first
                    for transcript in transcripts:
                        seq = transcriptome_dict.get(transcript, f"{col.RED}✘✘✘{col.ENDCOL} (sequence not found)")
                        length = f" ({len(seq)})" if seq[0] in ['A', 'T', 'C', 'G'] else ""
                        print(_text_format(23, f"{transcript}{length}", seq, pos_key='top'))
            else:
                print("  Transcripts          {}".format(blank.join(val['transcript'])))

            print()

    ### output not found items
    if not_found:
        print("=== NOT FOUND ===\n", *[f" {i}\n" for i in not_found])


def _text_format(marg:int, key:str , text:str, pos_key='left'):
    """ Function doc """
    size = os.get_terminal_size(2)[0] - marg                        # text chunks size
    ltext = [text[i:i+size] for i in range(0, len(text), size)]     # chunks
    smarg = '\n' + marg * ' '                                       # delimiter for chunks

    if pos_key == 'top':
        return f"{smarg[1:]}{key}{smarg}{smarg.join(line.lstrip(' ') for line in ltext)}"
    else:
        return f"{key:<{marg}}{smarg.join(line.lstrip(' ') for line in ltext)}"
