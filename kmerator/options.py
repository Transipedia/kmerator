# usage.py

'''
TODO
 - ajouter une option permettant de faire un fasta d'une release de transcriptome
'''

"""
Options and arguments to bo passed on to the program.
"""

import os
import sys
import argparse
import info
import tempfile

from color import *


class EditConfig(argparse.Action):
    """
    https://stackoverflow.com/questions/8632354/python-argparse-custom-actions-with-additional-arguments-passed
    https://stackoverflow.com/questions/11001678/argparse-custom-action-with-no-argument
    """
    def __init__(self, option_strings, conf, *args, **kwargs):
        """ Function doc """
        super(EditConfig, self).__init__(option_strings=option_strings, *args, **kwargs)
        self.conf = conf

    def __call__(self, parser, namespace, values, option_strings):
        self.conf.edit()
        sys.exit()


def usage(conf):
    """
        Help function with argument parser.
    """
    parser = argparse.ArgumentParser(
        description=info.DOC,
        epilog = info.EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    exclusive = parser.add_mutually_exclusive_group(required=True)
    exclusive.add_argument('-s', '--selection',
                        help=(
                            "list of gene IDs (ENSG, gene Symbol or alias) or transcript IDs "
                            "(ENST) from which you want to extract specific kmers from. For genes, "
                            "kmerator search specific kmers along the gene. For transcripts, it "
                            "search specific kmers to the transcript. You can also give a file "
                            "with yours genes/transcripts separated by space, tab or newline. "
                            "If you want to use your own unannotated sequences, you "
                            "must give your fasta file with --fasta-file option."
                            ),
                        nargs='+',
                        )
    exclusive.add_argument('-f', '--fasta-file',
                        help=(
                            "Use this option when yours sequences are unannonated or provided "
                            "by a annotation file external from Ensembl. Otherwise, "
                            "use --selection option."
                            ),
                        )
    parser.add_argument('-d', '--datadir',
                        help=("Storage directory for kmerator datasets."
                              "We recommend to set this parameter by editing the configuration "
                              f"file ({info.APPNAME} --edit)"
                              ),
                        )
    parser.add_argument('-g', '--genome',
                        help=(
                            "Genome jellyfish index (.jf) to use for k-mers requests."
                            ),
                        required=True,
                        )
    parser.add_argument('-S', '--specie',  ### replace -a --appris in julia version
                        help=(
                            "indicate a specie referenced in Ensembl, to help, follow the link "
                            "https://rest.ensembl.org/documentation/info/species. You can use "
                            "the 'name', the 'display_name' or any 'alias'. For example human or "
                            "homo_sapiens are valid (default: human)."
                            ),
                        default='human',
                        )
    parser.add_argument('-k', '--kmer-length',
                        type=int,
                        help="k-mer length that you want to use (default 31).",
                        default=31,
                        )
    parser.add_argument('-r', '--release',
                        help="release of transcriptome (default: last).",
                        default="last",
                        )
    parser.add_argument('--chimera',
                        action='store_true',
                        help="Only with '--fasta-file' option.",
                        )
    parser.add_argument('--stringent',
                        action='store_true',
                        help=(
                            "Only for genes with '--selection' option: use this option if you want "
                            "to select gene-specific k-mers present in ALL transcripts for your "
                            "gene. If false, a k-mer is considered as gene-specific if present "
                            "in at least one isoform of your gene of interest."
                            ),
                        )
    parser.add_argument('--max-on-transcriptome',
                        type=float,
                        help=argparse.SUPPRESS,
                        default=0,
                        )
    parser.add_argument('-o', '--output',
                        help="output directory, created if not exists (default: 'output')",
                        default='output',
                        )
    parser.add_argument('-t', '--thread',
                        type=int,
                        help="run n process simultaneously (default: 1)",
                        default=1,
                        )
    parser.add_argument('--tmpdir',
                        help="directory to temporary file (default: /tmp/kmerator_<random>)",
                        default=None,
                        )
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        help="Show more details while Kmerator is running.",
                        )
    parser.add_argument('--keep',
                        action='store_true',
                        help=("keep intermediate files (sequences, indexes, separate kmers and "
                            "contigs files)."
                            ),
                        )
    parser.add_argument('-y', '--yes',
                        action='store_true',
                        help=("assumes 'yes' as the prompt answer, run non-interactively."),
                        )
    parser.add_argument('-e', '--edit-config',
                        action=EditConfig,
                        nargs=False,
                        conf=conf,
                        help='Edit config file',
                        )
    exclusive.add_argument('-l', '--list-dataset', '--list-datasets',
                        help=("list the local datasets (based on the datadir option)."),
                        action='store_true',
                        )
    exclusive.add_argument('--rm-dataset',
                        action='store_true',
                        help="remove a dataset, according with --specie and --release options",
                      )
    exclusive.add_argument('--mk-dataset',
                        action='store_true',
                        help="make a dataset, according with --specie and --release options",
                      )
    exclusive.add_argument('--last-avail', '--last-available',
                        action='store_true',
                        help="last release available on Ensembl",
                       )
    exclusive.add_argument('-u', '--update-dataset',
                        action="store_true",
                        help="builds a new dataset if a new version is found on Ensembl",
                       )
    parser.add_argument('-v', '--version',
                        action='version',
                        version=f'{parser.prog} v{info.VERSION}',
                        )

    ### Set default values define by Config class (using ConfigParser)
    default_conf = conf.content['CMD_ARGS']
    parser.set_defaults(**default_conf)
    ### Reset `required` attribute when provided from config file
    for action in parser._actions:
        if action.dest in default_conf:
            action.required = False
            ### handle booleans cause ConfigParse return a string.
            if action.__class__.__name__ == "_StoreTrueAction":
                action.default = True if action.default.lower() in ['true', 'yes', '1', 'on'] else False

    ### Go to 'usage()' without arguments or stdin
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    ### Check args
    checkup_args(args)

    return args


def checkup_args(args):
    """ checkup some options. """
    ### --datadir - at the first use, configure path of kmerator files
    if not args.datadir:
        mesg = ("\nAt the first use of kmerator, you need to configure a path for the files that "
               "kmerator will store: a fasta file (fa), a jellyfish file (jf) and a binary"
               " file (pkl) per release. Each release takes more than 2GB "
               "(2.2GB for release 108), so allow enough space.\n\n"
               f"Also, to avoid to set the '--genome' options each time you launch {info.APPNAME},"
               " you should specify the path of the jellyfish index of your reference genome.\n\n"
               "Note that other parameters are available in the configuration file, so that you "
               "don't have to enter some recurring options (such as the species "
               "for example).\n\n"
               f"Please, use the {BOLD}{DARKCYAN}--edit-config{ENDCOL} "
               f"(or {BOLD}{DARKCYAN}-e{ENDCOL}) option first to set the 'datadir' path parameter.\n")
        sys.exit(mesg)
    ### --datadir - check if path exists and right read is set
    if not os.path.isdir(args.datadir):
        sys.exit(f"{ERROR}Error: {args.datadir!r} not found or not a directory.{ENDCOL}")
    if not os.access(args.datadir, os.R_OK) or not os.access(args.datadir, os.X_OK):
        sys.exit(f"{ERROR}Error: insufficient rights in {args.datadir!r}.{ENDCOL}")
    ### --fasta-file - check if query fasta file is present
    if args.fasta_file and not os.path.isfile(args.fasta_file):
        sys.exit(f"{ERROR}Error: {args.fasta_file!r} not found.{ENDCOL}")
    ### --fasta-file - Check query fasta file extension
    if args.fasta_file and args.fasta_file.split('.')[-1] not in ("fa", "fasta"):
        sys.exit(f" {ERROR}Error: {os.path.basename(args.fasta_file)} " \
                f"Does not appears to be a fasta file.{ENDCOL}")
    ### --selection - Check Gene Select option
    if args.selection:
        for gene in args.selection:
            if gene.startswith('ENS') and '.' in gene:
                sys.exit(f"{ERROR}ENSEMBL annotations with version (point) like "
                        f"ENSTXXXXXXXX.1 is forbidden, just remove the '.1'."
                        )
        ### Define genes/transcripts provided when they are in a file
        if len(args.selection) == 1 and os.path.isfile(args.selection[0]):
            with open(args.selection[0]) as fh:
                args.selection = []
                for line in fh:
                    args.selection += line.split('#')[0].split()
    ### --genome - check jellifish genome
    # ~ if not args.genome[-3:] == '.jf':
        # ~ sys.exit(f"{ERROR}Error: file not a jellyfish index ({args.genome!r}).{ENDCOL}")
    ### --chimera level works only with --fasta-file option
    if args.chimera and not args.fasta_file:
        sys.exit(f"{ERROR}Error: '--chimera' needs '--fasta-file' option.{ENDCOL}")
    ### temporary directory
    if not args.keep:
        if not args.tmpdir:
            args.tmpdir = tempfile.mkdtemp(prefix="kmerator_")
        elif not os.path.isdir(args.tmpdir):
            sys.exit(f"{ERROR}Error: temporary directory not found ({args.tmpdir}).{ENDCOL}")
        else:
            args.tmpdir = tempfile.mkdtemp(prefix="kmerator_", dir=args.tmpdir)
    else:
        args.tmpdir = args.output


'''
def is_transcriptome_ensembl_file(first_row):
    """Checks if the line matches the Ensembl transcriptome file format"""
    line = first_row.split()
    try:
        if (not line[6].startswith('gene_symbol:') or
            not line[0].startswith('>ENST') or
            not line[3].startswith('gene:')):
            return False
    except:
        return False
    return True
'''
