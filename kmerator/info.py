# info.py

"""
Genral informations on to the program.
"""

APPNAME = "kmerator"
VERSION = "0.9.3-beta"
VERSION_DATASET = 1
SHORTDESC = "Find specific gene or transcript kmers. And more."
LICENCE = "GPL3"
AUTHOR = 'Sébastien RIQUIER, IRMB, Montpellier'
AUTHOR_EMAIL = "sebastien.riquier@ucd.ie"
CONTIBUTORS = [
    'Chloe BESSIERE chloe.bessiere@inserm.fr>'
    'Benoit GUIBERT <benoit.guibert@inserm.fr>',
]
DOC = f"""
-------------------------------------------------
                                   _
 ____  __.                        | |
|    |/ _| ___ _ _   __  ___  __ _| |_ ___  _ __
|      <  |  _` ` \ / _ \ '__/ _` | __/ _ \| '__|
|    |  \ | | | | ||  __/ |   (_| | || (_) | |
|____|__ \|_| |_| |_\___|_|  \__,_|\__\___/|_|
        \/
 Version: v{VERSION}
 Dependencies:
   - Jellyfish >= v2.0
-------------------------------------------------
"""

EXAMPLES = f"""
------------------------------------------------------
EXAMPLES:

Before all, remember that kmerator needs a jellyfish index of the genome.

Good idea before requests kmerator:
 kmerator -e             # Edit config file to set default options

Some requests:
 kmerator -s npm1 braf           # get specific kmers from NPM1 and BRAF genes
 kmerator -s genes.txt           # you can alse use a file with gene list (#: comment)
 kmerator -f file.fa             # give a fasta file from unannotated sequences

Maintains yours kmerator indexes
 kmerator -l                     # list local avalaible indexes
 kmerator --mk-dataset -r 101    # install dataset for release 101
 kmerator -u -S zebrafish        # update dataset if new release avalaible

Get info on some genes/transcript
 kmerator --info MLLT3 Braf      # info about the MLLT3 and BRAF genes
 kmerator --info MMLT3 --all     # extended info about the MMLT3 gene
 kmerator --info genes.txt       # genes/transcripts are in a file
"""
