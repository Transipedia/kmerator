# info.py

"""
Genral informations on to the program.
"""

APPNAME = "kmerator"
SHORTDESC = "Find specific gene or transcript kmers. And more."
LICENCE = "GPL3"
VERSION = "0.8.2-beta"
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
 kmerator -s npm1        # get specific kmers form NPM1 gene
 kmerator -s genes.txt   # you can also use a file with gene list
 kmerator -f file.fa     # give a fasta file fr unannotated sequences

Maintains yours kmerator indexes
 kmerator -l                   # list local avalaible indexes
 kmerator --mk-dataset -r 101  # install dataset for release 101
 kmerator -u -S zebrafish      # update dataset if new release avalaible
"""
