
import os
import sys
import subprocess
from configparser import ConfigParser, ParsingError

import info
import color

DEFAULT_CONFIG = """
[CMD_ARGS]
## --datadir option
## path to kmerator datasets, one per release:
##   - modified transcriptome file (pkl)
##   - jellyfish transcriptome file (jf))
##   - metadata files (pkl)
##   - report file (md)
# datadir = /path/to/kmerator/directory

## --release option
##  release number (ex: 109) or last (last find for last release)
## default: last
# release = last

## --genome option
## Path of genome jellyfish index
# genome = /index/jellyfish/GRCh38_with_MT.jf

## --specie option
## set a specie (default: human)
# specie = human

## --thread option
## number of thread (default: 1)
# thread = 1

## --kmer-length option
## set kmer length (default: 31)
# kmer_length = 31

## --max-on-transcriptome
## Only with --fasta-file option
# --max-on-transcriptome = 0

## --max-on-genome
## Only with --fasta-file option
# --max-on-genome = 1

## --stringent option
## Only with --selection option
# stringent = False

## output directory
# output = ./output

## --yes option
## assumes 'yes' as the prompt answer, run non-interactively
## default: False
# yes = True

## --keep option
## keep intermediate files (default: False)
# keep = False
"""

class Config:
    """
    attributes:
        - filename: name of configuration file
        - filepath: path of the configuration file
        - content: dict of configuration file
    methods:
        - edit(): edit configuration file
    """
    filename = 'config-v3.ini'

    def __init__(self, appname):
        """ Class initialiser """
        ### define location of config file
        self.filePath = self.get_configPath(appname)
        ### create config file if not exists
        if not os.path.isfile(self.filePath):
            self._set_default()
        ### create a dict of config file
        self.content = ConfigParser()
        try:
            self.content.read(self.filePath)
        except  ParsingError as err:
            for line, text in err.errors:
                print(f"{color.RED}❗Syntax error at line {line} ({text})")
            ithem = 'it' if len(err.errors) == 1 else 'them'
            print(f"Please correct {ithem} in {err.source!r} file.")
            sys.exit()


    def get_configPath(self, appname):
        """
        Find to config
        """
        ### Define config file location (user or root)
        if os.geteuid():
            configdir = os.path.join(os.path.expanduser('~'), '.config', appname)
        else:
            configdir = os.path.join('/etc', appname)
        configfile = os.path.join(configdir, self.filename)
        return configfile


    def _set_default(self):
        try:
            os.makedirs(os.path.dirname(self.filePath), exist_ok=True)
            with open(self.filePath, 'w') as fs:
                fs.write(DEFAULT_CONFIG)
        except PermissionError as err:
            sys.exit(err.msg)


    def edit(self):
        """ Edit configuration file """
        ### define editor
        editor = 'editor'
        if os.environ.get('XDG_SESSION_TYPE') == 'x11':
            editor = 'xdg-open'
        elif 'EDITOR' in os.environ:
            editor = os.environ['EDITOR']
        ### launch config edit
        subprocess.call([editor, self.filePath])
