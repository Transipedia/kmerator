#!/usr/bin/env python3

'''
- if 'last' is given as release (default), get last release of transcriptome and it set to args.release
- check if dataset is present for this release (3 files)
- if not present, ask to build it
  - download some files from Ensembl database (mysql) and build 'geneinfo.pkl' file
  - download transcriptome (cdna + ncran), filter it and build a jellyfish index
  - write transcriptome as pickle file

--------
Ensembl links
  - http://www.ensembl.org/info/docs/api/core/core_schema.html
  - https://ftp.ensembl.org/pub/current-mysql/
  - http://www.ensembl.org/info/docs/api/core/diagrams/Core.svg
Other links:
  - https://lists.ensembl.org/pipermail/dev_ensembl.org/2013-January/003357.html
'''


import sys
import os
import argparse
import requests
from bs4 import BeautifulSoup
import gzip
import pickle
import threading
import subprocess

import info
from color import *
from mk_geneinfo import GeneInfoBuilder
from mk_transcriptome import TranscriptomeBuilder
import exit


species = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    'zebrafish': 'danio_rerio',
    "horse": "equus_caballus",
    "hen" : "gallus_gallus",
    "c.elegans": "caenorhabditis_elegans",
    "droso": "drosophila_melanogaster",
    # "e.coli": "escherichia_coli", # not in Ensembl
}


def main():
    """ Function doc """
    args = usage()
    dataset = Dataset(args)
    if args.list_dataset:
        dataset.list()
    elif args.update_last:
        dataset.update_last()
    elif args.rm_dataset:
        dataset.remove()
    elif args.mk_dataset:
        dataset.make()
    elif args.last_avail:
        dataset.last_available()
    elif args.load:
        geneinfo_dict = dataset.load_geneinfo()
        transcriptome_dict = dataset.load_transcriptome()
        for type in ('gene', 'symbol', 'alias', 'transcript'):
            print(type, next(iter(geneinfo_dict[type].items())))
        print(next(iter(transcriptome_dict.items())))
        print("Assembly in geneinfo:", geneinfo_dict['assembly'])
        print("Assembly in args.assembly:", args.assembly)
        print("Relesase", args.release)


class Dataset:
    """
    Methods:
      - __init__(): if release == 'last', redefine release with last release number
      - is_present(): Ensure release files are presents
      - build(): build files for the release
      - load_geneinfo(): load geneinfo file as dict
    """

    def __init__(self, args):
        """ Class initialiser """
        self.args = args
        self.base_url = "http://ftp.ensembl.org/pub"
        self.attended = ['transcriptome.pkl', 'transcriptome.jf', 'geneinfo.pkl', 'report.md']
        # ~ self.assembly = None
        self.transcriptome_fa = None         # transcriptome fasta path
        self.transcriptome_pkl = None        # transcriptome pickle path
        self.transcriptome_jf = None         # transcriptome jellyfish path
        self.geneinfo_pkl = None             # geneinfo path
        self.report_md = None                # report path
        self.report = []                     # report
        self.ebl_releases = []               # all releases avalaible on Ensembl
        self.dataset = self.set_dataset_dict()   # {complete: {specie: {release:[k31,...]}},partial:..., other:...}
        ### Is args.specie an alternative name ?
        if self.args.specie.lower() in species:
            self.args.specie =  species[self.args.specie.lower()]
        if self.args.release == 'last' and not self.args.list_dataset:
            self.args.release = self.get_ebl_current_release()
            self.get_ebl_releases()
        ### check if dataset is locally present and assign variable for each file
        if args.datadir and not self.args.list_dataset:
            self.dataset_ok = self.dataset_here()

    
    def set_dataset_dict(self):
        ### general variables
        files = next(os.walk(self.args.datadir))[2]
        attended = ['transcriptome.pkl', 'transcriptome.jf', 'geneinfo.pkl', 'report.md']
        dataset = {"complete": {},  # set complete 
                   "partial": {},   # set uncomplete, miss some files
                   "other": [],}     # files not part of dataset
        releases = {}
        
        ### dict of releases found: {<specie>: {<release>:[<file1>, <file2>]} }
        for file in files:
            try:
                file_l = file.split('.')
                specie = ".".join(file_l[:2])
                release = int(file_l[2])
                releases.setdefault(specie, {}).setdefault(release, []).append(file)
            except (IndexError, ValueError):
                dataset["other"].append(file)

        ### classify the releases (complete or incomplete)
        for specie in releases:
            for release, file_list in releases[specie].items():
                founded_items = set(".".join(f.split(".")[-2:]) for f in file_list)
                founded_files = list(map(lambda v: v in founded_items, self.attended))
                if all(founded_files):
                    dataset['complete'].setdefault(specie, {})[int(release)] = []
                else:
                    dataset['partial'].setdefault(specie, []).append(int(release))

        ### find for different version of jellyfish indexes of transcriptome
        for specie, rels in dataset['complete'].items():
            for rel in rels:
                for file in releases[specie][rel]:
                    if file.endswith('.transcriptome.jf'):
                        k = file.split('.')[3] if file.split('.')[3].startswith('k') else 'k?'
                        dataset['complete'][specie][rel].append(k)
        
        return dataset


    def get_ebl_current_release(self):
        """get number of last available release on Ensembl"""
        url = os.path.join(self.base_url, "current_mysql")
        try:
            r = requests.get(url)
        except requests.exceptions.ConnectionError:
            print(f"{YELLOW} Error connecting to Ensembl.{ENDCOL}")
            return None
        if not r.ok:
            print(r)
            return None
        soup = BeautifulSoup(r.text, 'lxml')
        core = [ a.text for a in soup.findAll('a') if a.text.startswith(f"{self.args.specie}_core")]
        try:
            ebl_current_release = core[0].split('_')[-2]
        except IndexError:
            print(f"{YELLOW}ErrorSpecie: {self.args.specie} not found in Ensembl database.{ENDCOL}")
            exit.gracefully(self.args)
        return ebl_current_release


    def get_ebl_releases(self):
        """ get avalaible releases on Ensembl, limited to 90 """
        try:
            r = requests.get(self.base_url)
        except requests.exceptions.ConnectionError:
            ### If no connection to Ensembl, suggest to use last local Release
            print(f"{YELLOW} Error connecting to Ensembl.{ENDCOL}")
            local_releases = self.get_local_releases()
            if not local_releases:
                print(f"{YELLOW} No dataset found in {self.args.datadir!r}, exit.")
                exit.gracefully(self.args)
            ask = 'y'
            if not self.args.yes:
                ask = input(f" {YELLOW}Last local Release found in datasets is {max(local_releases)}. Continue with it? (Yn): {ENDCOL}") or 'y'
            if ask.lower() == 'y':
                self.ebl_releases = local_releases  # this is a lie :)
                return
            else:
                print("Aborted by user.")
                exit.gracefully(self.args)

        ### If an error occur connecting to Ensembl, use the last local release avalaible.
        if not r.ok:
            all_files = next(os.walk(self.args.datadir))[2]
            self.ebl_releases = {int(i.split('.')[2]) for i in all_files if i.split('.')[0] == self.args.specie}
            ask = 'y'
            if not self.args.yes:
                print(f"{YELLOW} Error connecting to Ensembl (code {r.status_code}).")
                if not self.ebl_releases:
                    print(f" No dataset found in {self.args.datadir!r}, exit.{ENDCOL}")
                    exit.gracefully(self.args)
                ask = input(f" {YELLOW}Do you want to continue with release {str(max(self.ebl_releases))!r}? (Yn): {ENDCOL}") or 'y'
            if ask.lower() == 'y':
                return
            else:
                print("Aborted by user.")
                exit.gracefully(self.args)

        soup = BeautifulSoup(r.text, 'lxml')
        self.ebl_releases = [int(a.text.split('-')[1].rstrip('/')) for a in soup.findAll('a') if a.text.startswith('release-')]
        self.ebl_releases = [n for n in self.ebl_releases if n >= 90]


    def get_local_releases(self):
        print("get_local_releases()", {int(i.split('.')[2]) for i in next(os.walk(self.args.datadir))[2] if i.split('.')[0] == self.args.specie} )
        return {int(i.split('.')[2]) for i in next(os.walk(self.args.datadir))[2] if i.split('.')[0] == self.args.specie}


    def load_geneinfo(self):
        """ Load geneinfo.pkl as dict"""
        ### if dataset not found, ask to install
        if not self.dataset_ok:
            self.build()
            self.dataset_ok = self.dataset_here()
        ### Load geneinfo
        with open(self.geneinfo_pkl, 'rb') as fic:
            geneinfo_dict = pickle.load(fic)
            self.args.assembly = geneinfo_dict['assembly']
        ### check version of dataset
        version_geneinfo = geneinfo_dict.get('version', 0)
        if version_geneinfo < info.VERSION_DATASET:
            print(f"{RED}\n ErrorVersion: application and dataset versions are not compatible.\n"
                   "               Please, reinstall the dataset:\n"
                  f"{YELLOW} {info.APPNAME} --mk-dataset -S {self.args.specie} -r {self.args.release}\n")
            exit.gracefully(self.args)
        return geneinfo_dict


    def load_transcriptome(self):
        """ Load transcriptome.pkl as dict"""
        ### if dataset not found, ask to install
        if not self.dataset_ok:
            self.build()                                # build dataset
            self.dataset = self.set_dataset_dict()      # re-compute datasets
            self.dataset_ok = self.dataset_here()
        ### Load transcriptome
        with open(self.transcriptome_pkl, 'rb') as fic:
            transcriptome_dict = pickle.load(fic)
        return transcriptome_dict


    def dataset_here(self):
        """
        - Check if dataset is present in datadir
        - assign files names to matching variables
        - jump to rename_jf() if the releases jellyfish file names are in a older version (without kmer size)
         """
        def is_release_found(self):
            ### check if the specified release is present in the datasets
            releases_found = []
            assembly = None
            for specie, releases in self.dataset['complete'].items():
                if specie.startswith(self.args.specie):
                    assembly = specie.split('.')[1]
                    releases_found = releases.get(int(self.args.release), [])
                    if f"k{self.args.kmer_length}" in releases_found: 
                        return releases_found, assembly, True
            else:
                return releases_found, assembly, False        
        releases_found, assembly, found = is_release_found(self)
        
        ### Handle older file format (version < 2.3.0)
        if 'k?' in releases_found:
            rename = input(
                        f"{YELLOW}The dataset has an unspecified kmer value in the jellyfish file "
                        "name, because it was created with an older version of kmerator. Kmerator "
                        f"can rename the files to the new format before continuing [Yn]:{ENDCOL} "
                        ).lower() or 'y'
            if rename == 'y':
                self.rename_jf()
                releases_found, assembly, found = is_release_found(self)
                
        ### Assign name to file attributes
        if found:
            basename = f"{self.args.specie}.{assembly}.{self.args.release}"
            pathbasename = os.path.join(self.args.datadir, basename)
            self.transcriptome_fa = f"{pathbasename}.transcriptome.fa"
            self.transcriptome_pkl = f"{pathbasename}.transcriptome.pkl" 
            self.transcriptome_jf = f"{pathbasename}.k{self.args.kmer_length}.transcriptome.jf"
            self.geneinfo_pkl = f"{pathbasename}.geneinfo.pkl"
            self.report_md = f"{pathbasename}.report.md"
            
        return found


    def define_dataset(self):
        """
        Define the names of the local files for the specified dataset
        """

        ### Request Ensembl to get info (assembly name, release number, list of releases)
        if not self.ebl_releases:
            self.get_ebl_releases()
        if not int(self.args.release) in self.ebl_releases or not self.args.release.isnumeric():
            print(f"{ERROR} Error: {self.args.release!r} is not a valid Release (valid releases range from 90 to {max(self.ebl_releases)}).")
            exit.gracefully(self.args)

        ### Get assembly name
        url = f"{self.base_url}/release-{self.args.release}/fasta/{self.args.specie}/cdna/"
        r = requests.get(url)
        if not r.ok:
            print(f"{ERROR} Error: not a valid url: {url!r}.\n check for release ({self.args.release!r}) and specie ({self.args.specie!r}).{ENDCOL}")
            exit.gracefully(self.args)
        soup = BeautifulSoup(r.text, 'lxml')
        self.args.assembly = [a.text.split('.')[1] for a in soup.findAll('a') if a.text.startswith(self.args.specie.capitalize())][0]

        ### Assign files names of dataset
        basename = f"{self.args.specie}.{self.args.assembly}.{self.args.release}"
        self.geneinfo_pkl = os.path.join(self.args.datadir, f"{basename}.geneinfo.pkl")
        self.transcriptome_pkl = os.path.join(self.args.datadir, f"{basename}.transcriptome.pkl")
        self.transcriptome_jf = os.path.join(self.args.datadir, f"{basename}.k{self.args.kmer_length}.transcriptome.jf")
        self.report_md = os.path.join(self.args.datadir, f"{basename}.report.md")


    def build(self):
        """
        if dataset of this release is not present in the local directory of kmerator (datadir),
        we must build him, downloading some files and rearange them.
        """
        ### Ask user to download files
        valid = 'y' if self.args.yes else input(f"Dataset for release {self.args.release} ({self.args.specie}) not found, install it? (Yn) ")
        if valid.lower() in ['n', 'no']:
            print("Exited by user.")
            exit.gracefully(self.args)

        ### Check target directory permissions
        if not os.access(self.args.datadir, os.W_OK):
            print(f"{RED}\n Error: write acces denied to datadir ({self.args.datadir}).{ENDCOL}")
            exit.gracefully(self.args)

        ### Define target files
        if not self.dataset_ok:
            self.define_dataset()

        ### build kmerator dataset for the specie/release specified (multithreaded)
        geneinfo = GeneInfoBuilder(self.args, self.base_url, self.geneinfo_pkl, self.report)
        geneinfo.get_meta()         # first step : get meta info
        chr = geneinfo.data['chr']

        # ~ geneinfo_dict = geneinfo.build()
        # ~ TranscriptomeBuilder(self.args, self.base_url, self.transcriptome_jf, self.report, chr)

        ### build kmerator dataset for the specie/release specified (multithreaded)
        th1 = threading.Thread(target=TranscriptomeBuilder,
                               args=(self.args, self.base_url, self.transcriptome_jf, self.report, chr))
        th2 = threading.Thread(target=geneinfo.build)
        th1.start()
        th2.start()
        th1.join()
        th2.join()
        

        ### write report
        with open(self.report_md, 'w') as fh:
            fh.write(f"# Kmerator files for {self.args.specie}, release {self.args.release}\n")
            fh.write('\n'.join(self.report))


    def list(self):
        """ List local releases """
        ### Show releases
        print(f"\n {YELLOW}Location of datasets:{ENDCOL} {self.args.datadir}")
        if self.dataset['complete']:
            to_print = f"\n {YELLOW}Datasets found:{ENDCOL}\n"
            for specie, releases in self.dataset['complete'].items():
                rels = sorted([(int(k),v) for k,v in releases.items()])
                for rel in rels:
                    to_print += f"  - {specie.split('.')[0]}: {rel[0]} ({', '.join(sorted(rel[1]))})\n"
            print(to_print)
        else:
            print(f"\n {YELLOW}No releases found{ENDCOL}\n")
        if self.dataset['partial']:
            print(f" {YELLOW}Incompletes datasets:{ENDCOL}")
            for specie, releases in self.dataset['partial'].items():
                print(f"  - {specie.split('.')[0]}: {', '.join([str(i) for i in releases])}{ENDCOL}")    
        print()

        ### Other file found in dataset location
        if self.dataset['other']:
            print(f" {YELLOW}Files not part of a dataset:{ENDCOL}")
            print(*[f"  - {i}" for i in self.dataset['other']], sep="\n")
        print()
        
        ### Find for dataset built by old version of kmerator, and propose to update the datasets
        found = False
        rename = 'n'
        for species, releases in self.dataset['complete'].items():
            if found: break
            for release, values in releases.items():
                if 'k?' in values:
                    found = True
                    rename = input(f"{YELLOW}\nSome releases has a kmer size not determined (k?) "
                            "due to previous version of kmerator. Do you want to automatically "
                            f"correct their size? [Yn]: {ENDCOL}").lower().strip() or 'y'
                    break
        if rename == 'y':
            self.rename_jf()

        ### exit
        exit.gracefully(self.args)


    def remove(self):
        """ remove a release """
        ### list dataset files this specie/release
        release_files = []
        jellyfish_files = []
        for file in next(os.walk(self.args.datadir))[2]:
            l_file = file.split('.')
            try:
                specie, assembly, release = l_file[:3]
                if self.args.specie == specie and self.args.release == release:
                    release_files.append(file)
                    if l_file[-1] == 'jf' and l_file[-2] == 'transcriptome':
                        jellyfish_files.append(file)
            except ValueError:
                continue
        ### if release empty
        if not release_files:
            print(f"Dataset not found for {self.args.specie!r}, release {self.args.release!r}.")
            exit.gracefully(self.args)
        ### if multiple jellyfish for the release, remove only the transcriptome file
        if len(jellyfish_files) > 1:
            for file in jellyfish_files:
                k = file.split('.')[3]
                if k.startswith('k'):
                    if int(k[1:]) == self.args.kmer_length:
                        release_files = [file]
                        break
                elif k == 'transcriptome':
                    release_files = [file]

        ### if not files to delete
        if not release_files:
            print(f"Dataset not found for {self.args.specie!r}, release {self.args.release!r}.")
            exit.gracefully(self.args)
        ### Ask to remove
        resp = 'y'
        if not self.args.yes:
            print(f"\nspecie: {self.args.specie} - release: {self.args.release}\n")
            # ~ print("Files to delete", *to_delete, sep='\n - ')
            print("Files to delete", *release_files, sep='\n - ')
            resp = input("\nDelete files (Yn): ") or 'y'
        ### Remove
        if not resp.lower() == 'n':
            for file in release_files:
                os.remove(os.path.join(self.args.datadir, file))
        else:
            print("\nAborted by user.  ")
        exit.gracefully(self.args)


    def make(self):
        """ Download files and make dataset """
        resp = 'y'
        if self.dataset_ok and not self.args.yes:
            resp = input(f"Release {self.args.release} already exists, erase? (Yn): ") or 'y'
        if resp.lower() == 'y':
            # ~ self.args.release = str(release)           # because build() using self.args.release
            self.args.yes = True
            print(f" 🧬 Build kmerator dataset for {self.args.specie}, release {self.args.release}, please wait...")
            self.build()

        exit.gracefully(self.args)


    def update_last(self):
        """ Function doc """
        if not self.ebl_releases:
            self.get_ebl_releases()
        self.args.release = self.get_ebl_current_release()

        if not self.dataset_here():
            ask = 'y'
            if not self.args.yes:
                ask = input(f"Dataset for {self.args.specie} will be updated to release {self.args.release}, continue ? (Yn): ") or 'y'
            if ask.lower() == 'y':
                ### Build Dataset
                self.make()
            else:
                sys.exit("Aborted by user.")
        else:
            print(f"The last release for {self.args.specie} is {self.args.release}, nothing to do.")
        exit.gracefully(self.args)


    def last_available(self):
        print(self.get_ebl_current_release())
        exit.gracefully(self.args)


    def rename_jf(self):
        """
        Migrate dataset to new version: 
        Before version 2.3.0, the jellyfish index file name did not containt the size of the kmers.
        rename_jf aims to rename file, adding kmer size in the name, according the 2.3.0 version
        """
        files = next(os.walk(self.args.datadir))[2]
        ### get k value
        get_k = lambda res,i: res[i + 1] 
        # ~ print(files)
        for file in files:
            if file.endswith(".transcriptome.jf"):
                file_l = file.split(".")
                k_expected = file_l[3]
                if len(file_l) == 5 and k_expected == "transcriptome":
                    ### jellyfish keep the command line used to build the jf file.
                    cmd_get = f"jellyfish info {os.path.join(self.args.datadir, file)} | head -1"
                    res = subprocess.check_output(cmd_get, shell=True, text=True).split(' ')
                    for opt in ('-m', '--mer-len'):
                        k = f"k{get_k(res, res.index(opt))}"
                        break
                    file_l.insert(3, k)
                    new_name = '.'.join(file_l)
                    cmd_set = f"mv {os.path.join(self.args.datadir, file)} {os.path.join(self.args.datadir, new_name)}"
                    try:
                        subprocess.check_output(cmd_set, shell=True)
                    except Exception as e:
                        sys.exit(e)
        self.dataset = self.set_dataset_dict()          


def usage():
    parser = argparse.ArgumentParser()
    exclusive = parser.add_mutually_exclusive_group(required=True)
    ### OPTION
    parser.add_argument('-S', '--specie',
                        help=(
                            "indicate a specie referenced in Ensembl, to help, follow the link "
                            "https://rest.ensembl.org/documentation/info/species. You can use "
                            "the 'name', the 'display_name' or any 'aliases'. For example human, "
                            "homo_sapiens or homsap are valid (default: human)."
                            ),
                        default='human',
                        )
    parser.add_argument('-r', '--release',
                        help="release of transcriptome (default: last).",
                        default="last",
                        )
    exclusive.add_argument('-d', '--datadir',
                        help=(
                            "Directory where kmerator file are stored. Files are:"
                            "\n - fasta file of modified transcriptome (fa)"
                            "\n - jellyfish of this transcriptome (jf)"
                            "\n - metadata file like gene-symbols or aliases, "
                            ),
                        # ~ required = True,
                        )
    parser.add_argument('-l', '--list-dataset', '--datasets',
                        action="store_true",
                        help="list local releases",
                       )
    parser.add_argument('--debug',
                        action="store_true",
                        help="Show more",
                       )
    parser.add_argument('--load',
                        action="store_true",
                        help="load dataset",
                       )
    parser.add_argument('--tmpdir',
                        help="temporary dir",
                        default="/tmp",
                       )
    parser.add_argument('--rm-dataset',
                        action="store_true",
                        help="remove a dataset, according with --specie and --release options",
                       )
    parser.add_argument('--mk-dataset',
                        action="store_true",
                        help="make a dataset, according with --specie and --release options",
                       )
    parser.add_argument('-t', '--thread',
                        type=int,
                        help="thread number (default: 1)",
                        default=1,
                       )
    parser.add_argument('-u', '--update-last',
                        action="store_true",
                        help="builds a new dataset if a new version is found on Ensembl",
                       )
    parser.add_argument('-k', '--kmer-length',
                        type=int,
                        help="kmer length (default: 31)",
                        default=31,
                       )
    exclusive.add_argument('--last-avail', '--last-available',
                        action='store_true',
                        help="last release available on Ensembl",
                       )
    parser.add_argument('--keep',
                        action='store_true',
                        help=("keep kmerator transcriptome as fasta format."
                            ),
                        )
    parser.add_argument('-y', '--yes',
                        action='store_true',
                        help=("assumes 'yes' as the prompt answer, run non-interactively."),
                        )
    return parser.parse_args()


if __name__ == "__main__":
    main()
