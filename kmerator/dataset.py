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
  - https://ftp.ensembl.org/pub/release-108/mysql/
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
    elif args.load:
        geneinfo_dict, transcriptome_dict = dataset.load()
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
        self.assembly = None
        self.ebl_releases = None
        self.transcriptome_fa = None    # transcriptome fasta path
        self.transcriptome_pkl = None   # transcriptome pickle path
        self.transcriptome_jf = None    # transcriptome jellyfish path
        self.geneinfo_pkl = None        # geneinfo path
        self.report_md = None           # report path
        self.report = []                # report
        self.ebl_releases = []          # all releases avalaible on Ensembl
        ### Is args.specie is an alternative name ?
        if self.args.specie.lower() in species:
            self.args.specie =  species[self.args.specie.lower()]
        if self.args.release == 'last' and not self.args.list_dataset:
            self.get_ebl_releases()
            self.args.release = str(max(self.ebl_releases))
        ### check if dataset is locally present and assign variable for each file
        self.dataset_ok = self.dataset_here()


    def get_ebl_releases(self):
        """ get avalable releases on Ensembl, limited to 90 """
        r = requests.get(self.base_url)
        ### If an error occur connecting to Ensembl, use the last local release avalaible.
        if not r.ok:
            all_files = next(os.walk(self.args.datadir))[2]
            self.ebl_releases = {int(i.split('.')[2]) for i in all_files if i.split('.')[0] == self.args.specie}
            ask = 'y'
            if not self.args.yes:
                print(f"{YELLOW} Error connecting to Ensembl (code {r.status_code}).")
                if not self.ebl_releases:
                    print(f" No dataset found in {self.args.datadir!r}, exit")
                    exit.gracefully(args)
                ask = input(f" Do you want to continue with release {str(max(self.ebl_releases))!r}? (Yn): {ENDCOL}") or 'y'
            if ask.lower() == 'y':
                return
            else:
                print("Aborted by user.")
                exit.gracefully(self.args)

        soup = BeautifulSoup(r.text, 'lxml')
        self.ebl_releases = [int(a.text.split('-')[1].rstrip('/')) for a in soup.findAll('a') if a.text.startswith('release-')]
        self.ebl_releases = [n for n in self.ebl_releases if n >= 90]


    def get_local_releases(self):
        return {int(i.split('.')[2]) for i in next(os.walk(self.args.datadir))[2] if i.split('.')[0] == self.args.specie}


    def load(self):
        """ Load dataset as dict (geneinfo and transcriptome)"""
        ### if dataset not found, ask to install
        if not self.dataset_ok:
            self.build()
        ### Load geneinfo and transcriptome
        with open(self.geneinfo_pkl, 'rb') as fic:
            geneinfo_dict = pickle.load(fic)
            self.args.assembly = geneinfo_dict['assembly']
        with open(self.transcriptome_pkl, 'rb') as fic:
            transcriptome_dict = pickle.load(fic)
        return geneinfo_dict, transcriptome_dict


    def dataset_here(self):
        """
        - Check if dataset is present in datadir
        - assign files names to matching variables
        """
        attended  = ['transcriptome_pkl', 'transcriptome_jf', 'geneinfo_pkl', 'report_md']
        found = 0
        all_files = next(os.walk(self.args.datadir))[2]
        ### dict of releases : {specie: {release:[file1, file2]} }
        for file in all_files:
            try:
                file = file.split('.')
                specie = file[0]
                assembly = file[1]
                release = file[2]
                item = f"{file[3]}_{file[4]}"
                ext = f"{file[3]}.{file[4]}"
            except IndexError:
                continue
            if specie == self.args.specie and release == self.args.release and item in attended:
                found += 1
                basename = f"{specie}.{assembly}.{release}.{ext}"
                setattr(self, item, os.path.join(self.args.datadir, basename))
        return True if found == len(attended) else False


    def define_dataset(self):
        """
        Define the names of the local files for the dataset
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
        self.transcriptome_jf = os.path.join(self.args.datadir, f"{basename}.transcriptome.jf")
        self.report_md = os.path.join(self.args.datadir, f"{basename}.report.md")


    def build(self):
        """
        if dataset of this release is not present in the local directory of kmerator (datadir),
        we must build him, downloading some files and rearange them.
        """
        ### Ask user to download files
        valid = 'y' if self.args.yes else input(f"Dataset for release {self.args.release} ({self.args.specie}) not found, intall it? (Yn) ")
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
        print(f" 🧬 Build kmerator dataset for {self.args.specie}, release {self.args.release}, please wait...")
        th1 = threading.Thread(target=TranscriptomeBuilder, args=(self.args, self.base_url, self.transcriptome_jf, self.report))
        th2 = threading.Thread(target=GeneInfoBuilder, args=(self.args, self.base_url, self.geneinfo_pkl, self.report))
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
        attended  = ['transcriptome.pkl', 'transcriptome.jf', 'geneinfo.pkl', 'report.md']
        releases = {}
        files = next(os.walk(self.args.datadir))[2]

        ### dict of releases : {<specie>: {<release>:[<file1>, <file2>]} }
        for file in files:
            try:
                file = file.split('.')
                specie = file[0]
                release = file[2]
                item = f"{file[3]}.{file[4]}"
                releases.setdefault(specie, {}).setdefault(release, []).append(item)
            except IndexError:
                print(f"{DEBUG} Error: {'.'.join(file)!r} improperly formated.{ENDCOL}")

        ### classify the releases (complete or incomplete)
        releases_ok = {}
        releases_ko = {}
        for specie in releases:
            for release in releases[specie]:
                if all(map(lambda v: v in releases[specie][release], attended )):
                    releases_ok.setdefault(specie, []).append(release)
                else:
                    releases_ko.setdefault(specie, []).append(release)

        if releases_ko:
            print(f"\n{YELLOW} Incompletes datasets:")
            for specie, releases in releases_ko.items():
                print(f"{YELLOW}   {specie}: {', '.join(releases)}.{ENDCOL}")

        ### Show releases
        print(f"\n Location of datasets: {self.args.datadir}")
        if releases_ok:
            print("\n Datasets found:", *[f"{k}: {', '.join([str(a) for a in sorted([int(i) for i in v])])}" for k,v in releases_ok.items()], sep="\n  - ")
            if releases_ko:
                print("\n Incompletes datasets:", *[f"{k}: {', '.join([str(a) for a in sorted([int(i) for i in v])])}" for k,v in releases_ko.items()], sep="\n  - ")
        else:
            print(f"\n No releases found")
        print()

        ### exit
        exit.gracefully(self.args)


    def remove(self):
        """ remove a release """
        ### list dataset file for this specie/release
        to_delete = []
        for file in next(os.walk(self.args.datadir))[2]:
            l_file = file.split('.')
            try:
                specie, assembly, release = l_file[:3]
                if self.args.specie == specie and self.args.release == release:
                    to_delete.append(file)
            except ValueError:
                continue
        ### if not files to delete
        if not to_delete:
            print(f"Dataset not found for {self.args.specie!r}, release {self.args.release!r}.")
            exit.gracefully(self.args)
        ### Ask to remove
        resp = 'y'
        if not self.args.yes:
            print(f"\nspecie: {self.args.specie} - release: {self.args.release}\n")
            print("Files to delete", *to_delete, sep='\n - ')
            resp = input("\nDelete files (Yn): ") or 'y'
        ### Remove
        if not resp.lower() == 'n':
            for file in to_delete:
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
            print(f"Install {self.args.specie}, releaase {self.args.release}, please wait...")
            self.build()

        exit.gracefully(self.args)


    def update_last(self):
        """ Function doc """
        if not self.ebl_releases:
            self.get_ebl_releases()
        last_ebl_release = max(self.ebl_releases)
        self.args.release = str(last_ebl_release)

        if not self.dataset_here():
            ask = 'y'
            if not self.args.yes:
                ask = input(f"Dataset for {self.args.specie} will be updated to release {last_ebl_release}, continue ? (Yn): ") or 'y'
            if ask.lower() == 'y':
                ### Build Dataset
                self.args.release = str(last_ebl_release)
                self.make()
            else:
                sys.exit("Aborted by user.")
        else:
            print(f"The last release for {self.args.specie} is {last_ebl_release}, nothing to do.")
        exit.gracefully(self.args)


def usage():
    parser = argparse.ArgumentParser()
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
    parser.add_argument('-d', '--datadir',
                        help=(
                            "Directory where kmerator file are stored. Files are:"
                            "\n - fasta file of modified transcriptome (fa)"
                            "\n - jellyfish of this transcriptome (jf)"
                            "\n - metadata file like gene-symbols or aliases, "
                            ),
                        required = True,
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