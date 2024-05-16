import os
import requests
from bs4 import BeautifulSoup
import gzip
import pickle

import info
from color import *


class GeneInfoBuilder:
    """
    get files on Ensembl (mysql):
      - gene.txt
      - transcript.txt
      - xref.txt  (= gene-symbols)
      - external-synonym.txt (= aliases)
    merge info from theses files, and write has pickle file
    -----
    Example of result whith HSTF1 (alias):

    assembly: GRCh38
    ensg:
        ENSG00000185122:
            gene-symbol: HSF1
            canonical: ENST00000528838
            nb_transcript: 14
            chr: 11
            start: 144291591
            end: 144314720
            strand: 1
            aliases:
                - HSTF1
            enst:
                id: ENST00000528838
                seq: ATCG...
            enst:
                id:ENST00000400780
                seq: GCTA...
            ...
    transcript:
        ENST00000528838: ENSG00000185122
        ENST00000400780: ENSG00000185122
    symbols:
        HSF1: ENSG00000185122
        FGF4: ENSG00000075388
    aliases:
        HSTF1:
            - ENSG00000185122
            - ENSG00000075388
    -----
    http://www.ensembl.org/info/docs/api/core/core_schema.html#gene
    https://ftp.ensembl.org/pub/release-106/mysql/homo_sapiens_core_106_38/
    (http://www.ensembl.org/info/docs/api/core/diagrams/Core.svg)
    """

    def __init__(self, args, url, outfile, report):
        """ Class initialiser """
        if args.debug: print(f"{DEBUG}Build geneinfo, please wait...{ENDCOL}")
        self.args = args
        self.outfile = outfile
        self.report = report
        ### internal variables
        self.chr_dict = {}                      # ID:name for chromosomes
        ### get base url of files to download
        url = f"{url}/release-{args.release}/mysql/"
        r = requests.get(url)
        soup = BeautifulSoup(r.text, 'lxml')
        dir = [a.text for a in soup.findAll('a') if a.text.startswith(f"{args.specie}_core_{args.release}")][-1]
        self.url = url + dir
        ### dicts (refs: content downloaded, data: dict of results)
        self.refs = {"gene": {}, "xref": {}, "synonym": {}, "transcript": {}}
        self.data = {"gene": {}, "symbol": {}, "alias": {}, "transcript": {}, "assembly": "",
                     "chr": [], "version": info.VERSION_DATASET}


    def get_meta(self):
        '''
        2. get seq_region.txt to make a list of region_id (chromosomes + MT)
        1. get Assembly version of the specie
        AIMS - filter alternates chromosomes (scaffold, CHR_)
          - $0 = region_id
          - $1 = chromosome (ex: 16, X)
        push chromosomes list in geneinfo_dict['chr']
        push Assembly version in geneinfo_dic['assembly']

        https://lists.ensembl.org/pipermail/dev_ensembl.org/2013-January/003357.html
        '''
        if __name__ == "__main__": print("Get meta info, please wait...")

        '''
        CHROMOSOME LIST & ASSEMBLY
        https://lists.ensembl.org/pipermail/dev_ensembl.org/2013-January/003357.html
        '''
        ### get attrib_type file and keep id when second field == "karyotype_rank"
        attrib_type_url = f"{self.url}attrib_type.txt.gz"
        attrib_type_r = requests.get(attrib_type_url)
        attrib_type_str = gzip.decompress(attrib_type_r.content).decode()
        attrib_type_id = None
        for line in attrib_type_str.splitlines():
            line = line.split()
            if line[1] == "karyotype_rank":
                attrib_type_id = line[0]
                break
        del attrib_type_r, attrib_type_str

        ### get seq_region_attrib file and keep list of id when second field == attrib_type_id
        seq_region_ids = []
        seq_region_attrib_url = f"{self.url}seq_region_attrib.txt.gz"
        seq_region_attrib_r = requests.get(seq_region_attrib_url)
        seq_region_attrib_str = gzip.decompress(seq_region_attrib_r.content).decode()
        for line in seq_region_attrib_str.splitlines():
            line = line.split()
            if line[1] == attrib_type_id:
                seq_region_ids.append(line[0])
        del seq_region_attrib_r, seq_region_attrib_str

        ### get seq_region file and keep list of id when second field == attrib_type_id
        seq_region_url = f"{self.url}seq_region.txt.gz"
        seq_region_r = requests.get(seq_region_url)
        seq_region_str = gzip.decompress(seq_region_r.content).decode()
        for line in seq_region_str.splitlines():
            line = line.split()
            if line[0] in seq_region_ids:
                self.chr_dict[line[0]] = line[1]
        del seq_region_r, seq_region_str
        ## Add chromosome list in geneinfo
        self.data['chr'] = list(self.chr_dict.values())

        ### get meta file and keep assembly.default
        meta_url = f"{self.url}meta.txt.gz"
        meta_r = requests.get(meta_url)
        meta_str = gzip.decompress(meta_r.content).decode()
        for line in meta_str.splitlines():
            line = line.split()
            if line[2] == "assembly.default":
                assembly = line[3].split('.')[0]
                self.data["assembly"] = assembly
                break
        del meta_r, meta_str


    def build(self):
        """ Function doc """
        ### Get ENSG
        self._set_refs()
        ### Add SYMBOLS
        self._add_symbols()
        ### Add ALIASES
        self._add_aliases()
        ### Add TRANSCRIPTS
        self._add_transcripts()
        ### Write as file
        self._write()
        ### Add to report
        self._to_report()

        if self.args.debug:
            print(DEBUG)
            print("url:", self.url)
            print("Gene Info")
            print(f" nb genes: {len(self.data['gene'])}")
            print(f" nb symbols: {len(self.data['symbol'])}")
            print(f" nb aliases: {len(self.data['alias'])}")
            print(f" nb transcripts: {len(self.data['transcript'])}")
            print(ENDCOL)


    def _set_refs(self):
        if __name__ == "__main__": print("Get gene.txt, please wait...")
        '''
        get gene.txt (gene_id, seq_region_id, xref_id, canonical_transcript_id, stable_id)
          key: line[0] = gene_id
          tupple:
            - $3  = seq_region_id
            - $4  = start
            - $5  = end
            - $6  = strand
            - $7  = xref_id
            - $11 = canonical_transcript_id
            - $12 = stable_id (ENSG)
        '''

        '''
        GENES
        '''
        gene_url = f"{self.url}gene.txt.gz"
        gene_r = requests.get(gene_url)
        gene_str = gzip.decompress(gene_r.content).decode()
        gene_id = {}
        gene_name = {}
        xrefs = {}
        # ~ print("seq_region_ids:", self.seq_region_ids)
        for line in gene_str.splitlines():
            line = line.split('\t')
            # ~ print("line:", line)
            # ~ if line[3] in self.seq_region_ids:   # item must be in a regular chromosom
            if line[3] in self.chr_dict:   # item must be in a regular chromosom
                gene_id[line[0]] = (line[3], line[7], line[11], line[12])
                # ~ print('line3:', line[3])
                # ~ print(f"self.seq_region_ids is a {type(self.seq_region_ids)}")
                gene_name[line[12]] = {'chr': self.chr_dict[line[3]],
                                       'start': int(line[4]),
                                       'end': int(line[5]),
                                       'strand': int(line[6]),
                                       'biotype': line[1],
                                       'desc': line[9],
                                       }
                xrefs.setdefault(line[7], []).append(line[0])

        del gene_r, gene_str

        self.refs["gene"] = gene_id
        self.refs["xref"] = xrefs
        self.data["gene"] = gene_name


    def _add_symbols(self):
        '''
        XREF - list of symbol-name
        kept  if xref_id is found in refs[xrefs] dict
          - add symbol names in data[gene] dict
          - populate data[symbol] dict (symbol-name: [ENSG01, ENSG02])
        '''
        if __name__ == "__main__": print("Get xrefs.txt, please wait...")

        xref_url = f"{self.url}xref.txt.gz"
        xref_r = requests.get(xref_url)
        xref_str = gzip.decompress(xref_r.content).decode()
        for line in xref_str.splitlines():
            line = line.split()
            if line[0] in self.refs["xref"]:
                for gene_id in self.refs["xref"][line[0]]:
                    ensg = self.refs["gene"][gene_id][3]
                    xref_id = line[0]
                    symbol = line[3]
                    self.data["gene"][ensg]["symbol"] = symbol                        # add to gene info
                    self.data["symbol"].setdefault(symbol.upper(), []).append(ensg)   # add to symbol index
        del xref_r, xref_str


    def _add_aliases(self):
        '''
        SYNONYM - list of aliases (from external sources)
        kept if xref_id is found on refs[xref] dict
          - add alias in data[gene] dict
          - populate data[alias] dict (alias: [symbol1, symbol2])
        '''
        if __name__ == "__main__": print("Get external_synonym.txt, please wait...")

        ### EXT_SYNONYMS - get external_synonym.txt (Aliases)
        ext_synonym_url = f"{self.url}external_synonym.txt.gz"
        ext_synonym_r = requests.get(ext_synonym_url)
        ext_synonyms_str = gzip.decompress(ext_synonym_r.content).decode()
        aliases = {}
        for line in ext_synonyms_str.splitlines():
            line = line.split('\t')
            if line[0] in self.refs["xref"]:
                alias = line[1]
                for gene_id in self.refs["xref"][line[0]]:
                    ensg = self.refs["gene"][gene_id][3]
                    symbol = self.data["gene"][ensg]["symbol"]
                    self.data["gene"][ensg].setdefault('aliases', []).append(alias)           # add to gene info
                    aliases.setdefault(alias.upper(), set()).add(ensg)                        # add to alias index
        self.data["alias"] = { k:list(v) for k,v in aliases.items()}


    def _add_transcripts(self):
        '''
        TRANSCRIPT - list of transcripts
        $0 = transcript_id
        $1 = gene_id
        $3 = region_id
        $13 = stable_id (ENST)
        '''
        if __name__ == "__main__": print("Get transcript.txt, please wait...")

        transcript_url = f"{self.url}transcript.txt.gz"
        transcript_r = requests.get(transcript_url)
        transcript_str = gzip.decompress(transcript_r.content).decode()
        for line in transcript_str.splitlines():
            line = line.split('\t')
            if line[1] in self.refs["gene"]:
                ensg_id = line[1]
                ensg = self.refs["gene"][ensg_id][3]
                enst = line[13]
                self.data["gene"][ensg].setdefault("transcript", []).append(enst)  # add to gene info
                self.data["transcript"][enst.upper()] = ensg                       # add to transcript index
                if line[0] == self.refs["gene"][ensg_id][2]:
                    self.data["gene"][ensg]["canonical"] = enst                    # add to gene info if canonical


    def _write(self):
        if "\\N" in self.data["symbol"]:
            del self.data["symbol"]["\\N"]

        # ~ with open(self.outfile, 'w', encoding='utf8') as outfile:
            # ~ yaml.dump(self.data, outfile, default_flow_style=False, allow_unicode=True)
        with open(self.outfile, 'wb') as fh:
            pickle.dump(self.data, fh)
            if self.args.debug: print(f"{DEBUG}File {os.path.basename(self.outfile)!r} created.{ENDCOL}")
        return


    def _to_report(self):
        """ Function doc """
        self.report.append("\n## Gene Info\n")
        self.report.append(f"- nb genes: {len(self.data['gene'])}")
        self.report.append(f"- nb symbols: {len(self.data['symbol'])}")
        self.report.append(f"- nb aliases: {len(self.data['alias'])}")
        self.report.append(f"- nb transcripts: {len(self.data['transcript'])}")
