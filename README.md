# Kmerator

## Prototype for decomposition of transcript or gene sequences and extraction of their specific k-mers

Kmerator is a prototype tool designed for the prediction of specific k-mers (also called tags) from input sequences, considering a reference genome and an ENSEMBL-like transcriptome. From these specific k-mers, it also outputs their corresponding specific contigs which are sequences of consecutive k-mers (overlapping length between k-mers must be k-1, otherwise, it's a new contig). Kmerator first uses Jellyfish [1] to create 2 requestable indexes from the reference genome and transcriptome, and second, decomposes your input transcript or gene sequences to count the occurences of each k-mer in the genome and transcriptome. Number of occurrences are then interpreted, in different manners, to select specific k-mer from your input. 

Kmerator strictly depends on a reference genome (fasta or jellyfish index format) and on an Ensembl fasta format transcriptome, you can find it there: https://www.ensembl.org/info/data/ftp/index.html. For a more complete k-mer filtering, we advice to merge the coding (cDNA) and non-coding (ncRNA) as one unique reference transcript.

Interpretation of occurrences counts differs depending on the ''level'' option : 

- gene 
- transcript 
- chimeras (must not be found in your reference genome and transcriptome)

If you want to use Kmerator on unannotated RNA signatures, you have to use the -u option. In this case, only "transcript" and "chimeras" levels are available.

## Dependencies

- Julia >= v1.4.2 (https://julialang.org/downloads/)
- Kmerator will automatically ask you to install missing required packages: ParallelDataTransfer, ArgParse, FastaIO, JSON and HTTP.
- Jellyfish >= 2.0


## Usage
```
kmerator.jl [--selection [SELECTION...]] [-f SEQUENCES] [--appris SPECIE] 
			[-u] [-s] [-v] -g GENOME -t TRANSCRIPTOME -l LEVEL [-o OUTPUT]
            [--length LENGTH] -p [NB-THREADS] [--threshold THRESHOLD]
```

## arguments
```
  --selection [SELECTION...]
 						list of gene IDs or transcript IDs (ENSG, ENST or gene 
 						Symbol) to select inside your fasta transcriptome file
 						and that you want to extract specific kmers from. If 
 						you want to use your own sequences, you can give your 
 						fasta file with --fasta_file option.                       
  -f, --fasta_file SEQUENCES
  						sequences in fasta format (genes or transcripts) that 
  						you want to extract specific kmers from. If you don't 
  						have your own sequences file, you can use a list of 
  						gene IDs or transcript IDs with --selection option.
  --appris SPECIE        
  						indicate: 'homo_sapiens', 'mus_musculus', 
  						'rattus_norvegicus','danio_rerio', 'sus_scrofa', or 
  						virtually any specie available in APPRIS database [2]. 
  						Genes have multiple possible transcripts called 
  						isoforms. This option selects the principal transcript 
  						defined in APPRIS database. If this option is used and 
  						there is no data available, no connection or no 
  						principal transcript (non-coding genes), the longer 
  						transcript is selected. 						
  -u, --unannotated     	
  						use this option if your provided sequences fasta file 
  						corresponds to annotations external from Ensembl. 
  						Otherwise, use ensembl annotations. 	
  -s, --stringent       
  						FOR GENE LEVEL ONLY: use this option if you want to 
  						select gene-specific k-mers present in ALL known 
  						transcripts for your gene. If false, a k-mer is 
  						considered as gene-specific if present in at least one 
  						isoform of your gene.               
  -v, --verbose         
  						if you want some updates while Kmerator is running.						
  -g, --genome GENOME   
  						genome fasta file or jellyfish index (.jf) to use for 
  						k-mers request.
  -t, --transcriptome TRANSCRIPTOME
                        	transcriptome fasta file (ENSEMBL fasta format ONLY) 
                        	to use for k-mers request and transcriptional variants 
                        	informations.
  -l, --level LEVEL     
  						use 'gene', 'transcript' or 'chimera' to extract 
  						specific k-mers at these different levels. 'chimera' 
  						option MUST be done with 'unannotated' (-u) option!
  -o, --output OUTPUT   
  						your output directory (default is current directory).
  --length LENGTH       
  						k-mer length that you want to use (type Int64, default 31).
  --threshold THRESHOLD
                        FOR GENE LEVEL ONLY: minimum fraction of annotated 
                        transcripts, from a given gene, containing this k-mer to 
                        keep it (default 0).
  -p, --procs NB-THREADS     
  						run n local processes.
  -h, --help            
  						show help message and exit.
```

## References

[1] Guillaume Marçais, Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers, Bioinformatics, Volume 27, Issue 6, 15 March 2011, Pages 764–770, https://doi.org/10.1093/bioinformatics/btr011
[2] Rodriguez JM, et al. Nucleic Acids Res. Database issue; 2017 Oct 23
