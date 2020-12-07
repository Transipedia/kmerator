#!/usr/bin/env julia

# AUTHOR : Sébastien RIQUIER, IRMB, Montpellier
# MODIF : Chloe BESSIERE

println("                                                       
 ____  __.                        | |            
|    |/ _| ___ _ __ _ __ ___  __ _| |_ ___  _ __   
|      <  | '_ ` _ \\/ _ \\ '__/ _` | __/ _ \\| '__|  
|    |  \\ | | | | ||  __/ |   (_| | || (_) | |     
|____|__ \\|_| |_| |_\\___|_|  \\__,_|\\__\\___/|_|     
        \\/                                              
-------------------------------------------------- 
      VERSION v0.2.1

Dependencies :
- Julia language >= v1.4.2
- Jellyfish >= v2.0
--------------------------------------------------")

## Pkg = Julia's builtin package manager (installing, updating and removing packages)
using Pkg

packages = ("ArgParse", "ParallelDataTransfer", 
            "HTTP", "JSON", "FastaIO")

## Evaluate and install package if missing
function check_installed_pkg()
    missing_pkg = []
    for mod in packages
        if ! haskey(Pkg.installed(), mod) 
        push!(missing_pkg, mod)
        end
    end
    if length(missing_pkg) != 0
        println("List of missing packages: ")
        for pkg in missing_pkg
            println("  $pkg")
        end
        print("Install automatically ? (y,n*): ")
        ack = chomp(readline(stdin))
        ## ask for automatically install missing packages
        if lowercase(ack) == "y"
            for pkg in missing_pkg
                println("Installing $pkg")
                Pkg.add(pkg)
            end
        else
            println("\nProgram aborted by user\n")
            exit()
        end
    end
end

check_installed_pkg()

using Distributed
using HTTP
using JSON
using ArgParse


## Parse arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--selection"
            help = """list of gene IDs or transcript IDs (ENSG, ENST or gene Symbol) to select 
                    inside your fasta transcriptome file and that you want to extract specific 
                    kmers from. If you want to use your own sequences, you can give your fasta 
                    file with --fasta_file option."""
            action = :append_arg 
            nargs = '*'
        "--fasta_file", "-f" 
            help = """sequences in fasta format (genes or transcripts) that you want to extract 
                    specific kmers from. If you don't have your own sequences file, you can use 
                    a list of gene IDs or transcript IDs with --selection option."""
        "--appris"
            help = """indicate: 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus',
                    'danio_rerio', 'sus_scrofa', or virtually any specie available in APPRIS 
                    database [Rodriguez JM, et al. Nucleic Acids Res. Database issue; 2017 Oct
                    23]. Genes have multiple possible transcripts called isoforms. 
                    This option selects the principal transcript defined in APPRIS database. 
                    If this option is used and there is no data available, no connection or no 
                    principal transcript (non-coding genes), the longer transcript is selected."""
            arg_type = String
        "--unannotated", "-u"
            action = :store_true
            help = """use this option if your provided sequences fasta file corresponds to 
                    annotations external from Ensembl. Otherwise, use ensembl annotations."""
        "--stringent", "-s"
            action = :store_true
            help = """FOR GENE LEVEL ONLY: use this option if you want to select gene-specific 
                    k-mers present in ALL known transcripts for your gene. If false, a k-mer is 
                    considered as gene-specific if present in at least one isoform of your gene."""
        "--verbose", "-v"
            action = :store_true
            help = """if you want some updates while Kmerator is running."""
        "--genome", "-g"
            help = "genome fasta file or jellyfish index (.jf) to use for k-mers request."
            required = true
        "--transcriptome", "-t"
            help = """transcriptome fasta file (ENSEMBL fasta format ONLY) to use for k-mers 
                    request and transcriptional variants informations."""
            required = true
        "--level", "-l"
            help = """use 'gene', 'transcript' or 'chimera' to extract specific 
                kmers at these different levels. 'chimera' option MUST BE DONE WITH 
                'unannotated' option!"""
            required = true
        "--output", "-o"
            help = """output directory (default is current directory)."""
            default = "output"
        "--length"
            help = """k-mer length that you want to use (default 31)."""
            default = 31
            arg_type = Int
        "--procs", "-p"
            help = """run n local processes"""
            arg_type = Int
        "--threshold"
            help = """FOR GENE LEVEL ONLY: minimum fraction of annotated transcripts, 
                for a given gene, containing this kmer to keep it (default 0)."""
            default = 0.0
            arg_type = Float64
#        "--contig", "-c"
#            action = :store_true
#            help = """if you want to keep ."""
    end
    return parse_args(s)
end


## Set variables with user parameters
function set_variables()
    ## fasta file
    global fastafile = parsed_args["fasta_file"] 
    ## Verbose
    global verbose_option = parsed_args["verbose"]
    ## Output
    global output = parsed_args["output"]
    if occursin(r"\/$", output)
        global output = replace(output, r"/$" => "")
    end
    ## selection
    global select_option = parsed_args["selection"]
    if isempty(select_option)
	select_option = nothing
    else
        select_option = select_option[1]
    end
    ## APPRIS
    global APPRIS_option = parsed_args["appris"]
    if APPRIS_option != nothing
        if verbose_option println("APPRIS selection of principal transcripts for $APPRIS_option" ) end
    end
    ## Unannotated option
    global unannotated_option = parsed_args["unannotated"]
    ## Stringent option
    global stringent_option = parsed_args["stringent"]
    ## Admission threshold option
    global admission_threshold = parsed_args["threshold"]
    ## Genome and transcriptome options
    global genome = parsed_args["genome"]
    global transcriptome = parsed_args["transcriptome"]
    ## Kmer length option
    global kmer_length = parsed_args["length"]
    ## Threads number
    if isa(parsed_args["procs"], Number)
        global nbthreads = parsed_args["procs"]
    else
        nbthreads = 1
    end
    ## Level of kmer specificity
    global level = parsed_args["level"]
end

function add_workers()
    if isa(parsed_args["procs"], Number) 
            addprocs(parsed_args["procs"])
    end
end

## Check up the conformity of the variables
function checkup_variables()
    ## test genome and transcriptome files
    if ! isfile(genome)
        println("\nFileError: genome file $genome not found\n")
        exit(2)
    elseif ! occursin(r".*\.(fa|fasta|fna|jf|jl)$", genome)
        println("\nFileTypeError: $genome is not a fasta or jellyfish file\n")
        exit(2)
    end
    if ! isfile(transcriptome)
        println("\nFileError: transcriptome file $transcriptome not found\n")
        exit(2)
    elseif ! occursin(r".*\.(fa|fasta)$", transcriptome)
        println("\nFileTypeError: $transcriptome is not a fasta file\n")
        exit(2)
    else 
        FastaReader(transcriptome) do fr
            for (desc, seq) in fr
                desc_array = split(desc)
                if length(desc_array)<7 || !occursin("gene_symbol:", desc_array[7])
                    println("\nFileError: $transcriptome not in ENSEMBL fasta format, use ENSEMBL transcriptome fasta file\n")
                    exit(2)
                elseif (!occursin("ENST", desc_array[1]))
                    println("\nFileError: $transcriptome not in ENSEMBL fasta format, use ENSEMBL transcriptome fasta file\n")
                    exit(2)
                elseif (!occursin("gene:ENSG", desc_array[4]))                   
                    println("\nFileError: $transcriptome not in ENSEMBL fasta format, use ENSEMBL transcriptome fasta file\n")
                    exit(2)
                end
                break
            end
        end
    end
    ## test fasta file and selection options   
    if select_option != nothing
        if fastafile != nothing
            println("\nInputError: you can't provide both gene list and sequences fasta file\n")
            exit(2)
	else
	    global fastafile = parsed_args["transcriptome"]
        end
    else
        if fastafile == nothing
            println("\nInputError: you must provide a gene list OR a sequences fasta file\n")
            exit(2)        
        elseif !isfile(fastafile)
            println("\nFileError: sequences fasta file not found ($fastafile)\n")
            exit(2)
        elseif !occursin(r".*\.(fa|fasta|fna)$", fastafile)
            println("\nFileTypeError: $fastafile is not a fasta file\n")
            exit(2)           
        end
    end
    if select_option != nothing
        for gene in select_option
            # ENSEMBL annotations ~ ENSTXXX.X, ENSGXXX.X (version) => forbidden into select option (remove .X)
            if occursin(".", gene) && occursin("ENS", gene)
                println("""\n ENSEMBL annotations with version (point) like ENSTXXXXXXXX.1 
                    is forbidden, just remove the '.1'.\n""")
            end
        end
    end
    ## APPRIS option works only with the gene annotated level
    if APPRIS_option != nothing && (level != "gene" || unannotated_option == true)
        println("\n ApprisError: APPRIS option works only with the gene annotated level\n")
        exit(2)
    end
    ## chimera level works only with unnanotated option
    if !unannotated_option && level == "chimera" 
        println("\n DependenceError: you must use chimera level with unannotated option\n")
        exit(2)
    end
    ## Show options if verbose option
    if verbose_option 
        println("
        \r OPTIONS
        \r ------------
        \r Output directory:           $output
        \r Genes selection option:     ", (select_option != nothing) ? (select_option, ", ") : "nothing", "
        \r APPRIS option:              ", repr(APPRIS_option),"
        \r Verbose option:             ", verbose_option ? "yes" : "no", "
        \r Unannotated option:         ", unannotated_option ? "yes" : "no", "
        \r Stringent option:           ", stringent_option ? "yes" : "no", "
        \r Admission threshold option: $admission_threshold
        \r Genome file:                $genome
        \r Transcriptome file:         $transcriptome
        \r Kmer length:                $kmer_length
        \r Input sequences file:      ", (fastafile != nothing) ? (fastafile) : "no sequences", "
        \r Level of kmer specificity:  $level
        \r Nb process:                 $nbthreads")
    end
end

function load_transcriptome(transcriptome_file)
    ! verbose_option ? print("\r") : println("\r ------------")
    print("\nCreating a dictionary of the transcriptome fasta... \n")
    transcriptome_dict = Dict{String,String}()
    FastaReader(transcriptome_file) do fr
        if unannotated_option == true
            for (desc, seq) in fr
                transcriptome_dict[desc] = seq
            end
        else
            for (desc, seq) in fr
                desc_array = split(desc)
                gene_name = split(desc_array[7], ":")[2]
                ensembl_transcript_name = split(desc_array[1],'.')[1]
                ensembl_gene_name = split(split(desc_array[4], ":")[2], '.')[1]
                transcriptome_dict["$gene_name:$ensembl_transcript_name"] = seq
            end
        end
    end
    if verbose_option
        println("\nFirst sequence of the transcriptome: ")
        for (name, seq) in transcriptome_dict
            println(name, "\n", seq)
            break
        end
        println("\nTranscripts dictionary done")
    end
    return transcriptome_dict
end


function run_jellyfish(genome, transcriptome_fa)
    jf_dir = "$output/jellyfish_indexes/$kmer_length"
    ## run jellyfish on transcriptome
    ! verbose_option ? print("\r") : println("\r ------------ \n")
    print("Running jellyfish on the transcriptome... ")
    jf_transcriptome = replace(basename(transcriptome_fa), ".fa" => ".jf")
    run(`mkdir -p $jf_dir`)
    global transcriptome = "$jf_dir/$(replace(basename(transcriptome_fa), ".fa" => ".jf"))"
    run(`jellyfish count -m $kmer_length -s 10000 -t $nbthreads -o $jf_dir/$jf_transcriptome $transcriptome_fa`)
    if verbose_option println("\nTranscriptome kmer index output: $jf_dir/$jf_transcriptome") end
    
    ## run jellyfish on genome if genome is fasta file
    if occursin(r".*\.(fa|fasta)$", genome)
        ! verbose_option ? print("\r") : println("\r ------------")
        println("Compute jellyfish on Genome...  ")
        jf_genome = replace(basename(genome), ".fa" => ".jf")
        run(`jellyfish count -C -m $kmer_length -s 10000 -t $nbthreads -o $jf_dir/$jf_genome $genome`)
        if verbose_option println("\nGenome kmer index output: $jf_dir/$jf_genome") end
    else
        println("Jellyfish genome index already provided")
        jf_genome = basename(genome)
        jf_dir = dirname(genome)
    end
    return(jf_genome, jf_dir)
end


function APPRIS_function(gene_ref)
    ! verbose_option ? print("\r") : println("\r ------------")
    print("Finding the principal isoform based on APPRIS database for $gene_ref... \n")
    ## Request to appris
    url = "http://apprisws.bioinfo.cnio.es/rest/exporter/id/$APPRIS_option/" * 
            "$gene_ref?methods=appris&format=json&sc=ensembl"
    if verbose_option println("\nAPPRIS url: $url") end
    ## make_API_call(url)

    function http_req(url)
        try
            r = HTTP.request("GET", url; verbose=0)
            res = String(r.body)
            res = JSON.Parser.parse(res)
            return res
        catch err
            if verbose_option println("ERROR:\n $err") end
            res = "NODATA"
            return res
        end
    end

    transcripts = http_req(url)
    if verbose_option && isa(transcripts, Array)
        for transcript in transcripts
            println("First transcript found:")
            for (key,value) in transcript
                println("  $key: $value")
            end
            break
        end
    end
    ## Finding the best isoform, find Principals
    principals = []
    for (i,value) in enumerate(transcripts)
        try 
            if occursin("PRINCIPAL:", transcripts[i]["reliability"])
                push!(principals, transcripts[i]["reliability"]) 
            end
        catch
        end
    end
    ## if no principals, return NODATA
    if isempty(principals)
        if verbose_option
            println("No principal isoform detected, this function will return " * 
                        "'NODATA' and the longest transcript will be selected")
        end
        return(transcripts)
    end
    if verbose_option println("principals: $principals \n") end
    ## Select best Principal (minimum level)
    levels = []
    for i in 1:length(principals)
        push!(levels, split(principals[i],":")[2]) 
    end
    levels = map(x-> parse(Int, x), levels)
    level = minimum(levels)
    ## if multiple transcripts with better Principal, select the biggest
    selected_transcripts = hcat(map(x -> try x["transcript_id"] 
                                         catch end, transcripts),
                                map(x -> try parse(Int, x["length_na"]) 
                                         catch end, transcripts),  
                                map(x -> try x["reliability"] 
                                         catch end, transcripts))
    selected_transcripts = selected_transcripts[map(x -> x != nothing, selected_transcripts[:, 3]), :]
    selected_transcripts = selected_transcripts[map(x -> occursin("PRINCIPAL:$level",x), selected_transcripts[:, 3]), :]
    max_length = maximum(selected_transcripts[:,2])
    selected_transcripts = selected_transcripts[map(x -> x == max_length, selected_transcripts[:, 2]), :]
    best_transcript = String( unique(selected_transcripts[:,1])[1])
    if verbose_option println("APPRIS result : $best_transcript \nAPPRIS function finished \n\r ------------ \n\n") end
    return(best_transcript)
end # end of APPRIS function


function find_longest_variant(gene_name, transcriptome_dict)
    if verbose_option println("\r ------------ \n\nFinding the longest variant for the gene $gene_name") end
    gene_name = replace(gene_name, "@SLASH@" => "/")
    variants_dict = filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict)
    nb_variants = length(variants_dict)
    if verbose_option println("\nnb variants : $nb_variants") end
    variants_lengths = Vector{Any}()
    for (k,v) in variants_dict
        push!(variants_lengths, length(v))
    end
    longest_variant = String(collect(keys(filter(v -> length(v.second) == maximum(variants_lengths), variants_dict)))[1])
    ## sometimes, gene name contains a slash! 
    gene_name = replace(gene_name, "/" => "@SLASH@")
    if verbose_option println("$gene_name : longest variant = $longest_variant \nFinding the longest variant done\n \r ------------ \n") end
    return(longest_variant)
end


function build_sequences() ## Creating individual sequence files from input fasta file or genes/transcripts list
    if select_option == nothing
	println("\r ------------ \n\nSplitting sequences of your fasta input file... \n")
    else   
	println("\r ------------ \n\nCreating sequences from your transcripts/genes list... \n")     
    end 
    run(`mkdir -p "$output/sequences"`)
    if !unannotated_option	## ensembl annotations only

        genes_already_processed = []
        genes_analysed = []
        FastaReader(fastafile) do fr
            for (desc, seq) in fr
                desc_array = split(desc)
                APPRIS_transcript = ""
                gene_name = replace(desc_array[7], "gene_symbol:" => "")
                ensembl_transcript_name = split(desc_array[1],'.')[1]
                ensembl_gene_name = split(replace(desc_array[4], "gene:" => ""), '.')[1]
                ## Transcript level
                if level == "transcript" 
                    if !isempty(select_option) && !(ensembl_transcript_name in select_option) 
                       continue
                    else 
                        if length("$seq") >= kmer_length 
                            if verbose_option println("$ensembl_transcript_name: sequence length >= $kmer_length => continue") end
                            gene_name = replace(gene_name, "/" => "@SLASH@") # some gene names can contain slash characters that break the processus
                            FastaWriter("$output/sequences/$gene_name:$ensembl_transcript_name.fa") do fwsequence
                                write(fwsequence, [">$gene_name:$ensembl_transcript_name", "$seq"])
                            end
                        else println("$ensembl_transcript_name: sequence length < $kmer_length => ignored") end
                    end 
                end  
                ## Gene level
                if level == "gene"
                    ## testing if gene has not been already processed
                    if isempty(select_option) || (gene_name in select_option || ensembl_gene_name in select_option) && !(gene_name in genes_already_processed) ## keeping all sequences corresponding to genes into select_option or take them all if select_option == false
                        if APPRIS_option != nothing && !(gene_name in genes_analysed)
                            APPRIS_transcript = APPRIS_function(ensembl_gene_name) 
                            push!(genes_analysed, gene_name)
			    # println("APPRIS selected variant : $APPRIS_transcript")
                            # écrire dans un fichier ?
                        end
                        if APPRIS_option == nothing || (ensembl_transcript_name == APPRIS_transcript) || (APPRIS_transcript == "NODATA" && "$gene_name:$ensembl_transcript_name" == find_longest_variant(gene_name,transcriptome_dict)) 
                            if verbose_option println("ensembl_transcript_name : $ensembl_transcript_name") end
                            if length("$seq") >= kmer_length 
                                if verbose_option println("$gene_name:$ensembl_transcript_name: sequence length >= $kmer_length => continue") end
                                gene_name = replace(gene_name, "/" => "@SLASH@")
                                FastaWriter("$output/sequences/$gene_name:$ensembl_transcript_name.fa") do fwsequence
                                    write(fwsequence, [">$gene_name:$ensembl_transcript_name", "$seq"])
                                    push!(genes_already_processed, gene_name)
                                end
                            else 
                                println("$gene_name: sequence length < $kmer_length => ignored")
                            end
                        end
                    end
                end
            end
            for i in unique(genes_analysed)
                ## it can happen that APPRIS contains obsolete IDs/error => longer variant selected for the concerned gene 
                if (replace(i, "/" => "@SLASH@") in genes_already_processed) == false
                    if verbose_option
                        println("""WarningAPPRIS : There is a problem with $i, 
                            \rAPPRIS function has not worked properly, sometimes 
                            \ran outdated APPRIS reference, in this case the longer 
                            \rtranscript is selected from the reference fasta""")
                    end
                    id = find_longest_variant(i, transcriptome_dict)
                    FastaWriter("$output/sequences/$(replace(id, "/" => "@SLASH@"))") do fwsequence
                        write(fwsequence, [">$id", string(transcriptome_dict["$id"])])
                    end
                end
            end
        end

    else  ## unannotated option
        FastaReader(fastafile) do fr
            for (desc, seq) in fr
                gene_name = desc
                if length("$seq") >= kmer_length
                    FastaWriter("$output/sequences/$desc.fa") do fwsequence
			if length("$desc")>79			
				a = desc[1:79]
				write(fwsequence, [">$a", "$seq"])
			else
				write(fwsequence, [">$desc", "$seq"])
			end
                    end
                else println("$desc: sequence length < $kmer_length => ignored") end
            end
        end
    end
    println("OK")
    if verbose_option
        println("Sequences:")
        for seq in (readdir("$output/sequences/"))
            println(" $seq")
        end
    end
    if verbose_option println("\nSequences splitting into individual files done \n \r ------------ \n") end
end


## MAIN
parsed_args = parse_commandline()				# parse arguments
set_variables()							# set variables with input arguments
add_workers()
@everywhere using Dates
@everywhere using DelimitedFiles
@everywhere using FastaIO
using  ParallelDataTransfer
checkup_variables()						# check variables
global transcriptome_dict = load_transcriptome(transcriptome) 		# transcriptome dictionary
jf_genome, jf_dir = run_jellyfish(genome, transcriptome) 	# genome.jf and dir_genomejf/
build_sequences() 						# create 1 file by requested sequence

@passobj 1 workers() unannotated_option
@passobj 1 workers() verbose_option
@passobj 1 workers() transcriptome_dict
@passobj 1 workers() genome
@passobj 1 workers() jf_genome
@passobj 1 workers() jf_dir
@passobj 1 workers() transcriptome
@passobj 1 workers() level
@passobj 1 workers() output
@passobj 1 workers() kmer_length
@passobj 1 workers() stringent_option
@passobj 1 workers() admission_threshold


## Extraction of specific kmers for one splitted fasta file (one sequence)
@everywhere function f(splitted_fasta_files)
    fasta_array = Array([])             	# specific kmers table
    fasta_contig_array = Array([])      	# specific contigs table
    myid() == 1 ? X = myid() : X = myid()-1 	# myid() get the id of the current process
    elapsed_time = 0
    step_time = 0
    kmers_analysed = 0
    replace_line = function(nbline::Int64, x::String) # function to replace the current line on the screen
        print(string("\u1b[$(nbline)F\u1b[2K"*"$x"*"\u1b[$(nbline)E"))
    end

    if ! unannotated_option             	# if annotated
	name = first(split("$splitted_fasta_files",".fa"))
        gene_name, transcript_name = split(name, ":")
        variants_dict = filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict)
        nb_variants = length(variants_dict)
        if verbose_option println("gene_name: $gene_name --- transcript_name: $transcript_name --- nb variants : $nb_variants") end
    else                                	# if unannotated
        gene_name = transcript_name = first(split("$splitted_fasta_files","."))
        if verbose_option println("gene_name = transcript_name: $transcript_name") end
    end
    ## Create tag-files names
    if level == "gene"
        tag_file = "$gene_name-$transcript_name-gene_specific.fa"
        contig_file = "$gene_name-$transcript_name-gene_specific_contigs.fa"
    elseif level == "transcript"
	    if unannotated_option
        	tag_file = "$gene_name-transcript_specific.fa"
        	contig_file = "$gene_name-transcript_specific_contigs.fa"
	    else
        	tag_file = "$gene_name-$transcript_name-transcript_specific.fa"
        	contig_file = "$gene_name-$transcript_name-transcript_specific_contigs.fa"
	    end
    elseif level == "chimera"
        tag_file = "$gene_name-chimera_specific.fa"
        contig_file = "$gene_name-chimera_specific_contigs.fa"
    end
    ## take the transcript sequence for jellyfish query
    sequence_fasta = "error if seen"
    FastaReader("$output/sequences/$splitted_fasta_files") do fr
        for (desc, seq) in fr
            sequence_fasta = seq
        end
    end
    if verbose_option remotecall(replace_line, 1 , X, "worker number $(myid()), for $level $gene_name, reporting for duty ! :)\n") end

    ## building kmercounts dictionary from jellyfish query on the genome
    println("Jellyfish query -s $output/sequences/$splitted_fasta_files $jf_dir/$jf_genome")
    kmercounts_genome = read(`jellyfish query -s "$output/sequences/$splitted_fasta_files" "$jf_dir/$jf_genome"`, String)
    if verbose_option remotecall(replace_line, 1 , X, "jellyfish query $splitted_fasta_files on $jf_dir/$jf_genome finished") end
    kmercounts_genome = split(kmercounts_genome, "\n")[1:end-1]
    kmercounts_genome_dict = Dict()
    for mer in kmercounts_genome
        mer = split(mer)
        seq = mer[1]
        kmercounts_genome_dict["$seq"] = mer[2]
    end
    ## building kmercounts dictionary from jellyfish query on the transcriptome
    kmercounts_transcriptome = read(`jellyfish query -s "$output/sequences/$splitted_fasta_files" "$transcriptome"`, String)
    if verbose_option remotecall(replace_line, 1 , X, "Jellyfish query on $transcriptome done") end
    kmercounts_transcriptome_dict = Dict()
    kmercounts_transcriptome = split(kmercounts_transcriptome, "\n")[1:end-1]
    for mer in kmercounts_transcriptome
        mer = split(mer)
        seq = mer[1]
        kmercounts_transcriptome_dict["$seq"] = mer[2]
    end

    ## initialization of count variables
    i = 0
    j = 1
    total_kmers = length(keys(kmercounts_transcriptome_dict))
    println("Total kmers = $total_kmers")

    ## creating a new dictionary with kmers and their first position in our query sequence
    kmer_starts = Dict()
    kmer_placed = 0
    for mer in keys(kmercounts_transcriptome_dict)
        kmer_placed = kmer_placed +1
        if verbose_option remotecall(replace_line, 1 , X, "kmers placed = $kmer_placed on $total_kmers \n") end
        kmer_start = minimum(collect(findfirst("$mer", "$sequence_fasta")))
        kmer_starts["$mer"] = kmer_start
    end
    kmer_starts_sorted = sort(collect(zip(values(kmer_starts),keys(kmer_starts)))) # array sorted by kmer position
    position_kmer_prev = first(kmer_starts_sorted[1])
    contig_string = ("") # initialize contig string

    for tuple in kmer_starts_sorted
        ## from the kmer/position sorted array, we extract sequence if specific (occurence ==1)
        mer = last(tuple) # kmer sequence
        position_kmer = first(tuple) # kmer position
        startt = time()
        kmers_analysed = kmers_analysed + 1
        per = round(kmers_analysed/total_kmers*100)
        elapsed_time = elapsed_time + step_time
        if elapsed_time == 0
            time_remaining = 0
        else
            time_remaining = Integer(ceil((total_kmers - kmers_analysed)/(kmers_analysed/elapsed_time)))
        end
        ptime = Dates.canonicalize(Dates.CompoundPeriod(Dates.Second(time_remaining)))
        ptime = time_remaining
        #remotecall(replace_line, 1 , X, "$level $gene_name : $kmers_analysed kmers analysed ( $per %), $i specifics,:remaining time -> $ptime sec.")
        
        if !any(x->x == "$mer", keys(kmercounts_genome_dict))
            d = Dict("A"=>"T", "C"=> "G", "G"=>"C", "T"=>"A")
            new_mer = reverse(replace(mer, r"[ACGT]{1}" => x->d[x]))  # "TGCA"
            genome_count = kmercounts_genome_dict["$new_mer"]
        else
            genome_count = kmercounts_genome_dict["$mer"]
        end
        transcriptome_count = parse(Float64, kmercounts_transcriptome_dict["$mer"])

        if level == "gene"
            ## if the kmer is present/unique or does not exist (splicing?) on the genome
            if parse(Int, genome_count) <= 1
                mer_regex = Regex(mer)
                variants_containing_this_kmer = keys(filter(v -> occursin(mer_regex, v.second), variants_dict))
                if stringent_option && Float64(transcriptome_count) == Float64(nb_variants) == Float64(length(variants_containing_this_kmer)) 
                    i = i+1
                    tmp = length(variants_containing_this_kmer)
                    push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)")
                    push!(fasta_array,"$mer")
                    # contigs 
                    if (position_kmer == Int64(position_kmer_prev+1)) 
                        contig_string = string(contig_string,mer[end])
                        position_kmer_prev = position_kmer
                    else
                        if contig_string != ""
                            push!(fasta_contig_array,">$gene_name-$transcript_name.contig$j")
                            push!(fasta_contig_array,"$contig_string")
                            j = j+1
                        end
                        contig_string = mer
                        position_kmer_prev = position_kmer
                    end                   
                elseif !stringent_option && Float64(transcriptome_count) == Float64(length(variants_containing_this_kmer)) && Float64(transcriptome_count) > Float64(nb_variants*admission_threshold)
                    i = i+1
                    tmp = length(variants_containing_this_kmer)
                    push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)" )
                    push!(fasta_array, "$mer")
                    # contigs 
                    if (position_kmer == Int64(position_kmer_prev+1)) 
                        contig_string = string(contig_string,mer[end])
                        position_kmer_prev = position_kmer
                    else
                        if contig_string != ""
                            push!(fasta_contig_array,">$gene_name-$transcript_name.contig$j")
                            push!(fasta_contig_array,"$contig_string")
                            j = j+1
                        end
                        contig_string = mer
                        position_kmer_prev = position_kmer
                    end                   
                end 
            end
        end

        if level == "transcript"
            if unannotated_option && Float64(transcriptome_count) == Float64(0) && parse(Int, genome_count) <= 1
                i = i+1
                push!(fasta_array,">$gene_name.kmer$i")
                push!(fasta_array,"$mer")
                # contigs 
                if (position_kmer == Int64(position_kmer_prev+1)) 
                    contig_string = string(contig_string,mer[end])
                    position_kmer_prev = position_kmer
                else
                    if contig_string != ""
                        push!(fasta_contig_array,">$gene_name.contig$j")
                        push!(fasta_contig_array,"$contig_string")
                        j = j+1
                    end
                    contig_string = mer
                    position_kmer_prev = position_kmer
                end

            elseif !unannotated_option && Float64(transcriptome_count) == Float64(1) && parse(Int, genome_count) <= 1
                i = i+1
                push!(fasta_array,">$gene_name-$transcript_name.kmer$i")
                push!(fasta_array,"$mer")
                # contigs 
                if (position_kmer == Int64(position_kmer_prev+1)) 
                    contig_string = string(contig_string,mer[end])
                    position_kmer_prev = position_kmer
                else
                    if contig_string != ""
                        push!(fasta_contig_array,">$gene_name-$transcript_name.contig$j")
                        push!(fasta_contig_array,"$contig_string")
                        j = j+1
                    end         
                    contig_string = mer
                    position_kmer_prev = position_kmer
                end
            end
        end

        if level == "chimera" && unannotated_option == true
            if Float64(transcriptome_count) == Float64(0) && parse(Int, genome_count) == 0
                i = i+1
                push!(fasta_array,">$gene_name.kmer$i")
                push!(fasta_array,"$mer")
                # contigs 
                if (position_kmer == Int64(position_kmer_prev+1)) 
                    contig_string = string(contig_string,mer[end])
                    position_kmer_prev = position_kmer
                else
                    if contig_string != ""
                        push!(fasta_contig_array,"$gene_name.contig$j")
                        push!(fasta_contig_array,"$contig_string")
                        j = j+1
                    end         
                    contig_string = mer
                    position_kmer_prev = position_kmer
                end
            end
        end
        step_time = Float64(time()-startt)
    end # end of kmer analysis
    
    # writting last contig
    if level == "gene" && contig_string != ""
        push!(fasta_contig_array,">$gene_name-$transcript_name.contig$j")
        push!(fasta_contig_array,"$contig_string")    
    elseif level == "transcript" && unannotated_option == false && contig_string != ""
        push!(fasta_contig_array,">$gene_name-$transcript_name.contig$j")
        push!(fasta_contig_array,"$contig_string")
    elseif (level == "chimera" || (level == "transcript" && unannotated_option == true)) && contig_string != ""
        push!(fasta_contig_array,">$gene_name.contig$j")
        push!(fasta_contig_array,"$contig_string")
    end

    ## write tag files
    run(`mkdir -p "$output/tags/$kmer_length"`)
    run(`mkdir -p "$output/contigs"`)
    writedlm("$output/tags/$kmer_length/$tag_file", reshape(fasta_array,length(fasta_array)), "\n")
    writedlm("$output/contigs/$contig_file", reshape(fasta_contig_array,length(fasta_contig_array)), "\n")
end #end of function

print("\rExtracting specific kmers                                  ")
for n in 1:nworkers()+1 println("") end
pmap(f, readdir("$output/sequences/"))

