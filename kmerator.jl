#!/usr/bin/env julia

# AUTHOR : Sébastien RIQUIER, IRMB, Montpellier

println("                                                       
 ____  __.                        | |            
|    |/ _| ___ _ __ _ __ ___  __ _| |_ ___  _ __   
|      <  | '_ ` _ \\/ _ \\ '__/ _` | __/ _ \\| '__|  
|    |  \\ | | | | ||  __/ |   (_| | || (_) | |     
|____|__ \\|_| |_| |_\\___|_|  \\__,_|\\__\\___/|_|     
        \\/                                              
-------------------------------------------------- 
      VERSION v0.2.0

Dependencies :
- Julia language >= v1.4.2
- Jellyfish >= v2.0
--------------------------------------------------")

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



function parse_commandline()
    ## parse argument
    s = ArgParseSettings()
    @add_arg_table s begin
        "--selection"
            help = "list of gene to select directly inside your fasta transcriptome file"
            action = :append_arg 
            nargs = '*'
        "--appris"
            help = """indicate: 'homo_sapiens', 'mus_musculus', 'rattus_norvegicus',
                'danio_rerio', 'sus_scrofa', or virtually any specie available on 
                APPRIS database (Rodriguez JM, Rodriguez-Rivas J, Domenico TD, 
                Vázquez J, Valencia A, and Tress ML. Nucleic Acids Res. Database 
                issue; 2017 Oct 23.). A gene have multiple possible sequence depending 
                of its variants. this option select the principal transcript of a gene, 
                based of APPRIS database, or by length if no data or no connection. 
                If 'APPRIS option and not available data or no principal transcript 
                (non-coding genes), use length instead"""
            arg_type = String
        "--unannotated", "-u"
            action = :store_true
            help = """activated if the provided initial fasta file correspond to an 
                annotation external from Ensembl. Otherwise, use ensembl fasta files !"""
        "--stringent", "-s"
            action = :store_true
            help = """if you think a gene-specific tag must be unique but present in 
                ALL KNOWN TRANSCRIPT, if false, the kmer is considered as gene-specific 
                is present only in variant of the corresponding gene, independently 
                of the number know variant containing the kmer"""
        "--verbose", "-v"
            action = :store_true
            help = "If you want this script talk too much"
        "--genome", "-g"
            help = "the genome fasta (.fa) or index by jellyfish for kmer request"
            #type = String
            required = true
        "--transcriptome", "-t"
            help = """the transcriptome fasta (.fa) (FASTA ONLY) for kmer request 
                and transcriptional variant informations"""
            #type = String
            required = true
        "--level", "-l"
            help = """type 'gene', 'transcript' or 'chimera' to extract specific 
                kmer at these different levels. 'chimera' option must be done with 
                'unannotated' option!."""
            required = true
        "--output", "-o"
            help = "directory of output"
            default = "output"
        "--length"
            help = "length required for the kmer generation"
            default = 31
            arg_type = Int
        "--procs", "-p"
            help = "Run n local processes"
            arg_type = Int
        "fasta_file" #argument, no option
            help = "the fasta input file"
            required = true
        "--threshold"
            help = """FOR GENE LEVEL ONLY : Minimum fraction of annotated transcripts 
                containing this kmer to admit it (default 0)"""
            default = 0.0
            arg_type = Float64
    end
    return parse_args(s)
end



function set_variables()
    ## Verbose
    global verbose_option = parsed_args["verbose"]
    ## Output
    global output = parsed_args["output"]
    if occursin(r"\/$", output)
        global output = replace(output, r"/$" => "")
    end
    ## selection
    global select_option = parsed_args["selection"]
    if !isempty(select_option)
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
    ## Genome option
    global genome = parsed_args["genome"]
    ## Transcriptome option
    global transcriptome = parsed_args["transcriptome"]
    ## kmer length option
    global kmer_length = parsed_args["length"]
    ## Fasta file option (input sequences file)
    global fastafile = parsed_args["fasta_file"]
    ## Threads number
    if isa(parsed_args["procs"], Number)
        global nbthreads = parsed_args["procs"]
    else
        nbthreads = 1
    end
    ## Level of kmer specificity
    global level = parsed_args["level"]
    ## Show options if verbose option
    if verbose_option println("
        \r OPTIONS
        \r -------
        \r Output directory:           $output
        \r Selection option, genes:    ", length(select_option)>0 ? (select_option, ", ") : "nothing", "
        \r APPRIS option:              ", repr(APPRIS_option),"
        \r Verbose option:             $verbose_option
        \r Unannotated option:         ", unannotated_option ? "yes" : "no", "
        \r Stringent option:           $stringent_option
        \r Admission threshold option: $admission_threshold
        \r Genome option:              $genome
        \r Transcriptome option:       $transcriptome
        \r Kmer length option:         $kmer_length
        \r Input sequences files:      $fastafile
        \r Level of kmer specificity:  $level
        \r nb process:                 $nbthreads")
    end
end


function add_workers()
    if isa(parsed_args["procs"], Number) 
            addprocs(parsed_args["procs"])
    end
end


function checkup_variables()
    ## test genome file
    if ! isfile(genome)
        println("\nFileError: genome not found (file: $genome)\n")
        exit(2)
    elseif ! occursin(r".*\.(fa|fasta|jf|jl)$", genome)
        println("\nFileTypeError: not a fasta or jellyfish file (file: $genome)\n")
        exit(2)
    end
    ## test transcriptome file
    if ! isfile(transcriptome)
        println("\nFileError: transcriptome not found (file: $transcriptome)\n")
        exit(2)
    elseif ! occursin(r".*\.(fa|fasta)$", transcriptome)
        println("\nFileTypeError: not a fasta file (file: $transcriptome)\n")
        exit(2)
    end
    ## test fasta file
    if ! isfile(fastafile)
        println("\nFileError: fasta file not found (file: $fastafile)\n")
        exit(2)
    elseif ! occursin(r".*\.(fa|fasta)$", fastafile)
        println("\nFileTypeError: not a fasta file (file: $fastafile)\n")
        exit(2)
    end
    ## test selection option
    if select_option != nothing
        for gene in select_option
            # ENSEMBL annotation are ENSTXXXX.XX or ENSGXXXXX.XX in ensembl annotation, 
            # to avoid bad recognition for selection option, 
            # putting annotation like ENSTXXXXX.X for selection is forbidden
            if occursin(".", gene) && occursin("ENS", gene)
                println("""\n ENSEMBL annotations with a point like ENSTXXXXXXXX.1 
                    is forbidden, just remove the '.1'.\n""")
            end
        end
    end
    ## APPRIS option works only with the gene annotated level
    if APPRIS_option != nothing && (level != "gene" || unannotated_option == true)
        println("\n ApprisEror: APPRIS option works only with the gene annotated level\n")
        exit(2)
    end
    ## verification of fasta ensembl annotation
    if ! occursin(r".*\.(fa|fasta)$", transcriptome)
        println("\n FileTypeError: $transcriptome seems not a fasta file\n") 
        exit(2)
    end
end


function load_transcriptome(transcriptome_file)
    ! verbose_option ? print("\r") : println("------------")
    print("\nCreate dictionary of transcriptome fasta ")
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
        println("First Sequence of transcriptome:")
        for (name, seq) in transcriptome_dict
            println(name, "\n", seq)
            break
        end
    end
    if verbose_option println("transcripts_dictionary finished") end
    return transcriptome_dict
end

function run_jellyfish(genome, transcriptome_fa)
    jf_dir = "$output/jellyfish_indexes/$kmer_length"
    
    ## run jellyfish on Transcriptome
    ! verbose_option ? print("\r") : println("------------")
    print("Compute jellyfish on Transcriptome       ")
    jf_transcriptome = replace(basename(transcriptome_fa), ".fa" => ".jf")
    run(`mkdir -p $jf_dir`)
    run(`jellyfish count -m $kmer_length -s 10000 -t $nbthreads -o $jf_dir/$jf_transcriptome $transcriptome_fa`)
    global transcriptome = "$output/jellyfish_indexes/$kmer_length/$(replace(basename(transcriptome_fa), ".fa" => ".jf"))"
    if verbose_option println("\nTranscriptome kmer index output: $jf_dir/$jf_transcriptome.") end
    
    ## run jellyfish on genome if genome is fasta file
    if occursin(r".*\.(fa|fasta)$", genome)
        ! verbose_option ? print("\r") : println("------------")
        println("Compute jellyfish on Genome              ")
        jf_genome = replace(basename(genome), ".fa" => ".jf")
        run(`jellyfish count -C -m $kmer_length -s 10000 -t $nbthreads -o $jf_dir/$jf_genome $genome`)
        if verbose_option println("\nGenome kmer index output: $jf_dir/$jf_genome.") end
        return(jf_genome, jf_dir)
    else
        println("Jellyfish index provided, nothing to do")

        jf_genome = basename(genome)
        jf_dir = dirname(genome)

        return(jf_genome, jf_dir)

    end
end


function APPRIS_function(gene_ref)
    ! verbose_option ? print("\r") : println("------------")
    print("Find best isoform on APPRIS database for $gene_ref ")
    #### Request to appris
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
    
    #### find best isoform
    ## find Principals
    principals = []
    for (i,value) in enumerate(transcripts)
        try 
            if occursin("PRINCIPAL:", transcripts[i]["reliability"])
                push!(principals, transcripts[i]["reliability"]) 
            end
        catch
        end
    end
    ## if not principals, return NODATA
    if isempty(principals)
        if verbose_option
            println("No principal isoforms detected, this function will return " * 
                        "'NODATA' and the longest transcript will be selected")
        end
        return(transcripts)
    end
    if verbose_option println("principals: $principals") end
    ## Select better Principal (minimum level)
    levels = []
    for i in 1:length(principals)
        push!(levels, split(principals[i],":")[2]) 
    end
    levels = map(x-> parse(Int, x), levels)
    level = minimum(levels)
    ## if multiple transcripts with better Principal, select biggest
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
    if verbose_option println("APPRIS result : $best_transcript") end
    if verbose_option println("APPRIS function finished -----") end
    return(best_transcript)
end # end of APPRIS function


function find_longest_variant(gene_name, transcriptome_dict)
    if verbose_option println("-----\nStarting find_longest_variant function") end
    gene_name = replace(gene_name, "@SLASH@" => "/")
    nb_variants = length(filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict))
    if verbose_option println("nb variants : $nb_variants") end
    variants_dict = filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict)
    variants_lengths = Vector{Any}()
    for (k,v) in variants_dict
        #if verbose_option println(k, ": ", v[1:35], "...") end # bug if sequence length < 35
        push!(variants_lengths, length(v))
    end
    if verbose_option println("variants lengths : $variants_lengths") end
    longest_variant = String(collect(keys(filter(v -> length(v.second) == maximum(variants_lengths), variants_dict)))[1])
    ## sometimes, gene name contain a slash! 
    gene_name = replace(gene_name, "/" => "@SLASH@")
    if verbose_option println("$gene_name : longest variant = $longest_variant") end
    if verbose_option println("find_longest_variant function finished") end
    return(longest_variant)
end


function build_sequences()
    if verbose_option println("----------\nBegining Build_sequences function") end
    ## split each sequence of fasta input file into individual fastas
    ## WARNING works only with ensembl fasta descriptions!
    if ! unannotated_option
        run(`mkdir -p "$output/sequences/$kmer_length"`)
        genes_already_processed = []
        genes_analysed = []
        FastaReader(fastafile) do fr
            for (desc, seq) in fr
                desc_array = split(desc)
                APPRIS_transcript = ""
                # if verbose_option println("desc_array: $desc_array") end
                gene_name = replace(desc_array[7], "gene_symbol:" => "")
                ensembl_transcript_name = split(desc_array[1],'.')[1]
                ensembl_gene_name = split(replace(desc_array[4], "gene:" => ""), '.')[1]
                # println("ensembl gene name : $ensembl_gene_name")
                ## engage this loop only if the gene have not been already processed (otherwise, the gene will be processed n times, n being the number of transcripts of the same genes in the reference fasta)
                if isempty(select_option) || (select_option != nothing && (gene_name in select_option || ensembl_transcript_name in select_option || ensembl_gene_name in select_option)) && (gene_name in genes_already_processed) == false
                    ## take all sequence corresponding to asked gene names (select option)
                    ## take all if select_option not provided
                    if APPRIS_option != nothing && (gene_name in genes_analysed) == false
                        APPRIS_transcript = APPRIS_function(ensembl_gene_name) 
                        # println("APPRIS selected variant : $APPRIS_transcript")
                        push!(genes_analysed, gene_name)
                    end
                    # println("$ensembl_transcript_name")
                    # println("marqueur ici : $gene_name")
                    if APPRIS_option == nothing || (APPRIS_option != nothing && ensembl_transcript_name == APPRIS_transcript) || (APPRIS_option != nothing && APPRIS_transcript == "NODATA" && "$gene_name:$ensembl_transcript_name" == find_longest_variant(gene_name,transcriptome_dict)) 
                        if verbose_option println("ensembl_transcript_name : $ensembl_transcript_name") end
                        if length("$seq") >= kmer_length 
                            if verbose_option println("$gene_name:$ensembl_transcript_name : good length, continue") end
                            gene_name = replace(gene_name, "/" => "@SLASH@") # some genes names can caontain some slash characters, break the processus
                            FastaWriter("$output/sequences/$kmer_length/$gene_name:$ensembl_transcript_name.fa") do fwsequence
                                write(fwsequence, [">$gene_name:$ensembl_transcript_name", "$seq"])
                                push!(genes_already_processed, gene_name)
                            end
                        else 
                            println("$gene_name:$ensembl_transcript_name : wrong length, ignored")
                        end
                    end
                end
            end
            # println("genes analysed : $genes_analysed")
            for i in unique(genes_analysed)
                # println("in genes already expressed : $(i in genes_already_processed)")
                ## it can happen that APPRIS contain obsolete id or other error, in this case the longer variant is selected for the concerned genes (passed in loop but have not provided a "sequence" fasta)
                if (replace(i, "/" => "@SLASH@") in genes_already_processed) == false
                    if verbose_option
                        println("""WarningAppris : There is a problem with $i, 
                            \rAPPRIS function have not worked properly, sometimes 
                            \ran outdated APPRIS reference, in this case the longer 
                            \rtranscript is selected from the reference fasta""")
                    end
                    id = find_longest_variant(i, transcriptome_dict)
                    # println(transcriptome_dict["$id"])
                    FastaWriter("$output/sequences/$kmer_length/$(replace(id, "/" => "@SLASH@"))") do fwsequence
                        write(fwsequence, [">$id", string(transcriptome_dict["$id"])])
                    end
                end
            end
        end
    ## if the input fasta is not a part of ensembl annotation (-u)
    else
        FastaReader("$fastafile") do fr
            for (desc, seq) in fr
                gene_name = desc
                if length("$seq") >= kmer_length
                    FastaWriter("$output/sequences/$kmer_length/$desc.fa") do fwsequence
                        write(fwsequence, [">$desc", "$seq"])
                    end
                else 
                    println("$desc ========> wrong length, ignored")
                end
            end
        end
    end
    if verbose_option
        println("Sequences:")
        for seq in (readdir("$output/sequences/$kmer_length/"))
            println(" $seq")
        end
    end
    if verbose_option println("Build_sequences function finished-----") end
end



parsed_args = parse_commandline()
set_variables()
add_workers()
@everywhere using Dates
@everywhere using DelimitedFiles
@everywhere using FastaIO
using  ParallelDataTransfer
checkup_variables()
transcriptome_dict = load_transcriptome(transcriptome)
jf_genome, jf_dir = run_jellyfish(genome, transcriptome)
build_sequences()


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

## function for extraction of specifics kmers for one sequence fasta file.
@everywhere function f(splitted_fasta_files)
    fasta_array = Array([])
    myid() == 1 ? X = myid() : X = myid()-1
    elapsed_time = 0
    step_time = 0
    kmers_analysed = 0
    replace_line = function(nbline::Int64, x::String)
        print(string("\u1b[$(nbline)F\u1b[2K"*"$x"*"\u1b[$(nbline)E"))
    end
    if ! unannotated_option          # if annotated
        gene_name, transcript_name = split(splitted_fasta_files, ":")
        nb_variants = length(filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict))
        if verbose_option println("gene_name: $gene_name --- transcript_name: $transcript_name --- nb variants : $nb_variants") end
        variants_dict = filter(k -> startswith(k.first, "$gene_name:"), transcriptome_dict)
    else                            # if unannotated
        gene_name = transcript_name = "$splitted_fasta_files"
        if verbose_option println("gene_name = transcript_name: $transcript_name") end
    end
    ## Determine tag-file 
    if level == "gene"
        tag_file = "$gene_name-$transcript_name-gene_specific.fa"
    elseif level == "transcript"
        tag_file = "$gene_name-$transcript_name-transcript_specific.fa"
    elseif level == "chimera"
        tag_file = "$gene_name-chimera_specific.fa"
    end
    # take the transcript sequence for jellyfish query
    sequence_fasta = "error if seen"
    FastaReader("$output/sequences/$kmer_length/$splitted_fasta_files") do fr
        for (desc, seq) in fr
            sequence_fasta = seq
        end
    end
    if verbose_option remotecall(replace_line, 1 , X, "worker number $(myid()), for $level $gene_name, reporting for duty ! :)\n") end
    ## build kmer count dictionnary from genome's jellyfish
    println("jellyfish query -s $output/sequences/$kmer_length/$splitted_fasta_files $jf_dir/$jf_genome")
    kmercounts_genome = read(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$jf_dir/$jf_genome"`, String)
    if verbose_option remotecall(replace_line, 1 , X, "jellyfish query $splitted_fasta_files on $jf_dir/$jf_genome finished") end
    kmercounts_genome = split(kmercounts_genome, "\n")[1:end-1]
    kmercounts_genome_dict = Dict()
    for mer in kmercounts_genome
        mer = split(mer)
        counts = mer[2]
        seq = mer[1]
        kmercounts_genome_dict["$seq"] = counts
    end
    ## build kmer count dictionnary from transcriptome
    kmercounts_transcriptome = read(`jellyfish query -s "$output/sequences/$kmer_length/$splitted_fasta_files" "$transcriptome"`, String)
    if verbose_option remotecall(replace_line, 1 , X, "jellyfish query on transcriptome.jf finished") end
    kmercounts_transcriptome_dict = Dict()
    kmercounts_transcriptome = split(kmercounts_transcriptome, "\n")[1:end-1]
    for mer in kmercounts_transcriptome
        mer = split(mer)
        counts = mer[2]
        seq = mer[1]
        kmercounts_transcriptome_dict["$seq"] = counts
    end
    ## initialisation of count variables
    i = 0
    total_kmers = length(keys(kmercounts_transcriptome_dict))
    
    ## create a new dictionary containing starts of kmers coordinates, to index
    ## kmers by their place on the sequence
    kmer_starts = Dict()
    kmer_placed = 0
    for mer in keys(kmercounts_transcriptome_dict)
        kmer_placed = kmer_placed +1
        if verbose_option remotecall(replace_line, 1 , X, "kmers placed = $kmer_placed on $total_kmers") end
        kmer_interval = findfirst("$mer", "$sequence_fasta")
        kmer_start = minimum(collect(kmer_interval))
        kmer_starts["$mer"] = kmer_start
    end
    
    for tuple in sort(collect(zip(values(kmer_starts),keys(kmer_starts))))
        ## for now, I can only have a sorted array of tuples, I have to extract 
        ## sequence for each tuple
        mer = last(tuple)
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
        remotecall(replace_line, 1 , X, "$level $gene_name : $kmers_analysed kmers analysed ( $per %)  ,$i specifics,:remaining time -> $ptime sec.")
        
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
            if genome_count == "1" || genome_count == "0" 
                mer_regex = Regex(mer)
                variants_containing_this_kmer = keys(filter(v -> occursin(mer_regex, v.second), variants_dict))
                if stringent_option == true && Float64(transcriptome_count) == Float64(nb_variants) == Float64(length(variants_containing_this_kmer)) 
                    i = i+1
                    tmp = length(variants_containing_this_kmer)
                    push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)")
                    push!(fasta_array,"$mer")
                elseif stringent_option == false && Float64(transcriptome_count) == Float64(length(variants_containing_this_kmer)) && Float64( transcriptome_count) > Float64( nb_variants*admission_threshold)
                    i = i+1
                    tmp = length(variants_containing_this_kmer)
                    push!(fasta_array,">$gene_name-$transcript_name.kmer$i ($tmp/$nb_variants)" )
                    push!(fasta_array, "$mer")
                elseif unannotated_option == true && Float64(transcriptome_count) == Float64(0)
                    i = i+1
                    push!(fasta_array,">$gene_name-$transcript_name.kmer$i" )
                    push!(fasta_array,"$mer")
                end # end of fasta writing function
            end
        end
        if level == "transcript"
            if unannotated_option == true && Float64(transcriptome_count) == Float64(0) && parse(Int, genome_count) <= 1
                i = i+1
                push!(fasta_array,">$gene_name.kmer$i")
                push!(fasta_array,"$mer")
            elseif Float64(transcriptome_count) == Float64(1) && parse(Int, genome_count) <= 1
                i = i+1
                push!(fasta_array,">$gene_name-$transcript_name.kmer$i")
                push!(fasta_array,"$mer")
            end
        end
        if level == "chimera" && unannotated_option == true
            if Float64(transcriptome_count) == Float64(0) && parse(Int, genome_count) == 0
                i = i+1
                push!(fasta_array,">$gene_name.kmer$i")
                push!(fasta_array,"$mer")
            end
        end
        step_time = Float64(time()-startt)
    end # end of kmer analysis
    
    ## write tag file
    run(`mkdir -p "$output/tags/$kmer_length"`)
    writedlm("$output/tags/$kmer_length/$tag_file",
    reshape(fasta_array,length(fasta_array)), "\n")
end #end of function

print("\rExtract specific kmers                                  ")
for n in 1:nworkers()+1 println("") end
pmap(f, readdir("$output/sequences/$kmer_length/"))



