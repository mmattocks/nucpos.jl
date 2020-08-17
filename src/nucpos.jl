module nucpos

using DataFrames, BioSequences, BioBackgroundModels, CSV, FASTX

function add_position_sequences!(df::DataFrame, genome_path::String, genome_idx_path::String)
    scaffold_seq_record_dict::Dict{String,BioSequences.LongSequence} = BioBackgroundModels.build_scaffold_seq_dict(genome_path, genome_idx_path)

    seqs=[BioSequences.LongSequence{DNAAlphabet{4}}() for i in 1:size(df,1)]
    df[!, :seq] .= seqs

    Threads.@threads for entry in eachrow(df)
        entry.seq=BioBackgroundModels.fetch_sequence(entry.chr, scaffold_seq_record_dict ,entry.start, entry.end, '+')
    end

end

function get_cluster!(df::DataFrame, fasta::String, arch::String)
    df[!, :cluster] .= zeros(Int64, size(df,1))
    art=CSV.read(arch, header=0)

    position_vec=Vector{Vector{Int64}}()
    reader=FASTA.Reader(open(fasta))
    for (n, entry) in enumerate(reader)
        scaffold = FASTA.identifier(entry)
        desc_array = split(FASTA.description(entry))
        pos_start = parse(Int64, desc_array[2])
        pos_end = parse(Int64, desc_array[4])
        seq = FASTA.sequence(entry)

        idx = filter(in(findall(chr->chr==scaffold, df.chr)), findall(start->start==pos_start, df.start))

        if length(idx)==1
            match=df[idx[1],:]
            @assert match.end == pos_end
            @assert match.seq==seq
            match.cluster=art.Column1[n]
        end
    end
end


function make_position_df(position_fasta::String)
    position_reader = FASTA.Reader(open((position_fasta),"r"))
    position_df = DataFrame(SeqID = String[], Start=Int64[], End=Int64[], Seq = LongSequence[])

    for entry in position_reader
        scaffold = FASTA.identifier(entry)

        if scaffold != "MT"
            desc_array = split(FASTA.description(entry))
            pos_start = parse(Int64, desc_array[2])
            pos_end = parse(Int64, desc_array[4])
            seq = FASTA.sequence(entry)

            if !hasambiguity(seq)
                push!(position_df, [scaffold, pos_start, pos_end, seq])
            end
        end
    end
    
    close(position_reader)
    return position_df
end

function observation_setup(position_df::DataFrame; order::Int64=0, symbol::Symbol=:Seq)
    order_seqs = BioBackgroundModels.get_order_n_seqs(Vector{LongSequence{DNAAlphabet{2}}}(position_df[!, symbol]),order)
    coded_seqs = BioBackgroundModels.code_seqs(order_seqs)

    return coded_seqs
end

function map_positions!(base_df, map_df)
    df_rows = size(base_df,1)
    base_df[!, :mapped_pos] = zeros(Int64, df_rows)
    base_df[!, :overlap_bp] = zeros(Int64, df_rows)
    base_df[!, :rel_overlap] = zeros(df_rows)

    chrgroups=groupby(base_df, :chr)

    for chr_df in chrgroups #split the df by scaffold so base numbers are local to scaffold
        base_chr = chr_df.chr[1]
        matched=false
        for map_chr_df in groupby(map_df, :chr)
            base_chr == map_chr_df.chr[1] && subDFmap!(chr_df, map_chr_df)
            matched=true
        end
        matched=false && (chr_df.mapped_pos=zeros(Int64,length(chr_df.start));
                          chr_df.overlap_bp=zeros(Int64,length(chr_df.start));
                          chr_df.rel_overlap=zeros(length(chr_df.start)))
    end

    return base_df
end

                function subDFmap!(chr_df::SubDataFrame, map_chr_df::SubDataFrame)
                    Threads.@threads for position in eachrow(chr_df)
                        search_start = position.start; search_end=position.end
                        start_idxs=findall(in(findall(ende->search_start<=ende,map_chr_df.end)), findall(start->search_start>=start,map_chr_df.start))
                        end_idxs=findall(in(findall(ende->search_end<=ende,map_chr_df.end)), findall(start->search_end>=start,map_chr_df.start))
                        mapped_idxs = unique(vcat(start_idxs, end_idxs))
                    
                        overlaps=Vector{Int64}()
                        for idx in mapped_idxs
                            #println(idx)
                            mapped_start=map_chr_df.start[idx]; mapped_end=map_chr_df.end[idx]
                            if (mapped_start==search_start && mapped_end==search_end) #complete position overlap (all positions are the same size)
                                push!(overlaps,(mapped_end-mapped_start+1))
                                #println("same")
                            elseif (mapped_start < search_start <= mapped_end) #5' end of the searched position overlaps with the mapped position
                                push!(overlaps,(mapped_end-search_start+1))
                                #println("5' overlap")
                            elseif (mapped_start <= search_end < mapped_end) #3' end of searched position overlaps with the mapped position
                                push!(overlaps,(search_end-mapped_start+1))
                                #println("3' overlap")
                            end
                        end

                        if length(mapped_idxs)>0 #1 or more mapped candidates: take the one with the most overlap
                            overlap, oidx = findmax(overlaps)
                            map_idx = mapped_idxs[oidx]
                            position.mapped_pos=map_idx
                            position.overlap_bp=overlap
                            position.rel_overlap=overlap / (search_end-search_start+1)
                        end
                    end
                end

end # module
