#= from DoViPv0.9 and up, the detection of the contig shape is done in the checkV integrated ("export_aggregated_int" function) 
and CheckV nonintegrated ("merge_nonintegrated" function) =#

function post_genomad(proj::ProjGenomad, parentD::String)
    #all
    if proj.ext_res == true
        df =  CSV.read("$(proj.ext_res_D)/$(proj.genomad_out_table_p.p)", DataFrame; delim = '\t', header = 1)
    else
        df =  CSV.read("$(parentD)/$(proj.genomad_out_table_p.p)", DataFrame; delim = '\t', header = 1)
    end

    df[!, :contig_name] = Vector{Union{Missing,String}}(missing, nrow(df)) # ! modifies in place the dataframe
    df[!, :predictor_genomad] = fill("yes", nrow(df))
    df[!, :virus_type_genomad] = Vector{Union{Missing,String}}(missing, nrow(df))
    df = rename!(df, :virus_score => :virus_score_genomad, :topology => :topology_genomad, :taxonomy => :taxonomy_genomad)
    
    for row in 1:nrow(df)
        splits = split(df[row, :seq_name], "|")
        df[row, :contig_name] = splits[1]
        if length(splits) >= 2 && occursin("provirus", splits[2])
            df[row, :virus_type_genomad] = "Integrated"
        else
            df[row, :virus_type_genomad] = "Non-integrated"
        end
    end
    
    df = select!(df, [:seq_name, :contig_name, :predictor_genomad, :virus_type_genomad, :length, :coordinates, :virus_score_genomad, :topology_genomad, :taxonomy_genomad])
 
    #nonintegrated
    nonintegrated_df = subset(df, :virus_type_genomad => a -> a .== "Non-integrated")
    nonintegrated_df = select!(nonintegrated_df, [:contig_name, :predictor_genomad, :length, :virus_score_genomad, :topology_genomad, :taxonomy_genomad])
     

    if !isempty(nonintegrated_df)
        CSV.write("$(parentD)/$(proj.postgenomad_nonintegrated_df.p)", nonintegrated_df, delim = '\t', header = true)
    end


    #integrated
    integrated_df = subset(df, :virus_type_genomad => a -> a .== "Integrated")
    integrated_df[!, :provirus_start] = Vector{Union{Missing, Int64}}(missing, nrow(integrated_df))
    integrated_df[!, :provirus_end] = Vector{Union{Missing, Int64}}(missing, nrow(integrated_df))

    for row in 1:nrow(integrated_df)
        splits = split(integrated_df[row, :coordinates], "-")
        integrated_df[row, :provirus_start] = parse(Int64, splits[1])
        integrated_df[row, :provirus_end] = parse(Int64, splits[2])
    end

    integrated_df = select!(integrated_df, [:seq_name, :contig_name, :predictor_genomad, :length, :provirus_start, :provirus_end, :virus_score_genomad, :taxonomy_genomad, :topology_genomad])

    if !isempty(integrated_df)
        CSV.write("$(parentD)/$(proj.postgenomad_integrated_df.p)", integrated_df, delim = '\t', header = true)
    end
 
    return nothing
end

function post_DVF(proj::ProjDVF, parentD::String)
    if proj.ext_res == true
        df =  CSV.read("$(proj.ext_res_D)/$(proj.modif_output_dvf_f.p)", DataFrame; delim = '\t', header = 1)
    else
        df =  CSV.read("$(parentD)/$(proj.modif_output_dvf_f.p)", DataFrame; delim = '\t', header = 1)
    end

    df = rename!(df, :name => :contig_name, :len => :length, :score => :score_dvf, :pvalue => :pvalue_dvf)
    df[!, :predictor_dvf] = fill("yes", nrow(df))

    if !isempty(df)
        CSV.write("$(parentD)/$(proj.postdvf_nonintegrated_df.p)", df, delim = '\t', header = true)
    end

    return nothing
end

# unussed function, just in case I will need it later
function get_VS2_circ_length(f_p::FnaP)
    to_split = get_contig_names(f_p)

    dict = Dict{String, Tuple{Int64, Int64}}()

    for i in eachindex(to_split)
        splittedA = split(to_split[i], " ")
        contig_name = splittedA[1]

        splitB = split(splittedA[3], "||")

        splitC = split(splitB[2], ":")
        contig_start = parse(Int64, splitC[2])

        splitD = split(splitB[3], ":")
        contig_end = parse(Int64, splitD[2])

        dict[contig_name] = (contig_start, contig_end)
    end

    return dict
end

function post_virSorter2(proj::ProjVirSorter2, parentD::String)
    #all
    if proj.ext_res == true
        dfscore =  CSV.read("$(proj.ext_res_D)/$(proj.vs2_viral_score_f.p)", DataFrame; delim = '\t', header = 1)
    else
        dfscore =  CSV.read("$(parentD)/$(proj.vs2_viral_score_f.p)", DataFrame; delim = '\t', header = 1)
    end
    select!(dfscore, [:seqname, :max_score, :max_score_group, :length, :hallmark])
    rename!(dfscore, :seqname => :seqname_new, :max_score => :max_score_virSorter2, :max_score_group => :max_score_group_virSorter2, :hallmark => :hallmark_virSorter2)

    if proj.ext_res == true
        dfboundary =  CSV.read("$(proj.ext_res_D)/$(proj.vs2_viral_boundary_f.p)", DataFrame; delim = '\t', header = 1)
    else
        dfboundary =  CSV.read("$(parentD)/$(proj.vs2_viral_boundary_f.p)", DataFrame; delim = '\t', header = 1)
    end
    select!(dfboundary, [:seqname, :seqname_new, :shape, :partial, :trim_bp_start, :trim_bp_end])
    rename!(dfboundary, :seqname => :contig_name, :shape => :shape_virSorter2)

    df = innerjoin(dfscore, dfboundary, on = :seqname_new)
    df[!, :predictor_virSorter2] = fill("yes", nrow(df))
 
    #nonintegrated
    nonintegrated_df = subset(df, :partial => a -> a .== 0)
    rename!(nonintegrated_df, :length => :length_virSorter2)
    select!(nonintegrated_df, [:contig_name, :predictor_virSorter2, :length_virSorter2, :max_score_virSorter2, :max_score_group_virSorter2, :hallmark_virSorter2, :shape_virSorter2])

    if !isempty(nonintegrated_df)
        CSV.write("$(parentD)/$(proj.postvs2_nonintegrated_df.p)", nonintegrated_df, delim = '\t', header = true)
    end

    #integrated
    integrated_df = subset(df, :partial => a -> a .== 1)
    
    rename!(integrated_df, :seqname_new => :seq_name, :trim_bp_start => :provirus_start, :trim_bp_end => :provirus_end)
    select!(integrated_df, [:seq_name, :contig_name, :predictor_virSorter2, :length, :max_score_virSorter2, :max_score_group_virSorter2, :hallmark_virSorter2, :shape_virSorter2, :provirus_start, :provirus_end])
    
    if !isempty(integrated_df)
        #save
        CSV.write("$(parentD)/$(proj.postvs2_integrated_df.p)", integrated_df, delim = '\t', header = true)
    end

    return nothing
end

function post_vibrant_nonintegrated!(proj::ProjVibrant, parentD::String)
    #nonintegrated viruses
    if proj.ext_res == true
        nonintegrated_df = DataFrame(contig_name = get_contig_names(FnaP("$(proj.ext_res_D)/$(proj.vib_out_nonintegrated_fna.p)")))
    else
        nonintegrated_df = DataFrame(contig_name = get_contig_names(FnaP("$(parentD)/$(proj.vib_out_nonintegrated_fna.p)")))
    end

    if !isempty(nonintegrated_df)
        nonintegrated_df[!, :predictor_vibrant] = fill("yes", nrow(nonintegrated_df))

        if proj.ext_res == true
            nonintegrated_df[!, :length] = get_contig_length(FnaP("$(proj.ext_res_D)/$(proj.vib_out_nonintegrated_fna.p)"))
        else
            nonintegrated_df[!, :length] = get_contig_length(FnaP("$(parentD)/$(proj.vib_out_nonintegrated_fna.p)"))
        end

        CSV.write("$(parentD)/$(proj.postvib_nonintegrated_df.p)", nonintegrated_df, delim = '\t', header = true)
    end
    return nothing
end

function post_vibrant_integrated!(proj::ProjVibrant, parentD::String)
    #integratederate viruses
    if proj.ext_res == true
        integrated_df = CSV.read("$(proj.ext_res_D)/$(proj.vib_out_integrated_tsv.p)", DataFrame; delim = '\t', header = 1)
    else
        integrated_df = CSV.read("$(parentD)/$(proj.vib_out_integrated_tsv.p)", DataFrame; delim = '\t', header = 1)
    end
    rename!(integrated_df, :fragment => :seq_name, :scaffold => :contig_name, Symbol("nucleotide start") => :provirus_start, 
                    Symbol("nucleotide stop") => :provirus_end, Symbol("nucleotide length") => :length)

    select!(integrated_df, Not(Symbol("protein start"), Symbol("protein stop"), Symbol("protein length")))
    
    integrated_df[!, :predictor_vibrant] = fill("yes", nrow(integrated_df))

    if !isempty(integrated_df)
        CSV.write("$(parentD)/$(proj.postvib_integrated_df.p)", integrated_df, delim = '\t', header = true)
    end

    return nothing
end

function post_viralVerify(proj::ProjViralVerify, parentD::String)
    # it only predicts non-integrated viruses

    if proj.ext_res == true
        vv_df = CSV.read("$(proj.ext_res_D)/$(proj.viralVerify_out_p.p)", DataFrame, silencewarnings=true; delim = ',', header = 1)
    else
        vv_df = CSV.read("$(parentD)/$(proj.viralVerify_out_p.p)", DataFrame, silencewarnings=true; delim = ',', header = 1)
    end

    if !isempty(vv_df)
        if eltype(vv_df[!, :Score]) != Float64
            vv_df[!, :score_viralVerify] = Vector{Union{Missing, Float64}}(missing, nrow(vv_df))
            for i in 1:nrow(vv_df)
                if(vv_df[i, :Score] == "-")
                    vv_df[i, :score_viralVerify] = 0
                else
                    vv_df[i, :score_viralVerify] = parse(Float64, vv_df[i, :Score])
                end
            end
        else 
            rename!(vv_df, :Score => :score_viralVerify)
        end

        rename!(vv_df, Symbol("Contig name") => :contig_name, :Length => :length, :Circular => :shape_viralVerify)
        subset!(vv_df, :score_viralVerify => a -> a .>= proj.score_threshold, :Prediction => a -> a .== "Virus")
        select!(vv_df, :contig_name, :length, :score_viralVerify, :shape_viralVerify)
        vv_df[!, :predictor_viralVerify] = fill("yes", nrow(vv_df))

        if !isempty(vv_df)
            CSV.write("$(parentD)/$(proj.postViralVerify_df.p)", vv_df, delim = '\t', header = true)
        end
    end

    return nothing
end

#region INTEGRATED - POST PREDICTION AND CHECKV

function export_aggregated_int(inref::FnaP, proj::ProjCheckVIntegrated, parentD::String)

    dfj = DataFrame()

    for i in 1:length(proj.input_dfs_2_aggregate)
        df = CSV.read("$(parentD)/$(proj.input_dfs_2_aggregate[i].p)", DataFrame; delim = '\t', header = 1)

        dfj = vcat(dfj, df, cols=:union)
    end

    dfj[!, :virus_type_DoViP] = fill("Integrated", nrow(dfj))
    dfj[!, :provirID1] = fill("", nrow(dfj))

    # determine the shape and ends of the contigs with INTEGARTED virus regions and join it to the out DF
    out_contig_fna_p = FnaP("$(parentD)/$(proj.output_aggregated_int_contigs_fna.p)")
    contigNameSel(FnaP("$(parentD)/$(inref.p)"), out_contig_fna_p, unique(dfj[!, :contig_name]))   
    shape_ends_DF = detectContigShapeEnds(out_contig_fna_p)
    dfj = leftjoin!(dfj, shape_ends_DF, on = [:contig_name => :contig_name])
    
    #region viral regions that wrap around cicrular contig ends as determined by VirSorter2

    # split prophages which wrap around contig ends (for those contigs ending in DTRs, as detected by VirSorter2 
    dfj[!, :splitprovir] = Vector{Union{Missing, String}}(missing, nrow(dfj))
    dfj[!, :r_provir_start] = Vector{Union{Missing, Int64}}(missing, nrow(dfj))
    dfj[!, :r_provir_end] = Vector{Union{Missing, Int64}}(missing, nrow(dfj))

    # mark prophages which wrap themselves around a circular contig - meaning that they are overlapping with the left DTR
    #if "shape_virSorter2" in names(dfj) 
        for i in 1:nrow(dfj)
            if dfj[i, :contig_shape] == "circular" #!ismissing(dfj[i, :shape_virSorter2]) && dfj[i, :shape_virSorter2] == "circular" && dfj[i, :contig_shape] == "circular"
                if dfj[i, :provirus_end] > dfj[i, :contig_trimmed_end]  
                    dfj[i, :r_provir_start] = 1
                    dfj[i, :r_provir_end] = dfj[i, :provirus_end] - dfj[i, :contig_trimmed_end] 
                    dfj[i, :provirus_end] = dfj[i, :contig_trimmed_end] 
                    dfj[i, :splitprovir] = dfj[i, :seq_name]
                    dfj[i, :length] = dfj[i, :r_provir_end] + (dfj[i, :provirus_end] - dfj[i, :provirus_start] + 1)
                end
            elseif "shape_virSorter2" in names(dfj) && !ismissing(dfj[i, :shape_virSorter2]) && dfj[i, :shape_virSorter2] == "circular" && dfj[i, :contig_shape] == "linear" #
                if dfj[i, :provirus_end] > dfj[i, :contig_full_end]  
                    dfj[i, :provirus_end] = dfj[i, :contig_full_end]  
                end
            end
        end
    #end
    
    #= find corresponding contig names for the split proviruses (as detected by VirSorter2) and save them as a Set
    splitprovir_set = collect(skipmissing(unique(dfj[!, :splitprovir])))
    
    for i in eachindex(splitprovir_set)
        splitprovir_set[i] = replace(splitprovir_set[i], Regex("\\|\\|\\d+_partial") => "")
    end
    splitprovir_set = Set(unique(splitprovir_set)) =#

    for i in 1:nrow(dfj)
        dfj[i, :provirID1] = "provirID1num$(i)"

        #= trim prophages predicted by the other predictors found on circular contigs that extend in the right direct terminal repeat
        end, contigs for which virSorter2 found proviruses wrapping around the ends of contigs #
        # check if the value in column seq_name is in the splitprovir_vec 
        if dfj[i, :contig_name] in splitprovir_set
            if dfj[i, :provirus_end] > dfj[i, :contig_trimmed_end]
                dfj[i, :provirus_end] = dfj[i, :contig_trimmed_end]
                dfj[i, :length] = dfj[i, :provirus_end] - dfj[i, :provirus_start] + 1
            end
        end =#

    end
    #endregion 

    # save proviruses based on two sets of coordinates
    out_fna_p = FnaP("$(parentD)/$(proj.output_aggregated_fna.p)")
    colnames = names(dfj)
    if ("r_provir_start" in colnames && "r_provir_end" in colnames)
        contigExtractRegions(FnaP("$(parentD)/$(inref.p)"), out_fna_p, dfj[!, :contig_name], 
                                dfj[!, :provirID1], dfj[!, :provirus_start], dfj[!, :provirus_end], dfj[!, :r_provir_start], dfj[!, :r_provir_end])
    else
        contigExtractRegions(FnaP("$(parentD)/$(inref.p)"), out_fna_p, dfj[!, :contig_name], 
                                dfj[!, :provirID1], dfj[!, :provirus_start], dfj[!, :provirus_end])
    end

    # save
    CSV.write("$(parentD)/$(proj.output_aggregated_df.p)", dfj, delim = '\t', header = true)

    return dfj
end

function return_provir_CheckV_coord(fna_p::String, parentD::String)
    tosplit = get_contig_names(FnaP("$(parentD)/$(fna_p)"))

    df = DataFrame(provir_names_checkV = Vector{Union{Missing, String}}(missing, length(tosplit)),
                    provirID1 = Vector{Union{Missing, String}}(missing, length(tosplit)),
                    provir_start_checkV = Vector{Union{Missing, Int64}}(missing, length(tosplit)),
                    provir_end_checkV = Vector{Union{Missing, Int64}}(missing, length(tosplit)))

    for i in eachindex(tosplit)
        splittedA = split(tosplit[i], " ")

        splittedA[1] = replace(splittedA[1], ">" => "")
        df[i, :provir_names_checkV] = splittedA[1]
        
        df[i, :provirID1] = split(splittedA[1], "_")[1]

        splittedB = split(splittedA[2], "/")[1]
        coord = split(splittedB, "-")
        df[i, :provir_start_checkV] = parse(Int64, coord[1])
        df[i, :provir_end_checkV] = parse(Int64, coord[2])
    end

    return(df)
end

function post_checkv1_integrated!(proj::ProjCheckVIntegrated, agreg_df::DataFrame, parentD::String)
    
     # find start and ends of prophages from checkV output
    if filesize("$(parentD)/$(proj.checkV1_out_provir_fna.p)") >0  #find if checkV detected any prophages at all
        df = return_provir_CheckV_coord(proj.checkV1_out_provir_fna.p, parentD)

        # merge checkV info with the agregated integrated_df
        #agregDFj = leftjoin!(agreg_df, df, on = [:provirID1])
        agregDFj = outerjoin(agreg_df, df, on = [:provirID1])    #hm, here I placed the df first, in the case when there are multiple prophages predicted by checkV. But this way, I'm loosing those proviruses without end contamination
    else
        agregDFj = agreg_df
    end
    
    agregDFj[!, :provir_start_cor] = Vector{Union{Missing, Int64}}(missing, nrow(agregDFj))
    agregDFj[!, :provir_end_cor] = Vector{Union{Missing, Int64}}(missing, nrow(agregDFj))
    agregDFj[!, :r_provir_start_cor] = Vector{Union{Missing, Int64}}(missing, nrow(agregDFj))
    agregDFj[!, :r_provir_end_cor] = Vector{Union{Missing, Int64}}(missing, nrow(agregDFj))

    for i in 1:nrow(agregDFj)
        if (("provir_start_checkV" in names(agregDFj) && ismissing(agregDFj[i, :provir_start_checkV]) == false) && 
            ("provir_end_checkV" in names(agregDFj) && ismissing(agregDFj[i, :provir_end_checkV]) == false))
            interm_start = (agregDFj[i, :provirus_start] + agregDFj[i, :provir_start_checkV] - 1)
            interm_end = (agregDFj[i, :provirus_start] + agregDFj[i, :provir_end_checkV] - 1)

            if ("splitprovir" in names(agregDFj) && ismissing(agregDFj[i, :splitprovir]) == false)  # this corrects the proviruses spanning the ends of circular contigs as predicted by virSorter2
                if interm_start < agregDFj[i, :contig_trimmed_end] && interm_end <= agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :provir_start_cor] = interm_start
                    agregDFj[i, :provir_end_cor] = interm_end
                    agregDFj[i, :splitprovir] = missing
                elseif interm_start > agregDFj[i, :contig_trimmed_end] && interm_end > agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :provir_start_cor] = interm_start - agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :provir_end_cor] = interm_end - agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :splitprovir] = missing
                elseif interm_start <= agregDFj[i, :contig_trimmed_end] && interm_end > agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :provir_start_cor] = interm_start
                    agregDFj[i, :provir_end_cor] = agregDFj[i, :contig_trimmed_end]
                    agregDFj[i, :r_provir_start_cor] = 1
                    agregDFj[i, :r_provir_end_cor] = interm_end - agregDFj[i, :contig_trimmed_end]
                end
            else
                agregDFj[i, :provir_start_cor] = interm_start
                agregDFj[i, :provir_end_cor] = interm_end
            end

        elseif (!("provir_start_checkV" in names(agregDFj)) || ismissing(agregDFj[i, :provir_start_checkV])) && 
            (!("provir_end_checkV" in names(agregDFj)) || ismissing(agregDFj[i, :provir_end_checkV]))
            agregDFj[i, :provir_start_cor] = agregDFj[i, :provirus_start]
            agregDFj[i, :provir_end_cor] = agregDFj[i, :provirus_end]

            if ("splitprovir" in names(agregDFj) && ismissing(agregDFj[i, :splitprovir]) == false)
                agregDFj[i, :r_provir_start_cor] = agregDFj[i, :r_provir_start]
                agregDFj[i, :r_provir_end_cor] = agregDFj[i, :r_provir_end]
            end
        end
    end
        
    CSV.write("$(parentD)/$(proj.postcheckv1_integrated_cor_df.p)", agregDFj, delim = '\t', header = true)

    return agregDFj
end

function find_consensus(vc::AbstractVector)
    cons = unique(vc) |> skipmissing |> collect |> sort
        
    if length(cons) == 1
        return cons[1]
    elseif length(cons) > 1
        cons = join(cons, ", ")
        return cons
    elseif length(cons) == 0
        return missing
    end
end

function group_proviruses!(df::DataFrame, predictors::Vector{Symbol}, merge_circ_proph::Bool)
    # df = deepcopy(dfj)

    df[!, :provir] = fill("", nrow(df))
    gdf = groupby(df, :contig_name)

    #region merge overlapping regions
    newdf = DataFrame(
        provirus_name = Vector{String}(),
        contig_name = Vector{String}(),
        provirus_start = Vector{Int64}(),
        provirus_end = Vector{Int64}(),
        length = Vector{Int64}()
        )

    #find overlapping viral regions (on the same contig) and merge them
    for g in gdf
        counter = 1
        sort!(g, :provir_start_cor)

        g[1, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"
        provir_start = g[1, :provir_start_cor]
        provir_end = g[1, :provir_end_cor]

        for row in 2:nrow(g)
            if g[row, :provir_start_cor] >= provir_start && g[row, :provir_end_cor] <= provir_end
                g[row, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"
            elseif g[row, :provir_start_cor] <= provir_start && g[row, :provir_end_cor] >= provir_end
                provir_start = g[row, :provir_start_cor]
                provir_end = g[row, :provir_end_cor]
                g[row, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"
            elseif g[row, :provir_start_cor] <= provir_start && g[row, :provir_end_cor] <= provir_end && g[row, :provir_end_cor] >= provir_start
                provir_start = g[row, :provir_start_cor]
                g[row, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"
            elseif g[row, :provir_start_cor] >= provir_start && g[row, :provir_end_cor] >= provir_end && g[row, :provir_start_cor] <= provir_end
                provir_end = g[row, :provir_end_cor]
                g[row, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"         
            else
                # data in newdf 
                push!(newdf, (provirus_name = "$(g[1, :contig_name])___provirus_$(counter)", 
                        contig_name = g[1, :contig_name], 
                        provirus_start = provir_start,
                        provirus_end = provir_end,
                        length = (provir_end - provir_start + 1)))
                
                # new provirus
                counter += 1
                provir_start = g[row, :provir_start_cor]
                provir_end = g[row, :provir_end_cor]
                g[row, :provir] = "$(g[1, :contig_name])___provirus_$(counter)"
            end
        end

        # data in newdf
        push!(newdf, (provirus_name = "$(g[1, :contig_name])___provirus_$(counter)", 
        contig_name = g[1, :contig_name], 
        provirus_start = provir_start,
        provirus_end = provir_end,
        length = (provir_end - provir_start + 1)))

    end

    #create new columns to transfer the info from the df with unmerged pro-viruses to the df with merged proviruses
    newdf[!, :taxonomy] = Vector{Union{Missing, String}}(missing, nrow(newdf))

    for pred in predictors
        newdf[!, pred] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    end
    newdf[!, :splitprovir] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :r_provirus_start] = Vector{Union{Missing, Int64}}(missing, nrow(newdf))
    newdf[!, :r_provirus_end] = Vector{Union{Missing, Int64}}(missing, nrow(newdf))
    newdf[!, :shape_virSorter2] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :topology_genomad]= Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :shape_viralVerify] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :contig_shape] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :contig_end_repeat_type] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :contig_shape_remarks] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :contig_end_repeat_size] = Vector{Union{Missing, Int64}}(missing, nrow(newdf))
    newdf[!, :contig_trimmed_end] = Vector{Union{Missing, Int64}}(missing, nrow(newdf))
    newdf[!, :contig_full_end] = Vector{Union{Missing, Int64}}(missing, nrow(newdf))
    newdf[!, :todel] = Vector{Union{Missing, Bool}}(missing, nrow(newdf))

    #bring consesus info in newdf, which contains the merged viral regions
    ggdf = groupby(df, :provir)
    for i in 1:length(ggdf)
        for j in 1:nrow(newdf)
            if ggdf[i][1, :provir] == newdf[j, :provirus_name]
                if "taxonomy_genomad" in names(ggdf[i])
                    newdf[j, :taxonomy] = find_consensus(ggdf[i][!, :taxonomy_genomad])
                end
                
                if "shape_virSorter2" in names(ggdf[i])
                    newdf[j, :shape_virSorter2] = find_consensus(ggdf[i][!, :shape_virSorter2])
                end

                if "shape_viralVerify" in names(ggdf[i])
                    newdf[j, :shape_viralVerify] = find_consensus(ggdf[i][!, :shape_viralVerify])
                end

                if "topology_genomad" in names(ggdf[i])
                    newdf[j, :topology_genomad] = find_consensus(ggdf[i][!, :topology_genomad])
                end

                newdf[j, :contig_shape] = find_consensus(ggdf[i][!, :contig_shape])
                newdf[j, :contig_end_repeat_type] = find_consensus(ggdf[i][!, :contig_end_repeat_type])
                newdf[j, :contig_shape_remarks] = find_consensus(ggdf[i][!, :contig_shape_remarks])
                newdf[j, :contig_end_repeat_size] = find_consensus(ggdf[i][!, :contig_end_repeat_size])
                newdf[j, :contig_trimmed_end] = find_consensus(ggdf[i][!, :contig_trimmed_end])
                newdf[j, :contig_full_end] = find_consensus(ggdf[i][!, :contig_full_end])

                if "splitprovir" in names(ggdf[i])
                    newdf[j, :splitprovir] = find_consensus(ggdf[i][!, :splitprovir])
                end

                for p in predictors
                    if string(p) in names(ggdf)
                        newdf[j, p] = find_consensus(ggdf[i][!, p])
                    end
                end
            end
        end
    end
    #endregion


    #region merge prophages spanning (or very close to) the ends of a circular contig
    if merge_circ_proph == true
        anewdf = DataFrame()
        #arow = 0
        gnewdf = groupby(newdf, :contig_name)
        for gn in gnewdf     
            shape = find_consensus(gn[!, :contig_shape]) 

            # I will merge only pro-viruses found at the end of circular contigs and overlapp 
            if (!ismissing(shape) && shape == "circular") #&& !ismissing(split_vs)
                sort!(gn, :provirus_start)
                lastrow = nrow(gn)
                if gn[lastrow, :splitprovir] == gn[1, :splitprovir]
                    if (gn[lastrow, :contig_trimmed_end] - gn[lastrow, :provirus_end]) <= 1 
                        spaceL = gn[lastrow, :contig_trimmed_end] - gn[lastrow, :provirus_end]                       
                        if gn[1, :provirus_start] <= 1
                            spaceR = gn[1, :provirus_start]
                            if (spaceR + spaceL) <= 2

                                contig_end_repeat_type = find_consensus(gn[!, :contig_end_repeat_type])
                                contig_shape_remarks = find_consensus(gn[!, :contig_shape_remarks])
                                contig_end_repeat_size = find_consensus(gn[!, :contig_end_repeat_size])

                                colnames = names(gn)
                                if "r_provirus_start" in colnames
                                    r_provirus_start = 1
                                else
                                    r_provirus_start = missing
                                end

                                if "r_provirus_end" in colnames
                                    r_provirus_end = gn[1, :provirus_end]
                                else
                                    r_provirus_end = missing
                                end
            
                                if :predictor_genomad in predictors
                                    taxonomy = find_consensus([gn[lastrow, :taxonomy], gn[1, :taxonomy]])
                                    predictor_genomad = find_consensus([gn[lastrow, :predictor_genomad], gn[1, :predictor_genomad]])
                                    shape_g = find_consensus(gn[!, :topology_genomad])
                                else
                                    taxonomy = missing
                                    predictor_genomad = missing
                                    shape_g = missing
                                end

                    
                                if "splitprovir" in names(gn)
                                    splitprovir = find_consensus([gn[lastrow, :splitprovir], gn[1, :splitprovir]])
                                else
                                    splitprovir = missing
                                end

                                if :predictor_virSorter2 in predictors
                                    predictor_virSorter2 = find_consensus([gn[lastrow, :predictor_virSorter2], gn[1, :predictor_virSorter2]])
                                    shape_vs = find_consensus([gn[lastrow, :shape_virSorter2], gn[1, :shape_virSorter2]])
                                else
                                    predictor_virSorter2 = missing
                                    shape_vs = missing
                                end
                                    
                                if :predictor_vibrant in predictors
                                    predictor_vibrant = find_consensus([gn[lastrow, :predictor_vibrant], gn[1, :predictor_vibrant]])
                                else
                                    predictor_vibrant = missing
                                end

                                if :predictor_dvf in predictors
                                    predictor_dvf = find_consensus([gn[lastrow, :predictor_dvf], gn[1, :predictor_dvf]])
                                else
                                    predictor_dvf = missing
                                end
                                        
                                if :predictor_viralVerify in predictors
                                    predictor_viralVerify = find_consensus([gn[lastrow, :predictor_viralVerify], gn[1, :predictor_viralVerify]])
                                    shape_vv = find_consensus(gn[!, :shape_viralVerify])
                                else
                                    predictor_viralVerify = missing
                                    shape_vv = missing
                                end
                    
                                push!(anewdf, 
                                        (provirus_name = "$(gn[lastrow, :provirus_name])_____$(gn[1, :provirus_name])",
                                        contig_name = gn[1, :contig_name],
                                        provirus_start = gn[lastrow, :provirus_start],
                                        provirus_end = gn[lastrow, :contig_trimmed_end],
                                        length = (gn[1, :provirus_end] + (gn[lastrow, :contig_trimmed_end] - gn[lastrow, :provirus_start] +  1)),
                                        r_provirus_start = r_provirus_start,
                                        r_provirus_end = r_provirus_end,
                                        taxonomy = taxonomy,
                                        predictor_genomad = predictor_genomad,
                                        predictor_virSorter2 = predictor_virSorter2,
                                        predictor_vibrant = predictor_vibrant,
                                        predictor_dvf = predictor_dvf,
                                        predictor_viralVerify = predictor_viralVerify,
                                        shape_virSorter2 = shape_vs,
                                        topology_genomad = shape_g,
                                        shape_viralVerify = shape_vv,
                                        contig_shape = shape,
                                        contig_end_repeat_type = contig_end_repeat_type, 
                                        contig_shape_remarks = contig_shape_remarks, 
                                        contig_end_repeat_size = contig_end_repeat_size,
                                        contig_trimmed_end = find_consensus([gn[lastrow, :contig_trimmed_end], gn[1, :contig_trimmed_end]]),
                                        contig_full_end = find_consensus([gn[lastrow, :contig_full_end], gn[1, :contig_full_end]]),
                                        splitprovir = splitprovir
                                        ); cols=:union)


                                gn[1, :todel] = true
                                gn[lastrow, :todel] = true
                            end
                        end
                    end
                end
            end
        end
        #arow = nothing

        newdf = vcat(newdf, anewdf; cols = :union)

        # remove rows marked as todel (merged contigs spanning the cicular contig ends)
        for i in nrow(newdf):-1:1
            if ismissing(newdf[i, :todel]) == false && newdf[i, :todel] == true
                deleteat!(newdf, i)
            end
        end
    end 
    select!(newdf, Not(:todel))

    #endregion

    return newdf
end

#=
consensus_df_p = "/mnt/cephfs1/projects/DoViP_benchmarking/test_dataset/ALL_43_genomes_v1_mergeproph-true/03_I-01_checkV_Integrated/03_I-01_03_merged_integrated.tsv"
consensus_df = CSV.read(consensus_df_p, DataFrame; delim = '\t', header = 1)

dfj_p = "/mnt/cephfs1/projects/DoViP_benchmarking/test_dataset/ALL_43_genomes_v1_mergeproph-true/03_I-01_checkV_Integrated/03_I-01_02_postcheckV1_integrated_withmergedprovirIDs_cordf.tsv"
dfj = CSV.read(dfj_p, DataFrame; delim = '\t', header = 1) =#

function bring_circ_provirus_names!(dfj::DataFrame, consensus_df::DataFrame)
    lprovir_v = Vector{String}()
    rprovir_v = Vector{String}()

    for i in 1:nrow(consensus_df)
        if occursin("_____", consensus_df[i, :provirus_name])
            sp = split(consensus_df[i, :provirus_name], "_____")
            push!(lprovir_v, sp[1])
            push!(rprovir_v, sp[2])
        end
    end

    dfj[!, :provirus_name] = fill("", nrow(dfj))
    for i in 1:nrow(dfj)
        for j in eachindex(lprovir_v)
            if dfj[i, :provir] == lprovir_v[j] || dfj[i, :provir] == rprovir_v[j]
                dfj[i, :provirus_name] = "$(lprovir_v[j])_____$(rprovir_v[j])"
            end
        end
    
        if dfj[i, :provirus_name] == ""
            dfj[i, :provirus_name] = dfj[i, :provir] 
        end
    end

    return dfj
end

#dfj =  bring_circ_provirus_names!(dfj, consensus_df)

function initialize_predcov_dict(consensus_df, name_col::Symbol, start_col::Symbol, end_col::Symbol, length_col::Symbol)

    cov_dict = Dict{String, CovStruct}()

    for i in 1:nrow(consensus_df)
        key = consensus_df[i, name_col]

        # left end
        lDF_length = consensus_df[i, end_col] - consensus_df[i, start_col] + 1
        lDF = DataFrame(coordinates = fill(0, lDF_length), 
                        pred_cov = fill(0, lDF_length))

        lcoords_v = consensus_df[i, start_col]:consensus_df[i, end_col]             
        for j in 1:nrow(lDF)
            lDF[j, :coordinates] = lcoords_v[j]
        end

        # right end
        if "r_provirus_start" in names(consensus_df) && !ismissing(consensus_df[i, :r_provirus_start]) && !ismissing(consensus_df[i, :r_provirus_end])
            rDF_length = consensus_df[i, :r_provirus_end] - consensus_df[i, :r_provirus_start] + 1
            rDF = DataFrame(coordinates = fill(0, rDF_length), 
                            pred_cov = fill(0, rDF_length))
            rcoords_v = consensus_df[i, :r_provirus_start]:consensus_df[i, :r_provirus_end]             
            for j in 1:nrow(rDF)
                rDF[j, :coordinates] = rcoords_v[j]
            end
        else
            rDF = missing
        end

        cov_dict[key] = CovStruct(lDF, rDF, consensus_df[i, length_col], 0.0, 0.0)
    end

    return cov_dict
end

#= calculate predictor coverage for the unMixed contigs
predictors = [:predictor_genomad, :predictor_virSorter2, :predictor_vibrant, :predictor_dvf, :predictor_viralVerify]
int_df_p = "/mnt/cephfs1/projects/DoViP_benchmarking/NCBI_dataset/outputs_relaxed_DVF_th/ALL_45_genomes_v1/04_M_Detection/04_M_00_Mixed_Int_Viruses.tsv"
int_df = CSV.read(int_df_p, DataFrame; delim = '\t', header = 1)
nonint_df_p = "/mnt/cephfs1/projects/DoViP_benchmarking/NCBI_dataset/outputs_relaxed_DVF_th/ALL_45_genomes_v1/04_M_Detection/04_M_00_Mixed_NonInt_Viruses.tsv"
nonint_df = CSV.read(nonint_df_p, DataFrame; delim = '\t', header = 1)

consensus_df_p = "/mnt/cephfs1/projects/DoViP_benchmarking/NCBI_dataset/outputs_relaxed_DVF_th/ALL_45_genomes_v1/04_M_Detection/04_M_01_resolved_NonInt_Viruses.tsv"
consensus_df = CSV.read(consensus_df_p, DataFrame; delim = '\t', header = 1)

dfj = prep_joint_mixed_df(int_df, nonint_df, predictors)

name_col = :contig_name
start_col = :virus_start
end_col = :virus_end
length_col = :virus_length =#

function calculate_predcov!(dfj::DataFrame, consensus_df::DataFrame, predictors::Vector{Symbol}, name_col::Symbol, start_col::Symbol, init_start_col::Symbol, end_col::Symbol, init_end_col::Symbol, init_length_col::Symbol)
    # initialize initialize_predcov_dict
    predcov_dict = initialize_predcov_dict(consensus_df, name_col, init_start_col, init_end_col, init_length_col)

    # use dfj to calculate the "initial predictor" coverage for each merged pro-virus
    gdfj = groupby(dfj, name_col)

    for i in 1:length(gdfj)

        key = gdfj[i][1, name_col]

        if nrow(gdfj[i]) == 1
            predcov_dict[key].final_averagecov = 1.0
            predcov_dict[key].final_standad_deviation_cov = 0
        else # load predictor coverage at each position of the prophage

            if !occursin("_____", key)           # for linear prophages
                shift =  predcov_dict[key].lDF[1, :coordinates] 
                for a in 1:nrow(gdfj[i])

                    # find if it has more than one predictos (e.g. when beeign transfered from teh non-integrated branch)
                    pc = 0
                    for p in predictors 
                        if string(p) in names(gdfj[i]) && !ismissing(gdfj[i][a, p]) && gdfj[i][a, p] == "yes"
                            pc += 1
                        end
                    end

                    # left arm
                    for j in gdfj[i][a, start_col]:gdfj[i][a, end_col]
                        row = j - shift + 1
                        predcov_dict[key].lDF[row, :pred_cov] += pc                                           # add the number of predictors for one proviral fragment
                    end
                end 
            else                                # for wrapped prophages
                for a in 1:nrow(gdfj[i])

                    pc = 0
                    for p in predictors 
                        if string(p) in names(gdfj[i]) && !ismissing(gdfj[i][a, p]) && gdfj[i][a, p] == "yes"
                            pc += 1
                        end
                    end

                    # right arm
                    if ismissing(gdfj[i][a, :r_provir_start_cor]) && ismissing(gdfj[i][a, :r_provir_end_cor])
                        shift =  predcov_dict[key].rDF[1, :coordinates] 
                        for j in gdfj[i][a, :provir_start_cor]:gdfj[i][a, :provir_end_cor]
                            row = j - shift + 1
                            predcov_dict[key].rDF[row, :pred_cov] += pc 
                        end
                    else
                        # left arm
                        shift =  predcov_dict[key].lDF[1, :coordinates] 
                        for j in gdfj[i][a, start_col]:gdfj[i][a, end_col]
                            row = j - shift + 1
                            predcov_dict[key].lDF[row, :pred_cov] += pc
                        end
                    end
                end
            end

            # calculate coverage

            if ismissing(predcov_dict[key].rDF)
                all_covs = predcov_dict[key].lDF[!, :pred_cov]
            else
                all_covs = vcat(predcov_dict[key].lDF[!, :pred_cov], predcov_dict[key].rDF[!, :pred_cov])
            end

            favcov = round(mean(all_covs), digits = 2)
            fstd_cov = round(std(all_covs), digits = 2)

            predcov_dict[key].final_averagecov = favcov
            predcov_dict[key].final_standad_deviation_cov = fstd_cov
        end
    end

    consensus_df[!, :predictor_average_coverage] = fill(0.0, nrow(consensus_df))
    consensus_df[!, :predictor_stddev_coverage] = fill(0.0, nrow(consensus_df))

    for i in 1:nrow(consensus_df)
        key = consensus_df[i, name_col]
        consensus_df[i, :predictor_average_coverage] = predcov_dict[key].final_averagecov
        consensus_df[i, :predictor_stddev_coverage] = predcov_dict[key].final_standad_deviation_cov
    end

    if name_col == :contig_name
        for i in 1:nrow(consensus_df)
            if consensus_df[i, :mixed] == ""
                pc = 0

                for p in predictors
                    if (string(p) in names(consensus_df)) && (ismissing(consensus_df[i,p]) == false) && (consensus_df[i, p] == "yes")
                        pc += 1
                    end
                end

                consensus_df[i, :predictor_average_coverage] = pc 
                consensus_df[i, :predictor_stddev_coverage] = 0
            end
        end
    end

    return consensus_df
end

#consensus_df = calculate_predcov!(dfj, consensus_df, predictors, :contig_name, :virus_start, :virus_start, :virus_end, :virus_end, :virus_length)
#consensus_df = calculate_predcov!(dfj, consensus_df, predictors, :provirus_name, :provir_start_cor, :provirus_start, :provir_end_cor, :provirus_end, :length)


function merge_integrated!(inref::FnaP, checkV::ProjCheckVIntegrated, dfj::DataFrame, parentD::String) # I need to chck if the ouput files exist and if they are empty

    #make a new row for the right part of the prophages spanning the left end of the contig (overllaping on the left DTR)
    for i in 1:nrow(dfj)
        if ("splitprovir" in names(dfj) && ismissing(dfj[i, :splitprovir]) == false)

            actual_pred = Symbol()
            for pred in checkV.predictors 
                if String(pred) in names(dfj) && !ismissing(dfj[i, pred])
                    actual_pred = pred
                end
            end

            newrow = Dict(:seq_name => dfj[i, :seq_name],
                            :contig_name => dfj[i, :contig_name],
                            :splitprovir => dfj[i, :seq_name],
                            :provir_start_cor => dfj[i, :r_provir_start_cor],
                            :provir_end_cor => dfj[i, :r_provir_end_cor],
                            actual_pred => "yes")

            push!(dfj, newrow; cols = :union)
        end
    end

    # merge overlapping proviruses
    consensus_df = group_proviruses!(dfj, checkV.predictors, checkV.merge_circ_proph)  # this is modifying the input df!!! I should save it again to the HDD
    consensus_df[!, :virus_type_DoViP] = fill("Integrated", nrow(consensus_df))

    # bring provirus_name in dfj (it is missing for circular proviruses)
    dfj = bring_circ_provirus_names!(dfj, consensus_df)

    # save dfj, which now contains the name of the merged provirus as well (it was modified by group_proviruses)
    CSV.write("$(parentD)/$(checkV.postcheckv1_integrated_cor_withmergedprovirIDs_df.p)", dfj, delim = '\t', header = true)

    # calculate predictor coverage
    consensus_df = calculate_predcov!(dfj, consensus_df, checkV.predictors, :provirus_name, :provir_start_cor, :provirus_start, :provir_end_cor, :provirus_end, :length)

    # save df with merged proviruses
    CSV.write("$(parentD)/$(checkV.merged_integrated_DF.p)", consensus_df, delim = '\t', header = true)
    contigExtractRegions(FnaP("$(parentD)/$(inref.p)"), FnaP("$(parentD)/$(checkV.merged_integrated_fna.p)"), consensus_df[!, :contig_name], 
                            consensus_df[!, :provirus_name], consensus_df[!, :provirus_start], consensus_df[!, :provirus_end], 
                            consensus_df[!, :r_provirus_start], consensus_df[!, :r_provirus_end])
    
    return consensus_df

end

function post_checkV2_integrated!(proj::ProjCheckVIntegrated, indf::DataFrame, parentD::String)
    checkv2_summary_df = CSV.read("$(parentD)/$(proj.checkV2_out_integrated_summary_df.p)", DataFrame; delim='\t', header=1)
    rename!(checkv2_summary_df, :provirus => :provirus_checkV, :proviral_length => :proviral_length_checkV, :gene_count => :gene_count_checkV, :viral_genes => :viral_genes_checkV, :host_genes => :host_genes_checkV, 
                        :checkv_quality => :checkv_quality_checkV, :miuvig_quality => :miuvig_quality_checkV, :completeness => :completeness_checkV, :completeness_method => :completeness_method_checkV, 
                        :contamination => :contamination_checkV, :kmer_freq => :kmer_freq_checkV, :warnings => :warnings_checkV)

    checkV2_complete_df = CSV.read("$(parentD)/$(proj.checkV2_out_integrated_complete_df.p)", DataFrame; delim = '\t', header = 1)
    select!(checkV2_complete_df, :contig_id, :contig_length, :prediction_type, :confidence_reason, :repeat_length, :repeat_count)
    rename!(checkV2_complete_df, :prediction_type => :provirus_end_type_checkV, :confidence_reason => :completeness_confidence_reason_checkV, :repeat_length => :provirus_repeat_length_checkV, :repeat_count => :provirus_repeat_count_checkV)

    jdf = leftjoin!(indf, checkv2_summary_df, on = [:provirus_name => :contig_id, :length => :contig_length])
    jdf = leftjoin!(jdf, checkV2_complete_df, on = [:provirus_name => :contig_id, :length => :contig_length])

    CSV.write("$(parentD)/$(proj.postcheckV2_integrated_df_p.p)", jdf, delim = '\t', header = true)

    return jdf
end

#endregion

#region NONINTEGRATED - POST PREDICTION AND CHECKV
function merge_nonintegrated(inref::FnaP, proj::ProjCheckVNonintegrated, parentD::String)   # I need to chck if the ouput files exist and if they are empty
    
    # enter here post-predictor data, without any DoViP info regarding contig shape or ends
    # agregate and dereplicate the viral contigs (the join operation will dereplicate them)
    dfj = CSV.read("$(parentD)/$(proj.input_dfs_2_aggregate[1].p)", DataFrame; delim = '\t', header = 1)

    for i in 2:length(proj.input_dfs_2_aggregate)
        df = CSV.read("$(parentD)/$(proj.input_dfs_2_aggregate[i].p)", DataFrame; delim = '\t', header = 1)

        if ("length" in names(dfj)) && ("length" in names(df))
            dfj = outerjoin(dfj, df, on = [:contig_name, :length])
        else
            dfj = outerjoin(dfj, df, on = [:contig_name])
        end
    end

    if nrow(dfj) > 0
        # save a fasta file with dereplicated viral contigs
        out_fna_p = FnaP("$(parentD)/$(proj.output_aggregated_fna.p)")
        contigNameSel(FnaP("$(parentD)/$(inref.p)"), out_fna_p, dfj[!, :contig_name])
        
        # determine contig shape and join it to the out DF
        shape_ends_DF = detectContigShapeEnds(out_fna_p)
        dfj = leftjoin!(dfj, shape_ends_DF, on = [:contig_name => :contig_name])
        
        # set start and end positions for viruses (for the moment, they match the contig start and end), give the virus a name (if circular, it will be different than the contig name)
        dfj[!, :virus_type_DoViP] = fill("Non-integrated", nrow(dfj))
        dfj[!, :virus_start] = fill(1, nrow(dfj))
        dfj[!, :virus_end] = fill(1, nrow(dfj))
        dfj[!, :virus_name] = fill("", nrow(dfj))
    
        if "length" ∉ names(dfj)   # ∉  means NOT included in   # ∈ means included in
            dfj[!, :length] = fill(0, nrow(dfj))
        end

        for i in 1:nrow(dfj)
            # I'm not cutting the viral contigs here, because they should enter CheckV fully. DTRs don't always signal the contig ends
            dfj[i, :virus_end] = dfj[i, :contig_full_end]
            dfj[i, :length] = dfj[i, :contig_full_end]
            dfj[i, :virus_name] = dfj[i, :contig_name]
            #=if dfj[i, :contig_shape] == "circular"
                dfj[i, :virus_name] = dfj[i, :contig_name]*"___circ"
            else
                dfj[i, :virus_name] = dfj[i, :contig_name]
            end =#
        end
        
        rename!(dfj, :length => :virus_length)


        # save outDF
        CSV.write("$(parentD)/$(proj.output_aggregated_df.p)", dfj, delim = '\t', header = true)
    end

    return dfj
end

function postcheckV_nonintegrated!(proj::ProjCheckVNonintegrated, mergeddf::DataFrame, parentD::String)
    # bring CheckV data 
    checkV_summary_df = CSV.read("$(parentD)/$(proj.checkV_out_nonintegrated_summary_df.p)", DataFrame; delim = '\t', header = 1)
    checkV_complete_df = CSV.read("$(parentD)/$(proj.checkV_out_nonintegrated_complete_df.p)", DataFrame; delim = '\t', header = 1)
    select!(checkV_complete_df, :contig_id, :contig_length, :prediction_type, :confidence_reason, :repeat_length, :repeat_count)
    rename!(checkV_complete_df, :prediction_type => :contig_end_type_checkV, :confidence_reason => :completeness_confidence_reason_checkV, :repeat_length => :repeat_length_checkV, :repeat_count => :repeat_count_checkV)
    
    #region detect if some of the NONIntegrated viruse contigs are actually prophages (results from CheckV), pull them out of mergeddf and place them in a separate int_df
    contigs_intbranch = Vector{Union{String, InlineString}}()
    c = 0
    for i in nrow(checkV_summary_df):-1:1
        if checkV_summary_df[i, :provirus] == "Yes"
            push!(contigs_intbranch, checkV_summary_df[i, :contig_id])
            deleteat!(checkV_summary_df, i)
            c += 1
        end
    end
    
    if c > 0 && isfile("$(parentD)/$(proj.checkV_out_provir_fna.p)")
        # make a new mergeddf only for the int contigs, and delete the correposnding rows from the NonInt mergeddf
        int_df = DataFrame()

        for i in nrow(mergeddf):-1:1
            if (mergeddf[i, :contig_name] in contigs_intbranch) || (mergeddf[i, :virus_name] in contigs_intbranch)
                push!(int_df, mergeddf[i, :]; promote = true)
                deleteat!(mergeddf, i)
            end
        end

        #bring joined df to a format matching the postcheckV1 df from INTEGRATED Branch
        int_df[!, :seq_name] = fill("", nrow(int_df))
        for i in 1:nrow(int_df)
            int_df[i, :seq_name] = "$(int_df[i, :virus_name])_provirNonInt_$i"
        end

        int_df = rename!(int_df, :virus_length => :length, :virus_start => :provirus_start, :virus_end => :provirus_end)
        int_df_names = names(int_df)

        if "length_virSorter2" in int_df_names
            select!(int_df, Not(:length_virSorter2))
        end
        int_df = select!(int_df, Not(:virus_name, :virus_type_DoViP, :contig_shape, :contig_end_repeat_type, :contig_shape_remarks, 
                                        :contig_end_repeat_size, :contig_trimmed_end, :contig_full_end))

        CSV.write("$(parentD)/$(proj.postcheckV_integrated_df_p.p)", int_df, delim = '\t', header = true)
    end
    #endregion

    #region if there are any NONINTEGRATED viruses left in mergeddf, join mergeddgf with checkVdf and commpletedf (to bring data about the completeness detection method for complete contigs) !!!!!! PLACE A CHECK AFTER postcheckV_nonintegrated, if there are NONINTEGRATED, to go to final thresholding
    if nrow(checkV_summary_df) > 0
        rename!(checkV_summary_df, :provirus => :provirus_checkV, :proviral_length => :proviral_length_checkV, :gene_count => :gene_count_checkV, :viral_genes => :viral_genes_checkV, :host_genes => :host_genes_checkV, 
                            :checkv_quality => :checkv_quality_checkV, :miuvig_quality => :miuvig_quality_checkV, :completeness => :completeness_checkV, :completeness_method => :completeness_method_checkV, 
                            :contamination => :contamination_checkV, :kmer_freq => :kmer_freq_checkV, :warnings => :warnings_checkV)

        mergeddf = leftjoin!(mergeddf, checkV_summary_df, on = [:virus_name => :contig_id, :virus_length => :contig_length])
        mergeddf = leftjoin!(mergeddf, checkV_complete_df, on = [:virus_name => :contig_id, :virus_length => :contig_length])
        
        # prepare mergedDF tp export a fasta file with the non-integrated viruses, trimmed contigs
        mergeddf[!, :trimm] = Vector{Bool}(undef, nrow(mergeddf))
        mergeddf[!, :contig_end_for_saving] = Vector{Union{Missing, Int64}}(missing, nrow(mergeddf))
        mergeddf[!, :virus_name_trimmed] = fill("", nrow(mergeddf))
        
        for i in 1:nrow(mergeddf)
            # no need to check completeness, because if DTR or ITR is present in this column, completeness is always 100
            if occursin("DTR", mergeddf[i, :completeness_method_checkV]) #|| occursin("ITR", mergeddf[i, :completeness_method_checkV])
                mergeddf[i, :trimm] = true
                mergeddf[i, :contig_end_for_saving] = mergeddf[i, :contig_trimmed_end]
                mergeddf[i, :virus_name_trimmed] = "$(mergeddf[i, :virus_name])_trimmed"
            else
                mergeddf[i, :trimm] = false
                mergeddf[i, :contig_end_for_saving] = mergeddf[i, :contig_full_end]
                mergeddf[i, :virus_name_trimmed] = mergeddf[i, :virus_name] 
            end
        end

        # remove the column contig_end_for_saving
        # select!(mergeddf, Not(:contig_end_for_saving))
        # save NONINTEGRATED df after checkV
        CSV.write("$(parentD)/$(proj.postcheckV_nonintegrated_df_p.p)", mergeddf, delim='\t', header = true)
    end
    return mergeddf
    #endregion

end

#endregion

function merge_postCheckV_phaTYP_nonIntegrated!(phatypdf_p::TableP, indf_p::TableP, outp::TableP, parentD::String)
    indf = CSV.read("$(parentD)/$(indf_p.p)", DataFrame; delim = '\t', header =1)
    phatyp_df = CSV.read("$(parentD)/$(phatypdf_p.p)", DataFrame, delim=',', header=1)
    phatyp_df = rename!(phatyp_df, :Accession => :virus_name, :Pred => :prediction_PhaTYP, :Score => :score_PhaTYP)
    
    if "Length" ∈ names(phatyp_df) 
        select!(phatyp_df, Not(:Length)) 
    end

    jdf = leftjoin!(indf, phatyp_df, on = [:virus_name])
    
    CSV.write("$(parentD)/$(outp.p)", jdf, delim = '\t', header = true)

    return jdf
end

"""
    transfer_predictors(nonintMixed_tsv::DataFrame, intMixed_tsv::DataFrame, predictors_v::Vector{Symbol})
    This function brings the predictor data from the mixedInt Df to the nonInt Df (Based on the idea that, if a contig was predicted as NONINTEGRATED by some initial predictors and confirmed by CheckV, 
        then its merging with overlapping INTEGRATED viral sequences results in the same contig, full length)
"""
function transfer_predictors(nonintMixed_tsv::DataFrame, intMixed_tsv::DataFrame, predictors_v::Vector{Symbol})
    nonintMixed_df = deepcopy(nonintMixed_tsv)

    nonintMixed_df[!, :mixed] = fill("", nrow(nonintMixed_df))
    nonintMixed_df[!, :check_mixed] = fill("", nrow(nonintMixed_df))
    nonintMixed_df[!, :orig_NonInt_predictors] = fill("", nrow(nonintMixed_df))
    
    for i in 1:nrow(nonintMixed_df)

        # identify original NonInt predictors
        orig_pred_v = Vector{String}()
        for pred in predictors_v
            if ("$pred" in names(nonintMixed_df) && !ismissing(nonintMixed_df[i, pred]) && nonintMixed_df[i, pred] == "yes")
                push!(orig_pred_v, "$(pred)") 
            end
        end
        nonintMixed_df[i, :orig_NonInt_predictors] = find_consensus(orig_pred_v)

        # transfer predictors from int_df
        for j in 1:nrow(intMixed_tsv)
            if intMixed_tsv[j, :contig_name] == nonintMixed_df[i, :contig_name]
                nonintMixed_df[i, :mixed] = "yes"
                for pred in predictors_v
                    #println("i is $i, j is $j and pred is $pred")
                    if ("$pred" in names(nonintMixed_df) && ismissing(nonintMixed_df[i, pred]) || nonintMixed_df[i, pred] == "") && ("$pred" in names(intMixed_tsv) && !ismissing(intMixed_tsv[j, pred]) && intMixed_tsv[j, pred]  == "yes")
                        nonintMixed_df[i, pred] = "yes"
                    end                            
                end
                # do not merge int-viruses with non-int viruses if the non-int are larger by more than 5 times than the int 
                if (nonintMixed_df[i, :virus_length]/intMixed_tsv[j, :length]) > 4
                    nonintMixed_df[i, :check_mixed] = "yes"
                end
            end
        end
    end

    return nonintMixed_df
end

function export_nonint_fna(df::DataFrame, infna::FnaP, fullfna_p::FnaP, trimfna_p::FnaP, parentD::String)
    # export a fasta file with the non-integrated viruses, untrimmed contigs
    contigNameSel(FnaP("$(parentD)/$(infna.p)"), FnaP("$(parentD)/$(fullfna_p.p)"), df[!, :virus_name])

    # export a fasta file with non-integrated viruses, DTR right end trimmed
    trim = 0
    for i in 1:nrow(df)
        if df[i, :trimm] == true
        trim += 1 
        end
    end

    if trim > 0
        contigExtractRegions(FnaP("$(parentD)/$(infna.p)"), FnaP("$(parentD)/$(trimfna_p.p)"), df[!, :virus_name], df[!, :virus_name_trimmed], df[!, :virus_start], df[!, :contig_end_for_saving])
    end

    return nothing
end

function prep_joint_mixed_df(int_df::DataFrame, nonint_df::DataFrame, predictors::Vector{Symbol})

    cols_nonint = vcat([:contig_name, :virus_name, :virus_start, :virus_end], predictors)
    mixedNonInt_df_sel = select(nonint_df, cols_nonint)

    cols_int = vcat([:contig_name, :provirus_name, :provirus_start, :provirus_end], predictors)
    mixedInt_df_sel = select(int_df, cols_int)
    rename!(mixedInt_df_sel, :provirus_name => :virus_name, :provirus_start => :virus_start, :provirus_end => :virus_end)

    mixed_jdf = vcat(mixedNonInt_df_sel, mixedInt_df_sel, cols = :setequal)
    
    return mixed_jdf
end

#=
projb_p = "/mnt/cephfs1/projects/DoViP_benchmarking/NCBI_dataset/outputs_relaxed_DVF_th_predcov_mixed/ALL_45_genomes_v1/sproj.binary"
projb = deserialize(projb_p)
parentD = "/mnt/cephfs1/projects/DoViP_benchmarking/NCBI_dataset/outputs_relaxed_DVF_th_predcov_mixed"
proj = projb.detect_mixed_viruses =#

function detect_mixed_virs(proj::ProjDetectMixedViruses, parentD::String)
    indf_int_o = nothing
    indf_nonint_o = nothing

    outNonIntDf = nothing
    outIntDf = nothing


    if isfile("$(parentD)/$(proj.inDf_Int.p)")
        indf_int_o = CSV.read("$(parentD)/$(proj.inDf_Int.p)", DataFrame; delim = '\t', header=1)
        mixedInt_df = deepcopy(indf_int_o)
    end

    if isfile("$(parentD)/$(proj.inDf_NonInt.p)")
        indf_nonint_o = CSV.read("$(parentD)/$(proj.inDf_NonInt.p)", DataFrame; delim = '\t', header=1)
        mixedNonInt_df = deepcopy(indf_nonint_o)
    end

    if isfile("$(parentD)/$(proj.inDf_Int.p)") && isfile("$(parentD)/$(proj.inDf_NonInt.p)")
        mixed_contigs = intersect(mixedNonInt_df[!, :contig_name], mixedInt_df[!, :contig_name])

        if length(mixed_contigs) != 0
            for i in nrow(mixedInt_df):-1:1
                if !(mixedInt_df[i, :contig_name] in mixed_contigs)
                    deleteat!(mixedInt_df, i)
                end
            end
            CSV.write("$(parentD)/$(proj.outDf_mixed_Int_p.p)", mixedInt_df, delim = '\t', header = true)

            for i in nrow(mixedNonInt_df):-1:1
                if !(mixedNonInt_df[i, :contig_name] in mixed_contigs)
                    deleteat!(mixedNonInt_df, i)
                end
            end
            CSV.write("$(parentD)/$(proj.outDf_mixed_nonInt_p.p)", mixedNonInt_df, delim = '\t', header = true)
          
            # transfer predictor info from the INTEGRATED to the NONINTEGRATED mixed virus contigs
            outNonIntDf = transfer_predictors(indf_nonint_o, mixedInt_df, proj.predictors)
            # calculate coverage of mixed viruses and non-integrated viruses
            mixed_jdf = prep_joint_mixed_df(mixedInt_df, mixedNonInt_df, proj.predictors)
            outNonIntDf = calculate_predcov!(mixed_jdf, outNonIntDf, proj.predictors, :contig_name, :virus_start, :virus_start, :virus_end, :virus_end, :virus_length)

            CSV.write("$(parentD)/$(proj.outDf_resolved_nonInt_p.p)", outNonIntDf, delim = '\t', header = true)


            #clean the int_df by removing proviruses found on contigs also predicted as nonintegrated
            outIntDf = indf_int_o
            for i in nrow(outIntDf):-1:1
                if outIntDf[i, :contig_name] in mixed_contigs
                    deleteat!(outIntDf, i)
                end
            end

            if nrow(outIntDf) > 0
                CSV.write("$(parentD)/$(proj.outDf_resolved_Int_p.p)", outIntDf, delim = '\t', header = true)
            end

        end
    end

    # export nonintegrated fna
    if !isnothing(outNonIntDf)
        export_nonint_fna(outNonIntDf, proj.inFna_NonInt, proj.outFna_nonint_p, proj.outFna_nonint_DTR_trimmed_p, parentD)
    else
        if !isnothing(indf_nonint_o)
            # add initial predictor coverage
            indf_nonint_o[!, :predictor_average_coverage] = fill(0, nrow(indf_nonint_o)) 
            indf_nonint_o[!, :predictor_stddev_coverage] = fill(0, nrow(indf_nonint_o)) 

            for i in 1:nrow(indf_nonint_o)
                pc = 0
                for p in proj.predictors
                    if (string(p) in names(indf_nonint_o)) && (ismissing(indf_nonint_o[i,p]) == false) && (indf_nonint_o[i, p] == "yes")
                        pc += 1
                    end
                end
                indf_nonint_o[i, :predictor_average_coverage] = pc 
            end
            
            # save .tsv
            CSV.write("$(parentD)/$(proj.outDf_resolved_nonInt_p.p)", indf_nonint_o, delim = '\t', header = true)
            export_nonint_fna(indf_nonint_o, proj.inFna_NonInt, proj.outFna_nonint_p, proj.outFna_nonint_DTR_trimmed_p, parentD)
        end
    end

    # export integrated fna
    if !isnothing(outIntDf) && nrow(outIntDf) > 0
        contigNameSel(FnaP("$(parentD)/$(proj.inFna_Int.p)"), FnaP("$(parentD)/$(proj.outFna_Int_p.p)"), outIntDf[!, :provirus_name]) 
    else
        if !isnothing(indf_int_o)
            CSV.write("$(parentD)/$(proj.outDf_resolved_Int_p.p)", indf_int_o, delim = '\t', header = true)
            contigNameSel(FnaP("$(parentD)/$(proj.inFna_Int.p)"), FnaP("$(parentD)/$(proj.outFna_Int_p.p)"), indf_int_o[!, :provirus_name]) 
        end
    end
    return nothing
end

"""
    post_genomadTax(genomad_out_df_p::TableP, prev_df_p::TableP, col::Symbol, parentD::String)
    Bring genomad taxonomy in the NonIntegrated and INTEGRATED viral sequences tables
"""
function post_genomadTax(proj::ProjGenomadTax, col::Symbol, parentD::String)
    genomad_out_df = CSV.read("$(parentD)/$(proj.genomadtax_out_table_p.p)", DataFrame; delim = '\t', header =1)
    select!(genomad_out_df, [:seq_name, :taxid, :lineage])
    rename!(genomad_out_df, :taxid => :genomadTax_taxid, :lineage => :genomadTax_lineage)
    
    prev_df = CSV.read("$(parentD)/$(proj.previous_df.p)", DataFrame; delim = '\t', header =1)

    leftjoin!(prev_df, genomad_out_df, on = [col => :seq_name])

    CSV.write("$(parentD)/$(proj.postgenomadtax_df.p)", prev_df; delim = '\t', header = true )

    return nothing
end

function apply_thresholds!(proj::FinalThresholding, name_col::Symbol, parentD::String, orderfun::Function, sample_name::String, sample_set::String)

    if isfile("$(parentD)/$(proj.inTsv.p)")
        tdf = CSV.read("$(parentD)/$(proj.inTsv.p)", DataFrame; delim = '\t', header = 1)
        tdf[!, :predictors_total] = Vector{Union{Missing, Int64}}(missing, nrow(tdf))
        
        println("\nInital predictors for this branch are $(proj.predictors)")

        for i in nrow(tdf):-1:1
            # calculate the number of predictors for storage in DF
            pc = 0
            for p in proj.predictors
                if (string(p) in names(tdf)) && (ismissing(tdf[i,p]) == false) && (tdf[i, p] == "yes")
                    pc += 1
                end
            end
            tdf[i, :predictors_total] = pc 
            
            #=
            if name_col == :provirus_name
                th = tdf[i, :predictor_average_coverage]
            elseif name_col == :virus_name
                th = pc 
            end =#
            
            th = tdf[i, :predictor_average_coverage]

            if tdf[i, :completeness_checkV] == "NA" && th < proj.th_num_predictors_CheckV_NA
                deleteat!(tdf, i)
            elseif tdf[i, :completeness_checkV] != "NA"
                if typeof(tdf[i, :completeness_checkV]) == Float64
                    completeness = tdf[i, :completeness_checkV]
                else
                    completeness = parse(Float64, tdf[i, :completeness_checkV])
                end
    
                host_genes = ((tdf[i, :host_genes_checkV] * 100) / tdf[i, :gene_count_checkV])
    
                if th <= 1.1  && host_genes >= 70
                    deleteat!(tdf, i)
                elseif occursin("AAI-based (high-confidence)", tdf[i, :completeness_method_checkV])
                    if (completeness < 90) && (th < proj.th_num_predictors_CheckV_AAIHighConf || completeness <= proj.th_completeness_CheckV_AAIHighConf)
                        deleteat!(tdf, i)
                    end
                elseif occursin("AAI-based (medium-confidence)", tdf[i, :completeness_method_checkV])
                    if th < proj.th_num_predictors_CheckV_AAIMediumConf || completeness <= proj.th_completeness_CheckV_AAIMediumConf
                        deleteat!(tdf, i)
                    end
                elseif occursin("HMM-based (lower-bound)", tdf[i, :completeness_method_checkV])
                    if th < proj.th_num_predictors_CheckV_HMM || completeness <= proj.th_completeness_CheckV_HMM
                        deleteat!(tdf, i)
                    end
                elseif occursin("DTR", tdf[i, :completeness_method_checkV]) || occursin("ITR", tdf[i, :completeness_method_checkV]) 
                        if occursin("AAI-based", tdf[i, :completeness_confidence_reason_checkV]) 
                            # I will not check completeness, because if DTR or ITR is pressent, completeness is always 100
                            if ismissing(proj.th_num_predictors_CheckV_DTR_ITR_AAI) == false && th < proj.th_num_predictors_CheckV_DTR_ITR_AAI
                                deleteat!(tdf, i)
                            end
                        elseif occursin("HMM-based", tdf[i, :completeness_confidence_reason_checkV])
                            if ismissing(proj.th_num_predictors_CheckV_DTR_ITR_HMM) == false && th < proj.th_num_predictors_CheckV_DTR_ITR_HMM
                                deleteat!(tdf, i)
                            end
                        end
                else
                    if th < 2.5
                        deleteat!(tdf, i)
                    end   
                end
            end
            
        end

        if nrow(tdf) > 0 
            if !ismissing(proj.inFna_trimmed_DTR) && isfile("$(parentD)/$(proj.inFna_trimmed_DTR.p)")
                contigNameSel(FnaP("$(parentD)/$(proj.inFna_trimmed_DTR.p)"), FnaP("$(parentD)/$(proj.outFna_trimmed_DTR.p)"), tdf[!, :virus_name_trimmed])
            end

            tdf = orderfun(tdf, sample_name, sample_set)

            CSV.write("$(parentD)/$(proj.outTsv.p)", tdf, delim = '\t', header = true)
            contigNameSel(FnaP("$(parentD)/$(proj.inFna.p)"), FnaP("$(parentD)/$(proj.outFnaP.p)"), tdf[!, name_col])
        end
    end

    return nothing
end

