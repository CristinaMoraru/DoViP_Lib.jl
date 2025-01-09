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
    df =  CSV.read("$(parentD)/$(proj.modif_output_dvf_f.p)", DataFrame; delim = '\t', header = 1)

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

    # mark prophages which wrap themselves around a circular contig (with DTR ends)
    if "shape_virSorter2" in names(dfj) 
        for i in 1:nrow(dfj)
            if !ismissing(dfj[i, :shape_virSorter2]) && dfj[i, :shape_virSorter2] == "circular" && dfj[i, :contig_shape] == "circular"
                if dfj[i, :provirus_end] > dfj[i, :contig_trimmed_end]  
                    dfj[i, :r_provir_start] = 1
                    dfj[i, :r_provir_end] = dfj[i, :provirus_end] - dfj[i, :contig_trimmed_end] 
                    dfj[i, :provirus_end] = dfj[i, :contig_trimmed_end] 
                    dfj[i, :splitprovir] = dfj[i, :seq_name]
                    dfj[i, :length] = dfj[i, :r_provir_end] + (dfj[i, :provirus_end] - dfj[i, :provirus_start] + 1)
                end
            elseif !ismissing(dfj[i, :shape_virSorter2]) && dfj[i, :shape_virSorter2] == "circular" && dfj[i, :contig_shape] == "linear"
                if dfj[i, :provirus_end] > dfj[i, :contig_full_end]  
                    dfj[i, :provirus_end] = dfj[i, :contig_full_end]  
                end
            end
        end
    end

    # find corresponding contig names for the split proviruses (as detected by VirSorter2) and save them as a Set
    splitprovir_set = collect(skipmissing(unique(dfj[!, :splitprovir])))
    
    for i in eachindex(splitprovir_set)
        splitprovir_set[i] = replace(splitprovir_set[i], Regex("\\|\\|\\d+_partial") => "")
    end
    splitprovir_set = Set(unique(splitprovir_set))

    for i in 1:nrow(dfj)
        dfj[i, :provirID1] = "provirID1num$(i)"

        #= trim prophages predicted by the other predictors found on circular contigs that extend in the right direct terminal repeat
        end, contigs for which virSorter2 found proviruses wrapping around the ends of contigs =#
        # check if the value in column seq_name is in the splitprovir_vec 
        if dfj[i, :contig_name] in splitprovir_set
            if dfj[i, :provirus_end] > dfj[i, :contig_trimmed_end]
                dfj[i, :provirus_end] = dfj[i, :contig_trimmed_end]
                dfj[i, :length] = dfj[i, :provirus_end] - dfj[i, :provirus_start] + 1
            end
        end

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
    cons = unique(vc) |> skipmissing |> collect
    
    if length(cons) == 1
        return cons[1]
    elseif length(cons) > 1
        return "conflicting"
    elseif length(cons) == 0
        return missing
    end
end

function group_proviruses!(df::DataFrame, predictors::Vector{Symbol})
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
    newdf[!, :predictor_genomad] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :predictor_virSorter2] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :predictor_vibrant] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :predictor_dvf] = Vector{Union{Missing, String}}(missing, nrow(newdf))
    newdf[!, :predictor_viralVerify] = Vector{Union{Missing, String}}(missing, nrow(newdf))
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
    anewdf = DataFrame()
    #arow = 0
    gnewdf = groupby(newdf, :contig_name)
    for gn in gnewdf     
        shape = find_consensus(gn[!, :contig_shape]) 
        if :predictor_virSorter2 in predictors
            shape_vs = find_consensus(gn[!, :shape_virSorter2])
        else
            shape_vs = missing
        end

        if "splitprovir" in names(gn)
            split_vs = find_consensus(gn[!, :splitprovir])
        else
            split_vs = missing
        end

        # I will merge only pro-viruses found on circular contigs (with DTRs) and found by VirSorter2 to wrapp around the contig ends
        if (!ismissing(shape) && shape == "circular") && (:predictor_virSorter2 in predictors && !ismissing(shape_vs) && shape_vs == "circular" && !ismissing(split_vs))
            sort!(gn, :provirus_start)
            lastrow = nrow(gn)
            if (gn[lastrow, :contig_trimmed_end] - gn[lastrow, :provirus_end]) <= 1000 
                spaceL = gn[lastrow, :contig_trimmed_end] - gn[lastrow, :provirus_end]                       
                if gn[1, :provirus_start] <= 1000
                    spaceR = gn[1, :provirus_start]
                    if (spaceR + spaceL) <= 1000

                        contig_end_repeat_type = find_consensus(gn[!, :contig_end_repeat_type])
                        contig_shape_remarks = find_consensus(gn[!, :contig_shape_remarks])
                        contig_end_repeat_size = find_consensus(gn[!, :contig_end_repeat_size])

                        colnames = names(gn)
                        if "r_provirus_start" in colnames
                            #anewdf[arow, :r_provirus_start] = 1
                            r_provirus_start = 1
                        else
                            r_provirus_start = missing
                        end

                        if "r_provirus_end" in colnames
                            #anewdf[arow, :r_provirus_end] = gn[1, :provirus_end]
                            r_provirus_end = gn[1, :provirus_end]
                        else
                            r_provirus_end = missing
                        end
      
                        if :predictor_genomad in predictors
                            #anewdf[arow, :taxonomy] = find_consensus([gn[lastrow, :taxonomy], gn[1, :taxonomy]])
                            taxonomy = find_consensus([gn[lastrow, :taxonomy], gn[1, :taxonomy]])
                            #anewdf[arow, :predictor_genomad] = find_consensus([gn[lastrow, :predictor_genomad], gn[1, :predictor_genomad]])
                            predictor_genomad = find_consensus([gn[lastrow, :predictor_genomad], gn[1, :predictor_genomad]])
                            shape_g = find_consensus(gn[!, :topology_genomad])
                        else
                            taxonomy = missing
                            predictor_genomad = missing
                            shape_g = missing
                        end

                        if :predictor_virSorter2 in predictors
                            #anewdf[arow, :predictor_virSorter2] = find_consensus([gn[lastrow, :predictor_virSorter2], gn[1, :predictor_virSorter2]])
                            predictor_virSorter2 = find_consensus([gn[lastrow, :predictor_virSorter2], gn[1, :predictor_virSorter2]])
                            #anewdf[arow, :shape_virSorter2] = shape_vs
                        else
                            predictor_virSorter2 = missing
                        end
                             
                        if :predictor_vibrant in predictors
                            #anewdf[arow, :predictor_vibrant] = find_consensus([gn[lastrow, :predictor_vibrant], gn[1, :predictor_vibrant]])
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
                                ); cols=:union)


                        gn[1, :todel] = true
                        gn[lastrow, :todel] = true
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

    select!(newdf, Not(:todel))

    #endregion

    return newdf
end

function merge_integrated!(inref::FnaP, checkV::ProjCheckVIntegrated, dfj::DataFrame, parentD::String) # I need to chck if the ouput files exist and if they are empty

    #make a new row for the right part of the circular prophages predicted by VS2
    for i in 1:nrow(dfj)
        if ("splitprovir" in names(dfj) && ismissing(dfj[i, :splitprovir]) == false)
            newrow = Dict(:seq_name => dfj[i, :seq_name],
                            :contig_name => dfj[i, :contig_name],
                            :splitprovir => dfj[i, :seq_name],
                            :provir_start_cor => dfj[i, :r_provir_start_cor],
                            :provir_end_cor => dfj[i, :r_provir_end_cor],
                            :predictor_virSorter2 => "yes")

            push!(dfj, newrow; cols = :union)
        end
    end

    # merge overlapping proviruses
    consensus_df = group_proviruses!(dfj, checkV.predictors)  # this is modifying the input df!!! I should save it again to the HDD
    consensus_df[!, :virus_type_DoViP] = fill("Integrated", nrow(consensus_df))

    # save dfj, which now contains the name of the merged provirus as well (it was modified by group_proviruses)
    CSV.write("$(parentD)/$(checkV.postcheckv1_integrated_cor_withmergedprovirIDs_df.p)", dfj, delim = '\t', header = true)

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
        
        # export a fasta file with the non-integrated viruses, untrimmed contigs
        contigNameSel(FnaP("$(parentD)/$(proj.output_aggregated_fna.p)"), FnaP("$(parentD)/$(proj.postcheckV_nonintegrated_fna.p)"), mergeddf[!, :virus_name])

        # export a fasta file with the non-integrated viruses, trimmed contigs
        mergeddf[!, :trimm] = Vector{Bool}(undef, nrow(mergeddf))
        mergeddf[!, :contig_end_for_saving] = Vector{Union{Missing, Int64}}(missing, nrow(mergeddf))
        mergeddf[!, :virus_name_trimmed] = fill("", nrow(mergeddf))
        
        trim = 0
        for i in 1:nrow(mergeddf)
            # no need to check completeness, because if DTR or ITR is present in this column, completeness is always 100
            if occursin("DTR", mergeddf[i, :completeness_method_checkV]) #|| occursin("ITR", mergeddf[i, :completeness_method_checkV])
                mergeddf[i, :trimm] = true
                mergeddf[i, :contig_end_for_saving] = mergeddf[i, :contig_trimmed_end]
                mergeddf[i, :virus_name_trimmed] = "$(mergeddf[i, :virus_name])_trimmed"
                trim += 1
            else
                mergeddf[i, :trimm] = false
                mergeddf[i, :contig_end_for_saving] = mergeddf[i, :contig_full_end]
                mergeddf[i, :virus_name_trimmed] = mergeddf[i, :virus_name] 
            end
        end

        if trim > 0
            contigExtractRegions(FnaP("$(parentD)/$(proj.output_aggregated_fna.p)"), FnaP("$(parentD)/$(proj.postcheckV_nonintegrated_fna_trimmed_DTR.p)"), mergeddf[!, :virus_name], mergeddf[!, :virus_name_trimmed], mergeddf[!, :virus_start], mergeddf[!, :contig_end_for_saving])
        end

        # remove the column contig_end_for_saving
        select!(mergeddf, Not(:contig_end_for_saving))
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

function detect_mixed_virs(proj::ProjDetectMixedViruses, parentD::String)
    indf_int = CSV.read("$(parentD)/$(proj.inDf_Int.p)", DataFrame; delim = '\t', header=1)
    indf_nonint = CSV.read("$(parentD)/$(proj.inDf_NonInt.p)", DataFrame; delim = '\t', header=1)

    mixed_contigs = intersect(indf_nonint[!, :contig_name], indf_int[!, :contig_name])

    if length(mixed_contigs) != 0
        for i in nrow(indf_int):-1:1
            if !(indf_int[i, :contig_name] in mixed_contigs)
                deleteat!(indf_int, i)
            end
        end

        for i in nrow(indf_nonint):-1:1
            if !(indf_nonint[i, :contig_name] in mixed_contigs)
                deleteat!(indf_nonint, i)
            end
        end

    CSV.write("$(parentD)/$(proj.outDf_Int.p)", indf_int, delim = '\t', header = true)
    CSV.write("$(parentD)/$(proj.outDf_NonInt.p)", indf_nonint, delim = '\t', header = true)
    end

    return nothing
end




function apply_thresholds!(proj::FinalThresholding, name_col::Symbol, parentD::String, orderfun::Function, sample_name::String, sample_set::String)

    if isfile("$(parentD)/$(proj.inTsv.p)")
        tdf = CSV.read("$(parentD)/$(proj.inTsv.p)", DataFrame; delim = '\t', header = 1)
        tdf[!, :predictors_total] = Vector{Union{Missing, Int64}}(missing, nrow(tdf))
        
        println("
Inital predictors for this branch are $(proj.predictors)")

        for i in nrow(tdf):-1:1
            pc = 0

            for p in proj.predictors
                if (string(p) in names(tdf)) && (ismissing(tdf[i,p]) == false) && (tdf[i, p] == "yes")
                    pc += 1
                end
            end

            tdf[i, :predictors_total] = pc
             
            if tdf[i, :completeness_checkV] == "NA" && pc < proj.th_num_predictors_CheckV_NA
                deleteat!(tdf, i)
            elseif tdf[i, :completeness_checkV] != "NA"
                if typeof(tdf[i, :completeness_checkV]) == Float64
                    completeness = tdf[i, :completeness_checkV]
                else
                    completeness = parse(Float64, tdf[i, :completeness_checkV])
                end
    
                host_genes = ((tdf[i, :host_genes_checkV] * 100) / tdf[i, :gene_count_checkV])
    
                if pc == 1 && host_genes >= 60
                    deleteat!(tdf, i)
                elseif occursin("AAI-based (high-confidence)", tdf[i, :completeness_method_checkV])
                    if (completeness < 90) && (pc < proj.th_num_predictors_CheckV_AAIHighConf || completeness <= proj.th_completeness_CheckV_AAIHighConf)
                        deleteat!(tdf, i)
                    end
                elseif occursin("AAI-based (medium-confidence)", tdf[i, :completeness_method_checkV])
                    if pc < proj.th_num_predictors_CheckV_AAIMediumConf || completeness <= proj.th_completeness_CheckV_AAIMediumConf
                        deleteat!(tdf, i)
                    end
                elseif occursin("HMM-based (lower-bound)", tdf[i, :completeness_method_checkV])
                    if pc < proj.th_num_predictors_CheckV_HMM || completeness <= proj.th_completeness_CheckV_HMM
                        deleteat!(tdf, i)
                    end
                elseif occursin("DTR", tdf[i, :completeness_method_checkV]) || occursin("ITR", tdf[i, :completeness_method_checkV]) 
                        if occursin("AAI-based", tdf[i, :completeness_confidence_reason_checkV]) 
                            # I will not check completeness, because if DTR or ITR is pressent, completeness is always 100
                            if ismissing(proj.th_num_predictors_CheckV_DTR_ITR_AAI) == false && pc < proj.th_num_predictors_CheckV_DTR_ITR_AAI
                                deleteat!(tdf, i)
                            end
                        elseif occursin("HMM-based", tdf[i, :completeness_confidence_reason_checkV])
                            if ismissing(proj.th_num_predictors_CheckV_DTR_ITR_HMM) == false && pc < proj.th_num_predictors_CheckV_DTR_ITR_HMM
                                deleteat!(tdf, i)
                            end
                        end
                else
                    if pc < 3
                        deleteat!(tdf, i)
                    end   
                end
            end
            
        end

        if !ismissing(proj.inFna_trimmed_DTR) && isfile("$(parentD)/$(proj.inFna_trimmed_DTR.p)")
            contigNameSel(FnaP("$(parentD)/$(proj.inFna_trimmed_DTR.p)"), FnaP("$(parentD)/$(proj.outFna_trimmed_DTR.p)"), tdf[!, :virus_name_trimmed])
        end

        tdf = orderfun(tdf, sample_name, sample_set)

        CSV.write("$(parentD)/$(proj.outTsv.p)", tdf, delim = '\t', header = true)
        contigNameSel(FnaP("$(parentD)/$(proj.inFna.p)"), FnaP("$(parentD)/$(proj.outFnaP.p)"), tdf[!, name_col])
    end

    return nothing
end
