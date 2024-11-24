export ProjSViP_fun

function set_contig_length(inref::FnaP, pd::String, sampleName::String, min_contig_length::Int64) 
    sp = "ALL_00"
    dir = "$(sampleName)/$(sp)_min_length_$(min_contig_length)"
    my_mkpath(["$(pd)/$(dir)"])    
    
    outref = "$(dir)/$(sampleName).fna" 

    proj = Proj_sel_contig_length(dir, min_contig_length, inref, FnaP(outref))

    return proj
end

set_shape_contigs() = missing

function set_predictors(dosteps::Dict{String, WorkflowStatus}) 
    predictors_all = Vector{Symbol}()
    predictors_Int = Vector{Symbol}()
    
    for (k, v) in dosteps
        if k == "genomad" #, "DVF", "virSorter2", "vibrant"]
            if v.signal == "use" || v.signal == "do" || v.signal == "use_external"
                push!(predictors_all, :predictor_genomad)
                push!(predictors_Int, :predictor_genomad)
            end
        elseif k == "DVF"
            if v.signal == "use" || v.signal == "do" || v.signal == "use_external"
                push!(predictors_all, :predictor_dvf)
            end
        elseif k == "virSorter2"
            if v.signal == "use" || v.signal == "do" || v.signal == "use_external"
                push!(predictors_all, :predictor_virSorter2)
                push!(predictors_Int, :predictor_virSorter2)
            end
        elseif k == "vibrant"
            if v.signal == "use" || v.signal == "do" || v.signal == "use_external"
                push!(predictors_all, :predictor_vibrant)
                push!(predictors_Int, :predictor_vibrant)
            end
        elseif  k == "viralVerify"
            if v.signal == "use" || v.signal == "do" || v.signal == "use_external"
                push!(predictors_all, :predictor_viralVerify)
            end
        end
    end
    return predictors_all, predictors_Int
end

function set_ProjGenomad(pd::String, sampleName::String, inref::FnaP, args::Vector{String}, signal::String)
    sp = "ALL-01a"  # step prefix
    genomad_D = "$(sampleName)/$(sp)_genomad"
    my_mkpath(["$(pd)/$(genomad_D)"])

    if signal == "use_external"
        ext_res = true
        ext_res_D = extract_inPaths(args, "genomad_res")
        out_genomad_D = ""
    else
        ext_res = false
        ext_res_D = nothing
        out_genomad_D = "$(genomad_D)/$(sp)_01_genomad_out"
        genomadLogs_predict_D = "$out_genomad_D/logs"
        my_mkpath(["$(pd)/$(genomadLogs_predict_D)"])

        genomad_env = extract_args(args, "genomad_env")
        genomadDB_p = extract_inPaths(args, "genomadDB_p")
        min_score = extract_args(args, "genomad_min_score", Float64, 0.7, 0.0, 1.0)
        genomad_sbatch_time = extract_args(args, "genomad_sbatch_time")
        genomad_cpus_per_task = extract_args(args, "genomad_cpus_per_task", Int64, 20, 1, 200)
        genomad_sbatch_mem = extract_args(args, "genomad_sbatch_mem")
    end

    genomad_out_tsv = "$out_genomad_D/$(sampleName)_summary/$(sampleName)_virus_summary.tsv" |> TableP
    genomad_out_fna = "$out_genomad_D/$(sampleName)_summary/$(sampleName)_virus.fna" |> FnaP

    postgenomad_genomad_nonintegrated_df = "$genomad_D/$(sp)_02_postgenomad_nonintegrated_df.tsv" |> TableP
    postgenomad_genomad_integrated_df = "$genomad_D/$(sp)_02_postgenomad_integrated_df.tsv" |> TableP
    #postgenomad_genomad_integrated_fnaf = "$genomad_D/$(sp)_02_postgenomad_integrated_fna.fna" |> FnaP
    todelete = ["$out_genomad_D/$(sampleName)_aggregated_classification", 
                "$out_genomad_D/$(sampleName)_annotate", 
                "$out_genomad_D/$(sampleName)_find_proviruses",
                "$out_genomad_D/$(sampleName)_marker_classification", 
                "$out_genomad_D/$(sampleName)_nn_classification",
                "$out_genomad_D/$(sampleName)_aggregated_classification.log", 
                "$out_genomad_D/$(sampleName)_annotate.log", 
                "$out_genomad_D/$(sampleName)_find_proviruses.log",
                "$out_genomad_D/$(sampleName)_marker_classification.log",
                "$out_genomad_D/$(sampleName)_nn_classification.log",]

    if signal == "use_external"
        genomad = nothing
    else
        genomad = WrapCmd(; cmd = RunGenomadCmd("end-to-end", inref, out_genomad_D, genomadDB_p, min_score), 
                            log_p = "$(genomadLogs_predict_D)/genomad.log", err_p = "$(genomadLogs_predict_D)/genomad.err", 
                            exit_p = "$(genomadLogs_predict_D)/genomad.exit", env = genomad_env,
                            sbatch_maxtime = genomad_sbatch_time, sbatch_cpuspertask = genomad_cpus_per_task, sbatch_mem = genomad_sbatch_mem)
    end            

    proj = ProjGenomad(genomad_D, ext_res, ext_res_D, genomad, genomad_out_tsv, genomad_out_fna, postgenomad_genomad_nonintegrated_df, postgenomad_genomad_integrated_df, todelete)

    return proj
end

function set_ProjDVF(pd::String, sampleName::String, inref::FnaP, args::Vector{String}, min_contig_length::Int64)
    sp = "ALL-01b"
    DVF_D = "$sampleName/$(sp)_DVF"
    outDVF = "$(DVF_D)/$(sp)_01_DVF_out"
    DVFLogs_D = "$outDVF/logs"
    my_mkpath(["$(pd)/$(DVFLogs_D)"])

    DVF_env = extract_args(args, "DVF_env")
    DVF_p = extract_args(args, "DVF_p")

    #resources
    DVF_sbatch_time = extract_args(args, "DVF_sbatch_time")
    DVF_cpus_per_task = extract_args(args, "DVF_cpus_per_task", Int64, 15, 1, 200)
    DVF_sbatch_mem = extract_args(args, "DVF_sbatch_mem")

    #contig len
    #DVF_minContigLen = extract_args(args, "DVF_minContigLen", Int64, 100, 1, 2099000)
    DVF_maxContigLen = extract_args(args, "DVF_maxContigLen", Int64, 2099000, 1, 2099000)

    #scores
    DVF_scoreTh = extract_args(args, "DVF_scoreTh", Float64, 0.7, 0.0, 1.0)
    DVF_pThreshold = extract_args(args, "DVF_pThreshold", Float64, 0.01, 0.0, 1.0)

    # input and output files
    cleaned_dvf_input_f = "$outDVF/$(sampleName)_sizeForDVF.fna" |> FnaP
    output_dvf_f = "$outDVF/$(basename(cleaned_dvf_input_f.p))_gt$(min_contig_length)bp_dvfpred.txt" |> TableP
    modif_output_dvf_f = "$outDVF/$(getFileName(output_dvf_f.p))_modif.txt" |> TableP
    postdfv_nonintegrated_fna = "$DVF_D/$(sp)_postdfv_nonintegrated.fna" |> FnaP
    postdvf_nonintegrated_df = "$DVF_D/$(sp)_postdfvf_nonintegrated_df.tsv" |> TableP

    todelete = [cleaned_dvf_input_f.p]


    runDVF = WrapCmd(; cmd = RunDVFCmd(DVF_p, cleaned_dvf_input_f, outDVF, min_contig_length, DVF_cpus_per_task), 
                    log_p = "$(DVFLogs_D)/DVF.log", err_p = "$(DVFLogs_D)/DVF.err", exit_p = "$(DVFLogs_D)/DVF.exit", env = DVF_env,
                    sbatch_maxtime = DVF_sbatch_time, sbatch_cpuspertask = DVF_cpus_per_task, sbatch_mem = DVF_sbatch_mem)
    
    proj = ProjDVF(DVF_D, inref, DVF_maxContigLen, runDVF, DVF_scoreTh, DVF_pThreshold, output_dvf_f, modif_output_dvf_f, 
                    postdfv_nonintegrated_fna, postdvf_nonintegrated_df, todelete)

    return proj
end

function set_ProjVirSorter2(pd::String, sampleName::String, inref::FnaP, args::Vector{String}, signal::String, min_contig_length::Int64)
    sp = "ALL-01c"
    virSorter2_D = "$sampleName/$(sp)_virSorter2"
    my_mkpath(["$(pd)/$(virSorter2_D)"])

    if signal == "use_external"
        ext_res = true
        ext_res_D = extract_inPaths(args, "virSorter2_res")
        out_vs_D = ""
    else
        ext_res = false
        ext_res_D = nothing
        out_vs_D = "$(virSorter2_D)/$(sp)_01_virSorter2_out"
        virSorter2Logs_D = "$virSorter2_D/logs"
        my_mkpath(["$(pd)/$(virSorter2Logs_D)"])

        virSorter2_env = extract_args(args, "virSorter2_env")
        virSorter2DB_p = extract_inPaths(args, "virSorter2DB_p")
        virSorter2_high_confidence_only = extract_args(args, "virSorter2_high_confidence_only", Bool, "false")
        min_score = extract_args(args, "virSorter2_min_score", Float64, 0.5, 0.0, 1.0)
        #min_length = extract_args(args, "virSorter2_min_length", Int64, 1000, 1, 2099000)
        virSorter2_sbatch_time = extract_args(args, "virSorter2_sbatch_time")  #(days-hours)
        virSorter2_cpus_per_task = extract_args(args, "virSorter2_cpus_per_task", Int64, 20, 1, 200)
        virSorter2_sbatch_mem = extract_args(args, "virSorter2_sbatch_mem")
    end

    vs2_score = "$out_vs_D/final-viral-score.tsv" |> TableP
    vs2_boundary = "$out_vs_D/final-viral-boundary.tsv" |> TableP
    vs2_cont_f = "$(out_vs_D)/final-viral-combined.fa" |> FnaP

    postvs2_nonintegrated_df = "$virSorter2_D/$(sp)_02_postVS2_nonintegrated_df.tsv" |> TableP
    postvs2_integrated_df = "$virSorter2_D/$(sp)_02_postVS2_integrated_df.tsv" |> TableP
    #postvs2_integrated_fnaf = "$VirSorter_D/$(sp)_02_postVS2_integrated_fna.fna" |> FnaP

    if signal == "use_external"
        virSorter2 = nothing
    else
        virSorter2 = WrapCmd(; cmd = RunVirSorterCmd(inref, out_vs_D, virSorter2DB_p, virSorter2_high_confidence_only, virSorter2_cpus_per_task, min_score, min_contig_length), 
                                log_p = "$(virSorter2Logs_D)/virSorter2.log", err_p = "$(virSorter2Logs_D)/virSorter2.err", 
                                exit_p = "$(virSorter2Logs_D)/virSorter2.exit", env = virSorter2_env,
                                sbatch_maxtime = virSorter2_sbatch_time, sbatch_cpuspertask = virSorter2_cpus_per_task, sbatch_mem = virSorter2_sbatch_mem)
    end
    
    proj = ProjVirSorter2(virSorter2_D, ext_res, ext_res_D, virSorter2, vs2_score, vs2_boundary, vs2_cont_f, postvs2_nonintegrated_df, postvs2_integrated_df)#, postvs2_integrated_fnaf)

    return proj
end

function set_ProjVIBRANT(pd::String, sampleName::String, inref::FnaP, args::Vector{String}, signal::String, min_contig_length::Int64)
    sp = "ALL-01d"
    vibrant_D = "$sampleName/$(sp)_vibrant"
    my_mkpath(["$(pd)/$(vibrant_D)"])

    if signal == "use_external"
        ext_res = true
        ext_res_D = extract_inPaths(args, "vibrant_res")
        out_vibrant = ""
    else
        ext_res = false
        ext_res_D = nothing
        out_vibrant = "$(vibrant_D)/$(sp)_01_vibrant_out"
        vibrantLogs_D = "$out_vibrant/logs"
        my_mkpath(["$(pd)/$(vibrantLogs_D)"])

        vibrant_env = extract_args(args, "vibrant_env")
        vibrant_p = extract_inPaths(args, "vibrant_p")
        vibrantDB = extract_inPaths(args, "vibrantDB_p")
        #min_length = extract_args(args, "vibrant_min_length", Int64, 1000, 1, 2099000)
        vibrant_sbatch_time = extract_args(args, "vibrant_sbatch_time")
        vibrant_cpus_per_task = extract_args(args, "vibrant_cpus_per_task", Int64, 20, 1, 200)
        vibrant_sbatch_mem = extract_args(args, "vibrant_sbatch_mem")
    end
    
    vib_out_integrated_tsv = "$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_results_$(sampleName)/VIBRANT_integrated_prophage_coordinates_$(sampleName).tsv" |> TableP
    vib_out_integrated_fna = "$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_phages_$(sampleName)/$(sampleName).phages_lysogenic.fna" |> FnaP
    vib_out_nonintegrated_fna = "$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_phages_$(sampleName)/$(sampleName).phages_lytic.fna" |> FnaP
    postvib_nonintegrated_df = "$vibrant_D/$(sp)_02_postvibrant_nonintegrated_df.tsv" |> TableP
    postvib_integrated_df = "$vibrant_D/$(sp)_02_postvibrant_integrated_df.tsv" |> TableP

    todelete = ["$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_figures_$(sampleName)",
                "$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_HMM_tables_parsed_$(sampleName)",
                "$out_vibrant/VIBRANT_$(sampleName)/VIBRANT_HMM_tables_unformatted_$(sampleName)",
                "$out_vibrant/VIBRANT_$(sampleName)/$(sampleName).prodigal.faa",
                "$out_vibrant/VIBRANT_$(sampleName)/$(sampleName).prodigal.ffn",
                "$out_vibrant/VIBRANT_$(sampleName)/$(sampleName).prodigal.gff",
                "$out_vibrant/VIBRANT_log_annotation_$(sampleName).log"]

    if signal == "use_external"
        vibrant = nothing
    else
        vibrant = WrapCmd(; cmd = RunVibrantCmd(vibrant_p, inref, out_vibrant, vibrantDB, vibrant_cpus_per_task, min_contig_length), 
                            log_p = "$(vibrantLogs_D)/vibrant.log", err_p = "$(vibrantLogs_D)/vibrant.err", 
                            exit_p = "$(vibrantLogs_D)/vibrant.exit", env = vibrant_env,
                            sbatch_maxtime = vibrant_sbatch_time, sbatch_cpuspertask = vibrant_cpus_per_task, sbatch_mem = vibrant_sbatch_mem)
    end

    proj = ProjVibrant(vibrant_D, ext_res, ext_res_D, vibrant, vib_out_integrated_tsv, vib_out_integrated_fna, vib_out_nonintegrated_fna, postvib_nonintegrated_df, postvib_integrated_df, todelete)

    return proj
end

function set_ProjViralVerify(pd::String, sampleName::String, inref::FnaP, args::Vector{String}, signal::String) #pd::String
    sp = "ALL-01e"
    vv_D = "$(sampleName)/$(sp)_viralVerify"     # relative folder
    my_mkpath(["$(pd)/$(vv_D)"])                 # absolute path

    if signal == "use_external"
        ext_res = true
        ext_res_D = extract_inPaths(args, "viralVerify_res")
        out_vv_D = ""
    else
        ext_res = false
        ext_res_D = nothing
        out_vv_D = "$(vv_D)/$(sp)_01_viralVerify_out"
        vvLogs_D = "$(out_vv_D)/Logs"
        my_mkpath(["$(pd)/$(vvLogs_D)"])

        vv_env = extract_args(args, "viralVerify_env")
        vv_p = extract_inPaths(args, "viralVerify_p")
        vv_DB = extract_inPaths(args, "viralVerifyDB_p")
        viralVerify_sbatch_time = extract_args(args, "viralVerify_sbatch_time")
        viralVerify_cpus_per_task = extract_args(args, "viralVerify_cpus_per_task", Int64, 20, 1, 200)
        viralVerify_sbatch_mem = extract_args(args, "viralVerify_sbatch_mem")
    end

    vv_th = extract_args(args, "viralVerfify_threshold", Float64, 7.0, 1.0, 1000.0)
    vv_out_tsv = "$(out_vv_D)/$(sampleName)_result_table.csv"  |> TableP
    postVV = "$(vv_D)/$(sp)_02_postViralVerify.tsv" |> TableP
    
    if signal == "use_external"
        vv = nothing
    else
        vv = WrapCmd(; cmd = RunViralVerify(vv_p, inref, out_vv_D, vv_DB, Int64(vv_th), viralVerify_cpus_per_task),
                    log_p = "$(vvLogs_D)/viraVerify.log", err_p = "$(vvLogs_D)/viraVerify.err", exit_p = "$(vvLogs_D)/viraVerify.exit", env = vv_env,
                    sbatch_maxtime = viralVerify_sbatch_time, sbatch_cpuspertask = viralVerify_cpus_per_task, sbatch_mem = viralVerify_sbatch_mem)
    end

    proj = ProjViralVerify(vv_D, ext_res, ext_res_D, vv, vv_out_tsv, vv_th, postVV)

    return proj
end


function set_ProjCheckVNonintegrated(pd::String, input_dfs_2_aggregate::Vector{TableP}, args::Vector{String}, sampleName::String)
    sp = "N-02"
    checkVNonintegrated_D = "$(sampleName)/$(sp)_checkV_Nonintegrated"
    outcheckV_D = "$(checkVNonintegrated_D)/$(sp)_01_checkV"
    checkVNonintegratedLogs_D = "$outcheckV_D/logs"
    my_mkpath(["$(pd)/$(checkVNonintegrated_D)", "$(pd)/$(outcheckV_D)", "$(pd)/$(checkVNonintegratedLogs_D)"])

    output_aggregated_df = "$(checkVNonintegrated_D)/$(sp)_00_aggregated_nonintegrated_virus_contigs.tsv" |> TableP 
    output_aggregated_fna = "$(checkVNonintegrated_D)/$(sp)_00_aggregated_nonintegrated_virus_contigs.fna" |> FnaP

    checkV_nonintegrated_summary_df = "$(outcheckV_D)/quality_summary.tsv" |> TableP
    checkV_out_nonintegrated_complete_df = "$(outcheckV_D)/complete_genomes.tsv" |> TableP
    checkV_out_provir_fna = "$(outcheckV_D)/proviruses.fna" |> FnaP
    postcheckV_nonintegrated_df_p = "$(checkVNonintegrated_D)/$(sp)_02_postcheckV_nonintegrated_df.tsv" |> TableP
    postcheckV_nonintegrated_fna = "$(checkVNonintegrated_D)/$(sp)_02_postcheckV_nonintegrated.fna" |> FnaP
    postcheckV_nonintegrated_fna_trimmed_DTR = "$(checkVNonintegrated_D)/$(sp)_02_postcheckV_nonintegrated_trimmed_DTR.fna" |> FnaP
    postcheckV_integrated_df_p = "$(checkVNonintegrated_D)/$(sp)_02_postcheckV_integrated.tsv" |> TableP

    checkV_env = extract_args(args, "checkv_env")
    checkVDB_p = extract_args(args, "checkvDB_p")

    # resources
    checkv_sbatch_time = extract_args(args, "checkv_sbatch_time")
    checkv_cpus_per_task = extract_args(args, "checkv_cpus_per_task", Int64, 2, 1, 200)
    checkv_sbatch_mem = extract_args(args, "checkv_sbatch_mem")

    proj = ProjCheckVNonintegrated(checkVNonintegrated_D, input_dfs_2_aggregate, output_aggregated_df, output_aggregated_fna,
                                 WrapCmd(; cmd = RunCheckVCmd(output_aggregated_fna, outcheckV_D, checkVDB_p, checkv_cpus_per_task), 
                                        log_p = "$(checkVNonintegratedLogs_D)/checkVNonintegrated.log", 
                                        err_p = "$(checkVNonintegratedLogs_D)/checkVNonintegrated.err", 
                                        exit_p = "$(checkVNonintegratedLogs_D)/checkVNonintegrated.exit", env = checkV_env,
                                        sbatch_maxtime = checkv_sbatch_time, sbatch_cpuspertask = checkv_cpus_per_task, sbatch_mem = checkv_sbatch_mem),
                                        checkV_nonintegrated_summary_df, checkV_out_nonintegrated_complete_df, checkV_out_provir_fna, 
                                        postcheckV_nonintegrated_df_p, postcheckV_nonintegrated_fna, postcheckV_nonintegrated_fna_trimmed_DTR, postcheckV_integrated_df_p) 
    return proj
end

function set_ProjCheckVIntegrated(pd::String, input_dfs_2_aggregate::Vector{TableP}, predictors::Vector{Symbol}, args::Vector{String}, sampleName::String)
    sp = "I-02"
    checkVIntegrated_D = "$(sampleName)/$(sp)_checkV_Integrated"
    
    checkV_env = extract_args(args, "checkv_env")
    checkVDB_p = extract_args(args, "checkvDB_p")

    # resources
    checkv_sbatch_time = extract_args(args, "checkv_sbatch_time")
    checkv_cpus_per_task = extract_args(args, "checkv_cpus_per_task", Int64, 2, 1, 200)
    checkv_sbatch_mem = extract_args(args, "checkv_sbatch_mem")

    output_aggregated_int_contigs_fna = "$(checkVIntegrated_D)/$(sp)_00_aggregated_integrated_virus_contigs.fna" |> FnaP
    output_aggregated_df = "$checkVIntegrated_D/$(sp)_00_aggregated_integrated_virus_All_regions.tsv" |> TableP 
    output_aggregated_fna = "$checkVIntegrated_D/$(sp)_00_aggregated_integrated_virus_All_regions.fna" |> FnaP

    # 1
    checkVIntegrated_D1 = "$(checkVIntegrated_D)/$(sp)_01_checkVIntegrated_1"
    checkVIntegratedLogs_D1 = "$(checkVIntegrated_D1)/logs"
    my_mkpath(["$(pd)/$(checkVIntegrated_D1)", "$(pd)/$(checkVIntegratedLogs_D1)"])

    
    checkV1 = WrapCmd(; cmd = RunCheckVCmd(output_aggregated_fna, checkVIntegrated_D1, checkVDB_p, checkv_cpus_per_task), 
                        log_p = "$(checkVIntegratedLogs_D1)/checkVIntegrated.log", err_p = "$(checkVIntegratedLogs_D1)/checkVIntegrated.err", 
                        exit_p = "$(checkVIntegratedLogs_D1)/checkVIntegrated.exit", env = checkV_env,
                        sbatch_maxtime = checkv_sbatch_time, sbatch_cpuspertask = checkv_cpus_per_task, sbatch_mem = checkv_sbatch_mem)

    checkV1_out_provir_fna = "$(checkVIntegrated_D1)/proviruses.fna" |> FnaP
    postcheckv1_integrated_cor_df = "$(checkVIntegrated_D)/$(sp)_02_postcheckV1_integrated_cordf.tsv" |> TableP
    postcheckv1_integrated_cor_withmergedprovirIDs_df = "$(checkVIntegrated_D)/$(sp)_02_postcheckV1_integrated_withmergedprovirIDs_cordf.tsv" |> TableP

    merged_integrated_DF = "$(checkVIntegrated_D)/$(sp)_03_merged_integrated.tsv" |> TableP
    merged_integrated_fna = "$(checkVIntegrated_D)/$(sp)_03_merged_integrated.fna" |> FnaP

    # 2
    checkVIntegrated_D2 = "$(checkVIntegrated_D)/$(sp)_04_checkVIntegrated_2"
    checkVIntegratedLogs_D2 = "$(checkVIntegrated_D2)/logs"

    my_mkpath(["$(pd)/$(checkVIntegrated_D2)", "$(pd)/$(checkVIntegratedLogs_D2)"])

    checkV2 = WrapCmd(; cmd = RunCheckVCmd(merged_integrated_fna, checkVIntegrated_D2, checkVDB_p, checkv_cpus_per_task), 
                                log_p = "$(checkVIntegratedLogs_D2)/checkVIntegrated.log", err_p = "$(checkVIntegratedLogs_D2)/checkVIntegrated.err", 
                                exit_p = "$(checkVIntegratedLogs_D2)/checkVIntegrated.exit", env = checkV_env,
                                sbatch_maxtime = checkv_sbatch_time, sbatch_cpuspertask = checkv_cpus_per_task, sbatch_mem = checkv_sbatch_mem)

    checkV2_out_integrated_df = "$(checkVIntegrated_D2)/quality_summary.tsv" |> TableP

    postcheckV2_integrated_df_p = "$(checkVIntegrated_D)/$(sp)_05_postcheckV2_integrated_df_p.tsv" |> TableP

    proj = ProjCheckVIntegrated(
                    checkVIntegrated_D,
                    input_dfs_2_aggregate,
                    #input_fnas_2_aggregate,
                    output_aggregated_int_contigs_fna,
                    output_aggregated_df,
                    output_aggregated_fna,
                    checkV1, 
                    checkV1_out_provir_fna,
                    postcheckv1_integrated_cor_df,
                    postcheckv1_integrated_cor_withmergedprovirIDs_df,
                    predictors,
                    merged_integrated_DF,
                    merged_integrated_fna,
                    checkV2,
                    checkV2_out_integrated_df, 
                    postcheckV2_integrated_df_p)

    return proj
end

function set_ProjPhaTYPNonintegrated(pd::String, checkV_NonIntegrated::ProjCheckVNonintegrated, args::Vector{String}, sampleName::String, min_contig_length::Int64)
    sp = "N-03"
    phaTYPNonintegrated_D = "$sampleName/$(sp)_phaTYPNonintegrated"
    out_D = "$(phaTYPNonintegrated_D)/$(sp)_01_phaTYP_out"
    phaTYPNonintegratedLogs_D = "$phaTYPNonintegrated_D/logs"
    my_mkpath(["$(pd)/$(phaTYPNonintegrated_D)", "$(pd)/$(out_D)", "$(pd)/$(phaTYPNonintegratedLogs_D)"])

    phaTYP_env = extract_args(args, "phaTYP_env")
    phaTYP_p = extract_args(args, "phaTYP_p")
    phaTYPdb_p = extract_args(args, "phaTYPdb_p")
    phaTYPparam_p = extract_args(args, "phaTYPparam_p")
    res_suffix = "results"

    # resources
    phaTYP_sbatch_time = extract_args(args, "phaTYP_sbatch_time")
    phaTYP_cpus_per_task = extract_args(args, "phaTYP_cpus_per_task", Int64, 2, 1, 200)
    phaTYP_sbatch_mem = extract_args(args, "phaTYP_sbatch_mem")

    #contig len
    #minContigLen = extract_args(args, "phaTYP_minlength", Int64, 100, 1, 2099000)

    phaTYPNonintegrated_out_df = "$(out_D)/$(res_suffix)/phatyp_prediction.csv" |> TableP
    mergedPostCheckV_PhaTYP_p = "$(phaTYPNonintegrated_D)/$(sp)_02_NonINT_mergedPostCheckV_PhaTYP.tsv" |> TableP

    proj = ProjPhaTYP(phaTYPNonintegrated_D, WrapCmd(; cmd = RunPhaTYPCmd(phaTYP_p, checkV_NonIntegrated.postcheckV_nonintegrated_fna, phaTYPdb_p, phaTYPparam_p,
                                                        out_D, res_suffix, min_contig_length, phaTYP_cpus_per_task), 
                                                        log_p = "$(phaTYPNonintegratedLogs_D)/phaTYP_NonIntegrated.log", err_p = "$(phaTYPNonintegratedLogs_D)/phaTYP_NonIntegrated.err", 
                                                        exit_p = "$(phaTYPNonintegratedLogs_D)/phaTYP_NonIntegrated.exit", env = phaTYP_env,
                                                        sbatch_maxtime = phaTYP_sbatch_time, sbatch_cpuspertask = phaTYP_cpus_per_task, sbatch_mem = phaTYP_sbatch_mem),
                    phaTYPNonintegrated_out_df, mergedPostCheckV_PhaTYP_p)

    return proj
end

#=
function set_ProjPhaTYPIntegrated(pd::String, checkV_Integrated::ProjCheckVIntegrated, num_threads::Int64, args::Vector{String})
    sp = "0xI"
    phaTYP_integrated_D = "$pd/$(sp)_phaTYP_integrated"
    out_D = "$(phaTYP_integrated_D)/$(sp)_01_phaTYP_out"
    phaTYP_integratedLogs_D = "$phaTYP_integrated_D/logs"
    my_mkpath([phaTYP_integrated_D, out_D, phaTYP_integratedLogs_D])

    phaTYP_env = extract_args(args, "phaTYP_env")
    phaTYP_p = extract_args(args, "phaTYP_p")
    phaTYPdb_p = extract_args(args, "phaTYPdb_p")
    phaTYPparam_p = extract_args(args, "phaTYPparam_p")
    res_suffix = "results"

    #contig len
    DVF_minContigLen = extract_args(args, "DVF_minContigLen", Int64, 100, 1, 2099000)

    phaTYP_integrated_out_df = "$(out_D)/$(res_suffix)/phatyp_prediction.csv" |> TableP

    proj = ProjPhaTYP(phaTYP_integrated_D, WrapCmd(; cmd = RunPhaTYPCmd(phaTYP_p, checkV_Integrated.merged_integrated_fna, phaTYPdb_p, phaTYPparam_p,
                                        out_D, res_suffix, DVF_minContigLen, num_threads), 
            log_p = "$(phaTYP_integratedLogs_D)/phaTYP_integrated.log", err_p = "$(phaTYP_integratedLogs_D)/phaTYP_integrated.err", exit_p = "$(phaTYP_integratedLogs_D)/phaTYP_integrated.exit", env = phaTYP_env),
            phaTYP_integrated_out_df)

    return proj
end =#

function set_ProjDetectMixedViruses(pd::String, inDf_NonInt::TableP, inDf_Int::TableP, sampleName::String)
    sp = "M-02"
    mixed_D = "$(sampleName)/$(sp)_Detection"
    my_mkpath(["$(pd)/$(mixed_D)"])

    outDf_Int = "$(mixed_D)/$(sp)_Mixed_Int_Viruses.tsv" |> TableP
    #outFna_Int = "$(mixed_D)/$(sp)_Mixed_Int_Viruses.fna" |> FnaP
    outDf_NonInt = "$(mixed_D)/$(sp)_Mixed_NonInt_Viruses.tsv" |> TableP
    #outFna_NonInt = "$(mixed_D)/$(sp)_Mixed_NonInt_Viruses.fna" |> FnaP
    #outmixed_df = "$(mixed_D)/$(sp)_Mixed_All_Viruses.tsv" |> TableP
    #outmixed_fna = "$(mixed_D)/$(sp)_Mixed_all_Viruses.fna" |> FnaP

    proj = ProjDetectMixedViruses(mixed_D, inDf_Int, inDf_NonInt, outDf_Int, outDf_NonInt) 

    return proj
end

function set_ProjFinalThresholding(pd::String, args::Vector{String}, sp::String, inFna::FnaP, inFna_trimmed_DTR::Union{Missing, FnaP}, inTsv::TableP, predictors::Vector{Symbol}, sampleName::String, predictors2::Union{Vector{Symbol}, Missing} = missing)
    npd = "$(sampleName)/$(sp)_FinalThresholding"
    my_mkpath(["$(pd)/$(npd)"])

    if sp == "N-04"
        th_num_predictors_CheckV_NA = extract_args(args, "NONInt_th_num_predictors_CheckV_NA", Int64, 3, 1, 5)

        th_num_predictors_CheckV_AAIHighConf = extract_args(args, "NONInt_th_num_predictors_CheckV_AAIHighConf", Int64, 1, 1, 5)
        th_completeness_CheckV_AAIHighConf = extract_args(args, "NONInt_th_completeness_CheckV_AAIHighConf", Float64, 30.0, 1.0, 100.0)

        th_num_predictors_CheckV_AAIMediumConf = extract_args(args, "NONInt_th_num_predictors_CheckV_AAIMediumConf", Int64, 2, 1, 5)
        th_completeness_CheckV_AAIMediumConf = extract_args(args, "NONInt_th_completeness_CheckV_AAIMediumConf", Float64, 10.0, 1.0, 100.0)

        th_num_predictors_CheckV_HMM = extract_args(args, "NONInt_th_num_predictors_CheckV_HMM", Int64, 2, 1, 5)
        th_completeness_CheckV_HMM = extract_args(args, "NONInt_th_completeness_CheckV_HMM", Float64, 10.0, 1.0, 100.0)

        th_num_predictors_CheckV_DTR_ITR_AAI = extract_args(args, "NONInt_th_num_predictors_CheckV_DTR_ITR_AAI", Int64, 1, 1, 5)
        th_num_predictors_CheckV_DTR_ITR_HMM = extract_args(args, "NONInt_th_num_predictors_CheckV_DTR_ITR_HMM", Int64, 1, 1, 5)

    elseif sp == "I-03"
        th_num_predictors_CheckV_NA = extract_args(args, "Int_th_num_predictors_CheckV_NA", Int64, 3, 1, 5)

        th_num_predictors_CheckV_AAIHighConf = extract_args(args, "Int_th_num_predictors_CheckV_AAIHighConf", Int64, 1, 1, 5)
        th_completeness_CheckV_AAIHighConf = extract_args(args, "Int_th_completeness_CheckV_AAIHighConf", Float64, 30.0, 1.0, 100.0)

        th_num_predictors_CheckV_AAIMediumConf = extract_args(args, "Int_th_num_predictors_CheckV_AAIMediumConf", Int64, 2, 1, 5)
        th_completeness_CheckV_AAIMediumConf = extract_args(args, "Int_th_completeness_CheckV_AAIMediumConf", Float64, 10.0, 1.0, 100.0)

        th_num_predictors_CheckV_HMM = extract_args(args, "Int_th_num_predictors_CheckV_HMM", Int64, 2, 1, 5)
        th_completeness_CheckV_HMM = extract_args(args, "Int_th_completeness_CheckV_HMM", Float64, 10.0, 1.0, 100.0)

        th_num_predictors_CheckV_DTR_ITR_AAI = missing
        th_num_predictors_CheckV_DTR_ITR_HMM = missing
    elseif sp == "M-03"
        #= old, I did not implement yet the final thresholding for mixed viruses
        th_num_predictors = extract_args(args, "th_num_predictors_Mixed_nonint", Int64, 1, 1, 4)
        th_num_predictors2 = extract_args(args, "th_num_predictors_Mixed_int", Int64, 1, 1, 4)
        th_checkV_completeness = extract_args(args, "th_checkV_completeness_Mixed_nonint", Float64, 30.0, 1.0, 100.0)
        th_checkV_completeness2 = extract_args(args, "th_checkV_completeness_Mixed_int", Float64, 30.0, 1.0, 100.0)
        th_checkV_contamination = extract_args(args, "th_checkV_contamination_Mixed_nonint", Float64, 50.0, 0.0, 100.0)
        th_checkV_contamination2 = extract_args(args, "th_checkV_contamination_Mixed_int", Float64, 50.0, 0.0, 100.0) =#
    end

    out_fna = "$(npd)/$(sp)_Final_Viruses.fna" |> FnaP
    if !ismissing(inFna_trimmed_DTR)     
        outFna_trimmed_DTR = "$(npd)/$(sp)_Final_Viruses_trimmed_DTR.fna" |> FnaP
    else
        outFna_trimmed_DTR = missing
    end
    out_tsv = "$(npd)/$(sp)_Final_Viruses_df.tsv" |> TableP

    proj = FinalThresholding(npd, inFna, inFna_trimmed_DTR, inTsv, predictors, predictors2, 
                            th_num_predictors_CheckV_NA,
                            th_num_predictors_CheckV_AAIHighConf, th_completeness_CheckV_AAIHighConf,
                            th_num_predictors_CheckV_AAIMediumConf, th_completeness_CheckV_AAIMediumConf,
                            th_num_predictors_CheckV_HMM, th_completeness_CheckV_HMM,
                            th_num_predictors_CheckV_DTR_ITR_AAI, th_num_predictors_CheckV_DTR_ITR_HMM,
                            out_fna, outFna_trimmed_DTR, out_tsv)

    return proj
end


function ProjSViP_fun(args::Vector{String})

    # num_threads
    num_threads = extract_args(args, "num_threads", Int64, 1, 1, 40)

    # continue or not
    cont = extract_args(args, "continue", Bool, "false")

    use_slurm = extract_args(args, "use_slurm", Bool, "false")

    stop_after_initial_predictors = extract_args(args, "stop_after_initial_predictors", Bool, "false")

    # inref and sample set
    inref = extract_inFiles(args, "inref", BioS_Gen.ALLOWED_EXT["FnaP"]) |> FnaP
    sampleName = getFileName(inref.p) 
    sample_set = extract_args(args, "sample_set")
    #println(typeof(sample_set))

    # min_contig_length 
    min_contig_length = extract_args(args, "min_contig_length", Int64, 1000, 1000, 10000000)

    # folders
    #pd_prefix = "$(extract_args(args, "pd_prefix"))/"
    pd = "$(extract_args(args, "pd_prefix"))/" #pd = "$(pd_prefix)/$(sampleName)"


    #region dosteps
    if cont == false
        #rm_mkpaths([pd]) 
        do_pd("$(pd)/$(sampleName)")

        dosteps = Dict(
            "sel_contig_length" => WorkflowStatus("do", "not_done"),
            #"shape_contigs" => WorkflowStatus("do", "not_done"),
            "genomad" => WorkflowStatus(
                setsignal(extract_args(args, "genomad_signal", ALLOWED_VALS_PROJ["signal"]), "GeNomad"), 
                "not_done"),
            "DVF" => WorkflowStatus(
                setsignal(extract_args(args, "DVF_signal", ("do", "dont", "use", "ignore", "remove")), "DVF"), 
                "not_done"),
            "virSorter2" => WorkflowStatus(
                setsignal(extract_args(args, "virSorter2_signal", ALLOWED_VALS_PROJ["signal"]), "virSorter2"),
                 "not_done"),
            "vibrant" => WorkflowStatus(
                setsignal(extract_args(args, "vibrant_signal", ALLOWED_VALS_PROJ["signal"]), "VIBRANT"),
                 "not_done"),
            "viralVerify" => WorkflowStatus(
                setsignal(extract_args(args, "viralVerify_signal", ALLOWED_VALS_PROJ["signal"]), "viralVerify"),
                "not_done"),
            "checkV_NonIntegrated" => WorkflowStatus("do", "not_done"),
            "checkV_Integrated" => WorkflowStatus("do", "not_done"),
            "phaTYP_nonintegrated" => WorkflowStatus("do", "not_done"),
            #"phaTYP_integrated" => WorkflowStatus("do", "not_done"),
            "detect_mixed_viruses" => WorkflowStatus("do", "not_done"),
            "final_thresholding_NonIntegrated" => WorkflowStatus("do", "not_done"),
            "final_thresholding_Integrated" => WorkflowStatus("do", "not_done")
            #"final_thresholding_Mixed" => WorkflowStatus("do", "not_done")
        )
        
        if (dosteps["genomad"].signal == "dont" &&  dosteps["DVF"].signal == "dont"
            && dosteps["virSorter2"].signal == "dont" && dosteps["vibrant"].signal == "dont"
            && dosteps["viralVerify"].signal == "dont")
             
            println("At least one predictor form the INITIAL PREDICTION MODULE needs to have its signal 'do'. Exiting DoviP.")
            exit()
        end

        proj = ProjSViP(use_slurm = use_slurm, continue_project = cont, pd = pd, inref = inref, sampleName = sampleName, sampleSet = sample_set, dosteps = dosteps)
    else
        
        #sampleName = "outMFasta_small"
        #inref="/home/cmoraru/MY_SPACE/WorkflowsJl/TESTS/test_in_files/DoViP/outMFasta_small.fasta" |> FnaP
        #pd = "/data3/CLM_projs/TEST_Workflows/DoViP_outMFasta_small"
        proj = load_proj("$(pd)/$(sampleName)/sproj.binary", pd, inref, sampleName)

        dosteps = proj.dosteps
        prev_dosteps = deepcopy(proj.dosteps)

        # at continue project, the shape of contigs will always be reused
        if min_contig_length != proj.contig_length.min_contig_length
            
            println("The minimum contig length can't be modified at project restart. The previous value for min_contig_length was $(proj.contig_length.min_contig_length). 
                    DoViP will exit now. To continue with the analysis, either set the min_contig_length to its previous value or start a completely new project.")

                    exit()
        else
            dosteps["sel_contig_length"].signal = "use"
            #dosteps["shape_contigs"].signal = "use"
        end
        
        # the signal for the predictors is given by the user and it is extracted anew at continue project (the user can continue the project by running other predictors)
        dosteps["genomad"].signal = setsignal(extract_args(args, "genomad_signal", ALLOWED_VALS_PROJ["signal"]), 
                                            "GeNomad"; cont = true, 
                                            oldsignal = prev_dosteps["genomad"].signal,
                                            progress = prev_dosteps["genomad"].progress)
        dosteps["DVF"].signal = setsignal(extract_args(args, "DVF_signal", ALLOWED_VALS_PROJ["signal"]),
                                        "DVF"; cont = true, 
                                        oldsignal = prev_dosteps["DVF"].signal,
                                        progress = prev_dosteps["DVF"].progress)
        dosteps["virSorter2"].signal = setsignal(extract_args(args, "virSorter2_signal", ALLOWED_VALS_PROJ["signal"]),
                                                "virSorter2"; cont = true, 
                                                oldsignal = prev_dosteps["virSorter2"].signal,
                                                progress = prev_dosteps["virSorter2"].progress)
        dosteps["vibrant"].signal = setsignal(extract_args(args, "vibrant_signal", ALLOWED_VALS_PROJ["signal"]),
                                            "VIBRANT"; cont = true, 
                                            oldsignal = prev_dosteps["vibrant"].signal,
                                            progress = prev_dosteps["vibrant"].progress)
        dosteps["viralVerify"].signal = setsignal(extract_args(args, "viralVerify_signal", ALLOWED_VALS_PROJ["signal"]),
                                                "viralVerify"; cont = true, 
                                                oldsignal = prev_dosteps["viralVerify"].signal,
                                                progress = prev_dosteps["viralVerify"].progress)

        # exit if no predictor is active                                    
        if (dosteps["genomad"].signal in ["ignore", "remove", "dont"] &&  dosteps["DVF"].signal in ["ignore", "remove", "dont"]
            && dosteps["virSorter2"].signal in ["ignore", "remove", "dont"] && dosteps["vibrant"].signal in ["ignore", "remove", "dont"] 
            && dosteps["viralVerify"].signal in ["ignore", "remove", "dont"])
                
            println("At least one predictor form the INITIAL PREDICTION MODULE needs to have its signal 'do', 'use' or 'use_external'. Exiting DoviP.")
            exit()
        end
        
        # determine is the results from the INITIAL PREDICTION MODULE and THE CONSENSUS PREDICTION STEPS can be re-used as they are
        touse = 0
        toignore = 0
        for p in ["genomad", "DVF", "virSorter2", "vibrant", "viralVerify"]
            if (dosteps[p].signal == "use" && prev_dosteps[p].signal in ["do", "use", "use_external"] && prev_dosteps[p].progress == "finished")
                touse += 1
            elseif (dosteps[p].signal in ["dont", "remove", "ignore"] && prev_dosteps[p].progress == "not_done")
                toignore += 1
            end
        end


        if (touse + toignore) == 5
            #=(dosteps["genomad"] == "use" && dosteps["DVF"] == "use" && dosteps["virSorter2"] == "use" && dosteps["vibrant"] == "use" 
            && prev_dosteps["genomad"] in ["do", "use"] && prev_dosteps["DVF"] in ["do", "use"] && prev_dosteps["virSorter2"] in ["do", "use"]
            && prev_dosteps["vibrant"] in ["do", "use"]) =#

            println("
            You have chosen to continue an already existing DoViP project by reusing the results from the INITIAL PREDICTION MODULE which were calculated in the previous run.
            The followiung rules apply: 
            
            In the case of the Initial Prediction step, all the previous results from the predictors will be used. Any newly given parameters related to this step will NOT be applied.
            In the case of the Consensus Prediction step, all the previous results from these steps will be used if their progress in the previous run was 'finished'.
            In the case of the Final Thresholding steps, all the previous results will be removed and new results will calculated using the current parameters.")

            #= I don't have to use the setstep_primary! function because the status of all primary steps (predictors) is finished 
            (it was checked by the setsignal function before setting the signal to 'use') and I'm re-using their results. =#
                     
            #= I'm setting manually the intermediary steps to "use". Here I have to use the setstep_intermediary! function to first check 
            their progress in the previous run and delete the results from dependent steps =#

            dosteps["checkV_NonIntegrated"].signal = "use"
            dosteps = setstep_intermediary!(proj, "checkV_NonIntegrated", dosteps, "$(pd)/$(proj.checkV_NonIntegrated.pd)",
                                            ["$(pd)/$(proj.phaTYP_nonintegrated.pd)", "$(pd)/$(proj.final_thresholding_NonIntegrated.pd)", 
                                            "$(pd)/$(proj.checkV_Integrated.pd)", "$(pd)/$(proj.final_thresholding_Integrated.pd)",
                                            "$(pd)/$(proj.detect_mixed_viruses.pd)"]; #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"
                                            pds2remove_use = ["$(pd)/$(proj.final_thresholding_NonIntegrated.pd)", "$(pd)/$(proj.final_thresholding_Integrated.pd)"], logfun = printProjSViP) #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"

            dosteps["checkV_Integrated"].signal = "use"
            dosteps = setstep_intermediary!(proj, "checkV_Integrated", dosteps, "$(pd)/$(proj.checkV_Integrated.pd)",
                                            ["$(pd)/$(proj.final_thresholding_Integrated.pd)", "$(pd)/$(proj.detect_mixed_viruses.pd)"];  #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"
                                            pds2remove_use = ["$(pd)/$(proj.final_thresholding_Integrated.pd)"], logfun = printProjSViP) #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"

            dosteps["phaTYP_nonintegrated"].signal = "use"
            dosteps = setstep_intermediary!(proj, "phaTYP_nonintegrated", dosteps, "$(pd)/$(proj.phaTYP_nonintegrated.pd)",
                                            ["$(pd)/$(proj.final_thresholding_NonIntegrated.pd)"]; 
                                            pds2remove_use = ["$(pd)/$(proj.final_thresholding_NonIntegrated.pd)"], logfun = printProjSViP)
            
            dosteps["detect_mixed_viruses"].signal = "use"
            dosteps = setstep_intermediary!(proj, "detect_mixed_viruses", dosteps, "$(pd)/$(proj.detect_mixed_viruses.pd)",
                                            ["$(pd)/$(proj.detect_mixed_viruses.pd)"]; #["$(pd)/$(proj.final_thresholding_Mixed.pd)"]
                                            pds2remove_use = ["$(pd)/$(proj.detect_mixed_viruses.pd)"], logfun = printProjSViP)   #["$(pd)/$(proj.final_thresholding_Mixed.pd)"]                  

        else 
            println("
            You have chosen to continue an already existing DoViP project by re-runing, removing or ignoring at least one of the 4 predictors in the INITIAL PREDICTION STEP.
            The followiung rules apply: 

            In the case of the Initial Prediction step, only the results from the predictors with signal 'do' or 'use' will be used. 
            The previous results from any predictor with signal 'do' will be removed and new results will be calculated, using the current predictor-specific parameters.
            The previous results of the post-processing step from any predictor with signal 'use_external' will be removed and the exteral results will be put through the post-processing step, using the current predictor-specific parameters.
            The previous results from any predictors with signal 'remove' will be removed and NO new results will be calculated.
            The preivous results from any predictos with signal 'ignore' will NOT be removed, but also will not be considered in the following steps.
            The previous results from any predictor with signal 'use' will be kept as they are. Any change of the predictor-specific parameters will not be considered.
            
            In the case of the Consensus Prediction and Final Thresholding steps, all the previous results will be removed and new results will calculated, using the current parameters.
            ")

            # if only one of the 4 predictors have must be redone or removed, then all post-prediction steps will be deleted and redone
            todel = vcat("$(pd)/$(proj.checkV_NonIntegrated.pd)", "$(pd)/$(proj.phaTYP_nonintegrated.pd)", "$(pd)/$(proj.final_thresholding_NonIntegrated.pd)",
                            "$(pd)/$(proj.checkV_Integrated.pd)", "$(pd)/$(proj.final_thresholding_Integrated.pd)",
                            "$(pd)/$(proj.detect_mixed_viruses.pd)") #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"
            
            if ismissing(proj.genomad) == false            
                dosteps = setstep_primary!(proj, "genomad", dosteps, "$(pd)/$(proj.genomad.pd)", todel; logfun = printProjSViP)
            else
                dosteps = setstep_primary!(proj, "genomad", dosteps, missing, todel; logfun = printProjSViP)
            end

            if ismissing(proj.DVF) == false
                dosteps = setstep_primary!(proj, "DVF", dosteps, "$(pd)/$(proj.DVF.pd)", todel; logfun = printProjSViP)
            else
                dosteps = setstep_primary!(proj, "DVF", dosteps, missing, todel; logfun = printProjSViP)
            end

            if ismissing(proj.virSorter2) == false
                dosteps = setstep_primary!(proj, "virSorter2", dosteps, "$(pd)/$(proj.virSorter2.pd)", todel; logfun = printProjSViP)
            else
                dosteps = setstep_primary!(proj, "virSorter2", dosteps, missing, todel; logfun = printProjSViP)
            end

            if ismissing(proj.vibrant) == false
                dosteps = setstep_primary!(proj, "vibrant", dosteps, "$(pd)/$(proj.vibrant.pd)", todel; logfun = printProjSViP)
            else
                dosteps = setstep_primary!(proj, "vibrant", dosteps, missing, todel; logfun = printProjSViP)
            end

            if ismissing(proj.viralVerify) == false
                dosteps = setstep_primary!(proj, "viralVerify", dosteps, "$(pd)/$(proj.viralVerify.pd)", todel; logfun = printProjSViP)
            else
                dosteps = setstep_primary!(proj, "viralVerify", dosteps, missing, todel; logfun = printProjSViP)
            end

           # I need to set the signal for the steps below to "do", becase they might be set to use form a previous run
            dosteps["checkV_NonIntegrated"].signal = "do"
            dosteps = setstep_intermediary!(proj, "checkV_NonIntegrated", dosteps, "$(pd)/$(proj.checkV_NonIntegrated.pd)",
                                            ["$(pd)/$(proj.phaTYP_nonintegrated.pd)", "$(pd)/$(proj.final_thresholding_NonIntegrated.pd)", 
                                            "$(pd)/$(proj.checkV_Integrated.pd)", "$(pd)/$(proj.final_thresholding_Integrated.pd)",
                                            "$(pd)/$(proj.detect_mixed_viruses.pd)"]; logfun = printProjSViP) #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"
            
            dosteps["checkV_Integrated"].signal = "do"
            dosteps = setstep_intermediary!(proj, "checkV_Integrated", dosteps, "$(pd)/$(proj.checkV_Integrated.pd)",
                                            ["$(pd)/$(proj.final_thresholding_Integrated.pd)", "$(pd)/$(proj.detect_mixed_viruses.pd)"]; logfun = printProjSViP) #, "$(pd)/$(proj.final_thresholding_Mixed.pd)"
            
            dosteps["phaTYP_nonintegrated"].signal = "do"
            dosteps = setstep_intermediary!(proj, "phaTYP_nonintegrated", dosteps, "$(pd)/$(proj.phaTYP_nonintegrated.pd)",
                                            ["$(pd)/$(proj.final_thresholding_NonIntegrated.pd)"]; logfun = printProjSViP) #proj.checkV_Integrated.pd, proj.final_thresholding_Integrated.pd,

            dosteps["detect_mixed_viruses"].signal = "do"
            dosteps = setstep_intermediary!(proj, "detect_mixed_viruses", dosteps, "$(pd)/$(proj.detect_mixed_viruses.pd)", 
                                            ["$(pd)/$(proj.detect_mixed_viruses.pd)"]; logfun = printProjSViP) #["$(pd)/$(proj.final_thresholding_Mixed.pd)"]
        end

        # these final steps allways have their signal set to do, because they will be recalculated at every project run
        dosteps = setstep_final!(proj, "final_thresholding_NonIntegrated", dosteps, "$(pd)/$(proj.final_thresholding_NonIntegrated.pd)"; logfun = printProjSViP)
        dosteps = setstep_final!(proj, "final_thresholding_Integrated", dosteps, "$(pd)/$(proj.final_thresholding_Integrated.pd)"; logfun = printProjSViP)
        #dosteps = setstep_final!(proj, "final_thresholding_Mixed", dosteps, "$(pd)/$(proj.final_thresholding_Mixed.pd)")

        touse = nothing
        toignore = nothing
    end
    #endregion

    #initialize set minimum contig length
    contig_length = initialize_step(dosteps, "sel_contig_length", set_contig_length, (inref, pd, sampleName, min_contig_length), proj.contig_length, cont)

    #initialize shape contigs
    #shape_contigs = initialize_step(dosteps, "shape_contigs", set_shape_contigs, (), proj.shape_contigs, cont)

    predictors_all, predictors_Int = set_predictors(dosteps)


    nonintegrated_viruses_to_agregate = Vector{TableP}()
    integrated_viruses_to_agregate_df = Vector{TableP}()
    #integrated_viruses_to_agregate_fna = Vector{FnaP}()

    #region predict viruses and prophages with geNomad
    genomad = initialize_step(dosteps, "genomad", set_ProjGenomad, (pd, sampleName, contig_length.outref, args, dosteps["genomad"].signal), proj.genomad, cont)
    if dosteps["genomad"].signal in ["do", "use", "use_external"]
        push!(nonintegrated_viruses_to_agregate, genomad.postgenomad_nonintegrated_df)
        push!(integrated_viruses_to_agregate_df, genomad.postgenomad_integrated_df)
    end
    #endregion
       
    #region viruses with DVF
    DVF = initialize_step(dosteps, "DVF", set_ProjDVF, (pd, sampleName, contig_length.outref, args, min_contig_length), proj.DVF, cont)
    if dosteps["DVF"].signal in ["do", "use"]
        push!(nonintegrated_viruses_to_agregate, DVF.postdvf_nonintegrated_df)
    end
    #endregion

    #region predict viruses and prophages with vibrant
    vibrant = initialize_step(dosteps, "vibrant", set_ProjVIBRANT, (pd, sampleName, contig_length.outref, args, dosteps["vibrant"].signal, min_contig_length), proj.vibrant, cont)
    if dosteps["vibrant"].signal in ["do", "use", "use_external"]
        push!(nonintegrated_viruses_to_agregate, vibrant.postvib_nonintegrated_df)
        push!(integrated_viruses_to_agregate_df, vibrant.postvib_integrated_df)
    end
    #endregion

    #region viralVerify
    viralVerify = initialize_step(dosteps, "viralVerify", set_ProjViralVerify, (pd, sampleName, contig_length.outref, args, dosteps["viralVerify"].signal), proj.viralVerify, cont)
    if dosteps["viralVerify"].signal in ["do", "use", "use_external"]
        push!(nonintegrated_viruses_to_agregate, viralVerify.postViralVerify_df)
    end
    #endregion

    #region predict viruses and prophages with virSorter2s
    virSorter2 = initialize_step(dosteps, "virSorter2", set_ProjVirSorter2, (pd, sampleName, contig_length.outref, args, dosteps["virSorter2"].signal, min_contig_length), proj.virSorter2, cont)
    if dosteps["virSorter2"].signal in ["do", "use", "use_external"]
        push!(nonintegrated_viruses_to_agregate, virSorter2.postvs2_nonintegrated_df)
        push!(integrated_viruses_to_agregate_df, virSorter2.postvs2_integrated_df)
    end
    #endregion


    #region Cenote-Taker3 for virus prediction and annotation?
    #endregion

    #region nonintegrated viruses
    if length(nonintegrated_viruses_to_agregate) >= 1
        checkV_NonIntegrated = initialize_step(dosteps, "checkV_NonIntegrated", set_ProjCheckVNonintegrated, (pd, nonintegrated_viruses_to_agregate, args, sampleName), proj.checkV_NonIntegrated, cont)

        push!(integrated_viruses_to_agregate_df, checkV_NonIntegrated.postcheckV_integrated_df_p)

        phaTYP_NonIntegrated = initialize_step(dosteps, "phaTYP_nonintegrated", set_ProjPhaTYPNonintegrated, (pd, checkV_NonIntegrated, args, sampleName, min_contig_length), proj.phaTYP_nonintegrated, cont)
        final_thresholding_NonIntegrated = initialize_step(dosteps, "final_thresholding_NonIntegrated", 
                                                            set_ProjFinalThresholding, 
                                                            (pd, args, "N-04", checkV_NonIntegrated.postcheckV_nonintegrated_fna, checkV_NonIntegrated.postcheckV_nonintegrated_fna_trimmed_DTR,
                                                            phaTYP_NonIntegrated.mergedPostCheckV_PhaTYP_p, predictors_all, sampleName), 
                                                            proj.final_thresholding_NonIntegrated, cont)
    else    
        println("No predictor for non-integrated viruses is activated. If that was not your intention, check your input parameters.")  
        checkV_NonIntegrated = missing
        phaTYP_NonIntegrated = missing
        final_thresholding_NonIntegrated = missing
    end
    #endregion

    #region integrated viruses
    if length(integrated_viruses_to_agregate_df) >= 1 #&& length(integrated_viruses_to_agregate_fna) >= 1
        checkV_Integrated = initialize_step(dosteps, "checkV_Integrated", set_ProjCheckVIntegrated, (pd, integrated_viruses_to_agregate_df, predictors_all, args, sampleName), proj.checkV_Integrated, cont)
        final_thresholding_Integrated = initialize_step(dosteps, "final_thresholding_Integrated", 
                                                        set_ProjFinalThresholding, 
                                                        (pd, args, "I-03", checkV_Integrated.merged_integrated_fna, missing, checkV_Integrated.postcheckV2_integrated_df_p, predictors_all, sampleName),
                                                        proj.final_thresholding_Integrated, cont)
    else     
        println("No predictor for integrated viruses is activated. If that was not your intention, check your input parameters.")  
        checkV_Integrated = missing
        final_thresholding_Integrated = missing
    end
    #endregion


    #region Mixed viruses
    if length(nonintegrated_viruses_to_agregate) >= 1 && length(integrated_viruses_to_agregate_df) >= 1
        detect_mixed_viruses = initialize_step(dosteps, "detect_mixed_viruses", set_ProjDetectMixedViruses, (pd, checkV_NonIntegrated.postcheckV_nonintegrated_df_p, checkV_Integrated.postcheckV2_integrated_df_p, sampleName), proj.detect_mixed_viruses, cont)
        #=finalize_thresholding_Mixed = initialize_step(dosteps, "final_thresholding_Mixed", 
                                                        set_ProjFinalThresholding, 
                                                        (pd, args, "M-03", detect_mixed_viruses.outmixed_fna, detect_mixed_viruses.outmixed_df, predictors_all, sampleName, predictors_Int),
                                                        proj.final_thresholding_Mixed, cont) =#
    else
        detect_mixed_viruses = missing
        #finalize_thresholding_Mixed = missing
    end

    #endregion

    proj = nothing
    sproj = ProjSViP("singleworkflow",
                    use_slurm,
                    cont, 
                    pd,
                    inref, 
                    sampleName, 
                    sample_set,
                    dosteps, 
                    stop_after_initial_predictors,
                    contig_length,
                    #shape_contigs,
                    genomad, 
                    DVF, 
                    virSorter2, 
                    vibrant,
                    viralVerify,
                    checkV_NonIntegrated, 
                    checkV_Integrated, 
                    phaTYP_NonIntegrated, 
                    #phaTYP_integrated, 
                    detect_mixed_viruses,
                    final_thresholding_NonIntegrated,
                    final_thresholding_Integrated,
                    #finalize_thresholding_Mixed
                    )
    serialize("$(pd)/$(sampleName)/sproj.binary", sproj)
    printProjSViP("$(pd)/$(sampleName)/project_parameters_and_status.txt", sproj)

    return sproj
end
