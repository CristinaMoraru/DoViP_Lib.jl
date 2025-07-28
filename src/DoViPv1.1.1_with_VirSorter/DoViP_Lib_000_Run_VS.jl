export run_workflow, run_workflowMDoViP

function run_workflow_genomad(proj::ProjGenomad, signal::String, parentD::String; sbatch::Bool = false)
    
    if signal == "do"
        do_cmd(proj.genomad, "genomad", true, parentD; sbatch = sbatch)
    end

    if signal in ["do", "use"]
        prefix = parentD
    elseif signal == "use_external"
        prefix = proj.ext_res_D
    end

    if (isfile("$(prefix)/$(proj.genomad_out_table_p.p)") && filesize("$(prefix)/$(proj.genomad_out_table_p.p)") > 0) &&
        (isfile("$(prefix)/$(proj.genomad_out_fnap.p)") && filesize("$(prefix)/$(proj.genomad_out_fnap.p)") > 0)

        println("\nPost-processing the geNomad output files.")
        post_genomad(proj, parentD)

    else
        println("\ngeNomad did not produce any output files.")
    end

    if signal == "do"
        rm_path(proj.todelete)
    end

    return nothing
end

function run_workflow_DVF(proj::ProjDVF, signal::String, parentD::String; sbatch::Bool = false)

    if signal == "do"
        println("\nPre-processing the DVF input file to split large contigs.")
        split_contigs = contigSplit(FnaP("$(parentD)/$(proj.input_f.p)"), FnaP("$(parentD)/$(proj.runDVF.cmd.input_f.p)"), proj.max_contig_len)

        do_cmd(proj.runDVF, "DeepVirusFinder", true, parentD; sbatch = sbatch)
    end

    if signal in ["do", "use"]
        if isfile("$(parentD)/$(proj.output_dvf_f.p)") && filesize("$(parentD)/$(proj.output_dvf_f.p)") > 0
            println("\nPost-processing the DVF output file.")
            dvfdf = modif_dvf_out(TableP("$(parentD)/$(proj.output_dvf_f.p)"), proj.scoreTh, proj.pThreshold, split_contigs)
            

            if !isempty(dvfdf)
                #println("Writing the modified DVF output file.")
                CSV.write("$(parentD)/$(proj.modif_output_dvf_f.p)", dvfdf, delim = '\t', header = true)

                #println("Post-processing for aggregation the DVF output files.")
                post_DVF(proj, parentD)
            else
                println("\nDVF did not produce any virus contigs.")
            end
        else
            println("\nDVF did not produce any output files.")
        end
    elseif signal == "use_external"
        if isfile("$(proj.ext_res_D)/$(proj.modif_output_dvf_f.p)") 
            post_DVF(proj, parentD)
        else
            println("\nDVF did not produce any virus contigs.")
        end
    end


    if signal == "do"
        rm_path(proj.todelete)
    end

    return nothing
end

function run_workflow_virSorter2(proj::ProjVirSorter2, signal::String, parentD::String; sbatch::Bool = false)

    if signal == "do"
        do_cmd(proj.virSorter2, "virSorter2", true, parentD; sbatch = sbatch)
    end

    if signal in ["do", "use"]
        prefix = parentD
    elseif signal == "use_external"
        prefix = proj.ext_res_D
    end

    if (isfile("$(prefix)/$(proj.vs2_viral_score_f.p)") && filesize("$(prefix)/$(proj.vs2_viral_score_f.p)") > 0) &&
        (isfile("$(prefix)/$(proj.vs2_viral_boundary_f.p)") && filesize("$(prefix)/$(proj.vs2_viral_boundary_f.p)") > 0) &&
        (isfile("$(prefix)/$(proj.vs2_viral_contigs_f.p)") && filesize("$(prefix)/$(proj.vs2_viral_contigs_f.p)") > 0)

        println("\nPost-processing the virSorter2 output files.")
        post_virSorter2(proj, parentD)
    else
        println("\nVirSorter2 did not produce any output files.")
    end

    return nothing
end

function run_workflow_virSorter(proj::ProjVirSorter, signal::String, parentD::String;)

    if signal == "use_external"
        prefix = proj.ext_res_D
    end

    if (isfile("$(prefix)/$(proj.vs_global_phage_sig.p)") && filesize("$(prefix)/$(proj.vs_global_phage_sig.p)") > 0)

        println("\nPost-processing the virSorter output files.")
        post_virSorter(proj, parentD)
    else
        println("\nVirSorter did not produce any output files.")
    end

    return nothing
end

function run_workflow_vibrant(proj::ProjVibrant, signal::String, parentD::String; sbatch::Bool = false) 
    if signal == "do"
        do_cmd(proj.vibrant, "vibrant", true, parentD; sbatch = sbatch)
    end

    if signal in ["do", "use"]
        prefix = parentD
    elseif signal == "use_external"
        prefix = proj.ext_res_D
    end

    if (isfile("$(prefix)/$(proj.vib_out_nonintegrated_fna.p)") && filesize("$(prefix)/$(proj.vib_out_nonintegrated_fna.p)") > 0) 
        println("\nPost-processing the VIBRANT output files for non-integrated viruses.")
        post_vibrant_nonintegrated!(proj, parentD)
    else
        println("\nVIBRANT did not produce any output files for non-integrated viruses.")
    end

    if (isfile("$(prefix)/$(proj.vib_out_integrated_fna.p)") && filesize("$(prefix)/$(proj.vib_out_integrated_fna.p)") > 0) &&
            (isfile("$(prefix)/$(proj.vib_out_integrated_tsv.p)") && filesize("$(prefix)/$(proj.vib_out_integrated_tsv.p)") > 0)
        
        println("\nPost-processing the VIBRANT output files for integrated viruses.")
        post_vibrant_integrated!(proj, parentD)
    else
        println("\nVIBRANT did not produce any output files for integrated viruses.")
    end

    if signal == "do"
        rm_path(proj.todelete)
    end

    return nothing
end

function run_workflow_viralVerify(proj::ProjViralVerify, signal::String, parentD::String; sbatch::Bool = false)
    if signal == "do"
        do_cmd(proj.viralVerify, "viralVerify", true, parentD; sbatch = sbatch)
    end

    if signal in ["do", "use"]
        prefix = parentD
    elseif signal == "use_external"
        prefix = proj.ext_res_D
    end

    if (isfile("$(prefix)/$(proj.viralVerify_out_p.p)") && filesize("$(prefix)/$(proj.viralVerify_out_p.p)") > 0)
        println("\nPost-processing the viralVerify output file.")
        post_viralVerify(proj, parentD)
    else
        println("\nViralVerify did not produce any output files.")
    end


    return nothing
end

function run_workflow_checkV_NonIntegrated!(proj::ProjCheckVNonintegrated, inref::FnaP, parentD::String; sbatch::Bool = false)
    println("\nAggregate NON-INTEGRATED VIRUSES")
    merged_df = merge_nonintegrated(inref, proj, parentD) 
    out_df = DataFrame()

    if nrow(merged_df) > 0
        do_cmd(proj.checkV, "checkV", true, parentD; sbatch = sbatch)
        println("\nCheckV - postprocessing")
        out_df = postcheckV_nonintegrated!(proj, merged_df, parentD)
    else
        println("There are no contigs left on the NON-INTEGRATED VIRUSES branch!")
    end

    return nothing
end

function run_workflow_checkV_Integrated!(proj::ProjCheckVIntegrated, inref::FnaP, parentD::String; sbatch::Bool = false) 
    println("\nAggregate INTEGRATED VIRUSES")
    integ_all_df = export_aggregated_int(inref, proj, parentD)
    
    do_cmd(proj.checkV1, "checkV for integrated viruses - 1", true, parentD; sbatch = sbatch)
    println("\nCheckV 1 - postprocessing")
    post_checkVint_df = post_checkv1_integrated!(proj, integ_all_df, parentD)
    
    println("\nMerge INTEGRATED VIRUSES")
    merged_int_df = merge_integrated!(inref, proj, post_checkVint_df, parentD)
    
    do_cmd(proj.checkV2, "checkV for integrated viruses - 2", true, parentD; sbatch = sbatch)
    println("CheckV 2 - postprocessing")
    post_checkV2_integrated!(proj, merged_int_df, parentD)

    return nothing
end

function run_workflow_phaTYP_nonintegrated(proj::ProjPhaTYP, parentD::String; sbatch::Bool = false) 
    cd(proj.phatyp.cmd.programD)
    do_cmd(proj.phatyp, "phaTYP", true, parentD; sbatch = sbatch)
    println("PhaTYP - postprocessing")
    merge_postCheckV_phaTYP_nonIntegrated!(proj.phatyp_out_df, proj.indf, proj.mergedPostCheckV_PhaTYP_p, parentD)
    cd("$(parentD)/$(proj.pd)")

    return nothing
end

function run_workflow_genomadTax(proj::ProjGenomadTax, parentD::String, col::Symbol; sbatch::Bool = false)
    
    do_cmd(proj.genomadtax, "geNomad annotate for taxonomy", true, parentD; sbatch = sbatch)
    println("geNomad annotate for taxonomy - postprocessing")

    post_genomadTax(proj, col, parentD)

    return nothing
end

function run_workflow(proj::ProjSViP)

    kw_args = Dict(:sbatch => proj.use_slurm)

    # select contigs with minimum lengthT
    println("
            Length based contig selection ")
    do_wfstep("sel_contig_length", proj, contigLengthSel, (proj.contig_length.min_contig_length, proj.contig_length.inref, FnaP("$(proj.pd)/$(proj.contig_length.outref.p)")), proj.contig_length; logfun = printProjSViP, splatkwargs = (replace_whitespaces = true, ))

    #= detect contig shape || I removed it from here, because with the new algorithm checking till the middle of the sequence and also for ITRs, it takes too long. I've moved it after the viral contigs merging step
    println("
            DETECT CONTIG SHAPE")
    proj.shape_contigs = do_wfstep("shape_contigs", proj, detectContigShapeEnds, (FnaP("$(proj.pd)/$(proj.contig_length.outref.p)"),), proj.shape_contigs; logfun = printProjSViP)
    =# 
    println("\n----------- START INITIAL PREDICTION MODULE")

    # genomad
    println("\n-----------   GENOMAD SUBMODULE")
    do_wfstep("genomad", proj, run_workflow_genomad, (proj.genomad, proj.dosteps["genomad"].signal, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

    println("\n-----------   DEEP VIRUS FINDER SUBMODULE")
    # DVF
    do_wfstep("DVF", proj, run_workflow_DVF, (proj.DVF, proj.dosteps["DVF"].signal, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

    println("\n-----------   VirSorter2 SUBMODULE")
    # virSorter2
    do_wfstep("virSorter2", proj, run_workflow_virSorter2, (proj.virSorter2, proj.dosteps["virSorter2"].signal, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

    println("\n-----------   VirSorter SUBMODULE")
    # virSorter
    do_wfstep("virSorter", proj, run_workflow_virSorter, (proj.virSorter, proj.dosteps["virSorter"].signal, proj.pd); logfun = printProjSViP)

    println("\n-----------   VIBRANT SUBMODULE")
    # vibrant
    do_wfstep("vibrant", proj, run_workflow_vibrant, (proj.vibrant, proj.dosteps["vibrant"].signal, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))
    
    println("\n-----------   ViralVerify SUBMODULE")
    # viralVerify
    do_wfstep("viralVerify", proj, run_workflow_viralVerify, (proj.viralVerify, proj.dosteps["viralVerify"].signal, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

    if proj.stop_after_initial_predictors == false

        println("\n----------- CONSENSUS PREDICTION MODULE")

        #region agregate and checkV NON-INTEGRATED viruses
        if ismissing(proj.checkV_NonIntegrated) == false
            println("\n-----------     Starting the NON-INTEGRATED VIRUSES branch.")

            for i in length(proj.checkV_NonIntegrated.input_dfs_2_aggregate):-1:1
                if isfile("$(proj.pd)/$(proj.checkV_NonIntegrated.input_dfs_2_aggregate[i].p)") == false
                    deleteat!(proj.checkV_NonIntegrated.input_dfs_2_aggregate, i)
                end
            end

            if length(proj.checkV_NonIntegrated.input_dfs_2_aggregate) > 0
                #println(proj.checkV_NonIntegrated.input_dfs_2_aggregate)
                do_wfstep("checkV_NonIntegrated", proj, run_workflow_checkV_NonIntegrated!, (proj.checkV_NonIntegrated, proj.contig_length.outref, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))
                              
                println("\n-----------                 Finished the agregating and checkV steps in the NON-INTEGRATED VIRUSES branch.")
            else
                println("\n            No non-integrated virus contigs to process.")
                println("\n-----------     Finished the NON-INTEGRATED VIRUSES branch.")
            end
        end
        #endregion

        #region agregate and checkV INTEGRATED viruses
        if ismissing(proj.checkV_Integrated) == false
            
            println("\n-----------     Starting the INTEGRATED VIRUSES branch")

            for i in length(proj.checkV_Integrated.input_dfs_2_aggregate):-1:1
                if isfile("$(proj.pd)/$(proj.checkV_Integrated.input_dfs_2_aggregate[i].p)") == false #|| isfile(proj.checkV_Integrated.input_fnas_2_aggregate[i].p) == false
                    deleteat!(proj.checkV_Integrated.input_dfs_2_aggregate, i)
                    #deleteat!(proj.checkV_Integrated.input_fnas_2_aggregate, i)
                end
            end

            if length(proj.checkV_Integrated.input_dfs_2_aggregate) > 0
                do_wfstep("checkV_Integrated", proj, run_workflow_checkV_Integrated!, (proj.checkV_Integrated, proj.contig_length.outref, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))
                println("\n-----------     Finished the merging of overlaping viral sequences and checkV steps in the INTEGRATED VIRUSES branch")
            else
                println("\n            No integrated virus contigs to process.")
            end
        end
        #endregion

        #region unMix viruses
        if isfile("$(proj.pd)/$(proj.checkV_Integrated.postcheckV2_integrated_df_p.p)") || isfile("$(proj.pd)/$(proj.checkV_NonIntegrated.postcheckV_nonintegrated_df_p.p)")
            println("\n-----------     Starting unMixing viruses branch")

            do_wfstep("detect_mixed_viruses", proj, detect_mixed_virs, (proj.detect_mixed_viruses, proj.pd); logfun = printProjSViP)
            println("\n-----------     Finished unMixing viruses branch")
        end
        #endregion


        #region PhaTYP and Final Thresholding nonIntegrated viruses
        # check for the existance of output files for the nonintegrated viruses in the MIXED folder
        if isfile("$(proj.pd)/$(proj.detect_mixed_viruses.outDf_resolved_nonInt_p.p)")
            println("\n-----------     Starting the rest of the NON-INTEGRATED VIRUSES branch.")
            do_wfstep("phaTYP_nonintegrated", proj, run_workflow_phaTYP_nonintegrated, (proj.phaTYP_nonintegrated, proj.pd); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))
            
            do_wfstep("genomadTax_NonInt", proj, run_workflow_genomadTax, (proj.genomadTax_NonInt, proj.pd, :virus_name); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

            println("\nFinal Thresholding")
            do_wfstep("final_thresholding_NonIntegrated", proj, apply_thresholds!, (proj.final_thresholding_NonIntegrated, :virus_name, proj.pd, order_NonIntDf!, proj.sampleName, proj.sampleSet); logfun = printProjSViP)

            println("\n-----------     Finished NON-INTEGRATED VIRUSES branch.")
        else
            println("\n There are no viruses left on the  NON-INTEGRATED VIRUSES Branch")
        end
        #endregion

        #region Final Thresholding integrated viruses
        if isfile("$(proj.pd)/$(proj.detect_mixed_viruses.outDf_resolved_Int_p.p)")
            println("\n-----------     Starting the rest of the INTEGRATED VIRUSES branch.")
            do_wfstep("genomadTax_Int", proj, run_workflow_genomadTax, (proj.genomadTax_Int, proj.pd, :provirus_name); logfun = printProjSViP, splatkwargs = (sbatch = proj.use_slurm,))

            println("\n Final Thresholding")                
            do_wfstep("final_thresholding_Integrated", proj, apply_thresholds!, (proj.final_thresholding_Integrated, :provirus_name, proj.pd, order_IntDf!, proj.sampleName, proj.sampleSet); logfun = printProjSViP)
            
            println("\n -----------     Finished the INTEGRATED VIRUSES branch.")
        else
            println("\n There are no viruses left on the  INTEGRATED VIRUSES Branch")
        end
        #end region

    end

    return nothing
end

function run_workflowMDoViP(mproj::ProjMultiWorkflow)
    for i in eachindex(mproj.allSingleWorkflows)
        if mproj.dosteps[i].progress in ["not_done", "running"]
            println("\n            ---------------------- Starting TO RUN the DoViP workflow for input file $i -------------------------")
            set2running!(i, mproj)

            run_workflow(mproj.allSingleWorkflows[i])

            set2finished!(i, mproj)

            println("\n            ---------------------- Finished RUNNING the DoViP workflow for input file $i -------------------------")
        else
            println("\n            ---------------------- Skipping TO RUN the DoViP workflow for input file $i, because its status in a previous run was $(mproj.dosteps[i].progress) -------------------------")
        end
    end

    return nothing
end