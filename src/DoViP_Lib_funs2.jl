
function printProjSViP(outp::String, proj::ProjSViP)

    input_nonint_table_files = join([proj.checkV_NonIntegrated.input_dfs_2_aggregate[i].p for i in eachindex(proj.checkV_NonIntegrated.input_dfs_2_aggregate)], "\n     ")
    input_int_table_files = join([proj.checkV_Integrated.input_dfs_2_aggregate[i].p for i in eachindex(proj.checkV_Integrated.input_dfs_2_aggregate)], "\n  ")

    #region genomad
    if ismissing(proj.genomad) == false
        if isnothing(proj.genomad.genomad) == false
            genomad_cmd = "        task: $(proj.genomad.genomad.cmd.task)
            input file: $(proj.genomad.genomad.cmd.input_f.p)
            output folder: $(proj.genomad.genomad.cmd.output_d)
            database: $(proj.genomad.genomad.cmd.database)
            minimum score: $(proj.genomad.genomad.cmd.min_score)
            conda environment: $(proj.genomad.genomad.env) "

            genomad_sbatch = "        sbatch --time: $(proj.genomad.genomad.sbatch_maxtime)
            sbatch --cpus-per-task: $(proj.genomad.genomad.sbatch_cpuspertask)
            sbatch --mem: $(proj.genomad.genomad.sbatch_mem)"
        else
            genomad_cmd = "     unknown, due to the use of external results"
            genomad_sbatch = "      unknown, due to the use of external results"
        end

        genomad_params = "subfolder: $(proj.genomad.pd)
    use external results: $(proj.genomad.ext_res)
    folder with external results: $(proj.genomad.ext_res_D)
    genomad command parameters:
    $(genomad_cmd)
    slurm options for genomad command:    
    $(genomad_sbatch)
    genomad main output table: $(proj.genomad.genomad_out_table_p.p)
    genomad main fasta file: $(proj.genomad.genomad_out_fnap.p)

    genomad outputs post-processed by DoViP:
        table file for NON-INTEGRATED viruses: $(proj.genomad.postgenomad_nonintegrated_df.p)
        table file for INTEGRATED viruses: $(proj.genomad.postgenomad_integrated_df.p)"

    else
        genomad_params = "  missing"
    end
    #endregion

    #region DVD
    if ismissing(proj.DVF.runDVF) == false
        DVF_params = "subfolder: $(proj.DVF.pd)
    input file for this step (to be pre-processed by DoViP before DVF): $(proj.DVF.input_f.p)
    maximum contig length allowed during pre-processing: $(proj.DVF.max_contig_len)
    DVF command parameters: 
        program path: $(proj.DVF.runDVF.cmd.program)
        input fasta: $(proj.DVF.runDVF.cmd.input_f.p)
        output folder: $(proj.DVF.runDVF.cmd.output_d)
        minimum contig length: $(proj.DVF.runDVF.cmd.min_contig_len)
        number of threads: $(proj.DVF.runDVF.cmd.num_threads)
        conda environment: $(proj.DVF.runDVF.env)
    slurm options for DVF command:
        sbatch --time: $(proj.DVF.runDVF.sbatch_maxtime)
        sbatch --cpus-per-task: $(proj.DVF.runDVF.sbatch_cpuspertask)
        sbatch --mem: $(proj.DVF.runDVF.sbatch_mem)
    main DFV output file: $(proj.DVF.output_dvf_f.p)
    post-DFV results thresholding:
        score threshold: $(proj.DVF.scoreTh)
        p-value threshold: $(proj.DVF.pThreshold)
    output file after thresholding: $(proj.DVF.modif_output_dvf_f.p)

    final DVF outputs after further DoViP postprocessing:
        table file: $(proj.DVF.postdvf_nonintegrated_df.p)"
    else
        DVF_params = "  missing"
    end
    #endregion

    #region virSorter2
    if ismissing(proj.virSorter2) == false
        if isnothing(proj.virSorter2.virSorter2) == false
            virSorter2_cmd = "        input file: $(proj.virSorter2.virSorter2.cmd.input_f.p)
            output folder: $(proj.virSorter2.virSorter2.cmd.output_d)
            database path: $(proj.virSorter2.virSorter2.cmd.db_p)
            high_confidence_only: $(proj.virSorter2.virSorter2.cmd.high_confidence_only)
            number of threads: $(proj.virSorter2.virSorter2.cmd.num_threads)
            minimum score: $(proj.virSorter2.virSorter2.cmd.min_score)
            minimum contig length: $(proj.virSorter2.virSorter2.cmd.min_length)
            conda environment: $(proj.virSorter2.virSorter2.env)"

            virSorter2_sbatch = "        sbatch --time: $(proj.virSorter2.virSorter2.sbatch_maxtime)
            sbatch --cpus-per-task: $(proj.virSorter2.virSorter2.sbatch_cpuspertask)
            sbatch --mem: $(proj.virSorter2.virSorter2.sbatch_mem)"
        else
            virSorter2_cmd = "      external results"
            virSorter2_sbatch = "       external results"
        end

        virSorter2_params = "subfolder: $(proj.virSorter2.pd)
    use external results: $(proj.virSorter2.ext_res)
    folder with external results: $(proj.virSorter2.ext_res_D)
    VirSorter2 command parameters:
    $(virSorter2_cmd)
    slurm options for VirSorter2 command:
    $(virSorter2_sbatch)
    main VirSorter2 outputs:
        table with viral scores: $(proj.virSorter2.vs2_viral_score_f.p) 
        table with virus boundary info: $(proj.virSorter2.vs2_viral_boundary_f.p)
    
    VirSorter2 outputs post-processed by DoViP:
        table file for NON-INTEGRATED viruses: $(proj.virSorter2.postvs2_nonintegrated_df.p)
        table file for INTEGRATED viruses: $(proj.virSorter2.postvs2_integrated_df.p)"
    else
        virSorter2_params = "   missing"
    end
    #endregion

    #region VIBRANT
    if ismissing(proj.vibrant) == false
        if isnothing(proj.vibrant.vibrant) == false
            vibrant_cmd = "        program path: $(proj.vibrant.vibrant.cmd.program)
            input file: $(proj.vibrant.vibrant.cmd.input_f.p)
            output folder: $(proj.vibrant.vibrant.cmd.output_d)
            database path: $(proj.vibrant.vibrant.cmd.database)
            number of threads: $(proj.vibrant.vibrant.cmd.num_threads)
            minimum contig length: $(proj.vibrant.vibrant.cmd.min_length)
            conda environment: $(proj.vibrant.vibrant.env)"

            vibrant_sbatch = "        sbatch --time: $(proj.vibrant.vibrant.sbatch_maxtime)
            sbatch --cpus-per-task: $(proj.vibrant.vibrant.sbatch_cpuspertask)
            sbatch --mem: $(proj.vibrant.vibrant.sbatch_mem)"
        else
            vibrant_cmd = "     external results"
            vibrant_sbatch = "      external results"
        end

        vibrant_params = "subfolder: $(proj.vibrant.pd)
    use external results: $(proj.vibrant.ext_res)
    folder with external results: $(proj.vibrant.ext_res_D)
    parameters for VIBRANT command:
    $(vibrant_cmd)
    slurm options for VIBRANT command:
    $(vibrant_sbatch)
    main VIBRANT outputs:
        table file for integrated viruses: $(proj.vibrant.vib_out_integrated_tsv.p)
        fasta file for integrated viruses: $(proj.vibrant.vib_out_integrated_fna.p)
        fasta file for non-integrated viruses: $(proj.vibrant.vib_out_nonintegrated_fna.p)
    
    VINRANT outputs post-processed by DoViP:
        table file for NON-INTEGRATED viruses: $(proj.vibrant.postvib_nonintegrated_df.p)
        table file for INTEGRATED viruses: $(proj.vibrant.postvib_integrated_df.p)"
    else
        vibrant_params = "  missing"
    end
    #endregion

    #region ViralVerify
    if ismissing(proj.viralVerify) == false
        if isnothing(proj.viralVerify.viralVerify) == false
            viralVerfify_cmd = "        program path: $(proj.viralVerify.viralVerify.cmd.program)
            input file: $(proj.viralVerify.viralVerify.cmd.input_f.p)
            output folder: $(proj.viralVerify.viralVerify.cmd.output_d)
            database path: $(proj.viralVerify.viralVerify.cmd.database)
            viralVerify threshold: $(proj.viralVerify.viralVerify.cmd.threshold)
            number of threads: $(proj.viralVerify.viralVerify.cmd.num_threads)
            conda environment: $(proj.viralVerify.viralVerify.env)"

            viralVerify_sbatch = "        sbatch --time: $(proj.viralVerify.viralVerify.sbatch_maxtime)
            sbatch --cpus-per-task: $(proj.viralVerify.viralVerify.sbatch_cpuspertask)
            sbatch --mem: $(proj.viralVerify.viralVerify.sbatch_mem)"
        else
            viralVerfify_cmd = "        external results"
            viralVerify_sbatch = "      external results"
        end

        viralVerify_params = "subfolder: $(proj.viralVerify.pd)
    use external results: $(proj.viralVerify.ext_res)
    folder with external results: $(proj.viralVerify.ext_res_D)
    parameters for ViralVerify command:
    $(viralVerfify_cmd)
    slurm options for ViraVerify command:
    $(viralVerify_sbatch)
    main ViralVerify outputs:
        table file: $(proj.viralVerify.viralVerify_out_p.p)
    
    ViralVerify outputs post-processed by DoViP:
        table file: $(proj.viralVerify.postViralVerify_df.p)"
    else
        viralVerify_params = "  missing"
    end
    #endregion

    #region CheckV-NonIntegrated
   if ismissing(proj.checkV_NonIntegrated) == false
        checkV_nonInt_params = "subfolder: $(proj.checkV_NonIntegrated.pd)
    input table files from the initial predictors to agregate and use further: $(input_nonint_table_files)
    table for virus contigs agregated from the initial predictors: $(proj.checkV_NonIntegrated.output_aggregated_df.p)
    fasta file with virus contigs agregated from the initial predictors: $(proj.checkV_NonIntegrated.output_aggregated_fna.p)
    parameters for the CheckV command:
        input file: $(proj.checkV_NonIntegrated.checkV.cmd.input_f.p)
        output folder: $(proj.checkV_NonIntegrated.checkV.cmd.output_d)
        database path: $(proj.checkV_NonIntegrated.checkV.cmd.database)
        num_threads: $(proj.checkV_NonIntegrated.checkV.cmd.num_threads)
        conda environment: $(proj.checkV_NonIntegrated.checkV.env)
    slurm options for CheckV command:
        sbatch --time: $(proj.checkV_NonIntegrated.checkV.sbatch_maxtime)
        sbatch --cpus-per-task: $(proj.checkV_NonIntegrated.checkV.sbatch_cpuspertask)
        sbatch --mem: $(proj.checkV_NonIntegrated.checkV.sbatch_mem)
    main CheckV outputs:
        summary table file: $(proj.checkV_NonIntegrated.checkV_out_nonintegrated_summary_df.p)
        complete table file: $(proj.checkV_NonIntegrated.checkV_out_nonintegrated_complete_df.p)
        fasta file for proviruses: $(proj.checkV_NonIntegrated.checkV_out_provir_fna.p)
    
    CheckV outputs post-processed by DoViP:
        table file for NON-INTEGRATED viruses: $(proj.checkV_NonIntegrated.postcheckV_nonintegrated_df_p.p)
        fasta file for NON-INTEGRATED viruses: $(proj.checkV_NonIntegrated.postcheckV_nonintegrated_fna.p)
        fasta file for NON-INTEGRATED viruses, with trimmed DTRs on the right contig side: $(proj.checkV_NonIntegrated.postcheckV_nonintegrated_fna_trimmed_DTR.p)"
   else
    checkV_nonInt_params = "    missing"
   end
    #endregion
    
    #region phaTYP-NonIntegrated
    if ismissing(proj.phaTYP_nonintegrated) == false
        phaTYP_params = "subfolder: $(proj.phaTYP_nonintegrated.pd)
    parameters for PhaTYP command:
        program path: $(proj.phaTYP_nonintegrated.phatyp.cmd.programD)
        input file: $(proj.phaTYP_nonintegrated.phatyp.cmd.input_f.p)
        database path: $(proj.phaTYP_nonintegrated.phatyp.cmd.database_d)
        parameters folder: $(proj.phaTYP_nonintegrated.phatyp.cmd.parameters_d)
        folder for temporary outputs: $(proj.phaTYP_nonintegrated.phatyp.cmd.outputtemp_d)
        output folder: $(proj.phaTYP_nonintegrated.phatyp.cmd.output_d)
        minimum contig length: $(proj.phaTYP_nonintegrated.phatyp.cmd.min_len)
        number of threads: $(proj.phaTYP_nonintegrated.phatyp.cmd.num_threads)
        conda environment: $(proj.phaTYP_nonintegrated.phatyp.env) 
    main outputs from PHATYP:
        table file: $(proj.phaTYP_nonintegrated.phatyp_out_df.p)
    
    Joined results from PhaTYP and CheckV for NON-INTEGRATED viruses:
        table file: $(proj.phaTYP_nonintegrated.mergedPostCheckV_PhaTYP_p.p)"
    else
        phaTYP_params = "   missing"
    end
    #endregion

    #region FinalThresholding Non-integrated
    if ismissing(proj.final_thresholding_NonIntegrated) == false
        final_th_nonint_params = "subfolder: $(proj.final_thresholding_NonIntegrated.pd)
    input fasta file: $(proj.final_thresholding_NonIntegrated.inFna.p)
    input fasta file with contigs, with trimmed DTRs on the right contig side: $(proj.final_thresholding_NonIntegrated.inFna_trimmed_DTR.p)
    input table file: $(proj.final_thresholding_NonIntegrated.inTsv.p)
    initial predictors: $(proj.final_thresholding_NonIntegrated.predictors)
    selection thresholds for viral contigs:
        - minimum number of predictors required if CheckV completeness is NA: $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_NA)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'AAIighConf': $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_AAIHighConf)
        - minimum completeness required when the method for determination of CheckV completeness is 'AAIighConf': $(proj.final_thresholding_NonIntegrated.th_completeness_CheckV_AAIHighConf)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'AAIMediumConf': $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_AAIMediumConf)
        - minimum completeness required when the method for determination of CheckV completeness is 'AAIMediumConf': $(proj.final_thresholding_NonIntegrated.th_completeness_CheckV_AAIMediumConf)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'CheckV_HMM': $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_HMM)
        - minimum completeness required when the method for determination of CheckV completeness is 'CheckV_HMM': $(proj.final_thresholding_NonIntegrated.th_completeness_CheckV_HMM)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'CheckV_DTR_ITR_AAI': $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_DTR_ITR_AAI)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'CheckV_DTR_ITR_HMM': $(proj.final_thresholding_NonIntegrated.th_num_predictors_CheckV_DTR_ITR_HMM)
    
    final outputs:
        - table file with selected viral contigs, their initial predictors, CheckV statistics, taxonomy and life style predictions: $(proj.final_thresholding_NonIntegrated.outTsv.p)
        - fasta file with selected viral contigs: $(proj.final_thresholding_NonIntegrated.outFnaP.p)
        - fasta file with selected viral contigs, with trimmed DTRs on the right contig side: $(proj.final_thresholding_NonIntegrated.outFna_trimmed_DTR.p)"
    else
        final_th_nonint_params = "  missing"
    end
    #endregion

    #region CheckV-Integrated
    if ismissing(proj.checkV_Integrated) == false
        checkV_int_params = "subfolder: $(proj.checkV_Integrated.pd)
    input table files from the initial predictors to agregate and use further: $(input_int_table_files)

    table for virus contigs agregated from the initial predictors: $(proj.checkV_Integrated.output_aggregated_df.p)
    fasta file with virus contigs containing potentially INTEGRATED virus region: $(proj.checkV_Integrated.output_aggregated_int_contigs_fna.p)
    fasta file with INTEGRATED viral regions agregated from the initial predictors: $(proj.checkV_Integrated.output_aggregated_fna.p)
    parameters for first CheckV (CheckV1) command:
        input file: $(proj.checkV_Integrated.checkV1.cmd.input_f.p)
        output folder: $(proj.checkV_Integrated.checkV1.cmd.output_d)
        database path: $(proj.checkV_Integrated.checkV1.cmd.database)
        num_threads: $(proj.checkV_Integrated.checkV1.cmd.num_threads)
        conda environment: $(proj.checkV_Integrated.checkV1.env)
    main CheckV1 outputs:
        fasta file for proviruses: $(proj.checkV_Integrated.checkV1_out_provir_fna.p)
    
    Outputs from the merging step of the overlaping integrated viruses (DoViP post-CheckV1 steps):
        intermediary files:
            table file with corected coordinates for agregated virus contigs: $(proj.checkV_Integrated.postcheckv1_integrated_cor_df.p)
            table file with corected coordinates and merged IDs for agregated virus contigs: $(proj.checkV_Integrated.postcheckv1_integrated_cor_withmergedprovirIDs_df.p)
        final files for merged integrated viruses:
            table file for merged viruses: $(proj.checkV_Integrated.merged_integrated_DF.p)
            fasta file for merged viruses: $(proj.checkV_Integrated.merged_integrated_fna.p)

    parameters for second CheckV (CheckV2) command:
        input file: $(proj.checkV_Integrated.checkV2.cmd.input_f.p)
        output folder: $(proj.checkV_Integrated.checkV2.cmd.output_d)
        database path: $(proj.checkV_Integrated.checkV2.cmd.database)
        num_threads: $(proj.checkV_Integrated.checkV2.cmd.num_threads)
        conda environment: $(proj.checkV_Integrated.checkV2.env)
    main CheckV2 outputs:
        table file with CheckV summary: $(proj.checkV_Integrated.checkV2_out_integrated_df.p)

    slurm options for CheckV commands:
        sbatch --time: $(proj.checkV_Integrated.checkV1.sbatch_maxtime)
        sbatch --cpus-per-task: $(proj.checkV_Integrated.checkV1.sbatch_cpuspertask)
        sbatch --mem: $(proj.checkV_Integrated.checkV1.sbatch_mem)

    Final report file from this step:
        table file (including merged overlaping viruses and their CheckV statistics): $(proj.checkV_Integrated.postcheckV2_integrated_df_p.p)"
    else
        checkV_int_params = "   missing"
    end
    #endregion

    #region FinalThresholding Integrated
    if ismissing(proj.final_thresholding_Integrated) == false
        final_th_int_params = "subfolder: $(proj.final_thresholding_Integrated.pd)
    input fasta file: $(proj.final_thresholding_Integrated.inFna.p)
    input table file: $(proj.final_thresholding_Integrated.inTsv.p)
    initial predictors: $(proj.final_thresholding_Integrated.predictors)
    selection thresholds for viral contigs:
        - minimum number of predictors required if CheckV completeness is NA: $(proj.final_thresholding_Integrated.th_num_predictors_CheckV_NA)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'AAIighConf': $(proj.final_thresholding_Integrated.th_num_predictors_CheckV_AAIHighConf)
        - minimum completeness required when the method for determination of CheckV completeness is 'AAIighConf': $(proj.final_thresholding_Integrated.th_completeness_CheckV_AAIHighConf)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'AAIMediumConf': $(proj.final_thresholding_Integrated.th_num_predictors_CheckV_AAIMediumConf)
        - minimum completeness required when the method for determination of CheckV completeness is 'AAIMediumConf': $(proj.final_thresholding_Integrated.th_completeness_CheckV_AAIMediumConf)
        - minimum number of predictors required if the method for determination of CheckV completeness is 'CheckV_HMM': $(proj.final_thresholding_Integrated.th_num_predictors_CheckV_HMM)
        - minimum completeness required when the method for determination of CheckV completeness is 'CheckV_HMM': $(proj.final_thresholding_Integrated.th_completeness_CheckV_HMM)
    final outputs:
        - table file with selected viral contigs, their initial predictors, CheckV statistics, taxonomy and life style predictions: $(proj.final_thresholding_Integrated.outTsv.p)
        - fasta file with selected viral contigs: $(proj.final_thresholding_Integrated.outFnaP.p)
        "
    else
        final_th_int_params = " missing"
    end
    #endregion

    #region DetectMixedViruses
    if ismissing(proj.detect_mixed_viruses) == false
        detect_mixed_viruses_params = "subfolder: $(proj.detect_mixed_viruses.pd)
    input table file with integrated viruses: $(proj.detect_mixed_viruses.inDf_Int.p)
    input table file with non-integrated viruses: $(proj.detect_mixed_viruses.inDf_NonInt.p)
    output table file with integrated viruses found also as non-integrated: $(proj.detect_mixed_viruses.outDf_Int.p)
    output table file with non-integrated viruses found also as integrated: $(proj.detect_mixed_viruses.outDf_NonInt.p)"
    else
        detect_mixed_viruses_params = " missing"
    end
    #endregion

    toprint = "
============================================================== DoViP parameters and status ========================================================

PROJECT TYPE: $(proj.projtype)
PROJECT PARENT FOLDER: $(proj.pd)

INPUT FASTA FILE: $(proj.inref.p)
SAMPLE NAME: $(proj.sampleName)

PROJECT FULL FOLDER: $(proj.pd)/$(proj.sampleName)

PERFORM CALCULATIONS USING SLURM: $(proj.use_slurm)
CONTINUE PREVIOUS PROJECT: $(proj.continue_project)
STOP AFTER INITIAL PREDICTORS STEP: $(proj.stop_after_initial_predictors)  

STATUS REPORT of individual steps
    select contig length:                signal is $(proj.dosteps["sel_contig_length"].signal), progress is $(proj.dosteps["sel_contig_length"].progress)
    
    => Initial prediction steps

    geNomad:                             signal is $(proj.dosteps["genomad"].signal), progress is $(proj.dosteps["genomad"].progress)
    DVF:                                 signal is $(proj.dosteps["DVF"].signal), progress is $(proj.dosteps["DVF"].progress)
    VirSorter2:                          signal is $(proj.dosteps["virSorter2"].signal), progress is $(proj.dosteps["virSorter2"].progress)
    VIBRANT:                             signal is $(proj.dosteps["vibrant"].signal), progress is $(proj.dosteps["vibrant"].progress)
    ViralVerify:                         signal is $(proj.dosteps["viralVerify"].signal), progress is $(proj.dosteps["viralVerify"].progress)

    => Consensus prediction steps

    NON-INTEGRATED VIRUSES BRANCH
    checkV_NonIntegrated:                signal is $(proj.dosteps["checkV_NonIntegrated"].signal), progress is $(proj.dosteps["checkV_NonIntegrated"].progress)
    PhaTYP_nonintegrated:                signal is $(proj.dosteps["phaTYP_nonintegrated"].signal), progress is $(proj.dosteps["phaTYP_nonintegrated"].progress)
    Final_thresholding_NonIntegrated:    signal is $(proj.dosteps["final_thresholding_NonIntegrated"].signal), progress is $(proj.dosteps["final_thresholding_NonIntegrated"].progress)
        
    INTEGRATED VIRUSES BRANCH
    checkV_Integrated:                   signal is $(proj.dosteps["checkV_Integrated"].signal), progress is $(proj.dosteps["checkV_Integrated"].progress)
    Final_thresholding_Integrated:       signal is $(proj.dosteps["final_thresholding_Integrated"].signal), progress is $(proj.dosteps["final_thresholding_Integrated"].progress)
        
    MIXED VIRUSES BRANCH
    Detect_mixed_viruses:                signal is $(proj.dosteps["detect_mixed_viruses"].signal), progress is $(proj.dosteps["detect_mixed_viruses"].progress)
    
   

PARAMETERS FOR 'SELECT CONTIG LENGTH'
    subfolder: $(proj.contig_length.pd)
    minimum contig length: $(proj.contig_length.min_contig_length)
    input fasta file: $(proj.contig_length.inref.p)
    output fasta file: $(proj.contig_length.outref.p)

==========> INITIAL PREDICTION STEPS:

PARAMETERS FOR geNomad
    $(genomad_params)

PARAMETERS FOR DeepVirFinder
    $(DVF_params)

PARAMETERS FOR VirSorter2
    $(virSorter2_params)

PARAMETERS FOR VIBRANT
    $(vibrant_params)

PARAMETERS for ViraVerify
    $(viralVerify_params)      

==========> CONSENSUS PREDICTION STEPS

=> NON-INTEGRATED VIRUSES BRANCH

PARAMETERS FOR CheckV processing steps
    $(checkV_nonInt_params)


PARAMETERS FOR PHATYP
    $(phaTYP_params)

PARAMETERS FOR FINAL SELECTION OF NON-INTEGRATED VIRUSES
    $(final_th_nonint_params)

=> INTEGRATED VIRUSES BRANCH

PARAMETERS FOR CheckV processing steps
    $(checkV_int_params)

PARAMETERS FOR FINAL SELECTION OF INTEGRATED VIRUSES
    $(final_th_int_params)


=> MIXED (UNRESOLVED IF INTEGRATED OR NON-INTEGRATED) VIRUSES BRANCH

PARAMETERS FOR DETECTING MIXED VIRUSES
    $(detect_mixed_viruses_params)

============================================================== DoViP parameters and status ======================================================== 
"

    open(outp, "w") do file
        write(file, toprint)
    end    

    return nothing
end


function order_NonIntDf!(df::DataFrame, sample_name::String, sample_set::String)
    df[!, :sample_name] = fill(sample_name, nrow(df))
    df[!, :sample_set] = fill(sample_set, nrow(df))
    
    column_names = names(df)
    column_order = [:sample_set, :sample_name, 
                    :contig_name, :contig_shape, :contig_shape_remarks, :contig_end_repeat_type, :contig_end_repeat_size] 

    if "predictor_genomad" in column_names
        push!(column_order, :topology_genomad)
    end

    if "predictor_virSorter2" in column_names
        push!(column_order, :shape_virSorter2)
    end

    append!(column_order, [:contig_trimmed_end, :contig_full_end, 
                            :virus_name, :virus_length, :virus_start, :virus_end])

    if "predictor_genomad" in column_names
        push!(column_order, :predictor_genomad)
    end

    if "predictor_dvf" in column_names
        push!(column_order, :predictor_dvf)
    end

    if "predictor_virSorter2" in column_names
        push!(column_order, :predictor_virSorter2)
    end

    if "predictor_vibrant" in column_names
        push!(column_order, :predictor_vibrant)
    end

    if "predictor_viralVerify" in column_names
        push!(column_order, :predictor_viralVerify)
    end

    if "predictor_genomad" in column_names
        push!(column_order, :virus_score_genomad)
    end

    if "predictor_dvf" in column_names
        append!(column_order, [:score_dvf, :pvalue_dvf])
    end

    if "prediction_long_contig_dvf" in column_names
        push!(column_order, :prediction_long_contig_dvf)
    end

    if "predictor_virSorter2" in column_names
        append!(column_order, [:max_score_virSorter2, :hallmark_virSorter2])
    end
    
    if "predictor_viralVerify" in column_names
        append!(column_order, [:score_viralVerify, :shape_viralVerify])
    end

    append!(column_order, [:predictors_total, :completeness_checkV, :completeness_method_checkV, :contamination_checkV, :kmer_freq_checkV, :warnings_checkV, :gene_count_checkV, :viral_genes_checkV, :host_genes_checkV, :checkv_quality_checkV, :miuvig_quality_checkV])
                            
    if "predictor_genomad" in column_names
        push!(column_order, :taxonomy_genomad)
    end

    if "predictor_virSorter2" in column_names
        push!(column_order, :max_score_group_virSorter2)
    end
    
    append!(column_order, [:virus_type_DoViP, :prediction_PhaTYP, :score_PhaTYP])

    df = select!(df, column_order)

    return df
end


function order_IntDf!(df::DataFrame, sample_name::String, sample_set::String)
    df[!, :sample_name] = fill(sample_name, nrow(df))
    df[!, :sample_set] = fill(sample_set, nrow(df))
    
    df = rename!(df, :length => :provirus_length, :taxonomy => :taxonomy_genomad)
    column_names = names(df)

    df[!, :proviral_fraction] = Vector{Union{Missing, Float64}}(missing, nrow(df))
    for i in 1:nrow(df)
        df[i, :proviral_fraction] = (df[i, :provirus_length]*100 / df[i, :contig_trimmed_end])
    end

    column_order = [:sample_set, :sample_name, 
                    :contig_name, :contig_shape, :contig_shape_remarks, :contig_end_repeat_type, :contig_end_repeat_size] 
    
    if "predictor_genomad" in column_names
        append!(column_order, [:topology_genomad])
    end

    if "predictor_virSorter2" in column_names
        push!(column_order, :shape_virSorter2)
    end

    append!(column_order, [:contig_trimmed_end, :contig_full_end, 
                            :provirus_name, :provirus_length, :provirus_start, :provirus_end, :r_provirus_start, :r_provirus_end, :proviral_fraction])

    if "predictor_genomad" in column_names
        push!(column_order, :predictor_genomad)
    end

    if "predictor_dvf" in column_names
        push!(column_order, :predictor_dvf)
    end

    if "predictor_virSorter2" in column_names
        push!(column_order, :predictor_virSorter2)
    end

    if "predictor_vibrant" in column_names
        push!(column_order, :predictor_vibrant)
    end

    if "predictor_viralVerify" in column_names
        push!(column_order, :predictor_viralVerify)
    end

    append!(column_order, [:predictors_total, :completeness_checkV, :completeness_method_checkV, :contamination_checkV, :kmer_freq_checkV, :warnings_checkV, :gene_count_checkV, :viral_genes_checkV, :host_genes_checkV, :checkv_quality_checkV, :miuvig_quality_checkV])
    
    if "predictor_genomad" in column_names
        push!(column_order, :taxonomy_genomad)
    end

    push!(column_order, :virus_type_DoViP)

    df = select!(df, column_order)

    return df
end