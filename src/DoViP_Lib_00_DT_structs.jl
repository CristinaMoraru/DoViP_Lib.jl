
struct Proj_sel_contig_length <: BioinfProj
    pd::String
    min_contig_length::Int64
    inref::FnaP
    outref::FnaP
end


struct ProjGenomad <: BioinfProj
    pd::String
    ext_res::Bool
    ext_res_D::Union{Nothing, String}
    genomad::Union{Nothing, WrapCmd{RunGenomadCmd}}
    genomad_out_table_p::TableP
    genomad_out_fnap::FnaP
    postgenomad_nonintegrated_df::TableP
    postgenomad_integrated_df::TableP
    #postgenomad_integrated_fnaf::FnaP
    todelete::Vector{String}
end

struct ProjDVF <: BioinfProj
    pd::String
    input_f::FnaP
    max_contig_len::Int64
    runDVF::WrapCmd{RunDVFCmd}
    scoreTh::Float64
    pThreshold::Float64
    output_dvf_f::TableP
    modif_output_dvf_f::TableP
    postdfv_nonintegrated_fna::FnaP
    postdvf_nonintegrated_df::TableP
    todelete::Vector{String}
end

struct ProjVirSorter2 <: BioinfProj
    pd::String
    ext_res::Bool
    ext_res_D::Union{Nothing, String}
    virSorter2::Union{Nothing, WrapCmd{RunVirSorterCmd}}
    vs2_viral_score_f::TableP
    vs2_viral_boundary_f::TableP
    vs2_viral_contigs_f::FnaP
    postvs2_nonintegrated_df::TableP
    postvs2_integrated_df::TableP
end

struct ProjVibrant <: BioinfProj
    pd::String
    ext_res::Bool
    ext_res_D::Union{Nothing, String}
    vibrant::Union{Nothing, WrapCmd{RunVibrantCmd}}
    vib_out_integrated_tsv::TableP
    vib_out_integrated_fna::FnaP
    vib_out_nonintegrated_fna::FnaP
    postvib_nonintegrated_df::TableP
    postvib_integrated_df::TableP
    todelete::Vector{String}
end

struct ProjViralVerify <: BioinfProj
    pd::String
    ext_res::Bool
    ext_res_D::Union{Nothing, String}
    viralVerify::Union{Nothing, WrapCmd{RunViralVerify}}
    viralVerify_out_p::TableP
    score_threshold::Float64
    postViralVerify_df::TableP
end


#=struct ProjCenoteTaker3 <: BioinfProj
    pd::String
    ct3::WrapCmd{RunCenoteTaker3}
    ct3_out
end =#

Base.@kwdef mutable struct ProjCheckVNonintegrated <: BioinfProj
    pd::String
    input_dfs_2_aggregate::Vector{TableP}
    output_aggregated_df::TableP
    output_aggregated_fna::FnaP
    checkV::WrapCmd{RunCheckVCmd}
    checkV_out_nonintegrated_summary_df::TableP
    checkV_out_nonintegrated_complete_df::TableP
    checkV_out_provir_fna::FnaP
    #postcheckV_nonintegrated_df::Union{Missing, DataFrame} = missing
    postcheckV_nonintegrated_df_p::TableP
    postcheckV_nonintegrated_fna::FnaP
    postcheckV_nonintegrated_fna_trimmed_DTR::FnaP
    postcheckV_integrated_df_p::TableP
end

Base.@kwdef mutable struct ProjCheckVIntegrated <: BioinfProj
    pd::String
    input_dfs_2_aggregate::Vector{TableP}
    #input_fnas_2_aggregate::Vector{FnaP}
    output_aggregated_int_contigs_fna::FnaP
    output_aggregated_df::TableP
    output_aggregated_fna::FnaP
    checkV1::WrapCmd{RunCheckVCmd}
    checkV1_out_provir_fna::FnaP
    postcheckv1_integrated_cor_df::TableP
    postcheckv1_integrated_cor_withmergedprovirIDs_df::TableP
    #integrated_cor_fnaf::FnaP
    predictors::Vector{Symbol}
    merged_integrated_DF::TableP
    merged_integrated_fna::FnaP
    checkV2::WrapCmd{RunCheckVCmd}
    checkV2_out_integrated_summary_df::TableP
    checkV2_out_integrated_complete_df::TableP
    #postcheckV2_integrated_df::Union{Missing, DataFrame} = missing
    postcheckV2_integrated_df_p::TableP
end

Base.@kwdef mutable struct ProjPhaTYP <: BioinfProj
    pd::String
    phatyp::WrapCmd{RunPhaTYPCmd}
    phatyp_out_df::TableP
    #mergedPostCheckV_PhaTYP::Union{Missing, DataFrame} = missing
    mergedPostCheckV_PhaTYP_p::TableP
end

struct ProjDetectMixedViruses <: BioinfProj
    pd::String
    inDf_Int::TableP
    inDf_NonInt::TableP
    outDf_Int::TableP
    #outFna_Int::FnaP
    outDf_NonInt::TableP
    #outFnaP_NonInt::FnaP
    #outmixed_df::TableP
    #outmixed_fna::FnaP
end

Base.@kwdef mutable struct FinalThresholding <: BioinfProj
    pd::String
    inFna::FnaP
    inFna_trimmed_DTR::Union{Missing, FnaP} = missing
    inTsv::TableP
    predictors::Vector{Symbol}
    predictors2::Union{Missing, Vector{Symbol}} = missing
    th_num_predictors_CheckV_NA::Int64
    th_num_predictors_CheckV_AAIHighConf::Int64
    th_completeness_CheckV_AAIHighConf::Float64
    th_num_predictors_CheckV_AAIMediumConf::Int64
    th_completeness_CheckV_AAIMediumConf::Float64
    th_num_predictors_CheckV_HMM::Int64
    th_completeness_CheckV_HMM::Float64
    th_num_predictors_CheckV_DTR_ITR_AAI::Union{Missing, Int64} = missing
    th_num_predictors_CheckV_DTR_ITR_HMM::Union{Missing, Int64} = missing

    outFnaP::FnaP
    outFna_trimmed_DTR::Union{Missing, FnaP} = missing
    outTsv::TableP
end

Base.@kwdef mutable struct ProjSViP <: BioinfSProj
    projtype::String = "singleworkflow"
    use_slurm::Bool
    continue_project::Bool
    pd::String
    inref::FnaP
    sampleName::String
    sampleSet::String
    dosteps::Dict{String, WorkflowStatus}
    stop_after_initial_predictors::Bool = false
    contig_length::Union{Missing, Proj_sel_contig_length} = missing
    #shape_contigs::Union{DataFrame, Missing} = missing
    genomad::Union{Missing, ProjGenomad} = missing
    DVF::Union{Missing, ProjDVF} = missing
    virSorter2::Union{Missing, ProjVirSorter2} = missing
    vibrant::Union{Missing, ProjVibrant} = missing
    viralVerify::Union{Missing, ProjViralVerify} = missing
    checkV_NonIntegrated::Union{Missing, ProjCheckVNonintegrated} = missing
    checkV_Integrated::Union{Missing, ProjCheckVIntegrated} = missing
    phaTYP_nonintegrated::Union{Missing, ProjPhaTYP} = missing
    #phaTYP_integrated::Union{Missing, ProjPhaTYP} = missing
    detect_mixed_viruses::Union{Missing, ProjDetectMixedViruses} = missing
    final_thresholding_NonIntegrated::Union{Missing, FinalThresholding} = missing
    final_thresholding_Integrated::Union{Missing, FinalThresholding} = missing
    #final_thresholding_Mixed::Union{Missing, FinalThresholding} = missing
end

