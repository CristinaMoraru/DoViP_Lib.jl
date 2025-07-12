module DoViP_Lib

using Revise
using CSV
using DataFrames
using Statistics
using Serialization

using BioS_Gen
using BioS_ProjsWFs
import BioS_ProjsWFs.run_workflow
import BioS_ProjsWFs.run_workflow!
using BioS_SeqFuns
using BioS_ExtTools

include("DoViP_Lib_00_DT_structs.jl")
include("DoViP_Lib_funs2.jl")
include("DoViP_Lib_00_DT_MainProjConstructor.jl")
include("DoViP_Lib_000_Run.jl")
include("DoViP_Lib_funs1.jl")

end # module DoViP_Lib


#=
#region to activate only for testing
args = [
    "projtype=multipleworkflow",
    "spd=/mnt/ASBRU/Projects/RESIST_phaseI_viruses/TEMP-TRANSFER/debug_dovipv09/one_meta",
    "allrefs_params=/mnt/ASBRU/Projects/RESIST_phaseI_viruses/TEMP-TRANSFER/debug_dovipv09/CLM_Julia_EMC_DoViPv0.9_b1-2.tsv",
    "continue=true",
] =#
#=
args = [
    "inref=/mnt/cephfs1/projects/DoViP_benchmarking/test_dataset/inputs_tests/P_2_5_SV21_Ref6_C07_MG_P_2_5_SE_S13_347_length_46575_cov_34.fna",   #/mnt/cephfs1/projects/Prospectomics/Prospectomics_Main_Feb2024/virus_doVIP.a/assembly.d/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_min1000.fasta", 
    "pd_prefix=/mnt/cephfs1/projects/DoViP_benchmarking/test_dataset/outputs",
    "projtype=singleworkflow",
    "sample_set=DoViP_test_Prospectomics_dataset",
    "use_slurm=false",
    "continue=true",
    "stop_after_initial_predictors=false",
    "min_contig_length=1000",
    "merge_circ_proph=true",
    "user=CristinaM",
    
    # genomad related parameters
    "genomad_signal=use",
    "genomad_res=/mnt/cephfs1/projects/Prospectomics/Prospectomics_Main_Feb2024/virus_doVIP.a/virus.a/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_genomad1.10", 
    "genomad_env=conda_genomad_v1.11",
    "genomadDB_p=/software/conda/databases/genomad_v1.11/genomad_db",
    "genomad_min_score=0.7",
    "genomad_sbatch_time=2-0",
    "genomad_cpus_per_task=25",
    "genomad_sbatch_mem=20G",
    
    # DVF related parameters
    "DVF_signal=do",
    "dvf_res=/mnt/cephfs1/projects/Prospectomics/Prospectomics_DoViP_viruspredictions/analysis_round5/b2/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_min1000/01b_ALL_DVF/01b_ALL_01_DVF_out",
    "DVF_env=conda_DVF",
    "DVF_maxContigLen=2099000",
    "DVF_scoreTh=0.7",
    "DVF_pThreshold=0.01",
    "DVF_p=/home/conda/software/DeepVirFinder/dvf.py",
    "DVF_sbatch_time=2-0",
    "DVF_cpus_per_task=15",
    "DVF_sbatch_mem=20G",
    
    # virSorter2 related parameters
    "virSorter2_signal=use",
    "virSorter2_res=/mnt/cephfs1/projects/Prospectomics/Prospectomics_DoViP_viruspredictions/analysis_round2/b2/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_min1000/ALL-01c_virSorter2/ALL-01c_01_virSorter2_out", 
    "virSorter2_env=conda_virsorter2",
    "virSorter2DB_p=/mnt/XIO_3/data1/virsorter-data/virsorter2/db/", 
    "virSorter2_high_confidence_only=false",
    "virSorter2_min_score=0.5",
    "virSorter2_sbatch_time=4-0",
    "virSorter2_cpus_per_task=25",
    "virSorter2_sbatch_mem=20G",
    
    # VIBRANT related parameters
    "vibrant_signal=use",
    "vibrant_res=/mnt/cephfs1/projects/Prospectomics/Prospectomics_Main_Feb2024/virus_doVIP.a/virus.a/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_VIBRANT1.2.1", 
    "vibrant_env=conda_VIBRANT",
    "vibrant_p=/home/conda/software/VIBRANT/VIBRANT_run.py",
    "vibrantDB_p=/home/conda/software/VIBRANT/databases/",
    "vibrant_sbatch_time=2-0",
    "vibrant_cpus_per_task=25",
    "vibrant_sbatch_mem=20G",
    
    # viralVerify related parameters
    "viralVerify_signal=use",
    "viralVerify_res=/mnt/cephfs1/projects/Prospectomics/Prospectomics_DoViP_viruspredictions/analysis_round2/b2/P_2-5_SV21-Ref6-C07_MG_P_2-5_SE_S13_min1000/ALL-01e_viralVerify/ALL-01e_01_viralVerify_out", 
    "viralVerify_env=conda_viralVerify",
    "viralVerifyDB_p=/software/conda/conda_viralVerify/DB/viralverifyDB_nbc_hmms.hmm",
    "viralVerify_p=/software/conda/conda_viralVerify/viralVerify/bin/viralverify",
    "viralVerfify_threshold=7",
    "viralVerify_sbatch_time=2-0",
    "viralVerify_cpus_per_task=25",
    "viralVerify_sbatch_mem=20G",
    
    #checkV related parameters
    "checkv_env=conda_checkV",
    "checkvDB_p=/mnt/XIO_3/data1/CheckV/checkv-db-v1.5/",
    "checkv_sbatch_time=2-0",
    "checkv_cpus_per_task=2",
    "checkv_sbatch_mem=20G", 
    
    #PhaTYP related parameters
    "phaTYP_env=conda_PhaBOX",
    "phaTYP_p=/software/conda/soft/PhaBOX/",
    "phaTYPdb_p=/software/conda/soft/PhaBOX/database/",
    "phaTYPparam_p=/software/conda/soft/PhaBOX/parameters",
    "phaTYP_sbatch_time=2-0",
    "phaTYP_cpus_per_task=2",
    "phaTYP_sbatch_mem=20G",  

    # genomadTax
    "genomadtax_sbatch_time=2-0",
    "genomadtax_cpus_per_task=15",
    "genomadtax_sbatch_mem=20",
    
    # Final thresholding related parameters for NON-INTEGRATED
    "NONInt_th_num_predictors_CheckV_NA=3",
    "NONInt_th_num_predictors_CheckV_AAIHighConf=1",
    "NONInt_th_completeness_CheckV_AAIHighConf=20",
    "NONInt_th_num_predictors_CheckV_AAIMediumConf=2",
    "NONInt_th_completeness_CheckV_AAIMediumConf=20",
    "NONInt_th_num_predictors_CheckV_HMM=2",
    "NONInt_th_completeness_CheckV_HMM=20",
    "NONInt_th_num_predictors_CheckV_DTR_ITR_AAI=1",
    "NONInt_th_num_predictors_CheckV_DTR_ITR_HMM=1",
    
    
    # Final thresholding related parameters for INTEGRATED
    "Int_th_num_predictors_CheckV_NA=2.5",
    "Int_th_num_predictors_CheckV_AAIHighConf=1",
    "Int_th_completeness_CheckV_AAIHighConf=20",
    "Int_th_num_predictors_CheckV_AAIMediumConf=1.5",
    "Int_th_completeness_CheckV_AAIMediumConf=20",
    "Int_th_num_predictors_CheckV_HMM=1.5",
    "Int_th_completeness_CheckV_HMM=20"
]


println("Start DoViP!")
proj = initialize_workflow(args, ProjSViP_fun)

if proj.projtype == "singleworkflow"
    run_workflow(proj)
elseif proj.projtype == "multipleworkflow"
    run_workflowMDoViP(proj)
end

println("DoVip is done!")
#endregion 
#
end =#