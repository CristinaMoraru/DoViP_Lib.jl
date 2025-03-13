module DoViP_Lib

using Revise
using CSV
using DataFrames
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
] #

args = [
    "inref=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/in/3300029793_withISs.fna", 
    "pd_prefix=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/out/",
    "projtype=singleworkflow",
    "sample_set=EMC",
    "use_slurm=false",
    "continue=false",
    "stop_after_initial_predictors=false",
    "min_contig_length=1000",
    
    # genomad related parameters
    "genomad_signal=use_external",
    "genomad_res=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/metagenomes_errors_apply_th/3300029793_withISs/ALL-01a_genomad/ALL-01a_01_genomad_out", 
    "genomad_env=conda_genomad",
    "genomadDB_p=path/to/genomad_db/genomad_db/",
    "genomad_min_score=0.7",
    "genomad_sbatch_time=2-0",
    "genomad_cpus_per_task=25",
    "genomad_sbatch_mem=20G",
    
    # DVF related parameters
    "DVF_signal=do",
    "DVF_env=conda_DVF",
    "DVF_maxContigLen=2099000",
    "DVF_scoreTh=0.7",
    "DVF_pThreshold=0.01",
    "DVF_p=/home/conda/software/DeepVirFinder/dvf.py",
    "DVF_sbatch_time=2-0",
    "DVF_cpus_per_task=15",
    "DVF_sbatch_mem=20G",
    
    # virSorter2 related parameters
    "virSorter2_signal=use_external",
    "virSorter2_res=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/metagenomes_errors_apply_th/3300029793_withISs/ALL-01c_virSorter2/ALL-01c_01_virSorter2_out", 
    "virSorter2_env=conda_virsorter2",
    "virSorter2DB_p=path/to/virsorter-data/virsorter2/db/", 
    "virSorter2_high_confidence_only=false",
    "virSorter2_min_score=0.5",
    "virSorter2_sbatch_time=4-0",
    "virSorter2_cpus_per_task=25",
    "virSorter2_sbatch_mem=20G",
    
    # VIBRANT related parameters
    "vibrant_signal=use_external",
    "vibrant_res=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/metagenomes_errors_apply_th/3300029793_withISs/ALL-01d_vibrant/ALL-01d_01_vibrant_out", 
    "vibrant_env=conda_VIBRANT",
    "vibrant_p=path/to/VIBRANT/VIBRANT_run.py",
    "vibrantDB_p=/path/to/VIBRANT/databases/",
    "vibrant_sbatch_time=2-0",
    "vibrant_cpus_per_task=25",
    "vibrant_sbatch_mem=20G",
    
    # viralVerify related parameters
    "viralVerify_signal=use_external",
    "viralVerify_res=/mnt/ASBRU/Projects/ProkIS_viruses/JGI_EMC/TEMP/metagenomes_errors_apply_th/3300029793_withISs/ALL-01e_viralVerify/ALL-01e_01_viralVerify_out", 
    "viralVerify_env=conda_viralVerify",
    "viralVerifyDB_p=path/to/conda_viralVerify/DB/viralverifyDB_nbc_hmms.hmm",
    "viralVerify_p=path/to/conda/conda_viralVerify/viralVerify/bin/viralverify",
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
    "Int_th_num_predictors_CheckV_NA=3",
    "Int_th_num_predictors_CheckV_AAIHighConf=1",
    "Int_th_completeness_CheckV_AAIHighConf=20",
    "Int_th_num_predictors_CheckV_AAIMediumConf=2",
    "Int_th_completeness_CheckV_AAIMediumConf=20",
    "Int_th_num_predictors_CheckV_HMM=2",
    "Int_th_completeness_CheckV_HMM=20",

]

#
if "--help" in args
    println("DoViP - a workflow for virus prediction in metagenomes.")
    show_file_content("../README.md")
end

println("Start DoViP!")
proj = initialize_workflow(args, ProjSViP_fun)

if proj.projtype == "singleworkflow"
    run_workflow(proj)
elseif proj.projtype == "multipleworkflow"
    run_workflowMDoViP(proj)
end

println("DoVip is done!")
#endregion 
=#
