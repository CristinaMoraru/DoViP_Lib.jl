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

include("DoViP_Lib_00_DT_structs_VS.jl")
include("DoViP_Lib_funs2_VS.jl")
include("DoViP_Lib_00_DT_MainProjConstructor_VS.jl")
include("DoViP_Lib_000_Run_VS.jl")
include("DoViP_Lib_funs1_VS.jl")

#end # module DoViP_Lib



#=region to activate only for testing
args = [
    "projtype=multipleworkflow",
    "spd=/mnt/cephfs1/projects/Prospectomics/Prospectomics_DoViP_viruspredictions/a_r6_o/b3",
    "allrefs_params=/mnt/cephfs1/projects/Prospectomics/Prospectomics_DoViP_viruspredictions/params/round6_ODIN/CLM_Prospectomics_DoViPv1.1.1_b3.tsv",
    "continue=false",
] =#
#=
args = [
    "inref=/mnt/cephfs1/projects/DeepH/MSI/MSI_MultiflockONT_2024/polish.d/assembly.d/Polish_MSImulti_ONT_min1000.fasta", #/mnt/cephfs1/projects/DeepH/MSI/Hybrid_assembly/assembly.d/MSI_hybridspades3.15.5_min1000.fasta", 
    "pd_prefix=/mnt/cephfs1/projects/DeepH_virus_pred_DoViP/00_DoViPv1.0/analysis/test",
    "projtype=singleworkflow",
    "sample_set=DeepH_HiC",
    "use_slurm=false",
    "continue=true",
    "stop_after_initial_predictors=false",
    "min_contig_length=1000",
    "merge_circ_proph=true",
    "user=CristinaM",
    
    # genomad related parameters
    "genomad_signal=use",
    "genomad_res=/path/to/res", 
    "genomad_env=conda_genomad_v1.11",
    "genomadDB_p=/software/conda/databases/genomad_v1.11/genomad_db",
    "genomad_min_score=0.7",
    "genomad_sbatch_time=2-0",
    "genomad_cpus_per_task=25",
    "genomad_sbatch_mem=20G",
    
    # DVF related parameters
    "DVF_signal=use",
    "dvf_res=/path/to/res",
    "DVF_env=conda_DVF",
    "DVF_maxContigLen=2099000",
    "DVF_scoreTh=0.7",
    "DVF_pThreshold=0.005",
    "DVF_p=/home/conda/software/DeepVirFinder/dvf.py",
    "DVF_sbatch_time=2-0",
    "DVF_cpus_per_task=15",
    "DVF_sbatch_mem=20G",
    
    # virSorter2 related parameters
    "virSorter2_signal=use",
    "virSorter2_res=/path/to/res", 
    "virSorter2_env=conda_virsorter2",
    "virSorter2DB_p=/mnt/XIO_3/data1/virsorter-data/virsorter2/db/", 
    "virSorter2_high_confidence_only=false",
    "virSorter2_min_score=0.5",
    "virSorter2_sbatch_time=4-0",
    "virSorter2_cpus_per_task=25",
    "virSorter2_sbatch_mem=20G",

    # virSorter related parameters
    "virSorter_signal=use_external",
    "virSorter_res=/mnt/cephfs1/projects/DeepH/Hi_C_202503_MSI/Hi_C_Alti1_interaction/virsorter1/virsorter_out_polishedONT",  #/mnt/cephfs1/projects/DeepH/Hi_C_202503_MSI/Hi_C_Alti1_interaction/virsorter1/virsorter_out_polishedONT 
    
    # VIBRANT related parameters
    "vibrant_signal=do",
    "vibrant_res=/path/to/res", 
    "vibrant_env=conda_VIBRANT",
    "vibrant_p=/home/conda/software/VIBRANT/VIBRANT_run.py",
    "vibrantDB_p=/home/conda/software/VIBRANT/databases/",
    "vibrant_sbatch_time=2-0",
    "vibrant_cpus_per_task=25",
    "vibrant_sbatch_mem=20G",
    
    # viralVerify related parameters
    "viralVerify_signal=do",
    "viralVerify_res=/path/to/res", 
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
    "genomadtax_cpus_per_task=40",
    "genomadtax_sbatch_mem=20",
    
    # Final thresholding related parameters for NON-INTEGRATED
    "NONInt_th_num_predictors_CheckV_NA=1",
    "NONInt_th_num_predictors_CheckV_AAIHighConf=1",
    "NONInt_th_completeness_CheckV_AAIHighConf=20",
    "NONInt_th_num_predictors_CheckV_AAIMediumConf=1",
    "NONInt_th_completeness_CheckV_AAIMediumConf=20",
    "NONInt_th_num_predictors_CheckV_HMM=1",
    "NONInt_th_completeness_CheckV_HMM=20",
    "NONInt_th_num_predictors_CheckV_DTR_ITR_AAI=1",
    "NONInt_th_num_predictors_CheckV_DTR_ITR_HMM=1",
    
    
    # Final thresholding related parameters for INTEGRATED
    "Int_th_num_predictors_CheckV_NA=1.0",
    "Int_th_num_predictors_CheckV_AAIHighConf=1",
    "Int_th_completeness_CheckV_AAIHighConf=20",
    "Int_th_num_predictors_CheckV_AAIMediumConf=1.0",
    "Int_th_completeness_CheckV_AAIMediumConf=20",
    "Int_th_num_predictors_CheckV_HMM=1.0",
    "Int_th_completeness_CheckV_HMM=20"
] #


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