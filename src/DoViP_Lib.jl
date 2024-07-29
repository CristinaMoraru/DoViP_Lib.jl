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

#=ProjSViP - a single virus prediction workflow
args = [
    "projtype=singleworkflow",
    "pd_prefix=/data3/CLM_projs/TEST_Workflows/batch6", 
    "inref=/mnt/ceph/data/CLM_Nitrospira/ExStreamIMP/viruspred/in/ES22_IMP_S_C12_min1000.fasta",
    #"inref_topology=/home/cmoraru/MY_SPACE/WorkflowsJl/TESTS/test_in_files/DoViP/outMfasta_contig_topology.tsv",

    "continue=true",

    # should I add a general length cutoff for predicted viruses? we have the 1000 bases threshold from Lenni's pipeline
    # check individual tools for their length cutoffs

    "genomad_signal=use",
    "genomad_env=conda_genomad",
    "genomadDB_p=/mnt/XIO_3/data1/genomad_db/genomad_db/",
    "genomad_min_score=0.7",  #default in genomad is 0.7, value between 0 and 1. ################

    "DVF_signal=use",
    "DVF_env=conda_DVF",
    "DVF_minContigLen=2000",
    "DVF_maxContigLen=2099000", #2099000
    "DVF_scoreTh=0.7",
    "DVF_pThreshold=0.01",
    "DVF_p=/home/conda/software/DeepVirFinder/dvf.py",

    "virSorter2_signal=use",
    "virSorter2_env=conda_virSorter2",
    "virSorter2DB_p=/mnt/XIO_3/data1/virSorter2-data/virSorter2/db/",
    "virSorter2_high_confidence_only=false",
    "virSorter2_min_score=0.5",  #default in virSorter2 is 0.5, value between 0 and 1. ################
    "virSorter2_min_length=2000",  #default in virSorter2 is 0 ################

    "vibrant_signal=use",
    "vibrant_env=conda_VIBRANT",
    "vibrant_p=/home/conda/software/VIBRANT/VIBRANT_run.py",
    "vibrantDB_p=/home/conda/software/VIBRANT/databases/",
    "vibrant_min_length=2000",  #default and minimum in vibrant is 1000 ################

    #"cenote-taker3_signal=do",
    #"ct3_env=conda_Cenote-Taker3_v3.2.1",

    "checkv_env=conda_checkV",
    "checkvDB_p=/mnt/XIO_3/data1/CheckV/checkv-db-v1.5/",

    "phaTYP_env=conda_PhaBOX",
    "phaTYP_p=/software/conda/soft/PhaBOX/",
    "phaTYPdb_p=/software/conda/soft/PhaBOX/database/",
    "phaTYPparam_p=/software/conda/soft/PhaBOX/parameters",
    "phaTYP_minlength=1000",  #default in phaTYP is ? ################
                                                    #### In reality, below I don"t need the contamination threshold, because I'm cleaning all contaminated viral contigs after checkV and thus only those with contamination 0 remain
    "th_num_predictors_NonIntegrated=2", # default is 3, Integer
    "th_num_predictors_Integrated=2",  # default is 2
    #"th_num_predictors_Mixed=1",  # default is 1
    "th_checkV_completeness_NonInt=30",  # default is set to 30, Float
    "th_checkV_completeness_Int=30",
    #"th_checkV_completeness_Mixed=50",
    "th_checkV_contamination_NonInt=50",  # default is set to 30, Float
    "th_checkV_contamination_Int=10",
    #"th_checkV_contamination_Mixed=50",

    "num_threads=20"
] =#


#=
if "--help" in args
    println("DoViP - a workflow for virus prediction in metagenomes.")
end

println("Start DoViP!")
proj=initialize_workflow(args, ProjSViP_fun)
Threads.nthreads() = 20

if proj.projtype == "singleworkflow"
    run_workflow(proj)
elseif proj.projtype == "multipleworkflow"
    run_workflowMDoViP(proj)
end

println("DoVip is done!") =#

end # module DoViP_Lib


#= To manipulate the dosteps object, for reruning and testing
proj.dosteps["checkV_NonIntegrated"].signal = "do"
proj.dosteps["phaTYP_nonintegrated"].progress = "not_done"
serialize("$(proj.pd)/sproj.binary", proj)

proj.dosteps["final_thresholding_Integrated"].progress = "not_done"
serialize("$(proj.pd)/sproj.binary", proj)
=#