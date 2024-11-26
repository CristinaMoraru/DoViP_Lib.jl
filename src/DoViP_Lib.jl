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

#=region to activate only for testing
args = [
    "projtype=multipleworkflow",
    "spd=/mnt/ASBRU/Projects/RESIST_phaseI_viruses/TEMP-TRANSFER/debug_dovipv09/one_meta",
    "allrefs_params=/mnt/ASBRU/Projects/RESIST_phaseI_viruses/TEMP-TRANSFER/debug_dovipv09/CLM_Julia_EMC_DoViPv0.9_b1-2.tsv",
    "continue=true",
] #

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
#endregion =#

end # module DoViP_Lib

