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

