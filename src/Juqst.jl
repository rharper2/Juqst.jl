module Juqst

export 	Tableau,
		cnot,
		hadamard,
		X,
		Z,
		phase,
		measure,
		kets,
		reinitialise,
		initialise,
		initialize

using LinearAlgebra


mutable struct Tableau
	state::Array{Int32,2}
	qubits::Int32
	showRaw::Bool
	commands::Array{String,1}
	executeCommands::Array{Expr,1}
	trackCommands::Bool
end

include("Initial.jl")


export 	cliffordToTableau,
		tableauToClifford,
		decomposeState,
		drawCircuit,
		qiskitCircuit,
		getNumberOfCliffords,
		getNumberOfSymplecticCliffords,
		getNumberOfBitStringsCliffords,
		generateRawCliffords,
		makeFromCommand

include("Symplectic.jl")
include("open-systems.jl")

export mat,
       liou,
       choi_liou_involution,
       swap_involution,
       choi2liou,
       choiX2liou,
       choi2kraus,
       kraus2choi,
       kraus2liou,
       liou2choi,
       liou2choiX,
       liou2kraus,
       liou2pauliliou,
       pauliliou2liou,
       depol,
       istp,
       iscp,
       ischannel,
       isunital,
       nearestu,
       nicePrint,
       writemime,
	   makeSuper,
	   choi2chi

include("rchannels.jl")

export 	fidelity,
		unitarity,
		unitarityPercent,
		genChannel,
		randomFidelityNoise,
		randomPrepNoise,
		randomMeasureNoise,
		genChannelMap

include("marginal.jl")


export readInCSVFile,transformToFidelity,fitTheFidelities,convertAndProject,gibbsRandomField,marginalise
export getGrainedP,mutualInformation,relativeEntropy,conditionalMutualInfo,covarianceMatrix,JSD
export correlationMatrix,marginaliseFromRawData


end # module
