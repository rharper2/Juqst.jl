module Marginal

# Note all the graphical stuff and IBM specific stuff has been moved into a different file
# marginalDrawing.jl which should be included if you want to use it.



using LsqFit, Hadamard, DelimitedFiles, LinearAlgebra
using PyPlot, Statistics
import Base.show

export readInCSVFile,transformToFidelity,fitTheFidelities,convertAndProject,gibbsRandomField
export getGrainedP,mutualInformation,relativeEntropy,conditionalMutualInfo,covarianceMatrix,JSD
export correlationMatrix

⊗ = kron

""" 
readInCSVFile(filename)

## Arguments 
	Filename: full name of a csv delimited Filename

## File Format
	Assumes that each sequence appears on a new line and each line has the same number of entries.
	Assumes that it is a raw count for each of the possible measurement outcomes.
	Note it is important to know the order of the outcomes (this is significant).
	Although this can be overriden, it is a assumed that they are stored as:

```
			00000000
			00000001
			00000010
			00000011
			
			...
			11111110
			11111111

		i.e. lsb to the right.
```
Returns a full rows x col matrix depending on the number of entries of the file.
for a 16 qubit machine, where the data was taken over 6 sequences we will return a 6x655536 array of Ints.

"""
function readInCSVFile(filename)
	return readdlm(filename,',',Int64)
end

"""
	transformToFidelity(data)

##	Arguments:
	Data, as returned by readInCSVFile - see help for that function.
		
## Returns:
	The Hadamard transformed fidelities for each sequence as a list of lists.

"""
function transformToFidelity(data)
	splitMatrix = [data[i,:]/sum(data[i,:]) for i in 1:size(data,1)]
	return [ifwht_natural(x) for x in splitMatrix]
end


modelF(x, p) = p[1]*(p[2].^x)

"""
	fitTheFidelties(lengths, data; no =0)

## Arguments
-	`lengths: Array{Int64,1}`, e.g. collect(1:2:24) - the length at which each of the above sets of data was gathered.

-   `data: Array{Array{Float64,1},1}`, this is a list of observed probabilities at different sequences.

-   `no:` printed with the warning message if we fail to fit. Allowing identification of problem in "batch" fits.

Take the data, transform it  and see if we can fit it to an exponential decay.
This proves problematic where the data is a bit sparse, as the decays can be all over the place.
In general it is *extremely* helpful to have a 0 gate sequence as that helps to anchor the fits.
	
Here we check for convergence, and warn if it doesn't converge.

This enforces a cut-off to the data used for the fit. Specifically it gets rid of the 'tail' of data.

## Returns
	three things:
	-  the fitting parameters, 
	-  the length of the sequence that was used for each parameter.
	-  the indices of those that still failed to converge.

## Typical usage
```
	params,_ = fithTheFidelities(transformToFidelity(data),collect(1:2:22))
```
"""
function fitTheFidelities(lengths,matrix;no=0)
    revh = [ifwht_natural(x) for x in matrix]
    params = Array{Float64,1}[]
    dataCounted =[]
    failedToC = []
    _negCount =0
    for idx = 2:length(revh[1])
        extract = map(x->x[idx],revh)
        # The idea here is as per https://arxiv.org/abs/1901.00535 i.e. don't let the tail wag.
        # It doesn't make a large difference, but it does make a difference.
        p1 = (extract[1]*(1+1/16))/4
        lastData = findfirst(x->p1>x,extract)
        if lastData == nothing
            lastData = length(extract)
        end
        # two is probably sufficient, but at 3 for now.
        if lastData < 3
            lastData = 3
        end
        fit =[]
        try
            fit = curve_fit(modelF,lengths[1:lastData],extract[1:lastData],[0.8,0.8],upper=[1.0,1.0],lower=[0.01,0.01])
        catch 
            print("Fall back on $no, $idx\n")
            fit = curve_fit(modelF,lengths[1:lastData],extract[1:lastData],[0.8,0.4],upper=[1.0,1.0],lower=[0.01,0.01])
            
        end
        if !fit.converged
            push!(failedToC,idx)
        end
        push!(params,fit.param)
        push!(dataCounted,lastData)
    
    end
    if (_negCount > 0) print("$no: There were $(_negCount) with a negative first number.\n") end
    if (length(failedToC) > 0) print("$no: $(length(failedToC)) failed to converge.\n") end
    return params,dataCounted,failedToC
end        
   



"""
	Projects the P's onto the simplex where they are all >= 0
"""
function projectSimplex(ps)
    sortedSeq = reverse(sort(ps))
    # This is a bit crytic.
    #   a) Sort (greatest first)
    #   b) create a culmlative sum (cumsum)
    #   c) we want the last one where (Cum[k]-1)/k < sorted[k]
    #   d) we then reduce all probs by the last value we had in c, setting them to zero if they are negative.
    # (The reason why we divide by the length in c is because we subtract it more times the further into the index we are)
    τ=filter(x->x[1]<x[2],collect(zip((cumsum(sortedSeq).-1)./(1:length(ps)),sortedSeq)))[end][1]
    projected = map(x->max(x-τ,0),ps)
    return projected
end




"""
	Convert and project the fitted fidelities, to rerieve the 'p's representing the joint probabilities.

	Input: The parameters that have been fit.
	Output: The projected ps.
"""
function convertAndProject(params)
	ps = fwht_natural(map(x->x[2],params))
	pps = projectSimplex(ps)
	if !(isapprox(sum(pps),1))
		@warn "There may be a problem, the probabilities do not add to 1."
	end
	return pps
end


function getIndices(qs;dimension=16)
    sqs = sort(qs)
    indices = [2^(sqs[1]-1),2]
    for i in 2:length(sqs)
            push!(indices,2^(sqs[i]-sqs[i-1]-1))
            push!(indices,2)
    end
    push!(indices,2^(dimension-sqs[end]))
    return indices
end


""" 
    Returns the marginalised probability, given the (projected) probabilities in pps.
"""
function marginalise(q,pps)
    # Get indices sorts the entries
    dimension = 0
    try
    	dimension = Integer(log(2,length(pps)))
    catch
    	@warn "Error, the size of pps needs to be an integer power of 2"
    	return 
    end
    if 0 in q 
        print("Please index off 1, a 0 in the bits to marginalise over will cause grief")
        return nothing
    end
    x = reshape(pps,tuple((getIndices(q,dimension=dimension))...))
    # Below the permutation that will get us back.
    sorted = sort(q)
    # This is the reverse of sortperm is from sorted -> original
    permute = [findfirst(isequal(x),sorted) for x in q]
  
    # print(permute)
    # Sum over all the odd ones.
    marginalised = sum(x,dims=1:2:length(size(x)))
    # Work out how many 2 variables we have
    indices = length(q)
    # Shove them in a tuple (note splat so it works well with reshape)
    fullIndices = tuple([2 for i in 1:indices]...)
    # Might as well get it in a 2 dim array to return
    # Rows has to be even of course.
    rows = 2^(floor(Int,indices/2))
    if rows == 0
        rows = 1
    end
    cols = round.(Int,2^indices/rows)
    return reshape(permutedims(reshape(marginalised,fullIndices),permute),rows,cols)
end 


"""
gibbsRandomField(pps,constraints)

## Arguments
-	`pps: Array{Float64,1}` For example the output of convertAndProject
-	`constraints: Array{Int64,1}`: The division of the gibbs field. See below for example.

	Takes a joint probability and the gibbs variable constraints and returns a vector of reduced probability distributions.
	For example the constraints might be [[1,2,3,4],[3,4,5,6],[5,6,7,8]] over a field of 8 qubits
	This would return an array of 3 ϕ's each 16 long (2^4), from which the joint probability can be extracted on the assumption
	That, say, qubits 1,2 are independent of qubits 5-8.
    Note: the assumption inbuilt here is that the constraints end and start with qubits in common e.g. [1,2,->3,4],[3,4<-,5,6]
   

Note the other version of is one gibbsRandomField, where the constraints are constraints::Array{Array{Array{Any,1},1},1}
is probably (now) preferred - there you can specify the constraints more generally.
In particular as tuples of qubits and the qubits they are conditioned on.

"""
function gibbsRandomField(pps,constraints)
    overlaps = []
    for (idx,i) in enumerate(constraints[1:end-1])
        x = findfirst(x->x==i[end],constraints[idx+1])
        if x == nothing
            print("Error in finding overlaps the constraints need to overlap at least slightly\n")
            return nothing
        end
        push!(overlaps,x)
    end
    pxx = [marginalise(x,pps) for x in constraints]
    for i in 1:length(overlaps)
        pxx[i] = reshape(pxx[i],:,2^overlaps[i])
    end # reshape to correspand to the overlap.
    px =  [vec(marginalise(x,pps))' for x in [i[end-(overlaps[idx]-1):end] for (idx,i) in enumerate(constraints[1:end-1])]]
    # Map any zero entries in px to something non-zero or division will fail.
    # Not an issue as if px is zero, the corresponding pxx has to be zero. (and in this case by definition the probability is zero)
    pxNz = [map(x-> x == 0 ? 1e-8 : x,pxidx) for pxidx in px]
    ϕ = [vec(pxx[i]./pxNz[i]) for i in 1:(length(constraints)-1)]
    push!(ϕ, vec(marginalise(constraints[end],pps)))
    return ϕ
end

"""
    Helper function for the second form of gibbsRandomField. Takes a list of probability tuples and flattens it out.
    So [(1,2),(3,4)],[(5,),(6,7)]] becomes [[1,2,3,4],[5,6,7]] (read the first as p(1,2|3,4)p(5|6,7))
"""
function jointIfy(x)
    return [collect(Iterators.Flatten(y)) for y in x]
end

"""
gibbsRandomField(pps,gen_constraints)

## Arguments
-   `pps: Array{Float64,1}` For example the output of convertAndProject
-   `constraints: Array{Array{Array{Any,1},1},1}`: The division of the gibbs field. See below for example.

    Takes a joint probability and the variable constraints that are conditioned and returns a vector of reduced probability distributions.
    It is up to the user to make sure that the factorization is complete and makes sense.
    For example the constraints might be [[(1,2),(3,4)],[(3,),(4,5,6)],[(4,5,6,7,8),()]] over a field of 8 qubits
    This would return an array of 3 ϕ's the first 16 long (2^4), the second 2^4 long and the third 2^5 long from which the joint probability can be extracted on the assumption
    That, say, qubits 1,2 are independent of qubits 5-8.
"""
function gibbsRandomField(pps,constraints::Array{Array{Array{Any,1},1},1})
 # This is more flexible than the previous one, which took overlapping contraints, to try and work out
 # what the qubits were conditioned on, e.g. if give [[1,2,3,4],[3,4,5,6]], it would assume that you wanted:
 # p(1,2|3,4)p(3,4,5,6)
 
    
# Here we can expressly state the dependancies.
# It is expecting contraints as follows a list of [ [[bits],[conditioned on]],[[bits],[conditioned on]]]
# eg [[[1,2],[3,5]],[[3,4],[5]],[[5],[]]]
# means that we want p(1,2|3,5)p(3,4|5)p(5)
# Which will be returned as a 2^4, 2^2, 2^1 set of arrays.
    
    overlaps = [length(i[2]) for i in constraints]
    jointMarginals = jointIfy(constraints)
    
    pxx = [Marginal.marginalise(x,pps) for x in jointMarginals]
    for i in 1:length(overlaps)
        pxx[i] = reshape(pxx[i],:,2^overlaps[i])
    end # reshape to correspand to the overlap.
    
    px =  [x == [] ? 1 : vec(Marginal.marginalise([x...],pps))' for x in [i[2] for i in constraints]]
    # Map any zero entries in px to something non-zero or division will fail.
    # Not an issue as if px is zero, the corresponding pxx has to be zero. (and in this case by definition the probability is zero)
    pxNz = [map(x-> x == 0 ? 1e-8 : x,pxidx) for pxidx in px]
    ϕ = [vec(pxx[i]./pxNz[i]) for i in 1:(length(constraints))]

    # If we did just allow them to be zero, then it causes problems reconstructing the probability distribution.
    # Here we find where the reshaped vector e.g. a conditioned on b,c,d 
    # is equal to zero (it should be 1, i.e.the culmlative percentage chance of a being 0 or 1 must be 1.)
    # We set the chance of it being 0 to be the average of the non 'all zero' chances. Same with 1  etc.
    # Generalised where we have a joint distribution on the left hand size (reshapeFactor)
    for (cidx,toC) in enumerate(ϕ)
        reshapeFactor = 2^length(constraints[cidx][1])
        for (idx,p) in enumerate(sum(reshape(toC,reshapeFactor,:),dims=1))
            if p==0
                newVals = zeros(reshapeFactor)
                #display(reshape(toC,reshapeFactor,:))
                for i in 1:reshapeFactor
                    newVals[i] =  (mean(reshape(toC,reshapeFactor,:)[i,1:idx-1]))
                end
                for i in 1:reshapeFactor
                    reshape(ϕ[cidx],reshapeFactor,:)[i,idx]=newVals[i]                
                end
            end
        end
    end
    return ϕ
end


# The function below is probably worth a bit of comment. 
# uses the information in 'graining' to tells us which index the pattern we want resides in
# For instance say we want all bits zero apart from qubit 3 (which = 1)
# Then for ϕ[1], which represents bits [2,1,3,6] - we need the entry corresponding to decimal equiv of reverse(0,0,1,0)
# plus one, because 0 0 0 0 = first entry i.e. index 1.
# We will have length(ϕ) of these entries, this gives us the perVecIndex.
# Then we just read it out for ϕ, with a map, using foldl to multiply them for us.


"""
	For a given ϕ, the list of qubits we used to construct ϕ (graining) and the 
    Bit pattern, the probability of which we want to extract
"""
function getGrainedP(ϕ,tomatch,graining)
    vectorIndex = [[x&tomatch>0 ? 1 : 0 for x in (map(x->2^(x-1),i))] for i in graining]
    perVecIndex = map(x->parse(Int64,join(reverse(x)),base=2)+1,vectorIndex)
    return foldl(*,map(x->ϕ[x][perVecIndex[x]],1:length(ϕ)))
end

"""
    function mutualInformation(p1,p2,p)

    Returns the mutualInformation between probabilities at p1, p2 given probability distribution p.
    Makes use of the marginalise function to marginalise over probabilities

    \$MI(p1,p2) =\\sum\\limit_{x\\in p1,y\\in p2}\\left[p(x,y)\\log{\\frac{p(x,y)}{p(x),p(y)})\\right]\$


    Some gotchas:
        I return -1 if p1==p2 , this allows an identification of the qubit where I loop.
        if p(x,y) is 0, p(x) or p(y) is 0 then this is defined as a 'zero' part of the sum - even though the log is undefined.
        
"""
function mutualInformation(p1,p2,p)
    # Allow this for easy looping and identificaton of controlling qubit.
    #print("p1 in $p1\n")
    #print("p2 in $p2\n")
    for x in [p1...]
        if x in [p2...] 
            return -1
        end
    end
    p12 = Marginal.marginalise([p1...,p2...],p)
    #println(p12)
    p1 = sum(p12,dims=2)
    #print("p1 $p1\n")
    p2 = sum(p12,dims=1)
    #print("p2 $p2\n")
    logPart = (log2.(vec(p12)) - log2.(vec(p1⊗p2)))
    # Note that the log of 0 will give -Inf or NAN (if p1 or p2 is zero as well). 
    # The log can only be zero (or NAN) if the relevant
    # p12 entry is 0, which by definition means the MI contribution is 0, but of course -Inf*0 = Nan
    # So we zap the -Inf's (and Nan's) to zero which gives us the result we want.
    toReturn = vec(p12)'*map(x->x==-Inf||isnan(x) ? 0 : x,log2.(vec(p12)) - log2.(vec(p1⊗p2)))
    return toReturn
end



"""
	relativeEntropy(P, Q)

	Calculates the relative entropy between two joint probability distributions.
	D(P||Q) = Sum p_j*log(p_j/q_j)

	Undefined if any of the q_js are zero, unless the corresponding p_j is zero.
"""
function relativeEntropy(P,Q)
	# First of all assert that where p̃ is zero pps is also zero
    if !(length(filter(x->x!=0,P.*Q)) == length(filter(x->x!=0,P)))
            @warn "D(P||Q) - the distribution (Q) has zero values not matched by by zero values in (P) - this is undefined."
            return NaN
    end

    # We only need the values of P and Q where P is non zero
    nonZero = vec(map(x->x!=0,P))
	xp = P[nonZero]
	xq = Q[nonZero]
	pDq = xp./xq
	return xp'*log2.(pDq)
end

"""
    Jensen-Shannon Divergence - symmetric
"""
function JSD(dist1,dist2)
    m = 0.5.*(dist1.+dist2)
    return 0.5*(relativeEntropy(dist1,m)+relativeEntropy(dist2,m))
end

"""
    Given a distribution and a set of constraints.
    Reconstruct the distribution using the constraints.
    Then return the JSD of the reconstructed distribution and the original distribution.
"""
function reconstructedJS(distribution,constraints)
    return JSD(reconstruct(distribution,constraints),distribution)
end

"""
    Given a distribution and a set of constraints.
    Return the reconstruction of a distribution parameterized by the constraints.
"""
function reconstruct(dist,constraints)
    gibbs2ϕ = gibbsRandomField(dist,constraints);
    reconstructedPps2 = [getGrainedP(gibbs2ϕ,tomatch,constraints) for tomatch=0:(2^14-1)];
    @assert(isapprox(sum(reconstructedPps2),1))
    return reconstructedPps2
end

"""
    Note to self this will often be undefined for the types of distributions I am currently looking at
"""
function relEntropy(dist,constraints)
    return relativeEntropy(reconstruct(dist,constraints),dist)
end




""" 
    Helper function that given an xyz correctly indexes and pulls out the
    xyz from the joint probability distribution, jpXYZ
    and the xz and yz from jpXZ and jpYZ 
    then does the calculation, inside the sum for that x,y,z
"""
function getSummand(x,y,z,jpXYZ,jpXZ,jpYZ,jpZ)
    multipland = size(jpYZ)[1]
    pxyz = jpXYZ[multipland*(y-1)+x,z]
    if pxyz == 0 
        return 0
    end
    pz = jpZ[z]
    pxz = jpXZ[x,z]
    pyz = jpYZ[y,z]
    return pxyz*log2((pz*pxyz)/(pxz*pyz))
end




""" 
    Gives the conditional mutual information for the qubits 
    supplied in X,Y and Z where we want I(XY|Z)

	I(XY|Z) = Sum( P(x,y,z)log2( P(Z)P(X,Y,Z)/P(X,Z)P(Y,Z) )


    All the heavy lifting is in getSummand.
"""
function conditionalMutualInfo(X,Y,Z,p)
    XYZ = vcat(X,Y,Z)
    XZ = vcat(X,Z)
    YZ = vcat(Y,Z)
    xSize = 2^length(X)
    ySize = 2^length(Y)
    zSize = 2^length(Z)
    jpXYZ = reshape(marginalise(XYZ,p),:,zSize)
    jpXZ = reshape(marginalise(XZ,p),:,zSize)
    jpYZ = reshape(marginalise(YZ,p),:,zSize)
    jpZ = reshape(marginalise(Z,p),1,zSize)
    return sum([getSummand(x,y,z,jpXYZ,jpXZ,jpYZ,jpZ) for x=1:xSize,y=1:ySize,z=1:zSize])
end



"""
	covarianceMatrix(p)

## Covariance
	We can compute the covariance matrix between the 0/1 random variables representing no error / error. 
	If x is a column vector of bits representing an error pattern, then we can compute the matrix
	Expect_p[(x-μ) (x-μ)^T]
	where μ = Expect_p[x]. 

    If we set reverse to true, then we reverse the digits.
    The way IBM stores its digits we wouldn't normally want to reverse them, but in general
    this might not be true
"""
function covarianceMatrix(p;reverseDigits = false)
# Take a probability distribution on n qubits and compute the covariance matrix as above.
    n = convert(Int,log2(length(p)));
    v = digits.(0:(2^n-1),base=2, pad=n) # the 0/1 random variable of errors
    if reverseDigits
        v = reverse.(v)
    end
    v = convert(Array{Array{Int,1},1}, v);
    mu= sum(p.*v); # compute the mean
    v = hcat(v...).-mu; # center the random variable v
    covmat = p[1]*v[:,1]*v[:,1]';
    for i in 2:2^n
        covmat+=p[i]*v[:,i]*v[:,i]';
    end
    return covmat
end




function marginaliseFromRawData(rawData,constraints,lengths)
    overlaps = []
    for (idx,i) in enumerate(constraints[1:end-1])
        x = findfirst(x->x==i[end],constraints[idx+1])
        if x == nothing
            print("Error in finding overlaps the constraints need to overlap at least slightly\n")
            return nothing
        end
        push!(overlaps,x)
    end
    pxx=[]
    # Marginalise over measurements to get probabilities
    for x in constraints
        marginalMatrix = [Marginal.marginalise(x,ms) for ms in rawData]
        paramsM,_ = fitTheFidelities(lengths,marginalMatrix)
        pm =  fwht_natural(vcat([1],map(x->x[2],paramsM)))
        push!(pxx,vec(Marginal.projectSimplex(pm))');
    end
    for i in 1:length(overlaps)
        pxx[i] = reshape(pxx[i],:,2^overlaps[i])
    end # reshape to correspand to the overlap.
    px=[]
    # Here we can use our previous marginalisation to consistently marginalise the constraints
    for (idx,x) in enumerate(constraints[1:end-1])
        lx=convert(Int,log2(length(pxx[idx])))
        toM = [y for y in lx-overlaps[idx]+1:lx]
        push!(px,vec(Marginal.marginalise(toM,pxx[idx]))');
    end
    # Map any zero entries in px to something non-zero or division will fail.
    # Not an issue as if px is zero, the corresponding pxx has to be zero. (and in this case by definition the probability is zero)
    pxNz = [map(x-> x == 0 ? 1e-8 : x,pxidx) for pxidx in px]
    ϕ = [vec(pxx[i]./pxNz[i]) for i in 1:(length(constraints)-1)]
    # Marginalise over measurements to get probabilities
    push!(ϕ, vec(pxx[end]));
    return ϕ
end


# Take a probability distribution on n qubits and compute the covariance matrix as above.
function correlationMatrix(p)
    M = covarianceMatrix(p,reverseDigits=false)
    d = sqrt(inv(LinearAlgebra.Diagonal(M)));
    return d*M*d
end





end

