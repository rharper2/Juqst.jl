# Generate random channels -

# Copyright Robin Harper 2016-2020


using Distributions
using Random
using LinearAlgebra
using SparseArrays

# Note this is for reproducibility --- for multiple runs or spreadsheets you will need some
# other way to set.
#srand()
Random.seed!

_d=Normal()


#    Based on 1707.06926 (Rudnicki et al) and some others - see workbook for more complete description


"""
    returns the fidelity of the channel (in PauliLiouville basis)
"""
function fidelity(channel)
    diag = tr(channel)-1
    d = sqrt(size(channel)[1])
    p = diag/(d^2-1)
    return ((d-1)*p+1)/d
    #(1+1/3*(channel[2,2]+channel[3,3]+channel[4,4]))/2
end

"""
    returns the unitarity of the channel (PauliLiouville basis)
"""
function unitarity(noise)
    d = sqrt(size(noise)[1])
    return (1/(d^2-1))*tr(noise[2:end,2:end]'*noise[2:end,2:end])
end

"""
    A given fidelity has a minimum unitarity that it can have
    (so perfect fidelity = 1 = unitarity)
    This shows as a percentage how much the unitarity of the channel
    lies between the minimum and the maximum (always 1)
"""
function unitarityPercent(ch)
    fid = fidelity(ch)
    uni = unitarity(ch)
    d = sqrt(size(ch)[1])
    r = 1-fid
    minU = (1-(d*r)/(d-1))^2
    return (uni-minU)/(1-minU)
end

"""
    Simple function that takes a normal distribution and returns a random variable
    that is bounded (non-inclusively) by min and max.
"""
function genTruncated(min,max,average,sigma,normalDistribution)
    x=0.0
    while (x < min) || (x > max)
        x=rand(normalDistribution,1)[1]*sigma+average
    end
    return x
end

"""
    Given the eigenvalues λ return a
    vector of random τ, constrained by the λ
"""
function getτ(λ1,λ2,λ3)
    Maxτ1 = 1-abs(λ1)
    Maxτ2 = 1-abs(λ2)
    Maxτ3 = 1-abs(λ3)
    τ1 = (Maxτ1)*(2*rand()-1)
    τ2 = (Maxτ2)*(2*rand()-1)
    τ3 = (Maxτ3)*(2*rand()-1)
    return [τ1 τ2 τ3]
end

"""
    q(l) and Zη(t,l) See equations 20 and 21 of 1707.06926 (Rudnicki et al)
"""
function q(l)
 return (1+l[1]+l[2]+l[3])*(1 + l[1]-l[2]-l[3])*(1-l[1]+l[2]-l[3])*(1-l[1]-l[2]+l[3])
end

"""
    q(l) and Zη(t,l) See equations 20 and 21 of 1707.06926 (Rudnicki et al)
"""
function Zη(t,l)
    norm(t)^4-2*norm(t)^2-2*sum([l[i]^2*(2*(t[i]^2)-norm(t)^2) for i=1:3]) + q(l)
end


"""
    Main function to generate the 'internal' matrix of the channel i.e. the one
    in the form where the unital part consists of diagonal numbers, and the we have
    a non-unital part. This can then be acted on by unitaries either side (rotations)
    to generate an arbitarty CPTP channel
"""
function genSigma(min,max,average,sigma)
    λ1=genTruncated(min,max,average,sigma,_d)
    λ2=genTruncated(min,max,average,sigma,_d)
    λ3=genTruncated(min,max,average,sigma,_d)
    while 1+λ3 < abs(λ1+λ2) || 1-λ3 < abs(λ2-λ1)
        λ3=genTruncated(min,max,average,sigma,_d)
    end
    τ = getτ(λ1,λ2,λ3)
    while norm(τ)^2 > 1-sum([x^2 for x in [λ1 λ2 λ3]])+2*λ1*λ2*λ3 || Zη(τ,[λ1 λ2 λ3]) < 0
        τ = getτ(λ1,λ2,λ3)
    end
    return [1 0 0 0;τ[1] λ1 0 0;τ[2] 0 λ2 0;τ[3] 0 0 λ3]
end

"""
    genChannel(r1,sm,sv,r2)
    Rotate the sigma channel randomly to generate an arbitrary channel
    r1 represents a rotation. The elements of a vector of three normally distributed elements are multiplied by this and used as an XZX rotation. For high fidelity channels this will be small.
    sm, represents the diagonal channel. This is the mean. High fidelity channels will have this close to 1.
    sv, represents the variance in the diagonal channel.
    r2, like r1.

    For example randomFidelityNoise is defined as   genChannel(0.06,0.998,0.04,0.06)
    whereas the noisier randomMeasureNoise is       genChannel(0.05,0.98,0.15,0.05)
"""
function genChannel(r1,sm,sv,r2)
    createRotation(rand(_d,3)*r1) * genSigma(0.9,1,sm,sv) * createRotation(rand(_d,3)*r2)'
end

"""
    given a vector of three angles, generate a random rotation channel
    by using the vector to control a Z*X*Z rotation
"""
function createRotation(angles)
    rotation = Matrix{Float64}(I,4,4)
    rotation[2:4,2:4] =
        [cos(angles[1]) sin(angles[1]) 0;-sin(angles[1]) cos(angles[1]) 0;0 0 1] *
        [1 0 0;0 cos(angles[2]) sin(angles[2]);0 -sin(angles[2]) cos(angles[2])] *
        [cos(angles[3]) sin(angles[3]) 0;-sin(angles[3]) cos(angles[3]) 0;0 0 1]
    return rotation
end


"""
    Helper function to generate typical, high fidelity channel.
"""
function randomFidelityNoise()
    return genChannel(0.06,0.998,0.04,0.06)
end

"""
    Helper function to generate typical, lower fidelity channel.
    I use it to simulate a state preparation channel.
"""
function randomPrepNoise()
    return genChannel(0.05,0.985,0.15,0.05)
end

"""
    Helper function to generate typical, lower fidelity channel.
    I typically use it to simulate measurement error.
"""
function randomMeasureNoise()
    return genChannel(0.05,0.98,0.15,0.05)
end


"""
    function genChannelMap(f)
    Simple helper function - designed to help generate a map of a specific fidelity (supplied parameter).
    Just generates random channels until you get one of the correct fidelity \$\\pm 0.0001\$.
    In general though you would want to rewrite this to give genChannel
    parameters that will result in random channels that often bracket the desired fidelity
    May take a long time if the parameters are far from the types of channels generated
    by randomFidelityNoise.
"""
function genChannelMap(f)
    while true
        candidate = randomFidelityNoise()
        if abs(f-fidelity(candidate)) < 0.0001
            return candidate
        end
    end
end


pI=[1 0;0 1]
pZ=[1 0;0 -1]


# This can be a bit slow, because of the makeSuper on large qubits.
# I have done this better somewhere.
"""
    function add2QPauliNoise(originalMap,cbit,tbit;pauliToAdd = pZ)

    Adds some noise between two qubits in the noise map.
    You pass in an  noise map (super operator basis)
    the control bit and the target bit and optionally the axis of rotation for the noise.
    will return a noise map with a small amount of random noise added between these qubits.
    The noise itself is calculated as \$e^(i*\\frac{\\text{rand()}}{20}\\text{pauliToAdd})\$
"""
function add2QPauliNoise(originalMap,cbit,tbit;pauliToAdd = pZ)
    nqubitsInMap = round(Int,log2(sqrt(size(originalMap,1))))
    #print("Map is $nqubitsInMap qubits\n")
    ok = false
    ## So we are getting irritating (what I think are) rounding errors here, best way to avoid.
    zz = []
    zeroSelect = SparseArrays.sparse([1; 0]*[1 0])
    oneSelect  = SparseArrays.sparse([0;1]*[0 1])
    spI = SparseArrays.sparse(pI)
    thetaa = rand()/10
    thetab = rand()/10
    if cbit < tbit
            startC=cbit 
            endC=tbit 
             l1 = SparseArrays.sparse(exp(im*thetaa/2*pI))
             l2 = SparseArrays.sparse(exp(im*thetab/2*pauliToAdd))
             for i=cbit+1:tbit-1
                   l1=kron(spI,l1)
                   l2=kron(spI,l2)
             end
            cbitbit = kron(zeroSelect,l1)+kron(oneSelect,l2)
    else 
            startC=tbit 
            endC=cbit 
            l1 = exp(im*thetaa/2*pI)
            l2 = exp(im*thetab/2*pauliToAdd)
            for i=tbit+1:cbit-1 
                 l1=kron(l1,spI)
                 l2=kron(l2,spI)
            end
            cbitbit = kron(l1,zeroSelect)+kron(l2,oneSelect)
    end
    for i = 1:(min(cbit,tbit)-1)
            cbitbit = pI ⊗ cbitbit
    end
    for i = max(cbit,tbit)+1:nqubitsInMap
            cbitbit = cbitbit ⊗ pI
    end
    zz = makeSuper(Array(cbitbit))
    withZ = zz*originalMap
    return withZ
end



