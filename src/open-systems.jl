# Note this is an altered verion of the file from

#  http://marcusps.github.com
#  Specifically from  https://github.com/BBN-Q/QuantumInfo.jl/blob/master/src/open-systems.jl
#  Original authors: Blake Johnson and Marcus da Silva

# Altered by Robin Harper.

# Original Copyright Notice for QuantumInfo.jl stuff

#Copyright (c) 2014: Raytheon BBN Technologies.

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# With respect to changes made and stuff added:

# Copyright Robin Harper 2015-2020


## I needed to make a few changes and just did it locally. See below.
## I have also decoupled it from using the "Cliffords" package as I don't actually need that one.

#import Base.writemime

import Base.show

# Changed slightly from the version on BitBucket (i.e. the original author's version)
# I added a translate in the pauliliou conversion functions, so that the basis was IXYZ (rather than IXZY)
# I removed a division by d when converting from liou 2 choi (and vice versa) (slightly different conventions)
# The choi2liou uses a different involution - keeping as column stacked.
# There is also a dimensional factor depending on which convention we use.


mat( v::Vector, r=round(Int,sqrt(length(v))), c=round(Int,sqrt(length(v))) ) = reshape( v, r, c )

liou( left::T, right::T ) where T<:AbstractMatrix = kron( transpose(right), left )

## So conj, which equals transpose(m'), ie the liou of m is based on := m \rho m' = kron(conj(m),m) * vec(\rho)
liou( m::T ) where T<:AbstractMatrix = kron( conj(m), m )
using SparseArrays

# This violates a large number of programming principles, esp no repetition, but it should be clear.
# The purpose behind the repl style dots is in case you mistakenly try and print say a 1000x1000 Matrix
# At least on my computer this would hang Jupyter.
function latexArray(m)
  todisplay = 20
  if size(m)[1] < 20
    todisplay =size(m)[2]
  end
  start = string("\\left(\\begin{array}{$(join(["c" for i in 1: todisplay]))c}")
  if ndims(m) == 2
      (a,b) = size(m)
      if b > 20
        if a > 20
            for i = 1:6
              for j = 1:6
                  start = "$start $(m[i,j]) &"
              end
              if i==2
                start = "$start \\cdots &"
              else 
                start = "$start  &"
              end
              for j = (b-5):b-1
                  start = "$start $(m[i,j]) &"
              end
              for j = b:b
                  start = "$start $(m[i,j])"
              end
              start= "$start \\\\"
            end
            start = "$start & \\vdots & & & & & \\ddots & & & & & \\vdots\\\\ "
            for i = a-5:a
              for j = 1:6
                  start = "$start $(m[i,j]) &"
              end
              if i==a-1
                start = "$start \\cdots &"
              else 
                start = "$start  &"
              end
              for j = (b-5):b-1
                  start = "$start $(m[i,j]) &"
              end
              for j = b:b
                  start = "$start $(m[i,j])"
              end
              start= "$start \\\\"
            end
            start = "$start \\end{array}\\right)"
 
        else # a < 20 b > 20
            for i = 1:a
              for j = 1:6
                  start = "$start $(m[i,j]) &"
              end
              if i == 2 || i == a-1 
                start = "$start \\cdots &"
              else 
                start = "$start  &"
              end
              for j = (b-5):b-1
                  start = "$start $(m[i,j]) &"
              end
              for j = b:b
                  start = "$start $(m[i,j])"
              end
              start= "$start \\\\"
            end
            start = "$start \\end{array}\\right)"
         end
      else 
        if a > 20 # b < 20
          for i = 1:6
              for j = 1:b
                  start = "$start $(m[i,j]) &"
              end
              start= "$start \\\\"
          end
          start = "$start &  \\vdots "
                 
          for i = 1:b-3
              start = "$start & "
          end
          start = "$start \\vdots \\\\"
          for i = a-5:a
               for j = 1:b-1
                  start = "$start $(m[i,j]) &"
              end
              for j = b:b
                 start = "$start $(m[i,j]) "
              end 
              start= "$start \\\\"
            end
            start = "$start \\end{array}\\right)"
        else 
        for i = 1:size(m,1)
          for j = 1:size(m,2)-1
            start = "$start $(m[i,j]) &"
          end
          for j = size(m,2):size(m,2)
            start = "$start $(m[i,j])"
          end
          start= "$start \\\\"
        end
        start = "$start \\end{array}\\right)"
      end
    end
    return start
  end 
  return "NOT A MATRIX?" 
end 



function nicePrint(m)
  start = string("\\begin{equation*} ")
  start = start*latexArray(m)
  start = start*string("\\\\\\end{equation*}")
  return start
end

Base.show(stream, ::MIME"text/latex", x::AbstractMatrix) = write(stream, nicePrint(x))

#writemime(io, ::MIME"text/latex", x::Matrix) = write(io, nicePrint(x))

# NOTE USING the Rc involution. 1,3,2,4 or 4,2,1,3?
function choi_liou_involution( r::Matrix )
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [1,3,2,4] )
  reshape( rl, size(r) )
end


# RH changed again: this was 3, 4, 1 2 - changed to 2, 1, 4, 3 i.e. L(x \otimes y)->L(y \otimes x)
function swap_involution( r::Matrix )
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [2,1,4,3] )
  reshape( rl, size(r) )
end

"""
  choi2liou

  Takes a superoperator in choi matrix form and returns the liouville superoperator (computational basis)
"""
function choi2liou( r::Matrix  )
  choi_liou_involution( r )
end

"""
  liou2choi

  Take a liouville superoperator (computational basis) and return the choi matrix
"""
function liou2choi( r::Matrix )
  choi_liou_involution( r )
end


"""
  takes a liouville superoperator (computational basis) and returns the Chi matrix
  (basically the choi matrix in Pauli basis)
"""
function liou2choiX(r::Matrix)
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [3,1,4,2] )
  reshape( rl, size(r) )'
end

"""
  takes the Chi matrix (Choi matrix in Pauli basis) and return the liouville superoperator, computaitonal basis
"""
function choiX2liou(r::Matrix)
  d = round(Int, sqrt(size(r,1)) )
  rl = reshape( r, (d, d, d, d) )
  rl = permutedims( rl, [3,1,4,2] )
  reshape( rl, size(r) )
end

"""
  takes the choi matrix and returns the kraus operators.
   This stays the same, save that we need to take into account dimensional factors.
   Note slight negative eigenvalue (rounding can cause problems). Here I round to 15 digits.
"""
function choi2kraus(r::Matrix{T}) where T
  d = sqrt(size(r,1))
  (vals,vecs) = eigen( d*r )
  #vals = eigvals( sqrt(size(r,1))*r )
  kraus_ops = Matrix{T}[]
  for i in eachindex(vals)
    push!(kraus_ops, sqrt(round(vals[i],digits=15,RoundToZero))*mat(vecs[:,i]))
  end
  factor = tr(sum([k'*k for k in kraus_ops]))
  kraus_ops = kraus_ops/sqrt(factor/d)
  kraus_ops
end

"""
  Takes a vector of kraus operators (matrices) and 
  converts them to a liouville superoperator (computational basis)
"""
function kraus2liou( k::Vector{Matrix{T}}) where T
  l = zeros(T,map(x->x^2,size(k[1])))
  for i in 1:length(k)
    # The Adjoint is now lazy, so we need to copy it.
    l = l + liou(k[i],copy(k[i]'))
  end
  l
end

"""
  Takes a liouville superoperator (computational basis) and
  returns a vector of the equivalen Krawu vectors
"""
function liou2kraus( l::Matrix )
  choi2kraus( liou2choi( l ) )
end


""" 
  Takes a vector of kraus operators and returns the choi matrix
"""
function kraus2choi(k::Vector{Matrix{T}} ) where T
  liou2choi(kraus2liou(k))
  # Seems to fail for some kraus vectors.
  #c = zeros(T,map(x->x^2,size(k[1])))
  #for i in 1:length(k)
  #  c = vec(k[i])*vec(k[i])'
  #end
  #c/sqrt(size(c,1))
end


_translate(v) = [x==2 ? 3 : x== 3 ? 2 : x for x in v]

_Paulis = [ [1 0im;0 1],[0im 1;1 0],[0 -im;im 0],[1 0im;0 -1]]

_num2quat(n,l) = map(s->parse(Int,s),collect(string(n,base=4,pad=l)))

_toPauli(p) = reduce(kron,[_Paulis[x+1] for x in p])


"""
  Converts a liouville superoperator in the Pauli basis
  to one in the computational basis. Order of Paulis is I,X,Y and Z.
"""
function pauliliou2liou( m::Matrix )
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end
  dsq = size(m,1)
  res = zeros(ComplexF64,size(m))
  l = round(Int,log2(dsq)/2)
  for i=1:dsq
    for j=1:dsq
      res += m[i,j] * vec(_toPauli(_num2quat(i-1,l))) * vec(_toPauli(_num2quat(j-1,l)))' / sqrt(dsq)
    end
  end
  res
end

"""
  Converts a liouville superoperator in the computaional basis
  to one in the pauli basis. Order of Paulis is I,X,Y and Z.
"""
function liou2pauliliou( m::Matrix{T} ) where T
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end
  dsq = size(m,1)
  res = zeros(ComplexF64,size(m))
  l = round(Int,log2(dsq)/2)
  m = sparse(m)
  Pvecs = [sparse(vec(_toPauli(_num2quat(i-1,l)))) for i=1:dsq]
  Pvecd = [x' for x in Pvecs]
  sqr = sqrt(dsq)
  for i=1:dsq
    for j=1:dsq
      res[i,j] = Pvecd[i] * m * Pvecs[j]  / sqr
    end
  end
  return res
end


""" 
  Takes a choi matrix and returns the Chi Matrix 
  Note the basis here has changed from marcusps version,
  The rightmost paulis are varying the quickest.
"""
function choi2chi( m::Matrix{T} ) where T
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end

  dsq = size(m,1)
  dim = round(Int(log2(dsq)/2))
  pauliTranslate=vec(_toPauli(_num2quat(0,dim)))
  for i = 2:dsq
    pauliTranslate=hcat(pauliTranslate,vec(transpose(_toPauli(_num2quat(i-1,dim)))))
  end
  pauliTranslate'*m*pauliTranslate
end

function pauliliou2chi( m::Matrix{T}) where T
  choi2chi(liou2choi(pauliliou2liou(m)))
end

""" Takes a chi matrix and returns the choi Matrix 
  Note the basis here has changed from marcusps version,
  Note we normalise (divide by dsq). Some conventions differ.
  
"""
function chi2choi( m::Matrix{T} ) where T
  if size(m,1) != size(m,2)
    error("Only square matrices supported")
  elseif size(m,1) != 4^(floor(log2(size(m,1))/2))
    error("Only matrices with dimension 4^n supported.")
  end

  dsq = size(m,1)
  dim = round(Int(log2(dsq)/2))
  pauliTranslate=vec(_toPauli(_num2quat(0,dim)))
  for i = 2:dsq
    pauliTranslate=hcat(pauliTranslate,vec(transpose(_toPauli(_num2quat(i-1,dim)))))
  end
  (pauliTranslate*m*pauliTranslate')./dsq
end

#"""
#Returns a superoperator that replaces the input with a maximally
#mixed state with probability p, and leaves it unchanged with probability (1-p).
#"""
#function depol( d::Int, p=1.0 )
#  choi2liou( p * Matrix(I,d^2,d^2)/d^2 + (1-p) * projector(_max_entangled_state(d)) )
#end

"""
Given a superoperator, it extracts the closest superoperator (in Frobenius norm)
that is unital. The result may not be completely positive.
"""
function unitalproj( m::Matrix{T} ) where T
  d2 = size(m,1)
  d  = round(Int,sqrt(d2))
  id = projector(normalize(vec(eye(d))))
  id*m*id + (I-id)*m*(I-id)
end

# tweaked to increase eps slightly was getting a trivial false !cp right on the boundary.
"""
  iscp(m;tol)

  Returns whether is completely positivel.
  Assumes input is in liouville basis (not pauli-liouville - use pauliliou2liou)
  Converts to choi and checks the eigenvalues.
  -- tweaked to increase eps slightly was getting a trivial false !cp right on the boundary.
"""
function iscp(m; tol=0.0)
    evs = eigvals(liou2choi(m))
    tol = tol==0.0 ? 1.01*eps(abs.(one(eltype(m)))) : tol
    all(real(evs) .> -tol) && all(abs.(imag(evs)) .< tol)
end

"""
  istp(m;tol)

  Returns whether is trace perserving..
  Assumes input is in liouville basis (not pauli-liouville - use pauliliou2liou)
"""
function istp(m; tol=0.0)
    tol = tol==0.0 ? eps(abs.(one(eltype(m)))) : tol
    dsq = size(m,1)
    d = round(Int,sqrt(dsq))
    norm(m'*vec(eye(d))-vec(eye(d)),Inf) < tol
end

"""
  ischannel(m;tol)

  Returns whether is the liouville superoperator is trace perserving and completely positive.
  Assumes input is in liouville basis (not pauli-liouville - use pauliliou2liou)
"""
function ischannel(m; tol=0.0)
    #println(iscp(m,tol=tol))
    #println(istp(m,tol=tol))
    iscp(m,tol=tol) && istp(m,tol=tol)
end

"""
  unital(m;tol)

  Returns whether is the liouville superoperator is unital.
  Assumes input is in liouville basis (not pauli-liouville - use pauliliou2liou)
"""
function isunital(m; tol=0.0)
    tol = tol==0.0 ? eps(abs(one(eltype(m)))) : tol
    dsq = size(m,1)
    d = round(Int,sqrt(dsq))
    norm(m*vec(eye(d))-vec(eye(d)),Inf) < tol
end

"""
Computes the unitary CP map closest (interferometrically) to a given CP map.
See D. Oi, [Phys. Rev. Lett. 91, 067902 (2003)](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.91.067902)
"""
function nearestu(l)
    c = liou2choi(l)
    vals,vecs = eig(Hermitian(c))
    imax = indmax(vals)
    Λ = mat(vecs[:,imax])
    U,Σ,V = svd(Λ)
    W = U*V'
    return kron(conj(W),W)
end

"""
function makeSuper(u)

Takes an operator and turns into a superoperator. Pauli basis.
"""
function makeSuper(u)
    return round.(real(liou2pauliliou(liou(u))),digits=10)
end


# Some more stuff just to make life easier.
⊗ = kron

zbasis = [[1 0],[0 1]]
xbasis = map(x->1/sqrt(2) .* x,[[1 1],[1 -1]])
ybasis = map(x->1/sqrt(2) .* x,[[1 -im],[1 im]])

# This is an arbitrary ordering of our basis
basis_selection= [xbasis,ybasis,zbasis]


""" 
  pass in a qubit up/down state e.g. [0 0] for two qubits both in |0> state
  returns a vector of the qubits (in this case [1 0]⊗[1 0])
"""
function buildvec(ar)
    foldl(⊗,map(x->zbasis[x+1],ar))
end

"""
  pass in a qubit up/down state e.g. [0 0] for two qubits both in |0> state
  returns a density matrix of the qubits (in this case ([1 0]⊗[1 0])'*([1 0]⊗[1 0]))
"""
function buildDensity(ar)
    v = buildvec(ar)
    v'*v
end


"""
  The following are more genearlised versions of the above (which are in the computational basis)
  Here we can give an abrbitrary basis and get the vectors/density
"""
basisVec(bv) = foldl(⊗,map(x -> basis_selection[x[1]][x[2]+1],bv))
function build_basis_Density(ar)
    v = basisVec(ar)
    return v'*v
end
function build_basis_Density(basis,values)
    v = basisVec(zip(basis,values))
    return v'*v
end    

"""
    observeBasis(basis,Λ)
## Arguments
-   `basis: Array{Float64,1}` A list of the Paulis basis you want to oberve in. e.g. [1,2,1,3] would be four qubits
     representing [X,Y,X,Z]
-   `Λ: Array{Float64,2}`: The superoperator of the noise. Assumed appropriately sized (4^n x 4^n) , where n is the number
    of basis provided.


## Returns
   The 2^n probability distribution (of up/down measurements) that would be seen if you put those Paulis through the channel
   and tried to measure in the same basis.
"""
function observeBasis(basis,Λ::Array{Float64, 2}) 
    qubits = length(basis)
    ZeroBasisState = getSuperVec(build_basis_Density(basis,[0 for _ in 1:qubits]))./ √(2)^qubits
    genBasisStabs(basis,noQubits) = map(x->round.(getSuperVec(build_basis_Density(basis,x)),digits=15),[_num2bin(i,noQubits) for i=0:(2^noQubits-1)]);
    BasisStabilisers = [v ./ √(2)^qubits  for v in genBasisStabs(basis,qubits)] # √(2)^qubits as we don't normalise the supervecs by default
    obs = [V'*Λ*ZeroBasisState for V in BasisStabilisers]
    return obs
end

"""
    observeBasis(basis,Λ)
## Arguments
-   `basis: Array{Float64,1}` A list of the Paulis basis you want to oberve in. e.g. [1,2,1,3] would be four qubits
     representing [X,Y,X,Z]
-   `Λ: Array{Float64,1}`: In this case we assume a Pauli channel and the noise is represented as a vector of the Pauli eigenvalues (4^n) long.

## Returns
   The 2^n probability distribution (of up/down measurements) that would be seen if you put those Paulis through the channel
   and tried to measure in the same basis.
"""
function observeBasis(basis,Λ::Array{Float64, 1}) 
    qubits = length(basis)
    ZeroBasisState = getSuperVec(build_basis_Density(basis,[0 for _ in 1:qubits]))./ √(2)^qubits
    genBasisStabs(basis,noQubits) = map(x->round.(getSuperVec(build_basis_Density(basis,x)),digits=15),[_num2bin(i,noQubits) for i=0:(2^noQubits-1)]);
    BasisStabilisers = [v ./ √(2)^qubits  for v in genBasisStabs(basis,qubits)] # √(2)^qubits as we don't normalise the supervecs by default
    obs = [V'* (Λ .* ZeroBasisState) for V in BasisStabilisers]
    return obs
end


""" 
  pass in a density matrix, gives you the super vector that corresponds to it
"""
function getSuperVec(density)
    d = size(density,1)
    l = round(Int,log2(d))
    real(vec([(tr(_toPauli(_num2quat(i,l))*(density))) for i = 0:(d^2-1)]))
end


_num2bin(n,l) = map(s->parse(Int,s),collect(string(n,base=2,pad=l)))
"""
  returns an array of the 'measurement' operators that extract each of the z measurements
  e.g. for two qubits will return the II, IZ, ZI and ZZ measurements
"""
function genZs(noQubits)
    return map(x->getSuperVec(buildDensity(x)),[_num2bin(i,noQubits) for i=0:(2^noQubits-1)]);
end


""" 
  Helper function for labeling graphs, pass in the qubits get the string of Is and Zs
  e.g. with two you get ["00","01","10","11"]
"""
function genLabels(noQubits)
    return [string(n,base=2,pad=noQubits) for n = 0:(noQubits^2-1)];
end


"""
    measure(s,measureNoise, allZs)


  Helper function in open-systems.jl
  Pass in the state, the noise (as a superOperator) and the measurements - use allZs to generate)
  returns the mapped function (x->x'*measureNoise*s,allZs)
"""
function measure(s,measureNoise,allZs)
    map(x->x'*measureNoise*s,allZs)
end
using Combinatorics

pn=["I","X","Y","Z"]
ttn=["↑","↓"]

"""
  Helper function (chart axis)
  return all the paulis for a certain number of qubits as a list of strings
"""
 function genPauliLabels(noQubits)
    return  map(x->foldl(*,map(y->pn[parse(Int,y)+1],x)),[collect(string(i,base=4,pad=noQubits)) for i=0:(4^noQubits-1)])
end


"""
  Helper function (chart axis)
  return all the ups and downs for a certain number of qubits as a list of strings
"""
function genArrows(noQubits)
    return map(x->foldl(*,map(y->ttn[parse(Int,y)+1],x)),[collect(string(i,base=2,pad=noQubits)) for i=0:(2^noQubits-1)]) 
end

