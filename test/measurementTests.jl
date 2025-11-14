using Juqst

@testset "Distribution Measurement Tests" begin


qubits=3
x = sort(rand(4^qubits+1)).* 0.2
p = [x[i+1]- x[i] for i in 1:4^qubits ]
p[1]=1-sum(p)
projectSimplex(p)
p[1]+=1-sum(p)
@test isapprox(sum(p),1) # sums to one
@test filter(x->x<=0,p) == [] # no -ve or zero

# Eigenvalues


using LinearAlgebra

⊗ = kron
# full distribution, we actually want our Hadamard.
pT1 = [1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]
pT = foldl(⊗,[pT1 for _ in 1:qubits])
# round eigenvalues down to roughly machine precision
e = round.(pT*p,digits=14);
# check we haven't rounded too much.
@test all(isapprox.(1/4^qubits .* pT*e,p))

noise = diagm(e)

obs = observeBasis([3 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
pz_index = [foldl((x,y)->x<<2+y,map(x->x*3,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[pz_index[i]+1],ob_eigenvalues[i]) 
end

obs = observeBasis([2 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
py_index = [foldl((x,y)->x<<2+y,map(x->x*2,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[py_index[i]+1],ob_eigenvalues[i]) 
end

obs = observeBasis([1 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
px_index = [foldl((x,y)->x<<2+y,map(x->x*1,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[px_index[i]+1],ob_eigenvalues[i]) 
end


for _ in 1:10
  rbasis = rand(1:3,qubits)
  translate(basis,vector) = [x > 0 ? basis[ix] : 0 for (ix,x) in enumerate(vector)]

  obs = observeBasis(rbasis,noise)
  obse = observeBasis(rbasis,e)
  @test all(isapprox.(obs,obse))
  ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
  # generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
  p_index = [foldl((x,y)->x<<2+y,translate(rbasis,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
  for i in 1:2^qubits
      @test isapprox(e[p_index[i]+1],ob_eigenvalues[i]) 
  end
end


qubits=3
x = sort(rand(4^qubits+1)).* 0.2
p = [x[i+1]- x[i] for i in 1:4^qubits ]
p[1]=1-sum(p)
projectSimplex(p)
p[1]+=1-sum(p)
@test isapprox(sum(p),1) # sums to one
@test filter(x->x<=0,p) == [] # no -ve or zero

# Eigenvalues

using LinearAlgebra

⊗ = kron
# full distribution, we actually want our Hadamard.
pT1 = [1 1 1 1;1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]
pT = foldl(⊗,[pT1 for _ in 1:qubits])
# round eigenvalues down to roughly machine precision
e = round.(pT*p,digits=14);
# check we haven't rounded too much.
@test all(isapprox.(1/4^qubits .* pT*e,p))

noise = diagm(e)

obs = observeBasis([3 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
pz_index = [foldl((x,y)->x<<2+y,map(x->x*3,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[pz_index[i]+1],ob_eigenvalues[i]) 
end

obs = observeBasis([2 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
py_index = [foldl((x,y)->x<<2+y,map(x->x*2,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[py_index[i]+1],ob_eigenvalues[i]) 
end

obs = observeBasis([1 for _ in 1:qubits],noise)
ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
# generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
px_index = [foldl((x,y)->x<<2+y,map(x->x*1,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
for i in 1:2^qubits
    @test isapprox(e[px_index[i]+1],ob_eigenvalues[i]) 
end


for _ in 1:10
  rbasis = rand(1:3,qubits)
  translate(basis,vector) = [x > 0 ? basis[ix] : 0 for (ix,x) in enumerate(vector)]

  obs = observeBasis(rbasis,noise)
  ob_eigenvalues = ifwht_natural(obs); # Note the standard WHT as we are dealing with commuting observations.
  # generate indexes of Z Paulis, e..g IIII,IIIZ,IIZI,IIZZ etc
  p_index = [foldl((x,y)->x<<2+y,translate(rbasis,Juqst._num2bin(n,qubits))) for n in 0:2^qubits-1] 
  for i in 1:2^qubits
      @test isapprox(e[p_index[i]+1],ob_eigenvalues[i]) 
  end
end

end
