# This julia file builds on the Initial.julia
# It is based on the "How to efficiently select an arbitrary Clifford group element"
# Robert Koenig and John Smolin
# arXivv:1406.2170v1

# Copyright Robin Harper 2015-2020


# The premise here is that in order to generate a random Clifford, we can be sure of
# proper (Haar) randomness if either we generate all the cliffords for a certain quibit size
# and select one randomly OR if we can have a 1-1 mapping between the integers and the Cliffords
# and then we can just randomly select an integer. The latter is the one adopted in Koenig's paper.
# getNumberOfCliffords(n) returns the number of cliffords with n qubits.

# This ties into the Aaronsen/Gottesman paper becauase the generated cliffords are generated
# in such a way it is possible to create the tableau that would result IF the start myVector
# |000000> had been acted on by the chosen clifford unitary.
# We can then (decompose) the tableau using the algorithms specified in the Aaronson/Gottesman paper
# (see Initial.jl) to work out what combination of one-qubit (phase and hadmard) and two quibit
# (cnot) gates would create the unitary in the first place.

# The initial part of this file is a port of the python code in the Koenig/Smolin paper
# followed by an implementation of the decomposition algorithm in the Aaronson/Gottesman paper.
# There is also code to:
# -  Draw a decomposed circuit
# -  Display the kets of a |00..0> vector transformed by the circuit
# -  Form a 2^n 2^n complex matrix representing the "raw circuit" for use outside this formalism
# -  There is brute force method for finding the minimum set of gates to construct a given unitary
#    (even with only 3qubits this takes a tediously long time to run)
# -  There is the start of rationalise method, just now it eliminates chains of 4 phase gates or 2 hadamard gates.
#     Maybe the way to go here is to eliminate recognised patterns. [REF] uses a reverse tree building approach
#     which is interesting.
#
# Use is made of the ability of Julia to manipulate its own commands (a sort of hybrid functional/imperative)
# approach. With decomposition we are using arrays stored with the Tableau
# The first "commands" is a text represenetation of the commands needed to rebuild a decomposed end state.
# the second "executeCommands" contains the actual Julia instructions to rebuild the state.
# They are only populated when the appropriate decompositions have been made.
# This is better documented in the appropriate methods.

# The generation of the Clifford requires an element from the symplectic group
#include("Initial.jl")
#using Main.CHP



# Scratch variable in the global spce of Symplectic to be used with eval.
_t = Tableau(3)


# This just shows how big the groups get, Julia overflow is going to be a problem

function getNumberOfCliffords(n::BigInt)::BigInt
		return 2^(n^2+2*n)*prod([4^x-1 for x =1:n])
end

"""
    getNumberOfCliffords(n::Integer)::BigInt

returns the number of Cliffords for n qubits.
"""
function getNumberOfCliffords(n::Integer)::BigInt
		return getNumberOfCliffords(BigInt(n))
end


"""
    getNumberOfSymplecticCliffords(n::Integer)::BigInt

returns the number of Cliffords for n qubits, modulo the signs on each row
of the tableau. Note

getNumberOfCliffords(n) == getNumberOfSymplecticCliffords(n)*getNumerOfBitStringsCliffords(n)
"""
function getNumberOfSymplecticCliffords(n::BigInt)::BigInt
	return 2^(n^2)*prod([4^x-1 for x =1:n])
end

function getNumberOfSymplecticCliffords(n::Integer)::BigInt
	return getNumberOfSymplecticCliffords(BigInt(n))
end

"""
    getNumberOfBitStringsCliffords(n)

returns the number of bitstrings representing a +/- on each row for the Tableau.
Note:

getNumberOfCliffords(n) == getNumberOfSymplecticCliffords(n)*getNumerOfBitStringsCliffords(n)
"""
function getNumberOfBitStringsCliffords(n)
	return 2^(2*n)
end

"""
    cliffordToTableau(qubits,clifford,signs)

Converts an integer into a Clifford. You need to specify the number of qubits in the
tableau, the clifford number and the bit signs corresponding to the rows of the Tableau.
The commands used to create the Clifford can be found in Tableau.commands and Tableau.executeCommands.?
It is fine to specify numbers that are too large (they wrap), but the helper functions
```Julia
getNumberOfSymplecticCliffords(n::BigInt)::BigInt
getNumberOfBitStringsCliffords(n)
```
Will show how many cliffords there are for a certain number of qubits (Warning this is a very large number for n > 2)
"""
function cliffordToTableau(qubits,clifford,signs)
    t = Tableau(qubits)
	t.state = stabiliseSymp(symplectic(clifford,qubits),signs)
	decomposeState(t)
	return t
end

"""
    tableauToClifford(t::Tableau)

Returns the 'symplectic' number that corresponds to the clifford embedded in the Tableau.
Note that this is modulo the bit signs. I.e. the tableaus corresponding to e.g.
cliffordToTableau(1,3,3) and cliffordToTableau(1,3,4) will both return 3.
"""
function tableauToClifford(t::Tableau)
	symplecticinverse(decomposeStateToSymplectic(t.state))
end

"""
    decomposeState(tableau::Tableau,rationalise=true)

Decomposes the state stabilised by the Tableau using the Aaronsen/Gottesman method

That is it will decompose the "arbitrary" state into a series of CHP (cnot, hadamard and phase) gates.
This is unlikely to be the most concise decomposition. Unless rationalise is set to
false there will be some naive trimming of gates, basically self eliminating gates (e.g. 4 phase gates)
are removed. The commands (and the execCommands) in the tableau are wiped and recreated. Following the
decomposition the state is re-created using the new commands.
"""
function decomposeState(tableau::Tableau,rationalise=true)
	global _t
	tableau.commands = []
	tableau.executeCommands = []
	#addCommand(tableau,"output(state)",Expr(:call,:show,:_t))
    ss1=tableau.state
	while getState(ss1) < 11
		nextStep(tableau)
	end
	j=tableau.qubits
	addCommand(tableau,"initialise($j)",Expr(:call,:initialise,:_t))
	reverseCommands(tableau)
	removeRedundancy(tableau)
	_t = tableau
	tracking = _t.trackCommands
	_t.trackCommands = false
	for i in tableau.executeCommands
		eval(i)
	end
	_t.trackCommands = tracking

end

function drawCircuit(t::Tableau)
	commands = t.commands
    qasmCommands=[]
    labels=[]
	currentDir = pwd()
	for i = 1:size(commands,1)
		m = match(r"initialise\((.*)\)",commands[i])
		if (m!==nothing)
			for idx=1:t.qubits
    			#push!(qasmCommands,"qubit q$i")
                push!(labels,"q$idx")
			end
    	else
    		m=match(r"hadamard\((.*)\)",commands[i])
    		if (m!==nothing)
    			push!(qasmCommands,("H","q$(m.captures[1])"))
    		else
    			m=match(r"phase\((.*)\)",commands[i])
    			if (m!==nothing)
                    push!(qasmCommands,("P","q$(m.captures[1])"))
    			else
    				m=match(r"cnot\((.*),(.*)\)",commands[i])
    				if (m !== nothing)
    					push!(qasmCommands,("CNOT","q$(m.captures[2])","q$(m.captures[1])"))
    				end
    			end
    		end
    	end
    end
    return (qasmCommands,labels)
end

function qiskitCircuit(t::Tableau,name = "Circuit")
	commands = t.commands
    qasmCommands=["#Pass in the circuit and qubit register","def create$name(circs,qreg):"]

    labels=[]
	currentDir = pwd()
	for i = 1:size(commands,1)
		m = match(r"initialise\((.*)\)",commands[i])
		if (m !== nothing)
			for idx=1:Meta.parse(m.captures[1])
    			#push!(qasmCommands,"qubit q$i")
                push!(labels,"q$idx")
			end
    	else
    		m=match(r"hadamard\((.*)\)",commands[i])
    		if (m !== nothing)
    			push!(qasmCommands,"    circs.h(qreg[$(m.captures[1])])")
    		else
    			m=match(r"phase\((.*)\)",commands[i])
    			if (m !== nothing)
                    push!(qasmCommands,"    circs.s(qreg[$(m.captures[1])])")
    			else
    				m=match(r"cnot\((.*),(.*)\)",commands[i])
    				if (m !== nothing)
    					push!(qasmCommands,"    circs.cx(qreg[$(m.captures[1])],qreg[$(m.captures[2])])")
    				end
    			end
    		end
    	end
    end
	return join(qasmCommands,"\n")
end



# using Quantikz # This is ONLY for the display functionality. If it causes problems remove the following functions.

# function quantikzCircuit(t::Tableau)
# 	commands = t.commands
#     quantikzCommands=[]

#     labels=[]
# 	currentDir = pwd()
# 	for i = 1:size(commands,1)
# 		m = match(r"initialise\((.*)\)",commands[i])
# 		if (m !== nothing)
# 			for idx=1:Meta.parse(m.captures[1])
#     			#push!(quantikzCommands,"qubit q$i")
#                 push!(labels,"q$idx")
# 			end
#     	else
#     		m=match(r"hadamard\((.*)\)",commands[i])
#     		if (m !== nothing)
#                 push!(quantikzCommands,Quantikz.H(parse(Int,m.captures[1])))
#     		else
#     			m=match(r"phase\((.*)\)",commands[i])
#     			if (m !== nothing)
#                     push!(quantikzCommands,Quantikz.P(parse(Int,m.captures[1])))
#     			else
#     				m=match(r"cnot\((.*),(.*)\)",commands[i])
#     				if (m !== nothing)
#     					push!(quantikzCommands,Quantikz.CNOT(parse(Int,m.captures[1]),parse(Int,m.captures[2])))
#     				end
#     			end
#     		end
#     	end
#     end
# 	return quantikzCommands
# end


###################### Functions used to create Clifford tabelau.



#takes square arrays and places them in a larger array, as follows
#   m1 0
#   0 m2
function directsum(m1,m2)
  n1=size(m1)[1]
  n2=size(m2)[1]
  out = zeros(n1+n2,n1+n2)
  out[1:n1,1:n1]=m1
  out[(n1+1):end,(n1+1):end]=m2
  return out
end

# Helper function used in decomposing the clifford
# builds up an array of the commands used to call the hadamard/phase/cnot.
# We are using two globals here for the sake of simplicity.
# execute commands is the Julia code required to call the relevant functions.

function addCommand(t::Tableau,comToAdd,comm)
	t.commands = append!(t.commands,[comToAdd])
	t.executeCommands=append!(t.executeCommands,[comm]);
end



# Returns the symplectic inner product of two vectors.

function inner(v,w)
	t=0
	for i in 1:(size(v)[1]>>1)
		 t+=v[2*i-1]*w[2*i]
		 t+=w[2*i-1]*v[2*i]
    end
    return t%2
 end


 function transvection(k,v)
 	if (k==0) return mod.(v,2)
 	end
 	return mod.((v+inner(k,v)*k),2)
 end

 function int2bits(i,n)
 	out = zeros(n)
 	for j in 1:n
 		out[j]=i&1
 		i>>=1
 	end
 	return out
 end

 function findtransvection(x,y)
 	out = zeros(UInt8,2,size(x)[1])
 	if (x==y)
 		return out
 	end
 	if inner(x,y) ==1
 		out[1,:]=mod.((x+y),2)
 		return out
 	end
 	z=zeros(size(x)[1])
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1] + x[ii]) != 0) && ((y[ii-1]+y[ii])!=0)
 			# found the pair
 			z[ii-1]=(x[ii-1]+y[ii-1])%2
 			z[ii]=(x[ii]+y[ii])%2
 			if (z[ii-1]+z[ii])==0 #the same
 				z[ii]=1
 				if (x[ii-1] != x[ii])
 					z[ii-1]=1
 				end
 			end
 			out[1,:]=mod.((x+z),2)
 			out[2,:]=mod.((y+z),2)
 			return out
 		end
 	end
 	# ok we didn't find a pair
 	# need places where x is 00 and y does not
 	# and vice versa.
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1]+x[ii]) !=0) && ((y[ii-1]+y[ii]) == 0)
 			if x[ii-1]==x[ii]
 				z[ii]=1
 			else
 				z[ii]=x[ii-1]
 				z[ii-1]=x[ii]
 			end
 			break
 		end
 	end
 	for i in 1:(size(x)[1]>>1)
 		ii=2*i
 		if ((x[ii-1]+x[ii]) ==0) && ((y[ii-1]+y[ii]) !=0 )
 			if y[ii-1]==y[ii]
 				z[ii]=1
 			else
 				z[ii]=y[ii-1]
 				z[ii-1]=y[ii]
 			end
 			break
 		end
 	end
 	# (x+z)%2 doesn't work on arrays, so mod function instead.
 	out[1,:]=mod.((x+z),2)
 	out[2,:]=mod.((y+z),2)
 	return out
 end

# The steps are as set out on page 5 of the paper.
# We use bigInts because we can (and if n was greater than about 64 we would need to).
function symplectic(i,n)
	nn=2*n
	s=BigInt(2)^nn-1
	k=(i%s)+1
	i=div(i,s)
	# step 2
	# There is probably a better way to do this, but it works.
	f1=[x == '0' ? UInt8(0) : UInt8(1) for x  in reverse(string(k,base=2,pad=nn))]
	#step 3
	e1=zeros(UInt8,nn)
	#define the first basis
	e1[1]=1
	T=findtransvection(e1,f1)
	# step 4
	divisor = BigInt(2)^(nn-1)
	value = i%divisor
	bits = [x == '0' ? UInt8(0) : UInt8(1) for x  in reverse(string(value,base=2,pad=nn-1))]

	#step 5
	# Note that in Julia we need to expressly copy
	eprime=copy(e1)
	for j in  3:nn
		eprime[j]=bits[j-1]
	end
	h0=transvection(T[1:1,:]',eprime)
	h0=transvection(T[2:2,:]',h0)

	#step 6
	if bits[1]==1
		f1*=0
	end

	#step 7
	id2=[1 0;0 1]

	if n!=1
		g=directsum(id2,symplectic(i>>(nn-1),n-1))
	else
		g=id2
	end
	for j in 1:nn
		g[j,:]=transvection(T[1:1,:]',g[j:j,:]')
		g[j,:]=transvection(T[2:2,:]',g[j:j,:]')
		g[j,:]=transvection(h0,g[j:j,:]')
		g[j,:]=transvection(f1,g[j:j,:]')
	end
	return g
end

function numberofcosets(n)
	x = 2^(2*n-1)*(2^(2*n)-1)
	return x
end

function bits2int(b,nn)
	output = BigInt(0)
	tmp = BigInt(1)
	for j in 1:nn
		if b[j]==1
			output = output + tmp
		end
		tmp = tmp*2
	end
	return output
end

function symplecticinverse(gn)
	n = round(Int,size(gn)[1]/2)
	nn=2*n

	# step 1

	v = gn[1,:]
	w = gn[2,:]

	# step 2

	e1 = zeros(nn)
	e1[1] = 1
	T = findtransvection(v,e1)

	# step 3
	tw = copy(w)
	tw = transvection(T[1,:],tw)
	tw = transvection(T[2,:],tw)
	b = tw[1]
	h0 = zeros(nn)
	h0[1]=1
	h0[2]=0
	for j in 3:nn
		h0[j]=tw[j]
	end

	# step 4

	bb = zeros(nn-1)
	bb[1]=b
	for j in 3:nn
		bb[j-1] = tw[j]
	end

	zv = bits2int(v,nn)-1
	zw = bits2int(bb,nn-1)

	cvw = zw*(2^(2*n)-1)+zv

	if n == 1
		return cvw
	end

	gprime = copy(gn)
	if b == 0
		for j in 1:nn
			gprime[j,:] = transvection(T[2,:],transvection(T[1,:],gn[j,:]))
			gprime[j,:] = transvection(h0,gprime[j,:])
			gprime[j,:] = transvection(e1,gprime[j,:])
		end
	else
		for j in 1:nn
			gprime[j,:]=transvection(T[2,:],transvection(T[1,:],gn[j,:]))
			gprime[j,:]=transvection(h0,gprime[j,:])

		end
	end
	gnew = gprime[3:nn,3:nn]
	gnidx = symplecticinverse(gnew)*numberofcosets(n)+cvw
	return gnidx
end




#Splits up the symplectic into the alpha,beta,gamma and delta arrays that
#specify the action of the clifford unitary on the X and Z Paulis respectively

function parseSymplectic(symp)
	#note that here we have (so far) ignored that for a,b,c,d we need to multiply an r and s
	#then I am guessing we just use rj as bit string on rhs. ie r1 to rn
	# so with one bit, we have and additional 2*2 * 6
	# then with two bits its 2^2 * 2^2 = 16
	s2 = round(Int,(size(symp)[1]/2))
	a = zeros(s2,s2)
	b = zeros(s2,s2)
	c = zeros(s2,s2)
	d = zeros(s2,s2)
	for i=1:s2
		for j = 1:s2
			a[i,j]=symp[(i-1)*2+1,(j-1)*2+1]
			b[i,j]=symp[(i-1)*2+1,j*2]
			c[i,j]=symp[i*2,(j-1)*2+1]
			d[i,j]=symp[i*2,j*2]
		end
	end
	return (a,b,c,d)
end

# Another way of doing this (incomplete - just an aide memoir)
# reshape(collect(Iterators.flatten(collect(zip(a,b)))),3,6)



function decomposeStateToαβγδ(state)
    n = size(state)[1]
    boundary = round(Int,n/2)
    return (state[1:boundary,1:boundary],state[1:boundary,boundary+1:end],
        state[boundary+1:end,1:boundary],state[boundary+1:end,boundary+1:end])
end


function decomposeStateToSymplectic(state)
    n = round(Int,size(state)[1]/2)
    finalD = n*2
    sympy=zeros(finalD,finalD)
    (α,β,γ,δ)= decomposeStateToαβγδ(state)
    for i=1:n
        for j = 1:n
            sympy[(i-1)*2+1,(j-1)*2+1]= α[i,j]
            sympy[(i-1)*2+1,j*2]=β[i,j]
            sympy[i*2,(j-1)*2+1]= γ[i,j]
            sympy[i*2,j*2]= δ[i,j]
        end
    end
    return sympy
end



# Takes the symplectic, parses it and uses it to create an Aaronson/Gottesman tableau
#The bits specify the 'bit' pattern that controls the sign
#For example in the 1 qubit there are only 6 unique 'symplectic' patterns
#but each of the two rows in the "tableau" can have bit signs of 0 or 1, leading
# to 6*4 different cliffords.
#this is what the bits is.


function stabiliseSymp(symp,bits)
	(a,b,c,d)=parseSymplectic(symp)

	n=size(a)[1]
	top = hcat(a,b,[ (bits >> (x-1)) %2 for x=1:n])
	for co=1:n
		bits = div(bits,2)
	end
  	bottom = hcat(c,d,[ (bits >> (x-1)) %2 for x=1:n])
 	state = convert(Array{Int32,2},vcat(top,bottom))
 end


# Takes tableau state and works out what state we are in the the 11 step deomposition
# where we are taking an arbitrary tableau state and "decomposing" it with basic gates
# back to the initial |00...00> state.
function getState(state)
	as = state[:,1:end]
	# calculate the states, backwards, so we get the "most refined" first
	#define some useful comparison matrices
	# The tableau is split up as follows:
	# A | B
	# ______
	# C | D
	#
	n=div(size(state,1),2) # half the dimension of this 2n x (2n+1) matrix
	state_i = Matrix{Int32}(I,n,n)
	state_z = zeros(Int32,n,n)
	state_r = zeros(Int32,n,1)
	if (as==vcat(hcat(state_i,state_z,state_r),hcat(state_z,state_i,state_r)))
		return 11
	end
	#split up into four matrices
	A=as[1:n,1:n]
	B=as[1:n,n+1:2*n]
	C=as[n+1:2*n,1:n]
	D=as[n+1:2*n,n+1:2*n]
	R=as[1:end,end]
	R2=as[n+1:end,end:end]
	if (B==state_z && C==state_z && R2==state_r)
		return 10
	end
	if (C==state_z && A==B  && R2==state_r)
		return 9
	end
	if (A==state_i && C==state_z && D==state_i && R2==state_r)
		return 7
	end
	if (C==state_i && D==state_z && as[n+1:2*n,end:end] == state_r)
		return 6
	end
	if (D==state_z )
		return 5
	end
	if (C==D)
		return 4
	end
	if (C==state_i)
		return 2
	end
	if (rank(C) == n)
		return 1
	end
	return 0
end

# The next series of vectors each implement one of the manipulations specified in Aaronson/Gottesman
# Section VI Canonical form.

function getFullRank(tableau::Tableau)
	n=tableau.qubits # half the dimension of this 2n x (2n+1) matrix
	svec = tableau.state
	for ii=1:n
		for i=1:n
			r1 = rank(svec[n+1:2*n,1:n])
			if r1 == n
				return
			end
			hadamard(tableau,i,false)
			r2 = rank(svec[n+1:2*n,1:n])
			if (r2 == n )
				#println("hadamard(",i,")")
				#addCommand(tableau,"hadamard($i)",Expr(:call,:hadamard,:_t,i))
				return
			end
			if (r2 > r1)
				#println("hadamard(",i,")",Expr(:call,:hadamard,:svec,i))
				#addCommand(tableau,"hadamard($i)",Expr(:call,:hadamard,:_t,i))
				continue
			else
				hadamard(tableau,i,false)
			end
		end
		hadamard(tableau,ii,false)
	end
end

function getfirstOne(myVector)
	#println(myVector)
	for i=eachindex(myVector)
		if myVector[i]==1
			#println("return ",i)
			return i
		end
	end
	#print("******* THIS MAY BE AN ERROR BUT GET FIRST ONE COULDNT FIND A 1********\n")
	#print("$myVector\n")
	return 0
end


function makeCtheI(tableau::Tableau,offset=-1)
	n=tableau.qubits # half the dimension of this 2n x (2n+1) matrix
	svec = tableau.state
	sanity = 100;
	if offset==-1 offset=n # no offset means we assume its the "C" we are making Int32
	end
	allDone = false;
	while !allDone
		for i=1:n
			#print("Checking for $i\n")
			allDone = identifytheRow(tableau,offset,i)
			if allDone == false
				break
			end
		end
		sanity = sanity -1;
		if (sanity < 1)
			print("Stuck in a loop making C, I\n")
			exit()
		end
	end
	# now we need double phases to make r_n+1..r_2n to be zero
end

function zapPhase(tableau::Tableau,offset=-1)
	n=tableau.qubits # half the dimension of this 2n x (2n+1) matrix
	svec = tableau.state
	if offset==-1 offset=n # no offset means we assume its the "C" we are making Int32
	end
	for i=1:n
		if (svec[offset+i,end] ==1)
			phase(tableau,i,false)
			#addCommand(t,"phase($i)",Expr(:call,:phase,:_t,i))
			phase(tableau,i,false)
			#addCommand(t,"phase($i)",Expr(:call,:phase,:_t,i))
		end
	end
end

# What I mean is take the passed in row, and use cnots to make it an identity element.
function identifytheRow(tableau::Tableau,offset,i)
	n=tableau.qubits # half the dimension of this 2n x (2n+1) state matrix
	#Make sure the diagonal element is 1, we can only use after the element to cnot with
	#println("Checking Diagonal ",i,",",i)
	svec = tableau.state
	if svec[offset+i,i] ==0 # need to get this diagonal equal to 1
		#println("It wasn't one in the diagonal we have offset ",offset, " i ", i)
		t = getfirstOne(svec[offset+i,i+1:n])
		#println(svec[offset+i,i:n])
		#println("Got a t back of ",t)
		if t==0
			#println("Couldnt make this the identity - redoing rank\n")
			getFullRank(tableau)
			return false
		end
		t=t+i
		#println("cnot(",t,",",i,")")
		cnot(tableau,t,i,false)
		#addCommand(tableau,"cnot($t,$i)",Expr(:call,:cnot,:_t,t,i))
		#println(svec[offset+i,1:n])
	end

	for j=1:i-1
		#println("Check ",i,",",j)
		if svec[offset+i,j]!= 0
			#println("Wasnt zero")
			#addCommand(tableau,"cnot($i,$j)",Expr(:call,:cnot,:_t,i,j))
			cnot(tableau,i,j,false)
			#println(svec[offset+i,1:n])
		end
	end
	#println("Up to diagonal cleared")
	#println(svec[offset+1,1:n])
	# so the diagonal element is now 1, make the rest zero
	for j=i+1:n
		#println("Checking to make zero",i,",",j)
		if svec[offset+i,j]==1
			#println("It was one")
			#println("cnot(",i,",",j,")")
			cnot(tableau,i,j,false)
			#addCommand(tableau,"cnot($i,$j)",Expr(:call,:cnot,:_t,i,j))
			#println(svec)
		end
	end
	return true
end

function diagonaliseD(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		if svec[n+i,n+i]== 0
			phase(tableau,i,false)
			#println("phase(",i,")");
			#addCommand("phase($i)",Expr(:call,:phase,:_t,i))
		end
	end
end

function diagonaliseB(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		if svec[i,n+i]== 0
			phase(tableau,i,false)
			#println("phase(",i,")");
			#addCommand(tableau,"phase($i)",Expr(:call,:phase,:_t,i))
		end
	end
end

function cmdm(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		for j=1:n
			if svec[n+i,j] == svec[n+i,n+j]
				continue
			end
			if svec[n+i,j] == 0
				cnot(tableau,j,i,false)
				#println("cnot(",j,",",i,")")
				#addCommand("cnot($j,$i)",Expr(:call,:cnot,:_t,j,i))
			else
				phase(tableau,j,false)
				#println("phase(",j,")")
				#addCommand("phase($j)",Expr(:call,:phase,:_t,i))
			end
		end
	end
end
function ambm(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		for j=1:n
			if svec[i,j] == svec[i,n+j]
				continue
			end
			if svec[i,j] == 0
				cnot(tableau,j,i,false)
				#println("cnot(",j,",",i,")")
				#addCommand(tableau,"cnot($j,$i)",Expr(:call,:cnot,:_t,j,i))
			else
				phase(tableau,j,false)
				#println("phase(",j,")")
				#addCommand(tableau,"phase($j)",Expr(:call,:phase,:_t,j))
			end
		end
	end
end

function zapD(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		for j=1:n
			if svec[n+i,n+j]==1
				phase(tableau,j,false)
				#println("phase(",j,")")
				#addCommand("phase($j)",Expr(:call,:phase,:_t,j))
			end
		end
	end
end

function zapB(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		for j=1:n
			if svec[i,n+j]==1
				phase(tableau,j,false)
				#println("phase(",j,")")
				#addCommand(tableau,"phase($j)",Expr(:call,:phase,:_t,j))
			end
		end
	end
end

function hadamardHard(tableau::Tableau)
	n=tableau.qubits
	svec = tableau.state
	for i=1:n
		hadamard(tableau,i,false)
		#println("hadamard(",i,")")
		#addCommand(tableau,"hadamard($i)",Expr(:call,:hadamard,:_t,i))
	end
end

function makeItR0(tableau::Tableau,offset=-1)
	n=tableau.qubits
	svec = tableau.state
	if offset<0
		offset=n
	end
	for i=1:n
		if svec[offset+i,2*n+1]== 1
			phase(tableau,i,false)
			phase(tableau,i,false)
			#println("phase(",i,")")
			#println("phase(",i,")")
			#addCommand(tableau,"phase($i)",Expr(:call,:phase,:_t,i))
			#addCommand(tableau,"phase($i)",Expr(:call,:phase,:_t,i))
		end
	end
end





function nextStep(tableau::Tableau)
	n=tableau.qubits # half the dimension of this 2n x (2n+1) matrix
	svec = tableau.state

	currentState = getState(svec)
	#print(currentState)
	if (currentState == 0)
		getFullRank(tableau)
	elseif currentState == 1
		makeCtheI(tableau)
	elseif currentState < 4
		diagonaliseD(tableau)
		cmdm(tableau)
	elseif currentState == 4
		zapD(tableau)
		makeItR0(tableau)
	elseif currentState == 5
		makeCtheI(tableau)
		zapPhase(tableau)
	elseif currentState == 6
		hadamardHard(tableau)
	elseif currentState < 9
		diagonaliseB(tableau)
		ambm(tableau)
	elseif currentState == 9
		zapB(tableau)
		makeItR0(tableau,0)
	elseif currentState == 10
		makeCtheI(tableau,0) # 0 = actually its A
		zapPhase(tableau,0)
	elseif currentState == 11
		println("Done")
	end
	return
end

# maxGates is the number of gates we are allowed (an integer offset to an array)
# gates is the "stack" <- an array of gates we are applying.
# increment index 1, if its greater then maxGate, increment the one above it.
function incrementGate!(gates,maxGate,offset=1)
	if (size(gates,1) < offset )
		push!(gates,1)
		return
	end
	currentValue = gates[offset]
	currentValue = currentValue +1
	if (currentValue > maxGate)
		gates[offset]=1
		incrementGate!(gates,maxGate,offset+1)
	else
		gates[offset]=currentValue
	end
end


function bruteForceBreadthFirst(clifford)
	global svec
	count = 0
	n=div(size(clifford,1),2) # half the dimension of this 2n x (2n+1) matrix
	if (n > 4)
		println("This might take some time!, there are ",n*(n+1) + n + n, " different gates , we might need ", n*n/log(n), " of them so I *might* have to check ", n*n/log(n)^(n*(n+1) + n + n), "combos !")
		println("Ctrl-C could be your friend")
	end
	gates = (Expr)[]
	svec = []
	# If we have a n quibit system then there are a total of
	# n*(n+1) possible cnots, n phase gates and n hadamards we can apply
	for i = 1:n
	  for j = 1:n
	  	if i!=j
	  		push!(gates,Expr(:call,:cnot,:svec,i,j,false))
	  	end
	  end
	  push!(gates,draw)
	  push!(gates,Expr(:call,:phase,:svec,i,false))
	end
	# so we now have an array of all possible gates we can apply to the system
	# we want to search through the gates, we "should" need at most n^2 of them.
	# The process will be to have a stack of gates to apply, check them out
	# if we match getState(11) i.e. they reverted the clifford to the initial state
	# we have found the gates needed
	# otherwise we increment our "bottom" gate, and if we have gone though them, recursively increment the
	# gate above it (incrementGate! does this for us)
	state = 0
	maxGates = size(gates,1);
	gatesToApply=(Int32)[]
	if getState(clifford) == 11
		println("Very funny you don't need to do anything")
		return
	end
	currentAt = 1
	while (state!=11 && size(gatesToApply,1) < n*n*n*n*n*n*n*n*n*n)
		svec = copy(clifford)
		incrementGate!(gatesToApply,maxGates)
		count += 1
		if (size(gatesToApply,1) > currentAt)
			print("Trying: ")
				for i=axis(gatesToApply,1)
				print(gatesToApply[i]," ")
			end
			println("")
			currentAt+=1
		end
		for index in axis(gatesToApply,1)
			#println("Going to apply ",gates[gatesToApply[index]])
			eval(gates[gatesToApply[index]])
		end
		state = getState(svec)
	end

	if (state == 11)
		println("Found it after (only) $count permuations")
		println("Gates applied")
		reverse!(gatesToApply)
		for i in axis(gatesToApply,1)
			println(gates[gatesToApply[i]])
		end
	else
		println("Didn't find it!")
	end
	return gatesToApply
end

#coding for current bits
# 1 = phase, 2 = hadamard
# 2 + x  = cnot control (with target at x)
# 2 + n + x = cnot target (with control at x)

function checkGate(tableau::Tableau,index,currentBits,n,toDelete)
	commands = tableau.commands
	executeCommands = tableau.executeCommands

	checking = commands[index]
	m = match(r"Tableau\((.*)\)",checking)
	if (m !== nothing)
    	return
    end
    m=match(r"hadamard\((.*)\)",checking)
    if (m !== nothing)
    	bit = parse(Int,(m.captures[1]))
    	if size(currentBits[bit],1) > 0
    		if currentBits[bit][end][1] == 2
    			# two hadamards make a nothing
    			(been,gone) = pop!(currentBits[bit])
    			push!(toDelete,gone)
    			push!(toDelete,index)
			else
				push!(currentBits[bit],(2,index))
			end
		else
				push!(currentBits[bit],(2,index))
		end
		return
	end
	m=match(r"phase\((.*)\)",checking)
    if (m !== nothing) # Its a phase, we need 4 of these
    	bit = parse(Int,(m.captures[1]))
    	if size(currentBits[bit],1) > 2
    		if currentBits[bit][end][1] == 1 && currentBits[bit][end-1][1] == 1 && currentBits[bit][end-2][1] == 1
    			# four phases make a nothing.
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			(b,g) = pop!(currentBits[bit])
    			push!(toDelete,g)
    			push!(toDelete,index)
			else
				push!(currentBits[bit],(1,index))
			end
		else
				push!(currentBits[bit],(1,index))
		end
		return
	end
	m=match(r"cnot\((.*),(.*)\)",checking)
    if (m !== nothing)
    		cbit = parse(Int,(m.captures[1]))
    		tbit = parse(Int,(m.captures[2]))
    		# so just now we are going to catch two cnots in a row
    		# needs to match both bits.
    		cbitNo = 2+tbit
    		tbitNo = 2+n+cbit
    		if (size(currentBits[cbit],1) > 0 && size(currentBits[tbit],1) > 0 && currentBits[cbit][end][1] == cbitNo &&
    			currentBits[tbit][end][1] == tbitNo)
    			(b,g)=pop!(currentBits[cbit])
    			push!(toDelete,g)
    			(b,g)=pop!(currentBits[tbit])
    			# dont need to push the deletion of this gate twice.
    			push!(toDelete,index)
    		else
    			push!(currentBits[cbit],(cbitNo,index))
    			push!(currentBits[tbit],(tbitNo,index))
    		end
    end
end


function removeRedundancy(tableau::Tableau)
	n = tableau.qubits
	executeCommands = tableau.executeCommands
	commands = tableau.commands
	toDelete = (Int32)[]
	bitsOf = (Array{Tuple{Int32,Int32}})[]
	for i=1:n
    	push!(bitsOf,[])
	end
	for i=1:size(commands,1)
		checkGate(tableau,i,bitsOf,n,toDelete)
	end
	deleteat!(commands,sort!(toDelete))
	deleteat!(executeCommands,toDelete)
end



##I am going to use the symplectic method to generate all the cliffords
# - Noting that we have 24, therefore there are 6 'symplectics' and 2*2=4 different bit patterns
# - Decompose, internally generates the stabilised regime and then works out what pattern of Hadamards and Phase gates are necessary to generate this, these are stored in "commands". [Here we are in the one Qubit regime - so there are no CNOTS]
# - I use a simple loop to execute the commands using a phase gate and hadamard gate. If we were to repeat this > 1 qubit, we would need to blow up the phase gates etc (for example, if we were in 3 qubit regime and it was phase(2) i.e. phase on the second qubit, we would multiply by $I\otimes P \otimes I$
# - Better still we could define our gates (here sphase, and shadamard) to have errors and get an erroneous clifford

"""
    generateRawCliffords(;phaseGate=[1 0;0 im],hadmardGate=1/sqrt(2)*[1 1;1 -1])::Array{Complex{Float64},2}[]

returns the single qubit cliffords, generated using the passed in  phase and hadamard gates (perfect gates default)
"""
function generateRawCliffords(;phaseGate=[1 0;0 im],hadamardGate=1/sqrt(2)*[1 1;1 -1])
	straightClifford = Array{Complex{Float64},2}[]
	state=[1 0;0 1]
	for i=1:6
    	for j=1:4
        	tab = cliffordToTableau(1,i,j)
        	for t in tab.commands
           		m = match(r"initialise\((.*)\)",t)
            	if (m !== nothing)
                	state=[1 0;0 1]
            	else
                	m=match(r"phase\((.*)\)",t)
                	if (m !== nothing)
           	        	state=state*phaseGate
            	    else
                    	m=match(r"hadamard\((.*)\)",t)
                    	if (m !== nothing)
                        	state=state*hadamardGate
                    	end
                	end
            	end
        	end
        	push!(straightClifford,round.(state,digits=10))
    	end
    end
    return straightClifford
end


""" 
	drawCircuit(commands)
	Takes the command string and returns their qasm equivalent.
	Can then be passed into e.g. qiskit or something similar to draw them.
"""
function drawCircuit(commands)
	qasmCommands=[]
    labels=[]
	currentDir = pwd()
	for i in axis(commands,1)
		m = match(r"setup\((.*)\)",commands[i])
		if (m !== nothing)
			for idx=1:Meta.parse(m.captures[1])
    			#push!(qasmCommands,"qubit q$i")
                push!(labels,"q$idx")
			end
    	else
    		m=match(r"hadamard\((.*)\)",commands[i])
    		if (m !== nothing)
    			push!(qasmCommands,("H","q$(m.captures[1])"))
    		else
    			m=match(r"phase\((.*)\)",commands[i])
    			if (m !== nothing)
                    push!(qasmCommands,("P","q$(m.captures[1])"))
    			else
    				m=match(r"cnot\((.*),(.*)\)",commands[i])
    				if (m !== nothing)
    					push!(qasmCommands,("CNOT","q$(m.captures[2])","q$(m.captures[1])"))
    				end
    			end
    		end
    	end
    end
    return (qasmCommands,labels)
end

"""
	makeFromCommand(command)

	Takes the commands in a Tableau and uses them to build a matrix representing the Tableau (in the computational basis)
	Uses the original initialise command to determine the size of the system.
	Note: for many qubits this gets quite big (scales as \$2^n\$)
"""
function makeFromCommand(command)
	up=[1;0]
	down=[0;1]
	control=[1 0; 0 0]
	target= [0 0;0 1]
	id=[1 0;0 1]
	x=[0 1;1 0]
	maxBits=1
    currentState = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
    si = [1 0;0 1]
    phaseG = [1 0;0 im]
    sphase1=kron([1 0;0 im],si)
    sphase2=kron(si,[1 0;0 im])
    shadamard=1/sqrt(2)*[1 1;1 -1]
    cnot12 = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
    cnot21 = [1 0 0 0;0 0 0 1;0 0 1 0;0 1 0 0]
    for t in command
        m = match(r"initialise\((.*)\)",t)
        if (m !== nothing)
        	bit =parse(Int,(m.captures[1]))
        	maxBits=bit
            state = id;
            for i=2:bit
            	state = kron(state,id);
            end
            currentState = state
        else
           m=match(r"phase\((.*)\)",t)
           if (m !== nothing)
                bit = parse(Int,(m.captures[1]))
                tphase = phaseG
                for i=1:bit-1
                	tphase = kron(id,tphase)
				end
				for i=bit+1:maxBits
					tphase = kron(tphase,id)
				end
				currentState = currentState*tphase
            else
                m=match(r"hadamard\((.*)\)",t)
                if (m !== nothing)
                    bit = parse(Int,(m.captures[1]))
                    thad = shadamard
                    for i=1:bit-1
                    	thad=kron(id,thad)
                    end
                    for i=bit+1:maxBits
                    	thad=kron(thad,id)
                    end
                    currentState=currentState*thad
                else
                    m=match(r"cnot\((.*),(.*)\)",t)
                    if (m !== nothing)
                        cbit = parse(Int,(m.captures[1]))
                        tbit = parse(Int,(m.captures[2]))

                        # just doing this each way, its easier for my mind.
                        if cbit < tbit
                        	startC=cbit
                        	endC=tbit
                        	l1 = id
                        	l2  = x
                        	for i=cbit+1:tbit-1
                        		l1=kron(id,l1)
                        		l2=kron(id,l2)
                        	end
                        	cbitbit = kron(control,l1)+kron(target,l2)
                        else
                        	startC=tbit
                        	endC=cbit
                        	l1 = id
                        	l2  = x
                        	for i=tbit+1:cbit-1
                        		l1=kron(l1,id)
                        		l2=kron(l2,id)
                        	end
                        	cbitbit = kron(l1,control)+kron(l2,target)
                        end
                        for i = 1:startC-1
                        	cbitbit=kron(id,cbitbit)
                        end
                        for i=endC+1:maxBits
                        	cbitbit=kron(cbitbit,id)
                        end
                        currentState = currentState*cbitbit
                    end
                end
            end
        end
    end
    return  currentState
end

"""
	makeFromCommand(t::Tableau)

	Uses the commands stored in a Tableau and uses them to build a matrix representing the Tableau (in the computational basis)
	Uses the original initialise command to determine the size of the system.
	Note: for many qubits this gets quite big (scales as \$2^n\$)
"""
function makeFromCommand(t::Tableau)
	return makeFromCommand(t.commands)
end

"""
	makeInverseFromCommand(t::Tableau)

	Uses the commands stored in a Tableau and uses them to build a matrix representing the inverse of the Tableau (in the computational basis)
	Uses the original initialise command to determine the size of the system.
	This might be useful for instance if you have a Tableau represnting Pauli measurements and you wanted the Clifford that maps them to the computational basis.
	Note: for many qubits this gets quite big (scales as \$2^n\$)
"""
function makeInverseFromCommand(t::Tableau)
	return makeFromCommand(reverse(t.commands[2:end]))
end


function reverseCommands(t::Tableau)
	len = length(t.commands)
	@assert length(t.executeCommands)==len
	newCommands = []
	newExecute  = []
	for i = len:-1:1
		m=match(r"phase\((.*)\)",t.commands[i])
		if (m !== nothing)
			push!(newCommands,t.commands[i])
			push!(newCommands,t.commands[i])
			push!(newCommands,t.commands[i])
			push!(newExecute,t.executeCommands[i])
			push!(newExecute,t.executeCommands[i])
			push!(newExecute,t.executeCommands[i])
        else
        	push!(newCommands,t.commands[i])
			push!(newExecute,t.executeCommands[i])
        end
    end
	t.commands = newCommands
	t.executeCommands = newExecute
    return
end
