# This file contains all the basic functions needed to run a quantum simulator
# using the stabiliser formalism
# The format used is as described in
# Improved Simulation of Stabilizer Circuits
# Scott Aaronson and Daniel Gottesman
# arXiv:quant-ph/0406196v5

# It is, effectively, a Julia Port of the CHP program by Scott Aaronson
# The original CHP program can be found at http://www.scottaaronson.com/chp

# Primarily, I have tried to write for clarity rather than efficiency.
# (Although thanks to the algorithms by Aaronson/Gottesman its fast enough for hundreds of qubits)

# I am partially through rationalising some things for ease of REPL
# Currently  you

# include("Initial.l")
# using Main.CHP

# Then you can e.g.
# t = Tableau(4)
# hadamard(t,2,4) # Applies hadamard to qubits 2 and 4
# cnot(t,2,1)
# cnot(t,2,3)
# cnot(t,4,3)
# measure(t,1,2)
# print(kets(t))

# Copyright Robin Harper 2015-2020


"""
    Tableau(n::Integer)

 sets up and returns the a CHP Tableau for n qubits.

 This is based on the formalism given by:
 *Improved Simulation of Stabilizer Circuits*,
 Scott Aaronson and Daniel Gottesman,
 arXiv:quant-ph/0406196v5

 The initial tableau represents a |00...0\$\\rangle\$ ket in the stabiliser state
 This stabilises with "Z" and anti-stabilises with "X"
 So the tableau starts off with an Identity in the 'anti commute' top half in the X section 
 and identity in the commuting (bottom) half in the Z section.

 For the purposes of this port, the tableau is exactly replicated as per the paper
 i.e. the "state" (Tableau.state) is an Int32 array (used as a bit array)
 containing the following information.

```
 x11   .....  x1n | z11   ...      z1n | r1
  .    \\       .  |  .    \\          . |  .
  .     \\      .  |  .     \\         . |  .     Destabilisers
  .      \\     .  |  .      \\        . |  .
 xn1      \\   xnn | zn1      \\      znn| rn
 ______________________________________________
 x(n+1)1. x(n+1)n | z(n+1) ... z(n+1)n | r(n+1)
  .    \\       .  |  .      \\        . |  .
  .     \\      .  |  .       \\       . |  .     Stabilisers
  .      \\     .  |  .        \\      . |  .
 x(2n)1   \\x(2n)n | z(2n)1     \\ z(2n)n| r(2n)
```
Set Tableau.showRaw to true to see the underlying state as a matrix.


See also: [`cnot`](@ref), [`hadamard`](@ref), [`phase`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)

Related: [`cliffordToTableau`](@ref)
"""
function Tableau(x::Integer)
    t  = Tableau(setup(x),x,false,["initialise($x)"],[],true)
    t.executeCommands = [Expr(:call,:initialise,:_t)]
    return t
end


"""
    Tableau(cs::Array{String,1})

    creates a tableau from a string of commands.
    String parsing is not very robust, it assumes they are in exactly 
    the form stored internally in a tableau.
    e.g

    >t = Tableau(2)
    >hadamard(t,2)
    >phase(t,1)

    > t.commands
    3-element Array{String,1}:
     "initialise(2)"
     "hadamard(2)"  
     "phase(1)" 


     >x=Tableau(t.commands)
"""
function Tableau(cs::Array{String,1})
    m = match(r"initialise\((.*)\)",cs[1])
    if (m === nothing)
        @warn("Command string needs to start with an intialise(::Integer) command")
        return nothing
    end
    bits =parse(Int,(m.captures[1]))
    t = Tableau(Juqst.setup(bits),bits,false,["initialise($bits)"],[],true)
    t.executeCommands = [Expr(:call,:initialise,:_t)]
    for c in cs[2:end]
        m=match(r"phase\((.*)\)",c)
        if (m !== nothing)
           bit = parse(Int,(m.captures[1]))
           phase(t,bit)
        else
            m=match(r"hadamard\((.*)\)",c)
            if (m !== nothing)
               bit = parse(Int,(m.captures[1]))
               hadamard(t,bit)
            else
               m=match(r"cnot\((.*),(.*)\)",c)
               if (m !== nothing)
                  cbit = parse(Int,(m.captures[1]))
                  tbit = parse(Int,(m.captures[2]))
                  cnot(t,cbit,tbit)
                else
                    m=match(r"measure\((.*)\)",c)
                    if (m !== nothing)
                        bit = parse(Int,(m.captures[1]))
                        measure(t,bit)
                    end
                end

            end
        end
    end
    return t
end



"""
    reinitialise(t::Tableau)

Reinitialise the Tableau, setting it support the |000000> state.
Resets all commands associated with the Tableau.

See also: [`initialise`](@ref)
"""
function reinitialise(t::Tableau)
    initialise(t)
    t.trackCommands = true
    t.commands = ["initialise($(t.qubits))"]
    t.executeCommands = [Expr(:call,:initialise,:_t)]
    return t
end

function initialize(t::Tableau)
    initialise(t)
end


"""
    initialise(t::Tableau)

Initialises the Tableau, setting it support the |000000> state.
However, unlike reinitialse it does not change the commands associated with the Tableau.

See also: [`reinitialise`](@ref)
"""
function initialise(t::Tableau)
    t.state = setup(t.qubits)
    return t
end


Base.show(io::IO,z::Tableau) = print(io,pformat(z))
Base.show(stream, ::MIME"text/latex", x::Tableau) = write(stream, nicePrintTableauState(x))



"""
Used by `show` to 'text' print out the Tableau based on the state.
"""
function pformat(tab::Tableau)::String
    # prints the Heisnberg representation of the state array supplied
    state = tab.state

    toOutput::String = "\n"
    n=div(size(state,1),2)
    #  1 in a and a+n = Y
    #  else 1 in a = X
    #  else 1 in a+n=Z
    #  otherwise I
    #  last bit is sign. (exponent of -1).
    for row=1:n
      if state[row,2*n+1] == 1
        toOutput *= "-"
      else
        toOutput *= "+"
      end
      for col=1:n
         toOutput *= decode(state[row,col],state[row,col+n])
      end
      toOutput *= "\n"
    end
    for dot=1:n+1
        toOutput *= "-"
    end
    toOutput *= "\n"
    for row=n+1:2*n
      if state[row,2*n+1] == 1
          toOutput *= "-"
      else
          toOutput *= "+"
      end
      for col=1:n
          toOutput *= decode(state[row,col],state[row,col+n])
      end
      toOutput *= "\n"
    end
    return toOutput
end




"""
    cnot(t::Tableau,control::Integer,target::Integer,showOutput=false)

Performs a CNOT (controlled not) on the tableau, with the specified qubits
as control and target.

If showOutput is passed as true the state of the tableau subsequent to the cnot
is printed.

See also: [`hadamard`](@ref), [`phase`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function cnot(t::Tableau,control::Integer,target::Integer,showOutput=false)
  # From control a to taget b
  # You can supress the output by calling with the fourth parameter false
  # otherwise it will "output" the stabiliser state.
  # note the state is changed, (in Julia this should really be cnot!)
  state = t.state
  n = t.qubits
  endC=2*n+1

  if control < 1 || control > n
    print("Control qubit $(control) out of range of state which has $(n) qubits\n")
    throw(BoundsError())
    return
   end
  if target < 1 || target > n
    print("Target qubit $(target) out of range of state which has $(n) qubits\n")
    throw(BoundsError())
    return
  end
  if t.trackCommands
    append!(t.commands,["cnot($control,$target)"])
    append!(t.executeCommands,[Expr(:call,:cnot,:_t,control,target)])
  end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,control]
      xib = state[i,target]
      zia = state[i,control+n]
      zib = state[i,target+n]
      state[i,endC] = j_xor(ri,xia*zib*(j_xor(xib,j_xor(zia,1))))
      state[i,target]=j_xor(xib,xia)
      state[i,control+n]=j_xor(zia,zib)
      #print("a is $(a) b is $(b), xia is $(xia) xib is $(xib), setting $([i,b]) to $(state[i,b])\n")
   end
   if showOutput
       print(t)
   end
end

"""
    hadamard(t::Tableau,qubit::Integer,showOutput::Bool=false)

Performs a hadamard gate on the specified qubit of the tableau.

See also: [`cnot`](@ref), [`phase`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function hadamard(t::Tableau,qubit::Integer,showOutput::Bool=false)
  # Apply a hadamard gate to quibit a in the state, state.
  # again note it changes the state
  state = t.state
  n = t.qubits
  endC=2*n+1

  if qubit < 0 || qubit > n
    print("Qubit $(qubit) out of range of state (tableau only has $(n) qubits)\n")
    throw(BoundsError())
    return
  end
  if t.trackCommands
    append!(t.commands,["hadamard($qubit)"])
    append!(t.executeCommands,[Expr(:call,:hadamard,:_t,qubit)])
  end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,qubit]
      zia = state[i,qubit+n]
      state[i,endC]=j_xor(ri,xia*zia)
      state[i,qubit]=zia
      state[i,qubit+n]=xia
   end
   if showOutput
      print(t)
   end
end

"""
    hadamard(t::Tableau,qubits::Integer...;showOutput::Bool = false)

Peforms hadamards on **each** of the qubits supplied.

# Examples
```julia-repl
julia> hadamard(t,2,3,1)
```

See also: [`cnot`](@ref) [`phase`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function hadamard(t::Tableau,qubits::Integer...;showOutput::Bool = false)
    for i in qubits
        hadamard(t,i,false)
    end
    if showOutput
       print(t)
    end
end




"""
    phase(t::Tableau,qubit::Integer,showOutput::Bool=false)

Performs a phase gate on the specified qubit of the tableau.

See also: [`cnot`](@ref), [`hadamard`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function phase(t::Tableau,qubit::Integer,showOutput::Bool=false)
  # Apply the phase gate to qubit a
  # again changes the state (phase!)
  state = t.state
  n = t.qubits
  endC=2*n+1

  if qubit < 0 || qubit > n
    print("Qubit $(qubit) out of range of state (tableau only has $(n) qubits)\n")
    throw(BoundsError())
    return
   end
  # so we need two loops!
  if t.trackCommands
      append!(t.commands,["phase($qubit)"])
      append!(t.executeCommands,[Expr(:call,:phase,:_t,qubit)])
  end
  for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,qubit]
      zia = state[i,qubit+n]
      state[i,endC]=j_xor(ri,xia*zia)
   end
   for i = 1:(2*n)
      ri=state[i,endC]
      xia = state[i,qubit]
      zia = state[i,qubit+n]
      state[i,qubit+n]=j_xor(xia,zia)
   end
   if showOutput
      print(t)
   end
end


"""
    phase(t::Tableau,qubits::Integer...;showOutput::Bool = false)

Peforms phase on each of the qubits supplied.

# Examples
```julia-repl
julia> phase(t,2,3,1)
```

See also: [`cnot`](@ref) [`hadamard`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function phase(t::Tableau,qubits::Integer...;showOutput::Bool = false)
    for i in qubits
        phase(t,i,false)
    end
    if showOutput
       print(t)
    end
end



"""
    Z(t::Tableau,qubit::Integer,showOutput::Bool=false)

Convenience function to perform a 'Z' gate on the tableau. Performs two phase gates.

See also: [`cnot`](@ref), [`hadamard`](@ref), [`measure`](@ref), [`X`](@ref), [`phase`](@ref)
"""
function Z(t::Tableau,qubit::Integer,showOutput::Bool=false)
    phase(t,qubit,false)
    phase(t,qubit,false)
    if showOutput
       print(t)
    end
end

"""
    Z(t::Tableau,qubits::Integer...;showOutput::Bool = false)

Convenience function to perform a 'Z' gate on the tableau on each of the qubits supplied. Performs the gate by doing two phase gates.

# Examples
```julia-repl
julia> Z(t,2,3,1)
```

See also: [`cnot`](@ref) [`hadamard`](@ref), [`measure`](@ref), [`X`](@ref), [`Z`](@ref)
"""
function Z(t::Tableau,qubits::Integer...;showOutput::Bool = false)
    for i in qubits
        Z(t,i)
    end
    if showOutput
       print(t)
    end
end

"""
    X(t::Tableau,qubit::Integer,showOutput::Bool=false)

Convenience function to perform an 'X' gate on the tableau. Performs hadamard, two phase gates and a hadamard.

See also: [`cnot`](@ref), [`hadamard`](@ref), [`measure`](@ref), [`Z`](@ref), [`phase`](@ref)
"""
function X(t::Tableau,qubit::Integer,showOutput::Bool=false)
    hadamard(t,qubit,false)
    phase(t,qubit,false)
    phase(t,qubit,false)
    hadamard(t,qubit,false)
    if showOutput
       print(t)
    end
end

"""
    X(t::Tableau,qubits::Integer...;showOutput::Bool = false)

Convenience function to perform an 'X' gate on the tableau on each of the qubits supplied.
Performs the gate by doing a hadamard, two phase gates and a hadamard on each of the qubits.

# Examples
```julia-repl
julia> X(t,2,3,1)
```

See also: [`cnot`](@ref) [`hadamard`](@ref), [`measure`](@ref), [`phase`](@ref), [`Z`](@ref)
"""
function X(t::Tableau,qubits::Integer...;showOutput::Bool = false)
    for i in qubits
        X(t,i,false)
    end
    if showOutput
       print(t)
    end
end


"""
    measure(t::Tableau,qubit::Integer)

Performs a measurement on the supplied qubit. If the Tableau is not in a Z eigenstate, it
will (pseudo) randomly choose the result and collapse the Tableau approriately.

"""
function measure(t::Tableau,qubit::Integer)
  state = t.state
  n=t.qubits
  endC=2*n+1

  if qubit < 0 || qubit > n
     print("Qubit $(qubit) out of range of state (tableau only has $(n) qubits)\n")
     throw(BoundsError())
     return
  end
  if t.trackCommands
      append!(t.commands,["measure($qubit)"])
      append!(t.executeCommands,[Expr(:call,:measure,:_t,qubit)])
  end
  if sum(state[(n+1):end,qubit]) > 0
     # we have a one on the xia somewhere.
     # This makes it a random measurement
     p = n+argmax(state[(n+1):end,qubit]); # Find the first stabiliser with an X in this qubit.
     # set p-nth to equal p, (because of the measurement the stabiliser
     # has now become a destabiliser. I.e. it will anti-commute with the measurement, but
     # commutes with everything else).
     state[p-n,:]=state[p,:]
     # Now set up the measurment 'stabiliser' this can be just IIs with a Z in the appropriate
     # qubit (and the 'r' set to 1 or 0 depending on the random outcome)
     value = round.(Integer,rand()) # check distribution (did we roll a 1 or a 0)
     state[p,:]=zeros(2*n+1) # All of them zero
     state[p,endC]=value # Set the end approrpiately.
     state[p,n+qubit]=1 # Set the 'z' qubit 
     # Finally we need to make everything commute with the measurment.
     # we do this my right multiplying everything that doesn't commute with the de-stabiliser.
     for h=1:2n # h is what we use in rowmult!
       if (h!=p-n) && (state[h,qubit]==1) # We have an X (or Y) in the relevant qubit
         rowmult!(state,h,p-n)
       end
     end
    return value # the result of the measurement.
  else
      # print("Deterministic\n\n")
       #first set the 2n+1 row i.e. scratch space to be zero
       scratch=zeros(Int32,2n+1)' # hey 2n works dont you know
       state = vcat(state,scratch) # add it on.
       for i = 1:n
         if state[i,qubit]==1
           rowmult!(state,2n+1,i+n)
         end
       end
       value = state[2n+1,endC]
       state=state[1:end-1,:]
       return value
   end
end




"""
    measure(t::Tableau,qubits::Integer...;showOutput = false)

Performs a measurement on the supplied qubits. Results are reported in a dictionary of outcomes.

# Examples
```julia-repl
julia> using Main.CHP

julia> t = Tableau(4)

+XIII
+IXII
+IIXI
+IIIX
-----
+ZIII
+IZII
+IIZI
+IIIZ


julia> measure(t,1,2)
Dict{Any,Any}(2=>0,1=>0)
julia>
```
"""
function measure(t::Tableau,qubits::Integer...;showOutput = false)
    outcomes = Dict()
    for i in qubits
        outcomes[i] = measure(t,i)
    end
    if showOutput
        print(outcomes)
    end
    if showOutput
       print(t)
    end
    return outcomes
end


"""
    kets(tableau::Tableau)

returns a string showing the kets of the state stabilised by the tableau.
"""
function kets(originalTableau::Tableau)

    #Use gaussian Elimination to provide a minimal set of generators
    # containing X's and Y's in upper triangular form
    ss1 = gaussianElimination(originalTableau.state)

    state=ss1[1]
    n=ss1[2]

    start=seed(state,n)
    scratchState=copy(start) # the first nqubits are the destabilisers
    outputString = printBasisState(vec(scratchState))

    # So the idea is we have the stabilisers with an X in them in upper triangular form
    # Each X represents a "1" in the Ket
    # But we need to generate each of the states stabilised by these Paulis.
    # For the first row, the ones generated are where there is an X in the first row.
    # On the second row, the ones generated are the second row, and the second row multiplied by the first.
    # On the third row, its the third row, third times first, third times second and third times first and second
    # etc.
    # Aaronsen has a really neat algorithm to do this, but its not easily understood.
    # I have used this, slighly less efficient, but (I hope) more comprehenisible algorithm
    # For a particular row (say we are looking at our third row with an X in it (n 3)
    # loop through the numbers 4:7 (being 2^(3-1) to (2^3)-1)
    # In binary, this loops as 100, 101, 110, 111, So we can see This is applying
    # Row 3 with each possible combination of Rows 1 and 2.
    # Then we need to multiply the rows together if 2^(row-1) (i.e. row 2 - equates to bit 2) is not zero when
    # logically anded with the count.

    number=2^(n-1)
    nqubits=div(size(state,1),2)
    for count = 0:(2^(n-1)-2)
        c2=xor(count,(count+1))
        for rows=0:(n-2)
            if (2^(rows) & c2) > 0
               scratchState=rowmult(scratchState,state[nqubits+rows+1,:])
            end
        end
        # Note the use of vec here, just makes sure that the scratch state is a vectdor
        # This was only an issue for the |00> one.
        outputString *= printBasisState(vec(scratchState))
    end
    return outputString
end


"""
    nicePrintTableauState(t::Tableau)

Used by show to print out the latex equation of the current tableau for e.g. Jupyter notebook
    In theory, if you want the latex equation call this and post into tex!

"""
function nicePrintTableauState(t::Tableau)
  state = t.state
  n = t.qubits
  toOuptut = ""
  if !t.showRaw
      start = string("\\begin{equation*} \\left(\\begin{array}{l|$(join(["c" for i in 1:n]))}")
      for row=1:n
          if state[row,2*n+1] == 1
            toOutput = "-&"
          else
            toOutput = "+&"
          end
          start = "$start $toOutput &"
          for col=1:n
             start = "$start $(decode(state[row,col],state[row,col+n])) "
         if (col != n)
            start = "$start & "
          end
          end
          start = "$start \\\\\n"
        end
        start = "$start \\hline\n"
        for row=n+1:2*n
          if state[row,2*n+1] == 1
            toOutput = "-&"
          else
            toOutput = "+&"
          end
          start = "$start $toOutput &"
          for col=1:n
            start = "$start $(decode(state[row,col],state[row,col+n])) "
            if (col != n)
              start = "$start & "
            end

          end
          start = "$start \\\\\n"
        end
        start = "$start \\end{array}\\right)\\\\\\end{equation*}"
      return start
    else
        start = string("\\begin{equation*} \\begin{array}{l$(join(["c" for i =1:n]))|$(join(["c" for i =1:n]))|c}")
        toshift = round(Integer,n/2,RoundUp)
        start = "$start\n$(join(["&" for i =1:toshift]))\\textbf{X}$(join(["&" for i =1:(n-toshift)]))$(join(["&" for i =1:toshift]))\\textbf{Z}$(join(["&" for i =1:(n-toshift)]))\\\\\n"
        for row=1:n
          start = "$start D$(row): &"
          for col=1:(2n)
             start = "$start $(state[row,col]) &"
          end
          start = "$start $(state[row,2n+1]) \\\\\n"
        end
        start = "$start \\hline\n"
        for row=n+1:2*n
            start = "$start S$(row-n): &"
             for col=1:(2n)
             start = "$start $(state[row,col]) &"
          end
          start = "$start $(state[row,2n+1])\\\\\n"
        end
        start = "$start \\end{array}\\\\\\end{equation*}"
        return start
    end
end

# ============== Internal Functions ==========================
function setup(n)
  # sets up and returns the initial state vector for n qubits.
  # Heisenberg representation
  # We are representing a |00...0> ket in the stabiliser state
  # This stabilisises with "Z" and anti-stabilises with "X"
  # So the tableau starts off with an Identity in the Anti commute top X section
  # and identity in the commuting - Z section.

  # For the purposes of this port, the tableau is exactly replicated as drawn
  # i.e. the "state" is an Int32 array (should really be just binary)
  # of the form
  #
  #
  # x11   .....  x1n | z11   ...      z1n | r1
  #  .    \       .  |  .    \          . |  .
  #  .     \      .  |  .     \         . |  .     Destabilisers
  #  .      \     .  |  .      \        . |  .
  # xn1      \   xnn | zn1      \      znn| rn
  # ______________________________________________
  # x(n+1)1. x(n+1)n | z(n+1) ... z(n+1)n | r(n+1)
  #  .    \       .  |  .      \        . |  .
  #  .     \      .  |  .       \       . |  .     Stabilisers
  #  .      \     .  |  .        \      . |  .
  # x(2n)1   \x(2n)n | z(2n)1     \ z(2n)n| r(2n)

  top = hcat(Matrix{Int32}(I,n,n),zeros(Int32,n,n),zeros(Int32,n,1))
  bottom = hcat(zeros(Int32,n,n),Matrix{Int32}(I,n,n),zeros(Int32,n,1))
  state = vcat(top,bottom)
end

# This is just mod 2 arithmetic.
function j_xor(a,b)
  return mod.(a+b,2)
end




function rowswap!(alteredState,i,k)
      temp = alteredState[i,:]
      alteredState[i,:]=alteredState[k,:]
      alteredState[k,:]=temp
end

""" 
    cliffordPhase(i,k)

    Pass in two vectors representing rows in a Tableau.

    Returns the exponent i is raised to (0,1,2 or 3) if i is left multiplied by k.
    ---
    Note: Assumes that the 'r' part of i (in position 2n+1) has ALREADY been changed from
    an exponent of -1 (which is how it is in the Tableau) to an exponent of i
    (by multiplying it by 2). It also assumes k is still in the range of 0..1 i.e. an
    exponent of -1
"""
function cliffordPhase(i,k)
  # Returns the phase (0,1,2,3) when row i is left multiplied by row k
  # println("CliffordPhase $i,$k")
  e=0
  n=div(length(i),2)
  for j=1:n
    if k[j] == 1 && k[j+n]==0 # its an X
      if i[j] == 1 && i[j+n] == 1 # its a Y
        e+=1 # XY=iZ
      elseif   i[j] == 0 && i[j+n] == 1 # its a Z
        e+=-1 # XZ=-iY
      end
    elseif  k[j] == 1 && k[j+n]==1 # its a Y
      if i[j] == 0 && i[j+n] == 1       # its a Z
        e+=1 # YZ=iXZ
      elseif i[j] == 1 && i[j+n] == 0 # its a X
        e+=-1 # YX=-iZ
      end
    elseif  k[j] == 0 && k[j+n]==1 # its a Z
      if i[j] == 1 && i[j+n] == 0       # its a X
        e+=1 # ZX=iY
      elseif i[j] == 1 && i[j+n] == 1 # its a Y
        e+=-1 # ZY=-iX
      end
    end
  end

  e = e + i[2*n+1] + 2*k[2*n+1]
  e = e%4
  if (e < 0)
    e+=4
  end
  
  return e
end

"""
    rowmult!(stateArray,i,k)

    stateArray - the state of a Tableau
    i the (integer) row of the stateArray that will hold the result of i left multiplied by 
      row pointed to by k.
    k the (integer) row with which i is to be left multiplied.

    Note: assumes that both 'r' components of the Tableau rows are exponents of -1
    use _rowmult! if i has already been lifted up to be an exponent of i.
""" 
function rowmult!(stateArray,i,k)
    n=div(size(stateArray,1),2)
    # Note cliffordPhase assumes the first row has the 'r' component changed to represent
    # a phase of i, we do this here by multiplying it by 2
    # and then after the row multiplication dividing it by 2.
    stateArray[i,2n+1] *=2
    _rowmult!(stateArray,i,k)
    stateArray[i,2n+1] = div(stateArray[i,2n+1],2)
    
end

function _rowmult!(stateArray,i,k)
    n=div(size(stateArray,1),2)
    # print(stateArray[i,:],stateArray[k,:])
    phaseBack = cliffordPhase(stateArray[i,:],stateArray[k,:])
    if (phaseBack  == 1 || phaseBack  ==3)
        println("phaseBack $phaseBack WAS 1 OR 3")
    end
    stateArray[i,:] = xor.(stateArray[i,:],stateArray[k,:])
    stateArray[i,2*n+1] = phaseBack
end



"""
    rowmult(i,k)
    
    i,k two vectors, equivalent to rwos in a Tableau.
    Left multiplies i with k and returns the result.
    Note it assumes that the 'r' in i is an exponent of i
    Whereas the 'r' in k is an exponent of -1.
"""
function rowmult(i,k) # supply two vectors and get the multiplied vector out.
  n=div(length(i),2)
  temp = cliffordPhase(i,k)
  tempState = [xor(i[c],k[c]) for c =eachindex(i)] 
  tempState[length(i)]=temp
  return tempState
 end

#returns a tableau state and a size where gaussian elimination has
#been carried out on the supplied state, in order to
#provide a minimal set of generators containing X's and Y's in
#upper triangular form

function gaussianElimination(state)
  alteredState = copy(state)
  #println(alteredState)
  n=div(size(alteredState,1),2)
  i = n+1 # first row of commuting generators
  k=0
  k2=0
  j=0
  for j=1:n
    found = 0;

    for kidx=i:2*n
      # Find a generator containing X in the jth column
      if alteredState[kidx,j] == 1
        found = 1
        k = kidx
        break
      end
    end
    if found  == 1
      #swap row with the row pointed to by i
      if (i!=k)
        rowswap!(alteredState,i,k)
        rowswap!(alteredState,i-n,k-n) # swap the non-commutators as well
      end
      # then use the row to eliminate X from that bit in the rest of the tableau
      for k2=i+1:2*n
        if alteredState[k2,j] == 1
          #println("Before, $(k2) and j is $j and i is $i")
          #println(alteredState)
          #print("ROWMULTA $k,$i\n")
          rowmult!(alteredState,k2,i)
          #println(alteredState)
          #println("Doing $(i-n) and $(k2-n)")
          #print("ROWMULTB $(i-n) $(k2-n)\n")
          rowmult!(alteredState,i-n,k2-n)
          #println("Row mult $k2, $i $alteredState")
        end
        #println(alteredState)
      end
      #println(alteredState)
      i+=1
    end
  end
  gen = i - n
  # The first gen generators are X/Ys and in quasi upper triangular form.
  for j=1:n
    found = 0
    for kidx=i:2*n
      if alteredState[kidx,n+j] == 1 # we have found a Z in the jth bit
        k = kidx
        found = 1
        break
      end
    end
    if found == 1
      if (i!=k)
        rowswap!(alteredState,i,k)
        rowswap!(alteredState,i-n,k-n)
      end
      for k2=i+1:2*n
        if (alteredState[k2,j+n] ==1 ) # z in this "bit"
          #print("ROWMULTC $k2,$i\n")

          rowmult!(alteredState,k2,i)
          #print("ROWMULTD $(i-n),$(k2-n)\n")

          rowmult!(alteredState,i-n,k2-n)
        end
      end
      i+=1
    end

  end
  return (alteredState,gen)
end


function printBasisState(sS)
    toPrint = ""
    nqubits=div(size(sS,1),2)
    #println("\npBS: $sS and $nqubits and $(size(sS)) e=$(sS[2*nqubits+1])\n")
    e=sS[2*nqubits+1]
    for i=1:nqubits
      if sS[i]==1 && sS[nqubits+i]==1 # its a Y
        e = (e+1) % 4
      end
    end
    if (e==0) toPrint *= "+|"
    elseif e== 1 toPrint *= "+i|"
    elseif e==2 toPrint *= "-|"
    else toPrint *= "-i|"
    end
    for i=1:nqubits
      if sS[i]==1
        toPrint *= "1"
      else
        toPrint *= "0"
      end
    end
    toPrint *= "> \n"
    return toPrint
end

function seed(state,n)
  nqubits=div(size(state,1),2)
  scratchState = zeros(Int32,1,2*nqubits+1)
   for i=2*nqubits:-1:(nqubits+n)
    min = nqubits
    f = 2*state[i,2*nqubits+1]
    # count down through the qubits.
    for j=nqubits:-1:1
      # if the z part of this is 1
      if state[i,nqubits+j] == 1
        min = j
        if scratchState[1,j]==1
          f = (f+2)%4
        end
      end
    end
    if f==2
      #print("F was 2 so we need to flip the x");
      scratchState[1,min] = xor(scratchState[1,min],1)

    end
  end
  return scratchState
end



function decode(z,x)
    if z==1
        if x == 1
            return "Y"
        else
            return "X"
        end
    elseif x == 1
        return "Z"
    else
        return "I"
    end
end
