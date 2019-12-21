
@testset "Basic CHP operations" begin
		
			using Juqst

           t = Tableau(4)
           start = Tableau(4)
           @test size(t.state)==(8,9)
           @test t.state ==  [1  0  0  0  0  0  0  0  0
							 0  1  0  0  0  0  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  0  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  0  0  0  0  1  0  0
							 0  0  0  0  0  0  0  1  0]
           hadamard(t,2)
           @test t.state == [1  0  0  0  0  0  0  0  0
 							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  0  0  0  0
							 0  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  0
							 0  0  0  0  0  0  0  1  0]
		   phase(t,1)
		   @test t.state == [1  0  0  0  1  0  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  0  0  0  0
							 0  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  0
							 0  0  0  0  0  0  0  1  0]
			cnot(t,2,1)
			@test t.state == [1  0  0  0  1  1  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  1  0  0  0
							 1  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  0
							 0  0  0  0  0  0  0  1  0]
			X(t,3)
			@test t.state == [ 1  0  0  0  1  1  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  1  0  0  0
							 1  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  1
							 0  0  0  0  0  0  0  1  0]
			Z(t,4)
			@test t.state == [1  0  0  0  1  1  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  1
							 0  0  0  0  1  1  0  0  0
							 1  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  1
							 0  0  0  0  0  0  0  1  0]
			Z(t,3)
			@test t.state == [1  0  0  0  1  1  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  1  0  0  0  0  0  1
							 0  0  0  1  0  0  0  0  1
							 0  0  0  0  1  1  0  0  0
							 1  1  0  0  0  0  0  0  0
							 0  0  0  0  0  0  1  0  1
							 0  0  0  0  0  0  0  1  0]
			@test kets(t) == "+|0010> \n+|1110> \n"
			@test t.commands == ["initialise(t)"
								 "hadamard(2)"  
								 "phase(1)"     
								 "cnot(2,1)"    
								 "hadamard(3)"  
								 "phase(3)"     
								 "phase(3)"     
								 "hadamard(3)"  
								 "phase(4)"     
								 "phase(4)"     
								 "phase(3)"     
								 "phase(3)" ]
			t.trackCommands=false
			hadamard(t,4)
			@test t.commands == ["initialise(t)"
								 "hadamard(2)"  
								 "phase(1)"     
								 "cnot(2,1)"    
								 "hadamard(3)"  
								 "phase(3)"     
								 "phase(3)"     
								 "hadamard(3)"  
								 "phase(4)"     
								 "phase(4)"     
								 "phase(3)"     
								 "phase(3)" ]
			
			initialise(t)
			@test t.state ==  [1  0  0  0  0  0  0  0  0
							 0  1  0  0  0  0  0  0  0
							 0  0  1  0  0  0  0  0  0
							 0  0  0  1  0  0  0  0  0
							 0  0  0  0  1  0  0  0  0
							 0  0  0  0  0  1  0  0  0
							 0  0  0  0  0  0  1  0  0
							 0  0  0  0  0  0  0  1  0]
           @test t.commands == ["initialise(t)"
								 "hadamard(2)"  
								 "phase(1)"     
								 "cnot(2,1)"    
								 "hadamard(3)"  
								 "phase(3)"     
								 "phase(3)"     
								 "hadamard(3)"  
								 "phase(4)"     
								 "phase(4)"     
								 "phase(3)"     
								 "phase(3)" ]
			reinitialise(t)
			@test t.commands == ["initialise(t)"]

							    

    end;


    @testset "Basic Symplectic operations" begin
		
 			t = cliffordToTableau(3,3472,3)
 			@test t.state == [0  0  0  0  1  0  1
						 0  0  0  1  0  1  1
						 1  0  1  1  0  1  0
						 0  1  0  1  0  1  0
						 1  0  0  0  1  1  0
						 0  0  0  0  0  1  0]
			hadamard(t,3)
			cnot(t,2,1)
			cnot(t,1,2)
			cnot(t,1,3)
			@test tableauToClifford(t) == 1280169
			@test decomposeState(t)
			@test t.commands == ["initialise(3)"
									"cnot(3,2)"
									"phase(2)"
									"phase(2)"
									"phase(2)"
									"phase(3)"
									"phase(3)"
									"cnot(3,2)"
									"phase(3)"
									"phase(3)"
									"phase(3)"
									"phase(2)"
									"phase(2)"
									"phase(2)"
									"phase(1)"
									"phase(1)"
									"hadamard(3)"
									"hadamard(2)"
									"hadamard(1)"
									"cnot(2,1)"
									"phase(1)"
									"cnot(2,1)"
									"phase(2)"
									"phase(2)"
									"phase(2)"
									"phase(1)"
									"phase(1)"
									"phase(1)"
									"cnot(2,1)"]
			@test drawCircuit(t) == (Any[("CNOT", "q2", "q3"), ("P", "q2"), ("P", "q2"), ("P", "q2"), ("P", "q3"), ("P", "q3"), ("CNOT", "q2", "q3"), ("P", "q3"), ("P", "q3"), ("P", "q3"), ("P", "q2"), ("P", "q2"), ("P", "q2"), ("P", "q1"), ("P", "q1"), ("H", "q3"), ("H", "q2"), ("H", "q1"), ("CNOT", "q1", "q2"), ("P", "q1"), ("CNOT", "q1", "q2"), ("P", "q2"), ("P", "q2"), ("P", "q2"), ("P", "q1"), ("P", "q1"), ("P", "q1"), ("CNOT", "q1", "q2")], Any["q1", "q2", "q3"])
			@test qiskitCircuit(t) == "#Pass in the circuit and qubit register\ndef createCircuit(circs,qreg):\n    circs.cx(qreg[3],qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[3])\n    circs.s(qreg[3])\n    circs.cx(qreg[3],qreg[2])\n    circs.s(qreg[3])\n    circs.s(qreg[3])\n    circs.s(qreg[3])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[1])\n    circs.s(qreg[1])\n    circs.h(qreg[3])\n    circs.h(qreg[2])\n    circs.h(qreg[1])\n    circs.cx(qreg[2],qreg[1])\n    circs.s(qreg[1])\n    circs.cx(qreg[2],qreg[1])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[2])\n    circs.s(qreg[1])\n    circs.s(qreg[1])\n    circs.s(qreg[1])\n    circs.cx(qreg[2],qreg[1])"
			@test getNumberOfCliffords(3) == 92897280
			@test getNumberOfSymplecticCliffords(4) == 47377612800
			@test getNumberOfBitStringsCliffords(5) == 1024
			@test vec(generateRawCliffords(hadamardGate=1/sqrt(3).*([1 -1;-1 1]))[3] ) == vec([0.5773502692 + 0.0im 0.5773502692 + 0.0im 0.5773502692 + 0.0im  0.5773502692 + 0.0im])
			state = makeFromCommand(t.commands)
			


			end;