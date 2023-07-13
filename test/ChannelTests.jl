using Juqst

@testset "Random Channel" begin

testChannel = [1.0         0.0           0.0          0.0
               0.00101911  0.971286     -0.0706597   -0.000750742
               0.0126756   0.0711915     0.964136     0.00751068
               -0.0142784   0.000199791  -0.00735239   0.953949 ]

@test round(fidelity(testChannel),digits=10) == 0.9815618333
@test round(unitarity(testChannel),digits=10) == .9310485031
@test round(unitarityPercent(testChannel),digits=10) ==0.0475365736

testChannel = randomFidelityNoise()
@test size(testChannel) == (4,4)
@test iscp(pauliliou2liou(testChannel))
testChannel = randomPrepNoise()
@test size(testChannel) == (4,4)
@test iscp(pauliliou2liou(testChannel))
testChannel = randomMeasureNoise()
@test size(testChannel) == (4,4)
@test iscp(pauliliou2liou(testChannel))

end
@testset "Open Systems" begin

testChannel = [1.0         0.0           0.0          0.0
               0.00101911  0.971286     -0.0706597   -0.000750742
               0.0126756   0.0711915     0.964136     0.00751068
               -0.0142784   0.000199791  -0.00735239   0.953949 ]

t1 = pauliliou2liou(testChannel)
c_s = [0.969835+0.0im         9.98955e-5+0.0036762im   9.98955e-5-0.0036762im    0.0158863+0.0im
       0.000134184+0.0100931im     0.967711+0.0709256im     0.003575+0.0002659im  0.000884926+0.00258246im
       0.000134184-0.0100931im     0.003575-0.0002659im     0.967711-0.0709256im  0.000884926-0.00258246im
       0.0301647+0.0im        -9.98955e-5-0.0036762im  -9.98955e-5+0.0036762im     0.984114+0.0im  ]
@test isapprox(round.(t1,digits=6),round.(c_s,digits=6))
c_s =  [0.969835+0.0im        0.000134184-0.0100931im   9.98955e-5-0.0036762im      0.967711-0.0709256im
 0.000134184+0.0100931im    0.0301647+0.0im           0.003575+0.0002659im   -9.98955e-5+0.0036762im
  9.98955e-5+0.0036762im     0.003575-0.0002659im    0.0158863+0.0im         0.000884926-0.00258246im
    0.967711+0.0709256im  -9.98955e-5-0.0036762im  0.000884926+0.00258246im     0.984114+0.0im ]
@test isapprox(round.(liou2choi(t1),digits=6),round.(c_s,digits=6))
c_s =   [0.969835-0.0im         9.98955e-5-0.0036762im   0.000134184-0.0100931im     0.967711-0.0709256im
         9.98955e-5+0.0036762im    0.0158863-0.0im            0.003575-0.0002659im  0.000884926-0.00258246im
         0.000134184+0.0100931im     0.003575+0.0002659im     0.0301647-0.0im        -9.98955e-5+0.0036762im
         0.967711+0.0709256im  0.000884926+0.00258246im  -9.98955e-5-0.0036762im     0.984114-0.0im ]
@test isapprox(round.(liou2choiX(t1),digits=6),round.(c_s,digits=6))
@test isapprox(round.(testChannel,digits=6),round.(liou2pauliliou(choi2liou(liou2choi(pauliliou2liou(testChannel)))),digits=6))

c_s =   [3.88937+0.0im            0.00101911-0.0148631im   -0.0126756+0.000950533im    -0.0142784+0.141851im
        0.00101911+0.0148631im        0.053201+0.0im         -0.0005318-0.0142784im    -0.000550951+0.0126756im
        -0.0126756-0.000950533im    -0.0005318+0.0142784im     0.038901+0.0im           -0.00015829+0.00101911im
        -0.0142784-0.141851im     -0.000550951-0.0126756im  -0.00015829-0.00101911im       0.018527+0.0im  ]

@test isapprox(round.(c_s,digits=5),round.(choi2chi(liou2choi(pauliliou2liou(testChannel))),digits=5))
testChannel = liou2choi(pauliliou2liou(kron(randomFidelityNoise(),randomFidelityNoise())))
@test all(isapprox.(testChannel,chi2choi(choi2chi(testChannel))))

end
