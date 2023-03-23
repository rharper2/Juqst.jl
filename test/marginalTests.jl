@testset "Basic marginalisation" begin
    using Juqst

        test = [
        0.7648481684181914
        0.03128853038589197
        0.028471086302135486
        0.012768746026454392
        0.003189295233699477
        0.005692825953783805
        0.02086440951756285
        0.010380398979097155
        0.0012616008994507078
        0.01796575858112293
        0.004724338423184753
        0.017657967270810454
        0.02867187580848538
        0.008852146424179821
        0.022890207903508476
        0.020472643872441038 ]


        t1= [0.797970940359827 
           0.06379926134497853
           0.07695004214639156
           0.061279756148803044]
        
        t2 = [0.797970940359827
            0.07695004214639156
            0.06379926134497853
            0.061279756148803044]

        t3 = [ 0.7993051940429623
            0.07561578846325619
            0.07968100226427975
            0.04539801522950182]

        t4 = [ 0.7661097693176421
            0.03186117104218485
            0.049254288967014906
            0.014544972377963625
            0.033195424725320236
            0.043754617421071326
            0.030426713297264846
            0.03085304285153819]

        

        @test vec(marginalise([1,2],test)) == t1
        @test vec(marginalise([2,1],test)) == t2        
        @test vec(marginalise([3,1],test)) == t3
        @test vec(marginalise([3,1,2],test)) == t4
end
@testset "Basic projection, distance measures, information metrics" begin
    using Juqst
    test = [
        0.7648481684181914
        0.03128853038589197
        0.028471086302135486
        0.012768746026454392
        0.003189295233699477
        0.005692825953783805
        0.02086440951756285
        0.010380398979097155
        0.0012616008994507078
        0.01796575858112293
        0.004724338423184753
        0.017657967270810454
        0.02867187580848538
        0.008852146424179821
        0.022890207903508476
        0.020472643872441038 ]

        t5 = [0.8048635793760284
        0.029062499685849712
        0.026245055602093226
        0.01054271532641213
        0.0
        0.003466795253741543
        0.01863837881752059
        0.008154368279054892
        0.0
        0.01573972788108067
        0.0024983077231424907
        0.015431936570768192
        0.026445845108443118
        0.0
        0.020664177203466215
        0.018246613172398778]

        test2= copy(test)

        test2[5]=-0.03
        test2[14]=-0.0002
        test2[1] = 1-sum(test2[2:end])
        
        test3 = projectSimplex(test2)
        @test test3 == t5

        @test isapprox(JSD(test,test3),0.007975870587895048)
        @test isapprox(relativeEntropy(test3,test),0.024306286719978618)
        @test_logs (:warn,"D(P||Q) - the distribution (Q) has zero values not matched by by zero values in (P) - this is undefined.") relativeEntropy(test,test3)
        t6 = [ -1
        0.07868183380398325
        0.0427156565865602
        0.1053661547705379
        0.07868183380398325
       -1
        0.13395663894730434
        0.09665917504034713
        0.0427156565865602
        0.13395663894730434
       -1
        0.18387637145887897
        0.1053661547705379
        0.09665917504034713
        0.18387637145887897
       -1]

       test4= [mutualInformation(i,j,test) for i in 1:4 for j in 1:4]
       @test all(isapprox.(test4,t6))
        t7=[-1
        0.0822651509802061
        0.027222146071853907
        0.08834812333377204
        0.0822651509802061
       -1
        0.1469195222158155
        0.10151874888488262
        0.027222146071853907
        0.1469195222158155
       -1
        0.17229533207299907
        0.08834812333377204
        0.10151874888488262
        0.17229533207299907
       -1]

       @test all(isapprox.([mutualInformation(i,j,test3) for i in 1:4 for j in 1:4],t7))
       @test isapprox(conditionalMutualInfo([1,2],[4],[3],test),0.14194754202764726)
       @test isapprox(conditionalMutualInfo([1],[2],[4],test),0.03530292789372267)
       @test isapprox(conditionalMutualInfo([1],[3],[4],test), 0.049266193640797526)
       @test isapprox(conditionalMutualInfo([1],[2,3],[4],test),0.07804391541248586)
       @test isapprox(conditionalMutualInfo([1],[3,2],[4],test),0.07804391541248586)
       @test isapprox(conditionalMutualInfo([1],[3],[2,4],test),0.04274098751876294)
       @test isapprox(conditionalMutualInfo([1],[3],[4,2],test),0.04274098751876294)

       t8 = [1.0
       0.3852847467234717
       0.280484319588325
       0.4575652445864713
       0.38528474672347157
       1.0
       0.5141893352294187
       0.43136809515727703
       0.280484319588325
       0.5141893352294187
       1.0
       0.6178231101527025
       0.4575652445864713
       0.4313680951572769
       0.6178231101527025
       1.0]

       @test all(isapprox.(vec(correlationMatrix(test)),t8))

       t9=[ 0.9999999999999999
       0.4111294187181525
       0.22882547274592066
       0.43900921378758023
       0.4111294187181525
       1.0
       0.566215780260792
       0.4620411698575078
       0.22882547274592063
       0.5662157802607919
       1.0000000000000002
       0.6362755719951203
       0.43900921378758023
       0.4620411698575078
       0.6362755719951204
       0.9999999999999999]

       @test all(isapprox.(vec(correlationMatrix(test3)),t9))
end