# Probability Distribution Commands


## Reading in device data and fitting it

```@docs
 readInCSVFile
```


```@docs
transformToFidelity
```

```@docs
fitTheFidelities
```

```@docs
convertAndProject
```


## Projecting probabilities and marginalising probability distributions

```@docs
projectSimplex
```


```@docs
marginalise
```

## Measuring correlations and various other metrics

```@docs
covarianceMatrix
```


```@docs
correlationMatrix
```


```@docs
mutualInformation(p1,p2,p)
```


```@docs
 relativeEntropy(P,Q)
```


```@docs
JSD(dist1,dist2)
```

```@docs
conditionalMutualInfo(X,Y,Z,p)
```



## Gibbs Random Fields (one dimensional)


```@docs
gibbsRandomField
```

```@docs
gibbsRandomField(pps,constraints::Array{Array{Array{Any,1},1},1})
```


```@docs
getGrainedP(ϕ,tomatch,graining)
```


```@docs
reconstruct(dist,constraints)
```


```@docs
reconstructedJS(distribution,constraints)
```

```@docs
marginaliseFromRawData
```
