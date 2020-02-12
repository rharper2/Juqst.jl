# Probability Distribution Commands


## Creating and manipulating

```@docs
 readInCSVFile(filename)
```




```@docs
transformToFidelity(data)
```

```@docs
fitTheFidelities(lengths,matrix;no=0)
```

```@docs
projectSimplex(ps)
```


```@docs
convertAndProject(params)
```

```@docs
getIndices(qs;dimension=16)
```



```@docs
marginalise(q,pps)
```

```@docs
gibbsRandomField(pps,constraints)
```

```@docs
gibbsRandomField(pps,constraints::Array{Array{Array{Any,1},1},1})
```


```@docs
getGrainedP(ϕ,tomatch,graining)
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
reconstructedJS(distribution,constraints)
```


```@docs
reconstruct(dist,constraints)
```

```@docs
 relEntropy(dist,constraints)
```

```@docs
getSummand(x,y,z,jpXYZ,jpXZ,jpYZ,jpZ)
```

```@docs
conditionalMutualInfo(X,Y,Z,p)
```

```@docs
function covarianceMatrix(p;reverseDigits = false)
```

```@docs
marginaliseFromRawData(rawData,constraints,lengths)
```

```@docs
correlationMatrix(p)
```


