"""
    Current assumption is that the pairs are passed in so the first qubit is on top of or
    or to the left of the second qubit

    Designed to be called by hinton plots and plot the qubit to qubit interactions.
"""
function outlineQubits(qbitPairs,ax)
    colors = ["darkslategray","slategray"]
    for x in qbitPairs
        firstP = translate_14Q_toIndex(x[1])
        secondP = translate_14Q_toIndex(x[2])
        if firstP[1]<secondP[1] # To the left
            rect = plt.Rectangle([firstP[1]+0.55,firstP[2]-0.45], 1.9, 0.9, facecolor="slategray",edgecolor="darkslategray",linewidth=2)
        else # To the bottom
            rect = plt.Rectangle([firstP[1]+0.55,firstP[2]-0.45], 0.9, 1.9, facecolor="slategray",edgecolor="darkslategray",linewidth=2)
        end    
        ax.add_patch(rect)
        #rect = []
        #if x[1] < 6
        #    rect = plt[:Rectangle]([x[1]+0.55,0.55], 1.9, 0.9, facecolor="slategray",edgecolor="darkslategray",linewidth=2)
        #else 
        #    rect = plt[:Rectangle]([(13-x[1])+1.55,1.55], 1.9, 0.9, facecolor="slategray",edgecolor="darkslategray",linewidth=2)
        #end
        #ax[:add_patch](rect)
    end
end


"""
Adapted from demo algorithm found at 
https://matplotlib.org/gallery/specialty_plots/hinton_demo.html

No default on ax or max_weight - pass it in (allows different plots to keep it constant)

Modified to accept keyword, which allows for a fixed highlight of negative values.

## Arguments:
-   `matrix: Array{Float64,1}` of values to be plotted

-   `max_weight: Float` The value required to 'fill a box'

-   `ax: PyCall.PyObject`, the graphical axis (get from gca())

## Optional Arguments:


-   `addScale: Bool` If true will draw the scale of the box (from max_weight) to the rhs of the plot.
-   `highlightNegative: Bool` If true will print negative values as a qubit, value given by qubit and font by fontsize.
-   `fontsize: Float64` default is 12.
-   `qubit: Int` qubit to highlight i.e. the qubit wrt which all the other information is printed.

    
"""
function hinton(matrix,max_weight, ax;
                        errorsH=[],
                        errorsL=[],
                        qubitPairs=[],
                        highlightNegative = false,
                        fontsize = 12,
                        qubit=0,
                        addScale=false,
                        xoffset=0,
                        yoffset=0,
                        joined=[])
    #"""Draw Hinton diagram for visualizing a weight matrix."""
    if errorsH == []
        errorsH = deepcopy(matrix)
    end
    if errorsL == []
        errorsL = deepcopy(matrix)
    end
    matrix=transpose(matrix)
    errorsH=transpose(errorsH)
    errorsL = transpose(errorsL)
    ax.patch.set_facecolor("gray")
    ax.set_aspect("equal", "box")
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            color="lightgray"
            size = 1.0
            if matrix[x,y]!= -3
                rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                             facecolor="gray", edgecolor="lightgray")
                ax.add_patch(rect)
            end
        end
    end
    if addScale
        ax.add_patch(rect)
        rect = plt.Rectangle([size(matrix)[1]+xoffset+0.5,yoffset+1.3],0.1,0.1,
                facecolor = "gray")
    
        ax.add_patch(rect)
    
        plt.annotate(s="", xy=(size(matrix)[1]-1+0.7,1.5), xytext=(size(matrix)[1]-1+0.7,0.5), arrowprops=Dict([(:arrowstyle,"<->")]))
        plt.annotate( string(max_weight), 
                    xy=(size(matrix)[1]-1+0.8,0.8),
                    xytext=(size(matrix)[1]-1+0.8,0.8), 
                    va = "top", ha="left")

    end

    joined=[]
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            w= matrix[x,y]
            if w == -1
                 push!(joined,(x,y))
            end
        end
    end

  
    outlineQubits(qubitPairs,ax)
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            w= matrix[x,y]
            
            barsH = errorsH[x,y]
            barsL = errorsL[x,y]
            
            size = min(sqrt(abs(w) / max_weight),1)
            largestError = 0
            smallestError = 0
            if w < 0
                largestError = min(sqrt(abs(barsL)/max_weight),1)
                smallestError = min(sqrt(abs(barsH)/max_weight),1)
                if barsH > 0
                    smallestError = 0
                end
            else  
                largestError = min(sqrt(abs(barsH)/max_weight),1)
                smallestError = min(sqrt(abs(barsL)/max_weight),1)
                if barsL < 0 
                    smallestError = 0
                end
            end
            if largestError < size
                    #print("CHECK ($x,$y): largestError = $largestError, size = $size\n")
                largestError = size
            end
            if smallestError > size
                    #print("CHECK ($x,$y): smallestError = $smallestError, size = $size\n")
                smallestError = size
            end
            if w >= 0 
                
                color = "lightgray" 
                rect = plt.Rectangle([x - largestError / 2, y - largestError / 2], largestError, largestError,
                                 facecolor=color, edgecolor= color)
                ax.add_patch(rect)

                rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                            facecolor=color, edgecolor= "black",linestyle=":")
                ax.add_patch(rect)
                color = "white" 
                rect = plt.Rectangle([x - smallestError / 2, y - smallestError / 2], smallestError, smallestError,
                                            facecolor="white", edgecolor= "white")
                ax.add_patch(rect)


            else 
                if highlightNegative || w == -3 || w == -1
                    if w == -3
            #            rect = plt[:Rectangle]([x - size / 2, y - size / 2], size, size,
            #            facecolor="gray", edgecolor="gray")
            #            ax[:add_patch](rect)
    
                    elseif w == -2
                        color = "red"
                        plt.text(x,y,"X",color="red",fontsize=fontsize,
                                horizontalalignment="center",
                                verticalalignment="center"
                                )
                    else
                    
                        if length(joined) == 2 # need to generalise if we are putting in more than one 'pair'
                            (a,b),(c,d) = joined
                            if (b==d) # join on same plane
                                rect = plt.Rectangle([a ,b - 0.35],1,0.62,facecolor="blue",edgecolor="black")
                                ax.add_patch(rect)
                            else
                                rect = plt.Rectangle([a-0.33,b],0.62,1,facecolor="blue",edgecolor="black")
                                ax.add_patch(rect)
                            end
                        end


                        color = "blue"
                        size=0.7
                        circ = plt.Circle([x - 0.01  , y - 0.03  ], size/2,
                             facecolor=color, edgecolor= "black")
                        ax.add_patch(circ)
                        q=99
                        if y==1 
                            q = x-1
                        else
                            q = (15-x)
                        end
                        plt.text(x,y,"Q"*string(q),color="white",fontsize=fontsize,
                                horizontalalignment="center",
                                verticalalignment="center"
                                )
                        push!(joined,(x,y))
                    end
                else # just a plain old negative number
                    color = "darkslategrey" 
                    rect = plt.Rectangle([x - largestError / 2, y - largestError / 2], largestError, largestError,
                                 facecolor=color, edgecolor= color)
                    ax.add_patch(rect)

                    rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                            facecolor=color, edgecolor= "white",linestyle=":")
                    ax.add_patch(rect)
                        
                    color = "black" 
                    rect = plt.Rectangle([x - smallestError / 2, y - smallestError / 2], smallestError, smallestError,
                                            facecolor="black", edgecolor= "black")
                    ax.add_patch(rect)
                end
            end
            
        end
    end
    ax.autoscale_view()
    ax.invert_yaxis()
    
end

"""
    translateLocationd(data,qubit)

Rather specific IBMQX16 function. Basically remaps the qubits into the same order as they appear in the IBM 16 Qubit machine.
"""
function translateLocation(data,i)
    a = Array{Float64}(undef,2,8)
    a[2,1]=data[1]
    a[1,1]=data[2]
    a[1,2:8]= data[3:9]
    a[2,2:8]=reverse(data[10:16])
    
    return a
end




"""
    translateLocationd(data,qubit)

Rather specific IBMQX14 function. Basically remaps the qubits into the same order as they appear in the IBM 14 Qubit machine.
"""
function translate_14Q_Location(data,i)
    a = Array{Float64}(undef,2,8)
    a[1,1:7]= data[1:7]
    a[2,2:8]=reverse(data[8:14])
    a[2,1]=-3
    a[1,8]=-3
    return a
end


"""
Specific to a particular device. Translates the qubit to an x,y grid reference.
"""
function translate_14Q_toIndex(i)
   if i < 7
        return (i,1)
   end
   return (14-i,2)
end


""" 
    drawIBMHinton(titleText)

Another specific function (specific to Melbourne). Basically draws hinton plots for each of the qubits, given the mutualInformation Ps.
If you specify a save file it will save it.

"""
function drawIBM_14Q_Hinton(mutualPs,titleText;saveFile = "",titleFontSize=24)
    addScale = false
    fig = figure("Slightly larger",figsize=(18,6))
    scale = round(maximum(maximum.(mutualPs)/5),RoundUp,digits=3)
    fig[:suptitle](titleText*" Full box = $scale",fontsize=titleFontSize)
    ax = gca()
    ax[:set_facecolor]("gray")
    for (ix,i) in enumerate(1:4)
      subplot(4,4,ix)
      ax=gca()
     hinton(translate_14Q_Location(mutualPs[i],i),scale,ax,highlightNegative=true,qubit=i-1)
    end

    for (ix,i) in enumerate(5:7)
        subplot(4,4,4+ix)
        ax=gca()
        hinton(translate_14Q_Location(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1)
    end
    for (ix,i) in enumerate(14:-1:12)
        subplot(4,4,9+ix)
        ax=gca()
        hinton(translate_14Q_Location(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
            fontsize=i > 10 ? 9 : 12)
    end
    for (ix,i) in enumerate(11:-1:8)
        subplot(4,4,12+ix)
        ax=gca()
        hinton(translate_14Q_Location(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
            fontsize=i >10 ? 9 : 12)
    end
    fig[:subplots_adjust](top=1.3)
    plt[:tight_layout]()
    if saveFile != ""
        savefig(saveFile)
    end
end





""" 
    drawIBMHinton(titleText)

Another specific function (specific to IBMQX16). Basically draws hinton plots for each of the qubits, given the mutualInformation Ps.
If you specify a save file it will save it.

"""
function drawIBMHinton(mutualPs,titleText;saveFile = "",titleFontSize=24)
    addScale = false
    fig = figure("Slightly larger",figsize=(18,6))
    scale = round(maximum(maximum.(mutualPs)),RoundUp,digits=3)
    fig[:suptitle](titleText*" Full box = $scale",fontsize=titleFontSize)
    ax = gca()
    ax[:set_facecolor]("gray")
    for (ix,i) in enumerate(2:4)
      subplot(4,4,ix)
      ax=gca()
     hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1)
    end


    for (ix,i) in enumerate(5:5)
     subplot(4,4,4)
      ax=gca()
      hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,addScale=addScale)
    end
    for (ix,i) in enumerate(6:9)
        subplot(4,4,4+ix)
        ax=gca()
        hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1)
    end
    for (ix,i) in enumerate(vcat([1],16:-1:14))
        subplot(4,4,8+ix)
        ax=gca()
        hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
            fontsize=i > 10 ? 9 : 12)
    end
    for (ix,i) in enumerate(13:-1:10)
        subplot(4,4,12+ix)
        ax=gca()
        hinton(translateLocation(mutualPs[i]',i),scale,ax,highlightNegative=true,qubit=i-1,
            fontsize=i >10 ? 9 : 12)
    end
    fig[:subplots_adjust](top=1.3)
    plt[:tight_layout]()
    if saveFile != ""
        savefig(saveFile)
    end
end


"""
    Draw Hinton diagram for visualizing a weight matrix.
    
    Adapted from demo algorithm found at 
    https://matplotlib.org/gallery/specialty_plots/hinton_demo.html

    No default on ax or max_weight - pass it in (allows different plots to keep it constant)

    We appear to have suffered from parameter creep! There are a large number of adjustable parameters now, all with semi-sensible defaults

## Arguments:
-   `matrix: Array{Float64,1}` of values to be plotted
-   `max_weight: Float` The value required to 'fill a box'
-   `ax`, the graphical axis (get from gca())
-   `highCorr/lowCorr : Array{Float64,1}`- possible high/low values. Used for drawing error boxes.
-   `highlightNegative`: If true will print negative values as a qubit, value given by qubit and font by fontsize.
-   `addScale`: If true will draw the scale of the box (from max_weight) to the rhs of the plot.
-    `fontsize = 12`
-    `qubit=0`
-    `addScale=false`
-    `addAxis=false`
-    `addOneAxis=false`
-    `adjust = 0.5`
-    `noQ=false`
-    `ind=[]`
-    `showQubits=false`
-    `qlabels=[]`
-    `stagger=true`
-    `voffset=0`
-    `useBars=false`
-    `thickness=0.01`
-    `altColor = "gray"`
    
"""
function covhinton(matrix, max_weight, ax;
    highCorr= [],
    lowCorr=[],
    highlightNegative = false,
    fontsize = 12,
    qubit=0,
    addScale=false,
    addAxis=false,
    addOneAxis=false,
    adjust = 0.5,
    noQ=false,
    ind=[],
    showQubits=false,
    qlabels=[],
    stagger=true,
    voffset=0,
    useBars=false,
    thickness=0.01,
    altColor = "gray")
    matrix=transpose(matrix)
    if highCorr ==[]
        highCorr = matrix
    else
        highCorr=transpose(highCorr)
    end
    if lowCorr==[]
        lowCorr = matrix
    else
        lowCorr = transpose(lowCorr)
    end
    maxX = size(matrix)[1]
    maxY = size(matrix)[2]
    xoffset=10
    yoffset=10
    ax.axis("off")

    ax.patch.set_facecolor("lightgray")
    ax.set_aspect("equal", "box")
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())
    #xlim([0,size(matrix[1])+2])
    xlim([0,maxX+3])
    ylim([0,maxY+3])
    rect = plt.Rectangle([1.45,1.45],maxX+0.1,maxY+0.1,facecolor="black",edgecolor="black")
    ax.add_patch(rect)

    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            if x > y 
                color = altColor
            else
                color="gray"
            end
            size = 1.0
            rect = plt.Rectangle([x+1 - size / 2, y+1 - size / 2], size, size,
               facecolor=color, edgecolor="lightgray")
            ax.add_patch(rect)
        end
    end
    if addAxis
        if ind==[]
            ind = collect(0:13)
        end
        qstring = "Q"
        if noQ
            qstring=""
        end
        if qlabels == []
            qlabels=["$qstring$(ind[x])" for x =1:maxX]
        end
        for x in 1:maxX
            plt.annotate( "$(qlabels[x])", 
                xy=(2,0.2),
                xytext=(x+adjust,0.4+(stagger ? x%2 : 0)*0.4+voffset), 
                va = "top", ha="left",fontsize=fontsize)
            plt.annotate( "$(qlabels[x])", 
                xy=(1.5,0.2),
                xytext=(adjust-0.2,x+adjust+0.2), 
                va = "top", ha="left",fontsize=fontsize)
            if !(addOneAxis)
                plt.annotate( "$(qlabels[x])", 
                    xy=(2,0.2),
                    xytext=(maxX+1.7+adjust,x+adjust+0.2), 
                    va = "top", ha="left",fontsize=fontsize)
                plt.annotate( "$(qlabels[x])", 
                    xy=(2,0.2),
                    xytext=(x+adjust,maxY+2+(stagger ? x%2 : 0)*0.4+voffset), 
                    va = "top", ha="left",fontsize=fontsize)
            end
        end

    end
    if addScale
        ax.add_patch(rect)
        rect = plt.Rectangle([size(matrix)[1]+xoffset+0.5,yoffset+1.3],0.1,0.1,
            facecolor = "blue")
        
        ax.add_patch(rect)
        
        plt.annotate(s="", xy=(size(matrix)[1]+1.7,2.5), xytext=(size(matrix)[1]+1.7,1.5), arrowprops=Dict([(:arrowstyle,"<->")]))
        plt.annotate( string(max_weight), 
            xy=(size(matrix)[1]+1.8,1.8),
            xytext=(size(matrix)[1]+1.8,1.8), fontsize=fontsize,
            va = "top", ha="left")

    end
    for x in 1:size(matrix)[1]
        for y in 1:size(matrix)[2]
            w= matrix[x,y]
            if isapprox(w,1) && showQubits # assume its the qubit
                 color = "blue"
                size=0.7
                circ = plt.Circle([x+1 - 0.01  , y+1 - 0.03  ], size/2,
                   facecolor=color, edgecolor= "black")
                ax.add_patch(circ)
                q=x-1
                plt.text(x+1,y+1,"Q"*string(q),color="white",fontsize=fontsize,
                    horizontalalignment="center",
                    verticalalignment="center"
                   )

            else
                size = min(sqrt(abs(w) / max_weight),1)
                maxw = highCorr[x,y]
                minw = lowCorr[x,y]
                largestError = 0
                smallestError = 0
                if w < 0
                    largestError = min(sqrt(abs(minw)/max_weight),1)
                    smallestError = min(sqrt(abs(maxw)/max_weight),1)
                    if maxw > 0
                        if useBars
                            smallestError = -1*smallestError
                        else
                            smallestError = 0
                        end
                    end
                else  
                    largestError = min(sqrt(abs(maxw)/max_weight),1)
                    smallestError = min(sqrt(abs(minw)/max_weight),1)
                    if minw < 0 
                        if useBars
                            smallestError = -1*smallestError
                        else
                            smallestError = 0
                        end
                    end
                end
                if largestError < size
                    #print("CHECK ($x,$y): largestError = $largestError, size = $size\n")
                    largestError = size
                end
                if smallestError > size
                    #print("CHECK ($x,$y): smallestError = $smallestError, size = $size\n")
                    smallestError = size
                end
                        

                if w >= 0 
                    if !useBars
                        # large square light gray
                        # actual square light gray, black border
                        # small sqaure (if !=0) white
                            
                            
                        color = "lightgray" 
                        rect = plt.Rectangle([x+1 - largestError / 2, y+1 - largestError / 2], largestError, largestError,
                                 facecolor=color, edgecolor= color)
                        ax.add_patch(rect)

                        rect = plt.Rectangle([x+1 - size / 2, y+1 - size / 2], size, size,
                        facecolor=color, edgecolor= "black",linestyle=":")
                        ax.add_patch(rect)
                        color = "white" 
                        rect = plt.Rectangle([x+1 - smallestError / 2, y+1 - smallestError / 2], smallestError, smallestError,
                                            facecolor="white", edgecolor= "white")
                        ax.add_patch(rect)
                    else
                        color = "white" 
                        rect = plt.Rectangle([x+1 - size / 2, y+1 - size / 2], size, size,
                        facecolor=color, edgecolor= color)
                        ax.add_patch(rect)
                        # plot the error bar
                        if x!=y # not on the diagonals
                            color = "red" 
                            moresize = 0
                            if smallestError < 0 
                               moresize = size/2
                            end
                            rect = plt.Rectangle([x+1 + size/3 , y+1- smallestError/2],  thickness,-moresize-largestError/2+smallestError/2,
                                facecolor=color, edgecolor=color)
                            ax.add_patch(rect)
                #print("($x,$y): Size is $(size/2), going from $(-smallestError/2) to $(-moresize-largestError/2+smallestError/2)\n")
                        end
                    end 

                else 
                    if highlightNegative
                        color = "blue"
                        size=0.7
                        circ = plt.Circle([x+1 - 0.01  , y+1 - 0.03  ], size/2,
                        facecolor=color, edgecolor= "black")
                        ax.add_patch(circ)
                        if x==y
                            plt.text(x+1,y+1,"Q"*string(y-1),color="white",fontsize=fontsize,
                             horizontalalignment="center",
                             verticalalignment="center"
                            )
                        end
                        rect = plt.Rectangle([x+0.5- 0.01, y+0.5 - 0.01], 0.98,0.98,
                                 facecolor="slategrey", edgecolor= "darkslategrey")
                            ax.add_patch(rect)

                    else
                         if !useBars
                        # large square dim gray
                        # actual square dim gray, white border
                        # small sqaure (if !=0) white
                            
                            
                            color = "darkslategrey" 
                            rect = plt.Rectangle([x+1 - largestError / 2, y+1 - largestError / 2], largestError, largestError,
                                 facecolor=color, edgecolor= color)
                            ax.add_patch(rect)

                            rect = plt.Rectangle([x+1 - size / 2, y+1 - size / 2], size, size,
                            facecolor=color, edgecolor= "white",linestyle=":")
                            ax.add_patch(rect)
                            color = "black" 
                            rect = plt.Rectangle([x+1 - smallestError / 2, y+1 - smallestError / 2], smallestError, smallestError,
                                            facecolor="black", edgecolor= "black")
                            ax.add_patch(rect)
                        else
                            color = "black" 
                            rect = plt.Rectangle([x+1 - size / 2, y+1 - size / 2], size, size,
                            facecolor=color, edgecolor= color)
                            ax.add_patch(rect)
                            # plot the error bar
                            color = "red" 
                            moresize = 0
                            if smallestError < 0 
                               moresize = size/2
                            end
                            rect = plt.Rectangle([x+1 + size/3 , y+1- smallestError/2],  thickness,-moresize-largestError/2+smallestError/2,
                                facecolor=color, edgecolor=color)
                            ax.add_patch(rect)
                        end
                     end        
 

                end
            end 
        end #for y
    end # for x
    ax.autoscale_view()
    ax.invert_yaxis()

end









