## -----------------------------------------------------------------------------
## PKG: flowStats
## 
## FILE: pkg.HelperPlotFcns.R
##
## The functions in this file are called by the boxplot function in 
## pkg.GroupListBoxplot.R 
##
## -----------------------------------------------------------------------------

## --------------------------------------------------------------------------
## getPValueLabels
##
## Returns a text vector of labeled p-values, eg. (p < .01, p = .34, p = NA)
## --------------------------------------------------------------------------

getPValueLabels = function(pValues, testsRoundDigits=2, pvalueLabel="p"){
 
   numPValues = length(pValues)
   pvalueLabels = NA

   pvalMin = 1/(10^testsRoundDigits)
   pvalMax = 1-pvalMin

   for(i in 1:numPValues){

      if(is.na(pValues[i])){
            pvalueLabels[i] = paste(pvalueLabel, " = ", "NA", sep="")
      }
      else{
         if(pValues[i] < pvalMin){
            pvalueLabels[i] = paste(pvalueLabel, " < ", pvalMin, sep="")
         }
         else{
            if(pValues[i] > pvalMax){
               pvalueLabels[i] = paste(pvalueLabel, " > ", pvalMax, sep="")
            }
            else{
               pvalueLabels[i] = paste(pvalueLabel, " = ", pValues[i], sep="")
   }}}}

   return(pvalueLabels)

}

## ------------------------------------------------
## getMaxOfListData
##
## Returns the minimum value of the data
## found in a list of dataframes, excluding NA's
## ------------------------------------------------

getMaxOfListData = function(dataList){

   listLength = length(dataList)
   allData = NA

   for(i in 1:listLength){

      allData = rbind(allData,dataList[[i]])

   }   

   maxDataValue = max(allData, na.rm=T)

   return( maxDataValue )
}

## ------------------------------------------------
## geMinOfListData
##
## Returns the minimum value of the data
## found in a list of dataframes, excluding NA's
## ------------------------------------------------

getMinOfListData = function(dataList){

   listLength = length(dataList)
   allData = NA

   for(i in 1:listLength){

      allData = rbind(allData,dataList[[i]])

   }   

   minDataValue = min(allData, na.rm=T)

   return( minDataValue )
}

## ----------------------------------------------------------
## getNsOfListData(dataList)
##
## Returns: 
## 
## NsVector, a text vector: 
## one group example for 3 cats on x-axis: (3,4,5)
## two group example for 3 cats on x-axis: (3/3,4/3,5/1)
## ----------------------------------------------------------

getNsOfListData = function(dataList){

   numGroups = length(dataList)
   theNsAll = NA
   NsVector = NA

   for(i in 1:numGroups){

      dataGroup = dataList[[i]]

      ## theNs is a vector of the number of nonmissing values for each
      ## column of the dataGroup dataframe; i.e. for each category 
      ## to be plotted on the x-axis.  
      theNs = getNs(dataGroup)

      if(numGroups==1){NsVector=theNs}
      else{ ## Each row of theNsAll matrix stores theNs for a group
         theNsAll = rbind(theNsAll, theNs)            
      }
   }

   if(numGroups >= 2){ ## Insert / in between numbers

      ## Remove first row of NA's
      theNsAll = as.data.frame(theNsAll[-1,])  
      numCols = dim(dataGroup)[2]

      if(numGroups > 2){
         for(i in 1:numCols){ ## Add extra space around the /           
            NsVector[i] = paste(theNsAll[,i], collapse=" / ")
      }}
      else{
         if(numGroups == 2) { ## Use less space around the /
            for(i in 1:numCols){           
              NsVector[i] = paste(theNsAll[,i], collapse="/")
      }}}

   }
   return(NsVector)
}

## -------------------------------------------------
## getNs(dataframe) 
##
## -- called by getNsListOfData(dataList)
##
## Returns a vector of the number of nonmissing
## values in each column of the dataframe.
## --------------------------------------------------

getNs = function(theData){

   numCols = dim(theData)[2]
   theNs = NA

   for(i in 1:numCols){

         dataCol = theData[,i]
         theNs[i] = length(dataCol[!is.na(dataCol)])
      }

   return(theNs)
}

## ------------------------------------------------------------------------
## get2GroupPValues(dataframe1, dataframe2, pairs=FALSE, test="Wilcoxon")
##
## Returns a vector of p-values comparing two groups on one or more 
## categories.  
## ------------------------------------------------------------------------

get2GroupPValues = function(dataGroup1, dataGroup2, pairs=FALSE, test="Wilcoxon"){

   numCols = dim(dataGroup1)[2]
   pValues = NA

   for(i in 1:numCols){

      col1 = dataGroup1[,i]
      col2 = dataGroup2[,i]

      if(length(col1[!is.na(col1)]) > 0 && length(col2[!is.na(col2)]) > 0) {

         if(test=="Wilcoxon"){
            pValues[i] = wilcox.test(dataGroup1[,i], dataGroup2[,i], paired=pairs)$p.value }

      } 
      else{ pValues[i] = NA }         

    }      

   return(pValues)

}

## ---------------------------------------------------------------------------------
## getMultiGroupPValues
##
## Returns vector of pvalues for each column in each group's dataframe; i.e. for
## each category being plotted on the x-axis.
## ---------------------------------------------------------------------------------

getMultiGroupPValues = function(dataList, testType="Kruskal"){

   numCols = dim(dataList[[1]])[2]
   numGroups = length(dataList)
   pValues = NA

   if(testType=="Kruskal"){

      for(i in 1:numCols){

         catData = dataList[[1]][,i]

         for(j in 2:numGroups){
            catData = cbind(catData, dataList[[j]][,i])
         }

         catDataList = apply(catData,2,list)
         pValues[i] = kruskal.test(catDataList)$p.value                    
   }}

   return(pValues)

}

## -------------------------------------------------------------------
## getXAxisAtPoints
##
## Returns the spacing position for a set of groups to be plotted at 
## a given x-point.  For example, for plotting 2 groups at x=1, we call
## this function to get (-.2,.2) which we add to "x", so that group1
## will be plotted at 0.8 and group2 will be plotted at 1.2.  If we
## plan to compare groups at x=1,2,and 3, these spacings can be added
## to each x-value to get the x-positions for each group. 
## -------------------------------------------------------------------

getXAxisAtPoints = function(numGroups){

   # Experimentation may result in selecting different values below

   oneGroup = 0
   twoGroups = c(-.2,.2)
   threeGroups = c(-.3,0,.3)
   fourGroups = c(-.3,-.1,.1,.3)
   fiveGroups = c(-.3,-.15,0,.15,.3)
   sixGroups = c(-.25,-.15,-.05,.05,.15,.25)
   sevenGroups = c(-.225,-.15,-.075,0,.075,.15,.225)
   eightGroups = c(-.35,-.25,-.15,-.05,.05,.15,.25,.35)
   nineGroups = c(-.3,-.225,-.15,-.075,0,.075,.15,.225,.3)
   tenGroups = c(-.4,-.35,-.25,-.15,-.05,.05,.15,.25,.35,.4)
   elevenGroups = c(-.375,-.3,-.225,-.15,-.075,0,.075,.15,.225,.3,.375)

   groupPoints = list(oneGroup,twoGroups,threeGroups,fourGroups, fiveGroups,
      sixGroups, sevenGroups, eightGroups, nineGroups, tenGroups, elevenGroups)

   if(numGroups > 11 || numGroups < 1){

      cat("getXAxisAtPoints not defined for numGroups < 1 or > 11", "\n")
      cat("Returning missing value, NA", "\n")
      return(NA)
   }
   else{

      return( groupPoints[[numGroups]] )

   }

}

## -----------------------------------------------------
## computeYmaxBoxplot
##
## Returns list of: (ymax, yP, yN)
##
## ymaxData is the max of the data, excluding NA's.
## ymax is larger than ymaxData, depending on whether
## pvals and n's are being printed, or just one of the
## two, or neither.
##
## yP and yN are the y positions for the p-values and
## n's, respectively. 
## -----------------------------------------------------

computeYmaxBoxplot = function(dataList, ymaxData, testsTrue, printNs, addToYMax1, 
   addToYMax2){

      yP = yN = NA

      if(testsTrue && printNs){ 

         addToYMax = addToYMax2 

         #-- Print p-values above n's ------
         yP = ymaxData*(1+addToYMax2-.05) 
         yN = ymaxData*(1+addToYMax2-.10)
      }
      else{ if(testsTrue || printNs){ 

         addToYMax = addToYMax1 

         #-- Print only p-vals or n's, so subtract off same.
         yP = ymaxData*(1+addToYMax1-.05) 
         yN = ymaxData*(1+addToYMax1-.05)
          
      }
       else { addToYMax = 0 }      
      } 

      ymax = ymaxData*(1+addToYMax) 

      return(list(ymax, yP, yN))
}

 