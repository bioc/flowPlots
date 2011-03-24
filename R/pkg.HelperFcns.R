## ----------------------------------------------------------------------------
## Helper Functions 
##
## Functions in this file help to prepare group data for plotting.  The 
## original data might be of the type which could be stored in a stackedData
## class.  The data might be plotted using the boxplot function in
## pkg.GroupListBoxplot.R.  The data could also be plotted using ternaryplot() from
## the vcd pkg or by using barplot().    
## 
## ----------------------------------------------------------------------------

## -------------------------------------------------------------------
## makeDataList()
##
## This function makes a list of data, where each item of the list 
## contains data for a group.
##
## INPUTS:  
##   theData - dataframe containing data for 1 or more groups
##   groupVariableName - column in theData identifying group
##   columnsToKeep - vector specifying the column numbers in theData
##                to include in the list of data. 
##
## dataList = makeDataList(anssDF,"group", 10:13)
## -------------------------------------------------------------------

makeDataList = function(theData, groupVariableName, columnsToKeep){

   group = eval(parse(text=paste("theData$",groupVariableName,sep="")))   

   groupValuesUnique  = unique(group) 

   # group names in alphabetical order
   groupValues = groupValuesUnique[order(groupValuesUnique)]

   numGroups = length(groupValues)
   dataList = list(NA)
   
   for(i in 1:numGroups){

      groupValue = groupValues[i]

      subsetCmd = paste("theData[theData$",groupVariableName,"==groupValue,]",sep="")

      groupData = eval(parse(text=subsetCmd)) 

      dataList[[i]] = groupData[,columnsToKeep]     

   }

   return(dataList)

}

## ----------------------------------------------------------------------------------
## makeTernaryData()
##
## This function takes a dataframe of polyfunctional degree (pfd) data and 
## prepares a matrix with 3 columns (eg. PFD1,PFD2,PFD3, or PFD1-2, PFD3-4, PFD5-6) to 
## use as input to ternaryplot() in the vcd pkg.  
##   
## INPUTS: 
##   pfdData - dataframe containing pfd data as columns
##   columns1 - column(s) of pfd data to place in the first column of the matrix
##   columns2 - column(s) of pfd data to place in the second column of the matrix
##   columns3 - column(s) of pfd data to place in the third column of the matrix
##   columnNames  - (optional) vector of names for the three columns of the matrix
## 
## ternaryData = makeTernaryData(pfdData, columns1, columns2, columns3)
## ----------------------------------------------------------------------------------

makeTernaryData = function(pfdData, columns1=1, columns2=2, columns3=3, 
   columnNames=c("PFD=1","PFD=2","PFD=3")){

   numberRows = dim(pfdData)[1]

   ternaryData = matrix(NA, nrow=numberRows, ncol=3)
         
   ternaryData[,1] = ifelse(rep(length(columns1)==1, numberRows), 
     pfdData[,columns1], 
     apply(pfdData[,columns1], 1, sum))

   ternaryData[,2] = ifelse(rep(length(columns2)==1, numberRows), 
      pfdData[,columns2], 
      apply(pfdData[,columns2], 1, sum))

   ternaryData[,3] = ifelse(rep(length(columns3)==1, numberRows), 
      pfdData[,columns3], 
      apply(pfdData[,columns3], 1, sum))

   colnames(ternaryData) = columnNames

   return(ternaryData)

}

## --------------------------------------------------------------------------------
## makeBarplotData()
##
## This function takes a dataframe of profile data and prepares a matrix of 
## data for input to barplot().  
##
## INPUTS:
##   profileData - dataframe of profile data, such as cell %'s in different categories of cells
##   profileColumns - the columns of profileData to include in the barplot
##   groupVariableName - the column in the dataframe containing the group info
##    
## RETURNS: a matrix whose rows represent different profile categories and whose
##   columns represent different groups.  Each cell in the matrix contains the 
##   mean value for the group for a given profile category.
##
## barplotData = makeBarplotData(profileData, profileColumns, groupVariableName)
## --------------------------------------------------------------------------------

makeBarplotData = function(profileData, profileColumns, groupVariableName){

   groupValues = unique(eval(parse(text=paste("profileData$", groupVariableName, sep=""))))      
   groupValuesOrdered = groupValues[order(groupValues)]
   numGroups = length(groupValues)

   barplotData = matrix(NA, nrow=length(profileColumns), ncol=numGroups)

   for(i in 1:numGroups){

      groupValue = groupValuesOrdered[i]

      subsetCmd = paste("profileData[profileData$",groupVariableName,"==groupValue,]",sep="")
      groupData = eval(parse(text=subsetCmd)) 

      barplotData[,i] = apply(groupData[,profileColumns], 2, mean)

   }

   colnames(barplotData) = groupValuesOrdered
   rownames(barplotData) = colnames(profileData[,profileColumns])

   return(barplotData)

}

## --------------------------------------------------------------------------------------------
## computePFDGroupStatsList()
## 
## This function can be used with ternaryplot() to add PFD group stats to, say, the legend.
## The stats computed are group size (N), pfd group mean, and pfd group sd.
##
## INPUTS: 
##  groupPFDDataList - one list item per group, each list item contains a matrix of PFD
##       percentages; the rows are subjects, and the columns are pfd categories. 
##  pfdValues - vector of the PFD values that the columns in each matrix in the 
##       groupPFDDataList represent; eg. 1:3 for (PFD1,PFD2,PFD3).
##  numDigitsMean - return a mean rounded to this number of digits
##  numDigitsSD - return a standard deviation rounded to this number of digits.  
##
## computePFDGroupStatsList(groupPFDDataList, pfdValues=1:3, numDigitsMean=3, numDigitsSD=2)
##
## --------------------------------------------------------------------------------------------

weightedPFDMean = function(x, pfdValues){

  meanPFD = sum(pfdValues*(x/100))

  return(meanPFD)
}

computePFDGroupStatsList = function(groupPFDDataList,pfdValues=1:3,numDigitsMean=3,numDigitsSD=2){
 
   pfdGroupStatsList = list(NA)
   numGroups = length(groupPFDDataList)

   for(i in 1:numGroups){

      ## compute meanPFD for each subject; eg. for 30%PFD1,20%PFD2,50%PFD3,
      ## meanPFD = .3*1 + .2*2 + .5*3 = 2.2
      individMeanPFD = apply(groupPFDDataList[[i]], 1, weightedPFDMean, pfdValues)

      groupN = length(individMeanPFD)
      groupMean = format(mean(individMeanPFD), digits=numDigitsMean)
      groupSD = format(sd(individMeanPFD), digits=numDigitsSD)

      groupStats = NA
      groupStats = c(groupN, groupMean, groupSD)

      pfdGroupStatsList[[i]] = groupStats         
   }

   return(pfdGroupStatsList)
}

## -----------------------------------------------------------------------------
## legendPFDStatsGroupNames()
##
## This function returns a vector of groupNames of the form:
## "Adults (25) 1.5/.6", which represents the group name, number of subjects in
## the group, the pfd mean / pfd standard deviation.  
##
## INPUTS:
##   pfdGroupStatsList - a list of vectors containing the pfd group stats of 
##      group size, pfd group mean, and pfd group standard deviation.  
##   groupNames - a vector of group names, such as c("Adult", "Neonate").
## 
## legendPFDStatsGroupNames(pfdGroupStatsList,groupNames)
## -----------------------------------------------------------------------------

legendPFDStatsGroupNames = function(pfdGroupStatsList, groupNames){

   newGroupNames = NA
   numGroups = length(pfdGroupStatsList)

   for(i in 1:numGroups){

      groupN = pfdGroupStatsList[[i]][1]
      groupMean = pfdGroupStatsList[[i]][2]
      groupSD = pfdGroupStatsList[[i]][3]

      newGroupNames[i] = paste(groupNames[i], " (", groupN, "), ",
         groupMean, "/", groupSD, sep="")   
   }

   return(newGroupNames)
}




















