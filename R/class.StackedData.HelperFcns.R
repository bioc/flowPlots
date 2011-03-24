## --------------------------------------------------------------------------------------
## These functions are used by the methods in the StackedData class to compute the
## various types of data: profile percent, marginal, pfd, and pfd parts.
## --------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------
## numericFactor
## 
## This function converts the values of x to 'factor' and then to 'numeric'.  
##
## Inputs:
##     x - columns of a data frame consisting of the "byVar" columns of a stacked 
##         data set; for example, the columns: cell(=CD4), visit(=2), stim(=gag).          
##
## Returns - numeric values for x; examples:
# 
## x = c("a","b") becomes xn = numericFactor(x) = c(1,2), where xn is numeric.
## y = c(2,3,4) becomes yn = numericFactor(x) = c(1,2,3), where yn is numeric.
##  
## -----------------------------------------------------------------------------------

numericFactor = function(x){
   return(as.numeric(factor(x)))
}

## -----------------------------------------------------------------------------------
## assignCategoryGroups
##
## This function assigns a group number to each unique combination of byVars.  
##
## Inputs: 
##     sdSS - data frame of the "byVar" columns of a stacked data set; for 
##        example, the columns: cell(=CD4), visit(=2), concGroup(= 3), stim(=gag).  
##
## Returns: a list of 2 elements; the first element is a vector of the group number 
##   assignments for each unique byVar combination in sdSS.  The second element is
##   data frame of the byVars and a groups column indicating the group number 
##   assigned to each unique combination of byVars.   
## -----------------------------------------------------------------------------------

assignCategoryGroups = function(sdSS){

   ## Assign numeric values representing the levels for each column of the sdSS 
   ## data frame which contains the "byVar" data
   sdSSNumFactor = apply(sdSS, 2, numericFactor)
   
   numCols = dim(sdSS)[2]
   if(numCols > 1) {
      for(i in 2:numCols){
         ## Assign a unique number to each byVar combination
         ## Allowable Range for factor values for byVar 1 = [1-99], value 1-99
         ## Allowable Range for factor values for byVar 2 = [1-90], value 100-9000
         ## Allowable Range for factor values for byVar 3 = [1-90], value 10000-90000
         ## Allowable Range for factor values for byVar 4 = [1-90], value 100000-900000
         ## and so on.
         sdSSNumFactor[,i] = sdSSNumFactor[,i]*( 10^(2*(i-1)) )
      }
   }

   ## Assigns a unique group number to each of the possible combinations of "byVar" values   
   groups = as.numeric(factor(apply(sdSSNumFactor, 1, sum)))
   sdSS = cbind(sdSS, groups)

   uniqueGroups = unique(groups)
   categoryGroups = NA

   for(group in uniqueGroups){
      categoryValues = subset(sdSS, groups==group)              
      categoryGroups = rbind(categoryGroups, categoryValues[1,])             
   }

   # Remove meaningless row of NA's
   categoryGroups = categoryGroups[-1,]
   groupAssignments = list(groups, categoryGroups)
  
   return( groupAssignments )
}

## -------------------------------------------------------------------------------------
## createSummaryData
##
## This function is used to create all the types of summary data: profile data, 
## marginal data, pfd data, and pfdParts data.  The user specifies which type of 
## data to compute via the 'summaryDataFcnName' input parameter.
##
## Inputs: 
##    stackedData - data frame of "stacked" data
##    markers - matrix of 0/1's indicating +/- of marker presence in a combination,
##       one row per combination, cols=markers            
##    summaryDataFcnName - text string specifying the type of summary data to compute, 
##       eg. "markerPercentPositiveAnyMarker"
##    percentVarName - the variable of data to summarize  
##    idVarName - text string of id variable name
##    byVarNames - text vector of variables to use to subset the data; for example, 
##       want to summarize by stimulus, concentration Group, and cell type.
##    groupName - text string of the variable name specifying the group membership
##
## Returns: a data frame for profile data, marginal data, or pfd data.  Returns a list
##   for pfdParts data, where each item is the list is a data frame containing the 
##   the pfd components for a given pfd.  
## 
## -------------------------------------------------------------------------------------

createSummaryData = function(stackedData, markers, summaryDataFcnName, 
   percentVarName, idVarName, byVarNames, groupName)
{
   charData = charDataAll = NA
   resultsData = resultsDataAll = NA
   pfdNumber = NA

   numMarkerRows = dim(markers)[1] # Number of categories; i.e. 2^NumberOfMarkers

   ## Subset the data by variables like: cell type, stimulus, concentration group
   sdSS = subset(stackedData, select=byVarNames)

   ## assignCategoryGroups returns a list with two elements:
   ## assignCategoryGroups[[1]]   
   ## The first list element is an integer column representing groups of data specified by 
   ## the byVarNames variables; i.e. an example of a 'category group' is:
   ## visit=2, cellType = CD4, stim = gag, conc = low
   ## Add this column to the stacked data data frame
   ## assignCategoryGroups[[2]]
   ## The second list element is a data frame containing the demographic data with one
   ## row per "category group"
   groupAssignments = assignCategoryGroups(sdSS)
   categoryGroups = groupAssignments[[2]]
   catGroups = groupAssignments[[1]]   
   uniqueGroups = categoryGroups$groups      
   stackedData = cbind(stackedData, catGroups)
   numGroups = length(uniqueGroups)

   ## ------------------------------------------------------
   ## Compute the summary data by the "category groups"  
   ## ------------------------------------------------------
   for(i in 1:numGroups){      

      ssAll = subset(stackedData, catGroups == i) 

      subjectIds = eval(parse(text=paste("stackedData$",idVarName,sep="")))
      uniqueIds = unique(subjectIds)

      ## Vector of demographic data for this character group
      charData = categoryGroups[categoryGroups$groups==i,]

      ## -------------------------------------------------
      ## Compute the summary data by subject
      ## -------------------------------------------------      
      for( subjectId in uniqueIds ){   

         ## subset out the data for given subject
         sscmd = paste("subset(ssAll,", idVarName, "==subjectId)", sep="")
         ssID = eval(parse(text=sscmd))

         ## check if subject data has exactly the same number of rows as
         ## the marker matrix; if so, call the specified function to compute
         ## the summary data.
         if(dim(ssID)[1] == numMarkerRows){  # ok to proceed
            ## add sort var here? 
            resultsDataCmd = paste(summaryDataFcnName,"(ssID$", percentVarName,
               ",markers)", sep="")                        
            resultsData = eval(parse(text=resultsDataCmd))

            ## -----------------------------------------------------------------
            ## Build the summary results data 
            ##
            ## Check if returned data is a data frame or a list: pfdParts 
            ## returns a list, all the others return a data frame.
            ## -----------------------------------------------------------------

            if(is.data.frame(resultsData)){
               ## Build the complete set of results data by adding this subject's summary 
               ## data for the given "category group" to the results data structure
               resultsDataAll = rbind(resultsDataAll, resultsData)
            }
            else{ ## Special handling for results from computation of pfdParts data
               if(is.list(resultsData)){ ## i.e. a result from a call to pfdParts()   
                  resultsDataDF = resultsData[[1]]
                  resultsDataAll = rbind(resultsDataAll, resultsDataDF)
                  pfdNumber = resultsData[[2]]                                 
               }
            }

            groupCmd = paste("ssID$", groupName, "[1]", sep="")
            group = eval(parse(text=groupCmd))

            ## -----------------------------------------------------------------------
            ## Build the demographic data to accompany the summary data; i.e. things
            ## like: subject id, cell type, stimulus, concentration group, group
            ## -----------------------------------------------------------------------

            charDataSubject = cbind(charData, subjectId, group)
            charDataAll = rbind(charDataAll, charDataSubject)
         }         
         else{ ## Skip a subject if the subject doesn't have the expected rows of data
               ## Print the subject id to alert the user of the function.
            warning("\nSkipping id = ", subjectId, ".  Number of rows of subject data 
            does not match the number of rows of marker data.\n")
         }

       } # subjects: process
   }     # category groups: process                         

   ## Remove row of meaningless NA data
   resultsDataAll = resultsDataAll[-1,]
   charDataAll = charDataAll[-1,]

   ## ---------------------------------------------------------------------
   ## Combine the summary results data and the demographic (char) data
   ## 
   ## For pfdParts, return a list of data frames; for all others, return a 
   ##    single data frame
   ## ---------------------------------------------------------------------
   if(summaryDataFcnName != "pfdParts")
   {   
       resultsDataAll = cbind(resultsDataAll, charDataAll, stringsAsFactors=FALSE)
       rownames(resultsDataAll) = NULL

   }else{ ## pfdParts data: return a list of data frames, one data frame for each possible 
          ## PFD, except for the max PFD since there is only one possible combination for 
          ## the max PFD.  Here, build the list to return. 

      maxPFDNumber = max(pfdNumber)
      resultsDataList = list( rep(NA, maxPFDNumber) )
      colStart = 1

      for(pfdNum in 1:maxPFDNumber){

         pfdPartDF = NA
         numThisPFD = length(pfdNumber[pfdNumber==pfdNum])
         colEnd = colStart + (numThisPFD-1)
         pfdPartDF = resultsDataAll[,colStart:colEnd]
         ## Combine the summary results and demographic (char) data
         pfdPartDF = cbind(pfdPartDF, charDataAll)         
         rownames(pfdPartDF) = NULL

         resultsDataList[[pfdNum]] = pfdPartDF

         colStart = colEnd + 1
      }

      ## Here remove the data for the maxPFD
      ## Do it this way; i.e. add it above, but remove it here, in case later it's 
      ## decided that it's better to include it in the result
      indexMaxPFD = length(resultsDataList)
      resultsDataAll = resultsDataList[-indexMaxPFD]  
   }

   return(resultsDataAll)

}

## ----------------------------------------------------------------------------------------
## Create the profile percents data: profilePercents() 
##
## This function creates on a per subject level for a 'category group', such as:
## visit=2, cellType = CD4, stim = gag, conc = low, the profile percents data in a 
## horizontal fashion (rather than the "stacked" (vertical) fashion) for each of the 
## marker categories.  This function just creates a data frame with the input, percent,
## as a row, so that the percents are horizontal in a data frame.
##
## Inputs: 
##     percent - the cell percentages for a subject for a "category group", such as:
##         visit=2, cellType = CD4, stim = gag, conc = low.  The length of percents is
##         the same as the number of rows of markers.  The percent vector should be in
##         the same order as the markers matrix.  
##     markers -- matrix of 0 and 1's; rows = +/- categories
##        such as, TNF+IFNg+IL2+, cols = markers, such as
##        TNFa, IFNg, IL2.   
##
## Returns: a data frame of profile percent data.
##  
## profile percent data can look like:
##
## TNFa+IFNg+IL2+ TNFa+IFNg+IL2- TNFa+IFNg-IL2+ TNFa-IFNg+IL2+ TNFa+IFNg-IL2- TNFa-IFNg+IL2-
## 1.5            5.5            6.0             1.0           12             6
##
## TNFa-IFNg-IL2+ TNFa-IFNg-IL2-
## 1.0            67.0
## ----------------------------------------------------------------------------------------

profilePercents = function(percent, markers){

   ## Create a matrix identical to the marker matrix except that the 1's and 0's are
   ## replaced by +'s and -'s.  Then, use it below to create category labels.
   plusMinus = markers
   plusMinus[markers==1] = "+"
   plusMinus[markers==0] = "-"

   markerNames = colnames(markers)
   numRowsMarkers = dim(markers)[1]
   percentLabels = NA

   ## Create category labels for the profile percents; eg. TNFa+IFNg+IL2+   
   for(i in 1:numRowsMarkers){
      x=plusMinus[i,]
      percentLabels[i] = paste(paste(markerNames,x,sep=""),collapse="")
   }

   ## Set the percent data as a row in a data frame
   ## Use the category labels generated above as the column names
   profilePs = as.data.frame(t(as.data.frame(percent)))
   colnames(profilePs) = percentLabels

   return( profilePs )
}

## ----------------------------------------------------------------------------------------
## Create the Marginal Data: markerPercentPositive()
## 
## This function computes on a per subject level for a 'category group', such as:
## visit=2, cellType = CD4, stim = gag, conc = low, the marginal data for each of the
## markers. 
##
## Inputs: 
##     percent - the cell percentages for a subject for a "category group", such as:
##         visit=2, cellType = CD4, stim = gag, conc = low.  The length of percents is
##         the same as the number of rows of markers.  The percent vector should be in
##         the same order as the markers matrix.  
##     markers -- matrix of 0 and 1's; rows = +/- categories
##        such as, TNF+IFNg+IL2+, cols = markers, such as
##        TNFa, IFNg, IL2.   
##
## Returns: a data frame of marginal data
##
## Marginal data can look like this:
## TNFa IFNg IL2
## 6.3  3.5  1.1
##
## ----------------------------------------------------------------------------------------

markerPercentPositive = function(percent, markers){

   numMarkers = dim(markers)[2] # For markers TNFa,IFNg,IL2 equals 3
   numCats = dim(markers)[1]    # For 3 markers, equals 2^3 = 8

   markerPP = NA
   sumPercents = 0

   ## One marker at a time, compute the combined percent 
   for(i in 1:numMarkers){
      ## For the given marker, sum the percents where the marker is +;
      ## i.e. where the marker is 1 in the marker matrix
      for(j in 1:numCats){
         if(markers[j,i]==1){
            sumPercents = sumPercents + percent[j]
         }         
      }
      ## Set the marginal percent for current marker
      markerPP[i] = sumPercents
      sumPercents = 0     
   }

   markerPP = as.data.frame(t(as.data.frame(markerPP)))
   colnames(markerPP) = colnames(markers)

   return( markerPP )
}

## --------------------------------------------------------------------------------------
## Create the Marginal and AnyMarker Data: markerPercentPositiveAnyMarker
##
## This function computes on a per subject level for a 'category group', such as:
## visit=2, cellType = CD4, stim = gag, conc = low, the marginal data for each of the
## markers and the percent of all markers combined; i.e. the anyMarker percent.  This 
## function calls markerPercentPositive to compute the marginal data for each of the 
## markers and then sums that data to compute the anyMarker percent.
##
## Inputs: 
##     percent - the cell percentages for a subject for a "category group", such as:
##         visit=2, cellType = CD4, stim = gag, conc = low.  The length of percents is
##         the same as the number of rows of markers.  The percent vector should be in
##         the same order as the markers matrix.  
##     markers -- matrix of 0 and 1's; rows = +/- categories
##        such as, TNF+IFNg+IL2+, cols = markers, such as
##        TNFa, IFNg, IL2.   
##
## Returns: a data frame of marginal data
##
## Marginal data can look like this:
## TNFa IFNg IL2 AnyMarker
## 6.3  3.5  1.1 10.9
##
## --------------------------------------------------------------------------------------

markerPercentPositiveAnyMarker = function(percent, markers){

   ## Compute the individual marker percents 
   markerPP = markerPercentPositive(percent, markers)

   ## Compute the percent for all markers combined
   markerCount = apply(markers, 1, sum) 
   anyMarker = sum(percent[markerCount>0])

   ## Combine the individual marker percents with the all markers combined percent 
   markerPPAM = cbind(markerPP, anyMarker)

   return( markerPPAM )
}

## -------------------------------------------------------------------------------------
## Create the Polyfunctional Degree (PFD) Data
## 
## This function computes on a per subject level for a 'category group', such as:
## visit=2, cellType = CD4, stim = gag, conc = low., the polyfunctional degree data 
## (PFD).  The percents in the PFD data are based on reactive cells only.
##
## Inputs: 
##     percent - the cell percentages for a subject for a "category group", such as:
##         visit=2, cellType = CD4, stim = gag, conc = low.  The length of percents is
##         the same as the number of rows of markers.  The percent vector should be in
##         the same order as the markers matrix.  
##     markers -- matrix of 0 and 1's; rows = +/- categories
##        such as, TNF+IFNg+IL2+, cols = markers, such as
##        TNFa, IFNg, IL2.
##
## Returns: a data frame of pfd data.
##
## pfdData can look like this:
##
## PFD1 PFD2 PFD3 
## 70   20   10
##
## --------------------------------------------------------------------------------------

pfd = function(percent, markers){

   ## Compute the pfd number for each category; i.e. each row of the marker matrix
   pfdNumber = apply(markers, 1, sum, na.rm=T)

   ## Markers can be cytokines, such as: TNFa, IFNg, IL2
   numMarkers = dim(markers)[2]

   ## Sum percentages for each possible PFD in the range: [min=1, max=numMarkers]
   pfdData = NA   
   for(i in 1:numMarkers){
      pfdData[i] = sum(percent[pfdNumber==i])
   }

   ## The PFD percentages are based on reactive cells only, so revise via denominator 
   sumReactivePercent = sum(pfdData)
   pfdData = 100*pfdData/sumReactivePercent

   ## Prepare the final data frame to return   
   colLabels = paste("PFD", 1:numMarkers, sep="")   
   pfdData = as.data.frame(t(as.data.frame(pfdData)))
   colnames(pfdData) = colLabels   

   return( pfdData )
}

## ---------------------------------------------------------------------------------------
## Create the Composition of Polyfunctional Degree (PFD) Data; i.e. pfdParts data
## 
## This function computes on a per subject level the pfd parts data.  If there are
## 3 markers, pfd can range from 1 to 3.  When pfd=1, there are 3 possible ways for
## that to occur: marker1+marker2-marker3-, marker1-marker2+marker3-, 
## or marker1-marker2-marker3+.  The pfd parts data for pfd1 consists of the 
## percentages of pfd=1 cells in each of those categories, so that the three %'s
## add up to 100.  Similarly, for pfd=2.  For pfd=3, there is only one way to 
## achive this, that is, marker1+marker2+marker3+.   
## 
## Inputs:
##     percent - the cell percentages for a subject for a "category group", such as:
##         visit=2, cellType = CD4, stim = gag, conc = low.  The length of percents is
##         the same as the number of rows of markers.  The percent vector should be in
##         the same order as the markers matrix.  
##     markers -- matrix of 0 and 1's; rows = +/- categories
##        such as, TNF+IFNg+IL2+, cols = markers, such as
##        TNFa, IFNg, IL2.
##
## Returns: A list of a data frame of pfdPartsData and a vector indicating the pfd
##     number for each column of pfdPartsData.  
##
##  pfdPartsData can look like this:
##
##  TNFa IFNg IL2 TNFa:IFNg TNFa:IL2 IFNg:IL2 TNFa:IFNg:IL2
##  90   8    2   70        28       2        100
##
##  pfdNumber can look like this (providing the pfd for the above 7 columns):
##  c(1,1,1,2,2,2,3)  
## ---------------------------------------------------------------------------------------

pfdParts = function(percent, markers){

   ## Compute the pfd number for each category; i.e. each row of the marker matrix
   pfdNumber = apply(markers,1,sum,na.rm=T)

   ## For each PFD, sum the percents of all the possible combinations
   ## Use this as a denominator for each PFD, when computing the compositional 
   ## percentage for each pfd 'part'
   denom = NA
   for(i in 1:max(pfdNumber)){
      denom[i] = sum(percent[pfdNumber==i])
   }

   ## Compute the compositional percent for each pfd 'part'
   percentCat = NA
   indices = order(pfdNumber[pfdNumber!=0])
   ## Compute these in order; i.e. pfd=1's, pfd=2's, etc.
   for(index in indices){
      percentCat[index] = percent[index]/denom[pfdNumber[index]]
   }

   ## Put Percents on the 0-100 scale, rather than 0-1 scale.
   pfdPartsData = NA
   pfdPartsData = 100*percentCat[indices]
   ## Set cases of 0%/0% = NaN to NA
   pfdPartsData[is.nan(pfdPartsData)] = NA

   ## Create pfdParts names like: TNFa, IFNg, IL2, TNFa:IFNg, TNFa:IL2, etc.
   markerNames = colnames(markers)
   pfdPartsDataLabels = NA
   for(i in 1:length(indices)){ 
      pfdPartsDataLabels[i] = paste(markerNames[markers[indices[i],]!=0],collapse=":")
   }

   ## Vector indicating the pfd number for each column of the pfdPartsData.  
   pfdNumber = pfdNumber[indices]

   pfdPartsData = as.data.frame(t(as.data.frame(pfdPartsData)))
   colnames(pfdPartsData) = pfdPartsDataLabels   

   return( list(pfdPartsData, pfdNumber) )
}
