## -------------------------------------------------------------------------------------------------------
##  StackedData refers to data originating from an ICS Flow Cytometry experiment where the cytokine
##   combinations for a given cell type, stimulus, and concentration, say, are stacked; i.e. a subset
##   of such data could look like this:
##
##   id group stim concGroup cell percentAll count totalCount percentReactive            cytCombo
## a2004 adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6+IL12+IFNa+
## a2004 adult  LPS         3  mDC       0.43     3        700        0.940625 TNFa+IL6+IL12+IFNa-
## a2004 adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6+IL12-IFNa+
## a2004 adult  LPS         3  mDC      21.86   153        700       47.818750 TNFa+IL6+IL12-IFNa-
## a2004 adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6-IL12+IFNa+
## a2004 adult  LPS         3  mDC       0.29     2        700        0.634375 TNFa+IL6-IL12+IFNa-
## a2004 adult  LPS         3  mDC       0.00     0        700        0.000000 TNFa+IL6-IL12-IFNa+
## a2004 adult  LPS         3  mDC      19.71   138        700       43.115625 TNFa+IL6-IL12-IFNa-
##   
##   The StackedData class has 6 slots for data: 
##
##   stackedData - data frame, 
##   markers - matrix, 
##   profilePercents - data frame, 
##   marginalData - data frame, 
##   pfdData - data frame, 
##   pfdPartsData - list of data frames. 
##
##   The StackedData class has methods to compute: 
##
##   the marker matrix, 
##   profilePercents, 
##   marginalData, 
##   PFDData, 
##   PFDPartsData 
##
##   all from the stackedData. 
##
##   The class also includes getter and setter methods for each of the data slots.    
##
## -------------------------------------------------------------------------------------------------------

## ----------------------
## class = StackedData
## ----------------------

setClass("StackedData", representation( stackedData="data.frame", 
   profileData="data.frame", marginalData="data.frame", pfdData="data.frame", 
   pfdPartsData="list", markers="matrix" ) )

## -----------------
## get Methods
## -----------------

## -- Method Signatures --------------

setGeneric("stackedData", function(object) standardGeneric("stackedData") )
setGeneric("profileData", function(object) standardGeneric("profileData") )
setGeneric("marginalData", function(object) standardGeneric("marginalData") )
setGeneric("pfdData", function(object) standardGeneric("pfdData") )
setGeneric("pfdPartsData", function(object) standardGeneric("pfdPartsData") )
setGeneric("markers", function(object) standardGeneric("markers") )

## -- Method Defns --------------

setMethod("stackedData", signature(object = "StackedData"), 
   function(object) object@stackedData )

setMethod("profileData", signature(object = "StackedData"), 
   function(object) object@profileData )

setMethod("marginalData", signature(object = "StackedData"), 
   function(object) object@marginalData )

setMethod("pfdData", signature(object = "StackedData"), 
   function(object) object@pfdData )

setMethod("pfdPartsData", signature(object = "StackedData"), 
   function(object) object@pfdPartsData )

setMethod("markers", signature(object = "StackedData"), 
   function(object) object@markers )

## ------------------------------------
## set Methods
## ------------------------------------

## -- Method Signatures --------------

setGeneric("stackedData<-", 
   function(object, value) standardGeneric("stackedData<-") )

setGeneric("profileData<-", 
   function(object, value) standardGeneric("profileData<-") )

setGeneric("marginalData<-", 
   function(object, value) standardGeneric("marginalData<-") )

setGeneric("pfdData<-", 
   function(object, value) standardGeneric("pfdData<-") )

setGeneric("pfdPartsData<-", 
   function(object, value) standardGeneric("pfdPartsData<-") )

setGeneric("markers<-", 
   function(object, value) standardGeneric("markers<-") )

## -- Method Defns --------------

setReplaceMethod("stackedData","StackedData", 
   function(object, value){ 
      object@stackedData <- value
      return(object) 
   }
)

setReplaceMethod("profileData","StackedData",  
   function(object, value){ 
      object@profileData <- value
      return(object) 
   }
)

setReplaceMethod("marginalData","StackedData",  
   function(object, value){ 
      object@marginalData <- value
      return(object) 
   }
)

setReplaceMethod("pfdData","StackedData", 
   function(object, value){
      object@pfdData <- value
      return(object) 
   }
)

setReplaceMethod("pfdPartsData","StackedData", 
   function(object, value){
      object@pfdPartsData <- value 
      return(object)
   }
)

setReplaceMethod("markers","StackedData", 
   function(object, value){
      object@markers <- value
      return(object) 
   }
)

## -------------------------------------
## compute Methods
## -------------------------------------

## -- Method Signatures --------------

setGeneric("readStackedData", 
   function(fileName) standardGeneric("readStackedData") )

setGeneric("computeProfileData", 
   function(object, byVarNames, idVarName, percentVarName, groupVarName)
   standardGeneric("computeProfileData") )

setGeneric("computeMarginalData", 
   function(object, byVarNames, idVarName, percentVarName, groupVarName)
   standardGeneric("computeMarginalData") )

setGeneric("computePFDData", 
   function(object, byVarNames, idVarName, percentVarName, groupVarName) 
   standardGeneric("computePFDData") )

setGeneric("computePFDPartsData", 
   function(object, byVarNames, idVarName, percentVarName, groupVarName)
   standardGeneric("computePFDPartsData") )

setGeneric("computeMarkers", 
   function(markerNames, includeAllNegativeRow) standardGeneric("computeMarkers") )

## -- Method Defns --------------

setMethod("readStackedData", 

   signature(fileName = "character"), 

   function(fileName){
      stackedData = read.csv(fileName, stringsAsFactors=FALSE)
      return(stackedData)
   }
)

setMethod("computeProfileData", 

   signature(object = "StackedData", byVarNames="character", 
      idVarName="character", percentVarName="character", groupVarName="character"), 

   function(object, byVarNames=NA, idVarName=NA, percentVarName=NA, groupVarName=NA){

      fcnName = "profilePercents"

      stackedData = stackedData(object)
      markers = markers(object)

      profileData = createSummaryData(stackedData, markers, fcnName, percentVarName, 
         idVarName, byVarNames, groupVarName)
         
      return(profileData)

   } # fcn: end

) # method: end


setMethod("computeMarginalData", 

   signature(object = "StackedData", byVarNames="character", 
      idVarName="character", percentVarName="character", groupVarName="character"), 

   function(object, byVarNames=NA, idVarName=NA, percentVarName=NA, groupVarName=NA){

      fcnName = "markerPercentPositiveAnyMarker"

      stackedData = stackedData(object)
      markers = markers(object)

      marginalData = createSummaryData(stackedData, markers, fcnName, percentVarName, 
         idVarName, byVarNames, groupVarName)
         
      return(marginalData)

   } # fcn: end

) # method: end


setMethod("computePFDData", 

   signature(object = "StackedData", byVarNames="character", 
      idVarName="character", percentVarName="character", groupVarName = "character"), 

   function(object, byVarNames=NA, idVarName=NA, percentVarName=NA, groupVarName=NA){

      fcnName = "pfd"
      
      stackedData = stackedData(object)
      markers = markers(object)

      pfdData = createSummaryData(stackedData, markers, fcnName, percentVarName, idVarName, 
         byVarNames, groupVarName)
         
      return(pfdData)
     
   } # fcn: end

) # method: end


setMethod("computePFDPartsData", 
   
   signature(object = "StackedData", byVarNames="character", 
      idVarName="character", percentVarName="character", groupVarName="character"), 

   function(object, byVarNames=NA, idVarName=NA, percentVarName=NA, groupVarName=NA){
    
      fcnName = "pfdParts"

      stackedData = stackedData(object)
      markers = markers(object)

      pfdPartsData = createSummaryData(stackedData, markers, fcnName, percentVarName, 
         idVarName, byVarNames, groupVarName)
         
      return(pfdPartsData)
      
   } # fcn: end

) # method: end


setMethod("computeMarkers", 
   signature(markerNames = "character", includeAllNegativeRow = "logical"), 
   function(markerNames, includeAllNegativeRow=TRUE){

      numMarkers = length(markerNames)
      markers = matrix(NA, nrow=2^numMarkers, ncol=numMarkers)

      for(i in 1:numMarkers){
         markers[,i] = rep( rep(c(1,0),each=2^(numMarkers-i)), 2^(i-1) )
      }

      colnames(markers) = markerNames

      if(includeAllNegativeRow){

         return(markers)

      }
      else{
         indexLastRow = dim(markers)[1]
         return(markers[-indexLastRow,])
      }

   }
)

## -----------------------------------------------------------
## Show Method
## -----------------------------------------------------------

setMethod("show", "StackedData",  
  function(object){

     cat("Class: StackedData [package: flowPlots]\n")
     cat("\n")
     cat("Slot Dimensions:\n")
     cat("\n")
 
     slotDim=list( dim(object@stackedData), dim(object@profileData), dim(object@marginalData),  
        dim(object@pfdData), length(object@pfdPartsData), dim(object@markers) )

     names(slotDim) = c("stackedData","profileData","marginalData",
        "pfdData","pfdPartsData","Markers")

     print(slotDim)

}


)









