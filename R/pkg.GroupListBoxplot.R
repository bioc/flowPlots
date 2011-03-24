## ---------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------
##
## FUNCTION: GroupListBoxplot()
##
## PURPOSE:
## 
## Create a panel of boxplots, with one or more boxplots at each "x" point.
## The plots can be annotated with "n's" for each group, and p-values for tests
## comparing the groups at each "x" point.  
##
## DATA INPUT: 
##
## 1)dataList -- a list of dataframes, one dataframe for each group.  The columns of
## a group's dataframe represent each marker being examined; i.e. each "x" point. 
## 
## Eg., dataList = list(dataframeGroup1, dataframeGroup2, dataframeGroup3), and 
## dataframeGroup1 is a dataframe with: 
## 1)columns = percentPositiveTNFa, percentPositiveIL2, percentPositiveAnyMarker
## 2)rows = one row per subject for a given condition, eg. CD4, gag, visit 2.
##
## DEPENDS ON OTHER FUNCTIONS: 
##
## Functions defined in pkg.HelperPlotFcns.R
## 
## getMaxOfListData(), getMinOfListData(), getNsOfListData(), getNs(),
## getXAxisAtPoints(),    
## get2GroupPValues(), getMultiGroupPValues(), getPValueLabels()
##
## AUTHOR: Natalie Hawkins
## DATE: Sept 16, 2009
##
## NOTE: Theoretically handles up to 11 groups plotted at each "x" point.  At some point, 
## crowding problem possible.  The limitation in the code is found in setting the xGroup 
## variable which holds the "at" positions on the x-axis.  Also, groupColorVector has 11 
## values.  
##
## ---------------------------------------------------------------------------------------
## ---------------------------------------------------------------------------------------

GroupListBoxplot = function( 

   ## The INPUTS ----------:

   #---- The List of Data to Plot: 
   # eg., dataList[[1]] = group1DataFrame, dataList[[2]] = group2DataFrame, etc.
   dataList, 

   #---- The Boxplot Plot Region -----------------
   # ymax for the boxplot will be set by R's boxplot() fcn, unless set here. 
   # addToYmax2, this value is added to ymax if printing p-vals and N's, 
   # addToYmax1, this value is added to ymax if only printing p-vals or N's
   ymaxBoxplot=NA, addToYmax2=.2, addToYmax1=.1, 

   #---- The Boxes in the boxplots ---------------
   # Width of the box, Color of the box, print or suppress outlier points on the plot,
   # Vector of unique colors for a set of up to 11 groups.
   # Type and width of line around box, Type and width of line at the median,
   # Point char and size plotted at median, if desired. 
   boxWidth=.10, boxColor=8, boxOutliers=TRUE, 
   groupColorVector=c(2,4,5,6,7,9,3,10,11,8,1),
   boxlty=1, boxlwd=1, medlty=1, medlwd=3, medpch=NA, medcex=NA,

   #---- Legend ------------------
   # Print legend on plot, names to use in legend, legend X location, legend Y location,
   # point size for legend text, symbol or character to print next to legend names, colors 
   # of the points or lines next to the legend names, line type to print next to legend names,
   # width of line to print next to legend names, legend title, print points next to legend
   # names, print lines next to legend names
   legendInclude=TRUE, legendGroupNames=paste("Group ", 1:length(dataList), sep=""), 
   legendX=NA, legendY=NA, legendCEX=1, legendPCH=1, legendColors=1:length(legendGroupNames), 
   legendLTY=NA, legendLWD=NA, legendTitle=NA, legendPoints=TRUE, legendLines=FALSE,

   #---- The Points --------------
   # overlay the points on the boxplot, symbol or character to plot, size of point plotted, 
   # color of point plotted -- can be single value, vector, matrix, or list, amount to scatter
   # points around group x position  
   printPoints=TRUE, pointChar=1, pointCEX=1, pointColor=1:length(legendGroupNames), 
   pointJitter=.25,     

   #---- Main Title --------------
   # text for main title, size of main title, regular or bold font, margin line on
   # which to print main title - minimum is 0. 
   mainTitle="Boxplots", mainTitleCEX=1, mainTitleFont=1, mainTitleLine=1,

   #---- Stat Tests Title (can print under main title) --------
   # size of stat test title, regular or bold font, margin line on which to print title - minimum
   # is 0  
   testTitleCEX=1, testTitleFont=1, testTitleLine=0, 

   #---- X and Y Labels ----------
   # label for X axis, label for Y axis, size of x and y axes labels, regular or bold font
   # for X and Y axes labels
   xlabel="X Axis", ylabel="Y Axis", xylabelsCEX=1, xylabelsFont=1,

   #---- X Axis, X Margin Text (use for (0,0,0,1) cytokine combos, Y Axis --------
   # labels for ticks on x axis, size of labels on x ticks, regular or bold font for x tick
   # labels, horizontal or vertical x tick labels, text to place in x-axis margin, size of
   # text to place in x-axis margin, regular or bold font for text in x-axis margin, x position
   # for text in x-axis margin, size of y tick labels, regular or bold font for y axis tick
   # labels, horizontal or vertical y tick labels  
   xAxisLabels=NA, xAxisCEX=1, xAxisFont=1, xAxisRotation=1,
   xMtext="", xMtextCEX=1, xMtextFont=1, xAtMtext=0, 
   yAxisCEX=1, yAxisFont=1, yAxisRotation=1,

   #---- Box Around Entire Plot (draw separately from boxplot call) ---------
   # line width of box drawn around entire plot
   plotBoxLWD=1, 

   #---- Comparison Tests General ----
   # number of digits to report in p-value, size of p-value text, regular or bold font
   # for p-value text, y position for p-value text, label to use to precede numerical 
   # p-value
   testsRoundDigits=2, pCEX=1, pFont=1, yP=NA, pvalueLabel="p",

   #---- Group Comparison Tests --
   # include group comparison p-values in plots, are the groups paired?, group comparison
   # test to perform
   betweenGroupTestsCompute=TRUE, pairedGroups=FALSE, 
   # testType="Wilcoxon", # Maybe use this later, if add more test types.

   # --- N's -----------------------
   # include sample sizes on plots, size of the sample size text to print, 
   # regular or bold text for n's, y position for the sample size text
   printNs=TRUE, nCEX=1, nFont=1, yN=NA

   ## --- END of INPUTS --------------------------

   ){

   ## Start of Function ------------------------------------------------------

   ## --------------------------------------------------------------------------------
   ## Overview of Function Steps
   ##
   ## 00. Get data info: numGroups, grab group1 data, numSubjectsGroup1, numCats
   ## -2. Set box border colors
   ## -1. Set (ymin, ymax) for boxplot region, and yP and yN for p-vals and n's.
   ##  0. Get the group x-spacings for plotting the data points
   ##  1. Plot the boxplots
   ##  2. Plot the data points
   ##  3. Label the x and y Axes
   ##  4. Put box around plot
   ##  5. Add Legend
   ##  6. Add Main Title to Plot
   ##  7. Compute Between Group Tests
   ##  8. Add N's to Plot
   ## -------------------------------------------------------------------------------

   ## ----------------------------------------------------------------------------
   ## 00. Get data info: numGroups, grab group1 data, numSubjectsGroup1, numCats
   ## ----------------------------------------------------------------------------

   numGroups = length(dataList) 
   dataGroup1 = dataList[[1]]   
   numSubjectsGroup1 = dim(dataGroup1)[1]
   ## cats = categories: eg categories, TNFa,IFNg,IL2 
   numCats = dim(dataGroup1)[2]

   ## ---------------------------------------
   ## -2. Set box border colors
   ## ---------------------------------------

   if( length(boxColor)==numGroups ){ boxColorVector = boxColor }
   else{ boxColorVector = rep(boxColor[1], numGroups) }    
    
   ## ------------------------------------------------------------------------------
   ## -1. Set (ymin, ymax) for boxplot region, and yP and yN for p-vals and n's.
   ## -------------------------------------------------------------------------------

   ymin = getMinOfListData(dataList)
   
   ## The user can set the ymax with ymaxBoxplot, else ymax is computed here.                   
   ## ymax is determined by the max of the data, and how much is printed
   ## on the plot; i.e. 1)p-vals and n's 2)either p-vals or n's 3)neither.

   if( is.na(ymaxBoxplot) || is.na(yP) || is.na(yN) ){
      ymaxData = getMaxOfListData(dataList)
      ymaxResults = computeYmaxBoxplot(dataList, ymaxData, betweenGroupTestsCompute, 
         printNs, addToYmax1, addToYmax2)
   }

   if( is.na(ymaxBoxplot) ){ ymax = ymaxResults[[1]] }
   else { ymax = ymaxBoxplot } 
  
   if( is.na(yP) ){ yP = ymaxResults[[2]] }
   if( is.na(yN) ){ yN = ymaxResults[[3]] }

   ## ------------------------------------------------------------------------
   ## 0. Determine the x-points for each group for plotting the data points
   ## ------------------------------------------------------------------------

   ## Get the x-spacings for each group being plotted.

   xGroupPoints = getXAxisAtPoints(numGroups) 
   xGroup = matrix(NA, nrow=numGroups, ncol=numCats)

   for(i in 1:numGroups){

      ## Add the x-category-positions to the x-spacings (xGroupPoints) to get
      ## the x plotting positions for each group for the entire plot.

      xGroup[i,] = 1:numCats + xGroupPoints[i]
   }

   ## ------------------------------------------------------------------------------
   ## 1. Plot the Boxplots
   ## ------------------------------------------------------------------------------

   ## -------------------------
   ## Plot the first group
   ## -------------------------

   ## If legend, make space on the right for it by adding an empty x-category

   if(legendInclude){
      dataFirstGroup = cbind(dataGroup1, rep(NA, numSubjectsGroup1))
      xPointsFirstGroup = c(xGroup[1,], numCats+1)
   }
   else{
      dataFirstGroup = dataGroup1
      xPointsFirstGroup = xGroup[1,]
   }

   boxplot(dataFirstGroup, border=boxColorVector[1], ylab=ylabel, 
      xlab=xlabel, at=xPointsFirstGroup, boxwex=boxWidth, outline=boxOutliers,
      ylim=c(ymin,ymax), axes=F, cex.lab=xylabelsCEX, font.lab=xylabelsFont,
      boxlty=boxlty, boxlwd=boxlwd, medlty=medlty, medlwd=medlwd, medpch=medpch, medcex=medcex)

   ## ---------------------------------------------------
   ## Add the remaining groups, if any, to the boxplot 
   ## ---------------------------------------------------

   if(numGroups > 1){   
     for(i in 2:numGroups) {

       dataNextGroup = dataList[[i]]
    
       boxplot(dataNextGroup, border=boxColorVector[i], at=xGroup[i,], 
          add=TRUE, boxwex=boxWidth, axes=F, outline=boxOutliers,
          boxlty=boxlty, boxlwd=boxlwd, medlty=medlty, medlwd=medlwd, medpch=medpch, medcex=medcex)
    }
   }

   ## -----------------------------------------------------------------------------
   ## 2. Plot the Points
   ## -----------------------------------------------------------------------------

   if(printPoints){
      
      ## -------------------------------------
      ## All Points will have the same color
      ## -------------------------------------

      ## If pointColor, pointChar, or pointCEX is NA, set the Specific var to 1.
           
      if(is.vector(pointColor) && length(pointColor)==1 ) {
         if(!is.na(pointColor)){ pointColorSpecific = pointColor}
         else{ pointColorSpecific = 1 }
      }
      if(is.vector(pointChar) && length(pointChar)==1 ) {
         if(!is.na(pointChar)){ pointCharSpecific = pointChar }
         else{ pointCharSpecific = 1 }
      }
      if(is.vector(pointCEX) && length(pointCEX)==1 ) {
         if(!is.na(pointCEX)){ pointCEXSpecific = pointCEX }
         else{ pointCEXSpecific = 1 }
      }

      for(i in 1:numGroups){

         ## -----------------------------------------------------------------------------------
         ## All points w/i a group will have the same color, Groups can have different colors
         ## -----------------------------------------------------------------------------------

         if(is.vector(pointColor) && length(pointColor)>1) pointColorSpecific = pointColor[i]
         if(is.vector(pointChar) && length(pointChar)>1) pointCharSpecific = pointChar[i]
         if(is.vector(pointCEX) && length(pointCEX)>1) pointCEXSpecific = pointCEX[i]

         ## ------------------------------------------------------------------------------------
         ## All subjects within a group for a given x-category will have the same color,
         ## Points w/i a group can have different colors for each x-category, 
         ## Groups can have different colors 
         ## rows = groups, columns = x-category
         ## ------------------------------------------------------------------------------------

         if(is.matrix(pointColor)) pointColorSpecific = pointColor[i,]
         if(is.matrix(pointChar)) pointCharSpecific = pointChar[i,]
         if(is.matrix(pointCEX)) pointCEXSpecific = pointCEX[i,]

         dataGroup = dataList[[i]]
         numSubjectsGroup = dim(dataGroup)[1]

         for(j in 1:numSubjectsGroup){

            ## ---------------------------------------------------------------------
            ## Colors can be specified for subjects w/i an x-category w/i a group;
            ## In this case, a color can be specified for every point plotted, no
            ## color sharing here.  Here, it is assumed that the input is a list of
            ## matrices.  
            ## ---------------------------------------------------------------------

            if(is.list(pointColor)) pointColorSpecific = pointColor[[i]][j,]
            if(is.list(pointChar)) pointCharSpecific = pointChar[[i]][j,]
            if(is.list(pointCEX)) pointCEXSpecific = pointCEX[[i]][j,]

            points(jitter(xGroup[i,], pointJitter), dataGroup[j,], 
               pch=pointCharSpecific, col=pointColorSpecific, cex=pointCEXSpecific)

         }
      }
   }

   ## --------------------------------------------------------------------------------
   ## 3. Label the X and Y Axes 
   ## --------------------------------------------------------------------------------

   ## ----------------------------------
   ## X-Axis
   ## ----------------------------------
 
   ## -------------------------------------------------------------------
   ## Set x-axis tick mark labels
   ##
   ## If user does not supply labels for the tick marks on the x-axis,
   ## Label the x-axis ticks w/the column names from a group's dataframe
   ## eg. labels could be:  c("TNFa","IFNg","IL2")
   ## -------------------------------------------------------------------

   if(length(xAxisLabels)==1){if(is.na(xAxisLabels)){
      xAxisLabels = colnames(dataGroup1)
   }}    

   ## ---------------------------------------------------------------
   ## xAxisLabels are a vector, rather than matrix, print one line
   ## ---------------------------------------------------------------
   if(is.null(dim(xAxisLabels)) || dim(xAxisLabels)[1] == 1){
      axis(1, at=1:length(xAxisLabels), labels=xAxisLabels, cex.axis=xAxisCEX, 
         font.axis=xAxisFont, las=xAxisRotation)
   }

   ## ------------------------------------------------------------------------
   ## xAxisLabels are a matrix, print multiple label rows on x-axis 
   ## eg. printing 4 row binary combos for cytokine presence; 
   ## i.e. one combo = (0,0,0,1) 
   ## ------------------------------------------------------------------------
   else{
      numRowsLabels = dim(xAxisLabels)[1]
      numColsLabels = dim(xAxisLabels)[2]

      ## -------------------------------
      ## Set-up axis line w/tick marks
      ## -------------------------------
      axis(1, at=1:numCats, labels=rep("", numColsLabels-1), line=0, tick=T)

      ## ------------------------
      ## Print rows of labels
      ## ------------------------
      for(i in 1:numRowsLabels){

         axis(1, at=1:numCats, labels=xAxisLabels[i,2:numColsLabels], line=(i-1), 
            tick=F, cex.axis=xAxisCEX, font.axis=xAxisFont, las=xAxisRotation)

         ## -------------------------------------------------------------------
         ## Add an identifying label on left side for each row of axis labels 
         ## -------------------------------------------------------------------
         mtext(xAxisLabels[i,1], side=1, line=i, at=xAtMtext, adj=0, cex=xMtextCEX, 
            font=xMtextFont)
      }

      ## -----------------------------------------------------------------------------------
      ## Print x-label under the rows of x-axis labels at the tick marks
      ## (probably need to set xlabel=NA, used in first call to boxplot(), step 1 above.)
      ## -----------------------------------------------------------------------------------
      mtext(xMtext, side=1, line=numRowsLabels+2)
   }

   ## ----------------------------------
   ## Y-Axis 
   ## ----------------------------------

      axis(2, cex.axis=yAxisCEX, font.axis=yAxisFont, las=yAxisRotation)

   ## ------------------------------------------------------------------------
   ## 4. Box around entire plot
   ## ------------------------------------------------------------------------

   box(lwd=plotBoxLWD)

   ## ------------------------------------------------------------------------------
   ## 5. Add Legend
   ## ------------------------------------------------------------------------------

   if(legendInclude){

      if(is.na(legendX)){ xLegend = numCats + .75 } else { xLegend = legendX }
      if(is.na(legendY)){ yLegend = ymax-.1*ymax } else { yLegend = legendY } 

      if(legendPoints){

       if(is.na(legendTitle)){
         legend(xLegend, yLegend, legendGroupNames, pch=legendPCH, col=legendColors, 
            cex=legendCEX) }
       else{
          legend(xLegend, yLegend, legendGroupNames, pch=legendPCH, col=legendColors, 
            cex=legendCEX, title=legendTitle) }
      }

      if(legendLines){

       if(is.na(legendTitle)){
         legend(xLegend, yLegend, legendGroupNames, col=legendColors, lty=legendLTY, 
            lwd=legendLWD) }
       else{
         legend(xLegend, yLegend, legendGroupNames, col=legendColors, lty=legendLTY, 
            lwd=legendLWD, title=legendTitle) }
      }
   }

   ## -------------------------------------------------------------------------------
   ## 6. Add Main Title to Plot
   ## -------------------------------------------------------------------------------

   title(mainTitle, line=mainTitleLine, font=mainTitleFont, cex=mainTitleCEX)


   ## -------------------------------------------------------------------------------
   ## 7. Compute Between Group Tests 
   ## -------------------------------------------------------------------------------

   if(betweenGroupTestsCompute){

     if(numGroups == 2){
       
        ## ---------------------------------------------------------
        ## Compute 2 Group Comparisons -- default is Wilcoxon test
        ## ---------------------------------------------------------

           ## Set default testType to Wilcoxon, since it's the only option coded
           testType="Wilcoxon"

           dataGroup2 = dataList[[2]]

           twoGroupPValues = get2GroupPValues(dataGroup1, dataGroup2, 
              pairs=pairedGroups, testType)

           twoGroupPValues = round(twoGroupPValues, digits=testsRoundDigits)
       
        ## ------------------------------------------
        ## -- Create Text Labels to print on Plot
        ## ------------------------------------------
     
           ## pvalueLabels is a text vector containing result for each category
           ## being plotted on x-axis; example, (p < .01, p = .34, p = NA).

           pvalueLabels = getPValueLabels(twoGroupPValues, testsRoundDigits, pvalueLabel)

           # -- Print p-vals on the plot -----------------
           text(1:numCats, rep(yP, numCats), labels=pvalueLabels, cex=pCEX, font=pFont) 

           # -- Print test title on the plot -----------------
           if(!pairedGroups){testTitle = paste(testType, " Rank Sum Test (unpaired data)", sep="")} 
           else {testTitle = paste(testType, " Signed Rank Test (paired data)", sep="")}          
           mtext(testTitle, line=testTitleLine, cex=testTitleCEX, font=testTitleFont)

      } # if numGroups==2

      if(numGroups > 2){

         ## ---------------------------------------------------------------
         ## Compute Multi-Group Comparisons -- default is Kruskal-Wallis 
         ## ---------------------------------------------------------------
         
         ## Set default testType to Kruskal, since it's the only option coded
         testType="Kruskal"

         multiGroupPValues = getMultiGroupPValues(dataList, testType)
         multiGroupPValues = round(multiGroupPValues, digits=testsRoundDigits)

         ## pvalueLabels is a text vector containing result for each category
         ## being plotted on x-axis; example, (p<.01, p=.34, p=NA).

         pvalueLabels = getPValueLabels(multiGroupPValues, testsRoundDigits, pvalueLabel)

         ## ------------------------------------------
         ## -- Create Text Labels to print on Plot
         ## ------------------------------------------

         # -- Print p-vals on the plot ---------------------
         text(1:numCats, rep(yP, numCats), labels=pvalueLabels, cex=pCEX, font=pFont) 

         # -- Print test title on the plot -----------------
         if( testType=="Kruskal" ){ testTitle = "Kruskal-Wallis Rank Sum Test" }         
         mtext(testTitle, line=testTitleLine, cex=testTitleCEX, font=testTitleFont)
      }

   } # end If betweenGroupTests


   ## ---------------------------------------------------------------------------------
   ## 8. Add N's to Plot 
   ## ---------------------------------------------------------------------------------

   if(printNs){

      ## NsVector is a text vector: 
      ## one group example for 3 cats on x-axis: (3,4,5)
      ## two group example for 3 cats on x-axis: (3/3,4/3,5/1)

      NsVector = getNsOfListData(dataList)
      text(1:numCats, rep(yN,numCats), NsVector, cex=nCEX, font=nFont)

      ## Position the "n =" text. 

      if( numGroups > 1 ){ text(0.50, yN, "n =", cex=nCEX, font=nFont) }
                     else{ text(0.25, yN, "n =", cex=nCEX, font=nFont) }
   }

} # end GroupListBoxplot()

