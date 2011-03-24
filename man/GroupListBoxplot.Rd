\name{GroupListBoxplot}
\alias{GroupListBoxplot}
\title{A Boxplot Function With Embedded Statistical Tests for Comparing Groups}
\description{
   A function which creates boxplots side-by-side, with points overlaid, to compare groups.  
   Group sizes and p-values from tests comparing groups can be printed on the plot.  
}
\usage{
   GroupListBoxplot(dataList, ymaxBoxplot=NA,addToYmax2=.2,addToYmax1=.1,boxWidth=.10,boxColor=8,
    boxOutliers=TRUE,groupColorVector=c(2,4,5,6,7,9,3,10,11,8,1),boxlty=1,boxlwd=1, 
    medlty=1,medlwd=3,medpch=NA, medcex=NA, legendInclude=TRUE,
    legendGroupNames=paste("Group ", 1:length(dataList),sep=""),legendX=NA,legendY=NA,legendCEX=1,
    legendPCH=1,legendColors=1:length(legendGroupNames),legendLTY=NA,legendLWD=NA,legendTitle=NA,
    legendPoints=TRUE,legendLines=FALSE,printPoints=TRUE,pointChar=1,pointCEX=1,
    pointColor=1:length(legendGroupNames),pointJitter=.25,mainTitle="Boxplots", mainTitleCEX=1,
    mainTitleFont=1,mainTitleLine=1,testTitleCEX=1,testTitleFont=1,testTitleLine=0,xlabel="X Axis",
    ylabel="Y Axis",xylabelsCEX=1,xylabelsFont=1,xAxisLabels=NA,xAxisCEX=1, xAxisFont=1,
    xAxisRotation=1,xMtext="",xMtextCEX=1,xMtextFont=1,xAtMtext=0,yAxisCEX=1,yAxisFont=1,
    yAxisRotation=1,plotBoxLWD=1,testsRoundDigits=2,pCEX=1,pFont=1,yP=NA,pvalueLabel="p",
    betweenGroupTestsCompute=TRUE,pairedGroups=FALSE,printNs=TRUE,nCEX=1,nFont=1,yN=NA)
}
\arguments{
   \item{dataList}{A list of data frames.  Each data frame contains the data for a group.  Each 
     column in the data frame represents a time point or category point on the plot.}
   \item{ymaxBoxplot}{ymax for the boxplot will be set by R's boxplot() fcn, unless set here, 
     default=NA} 
   \item{addToYmax2}{this value is added to ymax if printing p-vals and N's, default=.2} 
   \item{addToYmax1}{this value is added to ymax if only printing p-vals or N's, default=.1}
   \item{boxWidth}{width of the box, default=.10} 
   \item{boxColor}{color of the box, default=8} 
   \item{boxOutliers}{print or suppress outlier points on the plot, default=TRUE} 
   \item{groupColorVector}{Vector of unique colors for a set of up to 11 groups,
     default=c(2,4,5,6,7,9,3,10,11,8,1)}
   \item{boxlty}{Box outline type, default=1}
   \item{boxlwd}{Box outline width, default=1}
   \item{medlty}{median line type, default=1}
   \item{medlwd}{median line width, default=3}
   \item{medpch}{median point character, default=NA}
   \item{medcex}{median point size expansion, default=NA}
   \item{legendInclude}{Print legend on plot, default=TRUE} 
   \item{legendGroupNames}{names to use in the legend items, default=paste("Group ", 1:length(dataList), 
     sep="")} 
   \item{legendX}{legend X location, default=NA} 
   \item{legendY}{legend Y location, default=NA}
   \item{legendCEX}{point size for legend text, default=1} 
   \item{legendPCH}{symbol or character to print next to legend names, default=1}
   \item{legendColors}{colors of the points or lines next to the legend names, 
     default=1:length(legendGroupNames)} 
   \item{legendLTY}{line type to print next to legend names, default=NA} 
   \item{legendLWD}{width of line to print next to legend names, default=NA}
   \item{legendTitle}{legend title, default=NA}
   \item{legendPoints}{print points next to legend names, default=TRUE} 
   \item{legendLines}{print lines next to legend names, default=FALSE}
   \item{printPoints}{overlay the points on the boxplot, default=TRUE}
   \item{pointChar}{symbol or character to plot, default=1} 
   \item{pointCEX}{size of point plotted, default=1} 
   \item{pointColor}{color of point plotted -- can be single value, vector, matrix, or list, 
     default=1:length(legendGroupNames)} 
   \item{pointJitter}{amount to scatter points around group x position, default=.25}
   \item{mainTitle}{text for main title, default="Boxplots"} 
   \item{mainTitleCEX}{size of main title, default=1} 
   \item{mainTitleFont}{regular or bold font, default=1}
   \item{mainTitleLine}{margin line on which to print main title - minimum is 0, default=1}
   \item{testTitleCEX}{size of stat test title, default=1}
   \item{testTitleFont}{regular or bold font, default=1}
   \item{testTitleLine}{margin line on which to print title - minimum is 0, default=0} 
   \item{xlabel}{label for X axis, default="X Axis"}
   \item{ylabel}{label for Y axis, default="Y Axis"} 
   \item{xylabelsCEX}{size of x and y axes labels, default=1}
   \item{xylabelsFont}{regular or bold font for X and Y axes labels, default=1}
   \item{xAxisLabels}{labels for ticks on x axis, default=NA}
   \item{xAxisCEX}{size of labels on x ticks, default=1}
   \item{xAxisFont}{regular or bold font for x tick labels, default=1} 
   \item{xAxisRotation}{horizontal or vertical x tick labels, default=1}
   \item{xMtext}{text to place in x-axis margin, default=""}
   \item{xMtextCEX}{size of text to place in x-axis margin, default=1} 
   \item{xMtextFont}{regular or bold font for text in x-axis margin, default=1} 
   \item{xAtMtext}{x position for text in x-axis margin, size of y tick labels, default=0} 
   \item{yAxisCEX}{size of y tick labels, default=1} 
   \item{yAxisFont}{regular or bold font for y axis tick labels, default=1} 
   \item{yAxisRotation}{horizontal or vertical y tick labels, default=1}
   \item{plotBoxLWD}{line width of box drawn around entire plot, default=1} 
   \item{testsRoundDigits}{number of digits to report in p-value, default=2} 
   \item{pCEX}{size of p-value text, default=1} 
   \item{pFont}{regular or bold font for p-value text, default=1} 
   \item{yP}{y position for p-value text, default=NA} 
   \item{pvalueLabel}{label to use to precede numerical p-value, default="p"}
   \item{betweenGroupTestsCompute}{include group comparison p-values in plots, default=TRUE}
   \item{pairedGroups}{are the groups paired?, default=FALSE}
   \item{printNs}{include sample sizes on plots, default=TRUE}
   \item{nCEX}{size of the sample size text to print, default=1} 
   \item{nFont}{regular or bold text for n's, default=1}
   \item{yN}{y position for the sample size text, default=NA}
   
}
\author{N. Hawkins, Fred Hutchinson Cancer Research Center, Seattle, WA}
\seealso{\code{\link{boxplot}} }
\examples{

# Create Sample dataList
group1DataFrame = as.data.frame(cbind(1:3,4:6))
group2DataFrame = as.data.frame(cbind(4:6,7:9))
dataList = list(group1DataFrame, group2DataFrame) 

# Make the plot
GroupListBoxplot(dataList, xlabel="Cytokine", ylabel="Percent of CD4 Cells", 
   xAxisLabels=c("IFNg","TNFa"), mainTitle="Compare Innate Immune Response", 
   legendGroupNames=c("Group 1","Group 2"))

## -- Adults vs. Neonates Data -----------------------------------

## Marginal Data boxplot

# Get the data
data(marginalDF)
marginalDataSubset = subset(marginalDF, stim=="LPS" & concGroup==3 & cell=="mDC")
dataList = makeDataList(marginalDataSubset, "group", 1:5) 

# Make the group boxplot of marginal data
GroupListBoxplot(dataList, xlabel="Cytokine", ylabel="Percent of All Cells", 
   xAxisLabels=c("TNFa","IL6","IL12","IFNa","AnyMarker"), 
   mainTitle="Stimulation = LPS and Concentration Group = 3 and Cell = mDC", 
   legendGroupNames=c("Adults","Neonates"), pointColor=c(2,4), testTitleCEX=.8, nCEX=.8, 
   pCEX=.8, legendColor=c(2,4), legendCEX=.7)

}





