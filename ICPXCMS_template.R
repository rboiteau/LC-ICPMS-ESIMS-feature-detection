#Rene Boiteau, LC-ICPMS-ESIMS Isotope matching algorithm analysis template

library(xcms)
source("ICPXCMS_functions5.8.R")
########################   Input settings


#Set parameters
timerange=c(200,2400) #time range in seconds
background<-500 #background threshold in cps
noise<-12 
peakwidth<-10 
minICPheight<-300 #defines minimum ICP height to be detected as a peak in the LCICPMS data.
########################   Load Data

Orbifile<-'LCESIMS_file.mzXML' #Can use any file format that can be imported by XCMS
ICPMSfile<-'LCICPMS_file.csv' # file must be a ';' delimited csv of time and intensity for each isotope 'M' of interest, with column titles in the format:'Time M' and 'M'
savefile<-paste(gsub(".mzXML","",Orbifile),"_isotopesearch.pdf",sep="") #by default, the pdf printout will be named after the LCESIMS file

########################   Initialize clean up and report function

metalreport <- function(isotopepattern,mzxcms,timerange,background,peakwidth,savefile,ICPdata,metalpeaks,metal,results,apodiff){
  print(paste('Results:',nrow(results)))
  cleanresults<-Removenoise4(results,background,peakwidth)
  print(paste("Cleaned results:",nrow(cleanresults)))
  if(nrow(cleanresults)>0){
    if(apodiff>0){ICPXCMSreport3(savefile,mzxcms,cleanresults,isotopepattern,timerange,ICPdata,metal,metalpeaks,apodiff)}
    else{ICPXCMSreport2(savefile,mzxcms,cleanresults,isotopepattern,timerange,ICPdata,metal,metalpeaks)}
    
  }
  return(cleanresults)
}

########################   Call functions


#Read in Orbitrap file
mzxcms <- xcmsRaw(Orbifile,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)

#Read in ICPMS file
ICPdata<-read.csv(ICPMSfile, header = TRUE, sep = ";", strip.white=TRUE)


# First, Align ICPMS and Orbitrap data. Use multiMSalign to align based on cyanocobalamin
#mzxcms2<-multiMSalign(mzxcms,ICPdata,timerange)

#To align with a constant offset instead, uncomment the following lines:
# offset <- 20 # define offset in seconds here.
# mzxcms2<-mzxcms
# mzxcms2@scantime<-mzxcms@scantime - offset

# Second, define LC-ICPMS peaks. 
##This function takes in an ICP file name ('ICPdata'), a string that matches the column headers within the ICP data file (e.g. '56Fe'), a time range in seconds ('timerange'), and a minimum ICP intensity 'minICPheight'. It returns a list of retention times with that peak height.

peaks_56Fe<-findICPpeaks(ICPdata,'56Fe',timerange,300)
peaks_63Cu<-findICPpeaks(ICPdata,'63Cu',timerange,1000)
peaks_60Ni<-findICPpeaks(ICPdata,'60Ni',timerange,minICPheight)
peaks_79Br<-findICPpeaks(ICPdata,'79Br',timerange,500)

# Third, hunt for isotopes patterns. 
##This function requires LCESIMS data as an xcmsRAW object (mzxcms2, which is created above), an isotope pattern file (e.g. 'Feisotope.csv'), and a time range in seconds ('timerange'). It searches each scan of the MS data for features matching the specified isotope pattern and returns a complete list.

results_fe<-isotopehunter8(mzxcms2,'Feisotope.csv',timerange)
results_cu<-isotopehunter8(mzxcms2,'Cuisotope.csv',timerange)
results_br<-isotopehunter8(mzxcms2,'Brisotope.csv',timerange)
results_ni<-isotopehunter8(mzxcms2,'Niisotope.csv',timerange)

# Fourth, clean up data and generate metal reports. This function requires objects defined above including background and peakwidth parameters (used for peak picking), a file name for the saved pdf (savefile), and finally an optional apo form mass difference (difference between the lightest isotope and the apo mass of interest), which is used to generate an EIC plot of the apo m/z of each feature in the final report. The function bins all of the results determined above by mass, finds chromatographic peakes within each bin (several peak picking algorithms can be used. The default is very lax), and returns a report of all potential isotope features.

cleanresults_fe<-metalreport('Feisotope.csv',mzxcms2,timerange,background,peakwidth,savefile,ICPdata,peaks_56Fe,'56Fe',results_fe,50.9161)
cleanresults_cu<-metalreport('Cuisotope.csv',mzxcms2,timerange,background,peakwidth,savefile,ICPdata,peaks_63Cu,'63Cu',results_cu,60.9139) 
cleanresults_br<-metalreport('Brisotope.csv',mzxcms2,timerange,background,peakwidth,savefile,ICPdata,peaks_79Br,'79Br',results_br,0) 
cleanresults_ni<-metalreport('Niisotope.csv',mzxcms2,timerange,background,peakwidth,savefile,ICPdata,peaks_60Ni,'60Ni',results_ni,55.9212) 

##