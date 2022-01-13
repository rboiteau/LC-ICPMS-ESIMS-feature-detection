#Rene Boiteau, 12/21/21 - File for Batch processing of mzXML files. 

#Updated v2.2 with theoretical isotope patterns in plots.

library(xcms)

########################   Input settings (Change settings here)

timerange=c(5*60,20*60) #time range in seconds
background<-500 #minimum signal intensity
minICPheight<-500 #defines minimum ICP height to be detected as a peak in the LCICPMS data.
peakwidth<-30 #expected peak width in seconds (larger number is less stringent for discarding false positives)

slope<-c(0.5,2) #minimum normalized slope cutoff.(Light/Heavy)
correlation<-0.7 #minimum correlation cutoff.

eicwidth<-0.003 #m/z width of extracted ion chromatogram (+/- Da)

#If using a constant offset instead of internal standard to align ICP and LC data:
offset <- -32 # define offset in seconds here.
filetype<-".mzML"

######################## Files are specified here

#Set  ESIMS file
mzdatafile<-'./161221_tutorial_soil.mzML' #Can use any file format that can be imported by XCMS
mzxcms <- xcmsRaw(mzdatafile,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)

#Set ICPMS file
ICPMSfile<-'./161220_tutorial_soil.csv' # file must be a ';' delimited csv of time and intensity for each isotope 'M' of interest, with column titles in the format:'Time M' and 'M'

#Set isotope pattern file and element to look for in ICPMS data
isotopepattern<-'./Cuisotope.csv'
element<-'63Cu'

#Next function will read in ICPMS file. Must be adapted to ICPMS file format (this is for iCAPq)
ICPread<-function(ICPMSfile,element){
  ICPdata<-read.csv(ICPMSfile, header = TRUE, sep = ",", strip.white=TRUE)
  ICP<-cbind(ICPdata[,paste('Time.',element,sep="")],ICPdata[,paste('X',element,sep="")])
  return(ICP)
  }



savefile<-paste(gsub(paste("^./(.*)",filetype,sep=""),"\\1",mzdatafile),"_",gsub("^./(.*).csv","\\1",isotopepattern),sep="") #by default, the pdf printout will be named after the LCESIMS file

##Read in all isotopes
all_elements<-read.table('Isotope_masses.txt',sep="\t",header=TRUE)
epattern<-c()
for(e in unique(all_elements$Name)){
  curr<-all_elements[which(all_elements$Name==e),]
  if(nrow(curr)>1){
    for(i in 2:nrow(curr)){
      diff<-curr$Mass[i]-curr$Mass[1]
      rat<-curr$Abund.[1]/curr$Abund.[i]
      epattern<-rbind(epattern,c(e,diff,rat))
    }  
  }
}
epattern=data.frame(epattern)
colnames(epattern)=c('element','diff','rat')

######################## Functions are defined here

###########Function to find isotope matches

isotopehunter8<-function(mzxcms,isotopefile,scantime){
  
  starttime<-Sys.time()
  
  if(scantime[1]=="all"){
    startscan<-1
    endscan<-length(mzxcms@scantime)
  }else{
    startscan<-which(mzxcms@scantime>scantime[1])[1]
    endscan<-tail(which(mzxcms@scantime<scantime[2]),1)
  }
  
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  
  uppermass<-pattern[-1,2]+pattern[-1,4]-pattern[1,2]
  lowermass<-pattern[-1,2]-pattern[-1,4]-pattern[1,2]
  upperratio<-pattern[-1,3]*pattern[-1,5]/pattern[1,3]
  lowerratio<-pattern[-1,3]/pattern[-1,5]/pattern[1,3]
  nisotope<-length(uppermass)
  
  columnnames<-c(sapply(pattern[,1],FUN=function(x) c(paste(x,'mass'),paste(x,'intensity'))))
  
  results<-matrix(0,ncol=(4+nisotope*2),nrow=1E6)
  r<-1
  colnames(results) = c('scan','time', columnnames)
  
  
  mzint<-20
  mzbuffer<-uppermass[length(uppermass)]
  
  for(i in startscan:endscan) {
    tscan<-getScan(mzxcms,i)
    tscan<-tscan[which(tscan[,2]>background),]
    
    #print(paste(i,'of',endscan))
    
    first<-tscan[1,1]
    last<-first+mzint
    
    while(last<tscan[nrow(tscan),1]){
      
      scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
      
      if(length(scan)>2){
        final<-length(which(scan[,1]<(last-mzbuffer)))
        for(j in 1:final){
          x<-scan[j,1]
          y<-scan[j,2]
          isotopes<-matrix(NA,ncol=nisotope)
          k<-1
          while(k<nisotope+1){
            isotopes[k]<-(which(
              scan[,1]>(x+lowermass[k]) & 
                scan[,1]<(x+uppermass[k]) & 
                scan[,2]>(y*lowerratio[k]) & 
                scan[,2]<(y*upperratio[k])
            )[1])
            if(is.na(isotopes[k])){k<-(nisotope+2)
            } else {k<-k+1
            }
          }
          if(k==(nisotope+1)){
            results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
            r<-r+1
          }
          
        }
      }
      first<-last-mzbuffer
      last<-last+mzint
    }
    scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
    if(length(scan)>2){
      for(j in 1:nrow(scan)){
        x<-scan[j,1]
        y<-scan[j,2]
        isotopes<-matrix(NA,ncol=nisotope)
        k<-1
        while(k<nisotope+1){
          isotopes[k]<-(which(
            scan[,1]>(x+lowermass[k]) & 
              scan[,1]<(x+uppermass[k]) & 
              scan[,2]>(y*lowerratio[k]) & 
              scan[,2]<(y*upperratio[k])
          )[1])
          if(is.na(isotopes[k])){k<-(nisotope+2)
          } else {k<-k+1
          }
        }
        if(k==(nisotope+1)){
          results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
          r<-r+1
        }
      }
    }
  }
  
  print(Sys.time()-starttime)
  
  return(results[1:(r-1),])
}

########Simplify list:

featurebinning<-function(results,mzxcms,background,peakwidth){
  
  masslist2<-c()

      #create index row for results
      results2<-results[order(results[,4],decreasing=TRUE),]
      
      #Remove elements that don't appear at least 2x
      results2<-results2[!isUnique(round(results2[,5],digits=3)),]
      
      #Remove duplicate masses, selecting max intensity scan
      results2<-results2[!duplicated(round(results2[,5],digits=3)),]
      
      #Calculate correlation coefficient for each element
      slope<-0
      rsquare<-0  
      
      if(class(results2)=='matrix'){
        if(nrow(results2)>1){
          for(m in 1:nrow(results2)){
            EIC1<-rawEIC(mzxcms,mzrange=c(-eicwidth,eicwidth)+results2[m,3],rtrange=c(-peakwidth,peakwidth)+results2[m,2])
            EIC2<-rawEIC(mzxcms,mzrange=c(-eicwidth,eicwidth)+results2[m,5],rtrange=c(-peakwidth,peakwidth)+results2[m,2])
            peakcorr<-lm(EIC1$intensity ~ EIC2$intensity)
            slope[m]=coef(peakcorr)[2]
            rsquare[m]=summary(peakcorr)$r.square
          }
        }
      }
      results3<-cbind(results2,slope,rsquare)
      
      return(data.frame(results3))
}

########### Function to align MS data based on cyanocobalamin peak in data
MSalign <- function(mzxcms,ICPMSfile,timerange){
  #Get retention time of Co peak (cobalamin) in ICPMS data
  ICP<-ICPread(ICPMSfile,'59Co')
  ICPComax<-ICP[which.max(ICP[,2]),1]
  
  #Get retention time (in seconds) of Co peak (cobalamin) in Orbitrap data
  B12profile<-rawEIC(mzxcms,c(678.285,678.295))
  B12profile<-do.call(cbind,B12profile)
  
  B12maxscan<-B12profile[which.max(B12profile[,2]),1]
  OrbiComax<-mzxcms@scantime[B12maxscan]
  
  offset<-ICPComax-OrbiComax
  
  #Use this to view data
  par(mfrow=c(2,1))
  plot(ICP[,1],ICP[,2],type='l',xlim=timerange)
  abline(v=ICPComax, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(678.285,678.295),rtrange=timerange)
  abline(v=ICPComax, lty=3,col='cyan')
  
  par(mfrow=c(2,1))
  plot(ICP[,1]-Offset,ICP[,2],type='l',xlim=timerange)
  abline(v=ICPComax-offset, lty=3,col='cyan')
  plotEIC(mzxcms,mzrange = c(678.285,678.295),rtrange=timerange)
  abline(v=ICPComax-offset, lty=3,col='cyan')
  title(paste('offset = ',offset))
  
  return(Offset)
}

########### Automatically find peaks in LC-ICPMS data
findICPpeaks<-function(ICPMSfile,element,timerange,minheight){
  ICP<-ICPread(ICPMSfile,element)
  ICP_base<-ksmooth(ICP[,1],ICP[,2],kernel="normal",bandwidth=120)
  ICP_peak<-ksmooth(ICP[,1],ICP[,2],kernel="normal",bandwidth=10)
  
  peaks<-which(diff(sign(diff(ICP_peak[[2]])))==-2)+1
  truepeaks<-peaks[which((ICP_peak[[2]][peaks]-ICP_base[[2]][peaks])>minheight)]
  
  plot(ICP[,1]/60,ICP[,2],type='l',
       xlim=timerange/60,
       xlab="Retention Time (min)",
       ylab=paste(element," Intensity (cps)"),
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1)
  title(paste('  LC-ICPMS peaks:', element))
  abline(v=ICP[truepeaks,1]/60, lty=3,col='black')
  
  
  return(ICP[truepeaks,1])
}

################# Plotting functions

#Isotope plotter for final report. Makes one plot with EIC's of every isotope found, scaled to the 'same' intensity.
isotopeplotter<-function(mzxcms,pattern,lowmass,timerange,maxscan){
  
  plotrange<-pattern[,2]-pattern[1,2]+lowmass
  isotoperatio<-pattern[,3]/pattern[1,3]
  namerange<-as.vector(pattern[,1])
  
  colors<-c('blue4','darkorange2','burlywood4','black','red')
  usedcolors<-colors[1:ncol(pattern)]
  
  maxy<-0
  times<-mzxcms@scantime
  
  for(i in 1:nrow(pattern)){
    mzrange<-c(-eicwidth,eicwidth)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    EIC1<-unlist(EIC[2])/isotoperatio[i]
    newmax<-max(EIC1[which(times>timerange[1]&times<timerange[2])])
    maxy<-max(c(maxy,newmax))
  }
  
  
  mzrange<-c(-eicwidth,eicwidth)+lowmass
  EIC<-rawEIC(mzxcms,mzrange)
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylim=c(0,maxy*1.1),
       ylab='Scaled Intensity',
       xlab='Retention Time (min)',
       col=usedcolors[1],
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  
  for(i in 2:ncol(pattern)){
    mzrange<-c(-eicwidth,eicwidth)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    lines(times,EIC[[2]]/isotoperatio[i],col=usedcolors[i],lwd=1)
  }
  
  abline(v=maxscan, lty=3,col='red')
  title(paste('EIC = ',toString(round(plotrange,digits=4))),line=0,adj=0,cex.main=1)
  legend('topright',bty="n",namerange,lwd=2,col=usedcolors,cex=0.7,horiz=TRUE)
  
}

#EIC plotter for final report.
EICplot<-function(mzxcms,mass,title,timerange,maxscan){
  mzrange<-c(-eicwidth,eicwidth)+mass
  EIC<-rawEIC(mzxcms,mzrange)
  times<-mzxcms@scantime
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylab='Intensity',
       xlab='Retention Time (min)',
       col='black',
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  title(paste(title,'EIC = ',round(mass,digits=4)),line=1,adj=0,cex.main=1)
  
  abline(v=maxscan, lty=3,col='red')
}

#Mass spectra plotter for final report.
MSplotter<-function(mzxcms,results,i,pattern){
  
  #determine correct isotope pattern
  theorymass<-pattern[,2]-pattern[1,2]+results[i,3]
  theoryratio<-pattern[,3]/pattern[1,3]*results[i,4]
  
  isotopes<-matrix(results[i,-1:-2],byrow=TRUE,ncol=2)
  massrange<-c(mean(isotopes[,1])-10,mean(isotopes[,1])+10)
  scanrange<-getScan(mzxcms,results[i,1],mzrange=massrange)
  
  plot(scanrange,type='h',ylim=c(0,max(scanrange[,2])*1.1),bty='n',ylab='Intensity',xlab='m/z')
  lines(theorymass,theoryratio,type='h',col='cornsilk3',lwd=3)
  lines(isotopes,type='h',col='red')
  title(paste('MS: ',toString(round(results[i,2]/60,digits=2)),' min'),line=1,adj=0,cex.main=1)
}

#ICPMS data plotter for final report.
ICPplotter<-function(ICPMSfile,element,timerange,offset){
  
  #Next line must be changed if ICP data format differs from iCAPq export.
  ICP<-ICPread(ICPMSfile,element)
  ICP[,1]<-ICP[,1]+offset
  ICP<-ICP[which(ICP[,1]>timerange[1] & ICP[,1]<timerange[2]),]
  
  plot(ICP[,1]/60,ICP[,2],type='l',
       xlab="Retention Time (min)",
       ylab=paste(element," Intensity (cps)"),
       bty='n',
       cex.axis=1,
       lwd=1,
       cex.lab=1)
  title(paste(element,'LC-ICPMS'),line=1,adj=0,cex.main=1)
  
}

#Final plot function for generating report for isotope pattern. Includes Apo form MS spectra as last panel
Reportplotter_ICP<-function(mzxcms,mass,results,isotopepattern,timerange,i,metalpeaks,ICPMSfile,element,offset){
  
  mzxcms@scantime<-mzxcms@scantime/60
  timerange_m<-timerange/60
  readpattern<-read.csv(isotopepattern, header = TRUE, sep = ",")
  pattern<-readpattern[which(readpattern[,6]=='Y'|readpattern[,6]=='N'),]
  apo<-readpattern[which(readpattern[,6]=='O'),]
  
  nplot<-3+nrow(apo)
  par(mfrow=c(nplot,1))
  par(mar=c(4,4,3,3))
  
  metalpeaks<-metalpeaks/60
  ICPplotter(ICPMSfile,element,timerange,offset)
  abline(v=metalpeaks, lty=3,col='black')
  title(savefile)
  
  
  isotopeplotter(mzxcms,pattern,mass,timerange_m,(results[i,2]/60))
  abline(v=metalpeaks, lty=3,col='black')
  
  #title(savefile,outer=TRUE)
  title(paste(' R-squared: ',toString(round(results$rsquare[i],digits=3))),line=-1,adj=0.0,cex.main=1)
  title(paste(' Slope: ',toString(round(results$slope[i],digits=3))),line=-2,adj=0.0,cex.main=1)
  title(paste(' d m/z: ',toString(round((results[i,5]-results[i,3]),digits=4))),line=-3,adj=0.0,cex.main=1)
  
  if(nrow(apo)>0){
    for(m in 1:nrow(apo)){
      EICplot(mzxcms,mass+apo[m,2],apo[m,1],timerange_m,(results[i,2]/60))
      abline(v=metalpeaks, lty=3,col='black')
    }
  }
  
  MSplotter(mzxcms,as.matrix(results[,1:(ncol(results)-2)]),i,pattern)
}

# Initialize Pipeline for algorithm search, correlation analysis, and report generation.

Isotope_pipeline <- function(isotopepattern,mzxcms,timerange,background,peakwidth,savefile){
  
  pattern<-read.csv(isotopepattern, header = TRUE, sep = ",")
  ratio<-pattern[2,3]/pattern[1,3]
  
  #Search for all matches to the defined isotope pattern
  results<-isotopehunter8(mzxcms,isotopepattern,timerange)
  print(paste('Results:',nrow(results)))
  
  
  if(class(results)=='matrix'){
  
    #Bin results by mass and remove patterns that only appear 1x:
    binnedresults<-featurebinning(results,mzxcms,background,peakwidth)
    binnedresults$slope<- binnedresults$slope*ratio
    
    #Filter clean results based on R-squared value (should be close to 1 for best hits)
    cleanresults<-binnedresults[which(binnedresults$rsquare>correlation),]
    
    #Filter clean results based on normalized slope (should be 1 if peak ratio really matches isotope pattern)
    cleanresults<-cleanresults[which(cleanresults$slope<slope[2] &
                                         cleanresults$slope>slope[1]),]
    
    print(paste("Cleaned results:",nrow(cleanresults)))
    
    #Generate report: Fix this cleanresults part at some point...Right now, wont generate report if just one hit...
  
    if(nrow(cleanresults)>0){
      
      ##### ICPMS Data Processing
      #Align ICPMS and Orbitrap data. Use multiMSalign to align based on cyanocobalamin
      #offset <- MSalign(mzxcms,ICPMSfile,timerange)
      
      
      #Define LC-ICPMS peaks. 
      ##This function takes in an ICP file name (ICPMSfile), a string that matches the column headers within the ICP data file (e.g. '56Fe'), a time range in seconds ('timerange'), and a minimum ICP intensity 'minICPheight'. It returns a list of retention times with that peak height.
      metalpeaks<-findICPpeaks(ICPMSfile,element,timerange,minICPheight)+offset
      title(paste('Offset = ',offset, 's'),line=1,adj=0,cex.main=1)
      
      #Will generate pdf w/ title savefile.
      pdf(paste(savefile,'.pdf'))
      
      #First, creates QC page:
      par(mfrow=c(3,1))
      #Histogram of Rsquared distribution
      hist(binnedresults$rsquare,breaks=seq(0,1,0.05),xlab='R-squared',main='R-squared Distribution')
      abline(v=correlation,col='red')
      #Histogram of Slope distribution
      hist(cleanresults$slope,breaks=seq(slope[1],slope[2],0.05),xlab='Slope (light/heavy)',main='Normalized Slope Distribution')
      abline(v=slope,col='red')
      
      ## Plot m/z and intensity ratio of all results vs. clean results..
      plot(results[,5]-results[,3],results[,4]/results[,6],col='lightgray',xlab='Mass difference',ylab='Intensity ratio (light/heavy)')
      points(binnedresults[,5]-binnedresults[,3],binnedresults[,4]/binnedresults[,6],col='black')
      points(cleanresults[,5]-cleanresults[,3],cleanresults[,4]/cleanresults[,6],col='red')
      legend('top',legend=c('All','Binned','Cleaned'),col=c('lightgray','black','red'),pch=1, xpd=TRUE, horiz=TRUE, inset=c(0,-.2),bty='n')
      points(epattern$diff,epattern$rat,col='blue')
      text(as.numeric(epattern$diff),as.numeric(epattern$rat),labels=c(epattern$element),cex=0.6,pos=1,col='blue')
      
      masslist<-cleanresults[,3]
      
      
      for(k in 1:length(masslist)){
        Reportplotter_ICP(mzxcms,masslist[k],cleanresults,isotopepattern,timerange,k,metalpeaks,ICPMSfile,element,offset)
        
      }
      
      dev.off()
      
      #write.table(cleanresults,paste(savefile,'.txt'),sep="\t")
    }
    return(cleanresults)
  }
}

######################### To run:


Isotope_pipeline(isotopepattern,mzxcms,timerange,background,peakwidth,savefile)

  
########### To run w/ additional conditions:
  
  
#Set isotope pattern file and element to look for in ICPMS data
isotopepattern<-'./Feisotope.csv'
element<-'56Fe'
  
savefile<-paste(gsub(paste("^./(.*)",filetype,sep=""),"\\1",mzdatafile),"_",gsub("^./(.*).csv","\\1",isotopepattern),sep="") #by default, the pdf printout will be named after the LCESIMS file
  
Isotope_pipeline(isotopepattern,mzxcms,timerange,background,peakwidth,savefile)

