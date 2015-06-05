#searchsources.r version="20150306"
#searchsources.r for Galaxy Workflow4Metabo 20150408
#version 1.0
#Author Guitton Yann #R SOAP for Golm Metabolome Web Service 2015 Y. Guitton IDEALG
#For Galaxy usage

golmsearch<-function(mspfilevar,percentvar=1,rivar=1500,riWindowvar=3000, columnvar="VAR5", maxhitsvar="all"){
    require(metaMS)
    DB<-read.msp(file=mspfilevar)
    ions<-DB


    #For Golm Metabolome search
    #attention not on Galaxy server!
    # dir.create("c:/Temp/GOLMsearchGalaxy/")
    # setwd("c:/Temp/GOLMsearchGalaxy/")

    #OPtion avec GET (httr) et xml parsing  
    library(httr)
    library(XML)
    ans<-NULL #for hwrite
    # pcgrp=48
 print("Step1")
 cat("Step1")
    for (pcgrp in 1:length(ions)){
        file=file.path(getwd(),paste("GOLMresult_spectra_",ions[[pcgrp]]$Name,".tsv", sep=""))
        percent=percentvar
        ri=rivar
        riWindow=riWindowvar
        column=columnvar
  print("Step2")
 cat("Step2")
        newpcgroup<-cbind(round(ions[[pcgrp]]$pspectrum[,"mz"],0),round((ions[[pcgrp]]$pspectrum[,"intensity"]*100)/max(ions[[pcgrp]]$pspectrum[,"intensity"]),1))
        spectrumpcgroup=paste0(as.vector(t(subset(newpcgroup,newpcgroup[,2]>percent))), collapse="%20", sep="")
        xmlRES<-GET(paste("http://gmd.mpimp-golm.mpg.de/webservices/wsLibrarySearch.asmx/LibrarySearch?ri=",ri,"&riWindow=",riWindow,"&AlkaneRetentionIndexGcColumnComposition=",column,"&spectrum=",spectrumpcgroup), sep="")
        xmlDOC<-content(xmlRES)
        xmlROOT<-xmlRoot(xmlDOC)
        a<-xmlToList(xmlDOC)
        if (a[[1]]=="success"){
            HIT<-which(names(a)=="Results")
            print (paste("The number of hits is ", length(HIT)))
            colnamesHIT<-names(a[HIT[1]]$Results)
            options(stringsAsFactors = FALSE)
            df <- data.frame(matrix(unlist(a[HIT]), nrow=length(HIT), byrow=T))
            options(stringsAsFactors = TRUE)
            colnames(df)<-colnamesHIT
            
            #keep a list with unique metaboliteID with the lowest euclideanDistance if several time the same metaboliteID
            #order by EuclideanDistance
            df<-df[order(as.numeric(df[,"EuclideanDistance"])),]
           
            
            #remove duplicate and keep median value for scores
                #get levels
                list_ids<-unique(df$analyteName)
                dfnew<-vector("list", length(list_ids))
                for (i in 1:length(list_ids)){
                  
                    tmp<-subset(df, df[,"analyteName"]==list_ids[i])
                    
                    if (dim(tmp)[1]>1){
                        
                        mini<-lapply(tmp[,4:9], FUN=function(x)min(as.numeric(x))) #can take median also?
                        mini<-lapply(mini, FUN=function(x)round(as.numeric(x),2))
                        
                        dfnew[[i]]<-unlist(c(tmp[1,1:3],unlist(mini),tmp[1,10:12]))
                    }else
                    {
                       
                       mini<-lapply(tmp[,4:9], FUN=function(x)round(as.numeric(x),2))
                       dfnew[[i]]<-unlist(c(tmp[1,1:3],unlist(mini),tmp[1,10:12]))
                    }
                    
                    
                }
                options(stringsAsFactors = FALSE)
                df1 <- data.frame(matrix(unlist(dfnew), nrow=length(list_ids), byrow=T))
                colnames(df1)<-colnamesHIT
                df1<-df1[order(as.numeric(df1[,"EuclideanDistance"])),]
                options(stringsAsFactors = TRUE)
                #reodrer columns
                df<-df1[,c(11,3,5:9,1,2,4,10,12)]
            
        }else
        {
            df<-"Empty"
        }
        if (maxhitsvar!="all"){
                df<-df[1:as.numeric(maxhitsvar),]
        }
        write.table(df, file=file, sep="\t", row.names=F)
        cname<-rep("",dim(df)[1])
        cname[1]<-ions[[pcgrp]]$Name
        ans<-rbind(ans,cbind(cname,df))
    }#end for pcgrp


    #next create HTML output
    #to do export as HTMLtable add link to MassBank , Kegg, Chebi?
    filename="GOLM_Result"
    require(hwriter)
    p=openPage(paste(filename,".html",sep=""))
    hwrite(ans, p,row.bgcolor='#ffdc98', row.names=F )
    closePage(p)
    #add pcgrp name for first row and then "" for other hits
    #next step add KEGG Chebi ids for all
}


massbanksearch<-function(mspfilevar, mspoption=TRUE, xsAnnotate=NULL){ #xsAnnotate option just for futur
#MassBank spectrumsearch for W4M
#Yann Guitton 
#Inputs
# either a xsAnnotate object with pcgrp or an msp file generated elswhere
#to do add option to search only some pcgrp not all

#option from xsAnnotate ou from msp  mspoption=TRUE or FALSE (if from xsAnnotate)
if (mspoption==TRUE){
    require(metaMS)
    DB<-read.msp(file=mspfilevar)
    ions<-DB
}else{
  #todo load Rdata from w4m and check if pcgrp are in
    ions<- vector("list", length(xsAnnotate@pspectra))
    for (i in 1:length(xsAnnotate@pspectra)){ions[[i]]<-getpspectra(xsAnnotate,i)[,c("mz","into")]}
}

#For MassBankSearch  not used in GAlaxy?
dir.create("./MBsearchGalaxy/")
setwd("./MBsearchGalaxy/")

#temporary file creation because idon't know how to change the perl code to use mass spectra (pcgrp ) as input
for (i in 1:length(ions)){
    write.table(ions[[i]]$pspectrum[,1:2], sep=" ", file=paste("spectra_",ions[[i]]$Name,".txt", sep=""), row.names=F, col.names=FALSE)
    #option change instGC hard coded by a dynamically generated file with choosen instrument 
    polarity="both"
    instfile="instGC"
    Res1 <- system(paste0('perl searchSpectrumMB.pl  \"', paste(file.path(getwd(),paste("spectra_",ions[[i]]$Name,".txt\"", sep="")), polarity, paste(instfile, sep=""), sep=" "), collapse=" " ), intern=TRUE)
    # unlink(paste(file.path(getwd(),paste("spectra_",i,".txt", sep=""))))
    write.table(Res1, file=file.path(getwd(),paste("MBresult_spectra_",ions[[i]]$Name,".txt", sep="")), sep="\t", row.names=F, col.names=FALSE)
 }
      

    
# Read MBresult files in order to get details about all hits
# Get MassBank Records
# %rec = &MassBank::getRecordInfo($msbk, @acc);
#not used in Galaxy?
setwd("./MBsearchGalaxy/")

options(stringsAsFactors = FALSE)
    Resfiles<-dir(getwd())[grep("MBresult_spectra_",dir(getwd()))]
    MD_IDs<-vector("list", length(ions))
   df1<-vector("list", length(ions))
    for (i in 1:length(Resfiles)){
        # print(paste("I= ",i))
        tmp<-read.table(Resfiles[i], sep="\t", h=F)
        tmp1<-vector("list", dim(tmp)[1])
        for (j in 1:dim(tmp)[1]){
            tmp1[[j]]<-strsplit(tmp[j,1], split="\t")[[1]]
        }
        df <- data.frame(matrix(unlist(tmp1), nrow=dim(tmp)[1], byrow=T))
        pcgroup<-strsplit(Resfiles[i], split="pcgroup_|.txt")[[1]][2]
        df1[[i]]<-cbind(pcgroup,df)
        colnames(df1[[i]])<-c("pcgroup","MBIDs","Name","Formula","Mass","Score")
        
    }
 options(stringsAsFactors = TRUE) 
 
 #Attention no column names in final files
 
 #Use webservice to get Record Info KEGG CHEBI NAME ?CAS?
Res1<-vector("list", length(ions))
for (l in 1: length(ions)){
    Res2<-vector("list", dim(df1[[l]])[1])
     for (k in 1:dim(df1[[l]])[1]){
        tmp<-system(paste0('perl getMBRecordInfo_single.pl ', df1[[l]][k,2], collapse=" " ), intern=TRUE)
        Res2[[k]]$Name<-tmp[grep("NAME_ASCII: ",tmp)]
        Res2[[k]]$Kegg<-tmp[grep("KEGG_IDs: ",tmp)]
        Res2[[k]]$Chebi<-tmp[grep("ChEBI_IDs: ",tmp)]
     } 
     cname<-unlist(lapply(Res2, FUN=function(x)strsplit(x$Name, split="NAME_ASCII: ")[[1]][2]))   
     ckegg<-unlist(lapply(Res2, FUN=function(x)if(length(x$Kegg)!=0) strsplit(x$Kegg, split="KEGG_IDs: ")[[1]][2] else ""))
     cchebi<-unlist(lapply(Res2, FUN=function(x)if(length(x$Chebi)!=0) strsplit(x$Chebi, split="ChEBI_IDs: ")[[1]][2] else ""))
     
     Res1[[l]]<-cbind(df1[[l]],cname,ckegg,cchebi)
}

#to do export as HTMLtable add link to MassBank , Kegg, Chebi?
filename="MassBank_Result"
require(hwriter)
p=openPage(paste(filename,".html",sep=""))
options(stringsAsFactors = FALSE)
ans<-NULL
for (k in 1:length(Res1)){
    ans<-rbind(ans,Res1[[k]])
}
options(stringsAsFactors = TRUE)
ans2<-ans[,1:9]
cunique<-unique(ans2[,1])
for (j in 1:length(unique(ans2[,1]))){
    ans2[grep(cunique[j],ans2[,1])[-1],1]<-""
}
hwrite(ans2, p,row.bgcolor='#ffdc98', row.names=F )
closePage(p)

}