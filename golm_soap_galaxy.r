#!/usr/local/public/bin/Rscript --vanilla --slave --no-site-file
# golm_soap_galaxy.r version="20150408"
#created by Yann GUITTON 

#Redirect all stdout to the log file
log_file=file("golmsearch.log", open = "wt")
sink(log_file)
sink(log_file, type = "out")


library(batch) #necessary for parseCommandArgs function

source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
}

#Import the different functions
source_local("searchsources.r")

listArguments = parseCommandArgs(evaluate=FALSE) #interpretation of arguments given in command line as an R list of objects
print(listArguments)


#saving the name of the function in a variable thefunction
thefunction = listArguments[["xfunction"]]
listArguments[["xfunction"]]=NULL #delete from the list of arguments


if (thefunction == "golmsearch") { #only CDF2RData use library
    library(metaMS)
	# gs<-golmsearch(mspfilevar=listArguments[["imagemsp"]],percentvar=1,rivar=listArguments[["rivar"]],riWindowvar=listArguments[["riwindowvar"]], columnvar=listArguments[["columnvar"]])
    cat(listArguments[["mspfilevar"]])
    print(listArguments[["mspfilevar"]])
    gs<-golmsearch(mspfilevar=listArguments[["mspfilevar"]])
}

#delete the parameters to avoid the passage to the next tool in .RData image
rm(listArguments)

#saving R data in .Rdata file to save the variables used in the present tool
save.image(paste(thefunction,"RData",sep="."))