#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source(paste0(Sys.getenv("HOPMA_PATH"),"/tools.R"))
if(args[4]=="combi"){
	mat=extractExcludingListFromColoredMapCombi(args[1],"CA",as.numeric(args[2]),as.numeric(args[3]),c("0.1","0.001"))
}
if(args[4]!="combi"){
	mat=extractExcludingListFromColoredMapCombi(args[1],"CA",as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]))
}
