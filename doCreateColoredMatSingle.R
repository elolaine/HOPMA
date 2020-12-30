#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

source(paste0(Sys.getenv("HOPMA_PATH"),"/tools.R"))
if(args[3]=="combi"){
	createColoredMat(args[1],"CA",as.numeric(args[2]),0.1,as.numeric(args[4]),as.numeric(args[5]))
	createColoredMat(args[1],"CA",as.numeric(args[2]),0.001,as.numeric(args[4]),as.numeric(args[5]))
}
if(args[3]!="combi"){
createColoredMat(args[1],"CA",as.numeric(args[2]),as.numeric(args[3]),as.numeric(args[4]),as.numeric(args[5]))
}

