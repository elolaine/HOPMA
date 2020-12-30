
###########################################################
#####################   convert distances #################
###########################################################

# scale a distance matrix between 0 and 1
# by default the max value (1) is given for distances lower than 3A
# and the min value (0) is given for distances beyond 15A
convertDist<-function(mat,valMin=3,valMax=15){

	mat[mat>valMax]=valMax
	mat[mat<valMin]=valMin
	return((valMax-mat)/(valMax-valMin))
}

convertProb<-function(x,valMin=3,valMax=15){
	return(valMax-x*(valMax-valMin))
}
###########################################################

###########################################################
#######################  read files #######################
###########################################################

# experimental (distance)
# the function reads a distance matrix pre-computed with another tool 
# and then scales it between 0 (very far away) and 1 (contact)
# dmin and dmax are the lower and upper bounds
getDistMat<-function(target,typeAt,dmin,dmax){
	return(convertDist(as.matrix(read.table(paste(target,"_",typeAt,".dist",sep=""))),dmin,dmax))
}

# from the CSV file of the PDB structure
# get the chain ids and resids 
getInfo<-function(target){
	fname=paste0(target,"_CA.csv")
	dat=read.table(fname,sep=",")
	myInfo=paste0(dat[,5],dat[,6])
#		gsub(",[[:space:]]","",toString(x)))})
	return(myInfo)
}

###########################################################

###########################################################
###############  get intervals in matrices  ###############
###########################################################

# return ranges from a vector of numbers
findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

# return ranges from a vector of numbers
# only ranges longer than a certan threshold are retained
findLargeIntRuns <- function(run,cutoff=2){
  rundiff = c(1, diff(run))
  difflist = split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x)>cutoff) paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

# get ranges of zero or negative values
getZeroRanges<-function(mat,z){
	n = dim(mat)[[1]]
	nam = c()
	res=list()
	k = 0 
	# for each line, get all ranges 
	if(z<2){tmp = apply(mat,1,f<-function(x){return(findIntRuns(which(x<=0)))})}
	else{tmp = apply(mat,1,f<-function(x){return(findLargeIntRuns(which(x<=0)))})}
	#print(sort(table(unlist(tmp))))
	# for each line
	for(i in 1:n){
		# for each range
		for(j in tmp[[i]]){
			# accumulate the line in the range data
			if(j %in% nam){res[[which(nam==j)]]=c(res[[which(nam==j)]],i)}
			else{k = k+1; nam=c(nam,j) ; res[[k]] = i}
		}
	}
	names(res)=nam
	if(z<2){res=sapply(res,findIntRuns)}
	else{res=sapply(res,findLargeIntRuns)}
	res[sapply(res, is.null)] <- NULL
	return(res)
}

# get the immediate neighbours of a cell in a matrix, provided that
# they are not verboten and excluding the cell itself
getNeighbours<-function(coord,n,verboten=c()){
	#print(coord)
	i = coord[1]
	j = coord[2]
	# immediate neighbours of i and j (x and y direction for the query cell)
	rangei = seq(max(1,i-1),min(n,i+1))
	rangej = seq(max(1,j-1),min(n,j+1))
	# initialize output list
	res = c()
	for(k in rangei){
		for(l in rangej){
			# if this is not the cell itself but some neighbour
			if(sum(c(k,l)==coord)<2){
				# if the neighbour is not verboten
				if(sum(verboten==toString(c(k,l)))==0){
					# then add the neighbour to the output list
					res=rbind(res,c(k,l))}
			}
		}
	}
	return(res)
}


getRange<-function(i,dmin,dmax,n){
	start1 = max(1,i-dmax)
	end1 = max(1,i-dmin)

	start2 = min(n,i+dmin)
	end2 = min(n,i+dmax)

	return(c(start1:end1,start2:end2))
}

###########################################################

###########################################################
###################  filter matrices  #####################
###########################################################

# given an input matrix and a selection of cells
# return a matrix with zeros in the background
# and ones on a subset of the selection 
# the subset is defined by the points which have 
# neighboring cells (at distance k) below a certain threshold 
# (thresh, zero by default)
getTbrmCarres<-function(mat,x,k,thresh=0){
	mat[x] = 0
	n = dim(mat)[[1]]
	newMat = matrix(0,nc=n,nr=n)
	sel = (apply(x,1,min)>k&apply(x,1,max)<=(n-k)&abs(x[,1]-x[,2])>(2*k))
	x = x[sel,]
	#print(dim(x))
	le = sum(sel)
	for(p in 1:le){
		submat = mat[(x[p,1]-k):(x[p,1]+k),(x[p,2]-k):(x[p,2]+k)]
		cond = sum(submat<=thresh)/length(submat) == 1
		if(cond){
			newMat[(x[p,1]-k):(x[p,1]+k),(x[p,2]-k):(x[p,2]+k)]=1
		}
	}
	return(newMat)
}

# This function binarizes a distance matrix
# the idea is to put 1 when the average value  
# mat is the converted distance matrix (from experimental starting structure)
# x is a 2D-array with x and y indices specifying a selection of cells in the matrix (lower triangle for example)
# k is the size of the window
# cutBot and cutTop are lower and upper bounds 
getTbrmCarresMoy<-function(mat,x,k,cutBot,cutTop=1){
	# get the total number of lines in the distance matrix
	n = dim(mat)[[1]]
	# initialize result matrix
	newMat = matrix(0,nc=n,nr=n)
	# select cells (i,j) within the input pre-selection (lowest triangle for example)
	# such that min(i,j)>=k and max(i,j)<=(n-k) and the |i-j| greater than 2k
	sel = (apply(x,1,min)>=k&apply(x,1,max)<=(n-k)&abs(x[,1]-x[,2])>=(2*k))
	x = x[sel,]
	#print(dim(x))
	# how many cells to we have?
	le = sum(sel)
	# for each cell (i,j) in the selection
	for(p in 1:le){
		# extract a submatrix ([i-k:i+k],[j-k,j+k]) (2D window of size k around the cell)
		submat = mat[(x[p,1]-k):(x[p,1]+k),(x[p,2]-k):(x[p,2]+k)]
		# test whether the submatrix has an average value between the two bounds
		cond = (mean(submat)>=cutBot)&&(mean(submat)<=cutTop)
		if(cond){
			newMat[(x[p,1]-k):(x[p,1]+k),(x[p,2]-k):(x[p,2]+k)]=1
		}
	}
	return(newMat)
}

###########################################################

###########################################################
###################  plot matrices  #####################
###########################################################

plotClusteredMat<-function(lClus,n){
	mat=matrix(0,nc=n,nr=n)
	for(i in 1:length(lClus)){
		mat[lClus[[i]]]=i
	}
	image(1:n,1:n,mat,col=c("white",1:length(lClus)))
}

###########################################################

###########################################################
###################  detect paches  #######################
###########################################################

# detect contiguous patches of 1s in a binarized matrix
# this is a very well known algorithm (found on Wikipedia)
# Connected component labeling, 2-pass algorithm
# https://en.wikipedia.org/wiki/Connected-component_labeling
# each patch is given an id (integer)
ConnCompLab2Pass<-function(mat){
	# select cells where there are 1s
	sel = which(mat==1,arr.ind = TRUE)
	# dimension of the matrix
	n = dim(mat)[[1]]
	# initialize output matrix with zeros
	newMat = matrix(0,nc=n,nr=n)
	c = 0
	k = 0
	lEqui = list()
	# for each selected cell
	for(i in 1:dim(sel)[[1]]){
		# get the indices of the cell
		myPoint = sel[i,]
		# check that it's really positive (this test is actually no necessary)
		if(mat[myPoint[1],myPoint[2]]>0){
			# get the immediate neighbours of the cell (excluding itself)
			myNeigh = getNeighbours(myPoint,n)
			# select the corresponding submatrix in the output matrix
			submat = newMat[myNeigh]
			# if there are only zeros in there
			if(sum(submat)==0){
				# attribute the point to a new patch 
				newMat[myPoint[1],myPoint[2]] = c+1
				# increase current patch id
				c = c+1
			}
			# otherwise
			else{
				# get the minimum patch id over the neighbourhood
				myMin = min(submat[submat>0])
				# attribute this patch id to the cell
				newMat[myPoint[1],myPoint[2]] = myMin
				# get all patch ids in there 
				vals = unique(submat[submat>0])
				# if there is nothing in the equivalence list
				p = length(lEqui)
				if(p==0){
					# then put these ids in the first element of the list
					lEqui[[1]] = vals
				}
				# otherwise
				else{
					# identify the element(s) of the list having some intersection with the ids
					myInd = which(sapply(lEqui,f<-function(x){return(length(intersect(x,vals))>0)}))
					# how many elements?
					nbInd = length(myInd)
					# zero, then just add a new equivalence class
					if(nbInd==0){
						lEqui[[p+1]] = vals
					}
					# otherwise 
					else{
						# get the first equivalence class
						myI = myInd[1]
						# and merge it with the discovered equivalence class
						lEqui[[myI]] = unique(c(lEqui[[myI]],vals))
						# if there are more than one
						if(nbInd>1){
							# for each of the remaining ones
							for(m in 2:nbInd){
								# merge them with the first equivalence class
								myJ = myInd[m]
								lEqui[[myI]] = unique(c(lEqui[[myI]],lEqui[[myJ]]))
								lEqui[[myJ]] = NULL
							}
						}
					}
				}
			}
		}
	}
	# for each equivalence class
	# (the id of the equivalence class will be the new patch id)
	nbG = length(lEqui)
	for(i in 1:nbG){
		# sort the patch ids
		vec = sort(lEqui[[i]])
		le = length(vec)
		# for each patch id
		for(j in 1:le){
			# replcae the patch id of the concerned cells by the class id
			# this ensures we have patch ids starting from 1 and increasing one by one.
			newMat[newMat==vec[j]] = i
		}
	}
	return(list(newMat,lEqui))
}


###########################################################

###########################################################
###################  main functions  ######################
###########################################################

# creates a colored matrix from a PDB file 
# target is the PDB file (code + chain combination), for example 1ex6A
# typeAt is the atom type, for example CA
# le determines the size of the window for the discretization
# cut is a lower bound for the distance matrix values 
createColoredMat<-function(target,typeAt,le,cut,dmin,dmax){
	# get a scaled distance matrix (0=far away, 1=contact)
	expMat = getDistMat(target,typeAt,dmin,dmax)
	n = dim(expMat)[[1]]
	# binarise the distance matrix
	matCarres = getTbrmCarresMoy(expMat,which(lower.tri(expMat),arr.ind=TRUE),le,cut)
	# identify patches
	matClass=ConnCompLab2Pass(matCarres)[[1]]
	# write out the colored matrix giving the definition of the patches (seen as "classes" here)
	write.table(matClass,paste0(target,"_",typeAt,"_",le,"_",cut,".class"))
	# plot the colored matrix
	pdf(paste0(target,"_",typeAt,"_",le,"_",cut,"_class.pdf")) 
	image(1:n,1:n,matClass,col=c("white",rainbow(max(matClass))),xlab="",ylab="")
	invisible(dev.off())
}



# identify ranges of 1s in a binary matrix
# (initially we have 1s and they are converted to 0s)
# this function is a bit complicated and probably not efficient
# the idea was just to get some convenient ranges of indices to give to NOLB
writeZeroRanges<-function(mat,k,namesPDB,fname){
	#print(c(sum(mat==1),sum(mat==0)))
	#print(which(mat==1,arr=TRUE))
	if(sum(mat==1)>0){
		# convert 1s into 0s and get the zero stretches
		# as a list of intervals
		indJ = getZeroRanges(1-mat,k)
		#print(indJ)
		# values are columns, names are lines
		indI = names(indJ)
		res = c()
		n = length(indI)
		# for each element in the list (interval in i)
		for(i in 1:n){
			# for each interval in j matching the interval in i
			for(j in 1:length(indJ[[i]])){
				# get the start and end coordinates
				startendI = as.numeric(strsplit(indI[i],"-")[[1]])
				startendJ = as.numeric(strsplit(indJ[[i]][j],"-")[[1]])

				x1 = namesPDB[startendI[1]]
				x2 = substr(namesPDB[startendI[2]],2,nchar(namesPDB[startendI[2]]))
				y1 = namesPDB[startendJ[1]]
				y2 = substr(namesPDB[startendJ[2]],2,nchar(namesPDB[startendJ[2]]))
				res = rbind(res,c(paste(x1,x2,sep="-"),paste(y1,y2,sep="-")))
			}
		}
		write.table(res,paste(fname,"excl",sep="."),quote=FALSE,row.names=FALSE,col.names = FALSE)
	}
}

# identify contacts to be removed from two colored matrices
# the matrices were built using two different thresholds for defining a contact
# the contacts to be removed are detected on both matrices and we make their union
# the specification of the contacts to be removed is written out in a NOLB-friendly format
# target is the PDB file (code + chain combination), for example 1ex6A
# typeAt is the atom type, for example CA
# cutoff is an arbitrary chosen value to select the smallest patches 
extractExcludingListFromColoredMapCombi<-function(target,typeAt,k,size=625,cutoffs){
	# get the chain ids and resids from the PDB file
	namesPDB = getInfo(target)
	lmatClass=list()
	N = length(cutoffs)
	# read colored matrices (2 thresholds are considered for defining a contact)
	for(i in 1:N){
		lmatClass[[i]] = as.matrix(read.table(paste0(target,"_",typeAt,"_",k,"_",cutoffs[i],".class")))}
	#lmatClass[[3]] = as.matrix(read.table(paste0(target,"_",typeAt,"_2_0.2.class")))
	n = dim(lmatClass[[1]])[[1]]
	newMat = matrix(0,nc=n,nr=n)
	newMats = list()
	# identify the small patches, in both colored matrices
	for(i in 1:N){
		newMats[[i]] = matrix(0,nc=n,nr=n)
		t = table(lmatClass[[i]])
		val = as.numeric(names(t)[t<=size])
		p = length(val)
		# put 1s in the output matrix where the contact should be removed/excluded
		if(p>0){
			for(j in 1:p){
				newMat[lmatClass[[i]]==val[j]]=1
				newMats[[i]][lmatClass[[i]]==val[j]]=1
			}
		}
	}
	# write out the selection of contacts to be removed (NOLB friendly format)
	for(i in 1:N){
		# the 2 parameter determines how long the detected ranges are
		writeZeroRanges(newMats[[i]],2,namesPDB,paste0(target,"_",typeAt,"_",k,"_",cutoffs[i],".class"))
	}
	writeZeroRanges(newMat,2,namesPDB,paste(target,typeAt,"combiClass",sep="_"))
	return(newMat)
}

extractExcludingList<-function(target,alter,typeAt,cutoff){
	namesPDB = getInfo(target)
	diffMat = computeDiffNegMat(target,alter,typeAt)
	sel=which(diffMat>cutoff,arr=TRUE)
	tbrm = t(apply(sel,1,f<-function(x){return(namesPDB[x])}))
	print(dim(tbrm))
	write.table(tbrm,paste(target,"_",typeAt,".excl",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
}

