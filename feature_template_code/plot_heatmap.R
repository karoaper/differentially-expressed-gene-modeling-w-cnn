do_mds_useggobi <- function(mydata,filename,title,classes,k=2,radius=1)
{
	
	d <- dist(mydata) # euclidean distances between the rows
#	fit <- cmdscale(d,eig=TRUE, k=k) # k is the number of dim
	fit <- princomp(mydata)
	

	fit # view results

#	x = fit$points[,1:k]
	x = fit$scores[,1:k]
	mydata = cbind(mydata,x)
#	print(dim(x))
	
	x.posclass = x[classes==0,]
	x.negclass = x[classes==1,]
	
	classlabs = unique(classes)
	
	mydata.posclass = mydata[classes==0,]
	mydata.negclass = mydata[classes==1,]
	
	
	c1 <- qtclust(mydata.posclass, radius=radius)
	c2 <- qtclust(mydata.negclass, radius=radius)
#	c3 <- qtclust(mydata,radius=radius)

	maxclass.pos = max(predict(c1)[!is.na(predict(c1))])
	maxclass.neg = max(predict(c2)[!is.na(predict(c2))])

	clusters = c(predict(c1),predict(c2))
	clusters[is.na(clusters)] = max(maxclass.pos,maxclass.neg)
	print(clusters)
	print(c(classes+1))

	g <- ggobi(mydata)	

	glyph_colour(g[1]) <- clusters
	glyph_type(g[1]) <- c(classes+1)	
}



plot_with_ggobi <- function(mydata,collabels,rowlabels,classes)
{
		
	mydata = matrix(rnorm(1000),nrow=10)
#	print(collabels)
#	print(rowlabels)
#	print(classes)

	library(rggobi)
#	colnames(mydata) = collabels
#	rownames(mydata) = rowlabels
#	print(mydata)

	mydata = data.frame(t(mydata[1:5,]))
	g <- ggobi(mydata)

	readline(prompt = "Pause. Press <Enter> to continue...")	
#	glyph_colour(g[1]) <- classes
#	glyph_type(g[1]) <- c(classes+1)		
}





do_pcatransform <- function(mydata,filename,title,classes,k=2,radius=1,line1,line2)
{	
#	d <- dist(mydata) # euclidean distances between the rows
#	fit <- cmdscale(d,eig=TRUE, k=k) # k is the number of dim
	fit <- princomp(mydata)
	fit # view results

#	x = fit$points[,1:k]
	x = fit$scores[,1:k]
	
	mydata = cbind(mydata,x)
#	print(dim(x))
	
	x.posclass = x[classes==0,]
	x.negclass = x[classes==1,]
	
	classlabs = unique(classes)
	
	mydata.posclass = mydata[classes==0,]
	mydata.negclass = mydata[classes==1,]
	
	
	line1_transformed = line1 %*% fit$loadings
	line2_transformed = line2 %*% fit$loadings

	line1_transformed = line1_transformed[,1:k]
	line2_transformed = line2_transformed[,1:k]
	
	c1 <- qtclust(mydata.posclass, radius=radius)
	c2 <- qtclust(mydata.negclass, radius=radius)
#	c3 <- qtclust(mydata,radius=radius)

	maxclass.pos = max(predict(c1)[!is.na(predict(c1))])
	maxclass.neg = max(predict(c2)[!is.na(predict(c2))])

	clusters = c(predict(c1),predict(c2))
	clusters[is.na(clusters)] = max(maxclass.pos,maxclass.neg)
	print(clusters)
	print(c(classes+1))

	if(filename != "")
	{
		pdf(filename,width=13,height=10)
	}
	
	opar <- par(c("mfrow","mar"))
	par(mfrow=c(2,1))


	maxclass = maxclass.neg + maxclass.pos
	mycolors = colors()[seq(367,657)][sample(length(seq(367,657)))[seq(1,maxclass)]]
	
#	print(mycolors)
	
	if(k==2)
	{
		require(car)
		plot(x[,1],x[,2],xlab="Coordinate 1", ylab="Coordinate 2",main='mds with class', col=classes+1,pch=20)
		plot(c(x.posclass[,1],x.negclass[,1]),c(x.posclass[,2],x.negclass[,2]),xlab="Coordinate 1", ylab="Coordinate 2",main='mds with clusters', col=c(mycolors[predict(c1)],mycolors[predict(c2)+maxclass.pos]),cex=1,pch=20)
	}
	else if(k==3)
{
		require(scatterplot3d)
		scatterplot3d(x[,1],x[,2],x[,3],xlab="Coordinate 1", ylab="Coordinate 2",zlab="Coordinate 3",main=title,color=classes+1,pch=20)
		scatterplot3d(c(x.posclass[,1],x.negclass[,1]),c(x.posclass[,2],x.negclass[,2]),c(x.posclass[,3],x.negclass[,3]),xlab="Coordinate 1", ylab="Coordinate 2",zlab="Coordinate 3",main=title,color=c(mycolors[predict(c1)],mycolors[predict(c2)+maxclass.pos]),pch=20)
	}
	
	par(opar)	
	if(filename != "")
	{
		dev.off()
	}
	
}






do_hclust <- function(data,k,rowinds,colinds,filename)
{
	d <- dist(data, method = "euclidean") # distance matrix
	fit <- hclust(d, method="complete")
	groups <- cutree(fit, k=k) # cut tree into 5 clusters

	draw_heatmaps(data,rowinds,colinds,4000,filename)
	groups
}

mycm.colors <- function (n, alpha = 1,lowcolor=0.7,highcolor=1,topsat=1,usezero=TRUE) 
{
	if ((n <- as.integer((n[1L]))) > 0) 
	{
		even.n <- n%%2 == 0
		if(even.n)
		{
			n = n + 1
			even.n <- FALSE
		}
		k <- n%/%2
		l1 <- k + 1 - even.n
		l2 <- n - k + even.n
		c(if (l1 > 0) hsv(h = lowcolor, s = seq.int(topsat, ifelse(even.n, topsat/k, 0), length.out = l1), v = 1, alpha = alpha),
		  if (l2 > 1) hsv(h = highcolor, s = seq.int(0, topsat, length.out = l2)[-1L], v = 1, alpha = alpha))
	}
	else character(0L)
}




do_qtclust <- function(x,k=6)
{
	require(flexclust)
	cl1 <- cclust(x, k)
	inds = predict(cl1,x)
#	print(parameters(cl1))
#	print(inds)
	c(parameters(cl1),inds)
}


#function to draw a heatmap, computing the gradient in a way that 0 takes the white color even if there may be more values above it than below it
draw_heatmaps <- function(clustRows,clustCols,x,n,rowlabels,filename,collabels=NULL,highcolor=1)
{
	

	numrows = dim(x)[1]
	numcols = dim(x)[2]
	

#	print(numrows)
#	print(numcols)
	
	if(is.null(collabels))
	{
		collabels = seq(1,numcols)
	}
#	print(collabels)
#	print(rowlabels)
	if(numrows == 1)
	{
		x = rbind(x,matrix(0,1,numcols))
	}
	
	if(numcols == 1)
	{
		x = cbind(x,matrix(0,numrows,1))
	}
	
	library(gplots)
	vals = c(x)
	
	pdf(filename,width=13,height=10)
	onlypositivevals = !any(vals < 0,na.rm=TRUE)
	
	ratio = max(vals,na.rm=TRUE)/abs(min(vals,na.rm=TRUE))
#	print(ratio)
	if(onlypositivevals)
	{
		Cols = gray(n:0 / n)
	}
	else
	{
		Cols = c(mycm.colors(n,highcolor=highcolor)[1:(n/2)],mycm.colors(n*ratio,highcolor=highcolor)[(n*ratio/2+1):(n*ratio)])
	}

#	print(Cols)
	hv <- heatmap(x, Colv=ifelse(clustCols,TRUE,NA), Rowv=ifelse(clustRows,TRUE,NA),col=Cols, scale="none",margins=c(5,10),xlab = "Kernel pixels", ylab= "",labCol=collabels,labRow=rowlabels)
	
	
	
	
	dev.off()
}



#function to draw a heatmap, computing the gradient in a way that 0 takes the white color even if there may be more values above it than below it
draw_several_heatmaps <- function(clustRows,clustCols,x,n,rowlabels,filename,highcolor=1)
{
	
	
	numrows = dim(x)[1]
	numcols = dim(x)[2]
	
	
	if(numrows == 1)
	{
		x = rbind(x,matrix(0,1,numcols))
	}
	
	if(numcols == 1)
	{
		x = cbind(x,matrix(0,numrows,1))
	}
	
	library(gplots)
	vals = c(x)
	
	pdf(filename,width=13,height=10)
	onlypositivevals = !any(vals < 0,na.rm=TRUE)
	
	ratio = max(vals,na.rm=TRUE)/abs(min(vals,na.rm=TRUE))
	
	if(onlypositivevals)
	{
		Cols = gray(n:0 / n)
	}
	else
	{
		Cols = c(mycm.colors(n,highcolor=highcolor)[1:(n/2)],mycm.colors(n*ratio,highcolor=highcolor)[(n*ratio/2+1):(n*ratio)])
	}
	
	hv <- heatmap(x, Colv=ifelse(clustCols,TRUE,NA), Rowv=ifelse(clustRows,TRUE,NA),col=Cols, scale="none",margins=c(5,10),xlab = "Kernel pixels", ylab= "",labRow=rowlabels)
	
	
	
	
	dev.off()
}



#function to draw a heatmap, computing the gradient in a way that 0 takes the white color even if there may be more values above it than below it
draw_heatmaps_custompalette <- function(clustRows,clustCols,x,n,rowlabels,filename)
{
	numrows = dim(x)[1]
	numcols = dim(x)[2]
	
	
	library(gplots)
		
	
	mycolors = readLines("/Users/khovsep/Dropbox/eblearn/demos/gene_expression/src/palette.dat")
	pdf(filename,width=10,height=10)
	
	
	image(seq(1,numcols),seq(1,numrows),z=t(x), zlim=range(1,256*256*256),col=mycolors, xlab = "Kernel pixels", ylab= "",axes=FALSE)
	axis(BELOW<-1, at=1:numcols, labels=seq(1,numcols), cex.axis=0.7)
	axis(LEFT <-2, at=1:length(rowlabels), labels=rowlabels, las= HORIZONTAL<-1,cex.axis=0.7)	
	dev.off()
}



get_palette <- function()
{
	x = array(0,256*256*256)
	for (i in seq(length(x)))
	{
		hexstr = as.character(as.hexmode(i),upper.case=TRUE)
		x[i] = paste("#",paste("",array(0,6-nchar(hexstr)),sep="",collapse=""),hexstr,sep="")
	}
	x
}
