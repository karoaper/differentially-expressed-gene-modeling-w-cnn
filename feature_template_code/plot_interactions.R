plot_interactions <- function(n,regionlabels,regions,pairs)
{
		
	
	
	
	
}


draw_G1 <- function(x,y,size)
{
	scalingfactor = 1
	symbols(x,y,stars=matrix(rep(c(1,0.5),5),nrow=1),inches=scalingfactor*size,add=TRUE)
}


draw_T1 <- function(x,y,size,Er4)
{
	scalingfactor = 1
	symbols(x,y,squares=c(1),inches=scalingfactor*size,add=TRUE)
	if Er4
		text(x,y,"+")
	else
		text(x,y,"-")
}

draw_G2 <- function(x,y,size,Er4)
{
	scalingfactor = 1
	symbols(x,y,circles=c(1),inches=scalingfactor*size,add=TRUE)
	if Er4
		text(x,y,"+")
	else
		text(x,y,"-")
}


draw_H3K4me1 <- function(x,y,size,Er4)
{
	scalingfactor = 1
	symbols(x,y,stars=matrix(c(1,2/3,1,2/3),nrow=1),inches=scalingfactor*size,add=TRUE)
	if Er4
		text(x,y,"+")
	else
		text(x,y,"-")
}

draw_H3K4me3 <- function(x,y,size,Er4)
{
	scalingfactor = 1
	symbols(x,y,stars=matrix(c(1,4/3,1,0),nrow=1),inches=scalingfactor*size,add=TRUE)
	if Er4
		text(x,y,"+")
	else
		text(x,y,"-")
}

draw_H3K4me3 <- function(x,y,size,Er4)
{
	scalingfactor = 1
	symbols(x,y,stars=matrix(c(1,4/3,1,0),nrow=1),inches=scalingfactor*size,add=TRUE)
	if Er4
		text(x,y,"+")
	else
		text(x,y,"-")
}