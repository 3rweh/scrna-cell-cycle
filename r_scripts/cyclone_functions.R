findChromosome <- function(geneInfoLine){
  chr = unique(as.array(strsplit(as.character(geneInfoLine[2]),";")[[1]]))
  if(length(chr) == 1){
    return (chr)
  }else{
    return (NA)
  }
}
findDir <- function(geneInfoLine){
  chr = unique(as.array(strsplit(as.character(geneInfoLine[5]),";")[[1]]))
  if(length(chr) == 1){
    return (chr)
  }else{
    return (NA)
  }
}


findStart <- function(geneInfoLine){
  starts = unique(as.array(strsplit(as.character(geneInfoLine[3]),";")[[1]]))
  start = min(starts)
  return (start)
}
findStop <- function(geneInfoLine){
  starts = unique(as.array(strsplit(as.character(geneInfoLine[4]),";")[[1]]))
  start = max(starts)
  return (start)
}


featureCount2genePred <- function(featureCount = exon){
  geneInfo = featureCount[,1:6]
  geneInfo$UCL <- apply(geneInfo, 1, findChromosome )
  geneInfo$Dir <- apply(geneInfo, 1, findDir )
  geneInfo$firstStart = apply(geneInfo, 1, findStart )
  geneInfo$lastStop = apply(geneInfo, 1, findStop )
  geneInfo$exonStart = paste(gsub(pattern = ";", replacement = ", ", geneInfo$Start),",", sep = "")
  geneInfo$exonStop = paste(gsub(pattern = ";", replacement = ", ", geneInfo$End),",", sep = "")
  geneInfo$nrOfExons = str_count(geneInfo$exonStop, ",")

  genePred = geneInfo[,c("Geneid", "UCL", "Dir", "firstStart","lastStop", "firstStart","lastStop", 
                         "nrOfExons", "exonStart", "exonStop")]
  colnames(genePred) = c("name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd",
                         "exonCount","exonStarts","exonEnds")
  return(genePred)
  
}

detect.genes <- function(x, cut=1) {
    length(which(x > cut))
}


corrdist <- function(x,method="pearson",transpose=T) {
	 if (transpose){
	    x<-t(x)	
	 }
         as.dist(1-cor(x,method=method))
}

corrdistS <- function(x,method="spearman",transpose=T) {
	 if (transpose){
	    x<-t(x)	
	 }	  
         as.dist(1-cor(t(x),method=method))
}

hclust.ward <-function(x) hclust(x, method="ward")


colplot.clust <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1,symbol="",...){
 ## modifiction of plclust for plotting hclust objects *in colour*!
 ## added option to put a symbol below each of the names
 ## Modified by Asa Bjorklund 2014
 ## from http://rafalab.jhsph.edu/batch/myplclust.R
 ## Copyright Eva KF Chan 2009
 ## Arguments:
 ##    hclust:    hclust object
 ##    lab:        a character vector of labels of the leaves of the tree
 ##    lab.col:    colour for the labels; NA=default device foreground colour
 ##    hang:     as in hclust & plclust
 ##    symbol:	 a vector of symbols to put below each leave
 ## Side effect:
 ##    A display of hierarchical cluster with coloured leaf labels.
 ## OBS! when changing cex, need to change hang to get labels at leaves

 if (length(symbol)>1) {
    lab<-paste(symbol,lab,sep="  ")
 }
 y <- rep(hclust$height,2)
 x <- as.numeric(hclust$merge)
 y <- y[which(x<0)]
 x <- x[which(x<0)]
 x <- abs(x)
 y <- y[order(x)]
 x <- x[order(x)]
 plot( hclust, labels=FALSE, hang=hang, ... )
 text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=lab[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
}


cv2.var.genes <- function(R,ERCC,plot=FALSE,cutP=0.1,main="CV2 vs mean",quant.fit=0.8,genes=NULL){
    # modified from Brenneke et al. by Asa Bjorlund 2014
    # input is: 
    # - R - matrix with RPKM/Counts for all genes (genes x samples)
    # - ERCC - matrix with RPKM/Counts for ERCCs  ( ercc x samples)
    # - plot=T/F if you want the plot or not
    # - cutP - P-value cutoff for selecting variable genes
    # - main - title of graph.
    # - quant.fit - cutoff for selecting spike-in genes to use for fit, 

    library(statmod)	      
    nSamp<-ncol(ERCC)
    
    # calculate for ERCC
    rem<-which(rowSums(ERCC)==0)
    if (length(rem)>0) { 
        ERCC<-ERCC[-rem,]
    }
    mE<-rowMeans(ERCC)
    vE<-apply(ERCC,1,var)
    cv2E<-vE/mE^2

   # for remaining genes.
   # remove all genes with sum<1
    remR<-which(rowSums(R)<1)
    if (length(remR)>0){
       R<-R[-remR,]
    }
    mR<-rowMeans(R)
    vR<-apply(R,1,var)
    cv2R<-vR/mR^2

    # creating the fit
    mForFit<- unname( quantile( mE[ which( cv2E > .3 ) ], quant.fit ) )
    useForFit <- mE >= mForFit
    fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/mE[useForFit] ), cv2E[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])

    # finding genes in R with variation over ercc
    minBiolDisp <- .5^2
    cv2th <- a0 + minBiolDisp + a0 * minBiolDisp
    psia1theta <- 1+a1
    testDenom <- ( mR * psia1theta + mR^2 * cv2th ) / ( 1 + cv2th/nSamp )
    
    p <- 1 - pchisq( vR * (nSamp-1) / testDenom, nSamp-1 )
    padj <- p.adjust( p, "BH" )
    sig <- padj < cutP
    sig[is.na(sig)] <- FALSE

    if (plot) {
         # Prepare the plot (scales, grid, labels, etc.)
       suppressWarnings(  plot( NULL, xaxt="n", yaxt="n",
             log="xy", xlim = range(mR), ylim = range(cv2R),
             xlab = "mean rpkm", ylab = "CV2" ))
        axis( 1, 10^(-4:5), c("0.0001", "0.001","0.01","0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
        axis( 2, 10^(-4:5), c("0.0001", "0.001","0.01","0.1", "1", "10", "100", "1000",expression(10^4), expression(10^5) ) )
        abline( h=10^(-3:4), v=10^(-5:5), col="#D0D0D0", lwd=2 )
        # Add the data points
        colS<-rep("red",nrow(R))
        colS[sig]<-"blue"
        points(mR,cv2R,pch=20,cex=0.2,col=colS)
        # add points for ERCC
        points(mE,cv2E,pch=20,cex=0.7,col="black")
        # Plot the fitted curve
        xg <- 10^seq( -2, 6, length.out=1000 )
        lines( xg, (a1)/xg + a0, col="#FF000080", lwd=3 )
        df <- nSamp - 1
        lines( xg, ( (a1)/xg + a0 ) * qchisq( .975, df ) / df,
              col="#FF000080", lwd=2, lty="dashed" )
        lines( xg, ( (a1)/xg + a0 ) * qchisq( .025, df ) / df,
              col="#FF000080", lwd=2, lty="dashed" )
        legend("topright",c("ERCC","non variable","variable"),col=c("black","red","blue"),pch=20,bty='n',cex=0.5)
	title(sprintf("%s\n%d var genes",main,length(which(sig))))
        if (!is.null(genes)){
           # remove genes that are not included in useForFit
           genes <- genes[genes %in% names(mR)]
           text(mR[genes],cv2R[genes],genes,cex=0.4)
        }

    }
    # convert the indices of sig to the original matrix before removing genes.
    idx.sig<-unname(sig)
    if (length(remR)>0){
      idx.sig<-rep(FALSE,nrow(R)+length(remR))
      idx.sig[-remR]<-sig
    }
    return(idx.sig)
}


calculate_overlap<-function(L,with.unique=TRUE, plot=FALSE){
   # takes a list with indices and creates a matrix with overlap between elements
   nL<-length(L)
   M<-mat.or.vec(nL,nL)
   nU<-mat.or.vec(nL,1)
   for (i in 1:nL){
      nU[i]<-length(setdiff(L[[i]],unlist(L[-i])))
      for (j in	 i:nL){
        M[i,j]<-length(intersect(L[[i]],L[[j]]))
      }
   }
   if (with.unique){
      M<-cbind(M,nU)
      colnames(M)<-c(names(L),"unique")
   }else {
      colnames(M)<-names(L)
   }
   rownames(M)<-names(L)

   if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=nL)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.5,mfrow=c(1,1)) 
      cexC = 0.2 + 0.5/log10(nL)
      h<-heatmap.2(M,cellnote=lab,scale="none",trace="none",density.info="none",notecex=0.7*cexC,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=F,cexRow=cexC,cexCol=cexC)
   }

   return(M)
}


overlap_phyper<-function(L,bg=length(unique(unlist(L))),with.unique=TRUE,plot=FALSE){
   # takes a list with indices and creates a matrix with overlap between elements
   # can also plot as a heatmap with coloring according to significance in overlap.
   # phyper test uses all entries in L as background if bg is not specified.

   nL<-length(L)
   M<-mat.or.vec(nL,nL)
   P<-mat.or.vec(nL,nL)
   P[,]<-1
   nU<-mat.or.vec(nL,1)
   for (i in 1:nL){
      nU[i]<-length(setdiff(L[[i]],unlist(L[-i])))
      for (j in  i:nL){
        M[i,j]<-length(intersect(L[[i]],L[[j]]))
	P[i,j]<-1-phyper(M[i,j],length(L[[i]]),bg-length(L[[i]]),length(L[[j]]))
	if (i==j) {P[i,j]<-NA}
      }
   }
   if (with.unique){
      M<-cbind(M,nU)
      colnames(M)<-c(names(L),"unique")
      P<-cbind(P,rep(1,length(nU)))
      colnames(P)<-c(names(L),"unique")
   }else {
      colnames(M)<-names(L)
   }

   rownames(M)<-names(L)
   rownames(P)<-names(L)
   if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=nL)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.5,mfrow=c(1,1))
      cexC = 0.2 + 0.2/log10(nL)
      notecex = 0.2 + 0.2/log10(nL)
      h<-heatmap.2(-log10(P+1e-16),cellnote=lab,scale="none",trace="none",density.info="none",notecex=notecex,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=T,cexRow=cexC,cexCol=cexC,key.xlab="-log10(p.value)",key.title='')
   }
   return(list(overlap=M,pval=P))
}

overlap_phyper2<-function(L1,L2,bg=length(unique(c(unlist(L1),unlist(L2)))),with.nTot=TRUE,plot=FALSE){
   # takes two list with indices and creates a matrix with overlap between elements
   # can also plot as a heatmap with coloring according to significance in overlap.
   # phyper test uses all entries in L as background if bg is not specified.

   nL1<-length(L1)
   nL2<-length(L2)
   M<-mat.or.vec(nL1,nL2)
   P<-mat.or.vec(nL1,nL2)
   P[,]<-1

   for (i in 1:nL1){
      for (j in  1:nL2){
        M[i,j]<-length(intersect(L1[[i]],L2[[j]]))
        P[i,j]<-1-phyper(M[i,j],length(L1[[i]]),bg-length(L1[[i]]),length(L2[[j]]))
      }
   }
   if (with.nTot){
      nT1<-unlist(lapply(L1,length))
      nT2<-unlist(lapply(L2,length))
      M<-cbind(M,nT1)
      colnames(M)<-c(names(L2),"Total")
      P<-cbind(P,rep(1,length(nT1)))
      colnames(P)<-c(names(L2),"Total")

      M<-rbind(M,c(nT2,bg))
      rownames(M)<-c(names(L1),"Total")
      P<-rbind(P,rep(1,length(nT1)+1))
      rownames(P)<-c(names(L1),"Total")
      nL1 = nL1+1
      nL2 = nL2+1

   }else {
      colnames(M)<-names(L2)
      rownames(M)<-names(L1)
   }

   if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=nL1)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.5,mfrow=c(1,1))
      cexC = 0.2 + 0.2/log10(nL1)
      notecex = 0.2 + 0.2/log10(nL1)
      h<-heatmap.2(-log10(P+1e-16),cellnote=lab,scale="none",trace="none",density.info="none",notecex=notecex,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=T,cexRow=cexC,cexCol=cexC,key.xlab="-log10(p.value)",key.title='')
   }
   return(list(overlap=M,pval=P))
}



plot.intersection <- function(v1,v2,plot=TRUE,...){
    # will plot a matrix with overlapping index positions for 2 vectors (defining group def in 2 different ways)	   
    t1<-table(v1) 	   
    t2<-table(v2)
    M<-mat.or.vec(length(t1)+1,length(t2)+1)
    for (i in 1:length(t1)){
    	c1<-which(v1==names(t1)[i])
    	for (j in 1:length(t2)){
	    c2<-which(v2==names(t2)[j])
	        M[i,j]<-length(intersect(c1,c2))
		}
    }
    rownames(M)<-c(names(t1),"total")
    colnames(M)<-c(names(t2),"total")
    M[,length(t2)+1]<-c(t1,0)
    M[length(t1)+1,]<-c(t2,0)
    # plot:
    if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=length(t1)+1)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.5,mfrow=c(1,1))
      maxlen<-max(length(t1),length(t2))
      cexC = 0.2 + 0.5/log10(maxlen)
      h<-heatmap.2(M,cellnote=lab,scale="none",trace="none",density.info="none",notecex=0.7,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=F,cexRow=cexC,cexCol=cexC,...)

    }
    return(M)
}


remove.uniquenames<-function(names){
        sub("\\.\\.\\d+$","",names)
}

change.duplicate.names<-function(names){
	# will give each gene one unique name, if some are duplicated
        t<-table(names)
        sel<-which(t>1)
        for (n in names(sel)){
            m<-which(names==n)
            for (i in 1:length(m)){
                names[m[i]] <- sprintf("%s..%d",names[m[i]],i)
            }
        }
        return(names)
}



plot.venn.ntot <- function(L,main="VENN"){
    # will make a venn diagram where total number in each group is displayed
    # input is a list containing the idfferent groups.
    library(gplots)
    v<-venn(L,show.plot=F)
    nTot<-length(unique(unlist(L)))
    nT<-unlist(lapply(L,length))
    names<-paste(names(L),as.character(nT),sep="\n")
    colnames(v)<-c("num",names)
    plot(v)
    title(sprintf("%s\n%d genes",main,nTot))
    return(v)
}

# make sets from a list of group memberships, can be strings or numbers
make_sets<-function(x){
	names<-names(table(x))
	l<-list()
	for (n in names){
	   l[[n]]<-which(x==n)    
	}
	return(l)
}



make_colors<-function(names,legend,colors){
	# make a color vector for each group based on sample names, a list of legend and colors 
	col<-mat.or.vec(1,length(names))
	for (i in 1:length(legend)){
	    col[names==legend[i]]<-colors[i]
	    }
	    return(col)
}
make.discrete.colors <- function(x,legend,col.def){
   # make a color vector for each group based on sample names, a list of legend and colors
   col<-mat.or.vec(1,length(x))
   for (i in 1:length(legend)){
       col[x==legend[i]]<-col.def[i]
   }
   return(list(cols=col,col.def=col.def[1:length(legend)],legend=legend))
}

 
convert.to.color <- function(x,colscale,col.range=NULL){
  x.range<-range(na.omit(x))
  by=0.1
  if (is.null(col.range)){
    by=10
    col.range<-seq(x.range[1],x.range[2],by=by)
  }
  col.def<-colscale(length(col.range))
  col.idx<-round((x-x.range[1])/by)+1
  col.idx[col.idx>length(col.range)]<-length(col.range)
  cols<-col.def[col.idx]
  return(list(cols=cols,col.def=col.def,col.range=col.range))
}

plot.color.bar<-function(col=col.log,crange=crange.log,main="log2(rpkm+1)"){
  barplot(rep(1,length(col)),col=col,border=NA,main=main,space=0,axes=F,cex.main=0.6)
  # round to closest 1,10,100,1000 etc.
  by.exp<-floor(log10(max(na.omit(crange))-min(na.omit(crange))))
  if (by.exp==1){
    by = 1
  }else {
    by<-10^(by.exp-1)
  }
  # if crange is not at even numbers.
  if (max(na.omit(crange)) > 100) {
    crange<-round(crange/10)*10
  }
  rem<-crange  %% by
  evens<-which(rem==0)
  if (length(evens) > 15) {
    rem<-crange  %% (by*5)
    evens<-which(rem==0)
  }
  lab<-crange[evens]
  axis(1,at=evens,labels=lab,cex.axis=0.5,las=2)
}

closest <- function(v,x){
  v[which.min(abs(v-x))]
}

newline.text <- function(x,len=35,split="[ _^;]"){
  n <- nchar(x)
  times <- floor(n/len)
  if (times==0){ return(x)}
    # find split places
    g <- gregexpr(split,x,perl=T)[[1]]
    for (t in 1:times){
        pos = closest(g,t*len)
        substr(x,pos,pos)<-"\n"
      }
    return(x)
}

invertList <- function(x) {
    xn <- unlist(unname(x), recursive=FALSE)
    xi <-rep(names(x),times=lengths(x))
    split(xi,xn)
}