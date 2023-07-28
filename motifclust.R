library(peakToGene)
library(dirfns)
library(BSgenome.Crobusta.HT.KY)
library(CrobustaTFs)
library(tfenrichr)
library(moreComplexHeatmap)
library(TFBSTools)
library(DESeq2)

coords <- import('motifpos.bed')
htky <- getFeatures('HT.Gene.gff3')

overlaps <- lapply(htky,getOverlaps,coords)

featpeaks <- sapply(overlaps, function(x) unique(x[,1]))
genepeak <- do.call(rbind,overlaps)
getUpset(featpeaks,'motiffeat')

damotifs <- read.csv('DAmotifs.csv',row.names=1)
damotifs$PeakID <- row.names(damotifs)
damotifs$KYID <- sub('-.*','',row.names(damotifs))

motifdat <- read.csv('motifdat.csv',row.names="KYID",na.strings='')

dat <- merge(damotifs,motifdat[,1:3],by.x="KYID",by.y='row.names')
row.names(dat) <- dat$PeakID

bymotif <- split(row.names(damotifs),damotifs$KYID)
bymotif <- bymotif[row.names(motifdat)]

byfamily <- sapply(split(bymotif,motifdat$Family_Name),unlist)
names(bymotif) <- motifdat$MotifID


foxfup <- row.names(damotifs)[damotifs$FoxF_KO>1&damotifs$padj<0.1]
foxfdown <- row.names(damotifs)[damotifs$FoxF_KO< -1&damotifs$padj<0.1]

tvcacc <- row.names(damotifs)[(damotifs$mesp_dnFGFR< -1 | damotifs$mesp_MekMut>1)&damotifs$padj<0.1]
atmacc <- row.names(damotifs)[(damotifs$mesp_dnFGFR>  1 | damotifs$mesp_MekMut< -1)&damotifs$padj<0.1]

asmacc <- row.names(damotifs)[(damotifs$handr_dnFGFR< -1 | damotifs$handr_MekMut>1)&damotifs$padj<0.1]
heartacc <- row.names(damotifs)[(damotifs$handr_dnFGFR>  1 | damotifs$handr_MekMut< -1)&damotifs$padj<0.1]

peakset <- list(TVC_accessible=tvcacc,ATM_accessible=atmacc,ASM_accessible=asmacc,heart_accessible=heartacc,FoxF_KO_open=foxfup,FoxF_KO_closed=foxfdown)
getUpset(peakset,'peaksets')

gene <- read.delim("TVC514_gene_features.tsv") 
clusts <- split(gene$Gene.ID,gene$gene_network_leiden)
clustmotifs <- lapply(clusts,function(x) lapply(overlaps,function(y) y$PeakID[y$GeneID%in%x]))
motifbyclust <- lapply(clustmotifs, function(x) unique(unlist(x)))

getUpset(motifbyclust,'clustmotifs')

getHyper <- function(x,y,bg){
	x <- intersect(x,bg)
	y <- intersect(y,bg)
	q <- length(intersect(x,y))
	m <- length(y)
	n <- length(setdiff(bg,y))
	k <- length(x)
	p <- 1-phyper(q,m,n,k)
	odds <- log2((q/k)/(m/length(bg)))
	return(c(p=p,log2OR=odds,ct=q))}

getHyper2 <- function(x,y,bg){
	x <- intersect(x,bg)
	y <- intersect(y,bg)
	mat <- matrix(c(length(intersect(x,y)),
		      length(setdiff(y,x)), 
		      length(setdiff(x,y)),
		      length(setdiff(bg,union(x,y)))),2)
	res <- fisher.test(mat)
	return(c(p=res$p,log2OR=log2(res$estimate[[1]]),ct=mat[1]))
}

hyperPlot <- function(lx,ly,bg=row.names(damotifs),file,
		      rowsel=quote(apply(p < 0.10,1,any)),
		      colsel=quote(apply(p < 0.10,2,any))){
  if(!dir.exists(mkdate(file,''))){
	require(parallel)
	#cl <- makeForkCluster(getOption('cl.cores',detectCores()-2))
	hyper <- lapply(lx,
			function(x) {
				sapply(ly, 
					  function(y) getHyper2(x,y,bg))
			})
	#stopCluster(cl)
	hyper <- sapply(
		c('log2OR', 'p', 'ct'),
		function(x) t(sapply(hyper, '[', x,)),
		simplify=F)
	hyper <- setNames(hyper,c('log2OR','p','ct'))
	dir.apply(hyper,file)
} else {
		hyper <- sapply(lrtab(mkdate(file,''),
				      sep='\t',
				      na.string='',
				      check.names=F),
				as.matrix,simplify = F)
  }

	colsel <- with(hyper,eval(colsel))
	rowsel <- with(hyper,eval(rowsel))

	if(any(colsel)&any(rowsel)){
		dotPscale(hyper$log2OR[rowsel,colsel],
			  hyper$p[rowsel,colsel],
			  hyper$ct[rowsel,colsel],
			  outl.name='motifs',
			  file=file,
			  append.date=T,
			  size.name="-log10(p)")
	}
}

stringent <- quote(apply(p < 0.05 & log2OR > 1,1,any))
top50 <- quote(order(apply(log2OR,1,max),decreasing = T)[1:50])

clustPlots <- function(l,file,bg=row.names(damotifs),peaksets=T,
		       sel=stringent){
	if(peaksets){
		hyperPlot(peakset,l,bg,paste0('peaksetHyper',file),T,colsel=T)
	}
	hyperPlot(byfamily,l,bg,paste0('familyHyper',file),colsel=T)
	hyperPlot(bymotif,l,bg,paste0('motifHyper',file),sel,colsel=T)
}

promoter <- lapply(clustmotifs,'[[','promoter')
genebody <- lapply(lapply(clustmotifs,'[',
			  c('intron','CDS',
			    'five_prime_UTR',
			    'three_prime_UTR')),
		   function(x) unique(unlist(x)))
distal <- lapply(lapply(clustmotifs,'[',
			  c('upstream','downstream')),
		   function(x) unique(unlist(x)))
tvcclust <- lapply(motifbyclust,intersect,
		   peakset$TVC_accessible)
tvcpromoter <- mapply(intersect,promoter,tvcclust)
tvcbody <- mapply(intersect,genebody,tvcclust)
tvcdistal <- mapply(intersect,distal,tvcclust)

batchyper <- function(lx,ly,bg){
	require(parallel)
	cl <- makeForkCluster(detectCores()-2)
	hyper <- parLapply(cl,lx,
		function(x) {
			sapply(ly,
			       function(y) getHyper2(x,y,bg))
		})
	stopCluster(cl)
	hyper <- sapply(
		c('log2OR', 'p', 'ct'),
		function(x) t(sapply(hyper, '[', x,)),
		simplify=F)
	hyper <- setNames(hyper,c('log2OR','p','ct'))
	return(hyper)
}

f <- function(x,y,bg,file,fdr=F,split=NULL,cutoff=0.05,width=6,height=10){
	enriched <- batchyper(x,y,bg)
	dir.apply(enriched,file)
	if(fdr) {
		d <- dim(enriched$p)
		enriched$p <- matrix(p.adjust(enriched$p,'fdr'),d[1],d[2])
		size.name <- expression(log[10](FDR))
	}else{
		size.name <- expression(log[10](P))
	}
	sel <- apply(enriched$p<cutoff,1,any)

	mat.name <- expression(log[2](OR))
	mat <- enriched$log2OR[sel,]
	mat[mat==Inf] <- max(mat[is.finite(mat)])
	mat[mat==-Inf] <- min(mat[is.finite(mat)])

	dotPscale(mat,enriched$p[sel,],
		  enriched$ct[sel,],
		  size.name=size.name,mat.name=mat.name,
		  split=split[sel],
		  file=paste0(file,'_dot'),append.date=T,row_title_rot=0)

	hm <- Heatmap(mat,name='log2(OR)',
		  split=split[sel],
		  row_title_rot=0)

	dir.pdf(file,width=width,height=height)
	draw(hm)
	dev.off()
}

y <- mapply(union,promoter,genebody)

f(bymotif,y,unlist(y),'accmotifs',T,motifdat$Family_Name)
f(bymotif,y,intersect(unlist(y),unlist(peakset)),'DAmotifs',F,motifdat$Family_Name,0.01,height=12)
f(bymotif,tvcclust, union(unlist(tvcpromoter), unlist(tvcbody)),'TVCmotifs',split=motifdat$Family_Name)

f(byfamily,y,unlist(y),'accfamilies',F,height=6,width=4)
f(byfamily,y,intersect(unlist(y),unlist(peakset)),'DAfamilies',F,height=6,width=4)
f(byfamily,tvcclust, union(unlist(tvcpromoter), unlist(tvcbody)),'TVCfamilies',height=4,width=4)

da <- do.call(rbind,lapply(names(tvcclust),function(x) cbind(dat[tvcclust[[x]],c(1,7:15)],clust=x)))
mat <- da[,2:6]

dir.pdf('tvcmotifs',height=8,width=5)
Heatmap(mat,show_row_names = F,right_annotation = rowAnnotation(df=da[,"Family_Name",drop=F]),row_split=da$clust)
dev.off()

hyperPlot(byfamily,peakset,row.names(damotifs),
	  'DAfamilyHyper')
hyperPlot(bymotif,peakset,row.names(damotifs),
	  'DAmotifsHyper',top50)

clustPlots(motifbyclust,'')
clustPlots(promoter,'Promoter')
clustPlots(genebody,"Genebody")
clustPlots(distal,"Distal")

clustPlots(tvcclust,"TVC",peaksets=F,sel=top50)
clustPlots(tvcpromoter,"TVCpromoter",peaksets=F,sel=top50)
clustPlots(tvcbody,"TVCgenebody",peaksets=F,sel=top50)
clustPlots(tvcdistal,"TVCdistal",peaksets=F,sel=top50)

clustPlots(motifbyclust,'_vClust',unlist(motifbyclust))
clustPlots(promoter,'Promoter_vClust',unlist(motifbyclust))
clustPlots(genebody,"Genebody_vClust",unlist(motifbyclust))

clustPlots(tvcclust,"TVC_vTVC",tvcacc,F)
clustPlots(tvcpromoter,"TVCpromoter_vTVC",tvcacc,F)
clustPlots(tvcbody,"TVCgenebody_vTVC",tvcacc,F)

clustPlots(tvcclust,"TVC_vClust",unlist(tvcclust),F)
clustPlots(tvcpromoter,"TVCpromoter_vClust",unlist(tvcpromoter),F)
clustPlots(tvcbody,"TVCgenebody_vClust",unlist(tvcbody),F)
