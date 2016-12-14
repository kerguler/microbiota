d<-read.csv("../data/Reads_OTU.csv",header=TRUE)

groups<-as.character(sapply(as.character(colnames(d[,c(-1,-2)])),function(a) {
    b<-strsplit(a,"[_.]")[[1]][1]
    if (substr(b,1,1)=="P") return("P") else return(b)
}))
groups<-factor(groups,levels=unique(groups),ordered=TRUE)

dd<-apply(d[,c(-1,-2)],2,function(a) a/sum(a))

dat<-data.frame(group=rep(groups,each=nrow(dd)),
                class=as.factor(rep(d[,1],length(groups))),
                freq=c(dd))

# Background information:
# http://www.coloss.org/beebook/I/gut-symbionts/2/2/5
#
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
# http://joey711.github.io/phyloseq-demo/unifrac.html
# https://joey711.github.io/phyloseq/import-data.html
# https://joey711.github.io/phyloseq/distance.html
library(phyloseq)

taxa<-strsplit(as.character(d[,1]),", ")
taxmat<-matrix(unlist(taxa),ncol=length(taxa[[1]]),byrow=TRUE)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:length(taxa[[1]])]

# Normalised sample abundances
OTU = otu_table(dd, taxa_are_rows = TRUE)
# Integer abundance values
OTU = otu_table(d[,c(-1,-2)], taxa_are_rows = TRUE)
# 
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
sampledata = sample_data(data.frame(
    Generation = groups,
    row.names = sample_names(physeq),
  stringsAsFactors=FALSE
))
physeq = merge_phyloseq(physeq,sampledata)

# These look better with normalised sample abundances
pdf("figures/barplot.pdf",width=24,height=6)
p<-plot_bar(physeq, fill="Phylum")
p <- p + scale_x_discrete(limits=sample_names(physeq))
p
dev.off()
# 
plot_heatmap(physeq, taxa.label="Phylum")
# 
# This requires integer abundance values
pdf("figures/richness.pdf",width=10,height=6)
plot_richness(physeq, x="Generation", color="Generation")
dev.off()

# This requires the phylogenetic distances (tree)
UniFrac(physeq, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)

# Principal Component Analysis
d<-read.csv("../data/Reads_OTU_Phylum.csv",header=TRUE)
dd<-apply(d[,c(-1,-2)],2,function(a) a/sum(a))

mat<-as.matrix(t(dd))
pr<-princomp(mat)
w1<-mat%*%pr$loadings[,1]
w2<-mat%*%pr$loadings[,2]

pdf("figures/PCA.pdf")
par(mar=c(5, 4, 0, 0)+0.1,mgp=c(3,1,0))
plot(w1,w2,t="p",pch=16,cex=2,col=rainbow(length(levels(groups)))[groups],xlab="PC1",ylab="PC2",frame=FALSE,asp=1)
text(w1,w2,as.character(colnames(dd)),cex=0.5)
dev.off()

###############################################################
