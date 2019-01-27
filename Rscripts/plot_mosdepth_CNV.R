library(ggplot2)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
name=""
if (length(args)==0) {
  stop("At least one argument must be supplied for GenomePrefix", call.=FALSE)
} else {
	name=args[1]
}

manualColors = c("dodgerblue2","red1","grey20")
bedwindows = read.table(sprintf("coverage/%s.mosdepth.10000bp.gg.tab.gz",name),
                        header=F)
colnames(bedwindows) = c("CHR","Start","End","Depth","Strain","Accession")
chrsizes = aggregate(End ~ CHR, data=bedwindows,max)
chrlist = subset(chrsizes$CHR,chrsizes$End > 50000)
chrlist = factor(chrlist)
print(chrlist)
d=bedwindows[bedwindows$CHR %in% chrlist, ]
d <- d[order(d$CHR, d$Start), ]
d$index = rep.int(seq_along(unique(d$CHR)),
                  times = tapply(d$Start,d$CHR,length)) 

d$pos=NA

nchr = length(unique(chrlist))
lastbase=0
ticks = NULL
minor = vector(,8)

for (i in 1:nchr ) {
    if (i ==1) {
        d[d$index==i, ]$pos = d[d$index==i, ]$Start
    } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
        lastbase = lastbase + max(d[d$index==(i-1),"Start"])
	      minor[i] = lastbase
        d[d$index == i,"Start"] =
             d[d$index == i,"Start"]-min(d[d$index==i,"Start"]) +1
        d[d$index == i,"End"] = lastbase
        d[d$index == i, "pos"] = d[d$index == i,"Start"] + lastbase
    }
}
ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
ticks
minorB <- tapply(d$End,d$index,max,probs=0.5)
minorB
minor
xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)

pdffile=sprintf("plots/%s.CovDepth_mosdepth.pdf",name)
pdf(pdffile,width=7,height=2.5)
Title=sprintf("Depth of sequence coverage %s",name)

p <- ggplot(d,
            aes(x=pos,y=Depth,color=d$CHR)) +
	        geom_vline(mapping=NULL, xintercept=minorB,alpha=0.5,size=0.1,colour='grey15')	+
    geom_point(alpha=0.8,size=0.4,shape=16) +
    scale_color_brewer(palette="RdYlBu",type="seq") +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(name="Chromosome", expand = c(0, 0),
                       breaks = ticks,                      
                       labels=(unique(d$CHR))) +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() + 
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) 

p


plot_strain <- function(strain,data) {
 l = subset(data,data$Strain == strain)
 Title=sprintf("Chr coverage plot for %s",strain)
 p <- ggplot(l,
            aes(x=pos,y=Depth,color=CHR))  +
    scale_colour_brewer(palette = "Set3") +
    geom_point(alpha=0.9,size=0.8,shape=16) +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(name="Chromosome", expand = c(0, 0),
                       breaks=ticks,
                       labels=(unique(d$CHR))) +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1))
}


plts <- lapply(unique(d$Strain),plot_strain,data=d)
pdf(sprintf("plots/%s.StrainPlot.pdf",name))
plts

plot_chrs <-function (chrom, data) {
  Title=sprintf("Chr%s depth of coverage",chrom)
  l <- subset(data,data$CHR==chrom)
  l$bp <- l$Start
  p<-ggplot(l,
            aes(x=bp,y=Depth,color=Strain)) +
    geom_point(alpha=0.7,size=0.8,shape=16) +
    # scale_color_brewer(palette="RdYlBu",type="seq") +
    labs(title=Title,xlab="Position",y="Normalized Read Depth") +
    scale_x_continuous(expand = c(0, 0), name="Position") +
    scale_y_continuous(name="Normalized Read Depth", expand = c(0, 0),
                       limits = c(0,3)) + theme_classic() +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1))
}
pdf(sprintf("plots/%s.ChrPlot_10kb.pdf",name))
plts <- lapply(chrlist,plot_chrs,data=d)
plts

