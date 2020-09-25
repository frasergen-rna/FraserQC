library(ggplot2)
library(reshape2)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

contentFile = args[1]

waterMark = 1
if (length(args) == 2){
	waterMark = args[2]
}

a <- read.table(contentFile, header=F, skip=1)

names(a) <- c("posi", "GC")

a$GC <- a$GC/sum(a$GC)

pdf(paste(contentFile, ".pdf", sep=""), width=12, height=7)

p0 <- ggplot(data=a)

p1 <- p0 + geom_point(aes(x=posi, y=GC), color=brewer.pal(9, "Set1")[1])

p2 <- p1 + geom_smooth(aes(x=posi, y=GC), se=F, method="loess", span=0.2, color=brewer.pal(9, "Set1")[2])

p3 <- p2 + scale_fill_brewer(palette="Set1") + theme_bw() + xlab("GC content (%)") + ylab("density") + theme(axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + ggtitle(paste("sequence GC content for ", strsplit(contentFile, "/")[[1]][length(strsplit(contentFile, "/")[[1]])-1])) + theme(plot.title = element_text(hjust = 0.5))

p_water <- annotate("text", x = Inf, y = -Inf, label = "Produced By Frasergen",  hjust=1.1, vjust=-1.1, col="black", cex=6, fontface = 'bold.italic', alpha = 0.8)

if (waterMark==0) {
	print(p3)
} else {
	print(p3 + p_water)
}

dev.off()
