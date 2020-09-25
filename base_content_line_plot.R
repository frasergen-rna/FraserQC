library(ggplot2)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)

contentFile = args[1]

waterMark = 1
if (length(args) == 2){
	waterMark = args[2]
}

a <- read.table(contentFile, header=F, skip=1)

names(a) <- c("posi", "G", "A", "T", "C")

b <- melt(a)

names(b) <- c("posi", "base", "content")

b$posi = factor(b$posi, level=a$posi)

b$base = factor(b$base, level=c("A", "T", "C", "G"))

pdf(paste(contentFile, ".pdf", sep=""), width=12, height=7)

p0 <- ggplot(data=b)

p1 <- p0 + geom_line(aes(x=posi, y=content, group=base, color=base), width=1.5)

p2 <- p1 + scale_color_brewer(palette="Set1") + scale_x_discrete() + theme_bw() + xlab("position in read(bp)") + ylab("base content (%)") + theme(axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15), axis.text.x = element_text(angle = 45, hjust = 1,size=11), axis.text.y = element_text(size=15))+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,10)) + ggtitle(paste("base content across all bases for ", strsplit(contentFile, "/")[[1]][length(strsplit(contentFile, "/")[[1]])-1])) + theme(plot.title = element_text(hjust = 0.5))

p_water <- annotate("text", x = Inf, y = -Inf, label = "Produced By Frasergen",  hjust=1.1, vjust=-1.1, col="black", cex=6, fontface = 'bold.italic', alpha = 0.8)

if (waterMark==0) {
	print(p2)
} else {
	print(p2 + p_water)
}

dev.off()
