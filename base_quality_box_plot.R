library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

qualityFile = args[1]

waterMark = 1
if (length(args) == 2){
	waterMark = args[2]
}

a <- read.table(qualityFile, header=F, skip=1)

names(a) <- c("base", "mean", "median", "median1", "median3", "min", "max")

a$base = factor(a$base, level=a$base)

pdf(paste(qualityFile, ".pdf", sep=""), width=12, height=7)

rects <- data.frame(y1=c(-Inf,20,30), y2=c(20,30,Inf), col=letters[1:3])

p0 <- ggplot() + geom_rect(data=rects, mapping=aes(xmin=-Inf, xmax=Inf, ymin=y1, ymax=y2, fill=col), alpha = 0.3)

p1 <- p0 + geom_errorbar(data=a, mapping=aes(x=base, ymin=min, ymax=max), width = 0.5)

p2 <- p1 + geom_boxplot(data=a, stat='identity', fill="yellow", mapping=aes(x=base, lower=median1, middle=median, upper=median3, ymin=min, ymax=max)) + guides(fill=FALSE)

p3 <- p2 + scale_fill_brewer(palette="Set1") + scale_y_continuous(limits = c(0,43)) + scale_x_discrete() + geom_line(data=a, mapping=aes(x=base, y=mean), color='red', group=1) + theme_bw() + xlab("position in read (bp)") + ylab("quality score") + theme(axis.title.x = element_text(face="bold", size=15), axis.title.y = element_text(face="bold", size=15), axis.text.x = element_text(angle = 45, hjust = 1,size=11), axis.text.y = element_text(size=15)) + ggtitle(paste("Quality scores across all bases for ", strsplit(qualityFile, "/")[[1]][length(strsplit(qualityFile, "/")[[1]])-1])) + theme(plot.title = element_text(hjust = 0.5))

p_water <- annotate("text", x = Inf, y = -Inf, label = "Produced By Frasergen",  hjust=1.1, vjust=-1.1, col="white", cex=6, fontface = 'bold.italic', alpha = 0.8)

if (waterMark==0) {
	print(p3)
} else {
	print(p3 + p_water)
}

dev.off()
