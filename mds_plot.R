args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]
mds.out <- cmdscale(as.dist(table), k = 2, eig = TRUE, add = TRUE)

if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t")
	colnames(table) <- c("sample", "name", "color")
} else {
	table <- data.frame(sample = tree$tip.label, name = tree$tip.label, color = "black")
}
table$color <- factor(table$color)

table <- cbind(table, mds.out$points[table[, 1], ])
colnames(table) <- c("sample", "name", "color", "x", "y")

xlab <- paste("MDS1 (", round(mds.out$eig[1] / sum(mds.out$eig) * 100, digits = 1), "%)", sep = "")
ylab <- paste("MDS2 (", round(mds.out$eig[2] / sum(mds.out$eig) * 100, digits = 1), "%)", sep = "")

library(ggplot2)
p <- ggplot(data = table, aes(x = x, y = y, label = name, color = color)) + geom_text(check_overlap = TRUE) + xlab(xlab) + ylab(ylab)
ggsave(filename = args[3], plot = p, width = as.numeric(args[4]), height = as.numeric(args[5]))
