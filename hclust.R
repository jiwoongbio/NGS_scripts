args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

table <- data.matrix(table)
if(identical(rownames(table), colnames(table)) && isSymmetric(table)) {
	hc <- hclust(as.dist(table))
} else {
	hc <- hclust(dist(t(table)))
}

tip.color <- NULL
tip.color[colnames(table)] <- "#000000"
if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
	colnames(table) <- c("sample", "name", "color")
	rownames(table) <- table$sample
	tip.color[hc$labels %in% table$sample] <- table[hc$labels[hc$labels %in% table$sample], "color"]
	hc$labels[hc$labels %in% table$sample] <- table[hc$labels[hc$labels %in% table$sample], "name"]
}

pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
par(mar = c(0, 0, 0, 0))

library(ape)
plot(as.phylo(hc), font = 1, tip.color = tip.color, underscore = TRUE)
dev.off()
