args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

tip.color <- NULL
tip.color[colnames(table)] <- "#000000"
if(args[2] != "") {
	color.table <- read.table(args[2], header = FALSE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
	tip.color[color.table[, 1]] <- color.table[, 2]
}

pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
par(mar = c(0, 0, 0, 0))

table <- data.matrix(table)
library(ape)
plot(nj(as.dist(table)), font = 1, tip.color = tip.color)
dev.off()
