args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
rownames(table) <- table[, 1]
table <- table[, -1]
table <- table[apply(table, 1, function(x) {length(unique(x))}) > 1, ]
if("scale" %in% args[-(1:2)]) {
	table <- t(scale(t(table)))
}

color <- colorRampPalette(c("blue", "white", "red"))(1000)
for(index in grep("^color=", args[-(1:2)])) {
	color <- unlist(strsplit(sub("^color=", "", args[-(1:2)][index]), ","))
}

to <- max(table)
from <- min(table)
if(max(table) > 0 && min(table) < 0) {
	to <- max(abs(table))
	from <- -to
}
buffer <- (to - from) / (length(color) - 1) / 2
breaks <- seq(from = from - buffer, to = to + buffer, length.out = length(color) + 1)
for(index in grep("^breaks=", args[-(1:2)])) {
	breaks <- as.numeric(unlist(strsplit(sub("^breaks=", "", args[-(1:2)][index]), ",")))
}

border_color <- "white"
for(index in grep("^border_color=", args[-(1:2)])) {
	border_color <- sub("^border_color=", "", args[-(1:2)][index])
}

cluster_rows <- ifelse("no_cluster_rows" %in% args[-(1:2)] || nrow(table) < 2, FALSE, TRUE)
cluster_cols <- ifelse("no_cluster_cols" %in% args[-(1:2)] || ncol(table) < 2, FALSE, TRUE)

cellwidth <- 20
for(index in grep("^cellwidth=", args[-(1:2)])) {
	cellwidth <- as.numeric(sub("^cellwidth=", "", args[-(1:2)][index]))
}

cellheight <- ifelse(nrow(table) < 15, round(150 / nrow(table)), 10)
for(index in grep("^cellheight=", args[-(1:2)])) {
	cellheight <- as.numeric(sub("^cellheight=", "", args[-(1:2)][index]))
}
show_rownames <- TRUE
if("squash" %in% args[-(1:2)]) {
	cellheight <- 1
	show_rownames <- FALSE
}

pdf(file = paste(dirname(args[2]), "Rplots.pdf", sep = "/"))
library(pheatmap)
pheatmap(table, color = color, breaks = breaks, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, show_rownames = show_rownames, cluster_rows = cluster_rows, cluster_cols = cluster_cols, filename = args[2])
dev.off()
