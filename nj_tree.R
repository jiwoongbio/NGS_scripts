args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

library(ape)
library(ggtree)
tree <- nj(as.dist(table))

if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
	colnames(table) <- c("sample", "name", "color")
} else {
	table <- data.frame(sample = tree$tip.label, name = tree$tip.label, color = "black")
}
scale_color_manual_values <- unique(table$color)
names(scale_color_manual_values) <- scale_color_manual_values
scale_color_manual_values[scale_color_manual_values == ""] <- NA

p <- ggtree(tree) + theme_tree2()
p <- p %<+% table + geom_tiplab(aes(label = name, color = color), show.legend = FALSE) + scale_color_manual(values = scale_color_manual_values)

if(length(args) == 5) {
	library(ggplot2)
	x.range <- layer_scales(p)$x$range$range
	p <- p + xlim(0, x.range[2] + (x.range[2] - x.range[1]) * 0.2)
} else {
	p <- p + xlim(as.numeric(args[6]), as.numeric(args[7]))
}
pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
plot(p)
dev.off()
