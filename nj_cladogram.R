args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

library(ape)
library(ggtree)
tree <- nj(as.dist(table))

if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t")
	colnames(table) <- c("sample", "name", "color")
} else {
	table <- data.frame(sample = tree$tip.label, name = tree$tip.label, color = "black")
}
table$color <- factor(table$color)

p <- ggtree(tree, layout = "circular", branch.length = "none")
p <- p %<+% table + geom_tiplab(aes(label = name, color = color), show.legend = FALSE) + scale_color_manual(values = levels(table$color))

library(ggplot2)
x.range <- layer_scales(p)$x$range$range

x.extend <- (x.range[2] - x.range[1]) * 0.2
p <- p + xlim(x.range[1] - x.extend, x.range[2] + x.extend)

pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
plot(p)
dev.off()
