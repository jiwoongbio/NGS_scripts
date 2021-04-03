args <- commandArgs(TRUE)

library(ggtree)
tree <- read.tree(args[1])

if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t")
	colnames(table) <- c("sample", "name", "color")
} else {
	table <- data.frame(sample = tree$tip.label, name = tree$tip.label, color = "black")
}
table$color <- factor(table$color)

p <- ggtree(tree) + theme_tree2()
p <- p %<+% table + geom_tiplab(aes(label = name, color = color), show.legend = FALSE) + scale_color_manual(values = levels(table$color))

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
