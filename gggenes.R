args <- commandArgs(TRUE)

library(ggplot2)
library(gggenes)

input.table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)

if(!("genome" %in% colnames(input.table))) {
	if("sample" %in% colnames(input.table)) {
		input.table$genome <- input.table$sample
	} else {
		input.table$genome <- input.table$chromosome
	}
}
if(!("fill" %in% colnames(input.table))) {
	input.table$fill <- input.table$gene
}
if(!("label" %in% colnames(input.table))) {
	input.table$label <- input.table$gene
}
if(!("forward" %in% colnames(input.table))) {
	input.table$forward <- input.table$strand == "+"
}

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
ggplot(input.table, aes(xmin = start, xmax = end, y = genome, fill = fill, label = label, forward = forward)) +
	geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
	geom_gene_label(align = "left") +
	facet_wrap(~ genome, scales = "free_y", ncol = 1) +
	scale_fill_brewer(palette = "Set3") +
	theme_genes()
dev.off()
