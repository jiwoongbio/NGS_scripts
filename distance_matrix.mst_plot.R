args <- commandArgs(TRUE)

d <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(d) <- d[, 1]
d <- d[, -1]

library(igraph)
set.seed(0)

g <- graph.adjacency(data.matrix(d), mode = "undirected", weighted = TRUE)
g_mst <- mst(g)

if(args[2] != "") {
	table <- read.table(args[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
	colnames(table) <- c("sample", "name", "color")
	rownames(table) <- table$sample
	V(g_mst)$color <- table[V(g_mst)$name, "color"]
	V(g_mst)$label.color <- table[V(g_mst)$name, "color"]
	V(g_mst)$name <- table[V(g_mst)$name, "name"]
}

edge.color <- colorRampPalette(c("#000000", "#E5E5E5"))(1000)[cut(E(g_mst)$weight, 1000)]

pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
plot(g_mst, vertex.size = as.numeric(args[6]), edge.color = edge.color, layout = layout_with_lgl(g_mst, maxiter = as.numeric(args[7])))
dev.off()
