args <- commandArgs(TRUE)

d <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(d) <- d[, 1]
d <- d[, -1]

library(igraph)
set.seed(0)

g <- graph.adjacency(data.matrix(d), mode = "undirected", weighted = TRUE)
g_mst <- mst(g)

write.table(as_edgelist(g_mst, names = TRUE), file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
