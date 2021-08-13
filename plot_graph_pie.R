args <- commandArgs(TRUE)

library(igraph)
set.seed(0)

edges <- read.table(args[1], header = FALSE, sep = "\t")
edges <- factor(as.vector(t(data.matrix(edges))))

names <- levels(edges)
edges <- as.numeric(edges)

table <- read.table(args[2], header = TRUE, sep = "\t")
rownames(table) <- table[, 1]
table <- table[names, -1]
rownames(table) <- 1:length(names)

g <- graph(edges)
x <- t(data.matrix(table))

pdf(file = args[3], width = as.numeric(args[4]), height = as.numeric(args[5]))
vertex.size <- as.numeric(args[6])
vertex.size <- rowSums(table) / max(rowSums(table)) * vertex.size
colors <- args[-(1:6)]
plot(g, vertex.shape = "pie", vertex.pie = split(x, rep(1:ncol(x), each = nrow(x))), vertex.pie.color = list(colors), vertex.size = vertex.size, vertex.label = names, vertex.label.dist = vertex.size ^ 0.5 / 2, vertex.label.degree = pi, layout = layout.reingold.tilford(g, root = 1))
legend("bottomright", legend = rownames(x), col = colors, pch=19)
dev.off()
