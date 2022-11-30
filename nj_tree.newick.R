args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

library(ape)

tree <- nj(as.dist(table))

write.tree(tree, file = args[2])
