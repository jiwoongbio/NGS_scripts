args <- commandArgs(TRUE)

table <- read.table(args[1], header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(table) <- table[, 1]
table <- table[, -1]

hc <- hclust(dist(table))
write.table(hc$labels[hc$order], file = args[2], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
