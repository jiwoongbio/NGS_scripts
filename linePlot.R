args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t")
colnames(table) <- c("x", "y")

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
plot(table, type = "l")
dev.off()
