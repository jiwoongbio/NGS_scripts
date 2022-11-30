args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t")
colnames(table) <- c("x", "y")

if(any(grep("\\.png$", args[2]))) {
	png(filename = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]), type = "cairo")
} else {
	pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
}
plot(table, type = "l")
dev.off()
