args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)

table <- table[, c(1, 2, 3, as.numeric(args[2]))]
colnames(table) <- c("chromosome", "start", "end", "y")

chromosome <- args[3]
start <- as.numeric(args[4])
end <- as.numeric(args[5])

table <- table[table$chromosome == chromosome & start <= table$start & table$end <= end, ]

xlim <- c(start, end)
ylim <- c(as.numeric(args[6]), as.numeric(args[7]))

pdf(file = args[8], width = as.numeric(args[9]), height = as.numeric(args[10]))
plot(x = NULL, y = NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "")
segments(table$start, table$y, table$end, table$y)
dev.off()
