args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
colnames(table) <- c("x", "y")
if(length(args) > 8) {
	levels <- args[-(1:8)]
} else {
	levels <- unique(table$x)
}
table$x <- factor(table$x, levels = levels)

library(ggpubr)
library(rstatix)

t <- table %>% t_test(y ~ x)
t <- as.data.frame(t)

step.increase <- as.numeric(args[7])
if(step.increase == 0) {
	t <- t[sapply(1:(length(levels) / 2) * 2, function(x) {which(t$group1 == levels[x - 1] & t$group2 == levels[x])}), ]
}

gap <- as.numeric(args[8])
ylim <- c(0, max(sapply(levels, function(x) {mean_sd(table$y[table$x == x])$ymax})))
ylim[2] <- ylim[1] + (ylim[2] - ylim[1]) / (1 - (nrow(t) - 1) * step.increase - gap)
t$y.position <- ylim[2] - (ylim[2] - ylim[1]) * (1:nrow(t) - 1) * step.increase

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
if("p.adj.signif" %in% colnames(t)) {
	ggbarplot(table, x = "x", y = "y", fill = "x", xlab = args[5], ylab = args[6], ylim = ylim, add = "mean_sd") + stat_pvalue_manual(t, label = "p.adj.signif", tip.length = 0.01) + theme(legend.position = "none")
} else {
	ggbarplot(table, x = "x", y = "y", fill = "x", xlab = args[5], ylab = args[6], ylim = ylim, add = "mean_sd") + stat_pvalue_manual(t, label = "p", tip.length = 0.01) + theme(legend.position = "none")
}
dev.off()
