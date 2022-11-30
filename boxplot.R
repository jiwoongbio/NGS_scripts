args <- commandArgs(TRUE)

table <- read.table(args[1], header = FALSE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
colnames(table) <- c("x", "y")
table$x <- factor(table$x, levels = args[-(1:6)])

library(ggpubr)
library(rstatix)

t <- table %>% t_test(y ~ x)

ylim <- c(min(table$y), max(table$y))
ylim[2] <- ylim[1] + (ylim[2] - ylim[1]) / (1 - nrow(t) * 0.05)
t$y.position <- ylim[2] - (ylim[2] - ylim[1]) * (1:nrow(t) - 1) * 0.05

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
ggboxplot(table, x = "x", y = "y", fill = "x", xlab = args[5], ylab = args[6], ylim = ylim) + stat_pvalue_manual(t, label = "p.adj.signif", tip.length = 0.01) + theme(legend.position = "none")
dev.off()
