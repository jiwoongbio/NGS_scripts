args <- commandArgs(TRUE)

d <- read.table(args[1], header = FALSE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
colnames(d) <- c("group", "time", "event")
if(length(args) > 4) {
	d$group <- factor(d$group, levels = args[-(1:4)])
} else {
	d$group <- factor(d$group)
}

library(survival)
library(survminer)
library(gridExtra)

fit <- survfit(Surv(time, event) ~ group, data = d, type = "kaplan-meier")
p <- ggsurvplot(fit, pval = TRUE, legend.title = "", legend.labs = levels(d$group))
t <- surv_median(fit)
t[, 1] <- sub("^group=", "", t[, 1])
colnames(t)[1] <- "group"
t$number <- table(d$group)[t[, 1]]
t <- t[, c("group", "number", colnames(t)[!(colnames(t) %in% c("group", "number"))])]

hr <- coxph(Surv(time, event) ~ group, data = d) %>% finalfit::fit2df()
hr <- hr$HR

pdf(file = args[2], width = as.numeric(args[3]), height = as.numeric(args[4]))
grid.arrange(p$plot, tableGrob(t, rows = NULL), tableGrob(paste0("HR = ", hr), rows = NULL), nrow = 3, heights = c(4, 1, 1))
dev.off()
