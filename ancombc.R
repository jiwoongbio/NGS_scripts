# Usage: Rscript ancombc.R read_count.tsv output.tsv alpha threads comma_separated_control_samples comma_separated_test_samples

# https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html
# https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html

library(phyloseq)
library(ANCOMBC)

args <- commandArgs(TRUE)

reads.file <- args[1]
output.file <- args[2]
alpha <- as.numeric(args[3])
n_cl <- as.numeric(args[4])
samples.list <- lapply(args[-(1:4)], function(x) {if(file.exists(x)) {scan(x, what = "character", sep = ",")} else {unlist(strsplit(x, ","))}})

reads <- read.table(reads.file, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
rownames(reads) <- reads[, 1]
reads <- reads[, -1]

colnames <- colnames(reads)
colnames(reads) <- 1:length(colnames)

colnumbers <- c()
conditions <- c()
for(index in 1:length(samples.list)) {
	for(colnumber in colnames(reads)[colnames %in% samples.list[[index]]]) {
		colnumbers <- c(colnumbers, colnumber)
		conditions <- c(conditions, index)
	}
}
reads <- reads[, colnumbers]

rownames <- rownames(reads)
rownumbers <- 1:length(rownames)
rownames(reads) <- rownumbers

sample_data <- sample_data(data.frame(conditions))
rownames(sample_data) <- colnumbers
phyloseq <- phyloseq(otu_table(reads, taxa_are_rows = TRUE), sample_data)

ancombc.out <- ancombc2(data = phyloseq, fix_formula = "conditions", p_adj_method = "fdr", alpha = alpha, n_cl = n_cl)

output <- ancombc.out$res[, grepl("_conditions$", colnames(ancombc.out$res))]
colnames(output) <- sub("_conditions$", "", colnames(output))
rownames(output) <- ancombc.out$res$taxon
output <- output[rownumbers, ]
rownames(output) <- rownames

write.table(cbind(feature = rownames(output), output), file = output.file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
