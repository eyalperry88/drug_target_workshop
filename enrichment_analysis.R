
all_genes <- read.table("D:\\Projects\\drug_target_workshop\\results\\b2h_aggregated_results.txt", sep="\t", header=FALSE)

sorted_genes = as.character(all_genes[[1]])

causal_genes <- read.table("D:\\Projects\\drug_target_workshop\\data\\AML_census_genes.txt", sep="\t", header=FALSE)
causal_genes <- as.character(causal_genes[[1]])

all <- 1:length(sorted_genes)
causal_indices <- c()

for (i in 1:length(sorted_genes)) {
  if (sorted_genes[i] %in% causal_genes) {
    causal_indices <- c(causal_indices, i)
  }
}

drug_targets <- read.table("D:\\Projects\\drug_target_workshop\\data\\AML_drug_targets.txt", sep="\t", header=FALSE)
drug_targets <- as.character(drug_targets[[1]])

dd_indices <- c()

for (i in 1:length(sorted_genes)) {
  if (sorted_genes[i] %in% drug_targets) {
    dd_indices <- c(dd_indices, i)
  }
}
# One directional
png(filename="D:\\Projects\\drug_target_workshop\\results\\enrich.png",units="in", width=9.5, height=6, res=300)
barcodeplot(all, index = causal_indices, index2 = dd_indices, quantiles=c(length(sorted_genes) / 10, 10000000),
            labels = c("", "Best"), col.bars=c("burlywood4", "firebrick4"))
dev.off()
