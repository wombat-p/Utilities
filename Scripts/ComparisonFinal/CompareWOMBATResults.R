library(tidyr)
library(grid)
library(eulerr)
library(reshape2)
library(matrixStats)
library(yaml)
library(CGPfunctions)
library(ggplot2)
library(gridExtra)
library(stringr)
library(gplots)
library(RColorBrewer)
library(UpSetR)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ReactomePA)

# Set folder for data set
# folder <- "~/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD020394/dev/"
folder <- "~/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD009815/0.9.1/"
#folder <- "~/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD005507/0.9.2/"
#folder <- "~/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD009203/0.9.2/"

# Set workflow names
workflows <- c("compomics","maxquant","proline", "tpp")

# Function to load results from JSON files
loadWorkflowResults <- function(folder, workflows) {
  results <- list()
  allnames <- NULL
  
  for (w in workflows) {
    r <- unlist(read_yaml(file.path(folder, paste0("benchmarks_", w, ".json"))))
    if (r[11] == 0)
      r <- r[-11:-12]
    results[[w]] <- r
    names(results[[w]]) <- make.unique(names(r))
    allnames <- append(allnames, names(results[[w]]))
  }
  
  allnames <- unique(allnames)
  return(list(results = results, allnames = allnames))
}

# Load results
resultsData <- loadWorkflowResults(folder, workflows)
results <- resultsData$results
allnames <- resultsData$allnames

# Read the 'ups_prots.csv' file
ups <- try(read.csv("../ups_prots.csv")[, 2])

# Create the benchmark matrix
createBenchmarkMatrix <- function(results, allnames, workflows) {
benchmatrix <- matrix(NA, ncol = length(allnames), nrow = length(workflows), 
                      dimnames = list(rows = workflows, cols = allnames))
benchmatrix <- as.data.frame(benchmatrix)
for (w in workflows) {
  r <- results[[w]]
  benchmatrix[w, names(r)] <- r
}
return(benchmatrix)
}

# Load and create benchmark matrix
benchmatrix <- createBenchmarkMatrix(results, allnames, workflows)

# Calculate the number of peptides per proteins
peps_per_prot <- benchmatrix[,grep("Functionality.Performance.Identification.PeptidesPerProtein.Freq", colnames(benchmatrix))]
for (i in 1:ncol(peps_per_prot))
peps_per_prot[,i] <- as.numeric(peps_per_prot[,i])
barplot(as.matrix(peps_per_prot), names.arg = 1:ncol(peps_per_prot), legend.text = workflows, 
      main="Peptides per protein distribution", beside=T, border=0, xlab="Number of peptides")

# colors
c_workflows <- brewer.pal(4, "Set1")
c_types <- c("#20a3ff", "#c44601", "#0035f0", "#f57600")

# Define a function to create lollipop charts
createLollipopChart <- function(dat, title, subtitle) {
dat$Metric <- gsub("[A-Z,a-z]*\\.","", dat$Metric)
print(dat)
g <- ggplot(dat, aes(x = Workflow, y = value, color = Metric)) +
  geom_segment(aes(xend = Workflow, yend = 0), size = 1, position = position_dodge(width = 0.5)) +
  geom_point(size = 5,  position = position_dodge(width = 0.5)) +
  geom_text(aes(label = round(value, max(0, round(log(1/min(value)))))), size = 3, vjust = -0.5, hjust=1, position = position_dodge(width = 0.9)) +
  labs(title = title, subtitle = subtitle, x = "") +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c_types) +
  theme_light() +
  theme(axis.text.x = element_text(colour = c_workflows, size = 12, hjust = 1, vjust = 0.5, angle = 90, face = "bold"))

return(g)
}

# Metrics for lollipop charts
metrics_list <- list(
list(sel = c("Functionality.Performance.Identification.PeptideNumber", "Functionality.Performance.Identification.ProteinNumber"),
     title = "A) Quantification depth", subtitle = "Fraction quantified in all samples"),
list(sel = c("Functionality.Performance.Identification.ProteinCoverage", "Functionality.Performance.Identification.PeptideCoverage"),
     title = "B) Coverage", subtitle = "Fraction quantified in all samples"),
list(sel = c("Functionality.Performance.Quantification.CVPeptides", "Functionality.Performance.Quantification.CVProteins", 
             "Functionality.Performance.Quantification.CorrelationPeptides", "Functionality.Performance.Quantification.CorrelationProteins"),
     title = "C) Variance and similarity", subtitle = "Pearson correlation and coefficient of variance"),
list(sel = c("Functionality.Performance.Quantification.DynamicPeptideRange", "Functionality.Performance.Quantification.DynamicProteinRange"),
     title = "D) Dynamic range", subtitle = "Difference between lowest and highest reported numbers")
)

# Create lollipop charts
prep_bmatrix <- melt(t(benchmatrix))
colnames(prep_bmatrix) <- c("Metric", "Workflow", "value")
prep_bmatrix$value <- as.numeric(prep_bmatrix$value)
lollipop_charts <- lapply(metrics_list, function(metric_info) {
sel <- metric_info$sel
dat <- prep_bmatrix[prep_bmatrix$Metric %in% sel, ]
createLollipopChart(dat, metric_info$title, metric_info$subtitle)
})

# Arrange lollipop charts in a grid
grid.arrange(grobs = lollipop_charts, ncol = 2)

# Quantative evaluation of proteins
proteins <- s_proteins <- NULL
for (w in workflows) {
p <- read.csv(paste0(folder, "/stand_prot_quant_merged",w,".csv"))
# remove contaminants
p <- p[!grepl("CON_",p$protein_group),]
# fix for compomics 
p <- p[!is.na(p$protein_group),]
pg <- p$protein_group
pg <- gsub("sp\\|","",pg)
pg <- gsub(",",";",pg)
pg <- gsub(" ","",pg)
pg <- gsub("\\|[A-Z,0-9]*_YEAST","",pg)
pg <- gsub("\\|[A-Z,0-9]*_HUMAN_UPS","",pg)
pg <- gsub("\\|[A-Z,0-9]*_HUMAN","",pg)
rownames(p) <- pg
p <- p[,2:ncol(p)]
colnames(p) <- paste(w, colnames(p))
s_proteins[[w]] <- p[, grep(" abundance_", colnames(p))]
if(length(proteins) > 0) {
  proteins <- merge(proteins, p, by.x=1, by.y=0, all=T)  
} else {
  proteins <- cbind(rownames(p), p)
}
# remove empty columns
proteins <- proteins[, colSums(!is.na(proteins)) > 0]
}

rownames(proteins) <- proteins[,1]
# Select columns containing 'abundance_'
qproteins <- proteins[, grep(" abundance_", colnames(proteins))]
qproteins <- as.matrix(qproteins)
# Filter out rows with 'REV_' and 'CON_' in row names
qproteins <- qproteins[!grepl("REV_|CON_", rownames(qproteins)),]

# Extract only the first differential abundance q-value column of each workflow
diff_proteins <- as.matrix(proteins[, grep("differential_abundance_qvalue", colnames(proteins))])
if (grepl("PXD009815", folder)) {
diff_proteins <- diff_proteins[, grep("10\\.amol",colnames(diff_proteins))]
diff_proteins <- diff_proteins[, grep("500\\.amol",colnames(diff_proteins))]
} else  {
ttt <- NULL
for (w in workflows) {
  ttt <- cbind(ttt, diff_proteins[,grep(w, colnames(diff_proteins))[1]])
}
diff_proteins <- ttt
}
colnames(diff_proteins) <- workflows

# Setting side colors for sample types and workflows
colnames(qproteins) <- gsub("\\.", "_", colnames(qproteins))
col_classes <- sub("_[0-9]*$","",colnames(qproteins))
sample_types <- sub(paste(paste0(workflows, " "), collapse="|"), "", col_classes)

# Remove unwanted strings from sample types (UPS data)
sample_types <- gsub("X|CT_mixture_QY_|_CN_UPS1_CV_Standards_Research_Group", "", sample_types)

# Conditional coloring of sample types
if(all(grepl("(amol)|(fmol)", sample_types))) {
n_sample_types <- paste0("abundance_", c(paste0(c(10, 50, 100, 250, 500), "_amol"),
                                         paste0(c(1, 5, 10, 25, 50), "_fmol")))
c_sample_types <- as.factor(sample_types)
n_sample_types <- unique(sort(sample_types))
} else {
c_sample_types <- as.factor(sample_types)
n_sample_types <- unique(sort(sample_types))
}

# Define sample type and workflow colors
c_sample_types <- colorpanel(length((n_sample_types)), "#115f9a", "#76c68f", "#d0f400")[c_sample_types]
workflow_types <- str_extract( col_classes, paste(workflows, collapse="|"))
c_workflow_types <- brewer.pal(4, "Set1")[as.factor(workflow_types)]

# Heatmap plot of presence/absence
heatmap.2(1-as.matrix(dist(t(!is.na(qproteins))/nrow(qproteins), method="manhattan")), scale="none", trace="none",
        col=colorpanel(50, "white", "black"), ColSideColors = c_sample_types, RowSideColors = c_workflow_types, 
        srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), 
        sepwidth=c(0,0), sepcolor = NA, key=T, keysize=1.5, labRow = "", labCol = "", xlab="Sample types", 
        ylab = "Workflows")
legend("topright", fill=colorpanel(length(unique(sample_types)), "#115f9a", "#76c68f", "#d0f400"),
     legend=n_sample_types, xpd=TRUE, cex=0.8, border=0)
legend("bottomleft", fill=brewer.pal(4, "Set1"), 
     legend=workflows, xpd=TRUE, cex = 0.8, border=0)

# Heatmap plot of correlation between quantitative values
heatmap.2(cor((qproteins), use="pairwise.complete.obs"), trace="none", scale="none",
        col=colorpanel(50, "red", "white", "blue"), ColSideColors = c_sample_types, RowSideColors = c_workflow_types, 
        srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), rowsep=0, colsep=0, 
        labRow = "", labCol = "", xlab="Sample types", 
        ylab = "Workflows")
legend("topright", fill=colorpanel(length(unique(sample_types)), "#115f9a", "#76c68f", "#d0f400"),
     legend=n_sample_types, xpd=TRUE, cex=0.8, border=0)
legend("bottomleft", fill=brewer.pal(4, "Set1"), 
     legend=workflows, xpd=TRUE, cex = 0.8, border=0)

### Now on statistical test results
# Replace NA values in diff_proteins with 1
diff_proteins[is.na(diff_proteins)] <- 1

# Filter diff_proteins for 5% FDR threshold
diff_proteins_5perc <- diff_proteins[rowSums(diff_proteins < 0.05) > 0,]

# Filter diff_proteins for 1% FDR threshold
diff_proteins_1perc <- diff_proteins[rowSums(diff_proteins < 0.01) > 0,]

# Extract workflow names from selected columns
workflow_names <- str_extract(colnames(diff_proteins_5perc), paste(workflows, collapse = "|"))

# Define colors for workflow types
c_workflow_types <- brewer.pal(4, "Set1")[as.factor(workflow_names)]

# Select relevant columns based on folder name
onediff <- diff_proteins[, ]
onediff <- onediff[, ]

#Calculate true FDR per and across workflows (only if ups proteins with ground thruth)
fdrs <- colSums(onediff[!grepl("ups", rownames(onediff)),] < 0.05) / colSums(onediff < 0.05)
workflow_fdr <- sum(rowMeans(onediff[!grepl("ups", rownames(onediff)),]) < 0.05) / sum(rowMeans(onediff) < 0.05)

# Calculate FDR for a range of threshold values
fdr_sum <- num_diff <- data.frame()
fdr_range <- 10^seq(-5, 0, 0.01)
for (i in fdr_range) {
fdr_sum <- rbind(fdr_sum, colSums(onediff[!grepl("ups", rownames(onediff)),] < i) / colSums(onediff < i))
num_diff <- rbind(num_diff, colSums(onediff < i))
}

# Create a plot for true FDR vs FDR thresholds (only for ground truth data)
plot(fdr_range, fdr_sum[,1], log="x", ylim=c(min(na.omit(fdr_sum)),1), col=brewer.pal(4, "Set1")[1], 
   type="l", lty=2, xlab="FDR threshold", ylab="True FDR")
lines(fdr_range, fdr_sum[,2], col=brewer.pal(4, "Set1")[2], type="l", lty=2)
lines(fdr_range, fdr_sum[,3], col=brewer.pal(4, "Set1")[3], type="l", lty=2)
lines(fdr_range, fdr_sum[,4], col=brewer.pal(4, "Set1")[4], type="l", lty=2)
legend("topleft", fill=brewer.pal(4, "Set1"), legend=workflows, border=NA)
legend(7e-6, 0.6, lty=c(1,2), legend=c("Proteins", "True FDR"))
par(new=T)
# Overlay the plot with number of proteins
plot(fdr_range, num_diff[,1], ylim=c(0,max(num_diff, na.rm=T)), type="l", 
   col=brewer.pal(4, "Set1")[1], axes=F, log="x", lwd=1, xlab="", ylab="")
axis(4, ylim=c(range(num_diff)))
mtext("Number of proteins", 4)
lines(fdr_range, num_diff[,2], type="l", col=brewer.pal(4, "Set1")[2], lwd=1)
lines(fdr_range, num_diff[,3], type="l", col=brewer.pal(4, "Set1")[3], lwd=1)
lines(fdr_range, num_diff[,4], type="l", col=brewer.pal(4, "Set1")[4], lwd=1)

# Pairs plot comparing FDR values
pairs(log10(diff_proteins), cex=.5, col="#00000099")

## UpSet plot comparing the performance
colnames(onediff) <- workflows
common_font_size <- 6
upset(as.data.frame(onediff<0.05)*1, sets = workflows, sets.x.label = "Differentially abundant proteins",
    sets.bar.color = c_workflows, main.bar.color = "#AA66AA", mainbar.y.label = "Proteins in subset")  
grid.text("Differentially abundant proteins, FDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=common_font_size))
pdf("Figure3DiffRegUPS_1.pdf", width=4, height=4)
upset(as.data.frame(onediff[!grepl("ups", rownames(onediff)), ]<0.05)*1, 
    sets = workflows, sets.x.label = "Total number\nper workflow",
    sets.bar.color = c_workflows, main.bar.color = "#AA6666", 
    mainbar.y.label = "Proteins in subset")  
grid.text("Differentially abundant yeast proteins\nFDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=common_font_size))
dev.off()
pdf("Figure3DiffRegUPS_2.pdf", width=4, height=4)
upset(as.data.frame(onediff[grepl("ups", rownames(onediff)), ]<0.05)*1, sets = workflows,
    sets.x.label = "Total number\nper workflow",
    sets.bar.color = c_workflows, main.bar.color = "#6666AA", 
    mainbar.y.label = "Proteins in subset")
grid.text("Differentially abundant UPS proteins\nFDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=common_font_size))
dev.off()
upset(as.data.frame(onediff[grepl("ups", rownames(onediff)), ]<0.05)*1, sets = workflows,
      sets.bar.color = c_workflows, main.bar.color = "#66AA66", mainbar.y.label = "Proteins in subset")
grid.text("Differentially abundant ups proteins, FDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=10))



## RUn clusterProfiler on differentially abundant proteins
# reduce to only unique proteins
reddiff <- onediff[!grepl(";|,", rownames(onediff)),]
# change to entrez gene ids
mapped_ids <- bitr(rownames(reddiff), "UNIPROT", "ENTREZID", OrgDb = org.Hs.eg.db)
reddiff <- reddiff[mapped_ids$UNIPROT, ]
rownames(reddiff) <- mapped_ids$ENTREZID
enriched <- list()
for (i in workflow_names) {
  enriched[[i]] <- rownames(reddiff[reddiff[,i]<0.05,])
}
enriched$ThreeOfFour <- rownames(reddiff[rowSums(reddiff<0.05) > 2, ])
cc <- compareCluster(enriched, fun = "enrichPathway", pvalueCutoff = 0.05)
dotplot(cc, showCategory = 15)  + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(), axis.text.y = element_text(size=9))

