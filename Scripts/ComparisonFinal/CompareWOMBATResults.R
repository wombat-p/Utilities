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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


workflows <- c("compomics","maxquant","proline", "tpp")

folder <- "~/devel/Bioinformatics/ELIXIR_EDAM/WOMBAT-P/WOMBAT-P_Processed/PXD020394/dev/"
results <- list()
allnames <- NULL

### pdf file?
# pdf(paste0(folder, "_plots.pdf"), width=10)

ups <- NULL
try(ups <- read.csv("../ups_prots.csv")[,2])

for (w in workflows) {
  r <- unlist(read_yaml(paste0(folder,"/benchmarks_",w,".json")))
  if (r[11] == 0)
    r <- r[-11:-12]
  results[[w]] <- r
  names(results[[w]]) <- make.unique(names(r))
  allnames <- append(allnames, names(results[[w]]))
  
}

allnames <- unique(allnames)

benchmatrix <- matrix(NA,ncol=length(allnames), nrow=4, dimnames=list(rows=workflows, cols=allnames))
benchmatrix <- as.data.frame(benchmatrix)
for (w in workflows) {
  r <- results[[w]]
  
  benchmatrix[w, names(r)] <- r
}


# barplot(as.numeric(benchmatrix$Functionality.Performance.Identification.PeptideNumber), names.arg = workflows, main="Number of peptides")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Identification.ProteinNumber), names.arg = workflows, main="Number of proteins")
# 
# barplot(as.numeric(benchmatrix$Functionality.Performance.Identification.ProteinGroupNumber), names.arg = workflows, main="Number of protein groups")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Identification.ProteinCoverage), names.arg = workflows, main="Coverage of protein identification across samples")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Identification.PeptideCoverage), names.arg = workflows, main="Coverage of peptide identification across samples")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Statistics.DifferentialRegulatedProteins5Perc), names.arg = workflows, main="Differentially regulated proteins 5%")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Statistics.DifferentialRegulatedProteins1Perc), names.arg = workflows, main="Differentially regulated proteins 1%")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Statistics.MissingProteinValues), names.arg = workflows, main="Missing protein values")

peps_per_prot <- benchmatrix[,grep("Functionality.Performance.Identification.PeptidesPerProtein.Freq", colnames(benchmatrix))]
for (i in 1:ncol(peps_per_prot))
  peps_per_prot[,i] <- as.numeric(peps_per_prot[,i])
barplot(as.matrix(peps_per_prot), names.arg = 1:ncol(peps_per_prot), legend.text = workflows, 
        main="Peptides per protein distribution", beside=T, border=0, xlab="Number of peptides")


# barplot(as.numeric(benchmatrix$Functionality.Performance.Quantification.CorrelationProteins), names.arg = workflows, main="Correlation between protein quant")
# barplot(as.numeric(benchmatrix$Functionality.Performance.Quantification.CorrelationPeptides), names.arg = workflows, main="Correlation between peptide quant")

## Create a lollipop chart
# colors
c_workflows <- brewer.pal(4, "Set1")
c_types <- c("#0073e6", "#c44601", "#2546f0", "#f57600")

sel <- c("Functionality.Performance.Identification.PeptideNumber", "Functionality.Performance.Identification.ProteinNumber") 
#       "Functionality.Performance.Identification.ProteinGroupNumber")
tttb <- melt(t(benchmatrix))
colnames(tttb) <- c("Metric", "Workflow", "value")
tttb$value <- as.numeric(tttb$value)
tttb <- tttb[!is.na(tttb$value),]
ttt <- tttb[tttb$Metric %in% sel, ]
ttt$Metric <- gsub("[A-Z,a-z]*\\.","", ttt$Metric)
g1 <- ggplot(ttt, aes(x = Workflow, y = value, color = Metric)) +
  # Add segments from the baseline to the value
  geom_segment(aes(xend = Workflow, yend = 0), size = 1,
               position = position_dodge(width = 0.5)) +
  # Add points at the value
  geom_point(size = 5,  position = position_dodge(width = 0.5)) +
  # Add labels for the value
  geom_text(aes(label = round(value, 0)), size=3, hjust = rep(c(1.5, -0.5), 4), 
            position = position_dodge(width = 0.5)) +
  # Add titles and caption
  labs(title = "A) Identification depth",
       subtitle = "Fraction identified in all samples",
       caption = "", x = "") +
  # Remove the legend title
  theme(legend.title = element_blank()) +  theme_light() + 
  # Specify the colors and labels for the color aesthetic
  scale_color_manual(values = c_types) +
  # Rotate the x labels 90 degrees
  theme(axis.text.x = element_text(colour=c_workflows, size = 12, hjust = 1, vjust = 0.2, angle = 90,  face = "bold"))
# g1 <- newggslopegraph(ttt, Workflow, value, Metric, RemoveMissing=T, Title = "A) Identification depth", 
#               DataLabelPadding = 0.2, WiderLabels = T,
#               SubTitle = "", Caption = "") + ylim(-1000, max(ttt$value)) + geom_hline(yintercept = c(0))

sel <- c("Functionality.Performance.Identification.ProteinCoverage", "Functionality.Performance.Identification.PeptideCoverage")
ttt <- tttb[tttb$Metric %in% sel, ]
ttt$Metric <- gsub("[A-Z,a-z]*\\.","", ttt$Metric)
g2 <- ggplot(ttt, aes(x = Workflow, y = value, color = Metric)) +
  # Add segments from the baseline to the value
  geom_segment(aes(xend = Workflow, yend = 0), size = 1,
               position = position_dodge(width = 0.5)) +
  # Add points at the value
  geom_point(size = 5,  position = position_dodge(width = 0.5)) +
  # Add labels for the value
  # geom_text(aes(label = round(value, 0)), size=3, hjust = rep(c(1.5, -0.5), 4), 
  #           position = position_dodge(width = 0.5)) +
  # Add titles and caption
  labs(title = "B) Coverage",
       subtitle = "Fraction identified in all samples",
       caption = "", x = "") +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Specify the colors and labels for the color aesthetic
  scale_color_manual(values = c_types) + theme_light() + 
  # Rotate the x labels 90 degrees
  theme(axis.text.x = element_text(colour=c_workflows, size = 12, hjust = 1, vjust = 0.2, angle = 90,  face = "bold"))


# g2 <- newggslopegraph(ttt, Workflow, value, Metric, RemoveMissing=T, Title = "B) Coverage", WiderLabels = T, 
#             SubTitle = "Fraction identified in all samples", Caption = "") + 
# ylim(0,1) + geom_hline(yintercept = c(0,1))

sel <- c("Functionality.Performance.Quantification.CVPeptides", "Functionality.Performance.Quantification.CVProteins", 
         "Functionality.Performance.Quantification.CorrelationPeptides", "Functionality.Performance.Quantification.CorrelationProteins")
ttt <- tttb[tttb$Metric %in% sel, ]
ttt$Metric <- gsub("[A-Z,a-z]*\\.","", ttt$Metric)
g3 <- ggplot(ttt, aes(x = Workflow, y = value, color = Metric)) +
  # Add segments from the baseline to the value
  geom_segment(aes(xend = Workflow, yend = 0), size = 1,
               position = position_dodge(width = 0.7)) +
  # Add points at the value
  geom_point(size = 5,  position = position_dodge(width = 0.7)) +
  # Add labels for the value
  # geom_text(aes(label = round(value, 0)), size=3, hjust = rep(c(1.5, -0.5), length(sel)*2), 
  #           position = position_dodge(width = 0.5)) +
  # Add titles and caption
  labs(title = "C) Variance and similarity",
       subtitle = "Pearson correlation and coefficient of variance",
       caption = "", x = "") +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Specify the colors and labels for the color aesthetic
  scale_color_manual(values = c_types) + theme_light() + 
  # Rotate the x labels 90 degrees
  theme(axis.text.x = element_text(colour=c_workflows, size = 12, hjust = 1, vjust = 0.2, angle = 90,  face = "bold"))

# g3 <- newggslopegraph(ttt, Workflow, value, Metric, RemoveMissing=T, Title = "C) Variance and similarity", WiderLabels = T, 
#                 SubTitle = "Pearsson correlation and coefficient of variance", Caption = "") +
# ylim(0,1) + geom_hline(yintercept = c(0,1))

sel <- c("Functionality.Performance.Quantification.DynamicPeptideRange", "Functionality.Performance.Quantification.DynamicProteinRange")
ttt <- tttb[tttb$Metric %in% sel, ]
ttt$Metric <- gsub("[A-Z,a-z]*\\.","", ttt$Metric)
g4 <- ggplot(ttt, aes(x = Workflow, y = value, color = Metric)) +
  # Add segments from the baseline to the value
  geom_segment(aes(xend = Workflow, yend = 0), size = 1,
               position = position_dodge(width = 0.5)) +
  # Add points at the value
  geom_point(size = 5,  position = position_dodge(width = 0.5)) +
  # Add labels for the value
  geom_text(aes(label = round(value, 0)), size=3, hjust = rep(c(1.7, -0.7), 4), 
            position = position_dodge(width = 0.5)) +
  # Add titles and caption
  labs(title = "D) Dynamic range",
       subtitle = "Difference between lowest and hights reported numbers",
       caption = "", x = "") +
  # Remove the legend title
  theme(legend.title = element_blank()) +
  # Specify the colors and labels for the color aesthetic
  scale_color_manual(values = c_types) + theme_light() + 
  # Rotate the x labels 90 degrees
  theme(axis.text.x = element_text(colour=c_workflows, size = 12, hjust = 1, vjust = 0.2, angle = 90,  face = "bold"))

# 
# g4 <- newggslopegraph(ttt, Workflow, value, Metric, RemoveMissing=T, Title = "D) Dynamic range", SubTitle = "", 
#                 Caption = "", WiderLabels = T) +
# geom_hline(yintercept = c(0))

grid.arrange(g1, g2, g3, g4, ncol=2)



proteins <- s_proteins <- NULL
for (w in workflows) {
  p <- read.csv(paste0(folder, "/stand_prot_quant_merged",w,".csv"))
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
}

rownames(proteins) <- proteins[,1]
qproteins <- proteins[, grep(" abundance_", colnames(proteins))]
qproteins <- as.matrix(qproteins)
qproteins <- qproteins[!grepl("REV_", rownames(qproteins)),]
qproteins <- qproteins[!grepl("CON_", rownames(qproteins)),]
qproteins <- t(t(qproteins) - colMedians(qproteins, na.rm=T))
# hist(rowSums(is.na(qproteins)), 50)
# image(as.matrix(qproteins))
# hist(colMeans(qproteins, na.rm=T))

diff_proteins <- as.matrix(proteins[, grep("differential_abundance_qvalue", colnames(proteins))])


# Setting side colors for sample types and workflows
colnames(qproteins) <- gsub("\\.", "_", colnames(qproteins))
col_classes <- sub("_[0-9]*$","",colnames(qproteins))
sample_types <- sub(paste(paste0(workflows, " "), collapse="|"), "", col_classes)
c_sample_types <- colorpanel(length(unique(sample_types)), "#115f9a", "#76c68f", "#d0f400")

# Only needed for PXD009815
sample_types <- gsub("X", "", sample_types)
sample_types <- gsub("CT_mixture_QY_", "", sample_types)
sample_types <- gsub("_CN_UPS1_CV_Standards_Research_Group", "", sample_types)
if(all(grepl("(amol)|(fmol)", sample_types))) {
  n_sample_types <- paste0("abundance_", c(paste0(c(10, 50, 100, 250, 500), "_amol"), 
                                           paste0(c(1, 5, 10, 25, 50), "_fmol")))
  names(c_sample_types) <- n_sample_types
  c_sample_types <- c_sample_types[sample_types]
} else {
  c_sample_types <- c_sample_types[as.factor(sample_types)]
  n_sample_types <- unique(sort(sample_types))
}



workflow_types <- str_extract( col_classes, paste(workflows, collapse="|"))
c_workflow_types <- brewer.pal(4, "Set1")[as.factor(workflow_types)]

heatmap.2(1-as.matrix(dist(t(!is.na(qproteins))/nrow(qproteins), method="manhattan")), scale="none", trace="none", 
          col=colorpanel(200, "white", "black"), ColSideColors = c_sample_types, RowSideColors = c_workflow_types, 
          srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), 
          sepwidth=c(0,0), sepcolor = NA, key=T, keysize=1.5, labRow = "", labCol = "", xlab="Sample types", 
          ylab = "Workflows")
legend("topright", fill=colorpanel(length(unique(sample_types)),"#115f9a", mid="#76c68f", "#d0f400"),
       legend=n_sample_types, xpd=TRUE, cex=0.8, border=0)
legend("bottomleft", fill=brewer.pal(4, "Set1"), 
       legend=workflows, xpd=TRUE, cex = 0.8, border=0)
heatmap.2(cor((qproteins), use="pairwise.complete.obs"), trace="none", scale="none",
          col=colorpanel(200, "red", "white", "blue"), ColSideColors = c_sample_types, RowSideColors = c_workflow_types, 
          srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), rowsep=0, colsep=0, 
          labRow = "", labCol = "", xlab="Sample types", 
          ylab = "Workflows")
legend("topright", fill=colorpanel(length(unique(sample_types)), "#115f9a", "#76c68f", "#d0f400"),
       legend=n_sample_types, xpd=TRUE, cex=0.8, border=0)
legend("bottomleft", fill=brewer.pal(4, "Set1"), 
       legend=workflows, xpd=TRUE, cex = 0.8, border=0)

# qproteins_ups <- qproteins[grep("ups", rownames(qproteins)), ]
# heatmap(is.na(as.matrix(qproteins_ups))*1)
# boxplot((qproteins_ups))

### Now on statistical test results
diff_proteins[is.na(diff_proteins)] <- 1
diff_proteins_5perc <- diff_proteins[rowSums(diff_proteins < 0.05) > 0,]
diff_proteins_1perc <- diff_proteins[rowSums(diff_proteins < 0.01) > 0,]
workflow_names <- str_extract( colnames(diff_proteins_5perc), paste(workflows, collapse="|"))
c_workflow_types <- brewer.pal(4, "Set1")[as.factor(workflow_names)]

if (grepl("PXD009815", folder)) {
  onediff <- diff_proteins[, grep("10\\.amol", colnames(diff_proteins_5perc))]
  onediff <- onediff[, grep("500\\.amol", colnames(onediff))]
} else {
  onediff <- diff_proteins[, ]
  onediff <- onediff[, ]
}
# FDR per workflow
fdrs <- colSums(onediff[!grepl("ups", rownames(onediff)),]<0.05) / colSums(onediff<0.05)
fdrs
# FDR combining workflows
sum(rowMeans(onediff[!grepl("ups", rownames(onediff)),]) < 0.05) / sum(rowMeans(onediff) < 0.05)
# now for continuous threshold
fdr_sum <- num_diff <- data.frame()
fdr_range <- 10^seq(-5,0,0.01)
for(i in fdr_range){
  fdr_sum <- rbind(fdr_sum, colSums(onediff[!grepl("ups", rownames(onediff)),]<i) / colSums(onediff<i))
  num_diff <- rbind(num_diff, colSums(onediff<i))
}
plot(fdr_range, fdr_sum[,1], log="x", ylim=c(min(na.omit(fdr_sum)),1), col=brewer.pal(4, "Set1")[1], 
     type="l", lty=2, xlab="FDR threshold", ylab="True FDR")
lines(fdr_range, fdr_sum[,2], col=brewer.pal(4, "Set1")[2], type="l", lty=2)
lines(fdr_range, fdr_sum[,3], col=brewer.pal(4, "Set1")[3], type="l", lty=2)
lines(fdr_range, fdr_sum[,4], col=brewer.pal(4, "Set1")[4], type="l", lty=2)
legend("topleft", fill=brewer.pal(4, "Set1"), legend=workflows, border=NA)
legend(7e-6, 0.6, lty=c(1,2), legend=c("Proteins", "True FDR"))
par(new=T)
plot(fdr_range, num_diff[,1], ylim=c(0,1000), type="l", 
     col=brewer.pal(4, "Set1")[1], axes=F, log="x", lwd=1, xlab="", ylab="")
axis(4, ylim=c(range(num_diff)))
mtext("Number of proteins", 4)
lines(fdr_range, num_diff[,2], type="l", col=brewer.pal(4, "Set1")[2], lwd=1)
lines(fdr_range, num_diff[,3], type="l", col=brewer.pal(4, "Set1")[3], lwd=1)
lines(fdr_range, num_diff[,4], type="l", col=brewer.pal(4, "Set1")[4], lwd=1)

heatmap.2((diff_proteins_5perc), scale="none", trace="none", 
          col=colorpanel(1000, "white", "black"), ColSideColors=c_workflow_types,
          RowSideColors = rainbow(2)[grepl("ups", rownames(diff_proteins_5perc))+1],
          srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), 
          sepwidth=c(0,0), sepcolor = NA)
legend(-5, 100, fill=brewer.pal(4, "Set1"), border=0, 
       legend=workflows, xpd=TRUE, cex = 0.7)

workflow_names <- str_extract( colnames(diff_proteins_1perc), paste(workflows, collapse="|"))
c_workflow_types <- brewer.pal(4, "Set1")[as.factor(workflow_names)]

heatmap.2((diff_proteins_1perc), scale="none", trace="none", 
          col=colorpanel(1000, "white", "black"), ColSideColors=c_workflow_types,
          RowSideColors = rainbow(2)[grepl("ups", rownames(diff_proteins_1perc))+1],
          srtCol = 45, srtRow = 45, cexRow=40/ncol(qproteins), cexCol=40/ncol(qproteins), 
          sepwidth=c(0,0), sepcolor = NA)
legend(-5, 100, fill=brewer.pal(4, "Set1"), border=0, 
       legend=workflows, xpd=TRUE, cex = 0.7)

## eulerr with automated annotation
# Input in the form of a named numeric vector
colnames(onediff) <- workflows
upset(as.data.frame(onediff<0.05)*1, sets = workflows, sets.x.label = "Differentially regulated proteins",
      sets.bar.color = c_workflows, main.bar.color = "#AA6666", mainbar.y.label = "Proteins in subset")  
grid.text("Differentially regulated proteins, FDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=12))
upset(as.data.frame(onediff[grepl("ups", rownames(onediff)), ]<0.05)*1, sets = workflows,
      sets.bar.color = c_workflows, main.bar.color = "#6666AA", mainbar.y.label = "Proteins in subset")
grid.text("Differentially regulated ups proteins, FDR < 0.05",x = 0.65, y=0.95, gp=gpar(fontsize=12))

# fit2 <- euler(onediff[grepl("ups", rownames(onediff)),] < 0.01, shape="circle", loss="square")
# plot(fit2, quantities=T, fill=brewer.pal(4, "Set1"))

## RUn clusterProfiler on differentially regulated proteins
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

## Not revelant
peptides <- s_peptides <- NULL
for (w in workflows) {
  p <- read.csv(paste0(folder, "/stand_pep_quant_merged",w,".csv"))
  p# fix for compomics 
  p <- p[!is.na(p$modified_peptide),]
  pg <- p$modified_peptide
  print(head(pg[grep("Oxid", pg)]))
  pg <- gsub("sp\\|","",pg)v
  pg <- gsub(",",";",pg)
  pg <- gsub(" ","",pg)
  pg <- gsub("\\|[A-Z,0-9]*_YEAST","",pg)
  pg <- gsub("\\|[A-Z,0-9]*_HUMAN_UPS","",pg)
  rownames(p) <- pg
  p <- p[,2:ncol(p)]
  colnames(p) <- paste(w, colnames(p))
  s_peptides[[w]] <- p[, grep(" abundance_", colnames(p))]
  if(length(peptides) > 0) {
    peptides <- merge(peptides, p, by.x=1, by.y=0, all=T)  
  } else {
    peptides <- cbind(rownames(p), p)
  }
}

rownames(peptides) <- peptides[,1]


dev.off()

