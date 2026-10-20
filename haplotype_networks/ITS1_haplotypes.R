library(tidyverse)
library(pegas)
library(Polychrome)
library(RColorBrewer)

ITS1 <- read.dna("final analyses/ITS1/ITS1_haps.fa", format = "fasta", as.character = FALSE)
new_labels <- labels(ITS1)
new_labels <- gsub("C\\.lav-S_", "C.lavaretus-Scotland_", new_labels)
new_labels <- gsub("C\\.lav-W_", "C.lavaretus-Wales_", new_labels)
new_labels <- gsub("C\\.lav-E_", "C.lavaretus-Europe_", new_labels)
new_labels <- gsub("C\\.lav-LD_", "C.lavaretus-LakeDistrict_", new_labels)
new_labels <- gsub("C\\.alb-LD_", "C.albula-LakeDistrict_", new_labels)
rownames(ITS1) <- new_labels
nrow(ITS1)
ITS1_char <- as.character(ITS1)
count_bases <- function(x) {
  sum(x != "-")
}
actual_lengths <- apply(ITS1_char, 1, count_bases)
summary(actual_lengths)
table(actual_lengths)
ITS1_clean <- ITS1[actual_lengths > 471, ]

labels_clean <- rownames(ITS1_clean)
split_list <- strsplit(labels_clean, "_")
split_matrix <- do.call(rbind, split_list)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)
ITS1_reg_cols <- c("C.lavaretus-LakeDistrict" = "green3", "C.lavaretus-Europe" = "tomato", "C.albula-LakeDistrict" = "steelblue2",
                   "C.lavaretus-Scotland"="goldenrod2", "C.lavaretus-Wales"="purple")

display.brewer.all(type = "qual")
pop_colours <- brewer.pal(n = length(unique(meta_df$Pop)), name = "Paired")
names(pop_colours) <- unique(meta_df$Pop)

meta <- data.frame(
  ID = split_matrix[, "ID"],
  Region = split_matrix[, "Region"],
  Pop = split_matrix[, "Pop"]
)
rownames(meta) <- labels_clean

ITS1_haps <- pegas::haplotype(ITS1_clean, prob = FALSE)
ITS1_haps
hapInfo <- stack(setNames(attr(ITS1_haps,"index"),rownames(ITS1_haps)))
names(hapInfo) <- c("index","haplotype")
head(hapInfo)
merged <- data.frame(cbind(hapInfo,meta[hapInfo$index,]))

net <- haploNet(ITS1_haps)
size <- attr(net, "freq")
adjusted_size <- sqrt(attr(net, "freq"))*2

ITS1_pie <- table(merged$haplotype,merged$Region)
ordered_cols <- ITS1_reg_cols[colnames(ITS1_pie)]

plot(net,size=adjusted_size, pie = ITS1_pie,bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL)
new_coords <- replot()
legend("bottomleft", colnames(ITS1_pie), fill = ordered_cols, title = "Regions")

ITS1_pop_pie <- table(merged$haplotype,merged$Pop)
ordered_cols <- pop_colours[colnames(ITS1_pop_pie)]

plot(net,size=adjusted_size, pie = ITS1_pop_pie, bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL)
new_coords <- replot()
legend("bottomleft", colnames(ITS1_pop_pie), fill = ordered_cols, title = "Populations")

hap_df <- as.data.frame(as.matrix(net))

pegas::hap.div(ITS1_clean, variance = TRUE)
pegas::nuc.div(ITS1_clean, pairwise.deletion = TRUE)

calb <- grepl("C.alb", labels_clean)
ITS1_calb <- ITS1[calb, ]
nrow(ITS1_calb)
hap.div(ITS1_calb, variance = TRUE)
nuc.div(ITS1_calb, pairwise.deletion = TRUE)

clav <- grepl("C.lav", labels_clean)
ITS1_clav <- ITS1[clav, ]
nrow(ITS1_clav)
hap.div(ITS1_clav, variance = TRUE)
nuc.div(ITS1_clav, pairwise.deletion = TRUE)

ITS1_mas <- read.dna("final analyses/ITS1/ITS1.fas", format = "fasta", as.character = FALSE)
labels_mas <- rownames(ITS1_mas)
split_list_mas <- strsplit(labels_mas, "_")
split_matrix <- do.call(rbind, split_list_mas)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)
all_regions <- unique(meta_df$Region)

calb <- grepl("C.alb", labels_mas)
ITS1_calb_mas <- ITS1_mas[calb, ]
clav <- grepl("C.lav", labels_mas)
ITS1_clav_mas <- ITS1_mas[clav, ]

dist_matrix <- as.matrix(dist.dna(ITS1_mas, model = "raw", pairwise.deletion = TRUE))
within_calb <- mean(dist_matrix[calb, calb][upper.tri(dist_matrix[calb, calb])])
within_clav <- mean(dist_matrix[clav, clav][upper.tri(dist_matrix[clav, clav])])
between_calb_clav <- mean(dist_matrix[calb, clav])
cat("Average p-distance within C.alb:", within_calb, "\n")
cat("Average p-distance within C.lav:", within_clav, "\n")
cat("Average p-distance between species:", between_calb_clav, "\n")

get_within_reg_dist <- function(reg_name, dist_mat, meta_info) {
  idx <- which(meta_info$Region == reg_name)
  if(length(idx) < 2) return(NA)
  sub_mat <- dist_mat[idx, idx]
  return(mean(sub_mat[upper.tri(sub_mat)]))
}
get_between_pop_stats <- function(pop1, pop2, dist_mat, meta_info) {
  idx1 <- which(meta_info$Region == pop1)
  idx2 <- which(meta_info$Region == pop2)

  sub_mat <- dist_mat[idx1, idx2]

  dist_mean <- mean(sub_mat)

  n <- length(sub_mat)
  dist_se <- sd(sub_mat) / sqrt(n)

  return(c(mean = dist_mean, se = dist_se))
}
within_distances <- sapply(all_regions, get_within_reg_dist, dist_mat = dist_matrix, meta_info = meta_df)
reg_combos <- combn(all_regions, 2)
pair_names <- apply(reg_combos, 2, paste, collapse = " vs ")
stats_list <- apply(reg_combos, 2, function(x) {
  get_between_pop_stats(x[1], x[2], dist_matrix, meta_df)
})
results_df <- as.data.frame(t(stats_list))
results_df$Comparison <- pair_names
results_df <- results_df[, c("Comparison", "Mean", "SE")]
print(results_df)
print(within_distances)

