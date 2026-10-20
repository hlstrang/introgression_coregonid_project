library(tidyverse)
library(pegas)
library(Polychrome)
library(RColorBrewer)

GH2c <- read.dna("final analyses/GH2c/GH2c_haps.fa", format = "fasta", as.character = FALSE)
new_labels <- labels(GH2c)
new_labels <- gsub("C\\.lav-S_", "C.lavaretus-Scotland_", new_labels)
new_labels <- gsub("C\\.lav-W_", "C.lavaretus-Wales_", new_labels)
new_labels <- gsub("C\\.lav-E_", "C.lavaretus-Europe_", new_labels)
new_labels <- gsub("C\\.lav-LD_", "C.lavaretus-LakeDistrict_", new_labels)
new_labels <- gsub("C\\.alb-LD_", "C.albula-LakeDistrict_", new_labels)
rownames(GH2c) <- new_labels
GH2c_char <- as.character(GH2c)
count_bases <- function(x) {
  sum(x != "-")
}
actual_lengths <- apply(GH2c_char, 1, count_bases)
summary(actual_lengths)
table(actual_lengths)

split_list <- strsplit((new_labels), "_")
split_matrix <- do.call(rbind, split_list)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)
GH2c_reg_cols <- c("C.lavaretus-LakeDistrict" = "green3", "C.lavaretus-Europe" = "tomato", "C.albula-LakeDistrict" = "steelblue2",
                   "C.lavaretus-Scotland"="goldenrod2", "C.lavaretus-Wales"="purple")

display.brewer.all(type = "qual")
pop_colours <- brewer.pal(n = length(unique(meta_df$Pop)), name = "Paired")
names(pop_colours) <- unique(meta_df$Pop)

meta <- data.frame(
  ID = split_matrix[, "ID"],
  Region = split_matrix[, "Region"],
  Pop = split_matrix[, "Pop"]
)
rownames(meta) <- new_labels

GH2c_haps <- pegas::haplotype(GH2c, prob = FALSE)
GH2c_haps
hapInfo <- stack(setNames(attr(GH2c_haps,"index"),rownames(GH2c_haps)))
names(hapInfo) <- c("index","haplotype")
head(hapInfo)
merged <- data.frame(cbind(hapInfo,meta[hapInfo$index,]))

net <- haploNet(GH2c_haps)
size <- attr(net, "freq")
adjusted_size <- sqrt(attr(net, "freq"))*2

GH2c_pie <- table(merged$haplotype,merged$Region)
ordered_cols <- GH2c_reg_cols[colnames(GH2c_pie)]

plot(net,size=size, pie = GH2c_pie,bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL)
new_coords <- replot()
legend("topright", colnames(GH2c_pie), fill = ordered_cols, title = "Regions")

GH2c_pop_pie <- table(merged$haplotype,merged$Pop)
ordered_cols <- pop_colours[colnames(GH2c_pop_pie)]

plot(net,size=size, pie = GH2c_pop_pie, bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL, xy = new_coords)
new_coords <- replot()
legend("bottomleft", colnames(GH2c_pop_pie), fill = ordered_cols, title = "Populations")

hap_df <- as.data.frame(as.matrix(net))
hap_df

pegas::hap.div(GH2c, variance = TRUE)
pegas::nuc.div(GH2c, pairwise.deletion = TRUE)

calb <- grepl("C.alb", labels_clean)
GH2c_calb <- GH2c[calb, ]
nrow(GH2c_calb)
hap.div(GH2c_calb, variance = TRUE)
nuc.div(GH2c_calb, pairwise.deletion = TRUE)

clav <- grepl("C.lav", labels_clean)
GH2c_clav <- GH2c[clav, ]
nrow(GH2c_clav)
hap.div(GH2c_clav, variance = TRUE)
nuc.div(GH2c_clav, pairwise.deletion = TRUE)

GH2c_mas <- read.dna("final analyses/GH2c/GH2c.fas", format = "fasta", as.character = FALSE)
labels_mas <- rownames(GH2c_mas)
split_list_mas <- strsplit(labels_mas, "_")
split_matrix <- do.call(rbind, split_list_mas)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)
all_regions <- unique(meta_df$Region)

calb <- grepl("C.alb", labels_mas)
GH2c_calb_mas <- GH2c_mas[calb, ]
clav <- grepl("C.lav", labels_mas)
GH2c_clav_mas <- GH2c_mas[clav, ]

dist_matrix <- as.matrix(dist.dna(GH2c_mas, model = "raw", pairwise.deletion = TRUE))
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
stats_matrix <- t(stats_list)
results_df <- data.frame(
  Comparison = pair_names,
  Mean = stats_matrix[, "mean"],
  SE = stats_matrix[, "se"],
  stringsAsFactors = FALSE
)

print(results_df)
print(within_distances)

