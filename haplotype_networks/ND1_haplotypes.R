library(tidyverse)
library(pegas)
library(Polychrome)
library(RColorBrewer)

ND1 <- read.dna("final analyses/ND1/ND1_haps.fa", format = "fasta", as.character = FALSE)
new_labels <- labels(ND1)
new_labels <- gsub("C\\.lav-S_", "C.lavaretus-Scotland_", new_labels)
new_labels <- gsub("C\\.lav-W_", "C.lavaretus-Wales_", new_labels)
new_labels <- gsub("C\\.lav-E_", "C.lavaretus-Europe_", new_labels)
new_labels <- gsub("C\\.lav-LD_", "C.lavaretus-LakeDistrict_", new_labels)
new_labels <- gsub("C\\.alb-LD_", "C.albula-LakeDistrict_", new_labels)
new_labels <- gsub("C\\.alb-E_", "C.albula-Europe_", new_labels)
rownames(ND1) <- new_labels
ND1_char <- as.character(ND1)
count_bases <- function(x) {
  sum(x != "-")
}
actual_lengths <- apply(ND1_char, 1, count_bases)
summary(actual_lengths)
table(actual_lengths)
short_seq_indices <- which(actual_lengths < mean(actual_lengths))
short_sequences <- new_labels[short_seq_indices]
short_seq_df <- data.frame(
  Sequence = short_sequences,
  Length = actual_lengths[short_seq_indices]
)
print(short_seq_df)
ND1_clean <- ND1[actual_lengths > 621, ]

ND1_matrix <- as.matrix(ND1_clean)
lengths <- apply(ND1_matrix, 1, length)
unique(lengths)

labels_clean <- rownames(ND1)
split_list <- strsplit(labels_clean, "_")
split_matrix <- do.call(rbind, split_list)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)

ND1_reg_cols <- c("C.lavaretus-LakeDistrict" = "green3", "C.lavaretus-Europe" = "tomato", "C.albula-LakeDistrict" = "steelblue2",
                   "C.lavaretus-Scotland"="goldenrod2", "C.lavaretus-Wales"="purple", "C.albula-Europe"="darkred")

display.brewer.all(type = "qual")
pop_colours <- brewer.pal(n = length(unique(meta_df$Pop)), name = "Paired")
names(pop_colours) <- unique(meta_df$Pop)

meta <- data.frame(
  ID = split_matrix[, "ID"],
  Region = split_matrix[, "Region"],
  Pop = split_matrix[, "Pop"]
)
rownames(meta) <- labels_clean

ND1_haps <- pegas::haplotype(ND1_clean, prob = FALSE)
ND1_haps
hapInfo <- stack(setNames(attr(ND1_haps,"index"),rownames(ND1_haps)))
names(hapInfo) <- c("index","haplotype")
head(hapInfo)
merged <- data.frame(cbind(hapInfo,meta[hapInfo$index,]))

net <- haploNet(ND1_haps)
size <- attr(net, "freq")
adjusted_size <- sqrt(attr(net, "freq"))*2

ND1_pie <- table(merged$haplotype,merged$Region)
ordered_cols <- ND1_reg_cols[colnames(ND1_pie)]

plot(net,size=size, pie = ND1_pie, bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL)
new_coords <- replot()
plot.new()
legend(
  x = "topleft",
  legend = colnames(ND1_pie),
  fill = ordered_cols,
  title = "Regions",
  xpd = NA,
  inset = 0.02
)

ND1_pop_pie <- table(merged$haplotype,merged$Pop)
ordered_cols <- pop_colours[colnames(ND1_pop_pie)]

plot(net,size=size, pie = ND1_pop_pie, bg=ordered_cols, labels = FALSE,
     scale.ratio = 1, threshold = NULL, xy = new_coords)
new_coords <- replot()
legend("topleft", colnames(ND1_pop_pie), fill = ordered_cols, title = "Populations")

hap_df <- as.data.frame(as.matrix(net))
hap_df

pegas::hap.div(ND1_clean, variance = TRUE)
pegas::nuc.div(ND1_clean, pairwise.deletion = TRUE)

calb <- grepl("C.alb", labels_clean)
ND1_calb <- ND1[calb, ]
nrow(ND1_calb)
hap.div(ND1_calb, variance = TRUE)
nuc.div(ND1_calb, pairwise.deletion = TRUE)

clav <- grepl("C.lav", labels_clean)
ND1_clav <- ND1[clav, ]
nrow(ND1_clav)
hap.div(ND1_clav, variance = TRUE)
nuc.div(ND1_clav, pairwise.deletion = TRUE)

ND1_mas <- read.dna("final analyses/ND1/ND1-final-jan.fas", format = "fasta", as.character = FALSE)
labels_mas <- rownames(ND1_mas)
split_list_mas <- strsplit(labels_mas, "_")
split_matrix <- do.call(rbind, split_list_mas)
colnames(split_matrix) <- c("ID", "Region", "Pop")
meta_df <- as.data.frame(split_matrix)
all_regions <- unique(meta_df$Region)

calb <- grepl("C.alb", labels)
ND1_calb <- ND1[calb, ]
clav <- grepl("C.lav", labels)
ND1_clav <- ND1[clav, ]

dist_matrix <- as.matrix(dist.dna(ND1, model = "raw", pairwise.deletion = TRUE))
within_calb <- mean(dist_matrix[calb, calb][upper.tri(dist_matrix[calb, calb])])
within_clav <- mean(dist_matrix[clav, clav][upper.tri(dist_matrix[clav, clav])])
between_calb_clav <- mean(dist_matrix[calb, clav], na.rm = TRUE)
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

  if(length(idx1) == 0 || length(idx2) == 0) return(c(mean = NA, se = NA))

  sub_mat <- dist_mat[idx1, idx2]

  dist_mean <- mean(sub_mat, na.rm = TRUE)

  valid_n <- sum(!is.na(sub_mat))

  if(valid_n < 1) return(c(mean = NA, se = NA))

  dist_se <- sd(sub_mat, na.rm = TRUE) / sqrt(valid_n)

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

