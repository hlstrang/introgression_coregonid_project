library(ggtree)
library(tidyverse)
library(ape)
library(pals)
library(ggtext)
library(tidytree)
library(ggforce)

tree <- read.tree("final analyses/ND1/ND1_apr.nwk")
tree <- root(tree, outgroup = "C._artedi")

tipcategories = read.csv("final analyses/ND1/ND1_tipcats.csv",
                         sep = ",",
                         col.names = c("sample", "region"),
                         header = FALSE,
                         stringsAsFactors = FALSE)
dd = as.data.frame(tipcategories)

region_colours <- c(
  "C. albula - Lake District"     = "skyblue3",
  "C. albula - Europe"            = "darkred",
  "outgroup"                      = "black",
  "C. lavaretus - Lake District"  = "green3",
  "C. lavaretus - Scotland"       = "orange",
  "C. lavaretus - Llyn Tegid"     = "purple",
  "C. lavaretus - Europe"         = "red"
)

coloured_labels <- setNames(
  paste0("<span style='color:", region_colours, ";'>", names(region_colours), "</span>"),
  names(region_colours)
)

max_depth <- max(ggtree(tree)$data$x, na.rm = TRUE)

p <- ggtree(tree) +
  xlim(0, 22) +
  geom_treescale(fontsize = 3)

p <- p %<+% dd +
  geom_tiplab(aes(color = factor(region), x = x),
              size = 3.5,
              align = FALSE,
              offset = 0.2) +
  scale_color_manual(
    name = "Region",
    values = region_colours,
    labels = coloured_labels,
    guide = guide_legend(override.aes = list(label = ""))
  ) +
  geom_nodelab(aes(label = label, x = x),
               size = 3, hjust = 1.25, vjust = -0.6,
               data = td_filter(!isTip & !is.na(as.numeric(label)) & as.numeric(label) >= 0.1)) +
  theme(
    legend.position = "right",
    legend.text = element_markdown(size = 12),
    legend.key = element_blank(),
    legend.title = element_text(
      hjust = 0.5,
      margin = margin(b = 5)
    ))
p
ggsave("final analyses/ND1/ND1_ggtree.png", width = 15, height = 15)

ND1_mas <- read.dna("final analyses/ND1/ND1_apr_tree.fa", format = "fasta", as.character = FALSE)
labels_mas <- dd[-1, ]
labels_mas <- labels_mas[-82, ]
meta_df <- labels_mas %>%
  separate_wider_delim(
    cols = region,
    delim = " - ",
    names = c("species", "location"),
    too_few = "align_start"
  )
all_regions <- unique(meta_df$location)

calb <- grepl("C. albula", meta_df$species)
clav <- grepl("C. lavaretus", meta_df$species)

calb_mat <- dist_matrix[calb, calb]
clav_mat <- dist_matrix[clav, clav]

dist_matrix <- as.matrix(dist.dna(ND1_mas, model = "raw", pairwise.deletion = TRUE))
within_calb       <- mean(calb_mat[upper.tri(calb_mat)], na.rm = TRUE)
within_clav       <- mean(clav_mat[upper.tri(clav_mat)], na.rm = TRUE)
between_calb_clav <- mean(dist_matrix[calb, clav], na.rm = TRUE)
cat("Average p-distance within C.alb:", within_calb, "\n")
cat("Average p-distance within C.lav:", within_clav, "\n")
cat("Average p-distance between species:", between_calb_clav, "\n")

meta_df <- meta_df %>%
  mutate(species_pop = paste(species, location, sep = "_"))
all_populations <- unique(meta_df$species_pop)

get_within_reg_dist <- function(reg_name, dist_mat, meta_info) {
  idx <- which(meta_info$species_pop == reg_name)
  if(length(idx) < 2) return(NA)
  sub_mat <- dist_mat[idx, idx]
  return(mean(sub_mat[upper.tri(sub_mat)], na.rm = TRUE))
}
get_between_pop_stats <- function(pop1, pop2, dist_mat, meta_info) {
  idx1 <- which(meta_info$species_pop == pop1)
  idx2 <- which(meta_info$species_pop== pop2)

  sub_mat <- dist_mat[idx1, idx2]

  dist_mean <- mean(sub_mat, na.rm = TRUE)
  sub_mat_clean <- sub_mat[!is.na(sub_mat)]
  n <- length(sub_mat_clean)

  dist_se <- sd(sub_mat_clean) / sqrt(n)

  return(c(mean = dist_mean, se = dist_se))
}
pop_combos <- combn(all_populations, 2)
within_distances <- sapply(all_populations, get_within_reg_dist, dist_mat = dist_matrix, meta_info = meta_df)
stats_list <- apply(pop_combos, 2, function(x) {
  get_between_pop_stats(x[1], x[2], dist_matrix, meta_df)
})
stats_matrix <- t(stats_list)
results_df <- data.frame(
  Pop1 = pop_combos[1, ],
  Pop2 = pop_combos[2, ],
  Mean = stats_matrix[, "mean"],
  SE   = stats_matrix[, "se"],
  stringsAsFactors = FALSE
)

results_df <- results_df %>%
  separate_wider_delim(Pop1, delim = "_", names = c("Species_1", "Location_1")) %>%
  separate_wider_delim(Pop2, delim = "_", names = c("Species_2", "Location_2"))

print(results_df)
print(within_distances)

clav_pops_idx <- grepl("C. lavaretus", names(within_distances))
within_distances[clav_pops_idx]
mean(within_distances[clav_pops_idx], na.rm = TRUE)
