library(driver)
library(rulesoflife)
library(ggplot2)
library(magrittr)
library(dplyr)
library(SpiecEasi)
library(Matrix)
library(propr)
library(ggraph)
library(igraph)
library(matrixsampling)
library(fido)
library(LaCroixColoR)
library(conflicted)
library(cowplot)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("alr", "driver")

source("COAT.R")

global_x_label <- "Proportion of joint zeros (overall samples)"

# ------------------------------------------------------------------------------
#   Plot effect of joint zeros in fido::basset
# ------------------------------------------------------------------------------

data <- load_data("ASV")

sname <- "ACA"

x <- data$counts[,data$metadata$sname == sname]
clr_x <- clr_array(x + 0.5, parts = 1)

z <- combn(1:(nrow(x)-1), m = 2)

joint0 <- sapply(1:ncol(z), function(i) {
  sum(x[z[1,i],] == 0 & x[z[2,i],] == 0)/length(x[z[1,i],])
})

either0 <- sapply(1:ncol(z), function(i) {
  sum(x[z[1,i],] == 0 | x[z[2,i],] == 0)/length(x[z[1,i],])
})

y <- readRDS(paste0("output/model_fits/asv_days90_diet25_scale1/MAP/", sname, ".rds"))
cormat <- cov2cor(y$Sigma[,,1])

clr_y <- fido::to_clr(y)

joint_cor <- sapply(1:ncol(z), function(i) {
  cormat[z[1,i],z[2,i]]
})

amboseli <- summarize_Sigmas("asv_days90_diet25_scale1")
scores <- apply(amboseli$rug, 2, function(x) calc_universality_score(x))
rho <- apply(amboseli$rug, 2, median)

plot_df_and <- data.frame(x = joint0,
                          y = joint_cor,
                          idx1 = amboseli$tax_idx1,
                          idx2 = amboseli$tax_idx2,
                          score = scores,
                          rho = rho)

plot_df_or <- data.frame(x = either0,
                         y = joint_cor,
                         idx1 = amboseli$tax_idx1,
                         idx2 = amboseli$tax_idx2,
                         score = scores,
                         rho = rho)

p1 <- ggplot(plot_df_and, aes(x = x, y = y, fill = score)) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(method = "loess", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 1) +
  labs(x = global_x_label, y = "ASV-ASV CLR correlation", fill = "Universality score") +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw()

# ggplot(plot_df_or, aes(x = x, y = y, fill = score)) +
#   geom_point(size = 2.5, shape = 21) +
#   geom_smooth(method = "loess", color = "black") +
#   labs(x = "frequency of 0-or-0", y = "ASV-ASV CLR correlation", fill = "Universality score") +
#   scale_fill_distiller(palette = "RdBu") +
#   theme_bw()

# ------------------------------------------------------------------------------
#   Plot effect of joint zeros when using proportionality
# ------------------------------------------------------------------------------

prop <- propr(t(x), metric = "rho")
prop_mat <- prop@matrix[1:125,1:125]
prop_rho <- prop_mat[lower.tri(prop_mat, diag = FALSE)]

prop <- propr(t(x), metric = "phi")
prop_mat <- prop@matrix[1:125,1:125]
prop_phi <- prop_mat[lower.tri(prop_mat, diag = FALSE)]

prop <- propr(t(x), metric = "phs")
prop_mat <- prop@matrix[1:125,1:125]
prop_phs <- prop_mat[lower.tri(prop_mat, diag = FALSE)]

plot_df_and2 <- data.frame(x = joint0,
                           y1 = prop_rho,
                           y2 = prop_phi,
                           y3 = prop_phs,
                           idx1 = amboseli$tax_idx1,
                           idx2 = amboseli$tax_idx2,
                           score = scores)

p2 <- ggplot(plot_df_and2, aes(x = x, y = y1, fill = score)) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(method = "loess", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 1) +
  labs(x = global_x_label, y = "ASV-ASV CLR proportionality (rho)", fill = "Universality score") +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw()

# ------------------------------------------------------------------------------
#   Plot effect of joint zeros when using COAT
# ------------------------------------------------------------------------------

coat_estimate <- coat(t(x + 0.5), soft = 1)$corr
coat_estimate <- coat_estimate[1:125,1:125]
coat_v <- coat_estimate[lower.tri(coat_estimate, diag = FALSE)]

plot_df_and2 <- plot_df_and
plot_df_and2$y <- coat_v
p3 <- ggplot(plot_df_and2, aes(x = x, y = y, fill = score)) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(method = "loess", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 1) +
  labs(x = global_x_label, y = "ASV-ASV CLR correlation", fill = "Universality score") +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw()

# ------------------------------------------------------------------------------
#   Plot effect of joint zeros when using SparCC
# ------------------------------------------------------------------------------

saved_scc <- "sparcc_output_saved.rds"
if(file.exists(saved_scc)) {
  scc <- readRDS(saved_scc)
} else {
  scc <- sparcc(t(x))
  saveRDS(scc, saved_scc)
}

sparcc_m <- scc$Cor[1:125,1:125]
sparcc_v <- sparcc_m[lower.tri(sparcc_m, diag = FALSE)]

plot_df_and <- data.frame(x = joint0,
                          y = sparcc_v,
                          idx1 = amboseli$tax_idx1,
                          idx2 = amboseli$tax_idx2,
                          score = scores)

p4 <- ggplot(plot_df_and, aes(x = x, y = y, fill = score)) +
  geom_point(size = 2.5, shape = 21) +
  geom_smooth(method = "loess", color = "black") +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 1) +
  labs(x = global_x_label, y = "ASV-ASV CLR correlation", fill = "Universality score") +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() +
  theme(legend.position = "bottom")

legend <- get_legend(p4)

p <- plot_grid(plotlist = list(p1 + labs(title = "fido::basset estimates") + theme(legend.position = "none"),
                               p2 + labs(title = "fido::basset estimates as proportionality (rho)") + theme(legend.position = "none"),
                               p3 + labs(title = "COAT estimates") + theme(legend.position = "none"),
                               p4 + labs(title = "SparCC estimates") + theme(legend.position = "none")), ncol = 2,
               labels = c("A", "B", "C", "D"))

ggsave(file.path("output", "figures", "Figure_1_Supplement_3.png"),
       plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.05)),
       units = "in",
       dpi = 200,
       height = 9,
       width = 9,
       bg = "white")

