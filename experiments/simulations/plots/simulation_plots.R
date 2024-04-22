library(tidyverse)
library(ggplot2)
library(pracma)
library(tidyverse)
theme_set(theme_bw(base_size = 14))

summ <- read_csv( "Fits/results_sparse_rrr_experiments.csv")

legend_order <- c("Oracle", "RRR-ADMM-opt", "FIT_SAR_CV",
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "Waaijenborg-CV", "Waaijenborg-Author", "SCCA_Parkhomenko")
labels_n <- c("Oracle",  
              "\nsparse-CCA-RRR\n(this paper)",
              "\nSAR\n(Wilms et al)",
              "\nSAR with BIC\n(Wilms et al)",
              "\nCCA-Lasso, permuted\n(Witten et al)",
              "\nCCA-Lasso\n(Witten et al)",
              "\nSSCA\n(Waaijenborg et al)",
              "\nSSCA with BIC\n(Waaijenborg et al)",
              "\nSCCA\n(Parkhomenko et al)")

cbPalette <- c("#000000",
               "#009E73",
               "#D55E00",
               "#CC79A7",
               "#0072B2",
               "#56B4E9",
               "#E69F00", 
               "#F0E442",
               "#999999")



summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))

ggplot(summ %>% filter(r_pca == 5, r==2,
                        nnzeros == 10, 
                        p1 <= 3000,
                        p2<= 50,
                        n == 500,
                        method %in% legend_order,
                        overlapping_amount==1),
aes(x=p1, y = distance_tot_mean, colour = method, fill = method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(p2~theta_strength, labeller = as_labeller(c(`10` = "q = 10",
                                                                         `30` = "q = 30",
                                                                         `50` = "q = 50",
                                                                         `80` = "q = 80",
                                                                         `100` = "n = 100",
                                                                         `200` = "n = 200",
                                                                         `300` = "n = 300",
                                                                         `500` = "n = 500",
                                                                         `high` = "High",
                                                                         `1000` = "n = 1,000",
                                                                         `2000` = "n = 2,000",
                                                                         `10000` = "n = 10,000",
                                                                         `medium` = "Medium",
                                                                         `low` = "Low"))) +
  xlab("Number of covariates in X (p)") + 
  ylab(expression("Subspace distance")) +
  ylim(0,2.8)+
  theme(legend.position = "top")+
  theme(legend.key.size = unit(0.1, 'cm'), legend.text = element_text(size=10), legend.title = element_blank())+
  guides(color=guide_legend(nrow=3))
ggsave("final_rrr_no_bars.pdf", width = 6.5, height = 7)

ggplot(summ %>% filter(r_pca == 5, r==2,
                        nnzeros==10, 
                        p2 == 30,
                        p1==500,
                        method %in% legend_order,
                        overlapping_amount==1),
aes(x=n, 
    y = distance_tot_q50, 
    colour = method,
    fill = method)) +
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
                    colour = NULL), alpha=0.2)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(p2 ~ theta_strength, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                          `30` = "p = 500, q = 30",
                                                                          `50` = "q = 50",
                                                                          `70` = "q = 70",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                          `500` = "n = 500",
                                                                          `high` = "High",
                                                                          `1000` = "n = 1,000",
                                                                          `2000` = "n = 2,000",
                                                                          `10000` = "n = 10,000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("Sample size (n)") + 
  ylab(expression("Subspace distance")) +
  labs(colour = "Method") + 
  scale_x_log10(breaks = c(100,1000,10000),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position = "none")
ggsave("final_sparse_rrr_function_of_n.pdf", width = 6.5, height = 2.4)


ggplot(summ %>% filter(r_pca == 5, n==500,
                       nnzeros==10, 
                       p2 == 30,
                       p1==500,
                       method %in% legend_order,
                       overlapping_amount==1),
       aes(x = r, 
           y = distance_tot_q50, 
           colour = method,
           fill = method)) +
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
                  colour = NULL), alpha=0.2)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(p2 ~ theta_strength, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
                                                                           `30` = "p = 500, q = 30",
                                                                           `50` = "q = 50",
                                                                           `70` = "q = 70",
                                                                           `100` = "n = 100",
                                                                           `200` = "n = 200",
                                                                           `300` = "n = 300",
                                                                           `500` = "n = 500",
                                                                           `high` = "High",
                                                                           `1000` = "n = 1,000",
                                                                           `2000` = "n = 2,000",
                                                                           `10000` = "n = 10,000",
                                                                           `medium` = "Medium",
                                                                           `low` = "Low"
                                                                           
  ))) +
  xlab("Number of CCA components (r)") + 
  ylab(expression("Subspace distance")) +
  labs(colour = "Method") + 
  theme(legend.position = "none")
ggsave("effect_of_r.pdf", width = 6.5, height = 2.4)



############## Groups #################

summ <- read_csv( "Fits/results_groups_sparse_rrr_experiments.csv")

legend_order <- c("Oracle", 
                  "RRR-ADMM-opt", 
                  "CVX-opt-group", 
                  "FIT_SAR_CV",
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "Waaijenborg-CV", "Waaijenborg-Author", "SCCA_Parkhomenko")
labels_n <- c("Oracle",  
              "\nsparse-CCA-RRR\n(this paper)",
              "\ngroup-CCA-RRR\n(this paper)",
              "\nSAR\n(Wilms et al)",
              "\nSAR with BIC\n(Wilms et al)",
              "\nCCA-Lasso, permuted\n(Witten et al)",
              "\nCCA-Lasso\n(Witten et al)",
              "\nSSCA\n(Waaijenborg et al)",
              "\nSSCA with BIC\n(Waaijenborg et al)",
              "\nSCCA\n(Parkhomenko et al)")

cbPalette <- c("#000000",
               "#009E73",
               "#999933",
               "#D55E00",
               "#CC79A7",
               "#0072B2",
               "#56B4E9",
               "#E69F00", 
               "#F0E442",
               "#999999")

summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))

ggplot(summ %>% filter(r_pca == 5, r==2,
                       nnzeros == 5, 
                       p1 <= 3000,
                       p2 <= 50,
                       n == 500,
                       method %in% legend_order),
       aes(x=p1, y = distance_tot_mean, colour = method, fill = method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  #geom_ribbon(aes(ymin=distance_tot_q975, ymax=distance_tot_q025, color = NULL), alpha = 0.2)+
  #geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
  #                  colour =method), width=0.05, alpha=0.5)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(p2~theta_strength, labeller = as_labeller(c(`10` = "q = 10",
                                                         `30` = "q = 30",
                                                         `50` = "q = 50",
                                                         `70` = "q = 70",
                                                         `100` = "n = 100",
                                                         `200` = "n = 200",
                                                         `300` = "n = 300",
                                                         `500` = "n = 500",
                                                         `high` = "High",
                                                         `1000` = "n = 1,000",
                                                         `2000` = "n = 2,000",
                                                         `10000` = "n = 10,000",
                                                         `medium` = "Medium",
                                                         `low` = "Low"))) +
  xlab("Number of covariates in X (p)") + 
  ylab(expression("Subspace distance")) +
  #scale_y_log10()+
  #scale_x_log10()+
  ylim(0, 2.8)+
  theme(legend.position = "top")+
  theme(legend.key.size = unit(0.1, 'cm'), legend.text = element_text(size=10), legend.title = element_blank())+
  guides(color=guide_legend(nrow=3))
ggsave("final_rrr_group_no_bars.pdf", width = 6.5, height = 7)


ggplot(summ %>% filter(r_pca == 5, r==2,
                       nnzeros == 5, 
                       p1 <= 3000,
                       p2 <= 50,
                       n == 500,
                       method %in% legend_order) %>%
         mutate(p = factor(p1)),
       aes(x=FNR_mean, y = FPR_mean, colour = method, fill = method, shape = p)) +
  geom_point()+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(p2~theta_strength, labeller = as_labeller(c(`10` = "q = 10",
                                                         `30` = "q = 30",
                                                         `50` = "q = 50",
                                                         `70` = "q = 70",
                                                         `100` = "n = 100",
                                                         `200` = "n = 200",
                                                         `300` = "n = 300",
                                                         `500` = "n = 500",
                                                         `high` = "High",
                                                         `1000` = "n = 1,000",
                                                         `2000` = "n = 2,000",
                                                         `10000` = "n = 10,000",
                                                         `medium` = "Medium",
                                                         `low` = "Low"))) +
  xlab("False Negative Rate (FNR)") + 
  ylab("False Positive Rate (FPR)") +
  theme(legend.position = "top")+
  guides(color = FALSE, fill = FALSE)+
  #theme(legend.key.size = unit(0.1, 'cm'), legend.text = element_text(size=12), legend.title = element_blank())+
  guides(shape = guide_legend(nrow = 1), size = "none")
ggsave("FPR_FNR_rrr_group_dots_shapes.pdf", width = 6.5, height = 5.5)


rbind(summ %>% filter(r_pca == 5, r==2,
                      nnzeros == 5, 
                      p1 <= 3000,
                      n == 500,
                      theta_strength == "medium",
                      method %in% legend_order) %>%
        select(p1, p2, method, FNR_mean) %>%
        mutate(rate = "FNR") %>%
        rename(value = FNR_mean),
      summ %>% filter(r_pca == 5, r==2,
                      nnzeros == 5, 
                      p1 <= 3000,
                      n == 500,
                      theta_strength == "medium",
                      method %in% legend_order) %>%
        select(p1, p2, method, FPR_mean) %>%
        mutate(rate = "FPR") %>%
        rename(value = FPR_mean)) %>%
  ggplot(aes(x=p1, y = value, colour = method, fill = method)) +
  geom_point(size = 2)+
  geom_line(linewidth=0.8)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(p2~rate, labeller = as_labeller(c(`10` = "q = 10",
                                               `30` = "q = 30",
                                               `50` = "q = 50",
                                               `70` = "q = 70",
                                               `100` = "n = 100",
                                               `200` = "n = 200",
                                               `300` = "n = 300",
                                               `500` = "n = 500",
                                               `high` = "High",
                                               `1000` = "n = 1,000",
                                               `2000` = "n = 2,000",
                                               `10000` = "n = 10,000",
                                               `medium` = "Medium",
                                               `low` = "Low",
                                               `FPR` = "FPR",
                                               `FNR` = "FNR"))) +
  xlab("Number of covariates in X (p)") + 
  ylab("Rate") +
  theme(legend.position = "top")+
  theme(legend.key.size = unit(1, 'cm'), legend.text = element_text(size=10), legend.title = element_blank())+
  guides(color=guide_legend(nrow=4))
ggsave("FNR_FPR_rrr_group_no_bars.pdf", width = 6.5, height = 9.7)


############Graphs###############

summ <- read_csv( "Fits/results_graph_sparse_rrr_experiments.csv")

legend_order <- c("Oracle", 
                  "RRR-ADMM-opt", 
                  "CVX-opt-graph", 
                  "FIT_SAR_CV",
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "Waaijenborg-CV", "Waaijenborg-Author", "SCCA_Parkhomenko")
labels_n <- c("Oracle",  
              "\nsparse-CCA-RRR\n(this paper)",
              "\ngraph-CCA-RRR\n(this paper)",
              "\nSAR\n(Wilms et al)",
              "\nSAR with BIC\n(Wilms et al)",
              "\nCCA-Lasso, permuted\n(Witten et al)",
              "\nCCA-Lasso\n(Witten et al)",
              "\nSSCA\n(Waaijenborg et al)",
              "\nSSCA with BIC\n(Waaijenborg et al)",
              "\nSCCA\n(Parkhomenko et al)")

cbPalette <- c("#000000",
               "#009E73",
               "#117733",
               "#D55E00",
               "#CC79A7",
               "#0072B2",
               "#56B4E9",
               "#E69F00", 
               "#F0E442",
               "#999999")

summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))


ggplot(summ %>% filter(r_pca == 5, r==2,
                       nnzeros == 5, 
                       p1 <= 2000,
                       p2 <= 50,
                       n == 500,
                       method %in% legend_order,
                       method != "Oracle"),
       aes(x=p1, y = distance_tot_mean, colour = method, fill = method)) +
  geom_point(size=2)+
  geom_line(linewidth=0.8)+
  scale_color_manual(values = cbPalette, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = cbPalette, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(p2~theta_strength, labeller = as_labeller(c(`10` = "q = 10",
                                                         `30` = "q = 30",
                                                         `50` = "q = 50",
                                                         `80` = "q = 80",
                                                         `100` = "n = 100",
                                                         `200` = "n = 200",
                                                         `300` = "n = 300",
                                                         `500` = "n = 500",
                                                         `high` = "High",
                                                         `1000` = "n = 1,000",
                                                         `2000` = "n = 2,000",
                                                         `10000` = "n = 10,000",
                                                         `medium` = "Medium",
                                                         `low` = "Low"))) +
  xlab("Number of covariates in X (p)") + 
  ylab(expression("Subspace distance")) +
  ylim(0, 2.8)+
  xlim(100, 1600)+
  theme(legend.position = "top")+
  theme(legend.key.size = unit(0.1, 'cm'), legend.text = element_text(size=10), legend.title = element_blank())+
  guides(color=guide_legend(nrow=3))
ggsave("final_rrr_graph_no_bars.pdf", width = 6.5, height = 7)
