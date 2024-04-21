library(tidyverse)
library(ggplot2)
library(pracma)
library(tidyverse)
theme_set(theme_bw(base_size = 14))

file_list <- list.files(path =  "experiments/simulations/results/group", 
                         pattern = "2024-group-newest_RRR_efficient_results14*", full.names = TRUE)
results <- bind_rows(lapply(file_list, read.csv))




legend_order <- c("Oracle",  "FIT_SAR_CV",
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "SCCA_Parkhomenko", "Waaijenborg-CV", "Waaijenborg-Author",
                  #"RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",
                  #"RRR-opt",    
                  "RRR-ADMM-opt", "CVX-opt-group"   )
my_colors <- c( "black", "red", "indianred4",
                "orange", "yellow", "chartreuse2",
                "burlywood2", "burlywood4",
                # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
                "navyblue", #"cyan", 
                "dodgerblue")

labels_n <-    c("Oracle",  "SAR CV (Wilms et al)",
                 "SAR BIC (Wilms et al)",
                 "Sparse CCA, permuted\n(Witten et al)",
                 "Sparse CCA with CV\n(Witten et al)",
                 "SCCA (Parkhomenko et al)", "Sparse CCA with CV\n(Waaijenborg et al)",
                 "Sparse CCA(Waaijenborg et al)",
                 # "RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",
                 "RRR-CCA (this paper)",
                 "group-RRR-CCA (this paper)")
theme_set(theme_bw(base_size = 18))


unique(summ$method)
colnames(results)
colnames(summ)
unique(summ$nnzeros)
unique(summ$r)
unique(summ$r_pca)
unique(summ$p1)
unique(summ$p2)
unique(summ$counts)

summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))


ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==5, 
                        method %in% legend_order
),
aes(x=p1, 
    y = time_mean, 
    colour =method)) +
  geom_point()+
  geom_line()+
  #geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
  #                  colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
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
                                                                         `low` = "Low"
                                                                         
  ))) +
  xlab("p") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==5, 
                        method %in% legend_order
),
aes(x=p1, 
    y = time_mean, 
    colour =method)) +
  geom_point()+
  geom_line()+
  #geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
  #                  colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
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
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Time (in seconds)")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(summ %>% filter( r_pca == 5, r==2,
                        nnzeros==5, 
                        method %in% legend_order
),
aes(x=p1, 
    y = FNR_mean, 
    colour =method)) +
  geom_point(size=1.6)+
  geom_line()+
  #geom_jitter(width = 0.2, height = 0) + # Jitter only in the x-direction
  geom_errorbar(aes(ymin=FNR_q25, ymax=FNR_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
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
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("False Negative Rate")) +
  labs(colour="Method") + 
  #scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


r_pca = 0
ggplot(summ %>% filter( r_pca == r_pca, 
                        nnzeros==5,
                        str_detect(method, "group-RRR"), !str_detect(method, "opt")),
       aes(x=lambda, 
       y = distance_tot_q50, 
       colour ="group")) +
  geom_point()+
  geom_line() +
  geom_point(data = summ %>% filter( r_pca == r_pca, 
                                      nnzeros==5,
                                      str_detect(method, "RRR"), 
                                      !str_detect(method, "opt"),
                                      !str_detect(method, "CVX")),
                     aes(x=lambda, 
                         y = distance_tot_q50, 
                         colour ="RRR"))+
  facet_grid(theta_strength~ p1 + n + nnzeros, scales = "free") +
  xlab("lambda") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

  geom_errorbar(aes(ymin=distance_tot_q25, ymax=distance_tot_q75,
                    colour =method), width=0.05)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) 


