library(tidyverse)
library(igraph)
library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(knitr)
library(kableExtra) 
library(ellipse)
library(viridis)
library(reshape2)
library(mclust)
library(combinat)

source('code/reduced_rank_regression.R')
source('code/eval.R')
source('code/graph_reduced_rank_regression.R')
source('code/alternative_methods/SAR.R')
source('code/alternative_methods/Parkhomenko.R')
source('code/alternative_methods/Witten_CrossValidation.R')
source('code/alternative_methods/Waaijenborg.R')
source('code/experiment_functions.R')

############# Prepare data ##############

data_presidents = read_csv("Data/1976-2020-president.csv") %>% 
  filter(year > 2007, is.na(candidate) == FALSE) %>%
  mutate(percentage_votes = candidatevotes/ totalvotes)

extract_format = function(name) {
  parts = strsplit(name, ", ")[[1]]
  last_name = parts[1]
  first_initial = substr(parts[2], 1, 2)
  return(paste(last_name, first_initial, sep = ", "))
}

X = data_presidents %>% 
  dplyr::select(year, state_po, candidate, percentage_votes) %>%
  group_by(year, state_po, candidate) %>%
  summarise(percentage_votes = sum(percentage_votes)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("year", "candidate"),
  names_from = state_po,
  values_from =  percentage_votes,
  values_fill = 0) 
X$candidate_simplified = sapply(X$candidate, extract_format)
X = X %>% dplyr::select(-candidate)

Y = read_csv("Data/politicians_positions.csv")
data = merge(X, Y, by = c("year", "candidate_simplified"))

state_names = colnames(X)[3:(ncol(X)-1)]
info = data[,1:2]
X = data[, 3:(ncol(data)-10)] %>% dplyr::select(-HI, -AK)
Y = data[, (ncol(data)-9):ncol(data)]

supp_info = data_presidents %>% 
  dplyr::select(year, candidate, party_detailed) %>% group_by(year, candidate, party_detailed) %>%
  mutate(n = n()) %>% ungroup() %>%
  distinct() %>% group_by(year, candidate) %>%
  slice(which.max(n)) %>% dplyr::select(1:3)
supp_info$candidate_simplified = sapply(supp_info$candidate, extract_format) 
supp_info = supp_info %>% data.frame() %>% dplyr::select(-candidate)
info = left_join(info, supp_info)

logit_smooth = function(p) {
  # Adjust values exactly equal to 0 or 1
  minX = 1e-3
  p[p == 0] <- minX
  log(p / (1 - p))
}

X = X %>% mutate(across(everything(), logit_smooth)) 
X = scale(as.matrix(X), center = T, scale = T)
Y = scale(as.matrix(Y), center = T, scale = T)
n = nrow(X)
p = ncol(X)
q = ncol(Y)

#### Download the adjacency matrix

get_edge_incidence = function(g, weight = 1){
  n_nodes = vcount(g)
  d_max = max(degree(g))
  #d_max = 1
  edges = data.frame(as_edgelist(g)) %>%
    arrange(X1, X2)
  Gamma = matrix(0, nrow(edges), n_nodes)
  # Make beta_v into a matrix
  names_st = unique(c(edges$X1, edges$X2))
  for (e in 1:nrow(edges)){
    ind1 = which( edges$X1[e] == names_st)
    ind2 = which( edges$X2[e] == names_st)
    Gamma[e, ind1] = weight
    Gamma[e,ind2] = - weight
  }
  return(Gamma)
}

A = read.table("Data/state_adjacency.csv", sep = ",", header = T, row.names = 1)
A = A[colnames(X), colnames(X)]

g = graph_from_adjacency_matrix(as.matrix(A), mode = "undirected")
plot(g,
     vertex.label = V(g)$name,  # Node labels
     vertex.size = 10,  # Adjust node size
     vertex.color = "skyblue",  # Node color
     edge.color = "grey",  # Edge color
     vertex.label.color = "black",  # Node label color
     edge.label.color = "red",  # Edge label color
     vertex.label.dist = 1.5,  # Distance of labels from nodes
     edge.label.cex = 0.8,  # Edge label size
     vertex.label.cex = 1.2  # Node label size
)
Gamma = get_edge_incidence(g, weight = 1)

#### Apply the algorithm

r = 2
set.seed(5)
#folds = createFolds(1:n, 6)
folds = createFolds(1:n, n) # for LOO
results = c()

lambdas = 10^seq(-3, 1, length.out = 30)
lambdas_comp = 10^seq(-3, 1, length.out = 30)

for (i in 1:length(folds)){
  test_ind = ifelse(i < length(folds), i + 1, (i+1) %% length(folds))
  val_ind = ifelse(i < length(folds) - 1, i + 2, (i+2) %% length(folds))
  test = folds[[test_ind]]
  val = folds[[val_ind]]
  train = -c(test, val)
  print(c(i, test_ind, val_ind))
  #for (lambda in lambdas_comp){  
    for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                     "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                     "SCCA_Parkhomenko")){
      print(paste0("Starting ", method))
      tryCatch({
        SCCA = additional_checks(#X[train,], Y[train,], 
                                  X[-val,], Y[-val,],
                                  S=NULL, 
                                  rank=r, kfolds = 17, 
                                  method.type = method,
                                  lambdax = lambdas_comp,
                                  lambday = 1e-4)
                                  #lambdax = lambda,
                                  #lambday = 10^seq(-4,-3, length.out = 5))
        results = rbind(results,
                             data.frame(method,
                                        lambda = 0,
                                        #lambda = lambda,
                                        fold = i,
                                        eval(X, Y, SCCA$u, SCCA$v, train, test, val)
                             ))
      }, error = function(e) {
        # Print the error message
        cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
        # Skip to the next iteration
      })
    }
  #}
  for (lambda in lambdas){
    RRR = CCA_graph_rrr(X[train,], Y[train,],
                        Gamma,
                        Sx=NULL, Sy=NULL, Sxy = NULL,
                        lambda = lambda,
                        Kx=NULL, r=r,
                        rho=1, niter=2 * 1e4,
                        do.scale = FALSE, lambda_Kx=0,
                        thresh=1e-6,
                        LW_Sy = FALSE)

    results = rbind(results,
                         data.frame(method = "graph-CCA-RRR",
                                    lambda,
                                    fold = i,
                                    eval(X, Y, RRR$U, RRR$V, train, test, val)
                         ))
    RRR = CCA_rrr(X[train,], Y[train,],
                  lambda = lambda,
                  Kx=NULL, r=r,
                  rho = 1, niter=2 * 1e4,
                  do.scale = FALSE, lambda_Kx = 0,
                  thresh=1e-6,
                  solver= "ADMM",
                  LW_Sy = FALSE)
    results = rbind(results,
                         data.frame(method = "sparse-CCA-RRR",
                                    lambda,
                                    fold = i,
                                    eval(X, Y, RRR$U, RRR$V, train, test, val)
                         ))
  }
}

saveRDS(results, "election_results_loo.rds")

###### Optimal lambda

results = readRDS("election_results_loo.rds")

ggplot(results, aes(lambda, test_mse, color = factor(fold)))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  facet_grid(method~comp, scale = "free")

summary_results = results %>%
  group_by(method, lambda) %>%
  summarise(train_cor = mean(train_cor, na.rm = T),
            train_mse = mean(train_mse, na.rm = T),
            test_cor = mean(test_cor, na.rm = T),
            test_mse = mean(test_mse, na.rm = T),
            test_sd = mean(test_sd, na.rm = T),
            val_cor = mean(val_cor, na.rm = T),
            val_mse = mean(val_mse, na.rm = T),
            val_sd = mean(val_sd, na.rm = T)) %>%
  ungroup()

ggplot(summary_results %>% filter(lambda < 0.7), aes(lambda, test_mse))+
  geom_point()+
  geom_line()+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~method, scale = "free")

summary_results %>% filter(lambda < 0.7) %>%
  group_by(method) %>%
  slice(which.min(test_mse)) %>% dplyr::select(val_sd)

lambda_opt = summary_results %>% filter(lambda < 0.3) %>%
  group_by(method) %>%
  slice(which.min(test_mse)) %>% 
  dplyr::select(method, lambda)

lambda_opt


############# Plot results

legend_order <- c("LIBERTARIAN", 
                  "REPUBLICAN", "GREEN", 
                  "INDEPENDENT", "DEMOCRAT", 
                  "CONSTITUTION PARTY")

my_colors <- c("#D55E00", 
               "red",  "#009E73", 
               "#999999", "dodgerblue", 
               "#332288")

labels_n <- c("Libertarian",
                 "Republican", "Green", "Independent",
                 "Democratic", "Constitution")

theme_set(theme_bw(base_size = 18))
ellipse.level = 0.95


plot_results = function(XU, U, V, name){
  
  df = data.frame(XU, name = info$candidate_simplified, party = info$party_detailed)
  df_green = df %>% filter(party == "GREEN") %>% dplyr::select(XU1, XU2)
  df_con = df %>% filter(party == "CONSTITUTION PARTY") %>% dplyr::select(XU1, XU2)

  ggplot(df, aes(x=XU1, y=XU2, color=party))+
    geom_point(size = 3)+
    ggrepel::geom_text_repel(aes(label = name,  colour = party), show.legend = FALSE, size = 3)+
    scale_color_manual(values = my_colors, breaks = legend_order,
                       labels = labels_n) +
    scale_fill_manual(values = my_colors, breaks = legend_order,
                       labels = labels_n) +
    stat_ellipse(level = ellipse.level)+
    geom_path(data = data.frame(ellipse(cov(df_green), 
                                      centre=colMeans(df_green),
                                      level = ellipse.level)),
              aes(x = XU1, y = XU2, color = "GREEN"))+
    xlab(expression(Xu[1]))+ 
    ylab(expression(Xu[2]))+
    labs(colour = "Party")+
    theme(legend.position = "none")#+
    #guides(color=guide_legend(nrow=2))
  ggsave(paste0("election-", name, "-var1-loo.pdf"), width = 4, height = 4)

  state_match = data_presidents %>% 
    dplyr::select(state_po, state) %>%
    mutate(state = tolower(state)) %>% 
    distinct() %>%
    right_join(data.frame(state_po = colnames(X)))
  map_data = map_data("state")
  dfU = merge(map_data, data.frame(U, state = state_match$state) %>% melt(id = c("state")), by.x = "region", by.y = "state")
  
  dfU %>% filter(variable == "U1") %>%
    mutate(variable = factor(variable,
                             levels = paste0("U", 1),
                             labels = c(expression(u[1]))))%>%
  ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = value),
                 color = "black") +
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    guides(x = "none", y = "none")+
    facet_grid(cols = vars(variable), scale = "free", labeller = label_parsed)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_gradient2(high = "#F0E442", low = "#332288")
  ggsave(paste0("election-", name, "-u1-loo.pdf"), width = 3.5, height = 2)
  
  dfU %>% filter(variable == "U2") %>%
    mutate(variable = factor(variable,
                             levels = paste0("U", 2),
                             labels = c(expression(u[2]))))%>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = value),
                 color = "black") +
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    guides(x = "none", y = "none")+
    facet_grid(cols = vars(variable), scale = "free", labeller = label_parsed)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_gradient2(high = "#F0E442", low = "#332288")
  ggsave(paste0("election-", name, "-u2-loo.pdf"), width = 3.5, height = 2)
  
  new_labels = c("Should Abortion\nRemain Legal?", 
                  "Should the Death Penalty\nBe Allowed?",
                  "Should Former Felons \nBe Allowed to Vote?",
                  "Should Federal Taxes \nBe Increased?" ,
                  "Should the US Expand \nIts Nuclear Power?",
                  "Are More Regulations \nOn Guns Needed?",
                  "Should the US Build a Fence\nAlong the Mexico Border?",
                  "Is Obamacare\nGood for America?",
                  "Are humans responsible \nfor global climate change?",
                  "Should the US tax \nCarbon Emissions?")
  dfV = data.frame(V, opinion = new_labels) %>% melt(id = "opinion")
  dfV %>% mutate(variable = factor(variable,
                                   levels = paste0("V", 1:r),
                                   labels = c(expression(v[1]), expression(v[2])))) %>%
  ggplot()+
    geom_bar(aes(value, opinion, fill = variable), stat = "identity")+
    xlab("Value")+
    ylab(NULL)+
    facet_grid(~variable, labeller = label_parsed)+
    theme(legend.position = "none")+
    scale_x_continuous(breaks = c(-0.5, 0, 0.5, 1))
  ggsave(paste0("election-", name, "-v-loo.pdf"), width = 7, height = 7)
}

RRR = CCA_graph_rrr(X, Y,
                    Gamma,
                    Sx=NULL, Sy=NULL, Sxy = NULL,
                    lambda = lambda_opt %>% filter(method == "graph-CCA-RRR") %>% pull(lambda),
                    Kx=NULL, r=r,
                    rho=1, niter=2 * 1e4,
                    do.scale = FALSE, lambda_Kx=0,
                    thresh=1e-6,
                    LW_Sy = FALSE)
XU = X %*% RRR$U
colnames(XU) = paste0("XU", 1:r)
U = RRR$U
colnames(U) = paste0("U", 1:r)
V = RRR$V
colnames(V) = paste0("V", 1:r)
plot_results(XU, U, V, "rrr")

SCCA = additional_checks(X, Y,
                        S=NULL, 
                        rank=r, kfolds = n, 
                        method.type = "FIT_SAR_CV",
                        lambdax = lambdas_comp,
                        lambday = 1e-4)
XU = X %*% SCCA$u
colnames(XU) = paste0("XU", 1:r)
U = SCCA$u
colnames(U) = paste0("U", 1:r)
V = SCCA$v
colnames(V) = paste0("V", 1:r)
plot_results(XU, U, V, "sar")

SCCA = additional_checks(X, Y,
                         S=NULL, 
                         rank=r, kfolds = n, 
                         method.type = "Witten_Perm",
                         lambdax = lambdas_comp,
                         lambday = 1e-4)

XU = X %*% SCCA$u
colnames(XU) = paste0("XU", 1:r)
U = SCCA$u
colnames(U) = paste0("U", 1:r)
V = SCCA$v
colnames(V) = paste0("V", 1:r)
plot_results(XU, U, V, "lasso")

######### Clustering accuracy

clus_acc = function(U){
  set.seed(0)
  XU = X %*% U
  colnames(XU) = paste0("XU", 1:r)
  df = data.frame(XU, party = info$party_detailed)
  df = df %>% filter(party %in% c("LIBERTARIAN", "REPUBLICAN", "DEMOCRAT", "GREEN", "CONSTITUTION PARTY"))
  
  model = Mclust(df %>% dplyr::select(XU1, XU2), G = 5:5)
  CM = confusionMatrix(factor(model$classification), factor(df$party, labels = 1:5))$table
  ca =  cm_accuracy(CM) 
  return(ca)
}

RRR = CCA_graph_rrr(X, Y,
                    Gamma,
                    Sx=NULL, Sy=NULL, Sxy = NULL,
                    lambda = lambda_opt %>% filter(method == "graph-CCA-RRR") %>% pull(lambda),
                    Kx=NULL, r=r,
                    rho=1, niter=2 * 1e4,
                    do.scale = FALSE, lambda_Kx=0,
                    thresh=1e-6,
                    LW_Sy = FALSE)
clustering_accuracy = data.frame(method = "graph-CCA-RRR", accuracy = clus_acc(RRR$U))

RRR = CCA_rrr(X, Y,
              lambda = lambda_opt %>% filter(method == "sparse-CCA-RRR") %>% pull(lambda),
              Kx=NULL, r=r,
              rho = 1, niter=2 * 1e4,
              do.scale = FALSE, lambda_Kx = 0,
              thresh=1e-6,
              solver= "ADMM",
              LW_Sy = FALSE)
clustering_accuracy = rbind(clustering_accuracy, data.frame(method = "sparse-CCA-RRR", accuracy = clus_acc(RRR$U)))


for(method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                 "SCCA_Parkhomenko")){
  print(paste0("Starting ", method))
  tryCatch({
    SCCA = additional_checks(X, Y,
                              S=NULL, 
                              rank=r, kfolds = n, 
                              method.type = method,
                              lambdax = lambdas_comp,
                              lambday = 1e-4)
    clustering_accuracy = rbind(clustering_accuracy, data.frame(method, accuracy = clus_acc(SCCA$u)))
  }, error = function(e) {
    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
  })
}


saveRDS(clustering_accuracy, "pol_clustering_accuracy.rds")

clustering_accuracy
