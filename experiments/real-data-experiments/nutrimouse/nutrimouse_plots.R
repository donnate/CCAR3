library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(caret)
library(knitr)
library(kableExtra)
library(mclust)
library(combinat)


source('code/eval.R')
source('code/reduced_rank_regression.R')
source('code/alternative_methods/SAR.R')
source('code/alternative_methods/Parkhomenko.R')
source('code/alternative_methods/Witten_CrossValidation.R')
source('code/alternative_methods/Waaijenborg.R')
source('code/experiment_functions.R')

##### Load data

data(nutrimouse)
X = as.matrix(nutrimouse$gene %>% mutate_all(~rank(.)))
Y = as.matrix(nutrimouse$lipid %>% mutate_all(~rank(.)))
X = scale(X, center = T, scale = T)
Y = scale(Y, center = T, scale = T)
n = nrow(X)
p = ncol(X)
q = ncol(Y)


#### Apply CCA

lambdas = 10^seq(-3, 0, length.out = 30)
lambdas_comp = 10^seq(-3, 3, length.out = 30)
r = 5
penalty="l21"
solver="rrr"
LW_Sy = TRUE

##### Cross validation

set.seed(0)

folds = createFolds(1:n, 8)
results = c()
for (i in  1:8){
  test_ind = ifelse(i < length(folds), i + 1, (i+1) %% length(folds))
  val_ind = ifelse(i < length(folds) - 1, i + 2, (i+2) %% length(folds))
  print(c(i, test_ind, val_ind))
  for (lambda in lambdas){
    test = folds[[test_ind]]
    val = folds[[val_ind]]
    train = -c(test, val)
    RRR = CCA_rrr(X[train,], Y[train,],  
                    lambda = lambda, 
                    Kx=NULL, r=r,
                    rho = 1, niter=2 * 1e4,
                    do.scale = FALSE, lambda_Kx=0,
                    thresh=1e-6,
                    solver= "ADMM",
                    LW_Sy = LW_Sy)
    
    results = rbind(results,
                         data.frame(method = "sparse-CCA-RRR",
                          lambda,
                          fold = i,
                          eval(X, Y, RRR$U, RRR$V, train, test, val)
                        ))
  }
  for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                   "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                   "SCCA_Parkhomenko")){
    print(paste0("Starting ", method))
    tryCatch({
      SCCA = additional_checks(X[-val,], Y[-val,],
                               S=NULL, 
                               rank=r, kfolds=8, 
                               method.type = method,
                               lambdax = lambdas_comp,
                               lambday =  0)
      results = rbind(results,
                           data.frame(method,
                                      lambda = 0,
                                      fold = i,
                                      eval(X, Y, SCCA$u, SCCA$v, train, test, val)
                           ))
    }, error = function(e) {
      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
    })
  }
}
saveRDS(results, "nutrimouse_results.rds")

##### Optimal lambda

results = readRDS("nutrimouse_results.rds")

summary_results = results %>% filter(lambda < 1) %>%
  group_by(method, lambda) %>%
  summarise(train_cor = mean(train_cor, na.rm = T),
            train_mse = mean(train_mse, na.rm = T),
            test_cor = mean(test_cor, na.rm = T),
            test_mse = mean(test_mse, na.rm = T),
            test_sd = mean(test_sd, na.rm = T),
            val_cor = mean(val_cor, na.rm = T),
            val_mse = mean(val_mse, na.rm = T),
            val_sd = mean(val_sd, na.rm = T)
            ) %>%
  ungroup()

ggplot(summary_results %>% filter(method == "sparse-CCA-RRR"), aes(lambda, test_mse))+
 geom_line()+
 geom_point()+
 scale_x_log10()+
 scale_y_log10()

summary_results %>% group_by(method) %>% 
  slice(which.min(test_mse)) %>% dplyr::select(val_cor, val_sd)

lambda_opt = summary_results %>% filter(method == "sparse-CCA-RRR") %>% 
  slice(which.min(test_mse)) %>% pull(lambda) %>% as.numeric()

lambda_opt

##### Projection plots

legend_order = c("lin", "sun", "fish",
                  "ref", "coc")
my_colors = c(  "red","orange",  "dodgerblue", 
                 "black", "brown")
labels_n =   c("LIN", "SUN", "FISH",
                "REF", "COC")
ellipse.level = 0.95

RRR = CCA_rrr(X, Y,  
                lambda = lambda_opt, 
                Kx = NULL, r=r,
                rho=1, niter=2 * 1e4,
                do.scale = FALSE, lambda_Kx=0,
                thresh=1e-6,
                solver= "ADMM",
                LW_Sy = LW_Sy)

XU = X %*% RRR$U
colnames(XU) = paste0("XU", 1:5)
YV = Y %*% RRR$V
colnames(YV) = paste0("YV", 1:5)
df = data.frame(XU, YV, gen =  nutrimouse$genotype, diet = nutrimouse$diet)
U = RRR$U
V = RRR$V
colnames(U) = paste0("U", 1:5)
colnames(V) = paste0("V", 1:5)

ggplot(df)+
  geom_point(aes(x=XU1, y=XU2), shape = 1, size = 2)+
  geom_point(aes(x=YV1, y=YV2), shape = 1, size = 2)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_segment(aes(x = XU1, y = XU2, xend = YV1, yend = YV2),
               arrow = arrow(length = unit(0.3, "cm")))+
  xlab("Canonical variate 1") + 
  ylab("Canonical variate 2")+
  labs(colour = "Diet", shape = "Genotype")
ggsave("nutrimouse-rrr-arrow.pdf", width = 4, height = 4)

ggplot(df, aes(x=XU1, y=XU2, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  scale_fill_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  # ggforce::geom_mark_ellipse(aes(fill = diet,
  #                                color = diet))+
  xlab(expression(Xu[1]))+ 
  ylab(expression(Xu[2]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
  #guides(color=guide_legend(nrow=2), shape = guide_legend(nrow=2))
ggsave("nutrimouse-rrr-var1.pdf", width = 4, height = 4)

ggplot(df, aes(x=XU3, y=XU4, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  xlab(expression(Xu[3]))+ 
  ylab(expression(Xu[4]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
ggsave("nutrimouse-rrr-var2.pdf", width = 4, height = 4)

ggplot(df, aes(x=XU1, y=XU5, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  xlab(expression(Xu[1]))+ 
  ylab(expression(Xu[5]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
ggsave("nutrimouse-rrr-var3.pdf", width = 4, height = 4)


##### Directions plots

dfU1 = data.frame(U, name = colnames(X))
dfU1 %>% reshape2::melt(id = "name") %>% filter(abs(value) > 0.1) %>%
  mutate(variable = factor(variable,
                           levels = paste0("U", 1:5),
                           labels = c(expression(u[1]), expression(u[2]), expression(u[3]), expression(u[4]), expression(u[5])))) %>% 
  ggplot() +
  geom_bar(aes(value, name, fill = variable), stat = "identity", show.legend = FALSE)+
  facet_grid(cols = vars(variable), labeller = label_parsed) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Genes", x = "Value")+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))
ggsave("nutrimouse-rrr-u.pdf", width = 10, height = 6.5)

dfV2 = data.frame(V, name = colnames(Y))
dfV2 %>% reshape2::melt(id = "name") %>% filter(abs(value) > 0.1) %>%
  mutate(variable = factor(variable,
                           levels = paste0("V", 1:5),
                           labels = c(expression(v[1]), expression(v[2]), expression(v[3]), expression(v[4]), expression(v[5])))) %>% 
  ggplot() +
  geom_bar(aes(value, name, fill = variable), stat = "identity", show.legend = FALSE)+
  facet_grid(cols = vars(variable), labeller = label_parsed) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Hepatic Fatty Acids", x = "Value")+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))+
  scale_x_continuous(breaks = seq(-0.3, 0.5, 0.2))
ggsave("nutrimouse-rrr-v.pdf", width = 10, height = 6)

SCCA = additional_checks(X, Y,  
                         S = NULL, 
                         rank=r, kfolds = 8, 
                         method.type =  "FIT_SAR_CV",
                         lambdax= lambdas_comp,
                         lambday = 0)

XU = X %*% SCCA$u
colnames(XU) = paste0("XU", 1:5)
YV = Y %*% SCCA$v
colnames(YV) = paste0("YV", 1:5)
df = data.frame(XU, YV, gen =  nutrimouse$genotype, diet = nutrimouse$diet)
U = SCCA$u
V = SCCA$v
colnames(U) = paste0("U", 1:5)
colnames(V) = paste0("V", 1:5)

ggplot(df)+
  geom_point(aes(x=XU1, y=XU2), shape = 1, size = 2)+
  geom_point(aes(x=YV1, y=YV2), shape = 1, size = 2)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_segment(aes(x = XU1, y = XU2, xend = YV1, yend = YV2),
               arrow = arrow(length = unit(0.3, "cm")))+
  xlab("Canonical variate 1") + 
  ylab("Canonical variate 2")+
  labs(colour = "Diet", shape = "Genotype")
ggsave("nutrimouse-scca-arrow.pdf", width = 4, height = 4)

ggplot(df, aes(x=XU1, y=XU2, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  xlab(expression(Xu[1]))+ 
  ylab(expression(Xu[2]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
ggsave("nutrimouse-scca-var1.pdf", width = 4, height = 4)

ggplot(df, aes(x=XU3, y=XU4, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  xlab(expression(Xu[3]))+ 
  ylab(expression(Xu[4]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
ggsave("nutrimouse-scca-var2.pdf", width =4, height = 4)

ggplot(df, aes(x=XU1, y=XU5, colour=diet))+
  geom_point(aes(shape=gen), size = 4)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  stat_ellipse(level = 0.95)+
  xlab(expression(Xu[1]))+ 
  ylab(expression(Xu[5]))+
  labs(colour = "Diet", shape = "Genotype")+
  theme(legend.position = "none")
ggsave("nutrimouse-scca-var3.pdf", width = 4, height = 4)


dfU2 = data.frame(U, name = colnames(X))
dfU2 %>% reshape2::melt(id = "name") %>% filter(abs(value) > 0.1) %>%
  mutate(variable = factor(variable,
                           levels = paste0("U", 1:5),
                           labels = c(expression(u[1]), expression(u[2]), expression(u[3]), expression(u[4]), expression(u[5])))) %>% 
  ggplot() +
  geom_bar(aes(value, name, fill = variable), stat = "identity", show.legend = FALSE)+
  facet_grid(cols = vars(variable), labeller = label_parsed) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Genes", x = "Value")+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))
ggsave("nutrimouse-scca-u.pdf", width = 10, height = 6.5)

dfV2 = data.frame(V, name = colnames(Y))
dfV2 %>% reshape2::melt(id = "name") %>% filter(abs(value) > 0.1) %>%
  mutate(variable = factor(variable,
                           levels = paste0("V", 1:5),
                           labels = c(expression(v[1]), expression(v[2]), expression(v[3]), expression(v[4]), expression(v[5])))) %>% 
  ggplot() +
  geom_bar(aes(value, name, fill = variable), stat = "identity", show.legend = FALSE)+
  facet_grid(cols = vars(variable), labeller = label_parsed) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Hepatic Fatty Acids", x = "Value")+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11))#+
ggsave("nutrimouse-scca-v.pdf", width = 10, height = 6)

##### Clustering accuracy

clus_acc = function(U){
  set.seed(0)
  XU = X %*% U
  df = data.frame(XU, gen =  nutrimouse$genotype, diet = nutrimouse$diet)
  test = createDataPartition(y = df$diet, p = .2, list = FALSE)
  LDA = train(diet ~ . -gen, data = df[-test,], method = "lda",
                     trControl = trainControl(method = "cv", number = 5),
                     preProcess = c("center", "scale"),
                     tuneLength = 5)
  predictions = predict(LDA, newdata = df[test,])
  ca1 = confusionMatrix(predictions, df[test,]$diet)$overall["Accuracy"]
  
  LDA = train(gen ~ . -diet, data = df[-test,], method = "lda",
              trControl = trainControl(method = "cv", number = 5),
              preProcess = c("center", "scale"),
              tuneLength = 5)
  predictions = predict(LDA, newdata = df[test,])
  ca2 = confusionMatrix(predictions, df[test,]$gen)$overall["Accuracy"]
  
  model = Mclust(XU, G = 5:5)
  CM = confusionMatrix(factor(model$classification), factor(as.numeric(df$diet)))$table
  ca3 =  cm_accuracy(CM) 
  
  model = Mclust(XU, G = 2:2)
  CM = confusionMatrix(factor(model$classification), factor(as.numeric(df$gen)))$table
  ca4 =  cm_accuracy(CM) 
  
  ca = data.frame(lda_diet = ca1, lda_gen = ca2, gauss_diet = ca3, gauss_gen = ca4)
  rownames(ca) = NULL
  return(ca)
}

RRR = CCA_rrr(X, Y,  
              lambda = lambda_opt, 
              Kx = NULL, r=r,
              rho=1, niter=2 * 1e4,
              do.scale = FALSE, lambda_Kx=0,
              thresh=1e-6,
              solver= "ADMM",
              LW_Sy = LW_Sy)
clustering_accuracy = data.frame(method = "sparse-CCA-RRR", clus_acc(RRR$U))


for(method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                 "SCCA_Parkhomenko")){
  print(paste0("Starting ", method))
  tryCatch({
    SCCA <- additional_checks(X, Y,
                              S=NULL, 
                              rank=r, kfolds=8, 
                              method.type = method,
                              lambdax = lambdas_comp,
                              lambday =  0)
    clustering_accuracy <- rbind(clustering_accuracy, data.frame(method, clus_acc(SCCA$u)))
  }, error = function(e) {
    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
  })
}

saveRDS(clustering_accuracy, "nutrimouse_clustering_accuracy_corrected.rds")

clustering_accuracy=readRDS("Fits/nutrimouse_clustering_accuracy_corrected.rds")
clustering_accuracy
