library(reshape2)
library(ggplot2)
library(plotly)

load("Fits/processed_data.RData")

ynames = colnames(Y) 
ylabels = c("Drive", "Funseeking", "Reward Response", "Total", "Distress", "Anhedonia", "Anxious Arousal", "Positive Affect", "Negative Affect")

gg_color_hue = function(n) {
  set.seed(1)
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gyri = sort(unique(positions$Gyrus))
colors0 = gg_color_hue(length(gyri))
topgyrilist= list()

filenames = c("RRR_res", "RRR_res_graph", "RRR_res_group", "SAR_res")
methods = c("rrr-sparse", "rrr-graph", "rrr-group", "sar")

for(k in 1:length(filenames)){
  method = methods[k]
  filename = filenames[k]
  load(paste0("Fits/", filename, ".RData"))
  
  V = Vhat_comp
  colnames(V) = paste0("V", 1:3)
  dfV = data.frame(V, question = colnames(Y)) %>%
    melt(id = c("question")) %>%
  mutate(variable = factor(variable,
                           levels = paste0("V", 1:3),
                           labels = c(expression(v[1]), expression(v[2]), expression(v[3]))),
         question = factor(question, 
                           levels = ynames, 
                           labels = ylabels))
  
  ggplot(dfV) +
    geom_bar(aes(value, question, fill = variable), stat = "identity", show.legend = FALSE)+
    facet_grid(cols = vars(variable), labeller = label_parsed) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(y = "Question", x = "Value")+
    theme(axis.text.x = element_text(size = 4))+
    scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5))+
    theme_bw()
  ggsave(paste0("brain-", method, "-rrr-v.pdf"), width = 4.5, height = 3)
  
  
  U = Uhat_comp
  colnames(U) = paste0("U", 1:3)
  dfU = data.frame(U, region = positions$Region, gyrus = positions$Gyrus) %>%
    mutate(U = abs(U1) + abs(U2) + abs(U3)) %>%
    mutate(gyrus = factor(gyrus, levels = gyri))
  
  dfUtop = dfU %>%
    slice_max(U, n = 20)
  eps = dfUtop %>% pull(U) %>% min()
  topgyri = unique(dfUtop$gyrus)
  topgyrilist[[method]] = topgyri
  colors = colors0[((1:30)*7) %% 30 + 1]
  colors[!gyri %in% topgyri] = "grey"
  
  dfU1 = rbind(dfUtop,
        dfU %>% filter(!gyrus %in% topgyri) %>% distinct(gyrus, .keep_all = TRUE) %>%
          mutate(U1 = 0, U2 = 0, U3 = 0)) %>%
    dplyr::select(-U) %>%
    mutate(gyrus = factor(gyrus, levels = gyri)) %>%
    arrange(gyrus)
  
  ma = dfU1 %>% select(U1, U2, U3) %>% max()
  mi = dfU1 %>% select(U1, U2, U3) %>% min()
  
  dfU2 = dfU1 %>%
    melt(id = c("region", "gyrus")) %>%
    mutate(variable = factor(variable,
                             levels = paste0("U", 1:3),
                             labels = c(expression(u[1]), expression(u[2]), expression(u[3]))),
           label = ifelse(abs(value) < 0.1, NA, region),
           region = factor(region, levels = dfU1$region),
           ticks = as.numeric(region),
           gyrus = factor(gyrus, levels = gyri)) 
  
  ticks = dfU2 %>% count(gyrus) %>% mutate(n = n/3, c = cumsum(n)) %>% arrange(gyrus)
  ticks = ticks %>% mutate(m = (c + c(0, ticks$c[-30]))/2)
  
  ggplot(dfU2) +
    geom_bar(aes(value, factor(ticks), fill = gyrus), stat = "identity", show.legend = FALSE)+
    facet_grid(cols = vars(variable), labeller = label_parsed) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = ticks$c+0.5, color = "grey", alpha = 0.3) + 
    geom_text(aes(-0.015 * sign(value), ticks, label = label, color = gyrus, hjust = 0.5 * (1+sign(value))), size = 1.6, color = "black")+
    labs(y = "Region", x = "Value")+
    scale_y_discrete(breaks = ceiling(ticks$m), labels = ticks$gyrus)+
    #xlim(mi-0.1, ma+0.1)+
    scale_fill_manual(values = colors)+
    #scale_color_manual(values = colors)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_text(colour = colors),
          axis.ticks.y = element_blank())
  ggsave(paste0("brain-", method, "-rrr-u.pdf"), width = 8, height = 5)
  
  
  XU = X %*% U
  colnames(XU) = c("Xu1", "Xu2", "Xu3")
  df = data.frame(XU)
  
  # for(i in 1:length((ynames))){
  #   plt = plot_ly(data = df, x = ~Xu1, y = ~Xu2, z = ~Xu3,
  #                color= Y[,i],
  #                colors = 'RdBu',
  #                type = 'scatter3d', mode = 'markers') %>% 
  #     colorbar(title = ylabels[i])
  #   htmlwidgets::saveWidget(as_widget(plt), paste0("brain-", method, "-rrr-xu-", ynames[i],".html"))
  #   #Sys.setenv("plotly_username" = "elenatuz")
  #   #Sys.setenv("plotly_api_key" = "6S7e4oo2erHqjHcTT5d6")
  #   #plotly_IMAGE(plt, format = "png", out_file = paste0("brain-", method, "-rrr-xu-", ynames[i],".png"), width = 650,
  #   #             height = 500)
  # }
  
  plt = data.frame(positions, u1 = U[,1]) %>%
    mutate(label = ifelse(abs(u1) < 0.2, NA, Region)) %>%
    plot_ly(x = ~x, y = ~y, z = ~z,
            color = ~u1,
            colors = 'Spectral',
            type = 'scatter3d', 
            mode = 'markers') %>%
    layout(scene = list(xaxis = list(title = "", showticklabels = FALSE), yaxis = list(title = "", showticklabels = FALSE), zaxis = list(title = "", showticklabels = FALSE)))
  htmlwidgets::saveWidget(as_widget(plt), paste0("brain-", method, "-rrr-u1.html"))
}


setdiff(topgyrilist[[1]],
topgyrilist[[4]])
intersect(topgyrilist[[1]], topgyrilist[[2]])

topgyrilist[[3]]
topgyrilist[[4]]


intersect(
  intersect(intersect(topgyrilist[[1]],
             topgyrilist[[2]]
),topgyrilist[[3]]), topgyrilist[[4]])

intersect(intersect(topgyrilist[[1]],
                      topgyrilist[[2]]
  ),topgyrilist[[3]])

