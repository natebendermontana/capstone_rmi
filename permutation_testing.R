
# Permutation testing for Capstone work with RMI

options(scipen=999) # force full notation not scientific
library(here)
library(ggplot2)
library(patchwork)
library(tidyverse)


get.prop.diff <- function(x){
  in_buffer <- sum(x$MINORPOP_intersect_count) / sum(x$ACSTOTPOP_intersect_count)
  out_buffer <- sum(x$out_mnr_count) / sum(x$out_pop_count)
  return(in_buffer - out_buffer)
}

permutation_test <- function(path, buffer_size){
  df <- read_tsv(path)
  
  # df <- na.omit(df[, c("MINORPOP_intersect_count", "ACSTOTPOP_intersect_count",
  #                      "out_mnr_count", "out_pop_count")])
  
  actual.value <- get.prop.diff(df) 
  
  n.sim <- 2000
  
  results <- tibble(statistic = c(actual.value,
                                  rep(NA,n.sim)))
  new.df <- df %>% 
    select(ACSTOTPOP, MINORPOP,
           MINORPOP_intersect_count, ACSTOTPOP_intersect_count,
           out_mnr_count, out_pop_count, intersect_prop)
  
  set.seed(20230331)
  
  for(i in 2:nrow(results)){
    new.df$intersect_prop <- sample(new.df$intersect_prop)
    
    # create new counts variables based on the sampled proportion
    new.df$ACSTOTPOP_intersect_count <- new.df$intersect_prop*new.df$ACSTOTPOP
    new.df$MINORPOP_intersect_count <- new.df$intersect_prop*new.df$MINORPOP
    new.df$out_pop_count <- new.df$ACSTOTPOP - new.df$ACSTOTPOP_intersect_count
    new.df$out_mnr_count <- new.df$MINORPOP - new.df$MINORPOP_intersect_count
    
    # Calc result statistic
    results$statistic[i] <- get.prop.diff(new.df)
  }
 
  return(list(results=results,actual.value=actual.value)) 
}

# permutation_test_plot <- function(data,buffer_size){
# 
#   actual.value <- data[[buffer_size]]$actual.value
#   data <- data[[buffer_size]]$results
#   
#   max_density <- max(data$statistic, na.rm = TRUE)
#   pval <- mean(data$statistic > actual.value)
#   pval_label <- paste0("p-value:\n", round(pval, 3))
#   act_label <- paste0("act. value:\n", round(actual.value, 3))
#   
# 
#   plot <- ggplot(data, 
#                  aes(x=statistic)) + 
#     geom_density() + 
#     geom_vline(xintercept=actual.value,color="red") +
#     #annotate("text", x=-.1, y=35, label=pval_label, color="red", size=3) + # add the label
#     #annotate("text", x=.12, y=35, label=act_label, color="red", size=3) + # add the label
#     #ylim(0, 50) +  # set y limits
#     xlim(-.3,.3)+
#     theme_minimal() + 
#     labs(x=NULL, y=paste0(buffer_size,"m"))  # remove x label and add y label
#   
#   return(plot)
# }  


# Combine the individual permutation tests
permutation_test_multi <- function(paths, buffer_sizes){
  plots <- list()
  for (i in seq_along(paths)){
    plots[[i]] <- permutation_test(paths[i], buffer_sizes[i])
  }
  # combine the plots using patchwork
  plot_combined <- wrap_plots(plots) + 
    plot_layout(ncol = 1) +
    plot_annotation(title = "Density Plots for Multiple Flare Buffers")
  return(plot_combined)
}

# Set up buffer sizes and file paths
buffer_sizes <- c(100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000)
# single size for permutation viz, if needed
#buffer_sizes <- 2000

file_names <- paste0("adjacent_df_",buffer_sizes,"m_forpermutation.tsv")
paths <- here("data",file_names)

buffer_sizes
paths

# Visualize density plot of single buffer size, if needed
permutation_test_plot2 <- function(data,buffer_size){
  
  actual.value <- data$actual.value
  data <- data$results
  
  #max_density <- max(data$statistic, na.rm = TRUE)
  #mean(abs(results$statistic) >= actual.value
  pval <- mean(abs(data$statistic) >= actual.value)
  pval_label <- paste0("p-value:\n", round(pval, 3))
  act_label <- paste0("act. value:\n", round(actual.value, 3))
  
  
  plot <- ggplot(data, 
                 aes(x=statistic)) + 
    geom_density() + 
    geom_vline(xintercept=actual.value,color="red") +
    #annotate("text", x=-.1, y=35, label=pval_label, color="red", size=3) + # add the label
    #annotate("text", x=.12, y=35, label=act_label, color="red", size=3) + # add the label
    #ylim(0, 50) +  # set y limits
    scale_x_continuous(labels=scales::percent,
                       breaks = seq(-0.3, 0.3, 0.1),
                       limits = c(-0.3, 0.3)) +
    #scale_y_discrete(labels = function(x) scales::comma(as.numeric(x))) +
    labs(y = "Density\n",
         x = "\n(In Buffer Minority %) - (Out Buffer Minority %)") +
    xlim(-.12,.12)+
    theme_minimal(base_size = 14) +
    # geom_label(aes(x = -0.2, y = 17.5, label = "Whites More Affected"), 
    #            size = 4, fill = "white", color = "black", label.padding = unit(0.3, "lines"),
    #            label.size = 0.15, hjust = 0) +
    # geom_label(aes(x = 0.15, y = 17.5, label = "BIPOC More Affected"), 
    #            size = 4, fill = "white", color = "black", label.padding = unit(0.3, "lines"),
    #            label.size = 0.15, hjust = 0)
   # expand_limits(y = c(0, 9.96)) # increase the y-axis limit to include the label elements
    #labs(x=NULL, y=paste0(buffer_size,"m"))  # remove x label and add y label
  
  print(pval)
  return(plot)
}   
test <- permutation_test(paths, buffer_sizes)
permutation_test_plot2(test, buffer_sizes)

test$actual.value


#####
# Create CI plot centered on actual values
results <- vector("list", length(buffer_sizes))

for(this_buffer in buffer_sizes){
#  if(exists("results")){
#    rm(results)
#  } 
  
  results[[this_buffer]] <- permutation_test(paths[buffer_sizes==this_buffer],this_buffer)
    
}

#TESTING
#results[100]
#test1 <- permutation_test(paths, buffer_sizes)
#test1

# Set up parameters and table for CI plot
ci_width <- 0.90
probs <- c((1-ci_width)/2,ci_width + (1-ci_width)/2)

for_plot <- tibble(
  buffer_size = buffer_sizes,
  average_effect = 0.0,
  lb = 0.0,
  ub = 0.0
)

for(i in 1: nrow(for_plot)){
  this_buf <- for_plot$buffer_size[i]
  for_plot$average_effect[i] <- results[[this_buf]]$actual.value
  for_plot$lb[i] <- quantile(results[[this_buf]]$results %>% pull(statistic),probs[1])
  for_plot$ub[i] <- quantile(results[[this_buf]]$results %>% pull(statistic),probs[2])
}

# Shift our CIs to be centered on the actual value
for_plot <- for_plot %>% 
  mutate(lb = lb + average_effect,
         ub = ub + average_effect,
         buff_fct = as_factor(buffer_size))

# Orig code - version with buffer size as a factor
# ggplot(for_plot,
#        aes(y=buff_fct,x=average_effect)) + 
#   geom_point() + 
#   geom_errorbarh(aes(xmin=lb,xmax=ub,y=buff_fct),height=0.1) + 
#   geom_vline(xintercept=0.0) + 
#   theme_minimal() +
#   
#   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
#   labs(y="Buffer Size (m)",
#        x="(In Buffer Minority %) - (Out Buffer Minority %)") +
#   geom_text(aes(x = -0.1,y=10,
#            label = "Whites More Affected"),
#            size=3) +
#   geom_text(aes(x = 0.1, y=10,
#                label = "BIPOC More Affected"),
#            size=3)

ggplot(for_plot, aes(y = buff_fct, x = average_effect)) +
  geom_point() + 
  geom_errorbarh(aes(xmin = lb, xmax = ub, y = buff_fct), height = 0.1) + 
  geom_vline(xintercept = 0.0) + 
  theme_minimal(base_size = 14) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  scale_x_continuous(labels=scales::percent) +
#  scale_y_discrete(labels=scales::comma) +
  scale_y_discrete(labels = function(x) scales::comma(as.numeric(x))) +
  labs(y = "Buffer Radius (m)\n",
       x = "\n(In Buffer Minority %) - (Out Buffer Minority %)")
  # geom_label(aes(x = -0.07, y = 9.8, label = "Whites More Affected"), 
  #            size = 3, fill = "white", color = "black", label.padding = unit(0.3, "lines"),
  #            label.size = 0.15, hjust = 1) +
  # geom_label(aes(x = 0.05, y = 9.8, label = "BIPOC More Affected"), 
  #            size = 3, fill = "white", color = "black", label.padding = unit(0.3, "lines"),
  #            label.size = 0.15, hjust = 0) +
  # expand_limits(y = c(0, 9.96)) # increase the y-axis limit to include the label elements





# version with buffer size as continuous
ggplot(for_plot,
       aes(y=buffer_size,x=average_effect)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin=lb,xmax=ub,y=buffer_size),height=0) + 
  geom_vline(xintercept=0.0) + 
  theme_minimal() + 
  labs(y="Buffer Size (m)",
       x="(In Buffer Minority %) - (Out Buffer Minority %)") +
  geom_text(aes(x = -0.15,y=4500,
                label = "Whites More Affected"),
            size=3) + 
  geom_text(aes(x = 0.1, y=4500,
                label = "BIPOC More Affected"),
            size=3) 


