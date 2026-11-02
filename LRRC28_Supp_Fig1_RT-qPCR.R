library(ggplot2)
library(ggpattern)
library(dplyr)
library(cowplot)

huh_pcr_data<-read.csv("Huh_pcr_result_4226_date_ordered.csv")  
#Experiment group 1 is 3.4 and experiment group 2 is 3.13 

# calculation of sd mean and sem from raw data
huh_pcr_stat <-huh_pcr_data %>%
  group_by(Experiment, Sample, Target) %>% # define the grouping order
  summarise(
    mean_rel = mean(Cq, na.rm = TRUE),                 
    sd_rel   = sd(Cq, na.rm = TRUE),                   
    n        = n(),         
    sem_rel  = sd_rel / sqrt(n),                              
    
    .groups = "drop"
  )


huh_summary_pcr <- huh_pcr_data %>%
  group_by(Experiment, Sample, Target) %>%
  summarise(mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  group_by(Experiment, Sample) %>%
  mutate(delta_Cq = mean_Cq - mean_Cq[Target == "EEF2"]) %>%
  filter(Target != "EEF2")  # filter EEF2
  

# 3. Calculating Delta Delta Cq and Rel. Expression 
huh_summary_pcr <- huh_summary_pcr %>%
  group_by(Experiment, Target) %>%
  mutate(
    delta_delta_Cq = delta_Cq - delta_Cq[grepl("Huh-7.5", Sample)],
    rel_conc = 2^(-delta_delta_Cq)
  )


# 4. sd calculation for the graph

huh_sd_error<- huh_pcr_stat %>%
  group_by(Experiment, Sample) %>%
  mutate (sd_pcr =sqrt( sd_rel^2 + sd_rel[Target == "EEF2"]^2))%>%
  filter(Target != "EEF2")

# 5. combinig the tables
huh_graph <- inner_join(huh_summary_pcr, huh_sd_error, by = c("Experiment", "Sample", "Target"))

huh_pcr_plot_2<-huh_graph %>%
  filter(rel_conc>10) %>%
  ggplot(aes(x = Sample, y = rel_conc)) +
  geom_point(color = "blue", size = 2) + 
 
  geom_errorbar(
    aes( ymin = 2^(-(delta_delta_Cq + sd_pcr)), ymax = 2^(-(delta_delta_Cq - sd_pcr))),
    width = 0.1,
    color = "black",
    linewidth = 0.7
  ) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  theme_bw()+
  labs(
    title = "Huh-7.5",
    x = "Sample",
    y = "Relative Expression to EEF2"
  ) +
  
  theme(  
    axis.text.x  = element_text(size = 13),     
    axis.text.y  = element_text(size = 13),      
    axis.title   = element_text(size = 13),     
    legend.position = "bottom",
    axis.title.x = element_text(
      angle = 0,     
      vjust =-1,    
      hjust = 0.5 ),    
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )
 
  huh_pcr_plot_2
shsy5y_pcr_data<-read.csv("shsy5y_pcr_data.csv")

#experiment 1 is dated back to 4.10.26 and experiment 2 is 4.16.26
shsy5y_pcr_stat <-shsy5y_pcr_data %>%
    group_by(Experiment, Sample, Target) %>% 
    summarise(
      mean_rel = mean(Cq, na.rm = TRUE),                   
      sd_rel   = sd(Cq, na.rm = TRUE),                   
      n        = n(),                                           
      sem_rel  = sd_rel / sqrt(n),                            
      
      .groups = "drop"
    )
  
  
shsy5y_summary_pcr <- shsy5y_pcr_data %>%
    group_by(Experiment, Sample, Target) %>%
    summarise(mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
    group_by(Experiment, Sample) %>%
    mutate(delta_Cq = mean_Cq - mean_Cq[Target == "EEF2"]) %>%
    filter(Target != "EEF2")

  
# 3. Calculating Delta Delta Cq and Rel. Expression 
shsy5y_summary_pcr <- shsy5y_summary_pcr %>%
    group_by(Experiment, Target) %>%
    mutate(
      delta_delta_Cq = delta_Cq - mean(delta_Cq[grepl("SHSY5Y", Sample)]),
      rel_conc = 2^(-delta_delta_Cq)
    )
  
  # 4. sd calculation for the graph
  
shsy5y_pcr_sd<- shsy5y_pcr_stat %>%
    group_by(Experiment, Sample) %>%
    mutate (sh_sd_pcr =sqrt( sd_rel^2 + sd_rel[Target == "EEF2"]^2))%>%
    filter(Target != "EEF2")
  # 5. combine the tables
shsy5y_graph <- inner_join(shsy5y_summary_pcr, shsy5y_pcr_sd, by = c("Experiment", "Sample", "Target"))
  

  
shsy5y_pcr_plot_2<-shsy5y_graph %>%
  filter(rel_conc>10) %>%
  filter(Sample !="28 MOI10 3.9") %>%
  ggplot(aes(x = Sample, y = rel_conc)) +
  geom_point(color = "blue", size = 2) + 

  geom_errorbar(
    aes(ymin = 2^(-(delta_delta_Cq + sh_sd_pcr)), ymax = 2^(-(delta_delta_Cq - sh_sd_pcr))),
    width = 0.1,
    color = "black",
    linewidth = 0.7
  ) +
  facet_grid(. ~ Experiment, scales = "free_x", space = "free_x") +
  
  theme_bw() +
  labs(
    title = "SH-SY5Y",
    x = "Sample",
    y = "Relative Expression to EEF2"
  ) +
  theme(
    axis.text.x  = element_text(size = 17),     
    axis.text.y  = element_text(size = 13),      
    axis.title   = element_text(size = 13),      
    strip.text   = element_text(size = 13, face = "bold"), 
    legend.position = "bottom",
    axis.title.x = element_text(
      angle = 0,    
      vjust =-1,      
      hjust = 0.5 ),  
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )
  shsy5y_pcr_plot
  
  

pcr_final_plot_2 <- plot_grid(huh_pcr_plot_2, shsy5y_pcr_plot_2,
                            nrow = 2,
                        labels = "AUTO")

ggsave("RTPCR_Figures_merged_dot.pdf", pcr_final_plot_2, width =16, height = 13)
