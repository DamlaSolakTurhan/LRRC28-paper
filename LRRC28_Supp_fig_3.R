library(readxl)
library(ggplot2)
HuhCompiled= read_xlsx("huh7.5_ocr_graph_data.xlsx")

#I obtained this data from Seahorse Wavepro desktop application by clicking on the graph. 
#For Huh-7.5 Replicate 1 I exported base1 corrected graph's data and for Replicate 2(dated:3.13) I exported base2 corrected graph's data.
class(HuhCompiled)
head(HuhCompiled)
dim(HuhCompiled)
HuhCompiled$Gene <- factor(HuhCompiled$Gene, levels = c("LRRC28", "LRRC42", "Huh-7.5"))
levels(HuhCompiled$Gene)
HuhCompiled$Gene

Huh7_5_34=ggplot(HuhCompiled[1:27,], aes(x = `Time (min)`, y = Value, color= Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% OCR", title = "Huh-7.5 Experiment 1")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth =  1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+
  theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
        legend.title = element_text(size = 10), 
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c( "Huh-7.5"= "#56B4E9", "LRRC28" = "#D55E00" , "LRRC42" = "#EEB479"))


Huh7_5_34

Huh_7_5_313=ggplot(HuhCompiled[28:54,], aes(x = `Time (min)`, y = Value, color= Gene))+
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0,200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% OCR", title = "Huh-7.5 Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+ theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                   legend.title = element_text(size = 10), 
                                                   legend.background = element_rect(fill = "transparent", color = NA),
                                                   legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c( "Huh-7.5"= "#56B4E9","LRRC28" = "#D55E00", "LRRC42" = "#EEB479" ))


Huh_7_5_313

SHSY5Y_Compiled= read_xlsx("shsy5y_ocr_graph_data.xlsx") 
#I obtained this data from Seahorse wavepro desktop application by clicking on the graph.
#SHSY5Y Replicate 2 (4.17) base1 and Replicate 1(4.10) base 2 was used for writing and figure generation
class(SHSY5Y_Compiled)
head(SHSY5Y_Compiled)
dim(SHSY5Y_Compiled)
unique(SHSY5Y_Compiled$Gene)
SHSY5Y_Compiled$Gene <- as.factor(SHSY5Y_Compiled$Gene)
#figures

SHSY5Y_4.10_MOI1=ggplot(SHSY5Y_Compiled[c(1:18,28:36),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2)   +             # Error bars
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 1")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 10), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI<1" = "#D55E00", "LRRC42 MOI<1"= "#EEB479"))




SHSY5Y_4.10_MOI1


SHSY5Y_4.17_MOI1=ggplot(SHSY5Y_Compiled[c(64:81,100:108),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0,1000)) +  
  scale_y_continuous(breaks = seq(0,1000, by = 200)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +# Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 10), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI<1" = "#D55E00", "LRRC42 MOI<1"= "#EEB479"))

SHSY5Y_4.17_MOI1


SHSY5Y_4.10_MOI1_NTD=ggplot(SHSY5Y_Compiled[c(1:18,46:54),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 1")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI<1" = "#D55E00", "LRRC NTD MOI<1"="#CC79A7"))


SHSY5Y_4.10_MOI10=ggplot(SHSY5Y_Compiled[c(1:9,19:27,37:45),], aes(x = `Time (min)`, y = Value, color= Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 1")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI10" = "#D55E00", "LRRC42 MOI10"= "#EEB479"))


SHSY5Y_4.10_MOI10_NTD=ggplot(SHSY5Y_Compiled[c(1:9,19:27,55:63),], aes(x = `Time (min)`, y = Value, color= Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 1")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI10" = "#D55E00", "LRRC NTD MOI10"= "#CC79A7"))



SHSY5Y_4.17_MOI1_NTD=ggplot(SHSY5Y_Compiled[c(64:81,118:126),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 1000)) +  
  scale_y_continuous(breaks = seq(0, 1000, by = 200)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   # Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI<1" = "#D55E00","LRRC NTD MOI<1"="#CC79A7"))

SHSY5Y_4.17_MOI10=ggplot(SHSY5Y_Compiled[c(64:72,82:90,109:117),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   # Error bars
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI10_4.3" = "#D55E00", "LRRC42 MOI10"= "#EEB479"))



SHSY5Y_4.17_MOI10_NTD=ggplot(SHSY5Y_Compiled[c(64:72,82:90,127:135),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +# Error bars
  coord_cartesian(ylim = c(0, 200)) +  
  scale_y_continuous(breaks = seq(0, 200, by = 20)) +   
  labs(x = "Time", y = "% OCR", title = "SH-SY5Y Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI10_4.3" = "#D55E00", "LRRC NTD MOI10"="#CC79A7"))

final_plot_sup <- plot_grid(Huh7_5_34,Huh_7_5_313, 
                            SHSY5Y_4.10_MOI1,SHSY5Y_4.17_MOI1,
                            SHSY5Y_4.10_MOI10, SHSY5Y_4.17_MOI10, 
                            SHSY5Y_4.10_MOI1_NTD, SHSY5Y_4.17_MOI1_NTD,
                            SHSY5Y_4.10_MOI10_NTD, SHSY5Y_4.17_MOI10_NTD,
                            labels = c("A","","B","","C","","D","","E",""),
                            nrow = 5)

ggsave("Supplementary_Fig3.pdf", final_plot_sup,  device = cairo_pdf, width = 15, height = 20)



