library(data.table)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library("magick")
library(readxl)
library(celda) # for decontX, it is included in celda


#rna_tissue_gtex.tsv.zip was dowloanded and renamed as rna_tissue_hpa_normalized_gtex.tsv
rna_gtex_normalized <- fread("rna_tissue_hpa_normalized_gtex.tsv", select = c("Gene name", "Tissue", "nTPM"))
LRRC28_gtex_data <- rna_gtex_normalized[`Gene name` == "LRRC28"]
class(LRRC28_gtex_data$nTPM)
str(LRRC28_gtex_data) 
LRRC28_gtex_data$Tissue <- fct_reorder(LRRC28_gtex_data$Tissue, LRRC28_gtex_data$nTPM, .fun = median)

LRRC28_Tissue_gtx_normalized_plot<- ggplot(LRRC28_gtex_data, aes(x = Tissue, y = nTPM)) +
  geom_col(fill = "steelblue") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        plot.title = element_text(hjust = 0.5, color = "black", size= 14 ),
        axis.text.y = element_text (size= 14, color= "black"),
        axis.title = element_text (size = 14, color = "black")) +
        labs(title = "Expression Profiles for LRRC28", x = "Tissue", y = "nTPM from GTEX Data")+
  coord_flip()
LRRC28_Tissue_gtx_normalized_plot



#18-34 young  35-49 middle 50+ old
samples_to_use <- list(
  Donor4  = list(gsm = "GSM8037606", age = 30, group = "young"),
  Donor5 = list(gsm = "GSM8037607", age = 32, group = "young"),
  Donor11 = list(gsm = "GSM8037608", age = 43, group = "middle"),
  Donor21 = list(gsm = "GSM8037609", age = 55, group = "old"),
  Donor7 = list(gsm = "GSM8037610", age = 38 , group = "middle"),
  Donor12 = list(gsm = "GSM8037611", age = 43, group = "middle"),
  Donor18 = list(gsm = "GSM8037612", age = 52, group = "old"),
  Donor20 = list(gsm = "GSM8037613", age = 55, group = "old"),
  Donor24 = list(gsm = "GSM8037614", age = 59, group = "old"),
  Donor16 = list(gsm = "GSM8037615", age = 46, group = "middle"),
  Donor17 = list(gsm = "GSM8037616", age = 48, group = "middle"),
  Donor22 = list(gsm = "GSM8037617", age =57 , group = "old"),
  Donor23 = list(gsm = "GSM8037618", age = 57, group = "old"),
  Donor15 = list(gsm = "GSM8037619", age = 44, group = "middle"),
  Donor25 = list(gsm = "GSM8037620", age = 59 , group = "old"),
  Donor6 = list(gsm = "GSM8037621", age = 33, group = "young"),
  Donor8 = list(gsm = "GSM8037622", age = 39 , group = "middle"),
  Donor9 = list(gsm = "GSM8037623", age = 42, group = "middle"),
  Donor10 = list(gsm = "GSM8037624", age = 43 , group = "middle"),
  Donor13 = list(gsm = "GSM8037625", age = 43, group = "middle"),
  Donor14 = list(gsm = "GSM8037626", age = 44 , group = "middle"),
  Donor19 = list(gsm = "GSM8037627", age = 53 , group = "old"),
  Donor35 = list(gsm = "GSM8037628", age = 69 , group = "old"),
  Donor30 = list(gsm = "GSM8037629", age = 66 , group = "old"),
  Donor34 = list(gsm = "GSM8037630", age = 67, group = "old")
  
)

data_dir <- "C:/Users/damla/Box/scRNA Seq/donors data"

seurat_list <- list()
for (donor in names(samples_to_use)) {
  gsm <- samples_to_use[[donor]]$gsm
  raw_counts <- Read10X(data.dir = file.path(data_dir, gsm))
  
  # Run DecontX to remove ambient RNA
  decont_obj <- decontX(x = raw_counts)
  
  # Create Seurat object from cleaned counts
  # filtering for min.cells and min.features
  obj <- CreateSeuratObject(counts = decont_obj$decontXcounts, 
                            project = donor, 
                            min.cells = 3, 
                            min.features = 200)
  
  # Add metadata
  obj$donor <- donor
  obj$gsm   <- gsm
  obj$age <- samples_to_use[[donor]]$age
  obj$group <- samples_to_use[[donor]]$group
  
  # Calculate Mitochondrial percentage since we will apply filter
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  

  obj <- subset(obj, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < 5500 & 
                  percent.mt < 15 & 
                  nCount_RNA > 800 & 
                  nCount_RNA < 17000)
  
  # Store in list
  seurat_list[[donor]] <- obj
}
ncol(merged) #72553
nrow(merged) #32572
table(merged@meta.data$donor)
#the average UMI for each donors
aggregate(nCount_RNA ~ donor, data = merged@meta.data, FUN = mean)
#average UMI for all donors
mean(merged$nCount_RNA)
mean(merged$nFeature_RNA)

# Merge into one object, keeping each sample as a separate layer.  
#add.cell.id function is adding prefix to the ids 
merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
merged 
table(merged$donor) #check the number of cells for each donor
head(colnames(merged)) #check whether cell prefixes are added
table(merged$donor, merged$group)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 5000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
ElbowPlot(merged, ndims=50)
#hormony integration
merged <- IntegrateLayers(
  object = merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

Layers(merged) #there are layers for each donor, we will join them together for further analysis
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

Layers(merged)
results <- list()
cluster_counts <- list()
dimplots <- list()   

for (dims in c(15,20)) {
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:dims)
  merged <- FindClusters(merged, resolution = 0.5)
  merged <- RunUMAP(merged, reduction = "harmony", dims = 1:dims)
  
  n_clusters <- length(unique(Idents(merged)))
  cluster_counts[[paste0("dims_", dims)]] <- n_clusters
  message("dims = ", dims, " → ", n_clusters, " clusters")
  
  # Labeled UMAP
  dimplots[[paste0("dims_", dims)]] <- DimPlot(merged, 
                                               label = TRUE, 
                                               repel = TRUE) +
    ggtitle(paste0("dims 1:", dims, " | ", n_clusters, " clusters"))
  
  # Feature plots
  results[[paste0("dims_", dims)]] <- FeaturePlot(merged,
                                                  features = c("TNP1", "VIM", "DDX4", "CD68"),
                                                  ncol = 2)
}

# View cluster counts
print(cluster_counts)

# View each result
results[["dims_15"]]
results[["dims_16"]]
results[["dims_17"]]
results[["dims_18"]]
results[["dims_19"]]
results[["dims_20"]]

# dim 1:18 was chosen
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:18)
merged <- FindClusters(merged, resolution = 0.35)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:18)


all_markers<- FindAllMarkers(object = merged, only.pos = TRUE,  min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox",    # default and fastest
                              max.cells.per.ident = 500  )


top_markers <- all_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = pct_diff)


canonical_markers_2 <- c(
  "DDX4",  #GERMLINE SPECIFIC MARKER
  "UTF1","ID4", "PIWIL4",   #SSC MARKER
  "KIT","MKI67","HMGN2", # Differentiating SPC MARKERS
  "SYCP3", "SPO11",   #upregulated from the leptotene to the zygotene stage
  "MEIOB","TOP2A", # early spermatocyte
  "CETN3" , "SPAG6",   # late primary spermatocyte
  "ZPBP", "ACRV1","SPACA1",
  "CCDC168" ,  # SPERMATID SPECIFIC GENES (S1-S2) # round spermatid
  "TEX29", #early stage of spermatids
  "NFKBIB","IQCF3","LELP1",  #s1-s4
  "TSSK6", "SPATA3", #elongating spermatid
  "SPEM1" , #late-stage spermatid markers
  "TNP1", "TNP2", #S2-S4 transition into PRM1 and PRM2
  "PRM2", #S2-S4 transition into PRM1 and PRM2
  "SOX9", "AMH",          # Sertoli
  "TSPAN33",  "WFDC1", "PTCH1","ACTA2", "MYH11",  # TPC
  "NOTCH3", "RGS5", #SMC
  "GJA4", #PERICYTES IN TESTIS HPA DATA
  "VWF", "CD34","WT1",  # endothelial
  "FCN3", "CLEC14A", # Capillary
  "ACKR1", #Vein
  "PODXL", "KDR", #vascular EC
  "DLK1", "IGF1",  "NEGR1", # Leydig
  "C1QC", "CD14",    # macrophage
  "CD3D", "CD69",     # lymphocyte
  "CD79A", "IGHM",     # B cell
  "LRRC28","LRRC42")

DotPlot(merged, features = canonical_markers_2 ,  cols = "PiYG" , scale.by = "size", dot.scale = 8) +
  RotatedAxis()


levels(merged)

my_order <- c("12","10","4","6","9","2","0","7","13","1","3","5","8","11")
# Apply the new factor levels to your Seurat object's active identity
Idents(merged) <- factor(Idents(merged), levels = my_order)

levels(merged)

# saveRDS(merged, "DecontX_before_merged_after_cluster_reorder.rds")
#merged<- readRDS("DecontX_before_merged_after_cluster_reorder.rds")


new.cluster.ids <- c("Undifferentiated Spermatogonia",
                     "Early Primary Spermatocytes" ,
                     "Late Primary Spermatocytes" ,
                     "Round Spermatids" ,
                     "Early Elongating Spermatids" ,
                     "Late Elongating Spermatids",
                     "Testicular Peritubular Cells",
                     "Smooth Muscle",
                     "Pericytes",
                     "Vascular Endothelial Cells - Venous",
                     "Vascular Endothelial Cells",
                     "Leydig",
                     "Macrophage",
                     "Lymphocyte")
names(new.cluster.ids) <- levels(merged)
# Apply the renaming to the Seurat object
merged <- RenameIdents(merged, new.cluster.ids)
# To reset the identities to the original cluster numbers stored in meta.data #Idents(merged) <- "seurat_clusters"

levels(merged)

merged_featureplot <- FeaturePlot(merged, features = c("TNP1","DDX4", "VIM"), ncol = 1)

LRRC28_FeaturePlot <-FeaturePlot(merged, features = c("LRRC28"), pt.size = 1.2)+
  theme(axis.title = element_text(color = "black", size=14),
        axis.text = element_text (color = "black",size = 14 ),
        plot.title = element_text(color = "black", size= 16)
        )

merged_dim_plot <- DimPlot(merged, reduction = "umap", label = TRUE, label.box = FALSE) +
  NoLegend()


table(Idents(merged)) # number of cells in each group

merged_dotplot <- DotPlot(merged, features = canonical_markers_2,   cols = "PiYG", scale.by = "size",dot.scale = 6) +
  theme(axis.title = element_text(color= "black", size=14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"), 
    axis.text.y = element_text(size = 14, color="black")               
  ) +
  coord_flip()


# Quick numeric summary — low nFeature = empty droplets, high = doublets
cluster_qc_whole <- merged@meta.data %>%
  summarise(
    n_cells       = n(),
    mean_nCount = mean(nCount_RNA),
    median_nCount = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    mean_nFeature = mean(nFeature_RNA),
    median_mt     = median(percent.mt),
    .groups = "drop"
  )


print(cluster_qc_whole)



Supplementary_Fig4 <- plot_grid(merged_dim_plot, merged_featureplot,
                          labels = "AUTO",
                          ncol = 2,
                          rel_widths = c(1, 0.7))

ggsave("LRRC28_scRNA_supp_fig.pdf", Supplementary_Fig4,  device = cairo_pdf,  width = 15, height = 10)



HuhCompiled= read.csv("Huh Base Corrected Compiled Data.csv") 
#I obtained this data from Seahorse Wavepro desktop application by clicking on the grap that drawn by the app. 
#For Huh3.4 I exported base1 corrected graph's data and for 3.13 I exported base2 corrected graph's data.
class(HuhCompiled)
head(HuhCompiled)
dim(HuhCompiled)
HuhCompiled$Gene <- factor(HuhCompiled$Gene, levels = c("LRRC28", "LRRC42", "Huh-7.5"))


Huh7_5_34=ggplot(HuhCompiled[1:27,], aes(x = Time..min., y = Value, color= Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 550)) +  
  scale_y_continuous(breaks = seq(0, 550, by = 100)) +   # <-- forces labels every 100
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   
  labs(x = "Time", y = "% ECAR", title = "Huh-7.5 Experiment 1")+
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

Huh_7_5_313=ggplot(HuhCompiled[28:54,], aes(x = Time..min., y = Value, color= Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 550)) +   
  scale_y_continuous(breaks = seq(0, 550, by = 100)) +   # <-- forces labels every 100
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   
  labs(x = "Time", y = "% ECAR", title = "Huh-7.5 Experiment 2")+
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

SHSY5Y_Compiled= read_xlsx("SH-SY5Y_data.xlsx") 
#I obtained this data from Seahorse wavepro desktop application by clicking on the grap that drawn by the app.
#SHSY5Y 4.17 base1 and 4.10 base  2 was used for writing and figure generation
class(SHSY5Y_Compiled)
head(SHSY5Y_Compiled)
dim(SHSY5Y_Compiled)
unique(SHSY5Y_Compiled$Gene)
SHSY5Y_Compiled$Gene <- as.factor(SHSY5Y_Compiled$Gene)
#figures
unique(SHSY5Y_Compiled$`Time (min)`)


SHSY5Y_4.10_MOI1=ggplot(SHSY5Y_Compiled[c(1:18,28:36),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2)   +             
  coord_cartesian(ylim = c(0, 900)) +   
  scale_y_continuous(breaks = seq(0, 900, by = 200)) +   # <-- forces labels every 200
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 1")+
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
  coord_cartesian(ylim = c(0,900)) +  
  scale_y_continuous(breaks = seq(0,900, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +# Error bars
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 2")+
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


my_image <- image_read("lrrc28_model0_chimeraX.png")
p_img <- ggdraw() + draw_image(my_image, scale = 1.1) 

spacer <- ggplot() + theme_void()

right_plot <- plot_grid(p_img,Huh7_5_34,Huh_7_5_313, SHSY5Y_4.10_MOI1, SHSY5Y_4.17_MOI1,
                        labels = c("A","B","","C",""),
                        ncol = 1,
                        rel_heights = c(1,1,1,1,1))


middle_fig <- plot_grid(LRRC28_Tissue_gtx_normalized_plot, LRRC28_FeaturePlot,
                       labels= c("D","E"),
                       ncol = 1,
                       nrow = 2,
                       rel_heights = c(1.5,0.5))

final_figure<- plot_grid (right_plot, middle_fig, merged_dotplot,
                        labels = c("","","F"),
                        ncol = 3)


ggsave("Figure_1.pdf", final_figure, device = cairo_pdf, width = 20, height = 20)


SHSY5Y_4.10_MOI1_NTD=ggplot(SHSY5Y_Compiled[c(1:18,46:54),], aes(x = `Time (min)`, y = Value, color= Gene, group = Gene)) +
  geom_line(linewidth=1)+
  geom_point(size=2)+
  coord_cartesian(ylim = c(0, 700)) +   
  scale_y_continuous(breaks = seq(0, 700, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 1")+
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
  coord_cartesian(ylim = c(0, 750)) +  
  scale_y_continuous(breaks = seq(0, 750, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 1")+
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
  coord_cartesian(ylim = c(0, 750)) +   
  scale_y_continuous(breaks = seq(0, 750, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                    # Error bars
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 1")+
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
  coord_cartesian(ylim = c(0, 700)) +   
  scale_y_continuous(breaks = seq(0, 700, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   # Error bars
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 2")+
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
  coord_cartesian(ylim = c(0, 750)) +   
  scale_y_continuous(breaks = seq(0, 750, by = 200)) +   # <-- forces labels every 200
  geom_errorbar(aes(ymin = Value - StandardDeviation, ymax = Value + StandardDeviation), 
                width = 0.2) +                   # Error bars
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 2")+
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
  coord_cartesian(ylim = c(0, 750)) +   
  scale_y_continuous(breaks = seq(0, 750, by = 200)) +   # <-- forces labels every 200
  labs(x = "Time", y = "% ECAR", title = "SH-SY5Y Experiment 2")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 14,hjust = 0.5),panel.background = element_rect(fill = "white"),panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        panel.spacing = unit(0.5, "lines"))+theme(legend.position =c(0.85,0.7),legend.text = element_text(size =14),
                                                  legend.title = element_text(size = 14), 
                                                  legend.background = element_rect(fill = "transparent", color = NA),
                                                  legend.key = element_rect(fill = "transparent", color = NA)) +
  scale_color_manual(values = c("Cont-SHSY5Y" = "#56B4E9", "LRRC28 MOI10_4.3" = "#D55E00", "LRRC NTD MOI10"="#CC79A7"))

final_plot_sup <- plot_grid(SHSY5Y_4.10_MOI10, SHSY5Y_4.17_MOI10, 
                            SHSY5Y_4.10_MOI1_NTD, SHSY5Y_4.17_MOI1_NTD,
                            SHSY5Y_4.10_MOI10_NTD, SHSY5Y_4.17_MOI10_NTD,
                            labels = c("A","","B","","C",""),
                            nrow = 3)
                          
ggsave("Supplementary_Fig2.pdf", final_plot_sup,  device = cairo_pdf, width = 15, height = 10)

