library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(SeuratDisk)
# Provide the paths to the 2 files below:
level2 = readRDS("/home/buysdb/bbio1_projects/grosveld_2024_2/data_phase2/level2_dCarlin_version_03.rds")
scars = read.csv2("/home/buysdb/bbio1_repos/sccarlinpy/data/full_dataset_calls.csv", sep = ",", row.names = 1, dec = ".")


for(day in c( "day7m")){ # "day0", "day4","day5m",
  fpath =  paste0("seurat_data_", day,".h5seurat");
  SaveH5Seurat(level2[[day]], filename = fpath, overwrite = TRUE)
  Convert(fpath, dest = "h5ad", overwrite = TRUE)
  
  write_csv2(   level2[[day]]@meta.data, paste0("seurat_meta_", day,".csv"))
}
output_base_path = './scar_plots'

FILTER_MODE = "STRICT" # RELAXED, STRICT, RELAXED_ALL_ENTROPIES

if(FILTER_MODE=="STRICT"){
  ## Filter parameters:
  min_cells_per_scar = 8 # This selects the minimum number of observations per scar
  doublet_threshold = 0.5 # Select a doublet threshold, 0.5 = default (matches ~92% with transcriptome doublets), a higher value will set the threshold less strict, eg 0.6 will select less doublets, a value >1 will disable the filter altogether
  entropy_threshold = 1.6
  allow_ambiguous_alignments = FALSE # Set to FALSE to only use calls which can be uniquely mapped to one long read alignment
} else if(FILTER_MODE=="RELAXED") {
  ## Filter parameters:
  min_cells_per_scar = 8 # This selects the minimum number of observations per scar
  doublet_threshold = 0.6 # Select a doublet threshold, 0.5 = default (matches ~92% with transcriptome doublets), a higher value will set the threshold less strict, eg 0.6 will select less doublets, a value >1 will disable the filter altogether
  entropy_threshold = 2
  allow_ambiguous_alignments = TRUE # Set to FALSE to only use calls which can be uniquely mapped to one long read alignment
} else if(FILTER_MODE=="RELAXED_ALL_ENTROPIES") {
  ## Filter parameters:
  min_cells_per_scar = 8 # This selects the minimum number of observations per scar
  doublet_threshold = 0.6 # Select a doublet threshold, 0.5 = default (matches ~92% with transcriptome doublets), a higher value will set the threshold less strict, eg 0.6 will select less doublets, a value >1 will disable the filter altogether
  entropy_threshold = 100
  allow_ambiguous_alignments = TRUE # Set to FALSE to only use calls which can be uniquely mapped to one long read alignment
} else {
  stop("Select FILTER_MODE = RELAXED_ALL_ENTROPIES, RELAXED, STRICT")
}


#day = "day0"
for(day in c("day0", "day4","day5m", "day7m")){

selected_reduction = 'umap_buys'
daydata = level2[[day]] # Select day of interest here
output_folder = file.path(output_base_path, FILTER_MODE,  day)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)



if(!allow_ambiguous_alignments){
  scars = scars[scars$n_unique_alignments>0,]
}


scars = scars[scars[, "scar_posterior_doublet"] <= doublet_threshold, ]

# Join this to the dataset
daydata@meta.data <- cbind(daydata@meta.data, scars[match(rownames(daydata@meta.data), rownames(scars)), ])
daydata@meta.data$scar <- factor(daydata@meta.data$scar)

# Pick some good scars, which are frequent and have low entropy across the days
# Count the frequency of each scar
freq <- table(daydata@meta.data$scar)

freq = freq[freq>min_cells_per_scar]

# Prepare dataframe which orders the scars based on entropy
df_entropy = daydata@meta.data[,c('scar',"scar_entropy")][!duplicated(daydata@meta.data$scar),]
df_entropy = df_entropy[!is.na(df_entropy$scar),]
rownames(df_entropy)=df_entropy$scar
df_entropy = df_entropy[ rownames(df_entropy) %in% rownames(freq),]
df_entropy = df_entropy[order(df_entropy$scar_entropy),]
df_entropy = df_entropy[(df_entropy$scar_entropy<entropy_threshold),]
# Plot the clusters:
DimPlot(daydata, group.by = "clusters_buys")

if(!selected_reduction %in% names(daydata@reductions)){
  local_selected_reduction = names(daydata@reductions)[-1]
} else {
  local_selected_reduction = selected_reduction
}

# Plots a subset of scars together in one plot:
##############
selected_scars = rownames( df_entropy[1:10,] ) # Select some scars to plot 
selected_scar_meta  = data.frame(daydata@meta.data) # Make a copy
selected_scar_meta$scar <- as.factor(ifelse(selected_scar_meta$scar %in% selected_scars, as.character(selected_scar_meta$scar), NA)) # Set all unselected scars as NA
selected_scar_meta = merge(selected_scar_meta, daydata@reductions[[local_selected_reduction]]@cell.embeddings) # Add umap information
gplot = ggplot(selected_scar_meta, aes(x=!!sym(colnames(daydata@reductions[[local_selected_reduction]]@cell.embeddings)[1]), 
                               y=!!sym(colnames(daydata@reductions[[local_selected_reduction]]@cell.embeddings)[2]), 
                               color=scar, shape=scar)) + scale_shape_manual(values = rep(c(15,16,17,18,19,20),20)) +
  geom_point(data = selected_scar_meta[is.na(selected_scar_meta$scar), ], color = "#AAAAAA99", shape="o") + 
  geom_point(data = selected_scar_meta[!is.na(selected_scar_meta$scar), ],  size=3, alpha=0.9) + theme(panel.background = element_rect(fill = "white")) 


ggsave(
  file.path(output_folder, '10_scars.png'),
  gplot,
  width = 6.25,
  height = 4.25,
  dpi = 300
)

###############




# Plot scars one by one:
n_colors = 9
dir.create(file.path(output_folder, 'per_scar'), showWarnings = FALSE, recursive = TRUE)

scar_colors <- brewer.pal(n_colors, "Set1")
for(i in 1: dim(df_entropy)[1]){ # Plot all scars above the thresholds
  
  selected_scar <- rownames(df_entropy)[i]

  
  print(paste0(c(selected_scar, df_entropy[selected_scar,]), sep=" "))
  scar_color <- scar_colors[ ( ((i-1)%%n_colors) )+1]
  daydata@meta.data$has_selected_scar = ( daydata@meta.data[, "scar"] == selected_scar )
  umap_tx = daydata@reductions[[local_selected_reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    cbind(has_selected_scar = daydata@meta.data$has_selected_scar)
  umap_tx$has_selected_scar = factor(umap_tx$has_selected_scar, levels = c("FALSE", "TRUE"))
  umap_tx = umap_tx[order(umap_tx$has_selected_scar, decreasing = FALSE), ]
  
  gplot = ggplot(umap_tx, aes(x=!!sym(colnames(daydata@reductions[[local_selected_reduction]]@cell.embeddings)[1]), 
                            y=!!sym(colnames(daydata@reductions[[local_selected_reduction]]@cell.embeddings)[2]), color=has_selected_scar)) + 
    geom_point(data = umap_tx[umap_tx$has_selected_scar == "FALSE", ], color = "#AAAAAA99") + 
    geom_point(data = umap_tx[umap_tx$has_selected_scar == "TRUE", ],  pch=21, color='black', size=3, fill = scar_color) + theme(panel.background = element_rect(fill = "white")) +
      ggtitle(paste(selected_scar, '\nN:', freq[selected_scar], 'cells'))
  
  
  
  ggsave(
    file.path(output_folder, 'per_scar', paste0(selected_scar, '.png')),
    gplot,
    width = 4.25,
    height = 4.25,
    dpi = 300
  )
  
}
}

