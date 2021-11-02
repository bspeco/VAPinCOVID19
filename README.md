# VAPinCOVID19
This the respository that contains the code used for bulk RNAseq, scRNAseq, and microbiome analysis.

# Bulk RNAseq
bulkDESeq.R: This R notebook can be used to reproduce differential gene expression and heatmap construction.

# scRNAseq 
PreprocessTrachAspirate_GITHUBnotations.ipynb: Data was processed and combined into adata ojbect  
scVI_GITHUBnotations.ipynb: scVI computation  
Clustering_lowclusterforDE_GITHUBannotations.ipynb: clustering (low resolution for entire mono/mac group - saved to take to R)  
Clustering_highres_GITHUBannotations.ipynb: higher resolution clustering, cell numbers, UMAPs  
R markdowns for pathway, volcano, and heatmaps: Rmarkdown_4GitHub.Rmd  
vap_classifier.ipynb: construction of the host-based classifier

# Microbiome Analyses
The following files contain code used for microbiome analyses.  
1_Format_Input_files.html is used to format the input files.
2_Coronavirus_percent.html is used to calculate the percent of reads assigned to SARS-CoV-2 in each sample.  
3_Calculate_Alpha_Beta.html is used to calculate the alpha and beta diversity.  
4_Heatmap_Top_Vap_Species_byVAPCAT.html is used to generate the heatmaps in Figure 6.   

# Host-based classifier
vap_classifier.ipynb: used to generate the host-based classifier to predict VAP.
