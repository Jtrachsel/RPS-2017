# resistant-starch

### This is the code associated with the submission to Frontiers in Immunology titled: 
## Dietary Resistant Potato Starch Alters Intestinal Microbial Communities and Their Metabolites, and Markers of Immune Regulation and Barrier Function in Swine  
All of the analysis code is in the file FS2.R.  Hopefully it is somewhat understandable.  Suggestions are welcome.  

Several of the functions used in this analysis are housed in my package "funfuns".  
This package can be installed by using the following commands:  
(if you don't have devtools run `install.packages('devtools')`)

`devtools::install_github('Jtrachsel/funfuns')`  
`library(funfuns)`  
  
All of the funcitons should have documentation with them that can be accessed using `?`. For example: `?NMDS_ellipse()`

Additionally you will need to install my package `ggscinames`  
Documentation is available at https://github.com/Jtrachsel/ggscinames  
Install using: `install.packages('grid')` and `devtools::install_github('Jtrachsel/ggscinames')`  


### File descriptions 
many of these files are written out by the FS2.R script itself just for my convenience.  I'm sure there's a better way of organizing this and hopefully I'll be able to get around to that one day. 

| file &nbsp; &nbsp; &nbsp;| description &nbsp; &nbsp; &nbsp; &nbsp;|
|---|	---|
| 16S_meta_forcorr.txt | 16S metadata, used to help build networks |	
|16S_shared_forcorr.txt|16s OTU table used to build networks|	
|but2.shared|but gene amplicon OTU table|	
|but_meta_forcorr.txt|but gene amplicon metadata to help build networks|	
|butmeta.txt|but gene amplicon metadata|	
|butreps_blastn.txt|but gene OTU representative sequence blastn results|	
|butreps_blastx2.txt|but gene OTU representative sequence blastx results|	
|but_shared_forcorr.txt|but gene amplicon OTU table for network building|	
|Cassidy_rel_expr.txt|qPCR ct values from Cassidy|	
|CD3-.txt.csv|CD3 -  cell types flow data|	
|CD3+.txt.csv|CD3+ cell types flow data|	
|cec_redo_R.txt|SCFA data from GC|	
|Cec_tissue_forcorr.txt|processed qPCR ct values for cecal tissue network|	
|cecum.shared|OTU table for cecal tissue network|	
|dif_ab_16s.txt|Differentially abundant 16S OTUs as determined by DEseq2|	
|dif_ab_but.txt|Differentially abundant but gene OTUs as determined by DEseq2|	
|fecalplasma_R.txt|fecal SCFA data from GC|	
|FS2_all_flow.txt|all flow data, kinda like an OTU table|	
|miscforcorr.txt|IgA and other data for network building|	
|V4.final.shared|16s OTU table output by mothur|	
|V4.final.taxonomy|16S SILVA taxonomy output by mothur|	
|V4.FS2.4200.shared|OTU table rarrefied to 4200 sequences per sample|	
|V4.metadata.txt|metadata for 16S samples|	
|vfas_orig.txt|SCFA data from GC|	
