# trDUHi_FVB


### Description of the repository


Comparison of the genes expressed in lower and upper first molar germs in DUHi and FVB 

RNA-Seq were made on twelve carefully dissected embryonic lower and upper first molar germs, from DUHi (embryo weight: 196, 219 and 239 mg) and FVB (195, 215 and 233 mg) strains. 

Reads were then mapped to the mouse genome using Kallisto ((Bray et al., 2016), version 0.44.0, options -l 200, -s 20). The reference cDNA sequences and annotation files for M. musculus are based on C57B6 strain. They were collected from Ensembl 88 (10 5129 cDNAs, (Zerbino et al., 2018), GRCm38). Reads were independently mapped to the FVB/NJ strain cDNAs, collected from Ensembl strains 94, using biomart (10 1520 cDNAs, strain FVB_NJ_v1, accession GCA_001624535.1). Tximport was used to import and summarize transcript-level estimates at gene level (version 1.6, (Soneson et al., 2016)). Differentially expressed genes were detected with DESeq2 (Love et al., 2014), version 1.18.1) with classical one-factor design, and using FDR significance threshold = 0.05. 

### Usage 

The R code allows to reproduce differential analysis between strains and between reference sequence (allowing to retrieve mapping artifacts) : RNASeq_DUHiFVB_Mapping.R  

This will automatically load the mapping data (Mapping folder, for mapping on classical mouse genome : Quant_Mus, and for mapping on FVB genome strain : Quant_FVB).  The annotation of the files is in "sample.txt". CaracteristicsERCC.txt is necessary to  study ERCC spike ins. Annotations retrieved from Ensembl is in Biomart_GRCm38.p5.txt, with corresponding ID in FVB corresp_FVB.txt.
