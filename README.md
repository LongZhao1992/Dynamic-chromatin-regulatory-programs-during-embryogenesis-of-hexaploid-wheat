# Dynamic-chromatin-regulatory-programs-during-embryogenesis-of-hexaploid-wheat
Here, we deposit the supporting code and data of the publication **"Dynamic chromatin regulatory programs during embryogenesis of hexaploid wheat"** (Zhao et al., 2022) for availability to the broad scientific community. The code and data deposited here are free for academic use. Any other types of use are prohibited.

## Abstract
### Background
Plant and animal embryogenesis have conserved and distinct features. Cell fate transitions occur during embryogenesis in both plants and animals. The epigenomic processes regulating plant embryogenesis remain largely elusive.
### Results
Here, we elucidate chromatin and transcriptomic dynamics during embryogenesis of the most cultivated crop, hexaploid wheat. Time-series analysis reveals stage-specific and proximal-distal distinct chromatin accessibility and dynamics concordant with transcriptome changes. Following fertilization, the remodeling kinetics of H3K4me3, H3K27ac and H3K27me3 differ from that in mammals, highlighting considerable species-specific epigenomic dynamics during zygotic genome activation. Polycomb repressive complex 2 (PRC2) mediated H3K27me3 deposition is important for embryo establishment. Later H3K27ac, H3K27me3, and chromatin accessibility undergo dramatic remodeling to establish a permissive chromatin environment facilitating the access of transcription factors to cis-elements for fate patterning. Embryonic maturation is characterized by increasing H3K27me3 and decreasing chromatin accessibility, which likely participates in restricting totipotency while preventing extensive
organogenesis. Finally, epigenomic signatures are correlated with biased expression among homeolog triads and divergent expression after polyploidization, revealing an epigenomic contributor to subgenome diversification in an allohexaploid genome.
### Conclusions
Collectively, we present an invaluable resource for comparative and mechanistic analysis of the epigenomic regulation of crop embryogenesis.

## Contents of this repository
* `metadata for ChIPseqSpikeInFree`	The directory contains the data used for ChIPseqSpikeInFree analysis.
* `0.data preprocess.sh`	The pipeline of CUT&Tag preprocess, including alignment and peak calling.
* `1.Fig1(PCA+tree).R`	The R script for PCA and clustering analysis in fig. 1. 
* `2.Fig3(chromHMM).R`	The R script for chromatin state analysis using ChromHMM result in fig.3.
* `3.Fig5(heatmap+network).R`	The R script for Heatmap and gene regulatory network analysis, including the module selection in fig.5A, PCA and pseudotime index calculation in fig. 5B, heatmap in fig. 5C and network construction.
* `4.ChIPseqSpikeInFree.R`	The R script for scale factor calculation in this project. R package ["ChIPseqSpikeInFree"](https://github.com/stjude/ChIPseqSpikeInFree) was used and the metadata for analysis was in the "metadata for ChIPseqSpikeInFree" directory.

Citation: Zhao, L., Yang, Y., Chen, J. et al. Dynamic chromatin regulatory programs during embryogenesis of hexaploid wheat. Genome Biol 24, 7 (2023). https://doi.org/10.1186/s13059-022-02844-2

