# Hotspot2Gene
Hotspot2Gene is a python package based on long reads WGS and RNA-seq. 

Hotspot2Gene consists of five integrated modules for detecting and characterizing structural variant (SV) hotspots in long-read cancer genomes. First, Structural Variant Calling generates raw SV callsets using long-read SV callers such as Sniffles2 and cuteSV. Next, Single-Sample Consensus SV Generation merges overlapping SVs within each sample to produce a consensus SV set. Cross-Sample Consensus SV Clustering then combines per-sample SVs to identify recurrent events shared across the cohort. SV Hotspot Detection applies Gaussian kernel density estimation (KDE) and shuffle-based statistical testing to pinpoint genomic regions with significant SV enrichment. Finally, SV–Expression Association links hotspots to nearby genes to identify candidate driver events.

<img width="612" height="300" alt="github" src="https://github.com/user-attachments/assets/2dade1d2-cedf-4da2-bc24-0c206dc43b87" />

