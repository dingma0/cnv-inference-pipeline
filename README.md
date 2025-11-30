# Copy Number Variation Inference Pipeline
Inference of copy number variation from scRNA-seq and scATAC-seq (TBA) data.

Author: Ding Ma

## Background
Copy-number variations (CNV) are a prominent feature of many cancer types, characterized by gains or losses of genomic regions [1]. These variations can reveal the deletion or amplification of cancer-related genes, subclonal structures, and genomic patterns associated with pathogenic phenotypes [2]. Traditionally, CNVs are profiled via bulk whole genome sequencing (WGS), which limits the resolution necessary to dissect intra-tumour heterogeneity [3]. To address this limitation, single-cell DNA sequencing (scDNA-seq) technologies was developed to enable characterization of CNVs on a single-cell basis. 

However, these methods remain constrained by high cost, low availability, and technical noises [4]. To overcome these challenges, a growing number of computational tools have been developed to infer CNVs from either scRNA-seq or single-cell ATAC sequencing (scATAC-seq) data [5]. 

CNV inference from RNA and ATAC data are not yet widely adopted in current single-cell analysis pipelines, as such, this project aims to enable reproducible workflows for popular CNV inference tools. This pipeline currently covers the basic workflow for Numbat, a tool for inferring allele-specific CNVs from scRNA-seq data: 

1. Matched WGS phasing (optional)
2. 1000Genomes Variant Reference Panel Preparation
3. SNP calling and phasing with cellsnp-lite and Eagle2
4. Numbat core algorithm


## References

[1] Ondrej Pös, Jan Radvanszky, Gergely Buglyó, Zuzana Pös, Diana Rusnakova, Bálint Nagy, and Tomas Szemes. DNA copy number variation: Main characteristics, evolutionary significance, and pathological aspects. 44(5):548–559. ISSN 23194170. doi: 10.1016/j.bj.2021.02.003. URL https://linkinghub.elsevier.com/retrieve/pii/S2319417021000093.

[2] Christopher D. Steele, Ammal Abbasi, S. M. Ashiqul Islam, Amy L. Bowes, Azhar Khandekar, Kerstin Haase, Shadi Hames-Fathi, Dolapo Ajayi, Annelien Verfaillie, Pawan Dhami, Alex McLatchie, Matt Lechner, Nicholas Light, Adam Shlien, David Malkin, Andrew Feber, Paula Proszek, Tom Lesluyes, Fredrik Mertens, Adrienne M. Flanagan, Maxime Tarabichi, Peter Van Loo, Ludmil B. Alexandrov, and Nischalan Pillay. Signatures of copy number alterations in human cancer. 606(7916):984–991. ISSN 0028-0836, 1476-4687. doi: 10.1038/s41586-022-04738-6. URL https://www.nature.com/articles/s41586-022-04738-6.

[3] David Lähnemann, Johannes Köster, Ewa Szczurek, Davis J. McCarthy, Stephanie C. Hicks, Mark D. Robinson, Catalina A. Vallejos, Kieran R. Campbell, Niko Beerenwinkel, Ahmed Mahfouz, Luca Pinello, Pavel Skums, Alexandros Stamatakis, Camille Stephan-Otto Attolini, Samuel Aparicio, Jasmijn Baaijens, Marleen Balvert, Buys De Barbanson, Antonio Cappuccio, Giacomo Corleone, Bas E. Dutilh, Maria Florescu, Victor Guryev, Rens Holmer, Katharina Jahn, Thamar Jessurun Lobo, Emma M. Keizer, Indu Khatri, Szymon M. Kielbasa, Jan O. Korbel, Alexey M. Kozlov, Tzu-Hao Kuo, Boudewijn P.F. Lelieveldt, Ion I. Mandoiu, John C. Marioni, Tobias Marschall, Felix Mölder, Amir Niknejad, Alicja R ˛aczkowska, Marcel Reinders, Jeroen De Ridder, Antoine-Emmanuel Saliba, Antonios Somarakis, Oliver Stegle, Fabian J. Theis, Huan Yang, Alex Zelikovsky, Alice C. McHardy, Benjamin J. Raphael, Sohrab P. Shah, and Alexander Schönhuth. Eleven grand challenges in single-cell data science. 21(1):31. ISSN 1474-760X. doi: 10.1186/s13059-020-1926-6. URL https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1926-6.

[4] Minfang Song, Shuai Ma, Gong Wang, Yukun Wang, Zhenzhen Yang, Bin Xie, Tongkun Guo, Xingxu Huang, and Liye Zhang. Benchmarking copy number aberrations inference tools using single-cell multi-omics datasets. 26(2):bbaf076. ISSN 1467-5463, 1477-4054. doi: 10.1093/ bib/bbaf076. URL https://academic.oup.com/bib/article/doi/10.1093/bib/bbaf076/8051529.

[5] Katharina T. Schmid, Aikaterini Symeonidi, Dmytro Hlushchenko, Maria L. Richter, Andréa E.Tijhuis, Floris Foijer, and Maria Colomé-Tatché. Benchmarking scRNA-seq copy number variation callers. 16(1):8777. ISSN 2041-1723. doi: 10.1038/s41467-025-62359-9. URL https://www.nature.com/articles/s41467-025-62359-9.

[6] Teng Gao, Ruslan Soldatov, Hirak Sarkar, Adam Kurkiewicz, Evan Biederstedt, Po-Ru Loh, and Peter V. Kharchenko. Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes. 41(3):417–426. ISSN 1087-0156, 1546-1696. doi: 10.1038/s41587-022-01468-y.URL https://www.nature.com/articles/s41587-022-01468-y.

[7] Flash-Frozen lymph node with B-cell lymphoma (14k sorted nuclei), Single Cell Multiome ATAC + Gene Expression Dataset, 10x Genomics (2021, May 3).

[8] PBMCs from a healthy donor (3 k, no cell sorting), Single Cell Multiome ATAC + Gene Expression Dataset, 10x Genomics (2021, May 3).