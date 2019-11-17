# anti_HBV_drug_repositioning_from_lincs
An NMF based method to conduct anti-HBV drug repositioning from the large scale cellular response data LINCS.
This is the project to implement anti-HBV drug repositioning from LINCS data and GEO HBV infection data.
Last page update: **17/11/2019**

# Prerequiries
1. Install [R](https://www.r-project.org) on your computer
2. Install the required R package [NMF](https://cran.r-project.org/web/packages/NMF/index.html), [R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html), [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) and [limma](http://www.bioconductor.org/packages/3.2/bioc/html/limma.html).
3. Download the project
    
    git clone https://github.com/hycyc/drug_repositioning_from_lincs
    
4. Download the [LINCS data](https://cbcl.ics.uci.edu/public_data/D-GEX/l1000_n1328098x22268.gctx), you need no less than 110GB free space on your computer.
5. Download the [annotation data for lincs](https://drive.google.com/file/d/19AlHVi2vv5T5hgvQ_uR6buygFQVHwKuD/view?usp=sharing) named inst.info.
6. Download the [GPL10558 platfrom annotation data](https://drive.google.com/file/d/12sjV2MJlPaTPNvc-ZFl540W1hoBqDMRv/view?usp=sharing)
7. Change all the '/your_path_to/' file path in the file nmf_drug.R to the directory on your computer and run nmf_drug.R to get the cluster results

# The pre-processed data
The zip file including annotation data for lincs, GPL10558 platfrom annotation data, original gene expression signatures, normalized signatures and non-negative data matrix are available to download at [Google Drive](https://drive.google.com/file/d/1OYRwsaV0wZepKc_nQVcHY7-tZmyMEH4T/view?usp=sharing), [Baidu SkyDrive](https://pan.baidu.com/s/11be_RWz8hcw7YAn34j2LGw) and [IEEE DataPort](https://ieee-dataport.org/documents/drug-repositioning-sitagliptin-promising-anti-hepatitis-b-virus-drug-candidate-silico).
