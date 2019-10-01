# anti_HBV_drug_repositioning_from_lincs
An NMF based method to conduct anti-HBV drug repositioning from the large scale cellular response data LINCS.
This is the project to implement anti-HBV drug repositioning from LINCS data and GEO HBV infection data.
Last page update: **23/09/2019**

# Prerequiries
1. Install [R](https://www.r-project.org) on your computer
2. Download the project
    
    git clone https://github.com/hycyc/drug_repositioning_from_lincs
    
3. Download the [LINCS data](https://cbcl.ics.uci.edu/public_data/D-GEX/l1000_n1328098x22268.gctx), you need no less than 110GB free space on your computer.
4. Download the [annotation data for lincs](https://drive.google.com/file/d/19AlHVi2vv5T5hgvQ_uR6buygFQVHwKuD/view?usp=sharing)
5. Download the [GPL10558 platfrom annotation data](https://drive.google.com/file/d/12sjV2MJlPaTPNvc-ZFl540W1hoBqDMRv/view?usp=sharing)
6. Install the required R packages in requirements.rtf
7. Change all the 'read-in' file path in the file nmf_drug.R to the directory on your computer and run nmf_drug.R to get the cluster results

# the processed data
The [original gene expression signatures](https://pan.baidu.com/s/1YBCNr_XJlgXd0fWEq7VnlQ), [normalized signatures](https://pan.baidu.com/s/1LVBLVPEyn8CYdLAoC-uC-g) and [non-negative data matrix](https://pan.baidu.com/s/1266VDblgPoVB7kA_YOk5Wg) are available to download.
