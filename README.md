# WheatCo-expression
build co-expression network for wheat

A paper named ["Integrating Co-Expression Networks with GWAS to Prioritize Causal Genes in Maize"](http://www.plantcell.org/content/early/2018/11/09/tpc.18.00299) was published on the plant cell journal. This paper combined GWAS data and co-expression data for key genes mining in maize. I want to use this method in wheat.

#### 1 install camoco
``` shell
conda create -n camoco_env python=3.6
source activate camoco_env
pip install numpy
pip install camoco
# camoco version: 0.6.1
source deactivate # exit environment
```
#### 2 Start
``` shell
# Start environment
source activate camoco_env 
# Buid a RefGen Object
# The gff3 file which only contained 107891 high confidence genes was download from Ensembl Plants release 41

camoco build-refgen --ID-attr gene_id Triticum_aestivum.IWGSC.41.gff3 "TaIWGSCv1.1" "IWGSCv1.1 from Chinese Spring" 1.1 "Triticum aestivum" 
```
#### Build a COB object (co-expression network)
```shell
# Here, co-expression network was constructed only used HC (high confidence) genes in wheat.

camoco build-cob --min-expr 0.5 --min-single-sample-expr 5 BCS_HC_genes_tpm.tsv TaBCSdevelopment "Development of Triticum aestivum L. cv. BCS" TaIWGSCv1.1
```