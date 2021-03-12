# Transcriptomics


- #### [Blast2GO tutorial](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=73)


> Use this command if you have InterProScan result (make sure you edit the annotation.prop file, change the value next to InterProScanImportParameters.inputFormat):

```bash

/usr/local/blast2go/blast2go_cli.run -properties annotation.prop -useobo go.obo -loadblast blastresults.xml -loadips50 ipsout.xml -mapping -annotation  -annex -statistics all -saveb2g myresult -saveannot myresult -savereport myresult -tempfolder ./ >& annotatelogfile &

```



***Very important file format problem***


The output file from BLAST2GO myresult.annot is formatted for Windows PC. It will cause problems if you use it on Linux computer. Make sure to convert the file to LINUX line ends for this file. The command is:

```bash

dos2unix myresult.annot


```

Reference paper
[Comparative transcriptome analysis of cultivated and wild seeds of Salvia hispanica (chia)](https://www.nature.com/articles/s41598-019-45895-5#Sec9)

[Seed Transcriptomics Analysis in Camellia oleifera Uncovers Genes Associated with Oil Content and Fatty Acid Composition](https://www.mdpi.com/1422-0067/19/1/118)

[Formation of two methylenedioxy bridges by a Sesamum CYP81Q protein yielding a furofuran lignan, ()-sesamin](https://www.pnas.org/content/pnas/103/26/10116.full.pdf)

[De novo Assembly of Leaf Transcriptome in the Medicinal Plant Andrographis paniculata](https://www.frontiersin.org/articles/10.3389/fpls.2016.01203/full)

[Transcriptome Dynamics during Black and White Sesame (Sesamum indicum L.) Seed Development and Identification of Candidate Genes Associated with Black Pigmentation](https://www.mdpi.com/2073-4425/11/12/1399)

[De novo transcriptome assembly, annotation and comparison of four ecological and evolutionary model salmonid fish species](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4379-x#Sec2)

[Data of de novo assembly and functional annotation of the leaf transcriptome of Impatiens balsamina](https://www.sciencedirect.com/science/article/pii/S2352340918315816#s0020)

[Transcriptome Analysis by RNA–Seq Reveals Genes Related to Plant Height in Two Sets of Parent-hybrid Combinations in Easter lily (Lilium longiflorum)](https://www.nature.com/articles/s41598-020-65909-x#Sec12)

[Transcriptome analysis of pecan seeds at different developing stages and identification of key genes involved in lipid metabolism](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195913)

[Transcriptome Analysis Reveals Key Seed-Development Genes in Common Buckwheat (Fagopyrum esculentum)](https://www.mdpi.com/1422-0067/20/17/4303)

[Transcriptome analysis of metabolic pathways associated with oil accumulation in developing seed kernels of Styrax tonkinensis, a woody biodiesel species](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-020-2327-4#Sec21)

[Genomic and Transcriptomic Analysis Identified Gene Clusters and Candidate Genes for Oil Content in Peanut (Arachis hypogaea L.)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6061501/pdf/11105_2018_Article_1088.pdf)

[Comparative transcriptome analysis of the genes involved in lipid biosynthesis pathway and regulation of oil body formation in Torreya grandis kernels](https://www.sciencedirect.com/science/article/pii/S0926669019310611#sec0010)

[Transcriptome Analysis of Acer truncatum Seeds Reveals Candidate Genes Related to Oil Biosynthesis and Fatty Acid Metabolism](https://journals.ashs.org/jashs/view/journals/jashs/aop/article-10.21273-JASHS05013-20/article-10.21273-JASHS05013-20.xml) 


[RNA-seq Transcriptome Analysis of Panax japonicus, and Its Comparison with Other Panax Species to Identify Potential Genes Involved in the Saponins Biosynthesis](https://www.frontiersin.org/articles/10.3389/fpls.2016.00481/full#h3)

[De Novo Assembly and Annotation of the Chinese Chive (Allium tuberosum Rottler ex Spr.) Transcriptome Using the Illumina Platform](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4512717/)




