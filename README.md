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
