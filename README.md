# Transcriptomics

# Qulaity control
[Reference site](https://genomics.sschmeier.com/ngs-qc/index.html#the-qc-process)

Create the channel for qc to install fastp fastqc multiqc install


* Step:1 
* Tool installation

```bash
# create env and install tools
$ conda create --yes -n qc fastp fastqc multiqc

# activate env
$ conda activate qc
```
* Step: 2
```bash

conda create --yes -n qc fastp fastqc multiqc

----------------------------

vi run_qc.sh


#!/bin/bash

set -e


#--Activate the qc environment

source activate qc


#-- Check the reads quality

cd ~/analysis/015_seed_transcriptomics/01.quality_preprocessing/raw

fastqc -t 32 -o ~/analysis/015_seed_transcriptomics/01.quality_preprocessing/raw *.fastq.gz &> log.fastqc &
```

#---Trimming

```bash

#!/bin/bash

set -e

for reads in *_1.fastq.gz

do 

base=$(basename $reads _1.fastq.gz)

fastp --detect_adapter_for_pe \
       --overrepresentation_analysis \
       --correction --cut_right --thread 2 \
       --html ../trim/${base}.fastp.html --json ../trim/${base}.fastp.json \
       -i ${base}_1.fastq.gz -I ${base}_2.fastq.gz \
       -o ../trim/${base}_R1.fastq.gz -O ../trim/${base}_R2.fastq.gz 

done
#----Run
bash run_fastp.sh &>log.fastp &
```

#---Let's run fastqc again on the trimmed data
```bash
cd ../trim/

fastqc -t 32 *

#--Let's run multiQC on both untrimmed and trimmed files

cd ..

multiqc raw trim 

source deactivate qc

#--Bye


$ bash run_qc.sh &> log.qc &
        
```

Note:

--detect_adapter_for_pe: Specifies that we are dealing with paired-end data.

--overrepresentation_analysis: Analyse the sequence collection for sequences that appear too often.

--correction: Will try to correct bases based on an overlap analysis of read1 and read2.

--cut_right: Will use quality trimming and scan the read from start to end in a window. If the quality in the window is below what is required, the window plus all sequence towards the end is discarded and the read is kept if its still long enough.

--thread: Specify how many concurrent threads the process can use.

--html and --json: We specify the location of some stat files.

-i data/anc_R1.fastq.gz -I data/anc_R2.fastq.gz: Specifies the two input read files

-o trimmed/anc_R1.fastq.gz -O trimmed/anc_R2.fastq.gz: Specifies the two desired output read files





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
* Transcriptome Quality Assessment
* [BUSCO](https://angus.readthedocs.io/en/2019/dammit_annotation.html#evaluation-with-busco)
* 
```bash


#!/bin/bash

set -e

cd /NABIC/HOME/senthil83/analysis/009_Assembly_quality_metrics_s_alatum

ln -fs /NABIC/HOME/senthil83/analysis/004_denovo_assembly_s_alatum_rna/s_alatum_trinity_out/Trinity.fasta .


echo "Activate the busco environment"

source activate busco_env

echo "BUSCO analysis for S_alatum" 

busco \
  -i Trinity.fasta \
  -o busco_S_alatum_embryophyta \
  -l  /NABIC/HOME/senthil83/analysis/009_Assembly_quality_metrics_s_alatum/eudicots_odb10.2020-09-10/eudicots_odb10 \
  -m transcriptome 


#deactivate busco environment
source deactivate busco_env


bash run.sh &> log &
```

#exit

* References

[Transcriptome Assembly Quality Assessment](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment)


#### Trinotate
[Prototype Ref1](https://informatics.fas.harvard.edu/trinotate-workflow-example-on-odyssey.html)
[Githubref2](https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required)


[A Tissue-Mapped Axolotl De Novo TranscriptomeEnables Identification of Limb Regeneration Factors](https://www.sciencedirect.com/science/article/pii/S2211124716317703?via%3Dihub#mmc9)


### Diamond
[DIAMOND](https://github.com/bbuchfink/diamond)

##Install DIAMOND BLASTX tool
https://anaconda.org/bioconda/diamond

##create specific environment for diamond
conda create --name diamond_env

source activate diamond_env

conda install -c bioconda diamond
conda install -c bioconda/label/cf201901 diamond

https://github.com/bbuchfink/diamond


# downloading the tool
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
# creating a diamond-formatted database file
./diamond makedb --in reference.fasta -d reference
# running a search in blastp mode
./diamond blastp -d reference -q queries.fasta -o matches.tsv
# running a search in blastx mode
./diamond blastx -d reference -q reads.fasta -o matches.tsv
# downloading and using a BLAST database
update_blastdb --decompress --blastdb_version 5 swissprot
./diamond blastp -d swissprot -q queries.fasta -o matches.tsv


##Making database###
```bash
/usr/bin/time -o output.txt -v diamond makedb -p 32 --in uniprot_sprot.fasta -d uniprot_sprot &>log &

/usr/bin/time -o outputswissprot.txt -v diamond makedb -p 32 --in swissprot.fasta -d swissprot &>log.swissprt &

/usr/bin/time -o outputnr.txt -v diamond makedb -p 32 --in nr.fasta -d nr &>log.nrdb &
```
##Run Blastx##
```bash
/usr/bin/time -o outputblastx.txt -v diamond blastx -p 32 -d uniprot_sprot -q alatum.transcripts.fasta -o alatum_matches.xml --outfmt 5 &> log.blastx.xml &

/usr/bin/time -o output_sprotblastx.txt -v diamond blastx -p 32 -d swissprot -q alatum.transcripts.fasta -o swissprot.alatum_matches.xml --outfmt 5 &> log.swisprot.blastx.xml &

/usr/bin/time -o output_nrblastx.txt -v diamond blastx -p 32 -sensitive -d nr -q alatum.transcripts.fasta -o nr.alatum_matches.xml --outfmt 5 &> log.nr.blastx.xml &

````

##### 
# CDD DATABASE ##
[NUCELOTIDE QUERY REF](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#NucleotideQuerySequence)
[Biostar_help](https://www.biostars.org/p/52309/)
[Biopython](https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/rpsblast/)

#####
# GO ENRICHMENT using R ##

[GOseq_R](https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-gene-set-testing.nb.html)


# AUGUSTUS GENE PREDICTION TOOL
[GENOMICS TUTORIAL](https://genomics.sschmeier.com/ngs-annotation/index.html)

# Using Sequence_IDS pooled from trasncript_dataset

$ seqtk subseq your.input.fasta the_header_of_interest_IDs.list > your_output.fasta

# Sequence Length distribution

$ cat file.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'

[How to generate sequence length distribution from Fasta file](https://www.biostars.org/p/421925/)

# PFAM
For PFAM annotations, I used hmmscan against the full database downloaded from the xfam page. To download the databank and create the hmmscan file:

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz #download

$ gunzip Pfam-A.hmm.gz #unzip the database

$ hmmpress Pfam-A.hmm #parse the database

And to run the scan against the query sequences:

$ hmmscan --tblout MySu_v1.PFAM.txt Pfam-A.hmm MySu01_v1.proteins.fasta

# InterproScan
[Download](http://www.ebi.ac.uk/interpro/download/)

[To install and run InterProScan](https://github.com/ebi-pf-team/interproscan)

[To Run](https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html)

```bash
./interproscan.sh -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,PROSITEPATTERNS,PROSITEPROFILES,SFLD,SMART,SUPERFAMILY,TIGRFAM -i /path/to/sequences.fasta

```

# SSR marker prediction_KRAIT #

[KRAIT](https://github.com/lmdu/krait/releases)

[Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design] (https://academic.oup.com/bioinformatics/article/34/4/681/4557187)

Du, L., Zhang, C., Liu, Q., Zhang, X. & Yue, B. Krait: an ultrafast tool for genome-wide survey of microsatellites and primer design. Bioinformatics 34, 681–683 (2018).

[Distribution patterns and variation analysis of simple sequence repeats in different genomic regions of bovid genomes
Download PDF](https://www.nature.com/articles/s41598-018-32286-5#Sec12)

[Genome-wide investigation of microsatellite polymorphism in coding region of the giant panda (Ailuropoda melanoleuca) genome: a resource for study of phenotype diversity and abnormal traits](https://link.springer.com/article/10.1007/s13364-019-00418-5?shared-article-renderer#Sec2)


# Reference paper ##

[Transcriptomic comparison reveals genetic variation potentially underlying seed developmental evolution of soybeans](https://academic.oup.com/jxb/article/69/21/5089/5064903#122474828)


[Methodology references---Identification of Glutathione S-Transferase Genes in Hami Melon (Cucumis melo var. saccharinus) and Their Expression Analysis Under Cold Stress](https://www.frontiersin.org/articles/10.3389/fpls.2021.672017/full#h4)

[EnTAP: Bringing Faster and Smarter Functional Annotation to Non-Model Eukaryotic Transcriptomes](https://www.biorxiv.org/content/biorxiv/early/2018/04/28/307868.full.pdf)

[A Comparison of Resources for the Annotation of a De Novo Assembled Transcriptome in the Molting Gland (Y-Organ) of the Blackback Land Crab, Gecarcinus lateralis](https://academic.oup.com/icb/article/56/6/1103/2647085)

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

[Seed Transcriptome Annotation Reveals Enhanced Expression of Genes Related to ROS Homeostasis and Ethylene Metabolism at Alternating Temperatures in Wild Cardoon](https://www.mdpi.com/2223-7747/9/9/1225/htm)


[Gene expression profiles that shape high and low oil content sesames](https://bmcgenomdata.biomedcentral.com/track/pdf/10.1186/s12863-019-0747-7.pdf)

[Transcriptome profiling of spike provides expression features of genes related to terpene biosynthesis in lavender](https://www.nature.com/articles/s41598-020-63950-4#Sec11)

[Metabolic engineering of fatty acid biosynthetic pathway in sesame (Sesamum indicum L.): assembling tools to develop nutritionally desirable sesame seed oil](https://link.springer.com/article/10.1007/s11101-015-9424-2)

[## unreplicated data_Comparative Characterization of the Leaf Tissue of Physalis alkekengi and Physalis peruviana Using RNA-seq and Metabolite Profiling](https://www.frontiersin.org/articles/10.3389/fpls.2016.01883/full#h3)

[Full-Length Transcriptome Sequencing and Comparative Transcriptome Analysis to Evaluate Drought and Salt Stress in Iris lactea var. chinensis](https://www.mdpi.com/2073-4425/12/3/434/htm#B23-genes-12-00434)

[Lamiales_ref_Functional characterization of the cytochrome P450 monooxygenase CYP71AU87 indicates a role in marrubiin biosynthesis in the medicinal plant Marrubium vulgare](https://bmcplantbiol.biomedcentral.com/articles/10.1186/s12870-019-1702-5)

[Genome sequencing of the important oilseed crop Sesamum indicumL (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-401#ref-CR37)

[Phylogenomic analysis of cytochrome P450 multigene family and their differential expression analysis in Solanum lycopersicum L. suggested tissue specific promoters](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5483-x#Sec10)
![img](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs12864-019-5483-x/MediaObjects/12864_2019_5483_Fig2_HTML.png?as=webp)

[Computational Identification and Systematic Classification of Novel Cytochrome P450 Genes in Salvia miltiorrhiza](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0115149)
![img](https://journals.plos.org/plosone/article/figure/image?size=inline&id=info:doi/10.1371/journal.pone.0115149.g002)



## Sesamin/Sesamolin/Lignan ###

[(+)‐Sesamin‐oxidising CYP92B14 shapes specialised lignan metabolism in sesame](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14989)

[Glycoside‐specific glycosyltransferases catalyze regio‐selective sequential glucosylations for a sesame lignan, sesaminol triglucoside](https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.14586)

[Formation of a Methylenedioxy Bridge in (+)-Epipinoresinol by CYP81Q3 Corroborates with Diastereomeric Specialization in Sesame Lignans](https://academic.oup.com/pcp/article/59/11/2278/5062626?login=true)

[Oxidative rearrangement of (+)-sesamin by CYP92B14 co-generates twin dietary lignans in sesame](https://www.nature.com/articles/s41467-017-02053-7)

[Formation of two methylenedioxy bridges by a Sesamum CYP81Q protein yielding a furofuran lignan, (+)-sesamin](https://www.pnas.org/content/103/26/10116.short)

[Sequential glucosylation of a furofuran lignan, (+)‐sesaminol, by Sesamum indicum UGT71A9 and UGT94D1 glucosyltransferases](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-313X.2008.03428.x) 

[Variations in the composition of sterols, tocopherols and lignans in seed oils from fourSesamum species](https://sci-hub.se/https://doi.org/10.1007/BF02541549)

[Lignan analysis in seed oils from fourSesamum species: Comparison of different chromatographic methods](https://aocs.onlinelibrary.wiley.com/doi/epdf/10.1007/BF02541548)

[Lignans of Sesame (Sesamum indicum L.): A Comprehensive Review](https://www.mdpi.com/1420-3049/26/4/883)

[Candidate genes involved in the biosynthesis of lignan in Schisandra chinensis fruit based on transcriptome and metabolomes analysis](https://www.sciencedirect.com/science/article/pii/S1875536420600073)

[The cytochrome P450 superfamily: Key players in plant development and defense](https://www.sciencedirect.com/science/article/pii/S2095311914609801)


[Identification of a binding protein for sesamin and characterization of its roles in plant growth](https://www.nature.com/articles/s41598-019-45003-7)



