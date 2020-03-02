# Metagenome assembly and binning

In this tutorial you'll learn how to inspect assemble metagenomic data and retrieve draft genomes from assembled metagenomes

We'll use a mock community of 20 bacteria sequenced using the Illumina HiSeq.
In reality the data were simulated using [InSilicoSeq](http://insilicoseq.readthedocs.io).

The 20 bacteria in the dataset were selected from the [Tara Ocean study](http://ocean-microbiome.embl.de/companion.html) that recovered 957 distinct Metagenome-assembled-genomes (or MAGs) that were previsouly unknown! (full list on [figshare](https://figshare.com/articles/TARA-NON-REDUNDANT-MAGs/4902923/1) )

## Getting the Data

```bash
mkdir -p ~/data
cd ~/data
curl -O -J -L https://osf.io/th9z6/download
curl -O -J -L https://osf.io/k6vme/download
chmod -w tara_reads_R*
```

## Quality Control

we'll use FastQC to check the quality of our data, as well as sickle for trimming the bad quality part of the reads.
If you need a refresher on how and why to check the quality of sequence data, please check the [Quality Control and Trimming](qc) tutorial

```bash
mkdir -p ~/results
cd ~/results
ln -s ~/data/tara_reads_* .
fastqc tara_reads_*.fastq.gz
```

!!! question
    What is the average read length? The average quality?

!!! question
    Compared to single genome sequencing, what graphs differ?


Now we'll trim the reads using sickle

```
sickle pe -f tara_reads_R1.fastq.gz -r tara_reads_R2.fastq.gz -t sanger \
    -o tara_trimmed_R1.fastq -p tara_trimmed_R2.fastq -s /dev/null
```

!!! question
    How many reads were trimmed?

## Assembly

Megahit will be used for the assembly.

```
megahit -1 tara_trimmed_R1.fastq -2 tara_trimmed_R2.fastq -o tara_assembly
```

the resulting assenmbly can be found under `tara_assembly/final.contigs.fa`.

!!! question
    How many contigs does this assembly contain? How small is the smallest contig?

## Binning

First we need to map the reads back against the assembly to get coverage information

```bash
ln -s tara_assembly/final.contigs.fa .
bowtie2-build final.contigs.fa final.contigs
bowtie2 -x final.contigs -1 tara_reads_R1.fastq.gz -2 tara_reads_R2.fastq.gz | \
    samtools view -bS -o tara_to_sort.bam
samtools sort tara_to_sort.bam -o tara.bam
samtools index tara.bam
```

then we run metabat

```bash
runMetaBat.sh -m 1500 final.contigs.fa tara.bam
mv final.contigs.fa.metabat-bins1500 metabat
```

!!! question
    How many bins did we obtain?

!!! extra credit:
    What's in the `final.contigs.fa.depth.txt` file?
    Can you use it with some other tools you know of to make a "coverage vs. GC contents" plot maybe?, and maybe then plot on it the result of th contig clustering?
    


## Checking the quality of the bins

The first time you run `checkm` you have to create the database

```bash
checkm data setRoot ~/.local/data/checkm
```

```bash
checkm lineage_wf -x fa metabat checkm/
checkm qa  checkm/lineage.ms checkm > checkm.txt
checkm qa_plot -x fa checkm metabat plots
```

and take a look at `plots/bin_qa_plot.png`

## Functional annotation

Now we should have a collection of MAGs that we can further analyze. The first step is to predict genes as right now we only have raw genomic sequences. We will use one of my all-time-favorites : `prokka`.

This tool does gene prediction as well as some decent and usefull annotations, and is actually quiet easy to run!

> Use `prokka` loaded with `module`
> Predict genes and annotate your MAGs!

`prokka` produces a number of output files that all kind of represent similar things. Mostly variants of FASTA-files, one with the genome again, one with the predicted proteins, one with the genes of the predicted proteins. Also it renames all the sequence with nicer IDs! Additionally a very useful file generated is a GFF-file, which gives more information about the annotations then just the names you can see in the FASTA-files.

The annotations of `prokka` are good but not very complete for environmental bacteria. Let's run an other tool I like a lot, eggNOGmapper. This is a bit heavier in computation and it is not on  `UPPMAX`, so you will have to [install](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2) it.

> Install eggNOG-mapper
> run it on at least one MAG (don't be greedy, it is not fast!)
> OPTIONAL : undestand the output of it ...


## Taxonomic annotation

We now know more about the genes your MAG contains, however we do not really know who we have?! `checkm` might have given us an indication but it is only approximative.

Taxonomic classification for full genomes is not always easy for MAGs, often the 16S gene is missing as it assembles badly, and which other genes to use to for taxonomy is not always evident. Typically marker genes, min-hashes or k-mer databases are used as reference. It is often problematic for environmental data as the databases are not biased into our direction! We will use a min-hash database I compiled specifically for this (based on other tools and the full-datasets) using a tool called [sourmash](https://sourmash.readthedocs.io/en/latest/).

> Install sourmash!
> run `sourmash compute -k 31 --scaled 10000  metabat/bin.*` to compute signatures
> Use the lca clasify function of sourmash with the database found here /proj/g2019027/2019_MG_course/dbs/gtdbtk_sourmash.lca.json  

This database was made by the representative genomes of the gtdb [database](https://gtdb.ecogenomic.org/). This database is used by the tool `gtdbtk`, it uses marker genes and loads of data. It is a bit heavy, and tricky to run/install but much more sensitive.

> Optional: install and run gtdbtk, you can use the database here `/proj/g2019027/2019_MG_course/dbs/gtdbtk`



## Further reading

* [Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life](https://www.nature.com/articles/s41564-017-0012-7)
* [The reconstruction of 2,631 draft metagenome-assembled genomes from the global oceans](https://www.nature.com/articles/sdata2017203)
