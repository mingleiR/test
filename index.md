
Multi-omics technologies have brought about a great advance in discovering microbial dark matter. However, how to utilize  metagenomic datasets to gain insight into the adaptive strategies based on along environmental gradients is still a limited and ambiguous step in the field of microbial ecology. Here, we developed a functional trait-based analytic framework, which integrates the abundance of metagenome-assembled genomes (MAGs) and community-weighted mean of functional genes to pinpoint the functional traits which covary with environmental gradients and improve our ability to understand the mechanisms driving microbial adaptation. 

The step-by-step procedure, including the codes and scripts for sequence and statistical analyses implemented in the trait-based framework is available below. The marker gene-based taxonomic profiling analysis using metagenomic assembly, genome-centric analysis, and functional annotation analysis. Please see our manuscript for the references of the bioinformatic tools called below. If the user has reconstructed a list of non-redundant MAGs from multiples samples, then only the abundance of these MAGs across samples and the functional traits for these MAGs are required to analyze the dataset using the trait-based framework. Us

The content list of the procedure is shown below.
1. Sequence analysis
  -  Clean metagenomic reads
  -  Assemble metagenomic reads and predict CDS
  -  Reconstruct metagenome-assembled genomes (MAGs)
  -  Identify species groups (SG) using rpS3 gene
  -  Screen the representative SG with MAGs
  -  Estimate functional traits for the representative SG

2. Statistical analysis
  -  Calculate the community-weighted means (CWM) of functional traits
  -  Identify the key functional traits through trait-environment relationship (TER)
  -  Calculate the index of environmental adaptation based on TER
  -  Fit the relationships of the species' index and its niche optima
  

## 1. Sequence analysis

### 1-1. Clean metagenomic reads

#### 1-1-1. Detect the potential adapter and contamination

```
fastq_f=$1
cores=$2

## reading fastq file in the command
fastqc --format fastq  --extract --threads $cores ${fastq_f}
```

#### 1-1-2. Trim the potential contamination sequence

According to the summary result from FastQC, the appropriate adapter set could be set in the Trimmomatic program.
```
sample="$1"
cores="$2"

r1=${sample}_1.fastq.gz
r2=${sample}_2.fastq.gz

# minimum length for trimmed reads (50 is okay for assembly)
min_len_thresh=50
# quality for a sliding window of 4 bases
quality_thresh=25

r1_p=${sample}_1.fq.gz
r1_s=${sample}_1_unpaired.fq.gz
r2_p=${sample}_2.fq.gz
r2_s=${sample}_2_unpaired.fq.gz
log=trimmomatic.log.${sample}
single=${sample}_unpaired.fq.gz
adapter="/apps/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa"

trimmomatic PE -threads $cores $r1 $r2 ${r1_p} ${r1_s} ${r2_p} ${r2_s} HEADCROP:2 ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:4:${quality_thresh} MINLEN:${min_len_thresh} >$log 2>&1
cat ${r1_s} ${r2_s} > $single

```

Note that it is necessary to check the quality of trimmed reads using FastQC to make sure that the contaminated adapters have been trimmed from the sequence. 


### 1-2. Assemble metagenomic reads and predict CDS

After cleaning the reads, MEGAHIT, a memory-efficient and fast metagenomic assembler, is selected to perform de novo assembly for each individual sample. Other popular metagenomic assembler of SPAdes, or the co-assembly strategy are also encouraged when the computational resource is sufficient.

```
sample="$1"
cores="$2"

r1=${sample}_1.fq.gz
r2=${sample}_2.fq.gz
single=${sample}_unpaired.fq.gz

# using default presets (meta-sensitive)
min_len=1000

megahit -1 $r1 -2 $r2 -r $single -o megahit_${sample} -m 0.95 -t $cores --min-contig-len ${min_len}

```

To keep track the fate of each contigs and facilitate downstream analysis, the header of sequence in each contigs file are formatted from `>k141_81535 flag=1 multi=8.0000 len=1108` to `>sample-name_contig0000001`, and the contigs file generated from megahit are renamed from `final.contigs.fa` to `megahit_sample-name_contigs.fa`.

The protein-coding genes are predicted by Prodigal for each metagenomic assembly, and the generated amino acid sequence is prepared for the downstream analysis, including the marker gene-based taxonomic profiling and the DAS-Tool binning program (See below).

```

sample_contig_f="$1"
prodigal -p meta -i ${sample_contig_f} -a ${sample_contig_f}.genes.faa -d ${sample_contig_f}.genes.fna

```

### 1-3. Reconstruct metagenome-assembled genomes (MAGs)

The DAS-Tool program is used to recover metagenome-assembled genomes (MAGs) from each assembly. Before running the DAS-Tool, the four widely-used binning tools were firstly used to reconstruct preliminary MAGs from each assembly, including MaxBin2, MetaBat1, MetaBat2 and CONCOCT. Then, these binning results were integrated by the DAS Tool, resulting in a final list of dereplicated MAGs. The quality of each MAG, including completeness and contamination were evaluated using CheckM, and the taxonomic assignment is performed by GDTB-tk.

#### 1-3-1. Align the reads against the contigs

For each of samples, The clean reads are aligned against the contigs from the samples using Bowtie 2. Then the alignment is sorted using Samtools for the next step.

```

sample="$1"
cores="$2"

refer="megahit_${sample}_contigs.fa"

r1="${sample}_1.fq.gz"
r2="${sample}_2.fq.gz"
single="${sample}_unpaired.fq.gz"

out="${sample}_sorted.bam"

bowtie2-build $refer $refer
bowtie2 --threads $cores -x $refer -1 $r1 -2 $r2 -U $single | samtools view -uT $refer - | samtools sort --reference $refer --threads $cores -O BAM -o $out -

```

#### 1-3-2. Perform initial binning using 4 binning algorithms 

The CONCOCT program is run as below.

```
sample="$1"
min_len_threshold="$2"
cores="$3"

refer="megahit_${sample}_contigs.fa"

read_len=150
bed="${refer}_split.bed"
split_refer="${refer}_split.fa"
depth="${sample}_coverage_table.tsv"
bam="${sample}_sorted.bam"
# output directory 
bin_out="concoct_${sample}"

cut_up_fasta.py $refer -c 10000 -o 0 --merge_last -b $bed > ${split_refer}

samtools index -@ $cores ${bam}
concoct_coverage_table.py $bed ${bam} > $depth

concoct --threads $cores --length_threshold ${min_len_threshold} --read_length ${read_len} --composition_file ${split_refer} --coverage_file $depth -b ${bin_out}/
merge_cutup_clustering.py ${bin_out}/clustering_gt${min_len_threshold}.csv > ${bin_out}/clustering_merged.csv
extract_fasta_bins.py $refer ${bin_out}/clustering_merged.csv --output_path ${bin_out}/

```

The METABAT2 program is run in two separate steps as below.

```
sample="$1"

refer="megahit_${sample}_contigs.fa"

##########################
bam="${sample}_sorted.bam"
out="${sample}_depth.txt"

jgi_summarize_bam_contig_depths --outputDepth $out ${bam}

```

```
sample="$1"
min_len_threshold="$2"
cores="$3"
refer="megahit_${sample}_contigs.fa"
out="${sample}_depth.txt"
method='metabat2'

## output directory and prefix of each binning MAGs
bin_out="${method}_${sample}/${method}"
metabat2 --numThreads $cores -i $refer -a $out -o ${bin_out} --minContig ${min_len_threshold}
```

The METABAT1 program is run as below. Note that the coverage depth generated from METABAT2 is re-used in the step. Although METABAT2 is recommended to replace METABAT1 by the author of the software, our prior test show that the binning result from METABAT1 are often selected as the best results for the test samples.

```
sample="$1"
min_len_threshold="$2"
cores="$3"

refer="megahit_${sample}_contigs.fa"
out="${sample}_depth.txt"

## output directory and prefix of each binning MAGs
method='metabat1'
bin_out="${method}_${sample}/${method}"

metabat1 --numThreads $cores -i $refer -a $out -o ${bin_out} --minContig ${min_len_threshold}

```

The MAXBIN2 program is run as below. Note that the coverage depth generated from METABAT2 is re-used in the step.

```
sample="$1"
min_len_threshold="$2"
cores="$3"

# reuse the depth file generated from MetaBat2
depth="${sample}_depth.txt"
refer="megahit_${sample}_contigs.fa"

# output directory (the prefix of each binning MAGs use the same name as the directory)
method='maxbin'
bin_dir="${method}_${sample}"
mkdir ${bin_dir}
cd ${bin_dir}
ln -s ../$refer .
cd ../
bin_out="${bin_dir}/${method}"

# reuse the depth file generated from MetaBat2
tail -n+2 $depth |cut -f1,3 > maxbin.${sample}.depth.txt
run_MaxBin.pl -abund maxbin.${sample}.depth.txt -contig $refer -out ${bin_out} -min_contig_length ${min_len_threshold} -thread $cores
```

The result of these four binning algorithms are summarized and exported into a file with tow columns, namely, contig ID and binning ID. The summary file is named as `megahit_sample-name_binning-name`.

#### 1-3-3. Select the optimum MAGs results

The DAS-Tool is run as below. Note that the protein sequence for each metagenomic assembly has been predicted above. To recover more genomes from each sample, the cutoff for the score in the program is adjusted to 0.3.

```
sample="$1"
cores="$2"
cwd=$(pwd)
refer="megahit_${sample}_contigs.fa"
refer_protein="${refer}.genes.faa"
out_dir="${cwd}/binning_${sample}"
mkdir ${out_dir}
ln -s $cwd/${refer_protein} ${out_dir}/${refer_protein}
out="${out_dir}/run1"
methods=('concoct' 'maxbin' 'metabat1' 'metabat2')
prefix="megahit_${sample}"

dastool_score_thresh=0.3

DAS_Tool -i ${prefix}_${methods[0]}.tsv,${prefix}_${methods[1]}.tsv,${prefix}_${methods[2]}.tsv,${prefix}_${methods[3]}.tsv  -l concoct,maxbin,metabat1,metabat2 -c $refer  -o $out --threads $cores --create_plots 0  --search_engine blast --score_threshold ${dastool_score_thresh} --proteins ${out_dir}/${refer_protein}

```

After running DAS-Tool, the nucleotide sequence for the MAGs is extracted and exported in the suffixed ".fa". Then, the MAGs are filtered based on the completeness threshold of 40%, which is retrieved in the output of the program suffixed with "_DASTool_scaffolds2bin.txt".


#### 1-3-4. Evaluate the quality of MAGs

The quality of MAGs is evaluated using CheckM, and only the MAG with the completeness >=70% and the contamination <=10% are subject to the downstream analysis.
```
genome_dir="$1"
cores=$2

suffix="fa"
output_dir="${genome_dir}.checkm"
output_sum="${genome_dir}.checkm_seq_summary.tsv"
output_tree="${genome_dir}.checkm_marker_gene_summary.tsv"
checkm lineage_wf -x $suffix  --threads $cores --pplacer_threads $cores --quiet ${genome_dir} ${output_dir}
checkm qa --threads $cores --quiet --tab_table -f ${output_sum}  --out_format 2 ${output_dir}/lineage.ms ${output_dir}
checkm tree_qa  --out_format 2  --tab_table  -f ${output_tree}  ${output_dir}
```

#### 1-3-5. Assign taxonomy for MAGs

The taxonomic assignment is performed by GTDB-tk.
```
genome_dir=$1
suffix="fna"
cores=$2
out_dir=${genome_dir}.gtdb.classify
gtdbtk classify_wf --extension $suffix --genome_dir ${genome_dir} --out_dir ${out_dir} --cpus $cores

```

### 1-4. Identify species groups (SG) using rpS3 gene

Taxonomic profiling of microbial communities is performed based on a conserved ribosomal protein rpS3 through a modified pipeline described previously (Diamond et al 2019). Not all microbial genomes are recovered from metagenomic sequencing due to insufficient sequencing, therefore, the approach of taxonomic profiling based on marker genes from metgenomic assembly could provide more comprehensively overview of microbial community compared to the genome-centric approach.


#### 1-4-1. Identify and cluster the rpS3 genes from contigs of each sample 

For each of the samples, the amino acid sequence is screened to identify the rpS3 genes, then these genes are clustered at the 99% similarity with USEARCH. The script for screening and the relevant HMM files are copied and modified the one from [Prof. Probst's github](https://github.com/AJProbst/rpS3_trckr).

```
script_dir="$1"
sample_contig_file="$2"
cores="$3"

mkdir -p intermediate_files

## ++++++++++++++ FUNCTIONS ++++++++++++++ ##
# 1. Extract rpS3 gene
extr_rp(){
  
  # retrive sample name from file
    sample_name=$(ls $sample | cut -f1 -d".")
    faa_file=${sample}.genes.faa
    echo -e "... working on sample ${sample_name} ..."
    for domain in Bacteria Archaea Eukaryotes; do
      # hmmsearch (query assembled gene, against three HMM files)
      hmmsearch --tblout /dev/stdout -o /dev/null --cpu $cores --notextw ${hmm}/${domain}_90_trimmed.hmm ${faa_file} | grep -v "^#" > ${sample_name}.${domain}.hmm_results.txt
            # extract the 1st field
      awk '{print $1}' ${sample_name}.$domain.hmm_results.txt >> ${sample_name}.$domain.all_hits.txt
      # extract the faa sequence for the query with hitted hmm
      pullseq -i ${faa_file} -n ${sample_name}.$domain.all_hits.txt > ${sample_name}.$domain.all_hits.fasta
      #      ruby ${tblx}/remove_linebreaks_from_fasta.rb -f ${sample_name}.$domain.all_hits.fasta > ${sample_name}.$domain.all_hits_woLB.fasta
      # format the faa sequence (one sequence one line)
      seqtk seq ${sample_name}.$domain.all_hits.fasta |cut -f1 -d" " > ${sample_name}.$domain.all_hits_woLB.fasta
      
      ## get length 
      for gene in $(cat ${sample_name}.$domain.all_hits.txt);do
          grep -w -A 1 $gene ${sample_name}.$domain.all_hits_woLB.fasta | sed 1d | wc;
      done | awk '{print $3}' > ${sample_name}.$domain.length.tmp

      ## get hmmsearch score  
      awk '{print $6}' ${sample_name}.$domain.hmm_results.txt > ${sample_name}.$domain.scores.tmp

      ## paste the name, length and score into one file (containing 3 column)
      paste ${sample_name}.$domain.all_hits.txt ${sample_name}.$domain.scores.tmp ${sample_name}.$domain.length.tmp > ${sample_name}.${domain}_hits.txt
  done
  
  
 # filter tables by score and length
    # BACTERIA
  awk '{ if ($2 >= 111) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Bacteria_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 > ${sample_name}_final_rpS3.ids
    # ARCHAEA
  awk '{ if ($2 >= 172) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Archaea_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 >> ${sample_name}_final_rpS3.ids
    # EUKARYOTES
  awk '{ if ($2 >= 175) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected1
  awk '{ if ($3 <= 450) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected2
  awk '{ if ($3 >= 120) print $1 }' ${sample_name}.Eukaryotes_hits.txt > ${sample_name}_selected3
  grep -w -f ${sample_name}_selected1 ${sample_name}_selected2 > ${sample_name}_crit1
  grep -w -f ${sample_name}_crit1 ${sample_name}_selected3 >> ${sample_name}_final_rpS3.ids
  # retrieve final set of rpS3 genes for sample
  cat ${sample_name}_final_rpS3.ids | pullseq -i ${faa_file} -l 10000 -N > ${sample_name}_final_rpS3.fa

  # clean up
  rm ${sample_name}.*.all_hits.fasta ${sample_name}.*.all_hits_woLB.fasta ${sample_name}.*.scores.tmp ${sample_name}.*.length.tmp ${sample_name}*_hits.txt ${sample_name}_selected1 ${sample_name}_selected2 ${sample_name}_selected3 ${sample_name}_crit1 ${sample_name}_final_rpS3.ids 
  echo -e "... done with ${sample_name} ..."
}
## ++++++++++++++ END OF FUNCTIONS ++++++++++++++ ##
```

#### 1-4-2. Retrieve the contigs representing the rpS3 cluster 

Among the metagenomic contigs containing a rpS3 gene belonging to a rpS3 cluster the one which is assigned to a MAG with the largest completeness is selected as the representative sequence for each SG. When multiple contigs corresponding to the rpS3 cluster are retrieved in MAGs, the longest one is selected. 

```
#!/usr/bin/env perl
use strict;
use warnings;

# Assume that the assembly file (one-line fasta) in the current directory 

# get  a list of rpS3 cluster in a directory
#my $target_dir = shift @ARGV;

my $target_dir = "rpS3_clusters";
$target_dir =~ s|/$||;
my @target_cluster = glob "$target_dir/rpS3_*";

# rename each rpS3 cluster as "SG000001"
my $count = 0;

# get the reference (metagenome contig)
my @target_ref = glob "*.fa";

my $rpS3_cluster_scaf_out = "final_rpS3_cluster_scaffolds.fna";
my $rpS3_cluster_scaf_table_out = "final_rpS3_cluster_scaffolds.tsv";
my $rpS3_cluster_faa = "final_rpS3_cluster_protein.faa";

open(OUT_scaf, ">", $rpS3_cluster_scaf_out) or die $!;
open(OUT_gene, ">", $rpS3_cluster_faa) or die $!;
open(OUT_table, ">", $rpS3_cluster_scaf_table_out) or die $!;
print OUT_table join("\t", "id", "cluster", "gene", "scaffold"), "\n";

# 1st time: scan reference/assembly contigs of all samples to extract their length
my %seqs_len = ();
foreach my $rr (@target_ref) {
    open(IN, "<", $rr) or die $!;
    while (<IN>) {
        chomp;
        if (/>(\S+)/) {
            my $header= $1;
            my $seq = <IN>;
            chomp $seq;
            $seqs_len{$header} = length($seq);
        }
    }
    close IN;
}

# extract (1) the longest scaffolds contain rpS3 gene,
#         (2) their corresponding genes.
my %longest_scafs = ();

foreach my $ss (@target_cluster) {
    $count++;
    my $new_id = sprintf("SG%05d", $count);
    my ($cluster_id) = $ss=~/\S+\/(\S+)/;

    open(IN_S, "<", $ss) or die $!;
    my $longest_scaf_id = "";
    my $longest_gene_id = "";
    my $longest_scaf_len = 0;
    my $last_gene_id = "";
    my $last_gene_seq = '';
    my %genes = ();
    
    while (<IN_S>) {
        chomp;
        if (/^>(\S+)\s+.*/) {
            if ( $last_gene_seq ne '') {
                $genes{$last_gene_id} = $last_gene_seq;
                $last_gene_seq = '';
            }
            $last_gene_id = $1;
            my ($scaf_id) = $last_gene_id=~/(\S+)_\d+/;
            my $scaf_len = $seqs_len{$scaf_id};
            # compare the length of scaffold, and kept the longest scaffold
            if ( $longest_scaf_len < $scaf_len ) {
                $longest_scaf_len = $scaf_len;
                $longest_scaf_id = $scaf_id;
                $longest_gene_id = $last_gene_id;
            }
        } else {
            $last_gene_seq .= $_;
        }
        
    } # end of while, cluster file
    close IN_S;
    # last gene seq
    if ( $last_gene_seq ne '') {
        $genes{$last_gene_id} = $last_gene_seq;
        $last_gene_seq = '';
    }
    
    $longest_scafs{$longest_scaf_id} = 1;

    # export the representative scaffold sequence    
    print OUT_gene ">", $longest_gene_id, "\n";
    print OUT_gene $genes{$longest_gene_id}, "\n";

    # export the mapping table 
    print OUT_table join("\t", $new_id, $cluster_id, $longest_gene_id, $longest_scaf_id), "\n";
    
}
# export the representative scaffold sequence

# 2nd time: scan reference/assembly contigs of all samples to extract the sequence of target scaf
foreach my $rr (@target_ref) {
    open(IN, "<", $rr) or die $!;
    while (<IN>) {
        chomp;
        if (/>(\S+)/) {
            my $header= $1;
            if ($longest_scafs{$header}) {
                my $seq = <IN>;
                print OUT_scaf ">", $header, "\n";
                print OUT_scaf $seq;
            }
        }
    }
    close IN;
}
close OUT_gene;
close OUT_table;
close OUT_scaf;

```

#### 1-4-3. Estimate the abundance of each rpS3 cluster

To get the abundance of each SG across samples, the clean reads are aligned against its representative sequence using Bowtie2, then the mapped reads with ≥ 99% sequence identity is filtered and counted using the ‘depth’ module of Samtools.The abundance of each SG in a sample is calculated as the total mapped bases on the representative sequence divided by the length of representative sequence and the total number of sequencing bases in the sample. 


```

ref_seq="$1"
cores="$2"

bowtie2-build --threads $cores --quiet ${ref_seq} ${ref_seq}

bowtie2 --threads $cores  -x ${ref_seq} -1 $read -2 $read2 | samtools view -F4 -uT ${ref_seq} - | tee ${ref_seq}__${id}_mapped.bam | samtools view -T ${ref_seq} -h -O SAM -o ${ref_seq}__${id}_mapped.sam

```

The alignment in the SAM file is filtered using a custom Perl script.

```

#!/usr/bin/env perl
use strict;
use warnings;

## filter the reads according to percent identity
my $ff = shift @ARGV;

my $read_pid_threshold = 0.99;

my $out = $ff;
$out =~ s/\.sam/_filter.sam/;
open(IN, "<", $ff) or die $!;
open(OUT, ">", $out) or die $!;
while (<IN>) {
    if (/^\@HD/ or /^\@SQ/ or /^\@PG/) {
        print OUT;
        next;
    }
    chomp;

    my ($num_mismatch) = $_=~/NM:i:(\d+)/;
    my @sp = split /\t/;
    my $read_len = length($sp[9]);
    my $pid = ($read_len - $num_mismatch)/$read_len;
    if ($pid >= $read_pid_threshold) {
        print OUT $_, "\n";
    }
}
close IN;
close OUT;
```

The total base mapped on the representative sequence for each SG is calculated using the 'depth' module of Samtools.

```
samtools sort --reference $genome --threads $cores -O BAM $ff  | samtools depth - > $out

```

The abundance of each SG across samples is calculated using a custom Perl script.

```

#!/usr/bin/env perl
use strict;
use warnings;

my $ref_f = shift @ARGV;
my $sample_depth_f = "sample_sequencing_depth.tsv";
my $constant = 1000;

## The depth of each site >= than 'threshold' is considered to be confident.
my $base_threshold = 0;
my $covered_threshold = 0;

my $output_f = "abundance.tsv";

open(IN_depth, "<", $sample_depth_f) or die $!;
open(OUT, ">", $output_f) or die $!;

# open the sample depth file
my %sample_depth = ();
<IN_depth>;
while (<IN_depth>) {
    chomp;
    my ($sample_id, $read_bases, $read_counts) = split /\t/;
    $sample_depth{$sample_id} = $read_bases;
}
close IN_depth;

## open the reference file in FASTA, then get the total length of each sequence
my %genome_seq = read_fasta($ref_f);
my %genome_size = map { $_ => length($genome_seq{$_}) } keys %genome_seq;

## Estimate the abundance from the depth files for each sample
my @target_files = glob '*.depth';
my %genome_abund = ();
my %genome_cov = ();
my %samples = ();

foreach my $ff (@target_files) {
    my ($ss) = $ff=~/\S+?__(\S+?)_mapped_filter.depth/;
    $samples{$ss} = 1;

    open(IN_D, "<", $ff) or die $!;
    my $total_bases = 0;
    my $covered_sites = 0;
    my $last_genome2 = '';
    # read the depth file for each sample
    while (<IN_D>) {
        chomp;
        my ($contig, $site, $bases) = split /\t/;

        # only consider the base passing the threshold
        if ($bases >= $base_threshold) {
            $total_bases+=$bases;
            $covered_sites++;
        }           
        
        my $gg2 = $contig;
        if ($last_genome2 eq '') {
            $last_genome2 = $gg2;
        } elsif ($last_genome2 ne $gg2) {
            my $gg_size = $genome_size{$last_genome2};
            my $ss_depth = $sample_depth{$ss};
            my $abund = sprintf("%.2f", ($total_bases*$constant)/($gg_size*$ss_depth) );
            my $cov = sprintf("%.2f", $covered_sites/$gg_size);
            $genome_abund{ $last_genome2 }{ $ss } = $abund;
            $genome_cov{ $last_genome2 }{ $ss } = $cov;

            # re-initialize the count
            $total_bases = 0;
                        $covered_sites = 0;
            
            $last_genome2 = $gg2;
        }

    } 
    close IN_D;
    # the last sequence
    if ($total_bases != 0) {
        my $gg_size = $genome_size{$last_genome2};
        my $ss_depth = $sample_depth{$ss};
        my $abund = sprintf("%.2f", ($constant*$total_bases)/($gg_size*$ss_depth) );
        my $cov = sprintf("%.2f", $covered_sites/$gg_size);
        $genome_abund{ $last_genome2 }{ $ss } = $abund;
        $genome_cov{ $last_genome2 }{ $ss } = $cov;
    }
    
}

# Print the table containing the abundance of each genome across mutliple samples
my @header = sort(keys %genome_size);
print OUT join("\t", 'sample', @header), "\n";

foreach my $ss (sort(keys %samples)) {
    print OUT $ss;
    foreach my $gg (sort(keys %genome_size)) {
        if ($genome_abund{ $gg }{ $ss }) {
            my $cov2 = $genome_cov{ $gg }{ $ss };
            if ($cov2 >= $covered_threshold) {
                print OUT "\t", $genome_abund{ $gg }{ $ss };
            } else {
                print OUT "\t", 0;
            }
        } else {
                    # No reads from the sample are mapped against the genome
            print OUT "\t", 0;
        }
    }
    print OUT "\n";
}

sub read_fasta {
    # the sequnce is either one-line or mutlipel-line format
    my $seq_f = shift @_;
    open(IN, "<", $seq_f) or die $!;
    my %seqs = ();
    my $last_seq = '';
    my $last_header = '';
    while (<IN>) {
        chomp;
        if (/^>(\S+)/) {
            if ($last_seq ne '') {
                $seqs{$last_header} = $last_seq;
                $last_seq = '';
            }
            $last_header = $1;
        } else {
            $last_seq .= $_;
        }
    }
    # the last sequence 
    if ($last_seq ne '') {
        $seqs{$last_header} = $last_seq;
    }
    close IN;
    return %seqs;
}

```

#### 1-4-4. Assign the taxonomy to SGs

Before determining the taxonomy of each SG, a local custom database is established, which consisting of the amino acid of rpS3 genes derived from available genomes the RefSeq prokaryotic genome database (downloading date: 2019-07, ~27000 genomes). After that, the sequence of each SG retrieved from the target samples are aligned against the datatbase using BlastP program. The taxonomy of the top one hit with the percent identity ≥ 50 is assigned to the query SG. 

```

query="$1"
rpS3_ref_db="$2"
cores="$3"
out=${query}.blastp

blastp -query $query -db ${rpS3_ref_db} -outfmt 6 -out $out  -evalue 1e-3 -max_target_seqs 50 -num_threads $cores

```

The taxonomy of SGs is further validated by the clustering pattern of phylogenetic tree of the rpS3 gene. The amino acid sequence of all representative rpS3 genes are aligned using MAFFT, trimmed using trimA. A maximum likelihood phylogenetic tree is built by FastTree. Then, manually inspection of each SG on the tree could validate the accuracy of taxonomy based on sequence alignment. The taxonomy of the SGs that have no hits in our custome reference database and branch deeply in microbial lineages in the rpS3 tree were designated as ‘Unassigned’.

```

seq="final_rpS3_cluster_protein.faa"
cores="$1"

einsi --thread $cores $seq > "${seq}.mafft"
trimal -in "${seq}.mafft" -out "${seq}.mafft.trimAl" -automated1
FastTree -log "${seq}.fastree.log" "${seq}.mafft.trimAl" > "${seq}.fastree.nwk"

```

### 1-5. Screen the representative SG with MAGs


The representative genomes for a subset of SGs were identified through the rpS3-containning contig sequence shared between the reconstructed MAGs and the SGs. Noted that not all genome of SGs were recovered using metagenome binning methods partly due to a great diversity and insufficient sequencing coverage for the samples. 


Two input files are required for the R script below. The first file contains the SG ID and contig ID generated when extracting the contig representing the SG from the metagenomic assembly (See Section 4-2 above), and the second contains the headers of contig sequence for all samples extracted using the `grep` command in Linux.

```
#!/usr/bin/env Rscript 

require(dplyr)
require(tidyr)

# get the input files from command line
cmd <- commandArgs(trailingOnly = TRUE)
rps3.list.f <- cmd[1]
mag.2.contig.f <- cmd[2]

# read the input files
mag.2.contig <- read.table(mag.2.contig.f, header=F,as.is=T,sep="\t")
colnames(mag.2.contig) <- c('mag', 'scaffold')
rps3.2.contig <- read.table(rps3.list.f, header=T,as.is=T,sep="\t")

# merge two files based on shared scaffold
rps3.2.contig %>% 
  left_join(mag.2.contig) %>% 
  drop_na() ->  rps3.2.contig

# check the consistency between rps3 and MAG,
#     remove the MAGs mapping two rpS3 genes on different contigs
#     '1 rps3 -- 1 MAG'
mag.list <- rps3.2.contig$mag
mag.replicate <- mag.list[duplicated(mag.list)]
rps3.2.contig %>% 
  filter(! mag %in% mag.replicate) %>% 
  select(mag, id, cluster, gene, scaffold) -> rps3.2.contig

write.table(rps3.2.contig, file= "mag_2_sg.tsv",row.names=F,quote=F,sep="\t")

```

### 1-6. Estimate functional traits for the representative SG


The protein-coding gene sequence of each representative genome is predicted by Prodigal and then annotated against KEGG database  and eggNOG database. The assignment of KOs to each gene is achieved using KOFamScan, which performs a homology search against a database of hidden Markov models with precomputed score thresholds for specific KOs. The annotation against eggNOG were performed using eggNOG mapper, then the COGs assigned to each gene is retrieved from the results. 


The KOFamScan is run as below. Before running the code below, the path to KO profile and KO list are modified accordingly.
```

query="$1"
cores="$2"

profile=/apps/kofamscan/1.3.0/db/profiles/prokaryote.hal
ko=/apps/kofamscan/1.3.0/db/ko_list

out=${query}.KOfam
tmpdir=${query}.tmpdir

exec_annotation -o $out --cpu $cores -p $profile -k $ko --tmp-dir $tmpdir  $query
rm -fr $tmpdir

```

The eggNOG mapper is run as below. Before running the code belwo, the path to the database of eggNOG is modified accordingly.

```
query="$1"
cores="$2"
eggnog_data_db="/app/mapper-data/"

# one-step annotation
emapper.py -i $query -o ${query}.out --cpu $cores

```

## 2. Statistical analysis

### 2-1. Calculate the community-weighted means (CWM) of functional traits

```
prefix <- 'marine'

# get the abundance of non-red MAGs corresponding to rps3 clusters
otu.f <- paste0('../Data/', prefix, '_species_abundance.tsv')
mag.f <- paste0('../Data/', prefix,'_rps3_2_mag.tsv')

otu <- read.table(otu.f , header=T,as.is=T,sep="\t")
rps3_2_mag <- read.table(mag.f, header=T,as.is=T,sep="\t")

# extract the sub-table and rename each species as MAG ID
otu %>% 
  gather(key= 'id', value= 'abund', -sample) %>% 
  right_join(rps3_2_mag[, c('mag', 'id')]) %>% 
  select(-id) %>% 
  spread(key= 'mag', value= 'abund') -> otu.sub

# covert into data matrix with row name
rownames(otu.sub) <- otu.sub[, 1]
otu.sub[, 1] <- NULL

# filter OTU
otu.sub <- otu.sub[, colSums(otu.sub)>0]

# 1. get the KO copy number for each MAG
ko.f <- paste0("../Data/", prefix, "_ko_per_genome.tsv")
ko.mag <- read.table(ko.f, header=T, sep="\t", as.is=T)
rownames(ko.mag) <- ko.mag[, 1]
ko.mag[, 1] <- NULL

ko.mag <- ko.mag[colnames(otu.sub), ]

# filter KO (at least present in 10% of MAGs)
ko.threshold <- 0
ko.mag <- ko.mag[, colSums(ko.mag>0)>ko.threshold]

# get the CWM-KO for each sample

# check consistent names of two datasets
sum(rownames(ko.mag) != colnames(otu.sub))

cwm.ko <- weimea::cwm(otu.sub, traits = ko.mag)

# save the result
cwm.f <- paste0("../Data/", prefix, "_cwm-ko.RData")
save(cwm.ko, file= cwm.f)

```


### 2-2. Identify the key functional traits

The trait and environment replationship (TER), that is, the relationships between CWM values of functional traits (KO or COG items) and environmental parameters, are examined by Spearman rank correlation coefficients. The significance was tested with a row- and column-based permutation test to control inflated Type I error in the statistical analyses of TERs.



```
prefix <- 'marine'

# get the CWM KO
ko.f <- paste0("../Data/", prefix, "_cwm-ko.RData")
load(ko.f)

# import environment factor for samples
env.f <- paste0("../Data/", prefix, "_environment.tsv")
env <- read.table(env.f, header=T, as.is=T, sep="\t")

# test the TER relationship using the max test
#   (the row- and column-permutation test)
cwm.ko.2.env.cor <- test_cwm(cwm.ko, env, 
                             method= 'cor', test= 'max', perm= 499, parallel = 8)

cwm.ko.2.env.cor$out %>% 
  mutate(var.pair= rownames(.)) %>% 
  mutate(ko= sub("(K\\d+).*", "\\1", var.pair)) %>% 
  select(ko, r, starts_with('P_')) %>% 
  rename(cor= r, pvalue= P_max) %>% 
  mutate(group= prefix) -> cwm.ko.env.cor

p.thresh.cor <- 0.05

cwm.ko.env.cor %>% 
  filter(pvalue <= p.thresh.cor) %>% 
  mutate(direc= if_else(cor<0, 'negative', 'positive')) -> cwm.ko.env.cor.signif

# assign the full KEGG annotation for each KO
ko.annot <- read.table("../Tables/kegg_ko_full_annotation.tsv", 
                       header=T,as.is=T,sep="\t", quote='')
cwm.ko.env.cor.signif %>% 
  arrange(ko) %>% 
  left_join(ko.annot) -> cwm.ko.env.cor.signif


```

After retrieving a list of KO significantly associated with environmental factors (like lake alkalinity, soil temperature and marine salinity), we could easily pinpoint the underlying functional traits that closely associated with their environmental distribution and adaptation of microbial species across gradients.


### 2-3. Calculate the index of environmental adaptation based on TER

We proposed An index of environmental adaptation (iEA) to quantify the adaptive capacity of a species to environmental stress by considering its key traits and their relationships with environments. The iEA is calculated as the mean value of the TER across the species’ traits weighted by the copy number of KOs.


Positive and negative iEA values indicate that a microbial species is dominated by alkaline-tolerant and alkaline-sensitive traits, respectively, while an iEA of zero suggests the two kinds of traits equally dominate 


```
prefix <- 'marine'

#  get the KO copy number for each MAG
ko.f <- paste0("../Data/", prefix, "_ko_per_genome.tsv")
ko.mag <- read.table(ko.f, header=T, sep="\t", as.is=T)

## TER (CWM-KO ~ gradient)
cor.f <- paste0("../Tables/", prefix, "_kegg_ko_2_env_corr_max_text.tsv")
cwm.ko.env.cor <- read.delim(cor.f, header=T, sep="\t", as.is=T)

p.thresh.cor <- 0.05
# using the KOs that significantly related to pH
cwm.ko.env.cor %>% 
  drop_na() %>% 
  rename(cor= r, p.value= P_max) %>% 
  filter(p.value < p.thresh.cor) %>% 
  select(ko, cor)  -> ko.env.index

# calculate the index considering the effect of the KO copy number
ko.mag %>% 
  gather(key= 'ko', value= 'copy', -id) %>% 
  right_join(ko.env.index) %>% 
  filter(copy !=0) %>% 
  group_by(id) %>%
  summarise(adapt.index= sum(cor*copy)/sum(copy)) -> mag.2.env.adapt


```


### 2-4. Fit the relationships of the species' index and its niche optima
  
The extent of adaptation index to the species optimum niche was accessed by the relationship between the index of a species and its pH optima through a linear regression model. The niche optima of an individual (i.e., pH or salinity) could reflect the optimum niche where the species is expected to be most abundant, and was calculated based on its abundance across samples and lake pH of samples where the species was identified.

```

## get optimal environment for each MAG
niche.f <- paste0("../Tables/", prefix, "_mag_env_optima.tsv")
mag.2.env <- read.table(niche.f, header=T, sep="\t", as.is=T)
mag.2.env %>% 
  select(mag, env.w.abund) %>% 
  rename(env= env.w.abund) -> mag.2.env

# merge KO adaptation index and niche optima together
mag.2.env.adapt %>%
  rename(mag= id) %>%
  right_join(mag.2.env) -> mag.2.env.adapt2

# fit the relationship
mag.2.env.adapt2 %>% 
  do(glance(lm(.$adapt.index ~ .$env))) -> mag.env.lm

p.thresh.lm <- 0.05
mag.env.lm %>% 
  mutate(signif= if_else(p.value < p.thresh.lm,
                         'signif', 'insignif')) %>%
  mutate(adj.r.squared = round(adj.r.squared, 2)) %>% 
  mutate(data= prefix)  -> mag.env.lm

mag.2.env.adapt2 %>% 
  mutate(data= prefix) %>% 
  left_join(mag.env.lm) %>%
  drop_na() -> mag.2.env.adapt2

## show the relationship using scatter
my.linetype <- c('signif'=1, 'insignif' = 2)

mag.2.env.adapt %>% 
  ggplot(aes(x= env, y= adapt.index)) +
  geom_point(size= 4, pch=21, fill= 'grey', stroke=0.1) +
  # add the fitted line
  geom_smooth(aes(linetype= signif),
              method= 'lm', se= T, size= 1) +
  scale_linetype_manual(name= 'Significance',
                        values= my.linetype) +
  
  labs(x= 'Niche optima', y= 'Adaptation index') +
  facet_wrap(.~direc, scales= 'free_y', nrow= 1) + 
  theme(legend.position = 'none', panel.grid.minor = element_blank())


```


You can use the [editor on GitHub](https://github.com/mingleiR/test/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

