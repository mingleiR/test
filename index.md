


## 1. Sequence analysis

1. Filter and trim metagenomic reads

1-1. detect the potential adapter and contamination


```
fastq_f=$1
cores=$2

## reading fastq file in the command
fastqc --format fastq  --extract --threads $cores ${fastq_f}
```


1-2. trim the adapter and cotamination


```
sample="$1"
cores="$2"

r1=${sample}_1.fastq.gz
r2=${sample}_2.fastq.gz

# minimum length for trimmed reads (50 is okay for assembly; 100 for SNP identification)
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



## Welcome to GitHub Pages test [test]

You can use the [editor on GitHub](https://github.com/mingleiR/test/edit/gh-pages/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/mingleiR/test/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
