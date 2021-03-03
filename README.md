## Finding Regions of Methylation in RRBS

Because we want to do entropy calculations in RRBS, we need a way of finding CpGs on the same read. We have already used `methclone` for this purpose, but want to test other existing programs as well. `methclone` only returns CpGs in groups of 4 and we would like variable amounts of CpGs to be identified as well.

### Setup

```bash
[kkeith]$ cd /mnt/data/data_kk
[kkeith]$ mkdir rrbs_region_finding
[kkeith]$ cd rrbs_region_finding
[kkeith]$ mkdir data
[kkeith]$ cd data
### link existing data to use in testing
[kkeith]$ ln -s /home/jjelinek/data/jj4/rrbs/data/rrbs35_novogene aging_data_batch1
[kkeith]$ ln -s /mnt/data/research_data/2019-12-19_novogene_rrbs46_aging/ aging_data_batch2
[kkeith]$ ln -s /mnt/data/data_kk/aging_rrbs/cord_blood_rrbs/ cord_blood_data
```
#### Reorganizing
*2021-03-03*

Put tools in their own folder and make into a git repository

```bash
[kkeith]$ cd /home/kkeith/data/rrbs_region_finding
[kkeith]$ mkdir tools
[kkeith]$ mv cgmaptools/ methclone_M13.txt.gz tools

### set up
[kkeith]$ git init
[kkeith]$ git add README.md
[kkeith]$ nano .gitignore
[kkeith]$ more .gitignore 
### folders to ignore
data/
tools

### files to ignore
.DS_Store
[kkeith]$ git add .gitignore
[kkeith]$ git commit -m 'first commit'
[kkeith]$ git branch -M main
[kkeith]$ git remote add origin git@github.com:coriell-research/rrbs_region_finding.git
[kkeith]$ git push -u origin main

### rearrange data folder because I want actual RRBS and amplicon data and index amplicon BAM file
[kkeith]$ cd data
[kkeith]$ [kkeith]$ mkdir rrbs
[kkeith]$ mv aging_data_batch* cord_blood_data rrbs/
[kkeith]$ mkdir amplicons
[kkeith]$ ln -s /mnt/data/research_data/2020-11-17_bis_amplicons_woonbok/processed_data/03_align/hg19/ amplicons/amps_round2_bams
[kkeith]$ mv M13_sorted.bam M13_rrbs_sorted.bam
[kkeith]$ mv M13_sorted.bam.bai M13_rrbs_sorted.bam.bai
[kkeith]$ samtools sort -t 8 amplicons/amps_round2_bams/M13_CKDL200165722-1a-AK845-GD06_HF5LTCCX2_L1_1_val_1_bismark_bt2_pe.bam -O BAM -o M13_amplicon_sorted.bam
[kkeith]$ samtools index -@ 4 M13_amplicon_sorted.bam
```

---

### Installing Programs

```bash
### 
[kkeith]$ cd /usr/local/programs
[kkeith]$ mkdir rrbs_region_finder/
[kkeith]$ cd rrbs_region_finder/
[kkeith]$ cd ~
### methpipe
[kkeith]$ git clone https://github.com/smithlabcode/methpipe.git
[kkeith]$ sudo mv methpipe /usr/local/programs/rrbs_region_finder/
[kkeith]$ cd /usr/local/programs/rrbs_region_finder/cgmaptools
[kkeith]$ sudo bash install.sh
### cgmaptools
[kkeith]$ git clone https://github.com/guoweilong/cgmaptools.git
[kkeith]$ sudo mv cgmaptools/ /usr/local/programs/rrbs_region_finder/
### methtuple
[kkeith]$ pip install methtuple
```

---

### Testing Programs

Try with one sample first

#### Set up tool comparison file
*2021-03-03*

```bash
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding
[kkeith]$ nano tool_features.txt
[kkeith]$ more tool_features.txt 
### File to list features of RRBS region finding tools as I test them
# tool = tool name
# multithread_available = Can the tool be multithreaded? no/yes/maybe, maybe for R programs that could
 maybe be made parallel
# number_steps = number of steps the program requires to give me a list of cpgs; does not include runn
ing my allele counting script or making a whitelist
tool	multithread_available	number_steps
methclone	no	1
cgmaptools	no	3
```

#### methclone

**OLD** Re-ran with new version and tracking time

```bash
[kkeith]$ tmux new -s rrbs_regions
[kkeith]$ pwd
mnt/data/data_kk/rrbs_region_finding
# methclone requires a sorted, indexed BAM file, so do that for my chosen file first
[kkeith]$ cd data
[kkeith]$ samtools sort -t 8 aging_data_batch1/bismark20hg19_pe/M13_CKDL190139057-1a-20_H3WM3BBXX_L1_1_val_1_bismark_bt2_pe.bam > M13_sorted.bam
[kkeith]$ samtools index -@ 4 M13_sorted.bam

### run methclone
[kkeith]$ methclone data/M13_sorted.bam data/M13_sorted.bam methclone_M13.txt.gz M13
```
<br>

*2021-03-03* **CORRECT** Re-ran with updated version of methClone, keeping track of how long it takes to run, and running on both the RRBS and amplicon version of the sample

```bash
[kkeith]$ cd /home/kkeith/data/rrbs_region_finding/tools
[kkeith]$ rm methclone_M13.txt.gz
[kkeith]$ mkdir methclone
[kkeith]$ cd methclone
[kkeith]$ tmux new -s regions
[kkeith]$ for i in ../../data/*.bam; do j=${i##*/}; echo -e "methclone\tstart\t$(date +'%F')\t$(date +'%T')\t${j/_sorted.bam/}" >> methclone_time.txt; methClone --f1=$i -s ${j/_sorted.bam/} -o ${j/_sorted.bam/}_methclone.txt.gz; echo -e "methclone\tend\t$(date +'%F')\t$(date +'%T')\t${j/_sorted.bam/}" >> methclone_time.txt; done
[detached (from session regions)]
```
<br>

### CGMapTools

```bash
[kkeith]$ tmux attach -t rrbs_regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding
[kkeith]$ mkdir cgmaptools
[kkeith]$ cd cgmaptools
```
**Step 1:** Convert BAM to CGMap file

```
### OLD, not this coode
[kkeith]$ tmux attach -t rrbs_regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding/cgmaptools
[kkeith]$ /usr/local/programs/rrbs_region_finder/cgmaptools/cgmaptools convert bam2cgmap -b ../data/M13_sorted.bam -g /mnt/data/gdata/human/hg19_hg37/hg19/hg19_soft_masked/hg19.fa -o cgmaptools_M13 --rmOverlap

### CORRECT
[kkeith]$ tmux attach -t regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding/tools/methclone
[kkeith]$ cd ../cgmaptools/

for i in ../../data/*.bam; do j=${i##*/}; echo -e "cgmaptools_convert\tstart\t$(date +'%F')\t$(date +'%T')\t${j/_sorted.bam/}" >> cgmaptools_time.txt; /usr/local/programs/rrbs_region_finder/cgmaptools/cgmaptools convert bam2cgmap --rmOverlap -b $i -g /mnt/data/gdata/human/hg19_hg37/hg19/hg19_soft_masked/hg19.fa -o 01_cgmaptools_${j/_sorted.bam/}; echo -e "cgmaptools_convert\tend\t$(date +'%F')\t$(date +'%T')\t${j/_sorted.bam/}" >> cgmaptools_time.txt; done
```
**Step 2:** Call SNPs; `CGmaptools` recommends the Bayesion mode (using flags `-m bayes --bayes-dynamicP`) since it's more accurate than the default binary mode. However since the documentation says binary mode is faster and I'm just trying to determine if the output of the tool is usable, I went with that.

```bash
[kkeith]$ tmux attach -t rrbs_regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding/cgmaptools
### rename existing files so I know what step they came from
[kkeith]$ for i in cgmaptools*; do mv $i 01_$i; done
### get VCF file
[kkeith]$ /usr/local/programs/rrbs_region_finder/cgmaptools/cgmaptools snv -i 01_cgmaptools_M13.ATCGmap.gz -v 02_cgmaptools_M13.vcf
```
**Step 3:** Call allele-specific methylation

```bash
[kkeith]$ tmux attach -t rrbs_regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding/cgmaptools
[kkeith]$ 

/usr/local/programs/rrbs_region_finder/cgmaptools/cgmaptools  asm -r /mnt/data/gdata/human/hg19_hg37/hg19/hg19_soft_masked/hg19.fa -b ../data/M13_sorted.bam -l 02_cgmaptools_M13.vcf -o 03_cgmaptools_M13_asmr.txt -t CG

```

