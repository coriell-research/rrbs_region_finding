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

#### methclone

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
[kkeith]$ tmux attach -t rrbs_regions
[kkeith]$ pwd
/home/kkeith/data/rrbs_region_finding/cgmaptools
[kkeith]$ /usr/local/programs/rrbs_region_finder/cgmaptools/cgmaptools convert bam2cgmap -b ../data/M13_sorted.bam -g /mnt/data/gdata/human/hg19_hg37/hg19/hg19_soft_masked/hg19.fa -o cgmaptools_M13 --rmOverlap
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

