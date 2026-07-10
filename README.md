# MErlin
_**M**ethylation-driven **E**xpression & **R**egulation **L**inkage in **I**nteracting **N**uclear-domains_

How is the repo structured?

```
  📁 v0.2-270326                             # version up to (the last) date (27th march 2026)
    > 📁 dev                                 # folder with all the scripts (development part)
        📄 bin.py                            # script for genomic binning + methylation count per each bin
        📄 mbr.py                            # script to parse and filter a bedMethyl file produced by modkit (ONT) 
        📄 multiomics_ml.py                  # script to apply a SHAP framework to our data (multi-omics)
        📄 xam.py                            # script to cross the methylation basecalling with genomic annotation
    > 📁 res                                 # folder with the results from bin.py and xam.py
        📁 bin_res
          📄 DCK546_basemods_w5000_5mC.bed   # BED file containing 5mC binned at 5Kbp
          📄 DCK546_basemods_w5000_m4C.bed   # BED file containing 4mC binned at 5Kbp
          📄 DCK546_basemods_w5000_m6A.bed   # BED file containing m6A binned at 5Kbp
        📁 xam_res
          📄 DCK546_basemods_w5000_5mC.tsv   # TSV file containing 5mC methylation basecalling
          📄 DCK546_basemods_w5000_m4C.tsv   # TSV file containing m4C methylation basecalling
          📄 DCK546_basemods_w5000_m6A.tsv   # TSV file containing m6A methylation basecalling
    > 📁 test                                # folder containg the data to be tested
        📁 data                              # data for mbr.py, xam.py and bin.py
        📁 shap-data                         # data for multiomics_ml.py
  ```

> [!IMPORTANT]
> The code was built using the (huge) ((thank you)) support of Claude Code.

## What's next?
Take a look to the data and check the script for a quick run. We can then meet up to align perhaps?


# Step by step

## 1 - Cross-reference the methylation basecalling (BED or GFF3 format) with genomic annotation
```
python3 xam.py -g SmBL225C.gff3 -m F5.gff3 -o results_xam —no-split
```
## 2 - Differential methylation distribution 
```
python3 dixam.py -g SmBL225C.gff3 -m F5.gff3 -M F9.gff3 -o results_dixam —no-split 
```
## 3 - Spatial distribution (OriC-Ter)
```
python3 spacem.py --gff SmBL225C.gff3 --fasta SmBL225C_reference.fasta --methylation F5.gff3 --motifs motifs.txt --outdir results_spacem --prefix F5
```
## 4 - Motifs analysis
```
python3 mematch.py --gff SmBL225C.gff3 --methylation F5.gff3 --motifs motifs.txt --fasta SmBL225C_reference.fasta --outdir results_mematch --prefix F5
```
## 5 - Superfancy image
```
to debug, something's wrong..
```
## 6 - Differential Expression and Methylation (DiEM) analysis
```
python3 diem.py --geneexp DESeq2_results_all_BL.csv --methylation SmBL225C_methylation_by_feature.tsv --outdir results_diem
```
