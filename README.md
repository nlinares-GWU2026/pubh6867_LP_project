# Tracking a European Lactase Persistence Variant Through Ancient DNA 

Final project for PUBH 6867 (Health Data Visualization), Nicole Linares, George Washington University.

This project uses published ancient DNA to follow one lactase persistence (LP) variant common in Europe, **rs4988235** (delivered allele **A**), from before farming to today. It asks two questions:

1. How did the frequency of the LP allele change across Europe, and did the change look the same in different regions?
2. Does the timing of the rise line up with the arrival of steppe-associated groups?

All figures and tables come from one script, `scripts/project_data_plot.R`. The data are from the Allen Ancient DNA Resource (AADR) v66.p1.

## What is in this repository

```
scripts/project_data_plot.R          All cleaning, analysis, figures, and tables (commented, run top to bottom)
data/lactase_extract_raw.txt         The one-SNP genotype extract (23,089) - included, so you do not need the 6.7 GB genotype file
figures/                             Figures used in the report (PNG)
figures/slides/                      Same figures re-exported at slide size
tables/                              Tables used in the report appendix (CSV)
environment.yml                      conda environment with PLINK and EIGENSOFT (only needed to rebuild the extract)
renv.lock, .Rprofile, renv/          R package versions (R 4.6.1)
Project_HDV.Rproj                    RStudio project (open this to work in the project root)
notes.md                             Running project notes, decisions, and instructor-feedback changes
AI_USAGE.md                          Documentation of how AI was used
```

The large AADR files are **not** in the repository (see `.gitignore`). They are free to download. 

## Two ways to reproduce

** Option A (recommended, quicker): Reproduce every figure and table from the included extract.**
You only need R and the AADR annotation file.

**Option B (optional, slower): Rebuild the extract from the raw AADR genotype files.**
Only needed if you want to check how `data/lactase_extract_raw.txt` was made. 

### Step 1. Download the AADR files (Harvard Dataverse)

AADR is maintained by the Reich Laboratory and is free, with no access request: <https://doi.org/10.7910/DVN/FFIDCW>

On the dataset page, download these files (the file names below are the ones this project used):

| File | Needed for |
|---|---|
| `v66.p1_1240K.aadr.PUB.anno` | **Options A and B** (dates, places, group labels) |
| `v66.p1_1240K.aadr.patch.PUB.geno` (about 6.7 GB) | Option B only |
| `v66.p1_1240K.aadr.patch.PUB.snp` | Option B only |
| `v66.p1_1240K.aadr.patch.PUB.ind` | Option B only |

Put `v66.p1_1240K.aadr.PUB.anno` in this repository's `data/` folder, so its path is `data/v66.p1_1240K.aadr.PUB.anno`. The script reads exactly that path. 

### Step 2 (Option B only). Rebuild the genotype extract

Run in a bash terminal (this project used WSL2 Ubuntu). Put the downloaded `.geno`, `.snp`, and `.ind` files in one working folder outside the repository.

```bash
# 1. Create the environment (PLINK 1.9 and EIGENSOFT)
conda env create -f environment.yml
conda activate aadr-lp
 
# 2. Tell convertf where the files are and where to write PLINK-format output
cat > convertf_lactase.par << 'EOF'
genotypename: v66.p1_1240K.aadr.patch.PUB.geno
snpname: v66.p1_1240K.aadr.patch.PUB.snp
indivname: v66.p1_1240K.aadr.patch.PUB.ind
outputformat: PACKEDPED
genotypeoutname: aadr_1240k.bed
snpoutname: aadr_1240k.bim
indivoutname: aadr_1240k.fam
hashcheck: NO
EOF

# 3. Convert the whole dataset to PLINK format (longer step; reads the 6.7 GB file)
convertf -p convertf_lactase.par 2>&1 | tee convertf_log.txt

# 4. Check the allele order: expect A1 = G (ancestral) and A2 = A (derived, LP allele)
grep -w "rs4988235" aadr_1240k.bim

# 5. Extract the one SNP. --keep-allele-order stops PLINK from silently swapping A1/A2.
plink --bfile aadr_1240k --snp rs4988235 --keep-allele-order --recode A --out lactase_extract

# 6. Sanity checks: 23,090 lines (header + 23,089); first row Loschbour.AG should have rs4988235_G = 2
wc -l lactase_extract.raw
head -3 lactase_extract.raw

# 7. Copy the result into the repository under the name the script expects
cp lactase_extract.raw <path-to-repository>/data/lactase_extract_raw.txt
```

Note on `convertf`: this project ran `convertf` from **AdmixTools built from source** (<https://github.com/DReichLab/AdmixTools>). The `eigensoft` package in `environment.yml` also provides a `convertf`. Use whichever one can read an AADR `.geno` file on your machine.

### Step 3. Run the analysis in R 

1. Open `Project_HDV.Rproj` in RStudio. This sets the working directory to the project root, which the script requires. The `.Rprofile` activates `renv` automatically.
2. Restore the R packages (versions are recorded in `renv.lock`. The project used is R 4.6.1):
```r
renv::restore()
```
3. Run the whole script from the project root:
```r
source("scripts/project_data_plot.R")
```

The script prints checks along the way (counts, sensitivity checks) and writes the files below. It creates the `figures/`, `figures/slides/`, and `tables/` folders if they do not exist.

## What the script produces

| Output | Report location | Script section | 
|---|---|---|
| `figures/q1_lactase_timeseries.png` | Figure 1 | Q1: time series |
| `figures/q1_regions.png` | Figure 2 | Q1: regional check |
| `figures/q2_steppe_timing.png` | Figure 3 | Q2: timing (steppe window) |
| `figures/q2_steppe_comparison.png` | Figure 4 | Q2: bar chart |
| `figures/q1_map.png` | Figure A1 | Q1: geographical map |
| `tables/regional_leave_one_out.csv` | Table A1 | Q1: regional check |
| `tables/steppe_restricted_comparison.csv` | Table A2 | Q2: instructor feedback check |
| `figures/slides/*.png` | Presentation slides | End of each figure section |

## Method Quick Summary

- **Genotype:** PLINK's `rs4988235_G` column counts the ancestral G allele. The script converts it to the derived LP allele with `A_count = 2 - G_count`.
- **One read per person:** ancient DNA is usually low coverage, so most ancient people contribute a single allele read (recorded as 0 or 2 copies of A). The analysis counts *people*, not DNA copies. The 95% ranges are Wilson intervals on the number of people.
- **Europe:** 40 countries (Russia and Turkey excluded). 7,711 people in total: 7,150 archaeological and 561 present-day.
- **Time:** 1,000-year periods in years before 1950 (BP). Everyone older than 10,000 BP is combined into one single pre-Neolithic group.
- **Steppe-associated:** AADR group labels: Yamnaya, Catacomb, and Corded Ware (129 people). Bell Beaker is kept separate (Olalde et al., 2018). Labels are archaeological names, not measured ancestry.
- **Checks:** Belgium sensitivity check for the most recent period, region-and-period checks without the largest country, and a steppe comparison restricted to countries with both groups.

 ## Citing the Data

 If you use AADR data, cite all three of the following (AADR's request):

 1. Mallick, S., Micco, A., Mah, M., Ringbauer, H., Lazaridis, I., Olalde, I., Patterson, N., & Reich, D. (2024). The Allen Ancient DNA Resource (AADR) a curated compendium of ancient human genomes. *Scientific Data*, 11(1), 182. <https://doi.org/10.1038/s41597-024-03031-7>
2. Mallick, S., & Reich, D. (2023). The Allen Ancient DNA Resource (AADR): A curated compendium of ancient human genomes. Harvard Dataverse, V14. <https://doi.org/10.7910/DVN/FFIDCW>
3. The original publication for each individual used. AADR compiles previously published data, and the source study for each person is recorded in the `.anno` file.
## Use of AI
 
AI assistance is documented in `AI_USAGE.md`.
