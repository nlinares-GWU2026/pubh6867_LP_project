## Entry 1: 
- **Tool:** Claude (Sonnet 5)
- **What it was used for:** It was used to help me pull and clean the AADR data, run PLINK conversion through WSL2, decode how to interpret the converted PLINK file, clean, merge, subset, and time-bin the data (referencing literature for time-binning), section out modern populations, restructure bins based on the results of the presence of LP based on the pre 10000 BP populations (would disrupt the signal), validate data aggregation, and confirm alignment with current literature.
- **Output:**
- Step 3 — write the convertf parameter file. Note this converts the entire dataset (all ~1.24 million SNPs) to PLINK format — convertf itself can't pull out just one SNP, only PLINK can do that filtering, in a step after this one. It's a one-time cost: once you have this PLINK-format copy, pulling out rs4988235 (or any other SNP you want later) becomes fast and cheap.
```bash
cd ~/data/aadr
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
```
- Now let's actually run it:
```bash
time /home/nelin/AdmixTools/bin/convertf -p convertf_lactase.par 2>&1 | tee convertf_log.txt
```
- A few notes on what this does and what to expect:

2>&1 | tee convertf_log.txt saves everything convertf prints (including any warnings) to a log file while still showing it on screen — useful since this produces a fair amount of output and you'll want a record of it if anything looks off.
time in front just reports how long the whole thing took once it's done.
This is converting the entire dataset (all ~1.24 million SNPs, ~20,000 individuals) from the packed AADR format into PLINK format, so it will take real time to run — possibly several minutes, possibly longer depending on your machine, since it's reading straight through a 6.7 GB file. That's normal; you don't need to interrupt it unless it's been sitting for a very long time with zero output or it exits with an error.
-Before pulling genotypes out, there's one more thing worth confirming: which allele PLINK will treat as "A1" in the extracted file, since that determines how we read the resulting numbers. Let's check that alongside activating the right environment.
```bash
conda activate popgen
plink --version
```
Then check how this SNP is recorded in the converted .bim file:
```bash
grep -w "rs4988235" aadr_1240k.bim
```
A .bim file has six columns: chromosome, SNP name, genetic position, physical position, allele 1 (A1), allele 2 (A2). We already know from the AADR .snp file that G is the reference/ancestral allele and A is the derived, lactase-persistence allele — I want to see whether PLINK's A1/A2 assignment here lines up with that G/A ordering, or flipped it. That matters because the next step (--recode A) counts copies of whichever allele is A1, and if it silently flipped orientation, your frequency trend would read backwards (declining instead of rising) without any error message to warn you.
Paste back both outputs and we'll write the extraction command against whichever allele is actually A1.
- Good — this confirms convertf preserved orientation cleanly: A1 = G (ancestral) and A2 = A (derived, lactase-persistence allele), matching exactly what the AADR .snp file told us. No flip happened during conversion.

That said, there's one more thing worth locking down before we extract, rather than trusting it'll stay that way by default: PLINK's --recode A counts copies of whatever is A1 — and by default, PLINK is technically allowed to re-sort A1/A2 by allele frequency during certain operations. We don't want that happening silently on this one SNP, since it would flip which allele the counts refer to. The fix is one flag: --keep-allele-order, which pins A1/A2 to exactly what's already in the .bim file — no surprises. 
```bash
plink --bfile aadr_1240k --snp rs4988235 --keep-allele-order --recode A --out lactase_extract
```
What this does: --bfile aadr_1240k points at your converted dataset, --snp rs4988235 pulls out just this one SNP (fast — no need to touch the full 1.2 million), and --recode A writes a simple text file with one row per individual and a genotype count column.

What to expect in the output: a file called lactase_extract.raw with columns FID IID PAT MAT SEX PHENOTYPE rs4988235_G — that last column name literally tells you it's counting G (the ancestral allele), with values 0, 1, or 2 per person (or NA for missing data). Since we care about the derived, lactase-persistence allele (A), you'll need one small transformation once this is in R: A_count = 2 - G_count per person. Worth writing that as a clearly commented line in your R script rather than doing it mentally — it's exactly the kind of silent-sign-flip mistake that's easy to introduce by accident later.
Run that and paste back: 
```bash
head -5 lactase_extract.raw
wc -l lactase_extract.raw
```
That'll confirm the column name matches what we expect and show roughly how many individuals actually have a called genotype here (missingness is common in ancient DNA, so don't be surprised if it's well under 23,089). Once that looks right, we'll bring it into R alongside the .anno file to get dates and regions attached.
- 5) This looks genuinely good — and there's a really nice sanity check hiding right in the head output, which is worth pointing out.
The numbers check out: 23,089 people extracted, matching your .fam file exactly; wc -l gives 23,090 (header + 23,089 data rows), consistent. One variant, as requested.
The sanity check: Look at those first four samples — Loschbour and the "Kou"/"Log" individuals. Loschbour in particular is a well-known Mesolithic hunter-gatherer genome from Luxembourg, and every one of those early samples shows rs4988235_G = 2 — homozygous G/G, meaning zero copies of the derived, lactase-persistence allele. That's exactly what the literature we discussed predicted: pre-Neolithic and early Neolithic Europeans essentially lacked this variant entirely. Seeing that pattern show up correctly in your very first few rows, unprompted, is a strong sign the whole chain — AADR → convertf → PLINK — preserved the data correctly rather than just producing well-formed but scrambled output.
One number worth noting, not worrying about: Total genotyping rate is 0.672311 — meaning about 67% of the 23,089 individuals actually have a called genotype here, and the rest are missing. That's normal for a single SNP in ancient DNA (coverage varies a lot sample to sample) — not a red flag — but it's a real number you'll want to state in your methods (something like "genotypes were available for N of 23,089 individuals at this locus") and you'll need to drop the missing ones in R rather than let them silently become zeros.
A couple of housekeeping notes: the PAT, MAT, SEX, and PHENOTYPE columns in that .raw file are just placeholder/default values PLINK always includes — not real data — so ignore them; you only need IID and rs4988235_G. And remember the flip we discussed: you'll compute A_count = 2 - G_count in R to get counts of the actual lactase-persistence allele.
Before we jump into R, there's one thing worth a 30-second check now rather than debugging a failed join later: does the IID format here (like Loschbour.AG) match how individuals are named in the .anno file? AADR sometimes uses a slightly different ID column or suffix convention between files. Quick check:
```bash
awk -F'\t' '$1=="Loschbour.AG"' v66.p1_1240K.aadr.PUB.anno 
```
If that returns a matching ID, we're clear to bring both files into R and join them. Paste what that shows and we'll move on — and you're very welcome, glad the explanations are landing well.
- That merge is exactly right — 15,523 rows (matching raw_clean precisely) and zero unmatched group_ids. Every genotyped individual found its metadata. That's a genuinely major milestone: you now have one table with genotype, date, location, and culture for every sample. Everything downstream is analysis, not data wrangling.
Before building the actual frequency-by-time/region result, there are two real decisions to make together — not mechanical steps, actual choices that affect what your results mean:
1. How to define "Europe" for filtering. Your merged table only has political_entity (country) and coordinates, no continent field. The simplest approach is filtering by a list of European country names from political_entity — readable and easy to justify in your methods. The catch: a few countries straddle continents in ways that need an explicit call (Russia, Turkey especially) — you'd want a short, stated rule like "Russian and Turkish samples are excluded" or "included only if west of X longitude," rather than an implicit one.
2. How to bin time into actual archaeological periods. date_mean_bp is just a continuous number of years — nothing currently sorts a sample into "Neolithic" vs. "Bronze Age" etc. This needs real boundary years, and rather than inventing round numbers, it's worth grounding them in the same ancient-DNA literature you've already been citing (Mathieson et al., Allentoft et al.), since those papers define the period cutoffs actually used in this specific research area — and citing that source will make your methods section more defensible.
Let me pull the actual boundary years those papers use before we lock anything in, rather than eyeballing it.Good news — the actual paper on this exact SNP (Segurel et al. 2020, which you'll already be citing) gives us a ready-made, directly citable convention rather than needing to invent one: they bin ancient samples into 1,000-year BP intervals (e.g., "the 5,000–4,000 BP period") rather than trying to pin down exact start/end years for "Neolithic" vs. "Bronze Age" by region. That's a real advantage for you, because archaeological period boundaries genuinely differ across Europe — the Bronze Age starts in different centuries in the Aegean vs. Scandinavia — so tying your bins to fixed calendar-year periods would require picking a specific regional chronology and defending it. Round 1,000-year BP bins sidestep that entirely and let you cite precedent directly: "time bins follow the convention used in Segurel et al. (2020)."
Two things this suggests for how to build it:
Bin as: floor(date_mean_bp / 1000) * 1000 in R — clean, defensible, and directly comparable to that paper's own reporting.
For the non-specialist audience the rubric wants, pure BP labels ("5000–4000 BP") won't mean much on sight — worth loosely annotating a couple of bins with era names for orientation (e.g., "~8000 BP (early farming spreads into Europe)"), without claiming those era boundaries are precise or universal.
The "what counts as Europe" question turned out to have a wrinkle worth flagging now rather than after you've built it: your Q1 (frequency rise within Europe) and Q2 (does it track steppe-culture arrival) may actually need slightly different geographic scopes. The Yamnaya steppe heartland — the source population for the "steppe migration" story — sits in what's now Russia, which most people would draw outside "Europe" proper. So a strict Europe-only filter is right for the Q1 trend line, but Q2 might need you to also pull in a few explicitly steppe-associated groups (by group_id, regardless of country) as a comparison point, rather than relying on the country filter alone.
For the country filter itself: simplest and most defensible is an explicit list of European country names matched against political_entity, with Russia and Turkey deliberately excluded and that exclusion stated outright in your methods, rather than left implicit.
Does that framing make sense — 1,000-year BP bins for the time axis, and a stated European-country list (minus Russia/Turkey) for Q1, with steppe-culture groups pulled in separately by name for Q2? Or would you rather simplify further and use named era bins instead of raw BP numbers, accepting the regional-boundary fuzziness as a stated limitation?
**MY OWN RESPONSE** *I initially wanted to include Russia and Turkey because they technically have parts of Europe in their borders, but upon investigation of the locations of where those samples were, they seemed to be mostly from Eastern parts of Russia (some near the Bering Strait/Alaska)*
- This is a solid, geographically sensible spread — good coverage from Iberia to the Baltics to the Balkans, with Hungary (1,127) and the UK (937) as your two largest groups, and predictably thin representation in small places like Malta (2) and Luxembourg (2). Worth a mental note for later: any breakdown by individual country will be shaky for those tiny-n places, but it won't hurt the Europe-wide aggregate we're about to build.
Now let's actually build the time-binned frequency table — this is the core empirical result for Q1.
```r
# Assign each individual to a 1,000-year BP bin (Segurel et al. 2020 convention)
merged_europe$bp_bin <- floor(merged_europe$date_mean_bp / 1000) * 1000

# Aggregate: sample size and allele frequency per bin
freq_by_bin <- aggregate(A_count ~ bp_bin, data = merged_europe,
                          FUN = function(x) c(n = length(x),
                                               freq = sum(x) / (2 * length(x))))
freq_by_bin <- do.call(data.frame, freq_by_bin)
names(freq_by_bin) <- c("bp_bin", "n", "freq")
```
```r
# BP counts DOWN toward the present, so sort by bp_bin descending
# to read the table chronologically: oldest first, most recent last
freq_by_bin <- freq_by_bin[order(-freq_by_bin$bp_bin), ]
freq_by_bin
```
Before trusting whatever that prints, there's one thing worth checking specifically: your research question runs "from the Neolithic through today," but AADR's ancient-DNA files sometimes also include present-day reference genomes bundled in alongside the truly ancient samples. If any exist here, they'd likely be dated 0 (or very close to it) and would land in your most-recent bin — which matters, because a bin mixing "a person who died 800 years ago" with "a living person genotyped in 2020" isn't really one coherent time period, and it's also your best opportunity to sanity-check the whole pipeline against a well-known real-world number (modern Northern Europeans run roughly 70–90% derived allele frequency).
```r
summary(merged_europe$date_mean_bp)
sum(merged_europe$date_mean_bp == 0)
sum(merged_europe$date_mean_bp < 100)
```
Paste back the binned frequency table and these three checks. If there's a real cluster at or near 0, we'll likely want to break "modern" out as its own explicit category rather than let it blend into the last ancient bin — but let's see the actual numbers before deciding.
This is a genuinely exciting result — before we get to the one thing that needs checking, look at the actual pattern: frequency sits at essentially zero from 53,000 BP all the way down through 5,000 BP, then climbs — 1.5% at 4000 BP, 5.4% at 3000 BP, 20.8% at 2000 BP, 27.9% at 1000 BP. That's not just "roughly what the literature says" — it's the specific shape (near-total absence through the Neolithic, rise starting around the Bronze/Iron Age transition) that Mathieson, Burger, and Segurel all independently reported. Your own extraction from raw AADR data is reproducing a published finding. That's a strong result to build your report around.
Now, the thing worth stopping on: that last bin (labeled 0, meaning "0–999 BP") has 1,298 individuals, and 561 of them are dated at exactly 0. A precise date of exactly zero for over 500 people isn't a coincidence — that's very likely AADR's placeholder for present-day reference genomes (modern population panels like Human Origins or 1000 Genomes samples bundled into the compendium for comparison), not archaeological burials with an actual measured date. If that's right, your "most recent" bin is currently a mix of two very different things: genuine post-medieval/early-modern burials (dated somewhere in the 1–999 BP range, which legitimately belong on your timeline) and modern living-population reference samples (which don't have an archaeological date at all — 0 is just a stand-in for "now"). Blending those into one bin muddies exactly the endpoint of your story.
Let's check what these actually are before deciding anything:
```r
zero_bp <- merged_europe[merged_europe$date_mean_bp == 0, ]
length(unique(zero_bp$group_id))
head(unique(zero_bp$group_id), 35)
```
If those group names look like ethnonyms without any archaeological/cultural qualifier (e.g., "Spanish", "French", "Basque" rather than something like "France_Medieval" or a named burial site), that confirms they're modern reference samples, and we'd want to split them into their own explicit "present-day" category rather than let them sit inside a BP bin at all — which actually strengthens your project, since a distinct modern endpoint is a natural, clean way to show "here's where this landed today."

- Let's split it out properly and use it as a sanity check against known modern lactose-tolerance geography — Northern Europeans typically run roughly 70–90% derived-allele frequency, Southern Europeans meaningfully lower (often 40–60%), so if your modern subset reproduces that north-south split, that's strong independent confirmation the whole pipeline is sound.
```r
# Separate genuinely modern reference samples from ancient archaeological samples
modern <- merged_europe[merged_europe$date_mean_bp == 0, ]
ancient <- merged_europe[merged_europe$date_mean_bp > 0, ]

# Overall modern frequency, as an endpoint
sum(modern$A_count) / (2 * nrow(modern))

# Rough north/south check using a few countries with decent n
table(modern$political_entity[modern$political_entity %in%
      c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")])
sum(modern$A_count[modern$political_entity %in% c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")]) /
  (2 * sum(modern$political_entity %in% c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")))

table(modern$political_entity[modern$political_entity %in% c("Italy", "Spain", "Greece")])
sum(modern$A_count[modern$political_entity %in% c("Italy", "Spain", "Greece")]) /
  (2 * sum(modern$political_entity %in% c("Italy", "Spain", "Greece")))

# Rebuild the archaeological-only time bins, now that modern samples are removed
ancient$bp_bin <- floor(ancient$date_mean_bp / 1000) * 1000
freq_by_bin <- aggregate(A_count ~ bp_bin, data = ancient,
                          FUN = function(x) c(n = length(x), freq = sum(x) / (2 * length(x))))
freq_by_bin <- do.call(data.frame, freq_by_bin)
names(freq_by_bin) <- c("bp_bin", "n", "freq")
freq_by_bin <- freq_by_bin[order(-freq_by_bin$bp_bin), ]
freq_by_bin
```
*pasted my code output*

And it passed, directionally: North (65.8%) well above South (24.7%) — exactly the pattern real-world lactose tolerance geography predicts (Northern Europeans have far higher lactase persistence than Southern Europeans). That's good.
Two honest caveats worth keeping in mind, not treating this as a precise measurement:
Your "Finland, Iceland, Norway, Sweden, UK" group is really "Finland (105) and UK (111) doing all the work," with Iceland and Norway at n=1 each contributing essentially nothing statistically. So 65.8% is really telling you about Finnish + British frequency, not a true pan-Nordic average — and it's Finland that's slightly pulling this down from the ~80-90% you sometimes see quoted for Scandinavia specifically, since Finland genuinely runs a bit lower than Sweden/Denmark for this allele.
Similarly, "Italy, Spain, Greece" is really Italy (164) and Spain (109) with Greece at n=4, contributing almost nothing.
None of that undermines the check — the direction is what mattered here, not the exact percentage — but it's the kind of caveat worth a one-line mention in your methods if you ever cite these specific numbers, rather than presenting them as precise national averages.
Two things worth working through before this becomes a chart — one is a real statistical problem, the other is a subtler interpretation issue that's easy to miss.
The statistical problem: look at bins older than about 10,000 BP — most have n of 1, 2, or 3. Bin 39000 shows freq = 1.0 — that's not "100% of the population carried this allele 39,000 years ago," it's one single person, who either did or didn't carry it, with nothing to average against. A line chart connecting these dots would show wild, meaningless swings between 0% and 100% for a period where you simply don't have enough data to say anything at all — and it would visually overwhelm the real signal starting around 8,000–10,000 BP where sample sizes actually become reasonable (33+ individuals and climbing). This is exactly the kind of thing that would undermine an "informed non-specialist" reader's trust if they noticed it, and it's also just not statistically honest to plot.
The natural fix, and it actually fits your research question well: your question is about "the Neolithic through today," so the deep pre-Neolithic era only needs to serve as a baseline — "before farming, this allele was essentially absent" — not a finely-resolved timeline. I'd collapse everything older than 10,000 BP into a single "Pre-Neolithic (>10,000 BP)" category, and keep individual 1,000-year bins from 10,000 BP down to the present, where you actually have the sample sizes to say something meaningful.
The subtler issue: notice that your most recent archaeological bin (0–999 BP, n=737) shows 48.0% — actually higher than the 43.4% Europe-wide modern average from the last step. That's not necessarily a contradiction, but it's a real trap worth naming: each bin's "Europe" isn't the same set of countries. Bin composition shifts depending on which sites happened to get sequenced — if the 0–999 BP bin happens to be dominated by, say, a large medieval England or Poland project, its frequency reflects those places, not a uniform Europe-wide sample. A single aggregated line risks implying a smoother, more uniform trend than the underlying data actually supports.
This doesn't need a fix right now — it's more something to be upfront about in your methods ("bin composition varies by which regions have been sequenced most heavily for each period") and something we should keep in mind when picking a chart type (a chart that can hint at sample composition, like point size scaled to n, would be more honest than a plain line).

- Good — collapsing it is the right call, and it's worth stating that reasoning directly in your methods section too: "bins with fewer than ~10 individuals were consolidated into a single pre-Neolithic baseline category to avoid over-interpreting small-sample noise" reads as careful methodology, not a shortcut. 

```r
# Collapse everything older than 10,000 BP into one baseline category;
# keep individual 1,000-year bins from 10,000 BP to the present, where
# sample sizes are actually large enough to be meaningful
ancient$bp_bin <- ifelse(ancient$date_mean_bp >= 10000,
                          "Pre-Neolithic (>10,000 BP)",
                          as.character(floor(ancient$date_mean_bp / 1000) * 1000))

freq_by_bin <- aggregate(A_count ~ bp_bin, data = ancient,
                          FUN = function(x) c(n = length(x), freq = sum(x) / (2 * length(x))))
freq_by_bin <- do.call(data.frame, freq_by_bin)
names(freq_by_bin) <- c("bp_bin", "n", "freq")

# Order rows chronologically: Pre-Neolithic first, then bins descending
# toward the present (this is just for display - it's a factor now, not a number)
bin_order <- c("Pre-Neolithic (>10,000 BP)", as.character(seq(9000, 0, by = -1000)))
freq_by_bin$bp_bin <- factor(freq_by_bin$bp_bin, levels = bin_order)
freq_by_bin <- freq_by_bin[order(freq_by_bin$bp_bin), ]
freq_by_bin

# Building modern data point in the same structure as freq_by_bin
modern_freq <- sum(modern$A_count) / (2 * nrow(modern))
modern_row <- data.frame(bp_bin = "Modern (Present)", n = nrow(modern), freq = modern_freq, stringsAsFactors = FALSE)

# Flag archaeological vs. modern so they can be styled differently
freq_by_bin$bp_bin <- as.character(freq_by_bin$bp_bin) # Drop factor temporarily to combine
freq_by_bin$period_type <- "Archaeological"
modern_row$period_type <- "Modern"

# Actually combine them - this step was missing
freq_by_bin_full <- rbind(freq_by_bin, modern_row)

# Re-apply chronological order with Modern as the final point
bin_order_full <- c(bin_order, "Modern (Present)")
freq_by_bin_full$bp_bin <- factor(freq_by_bin_full$bp_bin, levels = bin_order_full)
freq_by_bin_full <- freq_by_bin_full[order(freq_by_bin_full$bp_bin), ]
```
One thing to check once this prints: the pre-Neolithic bin should now have a much larger n (everything from your old 53000 down through 11000 BP bins combined) and a freq close to zero — that's the expected, honest baseline. If it's still small or oddly high, something didn't collapse correctly. 
This table is in great shape now — the collapse worked cleanly (n=93 in the baseline vs. those tiny single-digit bins before), and the story is intact: near-zero from the deep past through ~5,000 BP, then a clear, accelerating climb to 48% by the most recent archaeological period.
One small thing worth noticing, not fixing: the Pre-Neolithic baseline isn't exactly zero — it's 4.3% (roughly 8 allele copies out of ~186, across 93 people). That's actually a nice, honest detail rather than a problem. Segurel et al. specifically note a handful of derived-allele carriers turning up even before the frequency really takes off — rare, sporadic occurrences rather than a clean absence. You can genuinely say your baseline matches that same "very low but not literally zero" pattern reported in the literature, which is a stronger, more credible claim than "we found nothing."



## Entry 2:
- **Tool:** Claude (Sonnet 5)
- **What it was used for:** It was used to help me write the code to refine a plot that I created myself. It specifically added `scale_size = ... trans = sqrt` to scale the dots better on the graph and `theme(axis.text.x = element_text(angle = 45, hjust = 1))` upon my request to turn the x-axis labels to be slanted because I thought it would improve the readability of the plot. 
- **Output:** 
*BUILT ON MY OWN BESIDES ABOVE RECOMMENDATIONS*: `scale_size = ... trans = sqrt` AND `theme(axis.text.x = element_text(angle = 45, hjust = 1))`
Here's the base chart with point size scaled to n (using a square-root size transform, since n ranges from 34 to 2,259 — a huge spread that would make small bins invisible and the biggest bin swallow the plot without it):
```r
library(ggplot2)

bin_labels <- c(
  "Pre-Neolithic (> 10,000 BP)" = "Pre-Neolithic\n(>10,000 BP)",
  "9000" = "9,000-10,000\nBP",
  "8000" = "8,000-9,000\nBP",
  "7000" = "7,000-8,000\nBP",
  "6000" = "6,000-7,000\nBP",
  "5000" = "5,000-6,000\nBP",
  "4000" = "4,000-5,000\nBP",
  "3000" = "3,000-4,000\nBP",
  "2000" = "2,000-3,000\nBP",
  "1000" = "1,000-2,000\nBP",
  "0" = "0-1,000\nBP",
  "Modern (Present)" = "Present\nDay"
)

ggplot(freq_by_bin_full, aes(x = bp_bin, y = freq, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(size = n, color = period_type, shape = period_type), alpha = 0.85) +
  scale_size(range = c(5, 15), trans = "sqrt", name = "Sample size (n)") +
  scale_color_manual(values = c("Archaeological" = "steelblue", "Modern" = "firebrick"), name = "Sample type") +
  scale_shape_manual(values = c("Archaeological" = 16, "Modern" = 17), name = "Sample Type") +
  scale_x_discrete(labels = bin_labels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(
    title = "Rise of the Lactase Persistence Allele (rs4988235) in Europe",
    subtitle = "Frequency of the derived, milk-digesting allele across archaeological time periods",
    x = "Time Period",
    y = "Allele Frequency",
    caption = "Data: Allen Ancient DNA Resource (AADR) v66.p1.\nPoint size reflects sample size (n) per bin."
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(title = "Sample Type"),
         shape = guide_legend(title = "Sample Type"))
ggsave("figures/q1_lactase_timeseries.png", width = 10, height = 6, dpi = 300)
#############################
##### SENSITIVITY CHECK #####
#############################
# The 0-1000 BP bin is dominated by Belgium (225 of 737 individuals, ~30%).
# Checking whether this single country is driving the bin's frequency:
last_bin_countries <- ancient$political_entity[ancient$bp_bin == "0"]
sort(table(last_bin_countries), decreasing = TRUE)
# Result: 44.3% without Belgium vs. 48% with it - close to the modern estimate (43.4%),
#suggesting the apparent recent dip in the main chart is a sampling-composition artifact, not a real decline. 
  ```