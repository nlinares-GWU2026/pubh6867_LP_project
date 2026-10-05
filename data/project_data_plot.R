#########################
##### DATA CLEANING #####
#########################
raw <- read.table("data/lactase_extract_raw.txt", header = TRUE, stringsAsFactors = FALSE)

# PLINK's rs4988235_G column counts the ancestral G allele (0 "GG"/1 "AG"/2 "AA") (non LP)
# Want the derived, LP allele (A) instead so flip: 
# Calculates the derived, lactase-persistence A allele for each individual by subtracting the count of the ancestral (G) allele 
# from 2, converting 0 to 2, 1 to 1, and 2 to 0
raw$A_count <- 2 - raw$rs4988235_G 

# Drop individuals with no genotype data (NA values in A_count) and creates a simple
# dataframe only containing IID (individual ID and A_count)
raw_clean <- raw[!is.na(raw$A_count), c("IID", "A_count")]

# Checking
#nrow(raw_clean) # 15523
#head(raw_clean) # IID = Loschbour. AG, Kou01.SG, etc. | A_counts = 0 

# Read the annotation file
# Reads tab delimited file, treats the first row of the file as column titles, disables quote interpretation so that quote marks in the text 
# fields are read as literal characters rather than string delimiters, and ensures text columns are loaded as char vectors vs converted into categorical
anno <- read.delim("data/v66.p1_1240K.aadr.PUB.anno", header = TRUE, quote = "", stringsAsFactors =  FALSE)
# Checking
#nrow(anno) # 23089
#ncol(anno) # 49
#names(anno) # Header row embeds full documentation text for each column -> rename to avoid retyping

# Renaming
# Column 1 = Genetic ID (explains the .SG shotgun genome/.AG Agilent capture/.DG diploid genome)
# Column 11 = date mean in BP (years before 1950)
# Column 12 = date standard deviation in BP
# Column 13 = human readable calibrated date range
# Column 15 = Group.ID the archaeological group/culture label 
# Column 16-19 = locality, country, latitude, longitude labels for regional grouping 

names(anno)[1] <- "id"
names(anno)[11] <- "date_mean_bp"
names(anno)[12] <- "date_sd_bp"
names(anno)[13] <- "full_date"
names(anno)[15] <- "group_id"
names(anno)[16] <- "locality"
names(anno)[17] <- "political_entity"
names(anno)[18] <- "latitude"
names(anno)[19] <- "longitude"

# Subset of 9 specific metadata columns from the full anno df 
anno_sub <- anno[, c("id", "group_id", "locality", "political_entity",
                     "latitude", "longitude", "date_mean_bp", "date_sd_bp", "full_date")]
# Handling the latitude and longitude unknown values - 784 have a DATE but NO LOCATION - use in Q1 but NOT in Q2
anno_sub$latitude <- as.numeric(anno_sub$latitude)
anno_sub$longitude <- as.numeric(anno_sub$longitude)
#sum(is.na(anno_sub$latitude)) # 784 missing
#sum(is.na(anno_sub$longitude)) # 784 misisng

# Checking
#head(anno_sub)
#str(anno_sub)

# Rename the column "IID" in raw_clean to "id" so it matches the column name in anno_sub
names(raw_clean)[names(raw_clean) == "IID"] <- "id"
merged <- merge(raw_clean, anno_sub, by = "id", all.x = TRUE)

#Checking
#nrow(merged) # 15523 rows matches raw_clean
#sum(is.na(merged$group_id)) # 0 unmatched group ID's - all have metadata
#head(merged) # Now have one table "merged" with genotype, date, location, and culture for every sample

# Deciding on whether to include Russia/Turkey and if so, how to handle them _ See DISCUSSION for explanation
#russia_turkey <- merged[merged$political_entity %in% c("Russia", "Turkey"), c("political_entity", "locality", "longitude", "latitude")]
#table(merged$political_entity[merged$political_entity %in% c("Russia", "Turkey")])
#russia_turkey[order(russia_turkey$political_entity, russia_turkey$longitude), ]


# Binning
# Checking actual names of all countries included in AADR to establish "European" boundary
#sort(unique(merged$political_entity))

europe_countries <- c("Albania", "Austria", "Belgium", "Bulgaria", "Channel Islands", "Croatia",
                      "Crimea", "Czechia", "Denmark", "Estonia", "Faroe Islands", "Finland", "France", 
                      "Germany", "Gibraltar", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
                      "Latvia", "Lithuania", "Luxembourg", "Malta", "Moldova", "Montenegro", "Netherlands",
                      "North Macedonia", "Norway", "Poland", "Portugal", "Romania", "Serbia", "Slovakia",
                      "Slovenia", "Spain", "Sweden", "Switzerland", "Ukraine", "United Kingdom")
merged_europe <- merged[merged$political_entity %in% europe_countries, ]
nrow(merged_europe) #7711
table(merged_europe$political_entity) # Hungary (1127) and UK (937) biggest; Malta (2) and Luxemborg (2) smallest


# Assign each individual a 1,000 year bin (Segurel et al., 2020 convention) - Aged 4,350 BP -> 4.35 -> 4 -> 4,000
#merged_europe$bp_bin <- floor(merged_europe$date_mean_bp / 1000) * 1000

# Groups individuals by their 1000 year time bin (bp_bin) and calculates both sample size and allele frequency
#freq_by_bin <- aggregate(A_count ~ bp_bin, data = merged_europe,
#                       FUN = function(x) c(n = length(x), 
#                                           freq = sum(x) / (2 * length(x))))
#freq_by_bin <- do.call(data.frame, freq_by_bin)
#names(freq_by_bin) <- c("bp_bin", "n", "freq")
#freq_by_bin <- freq_by_bin[order(-freq_by_bin$bp_bin), ] # BP counts DOWN toward the present, sort descending for oldest first and most recent last 
#freq_by_bin

#summary(merged_europe$date_mean_bp)
#sum(merged_europe$date_mean_bp == 0) # 561
#sum(merged_europe$date_mean_bp < 100) #561
# Result showed last bin labeled "0" with 1298 individuals (exactly 561 at "present/0". This "0" bin means 0-999 BP or MODERN populations.
# This could be a mix of genuine post-medival/early-modern burials (1-999 BP) - which belong on the timeline AND modern reference samples (0 specifically)

# --> Check what these modern populations are:
zero_bp <- merged_europe[merged_europe$date_mean_bp == 0, ] # rows = 0, keep all columns
length(unique(zero_bp$group_id)) # 32 present day populations
head(unique(zero_bp$group_id), 35)

# Rebuilding archaeological-only time bin with modern samples are removed
modern <- merged_europe[merged_europe$date_mean_bp == 0, ]
ancient <- merged_europe[merged_europe$date_mean_bp > 0, ]

sum(modern$A_count) / (2*nrow(modern)) # Overall frequency, as an endpoint

# Checking if the modern population LP percentages match real-world percentages. 
#table(modern$political_entity[modern$political_entity %in%
#                                c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")])
#sum(modern$A_count[modern$political_entity %in% c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")]) /
#  (2 * sum(modern$political_entity %in% c("United Kingdom", "Finland", "Iceland", "Norway", "Sweden")))

#table(modern$political_entity[modern$political_entity %in% c("Italy", "Spain", "Greece")])
#sum(modern$A_count[modern$political_entity %in% c("Italy", "Spain", "Greece")]) /
#  (2 * sum(modern$political_entity %in% c("Italy", "Spain", "Greece")))

#ancient$bp_bin <- floor(ancient$date_mean_bp / 1000) * 1000
#freq_by_bin <- aggregate(A_count ~ bp_bin, data = ancient,
#                        FUN = function(x) c(n = length(x), freq = sum(x) / (2 * length(x)))) 
#freq_by_bin <- do.call(data.frame, freq_by_bin)
#names(freq_by_bin) <- c("bp_bin", "n", "freq")
#freq_by_bin <- freq_by_bin[order(-freq_by_bin$bp_bin), ]
#freq_by_bin

# Due to extremely small prevalence of LP in pre 10,000 BP - decided to collapse bins >10,000 BP into one category 
# In order to not disrupt signal for intepretation
# Keep individual 1,000-year bins from 10,000 BP to present where sample sizes are actually large enough to be meaningful
ancient$bp_bin <- ifelse(ancient$date_mean_bp >= 10000,
                         "Pre-Neolithic (> 10,000 BP)", as.character(floor(ancient$date_mean_bp / 1000) * 1000))
# Groups individuals by their 1000 year time bin (bp_bin) and calculates both sample size and allele frequency
freq_by_bin <- aggregate(A_count ~ bp_bin, data = ancient, # Groups the target variable (A_count) by the grouping factor (bp_bin)
                         FUN = function(x) c(n = length(x),  # Counts total number of individuals in that time bin
                                             freq = sum(x) / (2 * length(x)))) # totals the count of A alleles across all individuals in that bin and calculates the total number of alleles in that gene (out of 2)
freq_by_bin <- do.call(data.frame, freq_by_bin) # Flatten into data frame
names(freq_by_bin) <- c("bp_bin", "n", "freq")
bin_order <- c("Pre-Neolithic (> 10,000 BP)", as.character(seq(9000, 0, by = -1000))) # Order rows chronologically: Pre-Neo first -> descending towards present
freq_by_bin$bp_bin <- factor(freq_by_bin$bp_bin, levels = bin_order)
freq_by_bin <- freq_by_bin[order(freq_by_bin$bp_bin), ]
freq_by_bin

# Building modern data point in the same structure as freq_by_bin
modern_freq <- sum(modern$A_count) / (2 * nrow(modern))
modern_row <- data.frame(bp_bin = "Modern (Present)", n = nrow(modern), freq = modern_freq, stringsAsFactors = FALSE)

# Flag archaeological vs. modern so they can be styled differently later on
freq_by_bin$bp_bin <- as.character(freq_by_bin$bp_bin) # Drop factor temporarily to combine
freq_by_bin$period_type <- "Archaeological"
modern_row$period_type <- "Modern"
freq_by_bin_full <- rbind(freq_by_bin, modern_row) # Combine them

# Re-apply chronological order with Modern as the final point
bin_order_full <- c(bin_order, "Modern (Present)")
freq_by_bin_full$bp_bin <- factor(freq_by_bin_full$bp_bin, levels = bin_order_full)
freq_by_bin_full <- freq_by_bin_full[order(freq_by_bin_full$bp_bin), ]

##################################################################################### 
##### Q1 PLOT: HOW DID FREQ OF LP CHANGE ACROSS EUROPE FROM NEOLITHIC - MODERN? #####
#####################################################################################
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



#######################################################################################
##### Q2 PLOT: ARE STEPPE ASSOCIATED GROUPS ASSOCIATED WITH TIMING AND GEOGRAPHY? #####
#######################################################################################

# Finding which archaeological groups are actually in my data before labeling as "steppe-associated"
# Haak et al., 2015 and Allentoft et al., 2015 name specific archaeological cultures tied to Bronze Age steppe migration:
# Yamnaya (steppe population itself), Corded Ware, Bell Beaker (which both emerged in Europe after)
# Also a few cultures preceding Yamanaya west spread: Catacomb, Srubnaya, and Sintashta

unique_groups <- sort(unique(ancient$group_id))
length(unique_groups)

# Known steppe-associated culture names from aDNA lit
# (Haak et al., 2015 and Allentoft et al., 2015)
steppe_keywords <- c("Yamnaya", "Corded_Ware", "CordedWare", "Bell_Beaker", "BellBeaker",
                     "Catacomb", "Srubnaya", "Sintashta", "Afanasievo", "Potapovka",
                     "Poltavka", "Steppe")
steppe_candidates <- unique_groups[grepl(paste(steppe_keywords, collapse = "|"),
                                         unique_groups, ignore.case = TRUE)]
length(steppe_candidates) # Length was 96
steppe_candidates

######################################################################
##### Bell Beaker Mixed Ancestry & Quantifying Count Individuals #####
######################################################################

cand <- data.frame(group_id = steppe_candidates, stringsAsFactors = FALSE)

# Outliers caught first so "Yamnaya-o" lands in outlier, not core
cand$category <- ifelse(grepl("-o", cand$group_id), "Outlier (exclude)",
                        ifelse(grepl("Scythian", cand$group_id), "Scythian (exclude)",
                               ifelse(grepl("BellBeaker", cand$group_id), "Bell Beaker",
                                      ifelse(grepl("Yamnaya|Catacomb|CordedWare", cand$group_id),
                                             "Core steppe", "Other"))))
table(cand$category) # Groups per category
cand[cand$category == "Other", ] # Confirmed empty

# Individuals per category, and which time bins they go into
ancient_cand <- merge(ancient, cand, by = "group_id")
table(ancient_cand$category)
table(ancient_cand$category, ancient_cand$bp_bin)

# Any groups AADR flagged for exclusion
sum(grepl("^Ignore_", ancient$group_id))

