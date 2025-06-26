# MutSeqRData provides example data for MutSeqR package.
# Data is MutaMouse bonemarrow DNA sequenced using TwinStrand's Duplex
# Sequencing. Data was taken from LeBlanc et al., 2022 study.
# Adult MutaMouse males 9-14 weeks of age were exposed to three doses of
# Benzo[a]pyrene alongside vehicle controls for 28 days by oral gavage. 28
# days after the end of the exposure period, bone marrow was harvested from
# the femurs of euthanized animals. DNA was extracted using Qiagen's DNeasy
# Blood and Tissue kit. Libraries were built using TwinStrand Biosciences' Mouse
# Mutagenesis kit. The mouse mutagenesis panel consists of 20 2.4kb genomic
# targets spread across the mouse autosomes (48kb total). Libraries were
# sequenced on the NovaSeq 6000 at >10,000X depth. Demultiplexed FASTQ files
# were processed through the TwinStrand Biosciences Duplex Seq Mutagenesis App
# (v.3.20.1) hosted on DNAnexus platform. Pre-processing included extracting
# unique molecule identifier (UMI) sequences, correcting UMI sequences, raw
# read alignment, grouping reads by their UMI and strand defining element,
# quality trimming, error-correction of read groups by duplex consensus calling,
# consensus post-processing, re-alignment, and variant calling (described in
# Valentine et al., 2020). Pre-processing produced tabular and vcf files of
# mutation data for each mouse sample (24). The tabular and vcf files of mouse
# sample dna00996.1 were downloaded from DNAnexus and included in MutSeqRData
# to be used as example data for MutSeqR::import_mut_data() and
# MutSeqR::import_vcf_data respectively. The tabular mutation data for
# dna00996.1 was first loaded into R with read.delim(), and saved as an RDS
# file: Example import_mut_data. A copy of the tabular mutation data was
# included with altered column names to further showcase function utility:
# Example import_mut_data using Custom column names. The vcf file was included
# unaltered: Example import_vcf_data.
#
# The tabular mutation data files for the 24 mouse samples were downloaded from
# DNAnexus and imported into R using MutSeqR::import_mut_data() as shown in the
# following code. The resulting data frame was saved as an RDS file: Example
# mutation data.
library(MutSeqR)
sample_data <- data.frame(
  sample = c("dna00973.1", "dna00974.1", "dna00975.1", "dna00976.1",
             "dna00977.1", "dna00978.1", "dna00979.1", "dna00980.1",
             "dna00981.1", "dna00982.1", "dna00983.1", "dna00984.1",
             "dna00985.1", "dna00986.1", "dna00987.1", "dna00988.1",
             "dna00989.1", "dna00990.1", "dna00991.1", "dna00992.1",
             "dna00993.1", "dna00994.1", "dna00995.1", "dna00996.1"),
  dose_group = c(rep("Control", 6), rep("Low", 6), rep("Medium", 6),
                 rep("High", 6)),
  dose = c(rep(0, 6), rep(12.5, 6), rep(25, 6), rep(50, 6))
)
mutation_data <- import_mut_data(
  mut_file = "path/to/mut/files",
  is_0_based_mut = TRUE,
  sample_data = sample_data,
  regions = "TSpanel_mouse",
  genome = "mm10",
  species = "mouse"
)
saveRDS(mutation_data, file = "example_mutation_data.rds")

# mutation_data was run through MutSeqR::filter_mut()
# as demosntrated below. The resulting data frame was saved as an RDS
# file: Example mutation data filtered.
filtered_mutation_data <- filter_mut(
  mutation_data = mutation_data,
  vaf_cutoff = 0.01,
  snv_in_germ_mnv = TRUE,
  rm_abnormal_vaf = FALSE,
  custom_filter_col = "filter",
  custom_filter_val = "EndRepairFillInArtifact",
  custom_filter_rm = FALSE,
  regions = "TSpanel_mouse",
  regions_filter = "keep_within",
  allow_half_overlap = FALSE,
  rm_filtered_mut_from_depth = TRUE
)
saveRDS(filtered_mutation_data, file = "example_mutation_data_filtered.rds")

# The Precalculated Depth files were generated using MutSeqR::calculate_mf()
# as described below. They were saved as RDS files.
library(dplyr)
depth_6 <- calculate_mf(
  mutation_data = filtered_mutation_data,
  cols_to_group = "sample",
  subtype_resolution = "base_6",
  calculate_depth = TRUE,
  correct_depth = TRUE,
  correct_depth_by_indel_priority = TRUE
) %>%
  filter(normalized_ref != "N") %>%
  select(sample, normalized_ref, group_depth, subtype_depth)
saveRDS(depth_6, file = "precalc_depth_base_6_example.rds")

depth_12 <- calculate_mf(
  mutation_data = filtered_mutation_data,
  cols_to_group = "sample",
  subtype_resolution = "base_12",
  calculate_depth = TRUE,
  correct_depth = TRUE,
  correct_depth_by_indel_priority = TRUE
) %>%
  filter(short_ref != "N") %>%
  select(sample, short_ref, group_depth, subtype_depth)
saveRDS(depth_12, file = "precalc_depth_base_12_example.rds")

depth_96 <- calculate_mf(
  mutation_data = filtered_mutation_data,
  cols_to_group = "sample",
  subtype_resolution = "base_96",
  calculate_depth = TRUE,
  correct_depth = TRUE,
  correct_depth_by_indel_priority = TRUE
) %>%
  filter(normalized_context_with_mutation != "N") %>%
  select(sample, normalized_context_with_mutation, group_depth, subtype_depth)
saveRDS(depth_96, file = "precalc_depth_base_96_example.rds")

depth_192 <- calculate_mf(
  mutation_data = filtered_mutation_data,
  cols_to_group = "sample",
  subtype_resolution = "base_192",
  calculate_depth = TRUE,
  correct_depth = TRUE,
  correct_depth_by_indel_priority = TRUE
) %>%
  filter(context_with_mutation != "N") %>%
  select(sample, context_with_mutation, group_depth, subtype_depth)
saveRDS(depth_192, file = "precalc_depth_base_192_example.rds")


