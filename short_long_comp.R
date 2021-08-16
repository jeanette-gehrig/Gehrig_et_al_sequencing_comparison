# Install packages if not already installed
install.packages("vegan")
install.packages("RAM")
install.packages("janitor")
install.packages("ade4")
install.packages("reshape2")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("scales")
install.packages("BiocManager")
BiocManager::install("Maaslin2")
install.packages("pheatmap")

# Load package libraries
library(vegan)
library(RAM)
library(janitor)
library(ade4)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(BiocManager)
library(Maaslin2)
library(pheatmap)

# Change working directory to location where the Gehrig_et_al_sequencing_comparison repository was cloned locally
setwd("~/Desktop/Gehrig_et_al_sequencing_comparison")

# Read in metadata for amplicon sequencing data
amplicon_metadata<-read.table("input_files/amplicon_sequencing_metadata_for_R.txt", sep = "\t", header = T)

# Make metadata columns factors for plotting
amplicon_metadata$Participant<-as.factor(amplicon_metadata$Participant)
amplicon_metadata$Patient<-as.factor(amplicon_metadata$Patient)
amplicon_metadata$Timepoint<-as.factor(amplicon_metadata$Timepoint)
amplicon_metadata$DNA_extraction_method<-as.factor(amplicon_metadata$DNA_extraction_method)

# Subset out only StrainID samples (remove V3-V4 samples) from metadata table
strainid_metadata<-subset(amplicon_metadata,amplicon_metadata$DNA_extraction_method!="V3-V4")

# Read in StrainID strain-level taxonomy
st_strainid<-read.table("input_files/StrainID_strain_level_all_samples_taxonomy_table.txt", sep = "\t", header = T, check.names = F)

# Change last column name, since pcoa.plot (RAM) requires the last column to be "taxonomy"
colnames(st_strainid)[ncol(st_strainid)]<-"taxonomy"

# Make the rankID column the rownames, and remove the rankID column
rownames(st_strainid)<-st_strainid$rankID
st_strainid$rankID<-NULL

# Check to make sure order of samples in StrainID relative abundance table matches order of samples in metadata table
reorder_idx<-match(colnames(st_strainid)[1:ncol(st_strainid)-1],strainid_metadata$Sample)
strainid_metadata_ordered<-strainid_metadata[reorder_idx,]

# PCoA plot of StrainID strain data by participant
pcoa.plot(st_strainid, is.OTU=TRUE, meta=strainid_metadata_ordered, factors=c("Timepoint","Participant"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, ellipse = 2, main = NULL, file = NULL, ext = NULL, ggplot2 = FALSE, bw = FALSE)

# PCoA plot of StrainID strain data by DNA extraction method
pcoa.plot(st_strainid, is.OTU=TRUE, meta=strainid_metadata_ordered, factors=c("Participant","DNA_extraction_method"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, ellipse = 2, main = NULL, file = NULL, ext = NULL, ggplot2 = FALSE, bw = FALSE)

## Bray Curtis distance matrix for StrainID strain data

# Remove taxonomy column
st_strainid_notax<-st_strainid[,1:(ncol(st_strainid)-1)]

# Transpose
st_strainid_notax_t<-t(st_strainid_notax)

# Calculate Bray-Curtis distance matrix
braycurtis_st_strainid<-vegdist(st_strainid_notax_t, method = "bray")

# Calculate the PERMANOVA significance of grouping by Participant vs. DNA extraction method with Bray-Curtis distance matrix
adonis2(braycurtis_st_strainid~strainid_metadata_ordered$Participant,strainid_metadata_ordered, method = "bray")
adonis2(braycurtis_st_strainid~strainid_metadata_ordered$DNA_extraction_method,strainid_metadata_ordered, method = "bray")

# Get StrainID strain-level table with only the 12 samples selected
# Subset metadata table to get only StrainID replicates
strainid_metadata_ordered_sl<-subset(strainid_metadata_ordered, strainid_metadata_ordered$DNA_extraction_method=="Shoreline_Complete")
# For each StrainID replicate, get the sample with the higher number of reads
most_reads<-aggregate(strainid_metadata_ordered_sl$Reads,by=list(strainid_metadata_ordered_sl$SampleName),max)
most_reads_index<-match(most_reads$x,strainid_metadata_ordered_sl$Reads)
strainid_metadata_sl_most_reads<-strainid_metadata_ordered_sl[most_reads_index,]

# Get StrainID reads count table with only the StrainID replicate with the higher number of reads
strainid_most_reads_index<-match(strainid_metadata_sl_most_reads$Sample,colnames(st_strainid))
st_strainid_sl_replicate<-cbind(st_strainid[,strainid_most_reads_index],st_strainid$taxonomy)
colnames(st_strainid_sl_replicate)[ncol(st_strainid_sl_replicate)]<-"taxonomy"

# PCoA plot of StrainID strain data by participant, one replicate per sample
pcoa.plot(st_strainid_sl_replicate, is.OTU=TRUE, meta=strainid_metadata_sl_most_reads, factors=c("Patient","Participant"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, main = NULL, file = NULL, ext = NULL, ggplot2 = FALSE, bw = FALSE)

### Read in V3-V4 16S data, only the 12 samples, not rarefied, 10 read min across all samples from all cohorts, 2 sample minimum

v3_v4<-read.table("input_files/v3-v4_asv_table_non-rarefied.txt",sep = "\t", header = T, check.names = FALSE)

# Divide each sample by total number of reads in that sample, since sequencing depth varies widely across samples
v3_v4_percent<-adorn_percentages(v3_v4, denominator = "col")
# Verify that column sums are equal to 1
colSums(v3_v4_percent[,2:13])

# Make first column (ASV names) the rownames
rownames(v3_v4_percent)<-v3_v4_percent$ASV
v3_v4_percent$ASV<-NULL

# Get metadata for only v3-v4 data
v3_v4_metadata<-subset(amplicon_metadata, amplicon_metadata$DNA_extraction_method=="V3-V4")

# Check to make sure order of samples in v3-v4 relative abundance table matches order of samples in metadata table
reorder_idx<-match(colnames(v3_v4_percent)[1:ncol(v3_v4_percent)-1],v3_v4_metadata$Sample)
v3_v4_metadata_ordered<-v3_v4_metadata[reorder_idx,]

# PCoA plot - V3-V4
pcoa.plot(v3_v4_percent, is.OTU=TRUE, meta=v3_v4_metadata_ordered, factors=c("Patient","Participant"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, main = NULL, file = NULL, ext = NULL, height = 8, width = 10, ggplot2 = FALSE, bw = FALSE)

### Read in short-read metagenomics bracken data, minimum fractional abundance 0.0004
bracken<-read.table("input_files/short-read_metagenomics_bracken_taxonomy_min_abundance_0.0004.txt", sep = "\t", header = T, check.names = FALSE)

# Make first column (taxa) the rownames
rownames(bracken)<-bracken$Taxa
bracken$Taxa<-NULL

# Check to make sure order of samples in bracken relative abundance table matches order of samples in metadata table 
# Use the v3-v4 metadata table; the metadata is the same for short-read metagenomics
reorder_idx<-match(colnames(bracken)[1:ncol(bracken)-1],v3_v4_metadata$Sample)
v3_v4_metadata_ordered<-v3_v4_metadata[reorder_idx,]

# PCoA plot - bracken
pcoa.plot(bracken, is.OTU=TRUE, meta=v3_v4_metadata_ordered, factors=c("Patient","Participant"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, main = NULL, file = NULL, ext = NULL, height = 8, width = 10, ggplot2 = FALSE, bw = FALSE)

### Read in long-read metagenomics taxonomy data (GTDB normalized counts)

gtdb<-read.table("input_files/Protein-Taxa-GTDB-NormalizedCounts.txt", sep = "\t", header = T, check.names = F)

# Make first column (taxonomy) the rownames
rownames(gtdb)<-gtdb$Taxonomy
gtdb$Taxonomy<-NULL

# Check to make sure order of samples in bracken relative abundance table matches order of samples in metadata table 
# Use the v3-v4 metadata table; the metadata is the same for long-read metagenomics
reorder_idx<-match(colnames(gtdb)[1:ncol(gtdb)-1],v3_v4_metadata$Sample)
v3_v4_metadata_ordered<-v3_v4_metadata[reorder_idx,]

pcoa.plot(gtdb, is.OTU=TRUE, meta=v3_v4_metadata_ordered, factors=c("Patient","Participant"), rank=NULL, dist.method = "bray",  sample.labels = FALSE, top = 0, main = NULL, file = NULL, ext = NULL, height = 8, width = 10, ggplot2 = FALSE, bw = FALSE)

### Mantel tests to see if Bray Curtis distances between samples are correlated across methods ###

# Remove taxonomy columns (and the sample that is missing in the long-read metagenomics data: 3218_V6_D29) and transpose so samples are the rows and taxa are the columns
gtdb_notax_t<-t(gtdb[,1:ncol(gtdb)-1])
st_strainid_sl_replicate_notax_t<-t(st_strainid_sl_replicate[,1:(ncol(st_strainid_sl_replicate)-2)])
bracken_notax_t<-t(bracken[,1:(ncol(bracken)-2)])
v3_v4_percent_notax_t<-t(v3_v4_percent[,1:(ncol(v3_v4_percent)-2)])

# For StrainID matrix, make the sample names the same as other methods (remove replicate indication)
rownames(st_strainid_sl_replicate_notax_t)<-gsub('.{3}$','',rownames(st_strainid_sl_replicate_notax_t))

# Make distance matrices for all method's taxonomic profiles
gtdb_bray<-vegdist(gtdb_notax_t, method = "bray")
st_strainid_bray<-vegdist(st_strainid_sl_replicate_notax_t, method = "bray")
bracken_bray<-vegdist(bracken_notax_t, method = "bray")
v3_v4_bray<-vegdist(v3_v4_percent_notax_t, method = "bray")

# Mantel tests: measures correlation between two distance matrices (test for spatial autocorrelation)
mantel.rtest(st_strainid_bray, v3_v4_bray, nrepet = 9999)
mantel.rtest(st_strainid_bray, bracken_bray, nrepet = 9999)
mantel.rtest(v3_v4_bray, bracken_bray, nrepet = 9999)
mantel.rtest(v3_v4_bray, gtdb_bray, nrepet = 9999)
mantel.rtest(bracken_bray, gtdb_bray, nrepet = 9999)
mantel.rtest(st_strainid_bray, gtdb_bray, nrepet = 9999)

### Make stacked bar plots of different levels of taxonomy across all methods - Use only 11 samples shared across all methods ###

# Read in StrainID phylum table
strainid_phylum<-read.table("input_files/StrainID_phylum.txt", sep = "\t", header = T)
grep("3218_V6",colnames(strainid_phylum))

# Select phylum column and all sample columns except for sample 3218_V6 (missing in long read metagenomics) 
strainid_phylum_df<-strainid_phylum[,c(4,6:49)]

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
strainid_phylum_df_norm<-adorn_percentages(strainid_phylum_df, denominator = "col")
# Check to make sure column sums are 1
colSums(strainid_phylum_df_norm[,2:ncol(strainid_phylum_df_norm)])

# Add up total phyla across normalized samples
strainid_all_sample_phyla_totals<-rowSums(strainid_phylum_df_norm[,2:ncol(strainid_phylum_df_norm)])
strainid_phyla_totals_table<-as.data.frame(cbind(strainid_phylum_df_norm$phylum,strainid_all_sample_phyla_totals))

# Aggregate by phylum
strainid_phyla_sum<-aggregate(as.numeric(as.character(strainid_phyla_totals_table$strainid_all_sample_phyla_totals)),by = list(strainid_phyla_totals_table$V1), FUN=sum)
colnames(strainid_phyla_sum)<-c("phylum","StrainID")

## Get each sample's phyla sum separately for statistical tests
strainid_phylum_df_norm_collapsed<-aggregate.data.frame(strainid_phylum_df_norm[,2:ncol(strainid_phylum_df_norm)], by = list(strainid_phylum_df_norm$phylum), FUN = "sum")

## Read in V3-V4 percent abundance table with all levels of taxonomy
v3_v4_all_tax_levels<-read.table("input_files/V3-V4_percent_abundance_non-rarefied_gg_taxonomy_all_levels.txt", sep = "\t", header = T, check.names = F)

# Aggregate by phylum
v3_v4_phyla_sum<-aggregate(as.numeric(as.character(v3_v4_all_tax_levels$eleven_samp_sum)),by = list(v3_v4_all_tax_levels$phylum), FUN=sum)
colnames(v3_v4_phyla_sum)<-c("phylum","V3-V4")

# Get each sample's phyla sum and family sum separately for statistical tests
v3_v4_phylum_df_collapsed<-aggregate.data.frame(v3_v4_all_tax_levels[,2:12], by = list(v3_v4_all_tax_levels$phylum), FUN = sum)
v3_v4_family_df_collapsed<-aggregate.data.frame(v3_v4_all_tax_levels[,2:12], by = list(v3_v4_all_tax_levels$family), FUN = sum)

## Read in bracken taxonomy data with the plant and human fractional abundances removed, with all taxonomy levels
bracken_all_tax_levels<-read.table("input_files/short-read_metagenomics_bracken_taxonomy_min_abundance_0.0004_all_levels_taxonomy_contaminants_removed.txt", sep = "\t", header = T, check.names = F)

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
bracken_all_tax_levels_df<-bracken_all_tax_levels[,1:12]
bracken_all_tax_levels_df_norm<-adorn_percentages(bracken_all_tax_levels_df, denominator = "col")
# Check to make sure column sums are equal to 1
colSums(bracken_all_tax_levels_df_norm[,2:ncol(bracken_all_tax_levels_df_norm)])

# Add up total phyla across normalized samples
bracken_all_sample_phyla_totals<-rowSums(bracken_all_tax_levels_df_norm[,2:ncol(bracken_all_tax_levels_df_norm)])
bracken_phyla_totals_table<-as.data.frame(cbind(bracken_all_tax_levels$phylum,bracken_all_sample_phyla_totals))

# Aggregate by phylum
bracken_phyla_sum<-aggregate(as.numeric(as.character(bracken_phyla_totals_table$bracken_all_sample_phyla_totals)),by = list(bracken_phyla_totals_table$V1), FUN=sum)
colnames(bracken_phyla_sum)<-c("phylum","Short reads")

# Get each sample's phyla sum and family sum separately for statistical tests
bracken_phyla_totals_table_samples<-as.data.frame(cbind(bracken_all_tax_levels$phylum,bracken_all_tax_levels_df_norm))
colnames(bracken_phyla_totals_table_samples)[1]<-"phylum"
bracken_phylum_df_collapsed<-aggregate.data.frame(bracken_phyla_totals_table_samples[,3:13], by=list(bracken_phyla_totals_table_samples$phylum), FUN=sum)

bracken_family_totals_table_samples<-as.data.frame(cbind(bracken_all_tax_levels$family,bracken_all_tax_levels_df_norm))
colnames(bracken_family_totals_table_samples)[1]<-"family"
bracken_family_df_collapsed<-aggregate.data.frame(bracken_family_totals_table_samples[,3:13], by=list(bracken_family_totals_table_samples$family), FUN=sum)

## Read in long-read GTDB taxonomy data
gtdb_all_tax_levels<-read.table("input_files/Protein-Taxa-GTDB-NormalizedCounts_all_levels_taxonomy.txt", sep = "\t", header = T, check.names = F)

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
gtdb_all_tax_levels_df<-gtdb_all_tax_levels[,1:12]
gtdb_all_tax_levels_df_norm<-adorn_percentages(gtdb_all_tax_levels_df, denominator = "col")
# Check to make sure column sums are equal to 1
colSums(gtdb_all_tax_levels_df_norm[,2:ncol(gtdb_all_tax_levels_df_norm)])

# Add up total phyla across normalized samples
gtdb_all_sample_totals<-rowSums(gtdb_all_tax_levels_df_norm[,2:ncol(gtdb_all_tax_levels_df_norm)])
gtdb_phyla_totals_table<-as.data.frame(cbind(gtdb_all_tax_levels$phylum,gtdb_all_sample_totals))

# Aggregate by phylum
gtdb_phyla_sum<-aggregate(as.numeric(as.character(gtdb_phyla_totals_table$gtdb_all_sample_totals)),by = list(gtdb_phyla_totals_table$V1), FUN=sum)
colnames(gtdb_phyla_sum)<-c("phylum","Long reads")

# Get each sample's phyla sum and family sum separately for statistical tests
gtdb_phyla_family_samples<-as.data.frame(cbind(gtdb_all_tax_levels$phylum,gtdb_all_tax_levels$family, gtdb_all_tax_levels_df_norm))
colnames(gtdb_phyla_family_samples)[1:2]<-c("phylum","family")

gtdb_phylum_df_collapsed<-aggregate.data.frame(gtdb_phyla_family_samples[,4:14], by = list(gtdb_phyla_family_samples$phylum), FUN = sum)
sum(gtdb_phylum_df_collapsed$Siolta_1103_V2_D1)

gtdb_family_df_collapsed<-aggregate.data.frame(gtdb_phyla_family_samples[,4:14], by = list(gtdb_phyla_family_samples$family), FUN = sum)
sum(gtdb_family_df_collapsed$Siolta_1103_V2_D1)

## Combine all phylum sum tables

# Merge the phylum sums from StrainID and V3-V4
amplicon_phyla_sums<-merge(strainid_phyla_sum,v3_v4_phyla_sum,by="phylum", all = TRUE)

# Merge the phylum sums from metagenomics
metagenomics_phyla_sums<-merge(gtdb_phyla_sum,bracken_phyla_sum, by="phylum", all = TRUE)

# Merge the amplicon and metagenomics phyla sums
all_phyla_sums<-merge(amplicon_phyla_sums, metagenomics_phyla_sums, by="phylum", all = TRUE)

# Order by V3-V4 relative abundance, highest to lowest
all_phyla_sums_ordered<-all_phyla_sums[order(-all_phyla_sums$`V3-V4`),]

# Melt dataframe so the phylum counts are all in one column, with additional column of StrainID, V3-V4, long-read and short-read metagenomics
all_phyla_sums_melted<-melt(all_phyla_sums, measure.vars = c("StrainID","V3-V4","Long reads","Short reads"))
colnames(all_phyla_sums_melted)<-c("Phylum","Method","Relative abundance")

# Number of colors needed - 19
length(unique(all_phyla_sums_melted$Phylum))

# Specify 19 colors
phyla.19.colors<-c("#974661","#48A462","#FFB716","#4A72A6","#C28EA9","#E41A1C","#3E8E93","#5D995D","#7E6E85","#A35390","#D16948","#FF7F00","#FFF02D","#B97B2A","#B75F49","#DB728C","#EC83BA","#999999","#E1C62F")

# Change order of factors for plotting
all_phyla_sums_melted$Method<-factor(all_phyla_sums_melted$Method,levels=c("V3-V4","StrainID","Short reads","Long reads"))
levels(all_phyla_sums_melted$Method)

# Change the order of phyla for plotting
all_phyla_sums_melted$Phylum<-factor(all_phyla_sums_melted$Phylum, levels = all_phyla_sums_ordered$phylum)

# ggplot stacked barplot by phylum
ggplot(all_phyla_sums_melted,aes(fill=Phylum, y=`Relative abundance`, x=Method)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values = phyla.19.colors)

### Summarize by family level ###

## Read in StrainID family table
strainid_family<-read.table("input_files/StrainID_family.txt", sep = "\t", header = T, check.names = F)
grep("3218_V6",colnames(strainid_family))

# Get only the sample columns plus the first taxonomy column, no sample 3218_V6
strainid_family_df<-strainid_family[,c(3,6:49)]

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
strainid_family_df_norm<-adorn_percentages(strainid_family_df, denominator = "col")
# Check to make sure column sums are equal to 1
colSums(strainid_family_df_norm[,2:ncol(strainid_family_df_norm)])

# Add up total families across normalized samples
strainid_all_sample_family_totals<-rowSums(strainid_family_df_norm[,2:ncol(strainid_family_df_norm)])
strainid_family_totals_table<-as.data.frame(cbind(strainid_family_df_norm$taxon,strainid_all_sample_family_totals))

# Aggregate by family
strainid_family_sum<-aggregate(as.numeric(as.character(strainid_family_totals_table$strainid_all_sample_family_totals)),by = list(strainid_family_totals_table$V1), FUN=sum)
colnames(strainid_family_sum)<-c("family","StrainID")

# Get each sample's family sum separately for statistical tests
strainid_family_df_collapsed<-aggregate.data.frame(strainid_family_df_norm[,2:ncol(strainid_family_df_norm)], by = list(strainid_family_df_norm$taxon), FUN = sum)

## Summarize V3-V4 by family

# Aggregate by family
v3_v4_family_sum<-aggregate(as.numeric(as.character(v3_v4_all_tax_levels$eleven_samp_sum)),by = list(v3_v4_all_tax_levels$family), FUN=sum)
colnames(v3_v4_family_sum)<-c("family","V3-V4")

## Summarize Bracken short read metagenomic data by family

# Add up total families across normalized samples
bracken_all_sample_family_totals<-rowSums(bracken_all_tax_levels_df_norm[,2:ncol(bracken_all_tax_levels_df_norm)])
bracken_family_totals_table<-as.data.frame(cbind(bracken_all_tax_levels$family,bracken_all_sample_family_totals))

# Aggregate by family
bracken_family_sum<-aggregate(as.numeric(as.character(bracken_family_totals_table$bracken_all_sample_family_totals)),by = list(bracken_family_totals_table$V1), FUN=sum)
colnames(bracken_family_sum)<-c("family","Short reads")

## Summarize GTDB data by family

# Add up total family across normalized samples
gtdb_family_totals_table<-as.data.frame(cbind(gtdb_all_tax_levels$family,gtdb_all_sample_totals))

# Aggregate by family
gtdb_family_sum<-aggregate(as.numeric(as.character(gtdb_family_totals_table$gtdb_all_sample_totals)),by = list(gtdb_family_totals_table$V1), FUN=sum)
colnames(gtdb_family_sum)<-c("family","Long reads")

## Combine all tables

# Merge the family sums from StrainID and V3-V4
amplicon_family_sums<-merge(strainid_family_sum,v3_v4_family_sum,by="family", all = TRUE)

# Merge the family sums from metagenomics
metagenomics_family_sums<-merge(gtdb_family_sum,bracken_family_sum, by="family", all = TRUE)

# Merge the amplicon and metagenomics family sums
all_family_sums<-merge(amplicon_family_sums, metagenomics_family_sums, by="family", all = TRUE)

# Order by the relative abundance of V3-V4, highest to lowest
all_family_sums_ordered<-all_family_sums[order(-all_family_sums$`V3-V4`),]

## Melt dataframe so the family counts are all in one column, with additional column of StrainID, V3-V4, long-read and short-read metagenomics
all_family_sums_melted<-melt(all_family_sums_ordered, measure.vars = c("StrainID","V3-V4","Long reads","Short reads"))
colnames(all_family_sums_melted)<-c("Family","Method","Relative abundance")

# Number of colors needed - 163
length(unique(all_family_sums_melted$Family))

# Make color ramp palette with 163 colors
display.brewer.all()
Set1<-brewer.pal(9, "Set1")
Set1_ramp<-colorRampPalette(Set1)
fam.163.colors<-Set1_ramp(163)
set.seed(12)
fam.163.colors.rand<-sample(fam.163.colors)

# Show the 163 colors and their hex values
show_col(fam.163.colors)

# Select the first 24 colors for most abundant families
first_24_col<-c("#865070","#D77085","#3C899E","#FFCA1E","#599E59","#EA7520","#397CB6","#D8B52E","#FBF832","#D2232B","#6C866F","#B06A29","#CC8BAD","#4BAB52","#B0374A","#FFFD32","#C3655F","#526D9E","#F37FB9","#999999","#994EA0","#D4AD2D","#3D8B98","#FFC41B")
show_col(first_24_col)
length(unique(first_24_col))

# Get the remaining colors and combine with the first 24
remaining_col<-setdiff(fam.163.colors.rand,first_24_col)
col.163<-c(first_24_col,remaining_col)
length(col.163)
show_col(col.163)

# Make the Unknown family gray, and the Eubacteriacea family green
which(all_family_sums_ordered$family=="Unknown")
col.163[35]<-"#999999"
which(all_family_sums_ordered=="Eubacteriaceae")
col.163[31]<-"#449C72"

# Change order of factors for plotting methods
all_family_sums_melted$Method<-factor(all_family_sums_melted$Method,levels=c("V3-V4","StrainID","Short reads","Long reads"))
levels(all_family_sums_melted$Method)

# Change factor order for plotting families, with most abundant families at the top
all_family_sums_melted$Family<-factor(all_family_sums_melted$Family, levels = all_family_sums_ordered$family)
factor(all_family_sums_melted$Family)

## ggplot stacked barplot by family
ggplot(all_family_sums_melted,aes(fill=Family, y=`Relative abundance`, x=Method)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(breaks=c("Lachnospiraceae","Ruminococcaceae","Erysipelotrichaceae","Bifidobacteriaceae","Bacteroidaceae","Clostridiales unclassified","Akkermansiaceae","Veillonellaceae","Coriobacteriaceae","Clostridiaceae","Prevotellaceae","Streptococcaceae","Methanobacteriaceae","Porphyromonadaceae","Enterobacteriaceae","Rikenellaceae","Peptostreptococcaceae","Christensenellaceae","Eubacteriaceae","Unknown","Selenomonadaceae","Eggerthellaceae","Acidaminococcaceae","UBA1381"),values = col.163)

### Write out tables of phyla and family abundance by sample to see if levels across groups are significantly different ###

# Phylum level summaries by sample
write.table(strainid_phylum_df_norm_collapsed,"output_files/strainid_phylum_proportion_by_sample.txt", sep = "\t")
write.table(v3_v4_phylum_df_collapsed, "output_files/v3_v4_phylum_proportion_by_sample.txt", sep = "\t")
write.table(bracken_phylum_df_collapsed, "output_files/bracken_phylum_proportion_by_sample.txt", sep = "\t")
write.table(gtdb_phylum_df_collapsed, "output_files/gtdb_phylum_proportion_by_sample.txt", sep = "\t")

# Family level summaries by sample
write.table(strainid_family_df_collapsed, "output_files/straind_family_proportion_by_sample.txt", sep = "\t")
write.table(v3_v4_family_df_collapsed, "output_files/v3_v4_family_proportion_by_sample.txt", sep = "\t")
write.table(bracken_family_df_collapsed, "output_files/bracken_family_proportion_by_sample.txt", sep = "\t")
write.table(gtdb_family_df_collapsed, "output_files/gtdb_family_proportion_by_sample.txt", sep = "\t")

### Perform differential abundance at family level ###

## For StrainID family table, get only one replicate from each sample plus the first taxonomy column
col_list<-c("taxon",strainid_metadata_sl_most_reads$Sample)
strainid_family_df_12<-strainid_family[,col_list]

# Get the rows (families) with more than 1 read
strainid_family_df_12_nonzero<-strainid_family_df_12[rowSums(strainid_family_df_12[,2:ncol(strainid_family_df_12)])>1,]

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
strainid_family_df_12_norm<-adorn_percentages(strainid_family_df_12_nonzero, denominator = "col")
# Check to make sure column sums are equal to 1
colSums(strainid_family_df_12_norm[,2:ncol(strainid_family_df_12_norm)])

# Get each sample's family sum
strainid_family_df_12_collapsed<-aggregate.data.frame(strainid_family_df_12_norm[,2:ncol(strainid_family_df_12_norm)], by = list(strainid_family_df_12_norm$taxon), FUN = sum)

# Format family abundance table for maaslin2: features as columns, samples as rows
rownames(strainid_family_df_12_collapsed)<-strainid_family_df_12_collapsed$Group.1
strainid_family_df_12_collapsed$Group.1<-NULL
strainid_family_df_12_collapsed_t<-t(strainid_family_df_12_collapsed)

# Remove replicate designation from sample names
rownames(strainid_family_df_12_collapsed_t)<-gsub('.{3}$','',rownames(strainid_family_df_12_collapsed_t))

# Make general metadata table for maaslin2: features as columns, samples as rows
metadata_all<-v3_v4_metadata
metadata_all$DNA_extraction_method<-NULL
metadata_all$Reads<-NULL
metadata_all$SampleName<-NULL
rownames(metadata_all)<-metadata_all$Sample

# Run maaslin2 - StrainID family level
strainid_family_diff_abund<-Maaslin2(strainid_family_df_12_collapsed_t,metadata_all,'output_files/strainid_family_diff_abund_baseline_treatment', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, max_significance=0.2, cores = 2, min_prevalence = 0)

## Aggregate V3-V4 data by family
v3_v4_family_df_12_collapsed<-aggregate.data.frame(v3_v4_all_tax_levels[,2:13], by = list(v3_v4_all_tax_levels$family), FUN = sum)

# Format family abundance table for maaslin2: features as columns, samples as rows
rownames(v3_v4_family_df_12_collapsed)<-v3_v4_family_df_12_collapsed$Group.1 
v3_v4_family_df_12_collapsed$Group.1<-NULL
v3_v4_family_df_12_collapsed_t<-t(v3_v4_family_df_12_collapsed)

# Run maaslin2 - V3-V4 family level
v3_v4_family_diff_abund<-Maaslin2(v3_v4_family_df_12_collapsed_t,metadata_all,'output_files/v3_v4_family_diff_abund_baseline_treatment', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, max_significance=0.2, cores = 2, min_prevalence = 0)

## Aggregate Bracken data by family
bracken_family_df_12_collapsed<-aggregate.data.frame(bracken_all_tax_levels[,2:13], by = list(bracken_all_tax_levels$family), FUN = sum)
colSums(bracken_family_df_12_collapsed[,2:ncol(bracken_family_df_12_collapsed)])

# Divide each sample by total number of reads in that sample. First column not included with adorn_percentages
bracken_family_df_12_collapsed_norm<-adorn_percentages(bracken_family_df_12_collapsed, denominator = "col")
# Check to make sure column sums are equal to 1
colSums(bracken_family_df_12_collapsed_norm[,2:ncol(bracken_family_df_12_collapsed_norm)])

# Format family abundance table for maaslin2: features as columns, samples as rows
rownames(bracken_family_df_12_collapsed_norm)<-bracken_family_df_12_collapsed_norm$Group.1
bracken_family_df_12_collapsed_norm$Group.1<-NULL
bracken_family_df_12_collapsed_norm_t<-t(bracken_family_df_12_collapsed_norm)

# Run maaslin2 - Bracken family level
bracken_family_diff_abund<-Maaslin2(bracken_family_df_12_collapsed_norm_t,metadata_all,'output_files/bracken_family_diff_abund_baseline_treatment', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, max_significance=0.2, cores = 2, min_prevalence = 0)

## Format GTDB family abundance table for maaslin2: features as columns, samples as rows
rownames(gtdb_family_df_collapsed)<-gtdb_family_df_collapsed$Group.1
gtdb_family_df_collapsed$Group.1<-NULL
gtdb_family_df_collapsed_t<-t(gtdb_family_df_collapsed)
# Fix sample (row) names
rownames(gtdb_family_df_collapsed_t)<-gsub("Siolta_","",rownames(gtdb_family_df_collapsed_t))

# Run maaslin2 - GTDB family level
gtdb_family_diff_abund<-Maaslin2(gtdb_family_df_collapsed_t,metadata_all,'output_files/gtdb_family_diff_abund_baseline_treatment', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, max_significance=0.2, cores = 2, min_prevalence = 0)

### Heatmap of all families and the coefficients by method ###

# Read in file of concatenated all_results.tsv files output from maaslin2
diff_abund_fam_all<-read.table("input_files/Differentially_abundant_families_all_methods.txt", header = T, sep = "\t")

# Get only the families that were significantly differentially abundant with at least one method
diff_abund_fam_all_sig_gtdb<-diff_abund_fam_all[diff_abund_fam_all$pval<0.25&diff_abund_fam_all$method=="GTDB",]
diff_abund_fam_all_sig_other<-diff_abund_fam_all[diff_abund_fam_all$pval<0.2&diff_abund_fam_all$method!="GTDB",]
diff_abund_fam_all_sig<-rbind(diff_abund_fam_all_sig_gtdb,diff_abund_fam_all_sig_other)

# Make a matrix of family by method with the coefficients as values, for only the differentially abundant families
diff_abund_fam_all_sig_simp<-aggregate.data.frame(diff_abund_fam_all_sig$coef, by = list(diff_abund_fam_all_sig$feature,diff_abund_fam_all_sig$method), FUN = mean)
diff_abund_fam_all_sig_simp_dcast<-dcast(diff_abund_fam_all_sig_simp, Group.1 ~ Group.2)
rownames(diff_abund_fam_all_sig_simp_dcast)<-diff_abund_fam_all_sig_simp_dcast$Group.1
diff_abund_fam_all_sig_simp_dcast$Group.1<-NULL

# Get list of differentially abundant families (in at least one method)
diff_abund_fams<-rownames(diff_abund_fam_all_sig_simp_dcast)

# Get all coefficient and p-value data for all methods for these differentially abundant families
diff_abund_fam_all_methods<-subset(diff_abund_fam_all, diff_abund_fam_all$feature %in% diff_abund_fams)

# Make a matrix of family by method with the coefficients as values, with all methods' values for differentially abundant families 
diff_abund_fam_all_methods_simp<-aggregate.data.frame(diff_abund_fam_all_methods$coef, by = list(diff_abund_fam_all_methods$feature,diff_abund_fam_all_methods$method), FUN = mean)
diff_abund_fam_all_methods_simp_dcast<-dcast(diff_abund_fam_all_methods_simp, Group.1 ~ Group.2)
rownames(diff_abund_fam_all_methods_simp_dcast)<-diff_abund_fam_all_methods_simp_dcast$Group.1
diff_abund_fam_all_methods_simp_dcast$Group.1<-NULL

# Repeat, but make matrix of p-values
# Make a matrix of family by method with the coefficients as values, with all methods' values for differentially abundant families 
diff_abund_fam_all_methods_pval<-aggregate.data.frame(diff_abund_fam_all_methods$pval, by = list(diff_abund_fam_all_methods$feature,diff_abund_fam_all_methods$method), FUN = mean)
diff_abund_fam_all_methods_pval_dcast<-dcast(diff_abund_fam_all_methods_pval, Group.1 ~ Group.2)
rownames(diff_abund_fam_all_methods_pval_dcast)<-diff_abund_fam_all_methods_pval_dcast$Group.1
diff_abund_fam_all_methods_pval_dcast$Group.1<-NULL

# Convert significant p-values to asterisks
bracken_sig<-which(diff_abund_fam_all_methods_pval_dcast$Bracken<0.2)
bracken_nonsig<-which(diff_abund_fam_all_methods_pval_dcast$Bracken>0.2)

v3_v4_sig<-which(diff_abund_fam_all_methods_pval_dcast$`V3-V4`<0.2)
v3_v4_nonsig<-which(diff_abund_fam_all_methods_pval_dcast$`V3-V4`>0.2)

strainid_sig<-which(diff_abund_fam_all_methods_pval_dcast$StrainID<0.2)
strainid_nonsig<-which(diff_abund_fam_all_methods_pval_dcast$StrainID>0.2)

gtdb_sig<-which(diff_abund_fam_all_methods_pval_dcast$GTDB<0.25)
gtdb_nonsig<-which(diff_abund_fam_all_methods_pval_dcast$GTDB>0.25)

diff_abund_fam_all_methods_pval_dcast$Bracken[bracken_sig]<-"*"
diff_abund_fam_all_methods_pval_dcast$Bracken[bracken_nonsig]<-""

diff_abund_fam_all_methods_pval_dcast$`V3-V4`[v3_v4_sig]<-"*"
diff_abund_fam_all_methods_pval_dcast$`V3-V4`[v3_v4_nonsig]<-""

diff_abund_fam_all_methods_pval_dcast$StrainID[strainid_sig]<-"*"
diff_abund_fam_all_methods_pval_dcast$StrainID[strainid_nonsig]<-""

diff_abund_fam_all_methods_pval_dcast$GTDB[gtdb_sig]<-"*"
diff_abund_fam_all_methods_pval_dcast$GTDB[gtdb_nonsig]<-""

# Convert table of coefficients to a matrix
diff_abund_fam_all_methods_simp_dcast_mat<-as.matrix(diff_abund_fam_all_methods_simp_dcast)

# Replace period in family names with space
rownames(diff_abund_fam_all_methods_simp_dcast_mat)<-gsub("."," ",rownames(diff_abund_fam_all_methods_simp_dcast_mat),fixed = TRUE)
rownames(diff_abund_fam_all_methods_pval_dcast)<-gsub("."," ",rownames(diff_abund_fam_all_methods_pval_dcast),fixed = TRUE)

# Remove "NAs" from p value table
diff_abund_fam_all_methods_pval_dcast_no_NA<-diff_abund_fam_all_methods_pval_dcast
diff_abund_fam_all_methods_pval_dcast_no_NA[is.na(diff_abund_fam_all_methods_pval_dcast_no_NA)]<-""

# Make heatmap of differentially abundant families. Gray cells are NA, indicating that the family was not part of the method's dataset
pheatmap(diff_abund_fam_all_methods_simp_dcast_mat, display_numbers = diff_abund_fam_all_methods_pval_dcast_no_NA, cluster_rows=FALSE,cluster_cols = TRUE)

# Make heatmap of differentially abundant families. Cells marked NA indicate that the family was not part of the method's dataset
# Replace NAs with 0 to cluster
diff_abund_fam_all_methods_simp_dcast_mat_no_na<-diff_abund_fam_all_methods_simp_dcast_mat
diff_abund_fam_all_methods_simp_dcast_mat_no_na[is.na(diff_abund_fam_all_methods_simp_dcast_mat_no_na)]<-0

# Change color palette and make breaks so that 0 (or NA) is white
paletteLength<-50
myColor <- colorRampPalette(c("cornflowerblue","white","tomato"))(paletteLength)
myBreaks <- c(seq(min(diff_abund_fam_all_methods_simp_dcast_mat_no_na), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(diff_abund_fam_all_methods_simp_dcast_mat_no_na)/paletteLength, max(diff_abund_fam_all_methods_simp_dcast_mat_no_na), length.out=floor(paletteLength/2)))

pheatmap(diff_abund_fam_all_methods_simp_dcast_mat_no_na, display_numbers = diff_abund_fam_all_methods_pval_dcast,color = myColor, breaks=myBreaks)

### Plot differentially abundant families ###

## Merge metadata and StrainID family table
# Make matrix a dataframe
strainid_family_table<-as.data.frame(strainid_family_df_12_collapsed_t)
# Remove replicate designation from sample names
rownames(strainid_family_table)<-gsub('.{3}$','',rownames(strainid_family_table))
# Change "Timepoint" to "Time point" for graphs
metadata_all_plots<-metadata_all
colnames(metadata_all_plots)[6]<-"Time point"
# Add in Sample column to StrainID family table and merge with metadata to plot
strainid_family_table$Sample<-rownames(strainid_family_table)
strainid_family_table_meta<-merge(metadata_all_plots,strainid_family_table, by="Sample")

# Plot differentially abundant StrainID families for families that were differentially abundant with two or more methods
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Akkermansiaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Bacteroidaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Oscillospiraceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Christensenellaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Clostridiales Family XIII`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Streptococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Enterobacteriaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(strainid_family_table_meta, aes(x=`Time point`, y=`Micrococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))

## Merge metadata and V3-V4 family table
# Make matrix a dataframe
v3_v4_family_table<-as.data.frame(v3_v4_family_df_12_collapsed_t)
# Add in Sample column to V3-V4 family table and merge with metadata to plot
v3_v4_family_table$Sample<-rownames(v3_v4_family_table)
v3_v4_family_table_meta<-merge(metadata_all_plots,v3_v4_family_table, by="Sample")

# Plot differentially abundant V3-V4 families for families that were differentially abundant with two or more methods
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Akkermansiaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Bacteroidaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Christensenellaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Clostridiales Family XIII`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Streptococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Enterobacteriaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(v3_v4_family_table_meta, aes(x=`Time point`, y=`Micrococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))

## Merge metadata and short-read metagenomics (bracken) family table
# Make matrix a dataframe
bracken_family_table<-as.data.frame(bracken_family_df_12_collapsed_norm_t)
# Add in Sample column to bracken family table and merge with metadata to plot
bracken_family_table$Sample<-rownames(bracken_family_table)
bracken_family_table_meta<-merge(metadata_all_plots,bracken_family_table, by="Sample")

# Plot differentially abundant short-read metagenomics (bracken) families for families that were differentially abundant with two or more methods
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Akkermansiaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Bacteroidaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Oscillospiraceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Christensenellaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Clostridiales Family XIII`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Streptococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Enterobacteriaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(bracken_family_table_meta, aes(x=`Time point`, y=`Micrococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))

## Merge metadata and long-read metagenomics (GTDB) family table
# Make matrix a dataframe
gtdb_family_table<-as.data.frame(gtdb_family_df_collapsed_t)
# Add in Sample column to GTDB family table and merge with metadata to plot
gtdb_family_table$Sample<-rownames(gtdb_family_table)
gtdb_family_table_meta<-merge(metadata_all_plots,gtdb_family_table, by="Sample")

# Plot differentially abundant long-read metagenomics (GTDB) families for families that were differentially abundant with two or more methods
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Akkermansiaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Bacteroidaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Oscillospiraceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Clostridiales unclassified`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Streptococcaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(gtdb_family_table_meta, aes(x=`Time point`, y=`Enterobacteriaceae`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))


### Average the coverage of short metagenomic reads mapped to Strain 1 genome for coverage plots ###

# Read in sorted .paf files generated from mapping Participant 3 and Participant 4's metagenomic reads files to Strain 1's genome using minimap2
P3_treatment_Strain1_coverage<-read.table("input_files/P3_treatment_Strain1_mapped_sorted.txt", header = F, sep = "\t")
P4_treatment_Strain1_coverage<-read.table("input_files/P4_treatment_Strain1_mapped_sorted.txt", header = F, sep = "\t")

# Add column names
colnames(P3_treatment_Strain1_coverage)<-c("Reference", "Reference Start Position", "Coverage")
colnames(P4_treatment_Strain1_coverage)<-c("Reference", "Reference Start Position", "Coverage")

# Bin coverage and get the average of each 1000 bp bin, write out tables
P3_treatment_Strain1_coverage_1000bp_mean<-aggregate(P3_treatment_Strain1_coverage$Coverage, by=list(cut(P3_treatment_Strain1_coverage$`Reference Start Position`, seq(1000,2741000, by=1000))), mean)
write.table(output_files/P3_treatment_Strain1_coverage_1000bp_mean, "P3_treatment_Strain1_coverage_1000bp_mean.txt", sep = "\t")

P4_treatment_Strain1_coverage_1000bp_mean<-aggregate(P4_treatment_Strain1_coverage$Coverage, by=list(cut(P4_treatment_Strain1_coverage$`Reference Start Position`, seq(1000,2741000, by=1000))), mean)
write.table(output_files/P4_treatment_Strain1_coverage_1000bp_mean, "P4_treatment_Strain1_coverage_1000bp_mean.txt", sep = "\t")


### Functional microbiome profiles with metagenomic methods ###

## Read in humann2 pathway abundance data (short-read metagenomics)
humann2.pathabund<-read.table("input_files/combined_cohorts_12samp_humann2_pathabund_relab_unstratified.txt", sep = "\t", check.names = FALSE, header = TRUE)

# Make the rownames the pathway names
rownames(humann2.pathabund)<-humann2.pathabund$Pathway
humann2.pathabund$Pathway<-NULL

# Transpose pathway abundance table so samples are the rows and pathways are the columns
humann2.pathabund.t<-t(humann2.pathabund)

# Run differential abundance on humann2 pathway data
humann2.pathabund.diff.abund<-Maaslin2(humann2.pathabund.t,metadata_all,'output_files/humann2.pathabund.diff.abund.baseline.treatment.no.transform.TSS', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, transform = "NONE", max_significance=0.2, cores = 2, min_prevalence = 0)

## Read in SEED functional counts (long-read metagenomics)
seed.functional<-read.table("input_files/Protein-Functional-SEED-NormalizedCounts.siolta.mod.txt", sep = "\t", header = T, check.names = FALSE)

# Make first column of the SEED functions the rownames
rownames(seed.functional)<-seed.functional$Datasets
seed.functional$Datasets<-NULL

# Transpose SEED table so the samples are the rows and the SEED functions are the columns
seed.functional.t<-t(seed.functional)

# Run differential abundance on SEED functional data
seed.diff.abund<-Maaslin2(seed.functional.t,metadata_all,'output_files/seed.functional.diff.abund.baseline.treatment.no.transform.TSS', fixed_effects = c('Timepoint'), reference = 'Timepoint,baseline' ,random_effects = c('Patient'),standardize = FALSE, transform = "NONE", max_significance=0.2, cores = 2, min_prevalence = 0)

### Plot of differentially abundant pathways ####

# Combine humann2 pathway abundance table with metadata
humann2.pathabund.t.df<-as.data.frame(humann2.pathabund.t)
humann2.pathabund.t.df.samples<-cbind(rownames(humann2.pathabund.t.df),humann2.pathabund.t.df)
colnames(humann2.pathabund.t.df.samples)[1]<-"Sample"
humann2.functional.norm.meta<-merge(metadata_all_plots, humann2.pathabund.t.df.samples, by = "Sample")

# Plot differentially abundant pathways - humann2 pathways

ggplot(humann2.functional.norm.meta, aes(x=`Time point`, y=`PWY-7383: anaerobic energy metabolism (invertebrates. cytosol)`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(humann2.functional.norm.meta, aes(x=`Time point`, y=`HSERMETANA-PWY: L-methionine biosynthesis III`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(humann2.functional.norm.meta, aes(x=`Time point`, y=`DENOVOPURINE2-PWY: superpathway of purine nucleotides de novo biosynthesis II`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))

## Normalize SEED table (divide by total counts for each sample)
seed.functional.t.norm<-seed.functional.t/rowSums(seed.functional.t)
rowSums(seed.functional.t.norm)

# Combine SEED table with metadata
seed.functional.t.norm.df<-as.data.frame(seed.functional.t.norm)
seed.functional.t.norm.samples<-cbind(rownames(seed.functional.t.norm.df),seed.functional.t.norm.df)
colnames(seed.functional.t.norm.samples)[1]<-"Sample"
seed.functional.norm.meta<-merge(metadata_all_plots, seed.functional.t.norm.samples, by = "Sample")

# Plot differentially abundant pathways - SEED data
ggplot(seed.functional.norm.meta, aes(x=`Time point`, y=`Cluster containing Glutathione synthetase`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))
ggplot(seed.functional.norm.meta, aes(x=`Time point`, y=`Proline. 4-hydroxyproline uptake and utilization`)) + geom_boxplot() + geom_dotplot(binaxis='y',stackdir = 'center') + theme(axis.text = element_text(size = 18,colour = "black"), axis.title = element_text(size = 20,colour = "black",face = "bold"))


### Correlations between Picrust2 E.C. data and PacBio long-read E.C. data ###

picrust_ec<-read.table("input_files/pred_metagenome_unstrat.txt", sep = "\t", header = T, check.names = F)

pacbio_ec<-read.table("input_files/Functional-EC-absolute-read-counts-siolta-mod.txt", sep = "\t", header = T, check.names = F)

# Merge EC tables, keeping all ECs in both tables
ec_all<-merge(picrust_ec,pacbio_ec, by="function", all = T)

# Replace NAs with 0
ec_all[is.na(ec_all)]<-0

# Make function the rowname
rownames(ec_all)<-ec_all$`function`
ec_all$`function`<-NULL

# Remove rows (ECs) that have less than 1 count
ec_all_nonzero<-ec_all[rowSums(ec_all)>1,]

# Calculate Pearson correlations between all samples
rcorr_ec<-rcorr(as.matrix(ec_all_nonzero), type = c("pearson"))

## Flatten correlation matrix: from http://www.sthda.com/english/wiki/correlation-matrix-formatting-and-visualization

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

rrcor_ec_flat<-flattenCorrMatrix(rcorr_ec$r, rcorr_ec$P)

# Write out correlations to plot
write.table(rrcor_ec_flat,"output_files/Pearson_correlations_picrust_and_pacbio_ec.txt", sep = "\t")


### Make a grouped bar plot of the top 20 most represented genera of MAGs for short reads and long reads ###

# Read in file with the top 20 genera from MAGs from both short read and long read sequencing data. This file is derived from the summary excel file in the input_files folder, "All samples all genome bins", which contains all MAGs that passed the filtering criteria
top20_genera<-read.table("input_files/Top20_genera_MAGs.txt", sep = "\t", header = T, check.names = F)

# Melt dataframe and change column names
top20_genera_melted<-melt(top20_genera, measure.vars = c("Long reads","Short reads"))
colnames(top20_genera_melted)<-c("Genus","Method","Number of MAGs")

# Make color ramp palette with 20 colors
phyla.20.colors<-Set1_ramp(20)

# Change order of factors for plotting
# Sort from genera with greatest number of MAGs to least
top20_genera_melted$Genus<-factor(top20_genera_melted$Genus,levels=c("Ruminococcus","Blautia","Faecalibacterium","Bifidobacterium","Agathobacter","Gemmiger","Collinsella","Eubacterium","Dorea","Mediterraneibacter","Anaerobutyricum","Erysipelatoclostridium","Fusicatenibacter","Bacteroides","Coprococcus","Anaerostipes","Alistipes","Dialister","Roseburia","CAG-41 (o. Monoglobales)"))
# Plot short reads on left and long reads on right
top20_genera_melted$Method<-factor(top20_genera_melted$Method,levels=c("Short reads","Long reads"))

# Generate barplot
ggplot(top20_genera_melted,aes(fill=Genus, y=`Number of MAGs`, x=Method)) + geom_bar(position = "dodge", stat = "identity") + scale_fill_manual(values = phyla.20.colors)


