#Script for counting reads in features
directory = "/drives/ssd1/manuel/phaw/2021_analysis/phaw_mzt/rna-seq/zga_mzt_analysis/"
setwd(directory)

library(Rsubread)
#Call .bm sorted files polyA
files_directory <- "/drives/ssd1/manuel/phaw/2021_analysis/phaw_reference/bams"
bamFiles <- grep("S[1-9].bam$",list.files(files_directory,full.names = T),value=TRUE)

s1s9_gene_counts <- featureCounts(bamFiles,
                                     annot.ext =  "/drives/ssd1/manuel/phaw/2021_analysis/annotation/stringtie2/stringtie2_merge/stringtie2_merged_unigenes.gtf",
                                     isGTFAnnotationFile= TRUE,
                                     GTF.featureType = "transcript",
                                     GTF.attrType = "gene_id",
                                     useMetaFeatures= TRUE,
                                     allowMultiOverlap = FALSE,
                                     minOverlap = 2,
                                     strandSpecific = 2, #CHECK THIS PARAMETER
                                     isPairedEnd= TRUE,
                                     nthreads = 8)
s1s9_exon_counts <- featureCounts(bamFiles,
                                  annot.ext =  "/drives/ssd1/manuel/phaw/2021_analysis/annotation/stringtie2/stringtie2_merge/stringtie2_merged_unigenes.gtf",
                                  isGTFAnnotationFile= TRUE,
                                  GTF.featureType = "exon",
                                  GTF.attrType = "gene_id",
                                  useMetaFeatures= TRUE,
                                  allowMultiOverlap = FALSE,
                                  minOverlap = 2,
                                  strandSpecific = 2, #CHECK THIS PARAMETER
                                  isPairedEnd= TRUE,
                                  nthreads = 8)

s1s9_intron_counts = s1s9_gene_counts$counts - s1s9_exon_counts$counts
s1s9_intron_counts <- s1s9_intron_counts[,c(1,4,7,2,5,8,3,6,9)]
s1s9_intron_counts=pmax(s1s9_intron_counts,0)

#Call .bam sorted files Total
files_directory <- "/drives/ssd1/manuel/phaw/2021_analysis/phaw_reference/bams/embryo_totalRNA"
bamFiles <- grep("T[1-9].bam$",list.files(files_directory,full.names = T),value=TRUE)

t1t9_gene_counts <- featureCounts(bamFiles,
                                  annot.ext =  "/drives/ssd1/manuel/phaw/2021_analysis/annotation/stringtie2/stringtie2_merge/stringtie2_merged_unigenes.gtf",
                                  isGTFAnnotationFile= TRUE,
                                  GTF.featureType = "transcript",
                                  GTF.attrType = "gene_id",
                                  useMetaFeatures= TRUE,
                                  allowMultiOverlap = FALSE,
                                  minOverlap = 2,
                                  strandSpecific = 2, #CHECK THIS PARAMETER
                                  isPairedEnd= TRUE,
                                  nthreads = 8)

t1t9_exon_counts <- featureCounts(bamFiles,
                             annot.ext =  "/drives/ssd1/manuel/phaw/2021_analysis/annotation/stringtie2/stringtie2_merge/stringtie2_merged_unigenes.gtf",
                             isGTFAnnotationFile= TRUE,
                             GTF.featureType = "exon",
                             GTF.attrType = "gene_id",
                             useMetaFeatures= TRUE,
                             allowMultiOverlap = FALSE,
                             minOverlap = 2,
                             strandSpecific = 2, #CHECK THIS PARAMETER
                             isPairedEnd= TRUE,
                             nthreads = 8)

t1t9_intron_counts = t1t9_gene_counts$counts - t1t9_exon_counts$counts
t1t9_intron_counts <- t1t9_intron_counts[,c(1,4,7,2,5,8,3,6,9)]
t1t9_intron_counts=pmax(t1t9_intron_counts,0)

#write.table(gene_feature_counts$counts,"fc_counts.txt",quote = F,col.names = T,sep="\t",row.names = T)
count_matrix= read.table("fc_counts.txt",header = T,sep = "\t",stringsAsFactors = F)
