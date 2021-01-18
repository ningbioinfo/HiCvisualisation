## This script is used to generate integration-tracks plot showed in the manuscript "3DFAACTS-SNP: Using regulatory T cell-specific epigenomics data to uncover candidate mechanisms of Type-1 Diabetes (T1D) risk".
# Author: Ning Liu

library(Gviz)
library(GenomicInteractions)
library(coMET)
setwd('~/Documents/T1DPaper/PublicCodes/')

grepn <- function(pattern,x){
  y <- x[grep(pattern,x)] 
  return(y)
}

## This script is used to generate integration-tracks plot showed in the manuscript "3DFAACTS-SNP: Using regulatory T cell-specific epigenomics data to uncover candidate mechanisms of Type-1 Diabetes (T1D) risk".

# Here in this script, we aim at generating Figure S2 from Additional File 1.

# Before running this script, you should run the T1D_3DFAACTS-SNP_workflow.R to obtain the SNPs, enhancers, Hi-C data, ATAC-seq data, FOXP3 data that are used for plotting.

# first need to decide which Hi-C interaction should we include, because there is no point to plot all interactions and they can be very messy. In this case, we need Hi-C interactions that capturing our identified 3DFAACTS SNPs.

hic2plot <- interaction_intersecting_elements(interaction = treg_hic, gr_elements = t1d_3dfaacts_snp)

# then we need to output the Hi-C interactions into a file with bedpe format becuase the function from GenomicInteractions can only read from file.

hic2plot %>%
  mutate(Interactions = paste0("Interactions",seq_len(nrow(.))), int_count = 1, ex1 = 1) %>%
  write_delim('./hic2plot.bedpe', delim = '\t', col_names = F)

hic.track <- InteractionTrack(makeGenomicInteractionsFromFile('./hic2plot.bedpe',type = "bedpe", experiment_name = "HiC of SNPs"), name = "HiC Tregs")

displayPars(hic.track) = list(col.interactions="purple",
                              col.anchors.fill ="blue",
                              col.anchors.line = "black",
                              interaction.dimension="height",
                              interaction.measure ="counts",
                              plot.trans=FALSE,
                              plot.outside = TRUE,
                              col.outside="green",
                              anchor.height = 0.05,
                              collapse = FALSE)

# make ruler
gatrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE)

# load Super-enhancers

SE_files <- list.files("./PublicData4GvizPlotting/SE/") %>%
  grepn("^0.*SE.bed$",.)

SE_names <- c("B-cell.SE", "CD4-T.SE", "CD4mem-T.SE", "CD8-T-cell.SE", "Treg.SE", "CD4-T-cell.SE")

SE_grs <- SE_files %>%
  as.list() %>%
  lapply(function(x) import.bed(paste0("./PublicData4GvizPlotting/SE/",x)))

SE1 <- AnnotationTrack(SE_grs[[1]], name = SE_names[1], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)
SE2 <- AnnotationTrack(SE_grs[[2]], name = SE_names[2], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)
SE3 <- AnnotationTrack(SE_grs[[3]], name = SE_names[3], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)
SE4 <- AnnotationTrack(SE_grs[[4]], name = SE_names[4], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)
SE5 <- AnnotationTrack(SE_grs[[5]], name = SE_names[5], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)
SE6 <- AnnotationTrack(SE_grs[[6]], name = SE_names[6], width = 15, genome = "hg19", fill="#5D6D7E", col = NULL)

# load chromhmm

chromHMMid <- fread("./PublicData4GvizPlotting/chromHMM/chromhmm_ID.txt")

get_chromHMM_tracks <- function(chr, start, end){
  chromHMM_Treg <<- chromHMM_RoadMap(chr = chr, start = start, end = end, bedFilePath = file.path("./PublicData4GvizPlotting/chromHMM/E044_chromHMM.bed.gz"), featureDisplay = "all", colorcase='roadmap15', title = "Tregs")
  chromHMM_Th17 <<- chromHMM_RoadMap(chr = chr, start = start, end = end, bedFilePath = file.path("./PublicData4GvizPlotting/chromHMM/E042_chromHMM.bed.gz"), featureDisplay = "all", colorcase='roadmap15', title = "Th17")
  chromHMM_priT <<- chromHMM_RoadMap(chr = chr, start = start, end = end, bedFilePath = file.path("./PublicData4GvizPlotting/chromHMM/E034_chromHMM.bed.gz"), featureDisplay = "all", colorcase='roadmap15', title = "PrimaryT")
  chromHMM_priB <<- chromHMM_RoadMap(chr = chr, start = start, end = end, bedFilePath = file.path("./PublicData4GvizPlotting/chromHMM/E032_chromHMM.bed.gz"), featureDisplay = "all", colorcase='roadmap15', title = "PrimaryB")
  chromHMM_priHSC <<- chromHMM_RoadMap(chr = chr, start = start, end = end, bedFilePath = file.path("./PublicData4GvizPlotting/chromHMM/E035_chromHMM.bed.gz"), featureDisplay = "all", colorcase='roadmap15', title = "PrimaryHSC")
}

# set chromosome, start, end to plot, here we aim at generating Figure S2 from Additional File 1.

c <- "chr2"
s <- 204122714
e <- 204812714

# function to add other tracks
# in order to execute this function, you need to download the RNA-seq data (ENA Project id: PRJEB11844 and sample id: SQ_0350 and SQ_0351), which can also be found in the github repo.

get_track_elements <- function(chr, start, end){
  ideoTrack <<- IdeogramTrack(genome = "hg19", chromosome = chr, from = start ,to = end , showId = TRUE, bands = NULL)

  # gene
  txdb_hg19 <<- TxDb.Hsapiens.UCSC.hg19.knownGene
  knownGenes <<- GeneRegionTrack(txdb_hg19, genome="hg19", chromosome=chr, start = start, end = end, showId=TRUE, geneSymbol=TRUE, name="UCSC Transcripts", shape = "arrow", collapseTranscripts = F)
  if(length(gene(knownGenes)) != 0){
    sss <<- unlist(mapIds(org.Hs.eg.db, gene(knownGenes), "SYMBOL", "ENTREZID", multiVals = "first"))
    symbol(knownGenes) <<- sss[gene(knownGenes)]
  }

  # alignment

  bam1 <<- DataTrack(range="./PublicData4GvizPlotting/RNAseq/SQ_0350.bedgraph", genome="hg19", type=c("heatmap"), chromosome=chr, name="Th1 RNAseq", gradient = colorRampPalette(c("white","orange", "red"))(10)[1:10])
  bam2 <<- DataTrack(range="./PublicData4GvizPlotting/RNAseq/SQ_0351.bedgraph", genome="hg19", type=c("heatmap"), chromosome=chr, name="Treg RNAseq", gradient = colorRampPalette(c("white","orange", "red"))(10)[1:10])

  #atac
  dtrack.atac <<- DataTrack(atac_treg, chromosome = chr,
                            start = start, end = end, type="polygon",
                            name = "Tregs ATACseq", col = "orange")

  #foxp3

  dtrack.foxp3 <<- DataTrack(foxp3, chromosome = chr,
                             start = start, end = end,
                             type = "histogram", name = "FOXP3 ChIP-chip", fill.histogram = "blue", col.histogram = "blue")

  #  dtrack.ctcf <<- DataTrack(ctcf, chromosome = chr,
  #                            start = start, end = end,
  #                            type = "histogram", name = "CTCF NaiveT ChIP-seq", fill.histogram = "gray", col.histogram = "gray")

  #promoters
  promoters_ns <- promoters %>% subsetByOverlaps(GRanges(seqnames=chr, ranges=IRanges(start = start, end = end)))
  p.track <<- AnnotationTrack(promoters_ns, fill = "#B61456", col = NULL, name = "Promoters", stacking = "dense")

  #enhancers
  enhancers_ns <- enhancers %>% subsetByOverlaps(GRanges(seqnames=chr, ranges=IRanges(start = start, end = end)))
  e.track <<- AnnotationTrack(enhancers_ns, fill = "#80830C", col = NULL, name = "Enhancers")


  #, id = ifelse(risks$rsid %in% potential_causative_snps$rsid, risks$rsid, "")
  #risks
  snp.track <<- AnnotationTrack(t1d_3dfaacts_snp, name = "3DFAACTS-SNPs", col = NULL, fill = "red",shape = "box")

}

get_track_elements(c,s,e)
get_chromHMM_tracks(c,s,e)
pdf(paste0("./Figure_S2.pdf"), height = 8, width = 10)
plotTracks(trackList = c(ideoTrack, gatrack, knownGenes, bam1, bam2, hic.track, dtrack.atac, dtrack.foxp3, p.track, e.track, snp.track, SE1, SE2, SE3, SE4, SE5, SE6, chromHMM_Treg, chromHMM_Th17, chromHMM_priT, chromHMM_priB, chromHMM_priHSC),sizes = c(1,1,3,0.7,0.7,1.5,1,1,1,1.2,0.7,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),chromosome = c[i], from = s[i], to = e[i], cex.title = 0.65, rotation.title = 0, showAxis = FALSE,
           background.title = "white", title.width = 2.3,
           fontcolor.title = "black", col.border.title = "white")
dev.off()
