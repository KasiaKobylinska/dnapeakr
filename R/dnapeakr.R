# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# argument to znf_peaks w postaci tabelki (wczytane uprzednio przez read.table)
mapping <- function(znf_peaks, input_label, epsilon=5000) {
  result <- c()
  if (length(znf_peaks) == 0) {
    return (result)
  }
  ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  for (i in 1:dim(znf_peaks)[1]){
    # chr id: as numeric
    chromosome_id<-ifelse(
      stringi::stri_sub(znf_peaks[i, 1], 4)=="X",
      "X",
      as.numeric(stringi::stri_sub(znf_peaks[i, 1], 4)))
    st <- znf_peaks[i, 2] - epsilon
    en <- znf_peaks[i, 3] + epsilon
    gens <- biomaRt::getBM(attributes = c('ensembl_gene_id','chromosome_name','start_position', 'end_position'),
                 filters = c('chromosome_name','start','end'),
                 values = list(chromosome_id,st,en),
                 mart = ensembl)
    result <- rbind(result,
                 if (nrow(gens) == 0) {
                   data.frame(ensembl_gene_id=NA, chromosome_name=chromosome_id, start_position=st, end_position=en)
                   } else gens[1, ])
  }
  result$input <- input_label
  result <- cbind(znf_peaks, result)
  result <- result[,-8] # cutting out chromosome_name
  # and setting first 6 names of columns:
  colnames(result) <- c("chromosome", "start_peak", "end_peak", "MACS_peak", "peak",
                       "strand", "ensemble_gene_id", "start_position", "end_position", "input" )
  result$diff <- result$start_peak - result$start_position
  return (result)
}

read_human_genomic_file <- function(path) {
  COL_TYPES <- "idddciicccccici"
  COL_NAMES <- c("sw_score", "perc_div", "perc_del", "perc_ins", "query_sequence",
                 "pos_query_begin", "pos_query_end", "pos_query_left",
                 "strand", "matching_repeat", "repeat_class_family",
                 "pos_repeat_begin", "pos_repeat_end", "pos_repeat_left", "id")
  result <- readr::read_table2(path, skip=3, col_names=COL_NAMES, col_types=COL_TYPES)
  return (result)
}

append_genome_overlaps <- function(znf_peaks, genome) {
  grZNF <- GenomicRanges::GRanges(seqnames = znf_peaks$chromosome,
                   ranges = IRanges::IRanges(znf_peaks$start_peak, end = znf_peaks$end_peak),
                   strand = znf_peaks$strand
  )
  grGenome <- GenomicRanges::GRanges(seqnames = genome$query_sequence,
                    ranges = IRanges::IRanges(genome$pos_query_begin, end = genome$pos_query_end)
  )
  overlaps <- GenomicRanges::findOverlaps(grZNF, grGenome, select="first")
  # creating a copy of input data frame, cbind with 1 arg copies the data.frame
  appended <- cbind(znf_peaks)
  for (i in 1:length(overlaps)){
    number<-overlaps[i]
    appended$Start_seq[i] <- genome[number,6]
    appended$End_seq[i] <- genome[number,7]
    appended$matching_repeat[i] <- genome[number,10]
    appended$repeat_class_family[i] <- genome[number,11]
  }
  for (i in 1: dim(appended)[1]){
    appended[i,16] <- (min(as.numeric(appended[i,13]),appended[i,3])-max(as.numeric(appended[i,12]),appended[i,2]))
  }
  colnames(appended)[16]<-"length_of_the_overlap"
  appended[,12]<-as.numeric(appended[,12])
  appended[,13]<-as.numeric(appended[,13])
  appended[,14]<-as.character(appended[,14])
  appended[,15]<-as.character(appended[,15])
  return (appended)
}



