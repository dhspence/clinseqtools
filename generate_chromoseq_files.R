#!/usr/bin/R

suppressMessages(require(biomaRt))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))

args <- commandArgs(trailingOnly = T)

if (length(args) < 1){
  cat("generate_chromoseq_files <gene list>\nProvide a list of genes in a tab-delimited file:\n\tgene1\tgene2\t(for SVs)\n\tgene\t\t(for gene-level mutation)\n")
  quit()
}

infile <- args[1]

gene_clusters <- data.frame(
  chr = factor(c("7","7","14","14","14","22"),levels=c(seq(1,22,by=1),"X","Y","MT")),
  start = c(38240023,142299010,21621903,22422545,105586436,22026075),
  end = c(38368055,142813287,22552132,22466577,106879844,22922913),
  gene_name = c("TRG","TRB","TRA","TRD","IGH","IGL"),
  info = c("CS_svgenes|TRG|TRG|transcript|38240023|38368055|-1","CS_svgenes|TRB|TRB|transcript|142299010|142813287|1","CS_svgenes|TRA|TRA|transcript|21621903|22552132|1","CS_svgenes|TRD|TRD|transcript|22422545|22466577|1","CS_svgenes|IGH|IGH|transcript|105586436|106879844|-1","CS_svgenes|IGL|IGL|transcript|22026075|22922913|1")
)

canonical_chromosomes <- c(1:22, "X", "Y", "MT")  # Note: It's MT in Ensembl for mitochondrial DNA, not M

# Read and parse the text file
data <- read.table(infile, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene1","gene2"),fill=T)

# get rows with one gene, which are for small variant analysis
genelist <- data %>% filter(is.na(gene2) | gene2=='') %>% pull(gene1)

# get rows with 2 genes for SV analysis
svgenepairs <- data %>% filter(gene1!='' & gene2!='') %>% rowwise() %>% mutate(id=str_c(min(gene1,gene2),max(gene1,gene2),sep="-")) %>% dplyr::select(id,gene1,gene2) %>% distinct(id,.keep_all = T)

# Create gene list
svgenelist <- unique(c(svgenepairs$gene1, svgenepairs$gene2))

print(str_c("getting info for ",length(svgenelist)," SV genes."))

# Use biomaRt to query ensembl
ensembl <- useEnsembl('genes',dataset='hsapiens_gene_ensembl',version = 105)
genes <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 
                              'external_gene_name', 'gene_biotype', 'ensembl_gene_id', 
                              'ensembl_transcript_id', 'transcript_biotype',
                              'strand','genomic_coding_start','genomic_coding_end'),
               filters = 'external_gene_name',
               values = svgenelist,
               mart = ensembl)

notfound <- setdiff(svgenelist,genes$external_gene_name)
if (length(notfound) > 0){
  print(str_c(length(notfound)," SV genes not found:"))
  print(notfound)
}

genes$transcript_length <- genes$end_position-genes$start_position
genes$is_coding <- ifelse(genes$transcript_biotype=="protein_coding",1,0)

# Filter by protein coding and then by longest transcript length
trx <- genes[grep("_",genes$chromosome_name,invert=T), ]
trx <- trx %>% mutate(genomic_coding_start=ifelse(is.na(genomic_coding_start),start_position,genomic_coding_start),genomic_coding_end=ifelse(is.na(genomic_coding_end),end_position,genomic_coding_end))
trx <- merge(trx,trx %>% group_by(ensembl_transcript_id) %>% summarise(cds_start=min(genomic_coding_start,na.rm=T),cds_end=max(genomic_coding_end,na.rm=T)))

longest_transcripts <- trx[with(trx,order(-is_coding,-transcript_length)), ]
longest_transcripts <- longest_transcripts[!duplicated(longest_transcripts$ensembl_gene_id), ]
longest_transcripts <- longest_transcripts[grep("_",longest_transcripts$chromosome_name,invert=T),]
longest_transcripts <- longest_transcripts[longest_transcripts$chromosome_name %in% canonical_chromosomes,]

# Create output
output <- data.frame(chr = factor(longest_transcripts$chromosome_name,levels=c(seq(1,22,by=1),"X","Y","MT")),
                     start = longest_transcripts$start_position-1,
                     end = longest_transcripts$end_position,
                     gene_name = longest_transcripts$external_gene_name,
                     info = paste("CS_svgenes",
                                  longest_transcripts$ensembl_gene_id,
                                  longest_transcripts$ensembl_transcript_id,
                                  'transcript',
                                  longest_transcripts$cds_start,
                                  longest_transcripts$cds_end,
                                  longest_transcripts$strand,
                                  sep = "|"))

if (any(svgenelist %in% gene_clusters$gene_name)){
  x <- gene_clusters[which(gene_clusters$gene_name %in% svgenelist),]
  output <- rbind(output,x) %>% arrange(chr,start)
}

output <- output %>% arrange(chr,start) %>% mutate(chr=str_c('chr',chr)) %>% mutate(chr=ifelse(chr=="chrMT","chrM",chr))

# Write bed file
write.table(output, "chromoseq_sv_genes.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# write sv gene pair file
gene_strand <- output %>% mutate(strand=ifelse(str_detect(info,'-1'),'-','+')) %>% dplyr::select(gene_name,strand)
output2 <- merge(merge(svgenepairs,gene_strand,by.x='gene1','gene_name',all.x=T),gene_strand,by.x='gene2',by.y='gene_name',all.x=T) %>% dplyr::rename(strand1=strand.x,strand2=strand.y) %>% dplyr::select(gene1,strand1,gene2,strand2)

write.table(output2, "chromoseq_recurrent_sv.txt", sep = ",", row.names = FALSE, quote = FALSE, col.names = FALSE)

print("getting genes for gene-level analysis.")

genes <- getBM(attributes = c('chromosome_name','external_gene_name','ensembl_exon_id', 'gene_biotype', 'ensembl_gene_id', 
                              'ensembl_transcript_id', 'transcript_biotype',
                              'strand','genomic_coding_start','genomic_coding_end'),
               filters = c('external_gene_name','transcript_is_canonical','transcript_biotype'),
               values = list(genelist,TRUE,'protein_coding'),
               mart = ensembl)

notfound <- setdiff(genelist,genes$external_gene_name)
if (length(notfound) > 0){
  print(str_c(length(notfound)," genes not found:"))
  print(notfound)
}

trx <- genes[grep("_",genes$chromosome_name,invert=T), ] %>% filter(!is.na(genomic_coding_start))
trx <- merge(trx,trx %>% group_by(ensembl_transcript_id) %>% summarise(cds_start=min(genomic_coding_start,na.rm=T),cds_end=max(genomic_coding_end,na.rm=T)))
trx <- trx %>% arrange(ensembl_transcript_id,genomic_coding_start) %>% group_by(ensembl_transcript_id) %>% mutate(exon_number = ifelse(strand == 1, row_number(), n() - row_number() + 1)) %>% ungroup() %>% 
  mutate(exon_number = str_c('exon_',exon_number)) %>% arrange(chromosome_name,genomic_coding_start) %>% filter(!is.na(genomic_coding_start))

output <- data.frame(chr = factor(trx$chromosome_name,levels=c(seq(1,22,by=1),"X","Y","MT")),
                     start = trx$genomic_coding_start-3,
                     end = trx$genomic_coding_end+2,
                     gene_name = trx$external_gene_name,
                     info = paste("CS_genes",
                                  trx$ensembl_gene_id,
                                  trx$ensembl_transcript_id,
                                  trx$exon_number,
                                  trx$cds_start,
                                  trx$cds_end,
                                  trx$strand,
                                  sep = "|"))

output <- output %>% arrange(chr,start) %>% arrange(chr,start) %>% mutate(chr=str_c('chr',chr)) %>% mutate(chr=ifelse(chr=="chrMT","chrM",chr))
# Write bed file
write.table(output, "chromoseq_genes.bed", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
