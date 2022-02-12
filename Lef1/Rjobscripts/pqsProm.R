library(pqsfinder)

LefG4_prom_seq <- readDNAStringSet("Lef1GQ_H1_prom_int.fa", format = "fasta")
lefG4_prom_overlap_subset_bed <- read.table("lefG4_prom_overlap.bed", sep = "\t")
colnames(lefG4_prom_overlap_subset_bed) <- c("peakChr", "peakStart", "peakEnd", "peakLength",  "peakStrand", "name","score", "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", "geneId", "geneStrand", "transcriptId", "distanceToTSS")


#pqsFinder on promoter regions
winnersI <- vector()
numWinnersI <- 0
numWinnerSites <- 0
winnersIscore <-data.frame(score= integer(), OG_file_row=numeric())

for(ii in 1:length(LefG4_prom_seq)) {
  pqs <- pqsfinder(LefG4_prom_seq[[ii]], min_score = 52)
  dss <- as(pqs, "DNAStringSet")
  if (length(pqs@elementMetadata@listData[["score"]]) > 0) {
    numWinnerSites <- numWinnerSites +1
    numWinnersI <- numWinnersI + length(pqs@elementMetadata@listData[["score"]])
    winnersI[numWinnerSites] <- ii
    hit <- data.frame(pqs@elementMetadata@listData[["score"]], rep(ii, length(pqs@elementMetadata@listData[["score"]])))
    winnersIscore <-  rbind(winnersIscore, hit)
  }
  writeXStringSet(dss, file = "pred_GQ_prom.txt", format = "fasta", append = T)
}

colnames(winnersIscore) <- c("score", "og_df_row")
winnersIscore <- cbind(winnersIscore, seq(1:length(winnersIscore$score)))
colnames(winnersIscore) <- c("score", "og_df_row", "winnerNo")
winnersIscore <- winnersIscore[order(winnersIscore$winnerNo),]


prom_pqs_lef1G4_hits <- data.frame(matrix(ncol=ncol(lefG4_prom_overlap_subset_bed), nrow=length(winnersIscore)))
for(ii in 1:length(winnersIscore)) {
  prom_pqs_lef1G4_hits[ii,] <- lefG4_prom_overlap_subset_bed[winnersIscore[ii,2],]
}
colnames(prom_pqs_lef1G4_hits) <- c("peakChr", "peakStart", "peakEnd", "peakLength",  "peakStrand", "name","score", "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", "geneId", "geneStrand", "transcriptId", "distanceToTSS")
write.csv(prom_pqs_lef1G4_hits, "G4prom_lef1.txt", row.names = F, quote = F)


# take only unique ids
gene.ids <- unique(prom_pqs_lef1G4_hits$geneId)

# write names to an output file
write.table(gene.ids, file="lef1G4_regulated_genes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.csv(winnersIscore, "lef1GQ_promoters.csv", row.names = F)