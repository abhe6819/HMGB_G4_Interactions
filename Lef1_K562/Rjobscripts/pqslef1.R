library(pqsfinder)
library(stringi)
library(Biostrings)
dna <- readDNAStringSet("../lef1_293T_chipSeq.fa", format = "fasta")
print("No. Lef1 sequences:")
length(dna)
winners <- vector()
numWinners <- 0

for(ii in 1:length(dna)) {
  pqs <- pqsfinder(dna[[ii]], min_score = 25)
  dss <- as(pqs, "DNAStringSet")
  if (length(pqs@elementMetadata@listData[["score"]]) > 0) {
    numWinners <- numWinners +1
    winners[numWinners] <- ii
  }
}
print("Number of predicted GQ's: " )
print(numWinners)
