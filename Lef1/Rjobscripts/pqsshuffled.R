library(pqsfinder)
library(stringi)
library(Biostrings)
dna <- readDNAStringSet("../lef1_293T_chipSeq.fa", format = "fasta")
print("No. Lef1 sequences:")
length(dna)
set.seed(15)
shuffledDNA <- stri_rand_shuffle(dna)
shuffledDNA <- DNAStringSet(shuffledDNA)
print("No. shuffled sequences:")
length(shuffledDNA)

shuffledWinners <- vector()
numShuffledWinners <- 0

for(ii in 1:length(shuffledDNA)) {
  pqs <- pqsfinder(shuffledDNA[[ii]], min_score = 25)
  dss <- as(pqs, "DNAStringSet")
  if (length(pqs@elementMetadata@listData[["score"]]) > 0) {
    numShuffledWinners <- numShuffledWinners +1
    shuffledWinners[numShuffledWinners] <- ii
  }
}
print("Number of random GQ's: " )
print(numShuffledWinners)