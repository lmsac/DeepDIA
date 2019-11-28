library(readr)

detectability.file = list.files(pattern = '\\.prediction\\.detectability\\.csv$')[1]
detectability = read_csv(detectability.file)

detectabilityThreshold = 0.5

peptides = read_csv(paste0(
  '../',
  sub('\\.prediction\\.detectability\\.csv$', '.peptide.csv', detectability.file)
))

peptides$detectability = detectability$detectability[match(
  paste0(peptides$nTerminal, '.', peptides$sequence, '.', peptides$cTerminal), detectability$sequence
)]

peptides = subset(peptides, detectability >= detectabilityThreshold)

peptides = local({
  rank = order(peptides$detectability, decreasing = TRUE)
  filtered = duplicated(peptides$sequence[rank])
  indexes = rank[!filtered]
  indexes = sort(indexes)
  peptides = peptides[indexes, ]
})


write.csv(
  peptides,
  paste0(
    sub('\\.prediction\\.detectability\\.csv$', '', detectability.file),
    '_detectability',
    sprintf(
      '%02d%s', as.integer(detectabilityThreshold / 0.01), 
      gsub('^0(\\.00)?|0+$', '', sprintf('%f', detectabilityThreshold %% 0.01))
    ),
    '.peptide.csv'
  )
  row.names = FALSE
)

