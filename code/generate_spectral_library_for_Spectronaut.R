create.assays = function(peptides, ions, irt, precursorCharge) {
  assay = lapply(1:length(ions), function(i) {
    sequence = ions[[i]]$peptide
    precursorMz = (peptide.mass.1(sequence, carbamidomethyl = TRUE) + n.terminus.mass.monoisotopic['H'] * precursorCharge) / precursorCharge
    
    fragmentIntensity = c(
      ions[[i]]$ions$b1,
      ions[[i]]$ions$bn1,
      ions[[i]]$ions$bo1,
      ions[[i]]$ions$b2,
      ions[[i]]$ions$bn2,
      ions[[i]]$ions$bo2,
      ions[[i]]$ions$y1,
      ions[[i]]$ions$yn1,
      ions[[i]]$ions$yo1,
      ions[[i]]$ions$y2,
      ions[[i]]$ions$yn2,
      ions[[i]]$ions$yo2
    )
    fragmentType = c(
      rep('b', length(ions[[i]]$ions$b1)),
      rep('b', length(ions[[i]]$ions$bn1)),
      rep('b', length(ions[[i]]$ions$bo1)),
      rep('b', length(ions[[i]]$ions$b2)),
      rep('b', length(ions[[i]]$ions$bn2)),
      rep('b', length(ions[[i]]$ions$bo2)),
      rep('y', length(ions[[i]]$ions$y1)),
      rep('y', length(ions[[i]]$ions$yn1)),
      rep('y', length(ions[[i]]$ions$yo1)),
      rep('y', length(ions[[i]]$ions$y2)),
      rep('y', length(ions[[i]]$ions$yn2)),
      rep('y', length(ions[[i]]$ions$yo2))
    )
    fragmentNumber = c(
      1:length(ions[[i]]$ions$b1),
      1:length(ions[[i]]$ions$bn1),
      1:length(ions[[i]]$ions$bo1),
      1:length(ions[[i]]$ions$b2),
      1:length(ions[[i]]$ions$bn2),
      1:length(ions[[i]]$ions$bo2),
      1:length(ions[[i]]$ions$y1),
      1:length(ions[[i]]$ions$yn1),
      1:length(ions[[i]]$ions$yo1),
      1:length(ions[[i]]$ions$y2),
      1:length(ions[[i]]$ions$yn2),
      1:length(ions[[i]]$ions$yo2)
    )
    fragmentCharge = c(
      rep('1', length(ions[[i]]$ions$b1)),
      rep('1', length(ions[[i]]$ions$bn1)),
      rep('1', length(ions[[i]]$ions$bo1)),
      rep('2', length(ions[[i]]$ions$b2)),
      rep('2', length(ions[[i]]$ions$bn2)),
      rep('2', length(ions[[i]]$ions$bo2)),
      rep('1', length(ions[[i]]$ions$y1)),
      rep('1', length(ions[[i]]$ions$yn1)),
      rep('1', length(ions[[i]]$ions$yo1)),
      rep('2', length(ions[[i]]$ions$y2)),
      rep('2', length(ions[[i]]$ions$yn2)),
      rep('2', length(ions[[i]]$ions$yo2))
    )
    fragmentLossType = c(
      rep('noloss', length(ions[[i]]$ions$b1)),
      rep('NH3', length(ions[[i]]$ions$bn1)),
      rep('H2O', length(ions[[i]]$ions$bo1)),
      rep('noloss', length(ions[[i]]$ions$b2)),
      rep('NH3', length(ions[[i]]$ions$bn2)),
      rep('H2O', length(ions[[i]]$ions$bo2)),
      rep('noloss', length(ions[[i]]$ions$y1)),
      rep('NH3', length(ions[[i]]$ions$yn1)),
      rep('H2O', length(ions[[i]]$ions$yo1)),
      rep('noloss', length(ions[[i]]$ions$y2)),
      rep('NH3', length(ions[[i]]$ions$yn2)),
      rep('H2O', length(ions[[i]]$ions$yo2))
    )
    fragmentMz = fragment.ions.mz.1(
      sequence, 
      type = c('b', 'y'), 
      charge = c(1, 2),
      loss = c('NH3', 'H2O'),
      carbamidomethyl = TRUE
    )[paste0(
      fragmentType,
      c('noloss' = '', 'NH3' = '*', 'H2O' = 'o')[fragmentLossType],
      fragmentNumber,
      '^', fragmentCharge
    )]
    fragment = data.frame(
      fragmentType, fragmentNumber, fragmentCharge, fragmentLossType,
      fragmentMz, fragmentIntensity,
      stringsAsFactors = FALSE
    )
    fragment = fragment[fragment$fragmentIntensity > 0, ]
    if (nrow(fragment) == 0)
      return(NULL)
    
    
    iRT = irt$irt[match(sequence, irt$sequence)]
    
    # proteinId = report$PG.ProteinAccessions[match(sequence, peptides$PEP.StrippedSequence)]
    proteinId = peptides$protein[match(sequence, peptides$sequence)]
    
    data.frame(
      PrecursorMz = precursorMz,
      FragmentMz = fragment$fragmentMz,
      iRT = iRT,
      #iRTSourceSpecific,
      RelativeFragmentIntensity = fragment$fragmentIntensity,
      StrippedSequence = sequence,
      PrecursorCharge = precursorCharge,
      FragmentType = fragment$fragmentType,
      FragmentNumber = fragment$fragmentNumber,
      FragmentCharge = fragment$fragmentCharge,
      ProteinId = proteinId,
      #ProteinName,
      #ProteinDescription,
      ModifiedSequence = paste0('_', gsub('C', 'C[Carbamidomethyl (C)]', sequence), '_'),
      #LabelModifiedSequence,
      FragmentLossType = fragment$fragmentLossType,
      #ExcludeFromAssay,
      #IsProteotypic,
      stringsAsFactors = FALSE
    )
  })
}



file = list.files(pattern = '\\.peptide\\.csv$')[1]

library(readr)
peptides = read_csv(file)
irt = read_csv(sub('\\.peptide\\.csv$', '.prediction.irt.csv', file))

charges = 2:3

assays = lapply(charges, function(charge) {
  ions = rjson::fromJSON(file = sub('\\.peptide\\.csv$', paste0('_charge', charge, '.prediction.ions.json'), file))
  assay = create.assays(peptides, ions, irt, precursorCharge = charge)
})

write.csv(
  do.call(rbind, do.call(c, assays)),
  sub('\\.peptide\\.csv$', '.prediction.library.csv', file), 
  row.names = FALSE
)




