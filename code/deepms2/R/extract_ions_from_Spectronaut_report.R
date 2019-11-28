library(readr)
report.file = list.files(pattern = '\\.FragmentReport\\.csv$')[1]
fragmentReport = read_csv(report.file)

psm.indexes = which(sapply(1:nrow(fragmentReport), function(i) {
  if (i == 1) {
    return(TRUE)
  }
  
  if (fragmentReport$PSM.MS2ScanNumber[i] != fragmentReport$PSM.MS2ScanNumber[i - 1] ||
    fragmentReport$PEP.StrippedSequence[i] != fragmentReport$PEP.StrippedSequence[i - 1] ||
    fragmentReport$PP.Charge[i] != fragmentReport$PP.Charge[i - 1] ||
    fragmentReport$PG.ProteinAccessions[i] != fragmentReport$PG.ProteinAccessions[i - 1] ||
    fragmentReport$R.FileName[i] != fragmentReport$R.FileName[i - 1]) {
    TRUE
  }
  else {
    FALSE
  }
}))

lapply(2:3, function(charge) {
  extracted.ions = lapply(1:length(psm.indexes), function(i) {
    if (fragmentReport$PP.Charge[psm.indexes[i]] != charge) {
      return(NULL)
    }
    
    if (i < length(psm.indexes)) {
      rowIndexes = psm.indexes[i]:psm.indexes[i + 1]
    }
    else {
      rowIndexes = psm.indexes[i]:nrow(fragmentReport)
    }
    
    sequence = fragmentReport$PEP.StrippedSequence[rowIndexes[1]]
    charge = fragmentReport$PP.Charge[rowIndexes[1]]
    
    ions = list(
      b1 = rep(0, nchar(sequence) - 1),
      bn1 = rep(0, nchar(sequence) - 1),
      bo1 = rep(0, nchar(sequence) - 1),
      b2 = rep(0, nchar(sequence) - 1),
      bn2 = rep(0, nchar(sequence) - 1),
      bo2 = rep(0, nchar(sequence) - 1),
      y1 = rep(0, nchar(sequence) - 1),
      yn1 = rep(0, nchar(sequence) - 1),
      yo1 = rep(0, nchar(sequence) - 1),
      y2 = rep(0, nchar(sequence) - 1),
      yn2 = rep(0, nchar(sequence) - 1),
      yo2 = rep(0, nchar(sequence) - 1)
    )
    
    lapply(rowIndexes, function(row) {
      type = paste0(
        fragmentReport$FI.FrgType[row],
        ifelse(
          fragmentReport$FI.LossType[row] == 'noloss', '',
          ifelse(
            fragmentReport$FI.LossType[row] == 'NH3', 'n',
            ifelse(
              fragmentReport$FI.LossType[row] == 'H2O', 'o', 
              ''
            )
          )
        ),
        fragmentReport$FI.Charge[row]
      )
      num = as.integer(fragmentReport$FI.FrgNum[row])
      intensity = as.numeric(fragmentReport$FI.Intensity[row])
      if (type %in% names(ions) && num > 0 && num <= nchar(sequence) - 1 && !is.na(intensity)) {
        ions[[type]][num] <<- intensity
      }
    })
    
    qvalue = fragmentReport$PSM.Qvalue[rowIndexes[1]]
    
    list(
      peptide = sequence,
      charge = charge,
      ions = ions,
      qvalue = qvalue
    )
  })
  
  extracted.ions = extracted.ions[!sapply(extracted.ions, is.null)]
  save(extracted.ions, file = paste0(sub('\\.FragmentReport\\.csv', '', report.file), '_charge', charge, '.ions.RData'))
})

