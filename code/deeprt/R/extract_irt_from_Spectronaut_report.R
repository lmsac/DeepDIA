library(readr)
report.file = list.files(pattern = '\\.FragmentReport\\.csv$')[1]
psmReport = read_csv(report.file)

irt = local({
  rank = order(psmReport$PSM.Qvalue, decreasing = FALSE)
  dup = duplicated(psmReport$PEP.StrippedSequence[rank])
  uni = rank[!dup]
  indexes = uni[order(uni)]
  
  irt = data.frame(
    sequence = psmReport$PEP.StrippedSequence[indexes],
    irt = psmReport$PP.iRTEmpirical[indexes],
    rt = psmReport$PP.EmpiricalRT[indexes],
    stringsAsFactors = FALSE
  )
})

write.csv(
  irt, 
  paste0(sub('\\.FragmentReport\\.csv', '', report.file), '.irt.csv'),
  row.names = FALSE
)
