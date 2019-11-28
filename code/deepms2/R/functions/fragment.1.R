fragment.ions.mz.1 = function(sequence, type = c('b', 'y'), charge = c(1, 2), 
                              modification = NULL, loss = c('NH3', 'H2O', 'H3PO4'),
                              carbamidomethyl = FALSE, monoisotopic = TRUE, ...) {
  if (nchar(sequence) <= 1) {
    stop('invalid sequence')
  }
  sequence = toupper(sequence)
  
  if (!is.null(modification) && !is.na(modification) && nchar(modification) > 0) {
    s = strsplit(modification, ';')[[1]]
    name = gsub('[()]', '', stringr::str_extract(s, '\\(.+\\)$'))
    pos = as.integer(gsub('[A-Z(]', '', stringr::str_extract(s, '^[A-Z]+[0-9]+\\(')))
    modification = list(pos = pos, name = name)
  }
  else {
    modification = NULL
  }
  if (carbamidomethyl) {
    cysteine = stringr::str_locate_all(sequence, 'C')[[1]][, 1]
    modification$pos = c(modification$pos, cysteine)
    modification$name = c(modification$name, rep('Carbamidomethyl', length(cysteine)))
  }
  
  new.type = c()
  if ('NH3' %in% loss) {
    new.type = c(new.type, paste0(type, '*'))
  }
  if ('H2O' %in% loss) {
    new.type = c(new.type, paste0(type, 'o'))
  }
  type = c(type, new.type)
  
  fragments = fragment.ions.mz(sequence, type = type, charge = charge, modification = modification, monoisotopic = monoisotopic, ...)
  
  if ('H3PO4' %in% loss) {
    h3po4 = (if (monoisotopic) neutral.loss.mass.monoisotopic else neutral.loss.mass.average)['H3PO4']
    
    phospho = modification$pos[grep('Phospho', modification$name)]
    if (length(phospho) > 0) {
      p.fragments = c()
      if ('b' %in% type) {
        first = min(phospho)
        if (first < nchar(sequence)) {
          p.fragments = c(p.fragments, unlist(lapply(charge, function(ch) {
            bp = fragments[paste0('b', first:(nchar(sequence) - 1), '^', ch)] - h3po4 / ch
            names(bp) = paste0('bp', first:(nchar(sequence) - 1), '^', ch)
            bp
          })))
        }
      }
      if ('y' %in% type) {
        last = max(phospho)
        if (last > 1) {
          p.fragments = c(p.fragments, unlist(lapply(charge, function(ch) {
            yp = fragments[paste0('y', (nchar(sequence) - last + 1):(nchar(sequence) - 1), '^', ch)] - h3po4 / ch
            names(yp) = paste0('yp', (nchar(sequence) - last + 1):(nchar(sequence) - 1), '^', ch)
            yp
          })))
        }
      }
      
      fragments = c(fragments, p.fragments)
    }
  }
  fragments
}

peptide.mass.1 = function(sequence, charge = NA, modification = NULL,
                          carbamidomethyl = FALSE, monoisotopic = TRUE, ...) {
  sequence = toupper(sequence)
  mass = peptide.mass(sequence = sequence, monoisotopic = monoisotopic, ...)
  
  if (carbamidomethyl) {
    cysteine = which(strsplit(sequence, '')[[1]] == 'C')
    sapply(cysteine, function(i) {
      label = paste0('C', i, '(Carbamidomethyl)')
      if (is.null(modification) || is.na(modification) || nchar(modification) == 0) {
        modification <<- label
      }
      else if (!grepl(label, modification, fixed = TRUE)) {
        modification <<- paste0(modification, ';', label)
      }
    })
  }
  
  if (!is.null(modification) && !is.na(modification) && nchar(modification) > 0) {
    s = strsplit(modification, ';')[[1]]
    name = gsub('[()]', '', stringr::str_extract(s, '\\(.+\\)$'))
    mass = mass + sum(sapply(name, function(m) (if (monoisotopic) modification.delta.monoisotopic else modification.delta.average)[m]))
  }
  
  if (!is.na(charge) && !is.null(charge) && charge != 0)
    mass = (mass + (if (monoisotopic) n.terminus.mass.monoisotopic['H'] else n.terminus.mass.average['H']) * charge) / charge
  as.numeric(mass)
}
