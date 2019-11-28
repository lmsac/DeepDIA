aa.residues.mass.monoisotopic = c(
  A =	71.037114,
  R =	156.101111,
  N =	114.042927,
  D =	115.026943,
  C =	103.009185,
  E =	129.042593,
  Q =	128.058578,
  G =	57.021464,
  H =	137.058912,
  I =	113.084064,
  L =	113.084064,
  K =	128.094963,
  M =	131.040485,
  F =	147.068414,
  P =	97.052764,
  S =	87.032028,
  T =	101.047679,
  # U =	150.95363,
  W =	186.079313,
  Y =	163.06332,
  V =	99.068414
)

aa.residues.mass.average = c(
  A =	71.0779,
  R =	156.1857,
  N =	114.1026,
  D =	115.0874,
  C =	103.1429,
  E =	129.114,
  Q =	128.1292,
  G =	57.0513,
  H =	137.1393,
  I =	113.1576,
  L =	113.1576,
  K =	128.1723,
  M =	131.1961,
  F =	147.1739,
  P =	97.1152,
  S =	87.0773,
  T =	101.1039,
  # U =	150.0379,
  W =	186.2099,
  Y =	163.1733,
  V =	99.1311
)

n.terminus.mass.monoisotopic = c(
  H = 1.007825
)

n.terminus.mass.average = c(
  H = 1.008
)

c.terminus.mass.monoisotopic = c(
  OH = 17.00274
)

c.terminus.mass.average = c(
  OH = 17.007
)

neutral.loss.mass.monoisotopic = c(
  H2O = 18.010565,
  NH3 = 17.026548,
  H3PO4 = 97.977
)

neutral.loss.mass.average = c(
  H2O = 18.015,
  NH3 = 17.031,
  H3PO4 = 97.995
)

modification.delta.monoisotopic = c(
  Acetyl = 42.010565, # 1
  Carbamidomethyl = 57.021464, # 4
  Carbamyl = 43.005814, # 5
  Carboxymethyl = 58.005479, # 6
  Deamidated = 0.984016, # 7
  `ICAT-D:2H(8)` = 450.275205, # 12
  `ICAT-D` = 442.224991, # 13
  `PEO-Iodoacetyl-LC-Biotin` = 414.193691, # 20
  Phospho = 79.966331, # 21
  Propionamide = 71.037114, # 24
  `Pyro-carbamidomethyl` = 39.994915, # 26
  `Glu->pyro-Glu` = -18.010565, # 27
  `Gln->pyro-Glu` = -17.026549, # 28
  Oxidation = 15.994915, #35
  Dimethyl = 28.031300, # 36
  Methylthio = 45.987721, # 39
  `Propionamide:2H(3)` = 74.055944, # 97
  `Propionamide:13C(3)` = 74.047178,
  `ICAT-C` = 227.126991, # 105
  `ICAT-C:13C(9)` = 236.157185, # 206
  iTRAQ4plex = 144.102063, # 214
  TMT6plex = 229.162932 # 737
)

modification.delta.average = c(
  Acetyl = 42.0367, # 1
  Carbamidomethyl = 57.0513, # 4
  Carbamyl = 43.0247, # 5
  Carboxymethyl = 58.0361, # 6
  Deamidated = 0.9848, # 7
  `ICAT-D:2H(8)` = 450.6221, # 12
  `ICAT-D` = 442.5728, # 13
  `PEO-Iodoacetyl-LC-Biotin` = 414.5196, # 20
  Phospho = 79.9799, # 21
  Propionamide = 71.0779, # 24
  `Pyro-carbamidomethyl` = 40.0208, # 26
  `Glu->pyro-Glu` = -18.0153, # 27
  `Gln->pyro-Glu` = -17.0305, # 28
  Oxidation = 15.9994, # 35
  Dimethyl = 28.0532, # 36
  Methylthio = 46.0916, # 39
  `Propionamide:2H(3)` = 74.0964, # 97
  `Propionamide:13C(3)` = 74.0558,
  `ICAT-C` = 227.2603, # 105
  `ICAT-C:13C(9)` = 236.1942, # 206
  iTRAQ4plex = 144.1544, # 214
  TMT6plex = 229.2634 # 737
)

peptide.mass = function(sequence, monoisotopic = FALSE, n.terminus = 'H', c.terminus = 'OH') {
  aa.residues.mass = if(monoisotopic) 
    aa.residues.mass.monoisotopic
  else
    aa.residues.mass.average
  n.terminus.mass = if(monoisotopic) 
    n.terminus.mass.monoisotopic
  else
    n.terminus.mass.average
  c.terminus.mass = if(monoisotopic) 
    c.terminus.mass.monoisotopic
  else
    c.terminus.mass.average
  
  n.mass = if(is.null(n.terminus))
    0
  else if(!(n.terminus %in% names(n.terminus.mass)))
    stop(paste('undefined n-terminus: ', n.terminus))
  else
    n.terminus.mass[n.terminus]
  
  c.mass = if(is.null(c.terminus))
    0
  else if(!(c.terminus %in% names(c.terminus.mass)))
    stop(paste('undefined c-terminus: ', c.terminus))
  else
    c.terminus.mass[c.terminus]
  
  result = sapply(strsplit(sequence, ''), function(aa.seq) {
    undefined.aa = aa.seq[!(aa.seq %in% names(aa.residues.mass))]
    if (length(undefined.aa) > 0) {
      stop(do.call(paste, as.list(c('undefined residue: ', undefined.aa))))
    }
    sum(aa.residues.mass[aa.seq]) + n.mass + c.mass
  })
  names(result) = NULL
  result
}

fragment.mass.neutral = function(
  sequence,
  type, 
  monoisotopic = FALSE,
  modification = NULL,
  n.terminus = 'H',
  c.terminus = 'OH'
) {
  
  fragment.mass.funcs = local({
    cho.mass.monoisotopic = 29.00274
    cho.mass.average = 29.018
    
    if (monoisotopic) {
      n.terminus.mass = n.terminus.mass.monoisotopic
      cho.mass = cho.mass.monoisotopic
      neutral.loss.mass = neutral.loss.mass.monoisotopic
      c.terminus.mass = c.terminus.mass.monoisotopic
    }
    else {
      n.terminus.mass = n.terminus.mass.average
      cho.mass = cho.mass.average
      neutral.loss.mass = neutral.loss.mass.average
      c.terminus.mass = c.terminus.mass.average
    }
    
    m = function(sequence)
      peptide.mass(
        sequence,
        monoisotopic = monoisotopic,
        n.terminus = NULL,
        c.terminus = NULL
      )
    
    n = function(n.terminus)
      if(is.null(n.terminus))
        0
    else if(!(n.terminus %in% names(n.terminus.mass)))
      stop(paste('undefined n-terminus: ', n.terminus))
    else
      n.terminus.mass[n.terminus]
    
    c.. = function(c.terminus)
      if(is.null(c.terminus))
        0
    else if(!(c.terminus %in% names(c.terminus.mass)))
      stop(paste('undefined c-terminus: ', c.terminus))
    else
      c.terminus.mass[c.terminus]
    
    cho = function() cho.mass
    nh3 = function() neutral.loss.mass['NH3']
    h2o = function() neutral.loss.mass['H2O']
    h = function() n.terminus.mass['H']
    nh2 = function() nh3() - h()
    co = function() cho() - h()
    
    a = function(sequence, n.terminus, c.terminus = NULL) 
      n(n.terminus) + m(sequence) - cho()
    a.nh3 = function(sequence, n.terminus, c.terminus = NULL)
      a(sequence, n.terminus) - nh3()
    a.h2o = function(sequence, n.terminus, c.terminus = NULL)
      a(sequence, n.terminus) - h2o()
    b = function(sequence, n.terminus, c.terminus = NULL)
      n(n.terminus) + m(sequence) - h()
    b.nh3 = function(sequence, n.terminus, c.terminus = NULL)
      b(sequence, n.terminus) - nh3()
    b.h2o = function(sequence, n.terminus, c.terminus = NULL)
      b(sequence, n.terminus) - h2o()
    c. = function(sequence, n.terminus, c.terminus = NULL)
      n(n.terminus) + m(sequence) + nh2()
    x = function(sequence, n.terminus = NULL, c.terminus)
      c..(c.terminus) + m(sequence) + co() - h()
    y = function(sequence, n.terminus = NULL, c.terminus)
      c..(c.terminus) + m(sequence) + h()
    y.nh3 = function(sequence, n.terminus = NULL, c.terminus)
      y(sequence, NULL, c.terminus) - nh3()
    y.h2o = function(sequence, n.terminus = NULL, c.terminus)
      y(sequence, NULL, c.terminus) - h2o()
    z = function(sequence, n.terminus = NULL, c.terminus)
      c..(c.terminus) + m(sequence) - nh2()
    
    fragment.mass.funcs = list(
      a = a,
      `a*` = a.nh3,
      ao = a.h2o,
      b = b,
      `b*` = b.nh3,
      bo = b.h2o,
      c = c.,
      x = x,
      y = y,
      `y*` = y.nh3,
      yo = y.h2o,
      z = z
    )
  })
  
  mod = local({
    modification.delta = if(monoisotopic)
      modification.delta.monoisotopic
    else
      modification.delta.average
    
    function(sequence, modification) {
      if(!is.null(modification) && length(modification$pos) > 0) {
        undefined.mod = modification$name[!(modification$name %in% names(modification.delta))]
        if (length(undefined.mod) > 0) {
          stop(do.call(paste, as.list(c('undefined modification: ', undefined.mod))))
        }
        
        result = sum(sapply(1:length(modification$pos), function(i) {
          if(modification$pos[i] <= nchar(sequence) && modification$pos[i] >= 0)
            modification.delta[modification$name[i]]
          else
            0
        }))
        names(result) = NULL
        result
      }
      else
        0
    }
  })
  
  local({
    if(!(type %in% names(fragment.mass.funcs)))
      stop(paste('undefined ion type: ', type))
    result = fragment.mass.funcs[[type]](sequence, n.terminus, c.terminus)
    result = result + mod(sequence, modification)
    names(result) = NULL
    result
  })
}

fragment.ions.mz = function(sequence, type = c('b', 'y'), charge = c(1, 2), monoisotopic = TRUE, modification = NULL, ...) {
  n.frag.seq = substring(sequence, 1, 1:(nchar(sequence) - 1))
  c.frag.seq = rev(substring(sequence, 2:(nchar(sequence)), nchar(sequence)))
  
  unlist(lapply(type, function(t) {
    frag.mass = if (t %in% c('a', 'a*', 'ao', 'b', 'b*', 'bo', 'c')) {
      frag.seq = n.frag.seq
      sapply(frag.seq, function(frag) {
        mass = fragment.mass.neutral(frag, type = t, monoisotopic = monoisotopic, modification, ...)
      })
    }
    else {
      frag.seq = c.frag.seq
      sapply(frag.seq, function(frag) {
        mods = if(is.null(modification))
          modification
        else {
          pos = modification$pos - nchar(sequence) + nchar(frag)
          name = modification$name[pos > 0]
          pos = pos[pos > 0]
          list(pos = pos, name = name)
        }
        mass = fragment.mass.neutral(frag, type = t, monoisotopic = monoisotopic, mods, ...)
      })
    }
    
    result = lapply(charge, function(ch) {
      mz = (frag.mass + ch) / ch
      names(mz) = paste0(t, 1:(nchar(sequence) - 1), '^', ch)
      mz
    })
  }))
}
