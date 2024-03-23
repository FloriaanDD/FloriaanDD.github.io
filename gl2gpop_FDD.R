
gl2gpop <- function(glthing, filename) {
 
  afile <- file( filename, open = 'w')

  # Header line
  cat( 'description\n', file = afile)
  
  # Locus names:  
  cat(adegenet::locNames( glthing), sep = '\n', file = afile)
  
  # Get "genotypes" as integers
  glint <- as.matrix( glthing) # rows=samples, cols=loci; 0, 1, 2 or NA
  # ... found out by table( c( testgl), useNA='always')
  
  glint[is.na( glint)] <- 3 # for now
  glint <- glint + 1
  glint <- t(glint)  #This is important- fixes the error that the function is making (reading matrix in wrong order)

  outputs <- c( '0303', '0304', '0404', '0000')
  
  each_pop <- adegenet::pop( glthing)
  
  # Per population:
  for (this_pop in adegenet::popNames(glthing)) {
    cat( 'pop\n', file = afile)
    
    # pop( glthing)==this_pop
    
    this_recode <- outputs[ glint[ ,each_pop == this_pop]] # Important- this has also changed- pops now assigned to columns, to match transposed matrix. 
    
    this_recode <- matrix( this_recode, ncol = adegenet::nLoc( glthing), byrow = TRUE)
    
    these_rows <- apply( this_recode, 1, paste, collapse = ' ') # one string per sample in this pop
    
    these_rows <- paste( adegenet::indNames( glthing)[ each_pop == this_pop], these_rows, sep = ' , ')
    
    cat( these_rows, sep = '\n', file = afile)
  }
  
  close( afile)
}