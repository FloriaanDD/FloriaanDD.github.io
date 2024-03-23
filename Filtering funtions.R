#Multiple SNPs on same fragment
gl.filter.secondaries.FDD <- function (gl) {
  x <- gl
  if (class(x) == "genlight") {
    cat("Filtering a genlight object","\n")
  }
  else {
    cat("Fatal Error: Specify either a genlight object","\n")
    stop()
  }
  cat("Total number of SNP loci:", nLoc(x), "\n")
  oldlength <- length(x@other$loc.metrics$AlleleID)
  a <- data.frame(x@other$loc.metrics$CloneID,x@other$loc.metrics$RepAvg,
                  x@other$loc.metrics$AvgPIC, x@other$loc.metrics$AvgCountRef,
                  x@other$loc.metrics$AvgCountSnp, x@other$loc.metrics$AlleleID)
  colnames(a) <- c("CloneID", "RepAvg", "AvgPIC","AvgCountRef", "AvgCountSnp", "AlleleID")
  a <- a[order (a$CloneID, -a$RepAvg, -a$AvgPIC, -a$AvgCountRef, -a$AvgCountSnp),]
  nodups <- a[!duplicated(a$CloneID),]
  nosecondaries <- nodups$AlleleID
  x <- x[, x@other$loc.metrics$AlleleID %in% nosecondaries, treatOther=TRUE] 
  #Subset loc.metrics
  x@other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics$AlleleID %in% nosecondaries ,] 
  newlength <- length(x@other$loc.metrics$AlleleID)
  cat("   Number of secondaries:", oldlength-newlength, "\n")
  cat("   Number of loci after secondaries removed:", newlength, "\n")
  return(x)
}

 #CallRate
gl.filter.missing.data.FDD <- function(gl, loc.lb, loc.hb, ind.lb, ind.hb, iterations) {
  x <- gl
  ind.ite <- seq( from=ind.lb, to=ind.hb, length=iterations)
  loc.ite <- seq( from=loc.lb, to=loc.hb, length=iterations)
 
  for( iloop in 1:iterations){
    x <- gl.filter.callrate( x, method = "ind", threshold=ind.ite[ iloop], recalc = FALSE, v = 0)
    x <- gl.filter.callrate( x, method = "loc", threshold=loc.ite[ iloop], recalc = FALSE, v = 0)
  }
  cat("Summary of filtered dataset\n")
  cat("  SNPs with CallRate >", loc.hb, ":", nLoc(x), "\n")
  cat("  Individuals with CallRate >", ind.hb, ":", nInd(x), "\n")
  cat("  No. of loci removed:", nLoc(gl) - nLoc(x), "\n")
  cat("  No. of individuals removed:", nInd(gl) - nInd(x), "\n")
  if((nInd(gl) - nInd(x))>0){
    cat(paste("  Removed individuals: ", gl$ind.names[!(gl$ind.names %in% x$ind.names)],
            "from populations",  gl$pop[!(gl$ind.names %in% x$ind.names)],
            "\n"))
  }
  return( x) 
}

#Individual heterozygosity
gl.report.het.FDD<-function (glall){
  x <- glall
  n0 <- nInd(x)
  ind.Het <-rowSums(as.matrix(x)==1, na.rm =TRUE)/nLoc(x)
	par(mfrow=c(1,1))
      plot(ind.Het, ylim=c(min(ind.Het)-0.05,max(ind.Het)+0.05), main="ind.Het=(nAB, na.rm=T)/nLoc(x)")
 }
 
gl.filter.het.FDD <- function (gl, LowerT, UpperT)
  {
  x <- gl
   if (class(x) == "genlight") {
     cat("Reporting for a genlight object\n")
  }
  n0 <- nInd(x)
   cat("Initial no. of individuals =", n0, "\n")
  ind.Het <-rowSums(as.matrix(x)==1, na.rm =TRUE)/nLoc(x)
  if (sum(ind.Het >= LowerT) == 0) 
     stop(paste("Minimum individual heterozygosity =", min(ind.Het), 
                ". Nominated threshold of", LowerT, "too stringent.\n No individuals remain.\n"))
   if (sum(ind.Het <= UpperT) == 0) 
    stop(paste("Maximum individual heterozygosity =", max(ind.Het), 
                ". Nominated threshold of", UpperT, "too stringent.\n No individuals remain.\n"))
   x2 <- x[ind.Het >= LowerT & ind.Het <= UpperT, ]
   if (class(x) == "genlight") {
     cat("Filtering a genlight object\n  no. of individuals deleted =", 
        (n0 - nInd(x2)), "\nIndividuals retained =", 
         nInd(x2), "\n")
   }
   if (any(ind.Het <= LowerT)) {
     x3 <- x[ind.Het <= LowerT, ]
    if (length(x3) > 0) {
      cat("\nList of individuals deleted because of low heterozygosity\n", indNames(x3), "\n")
      cat("\n   from populations\n", as.character(pop(x3)), "\n")
    }
  }
  if (any(ind.Het >= UpperT)) {
      x3 <- x[ind.Het >= UpperT, ]
      if (length(x3) > 0) {
      cat("\nList of individuals deleted because of high heterozygosity\n", indNames(x3), "\n")
      cat("\n   from populations\n", as.character(pop(x3)), "\n")
    }  
  }
  cat("\nSummary of filtered dataset\n")
  cat(paste(LowerT, "<  Individuals with heterozygosity <=", UpperT, "\n"))
  cat(paste("  No. of loci:", nLoc(x2), "\n"))
  cat(paste("  No. of individuals:", nInd(x2), "\n"))
  cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))), 
            "\n"))
  par(mfrow=c(1,1))
  plot(ind.Het, ylim=c(min(ind.Het)-0.05,max(ind.Het)+0.05), main="ind.Het=(nAB, na.rm=T)/nLoc(x)")
  abline(h=LowerT)
  abline(h=UpperT)
  return(x2)
}

#Minor allele frequency
gl.filter.maf.FDD <- function (gl, threshold=0.01) {
  x <- gl
  if (class(x) == "genlight") {
    cat("Filtering a genlight object","\n")
  }
  else {
    cat("Fatal Error: Specify either a genlight object","\n")
    stop()
  }
  cat("Total number of SNP loci:", nLoc(gl), "\n")
  #calculate the allele frequencies
  af <- (colMeans(as.matrix(x), na.rm = T)/2)
  #switch the frequencies around if major and minor allele have been genotyped incorrectly
  af <- ifelse(af > 0.5, 1 - af, af)
  maf <- data.frame(x@other$loc.metrics$AlleleID, af)
  maf <- maf[maf[,2] >= threshold,]
  x <- x[, x@other$loc.metrics$AlleleID %in% maf[,1], treatOther=TRUE] #doesn't subset loc.metrics
  x@other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics$AlleleID %in% maf[,1] ,] #subsetting loc.metrics
  cat("   Number of loci with MAF <", threshold, ":",nLoc(gl)-nLoc(x), "\n")
  cat("   Number of loci after filtering:", nLoc(x), "\n")
  cat("   Number of individuals:", nInd(x), "\n")
  cat("   Number of populations: ", nlevels(pop(x)),"\n")
  return(x)
}

#Count per SNP
gl.filter.lowcount.FDD <- function (gl, threshold = 10){
x <- gl
cat("Total number of SNP loci:", nLoc(x), "\n")
        x2 <- x[, x$other$loc.metrics["SumCount"] >= threshold]
        x2$other$loc.metrics <- x$other$loc.metrics[x$other$loc.metrics["SumCount"] >= 
            threshold, ]
    cat("No. of loci deleted =", (nLoc(x) - nLoc(x2)), "\n")
    cat(paste("  Read depth >=", threshold, "\n"))
    cat(paste("  No. of loci:", nLoc(x2), "\n"))
    cat(paste("  No. of individuals:", nInd(x2), "\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))), 
        "\n"))
return (x2)
}

gl.filter.highcount.FDD <- function (gl, threshold = 300){
  x <- gl
  cat("Total number of SNP loci:", nLoc(x), "\n")
  x2 <- x[, x$other$loc.metrics["SumCount"] <= threshold]
  x2$other$loc.metrics <- x$other$loc.metrics[x$other$loc.metrics["SumCount"] <= threshold, ]
  cat("No. of loci deleted =", (nLoc(x) - nLoc(x2)), "\n")
  cat(paste("  Read depth <=", threshold, "\n"))
  cat(paste("  No. of loci:", nLoc(x2), "\n"))
  cat(paste("  No. of individuals:", nInd(x2), "\n"))
  cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))),"\n"))
  return (x2)
}
	
	
	
	


gl.filter.loc.het.FDD <- function (glall, threshold = 0.6){
  x <- glall
  cat("Total number of SNP loci:", nLoc(x), "\n")
  loc.Het <-colSums(as.matrix(x)==1, na.rm =T)/nInd(x)

  het <- data.frame(x@other$loc.metrics$AlleleID, loc.Het)
  het <- het[het[,2] <= threshold,]
  x2 <- x[, x@other$loc.metrics$AlleleID %in% het[,1], treatOther=TRUE] #doesn't subset loc.metrics
  x2$other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics$AlleleID %in% het[,1] ,] #subsetting loc.metrics
  cat("No. of loci deleted =", (nLoc(x) - nLoc(x2)), "\n")
  cat(paste("  SNPs with heterozygosity <=", threshold, "\n"))
  cat(paste("  No. of loci:", nLoc(x2), "\n"))
  cat(paste("  No. of individuals:", nInd(x2), "\n"))
  cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))), 
            "\n"))
  plot(loc.Het)
  abline(h=threshold)
  return (x2)
}


#' Converts genlight objects to STRUCTURE formated files
#'
#' This function exports genlight objects to STRUCTURE formatted files (be aware there is a gl2faststruture version as well). It is based on the code provided by Lindsay Clark (originally for genind objects: [Lindsay happy to link to your github repo if you want]) and this function is basically a wrapper around his genind2structure function.
#' @param gl -- genlight containing lat longs  [required]
#' @param pops -- switch if population column should be added
#' @param outfile -- name (path) of the output shape file
#' @param outpath -- path of the output file. Default is to tempdir(). If to be saved in the current working directory change to "."
#' @param v -- verbosity: if v=0 no output, v=1 reports name and path of output file. default 1
#' @export
#' @author Lindsay V. Clark [your email], wrapper by Bernd Gruber 
#' @examples
#' \dontrun{
#' gl2structure(testset.gl)
#'}

gl2structure <- function(gl, pops=FALSE, outfile="gl.str",  outpath=getwd(), v=1){
  if(!"genlight" %in% class(gl)){
    warning("Function was designed for genlight objects.")
  }
  # internally convert to genind
  gi <- gl2gi_mvb(gl, v = 0)
  # get the max ploidy of the dataset
  pl <- max(gi@ploidy)
  # get the number of individuals
  S <- nInd(gi)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(gi), each=pl))
  # column of pop ids to write
  if(pops==F){
    popnums <- 1:nPop(gi)
    names(popnums) <- as.character(unique(pop(gi)))
    popcol <- rep(popnums[as.character(pop(gi))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- locNames(gi)
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=nLoc(gi), dimnames=list(NULL,loci)))
  # begin going through loci
  for(L in loci){
    thesegen <- gi@tab[,grep(paste(L, ".", sep=""), dimnames(gi@tab)[[2]], fixed=TRUE), drop=FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(gi)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  # export table
  write.table(tab, file=file.path(outpath, outfile), sep="\t", quote=FALSE, row.names=FALSE)
  if (v==1)  cat(paste("Structure file saved as:", outfile,"\nin folder:",outpath))
}


#' Converts a genlight object to genind object
#' 
#' @param gl -- a genlight object
#' @param v -- level of verbosity. v=0 is silent, v=1 returns more detailed output during conversion.
#' @return A genind object, with all slots filled.
#' @export
#' @author Bernd Gruber (glbugs@@aerg.canberra.edu.au)
#' @details this function uses a faster version of df2genind (from the adgegenet package)

gl2gi <- function(gl, v=1) {

  if (v==1) {
    cat("Start conversion....\n")
    ptm <- proc.time()[3]
    cat("Please note conversion of bigger data sets will take some time!\n" )
    cat("Once finished, we recommend to save the object using >save(object, file=\"object.rdata\")\n")
  }
  #convert to genind....
  x <- as.matrix(gl[,])
  if (v==1) pb <- txtProgressBar(min=0, max=1, style=3, initial=NA)

  for (i in 1:nrow(x))
    {
    for (ii in 1:ncol(x))
      {

      inp <- x[i,ii]
      if (!is.na(inp))
        {
        if (inp==0) x[i,ii] <- "A/A" else if (inp==1) x[i,ii] <- "A/B" else if (inp==2) x[i,ii] <- "B/B"
        }
      }
  if (v==1)   setTxtProgressBar(pb, i/nrow(x))
    }
    
  if (v==1) {
    cat("\nMatrix converted.. Prepare genind object...\n")
    close(pb)
  }
  gen<-df2genind(x[,], sep="/", ncode=1, ind.names=gl@ind.names, pop = gl@pop, ploidy=2)#, probar=probar)
  gen@other <- gl@other

  if (v==1)cat(paste("Finished! Took", round(proc.time()[3]-ptm),"seconds.\n") )
  gen
}

gl2gi_mvb <- function( gl, v=1) { # eg glss from file "glss1-giss1.rdata"
  x <- as.matrix( gl)
  
  xch <- matrix( '', nrow( x), ncol( x))
  
  xch[ is.na( x)] <- 'NA'
  xch[ x==0] <- 'A/A'
  xch[ x==1] <- 'A/B'
  xch[ x==2] <- 'B/B'

  gen <- df2genind( xch, sep='/', ncode=1, ind.names=gl@ind.names, pop=gl@pop, ploidy=2)
  gen@other <- gl@other
  
gen
}

gl2gpop <- function ( glthing, filename) {
 
  afile <- file( filename, open='w')

  # Header line
  cat( 'description\n', file=afile)
  
  # Locus names:  
  cat( locNames( glthing), sep='\n', file=afile)
  
  # Get "genotypes" as integers
  glint <- as.matrix( glthing) # rows=samples, cols=loci; 0, 1, 2 or NA
  # ... found out by table( c( testgl), useNA='always')
  
  glint[ is.na( glint)] <- 3 # for now
  glint <- glint + 1
  glint <- t(glint)  #This is important- fixes the error that the function is making (reading matrix in wrong order)

  outputs <- c( '0303', '0304', '0404', '0000')
  
  each_pop <- pop( glthing)
  
  # Per population:
  for( this_pop in popNames( glthing)) {
    cat( 'pop\n', file=afile)
    
    # pop( glthing)==this_pop
    
    this_recode <- outputs[ glint[ ,each_pop==this_pop]] # Important- this has also changed- pops now assigned to columns, to match transposed matrix. 
    #this_recode <- outputs[ glint[ each_pop==this_pop,]]
    
    this_recode <- matrix( this_recode, ncol=nLoc( glthing), byrow=TRUE)
    
    these_rows <- apply( this_recode, 1, paste, collapse=' ') # one string per sample in this pop
    
    these_rows <- paste( indNames( glthing)[ each_pop==this_pop], these_rows, sep=' , ')
    
    cat( these_rows, sep='\n', file=afile)
  }
  
  close( afile)
}

#fuction for export to ADMIXTURE
gl2Adm <- function(gl, filename = "ADMIXTURE"){
  rep0map <- rep(0, times = adegenet::nLoc(gl))
  mapdf <- data.frame(rep0map, gl$loc.names, rep0map, rep0map)
  write.table(mapdf, file = paste0(filename,".map", sep = ""), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  rep0ped <- rep(0, times = adegenet::nInd(gl))
  rep9ped <- rep(as.numeric(-9), times = adegenet::nInd(gl))
  pop <- gsub(" ", "_", gl$pop)
  ind <- gsub(" ", "_", gl$ind.names)
  peddf <- data.frame(pop, ind, rep0ped, rep0ped, rep0ped, rep9ped)
  glmat <- as.matrix(gl)
  glmat[ is.na( glmat)] <- 3 # for now
  glmat <- glmat + 1
  glmat[glmat==1] <- '2 2'
  glmat[glmat==2] <- '1 2'
  glmat[glmat==3] <- '1 1'
  glmat[glmat==4] <- '0 0'
  peddf <- cbind(peddf, glmat)
  write.table(peddf, file = paste0(filename,".ped", sep = ""), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
}

gl2nhyb <-  function(gl, filename = "gl_NewHybrids.txt"){
  afile <- file( filename, open='w')
  #Header line
  cat( paste0("NumIndivs ",nInd(gl), "\n",
              "NumLoci ", nLoc(gl), "\n",
              "Digits 2\n",
              "Format  Lumped\n", sep = ""), file=afile)
  
  # Locus names:  
  cat( "LocusNames ", locNames( gl), sep=' ', file=afile)
  cat("\n", sep = "", file = afile)
  
  # Get "genotypes" as integers
  glint <- as.matrix(gl) # rows=samples, cols=loci; 0, 1, 2 or NA
  # ... found out by table( c( testgl), useNA='always')
  
  glint[ is.na( glint)] <- 3 # for now
  glint <- glint + 1
  outputs <- c( '0303', '0304', '0404', '0000')
  this_recode <- outputs[glint] 
  this_recode <- matrix( this_recode, ncol=nLoc( gl), byrow=TRUE)
  these_rows <- apply( this_recode, 1, paste, collapse=' ') # one string per sample in this pop
  these_rows <- paste( 1:nInd(gl), "n", indNames( gl), these_rows, sep=' ')
  cat( these_rows, sep='\n', file=afile)
  close( afile)
}

gl2arl <- function ( glthing, filename) {
  
  afile <- file( filename, open='w')
  
  # Header line
  cat( paste0('[Profile]\n
       Title="SNP data"\n
       NbSamples=', nPop(glthing),'\n\n',
      'GenotypicData=1\n
       GameticPhase=0\n
       RecessiveData=1\n
       RecessiveAllele="2"\n
       DataType=STANDARD\n
       LocusSeparator=WHITESPACE\n
       MissingData= "3"
       \n[Data]\n
       \n[[Samples]]\n', sep=""), file=afile)
  
  
  # Get "genotypes" as integers
  glint <- as.matrix( glthing) # rows=samples, cols=loci; 0, 1, 2 or NA
  # ... found out by table( c( testgl), useNA='always')
  
  glint[ is.na( glint)] <- 3 # for now
  
  each_pop <- pop( glthing)
  
  # Per population:
  for( this_pop in popNames( glthing)) {
    
    cat(paste('
    SampleName=','"', this_pop,'"',
    '\nSampleSize=',length(which(each_pop==this_pop)),
    '\nSampleData={\n', sep = ""), file=afile)
      
    
    # pop( glthing)==this_pop

    these_rows <- apply( glint, 1, paste, collapse=' ') # one string per sample in this pop
    
    these_rows <- paste( indNames(glthing)[ each_pop==this_pop],' 1 ', these_rows, sep='')
    
    cat( these_rows, '}', sep='\n', file=afile)
  }
  cat( paste('[[Structure]]\n'), file=afile)
  cat( paste('StructureName="SNP data"\n
             NbGroups=1\n
             Group={\n',
             sep = ""), file=afile)
  for( this_pop in popNames( glthing)) {
  cat( paste('"',this_pop,'"',
             '\n', sep = ""), file=afile)
  }
  cat( paste('}',
             sep = ""), file=afile)
  close( afile)
}


gl.filter.secondaries.kp <- function (gl) {
  x <- gl
  if (class(x) == "genlight") {
    cat("Filtering a genlight object","\n")
  }
  else {
    cat("Fatal Error: Specify either a genlight object","\n")
    stop()
  }
  cat("Total number of SNP loci:", nLoc(x), "\n")
  oldlength <- length(x@other$loc.metrics$AlleleID)
  a <- data.frame(x@other$loc.metrics$CloneID,x@other$loc.metrics$RepAvg, x@other$loc.metrics$AvgPIC, x@other$loc.metrics$AvgCountRef, x@other$loc.metrics$AvgCountSnp, x@other$loc.metrics$AlleleID)
  colnames(a) <- c("CloneID", "RepAvg", "AvgPIC","AvgCountRef", "AvgCountSnp", "AlleleID")
  a <- a[order (a$CloneID, -a$RepAvg, -a$AvgPIC, -a$AvgCountRef, -a$AvgCountSnp),]
  nodups <- a[!duplicated(a$CloneID),]
  nosecondaries <- nodups$AlleleID
  x <- x[, x@other$loc.metrics$AlleleID %in% nosecondaries, treatOther=TRUE] #doesn't subset loc.metrics
  x@other$loc.metrics <- x@other$loc.metrics[x@other$loc.metrics$AlleleID %in% nosecondaries ,] #subsetting loc.metrics
  newlength <- length(x@other$loc.metrics$AlleleID)
  cat("   Number of secondaries:", oldlength-newlength, "\n")
  cat("   Number of loci after secondaries removed:", newlength, "\n")
  return(x)
}


write.fasta <- function ( sexmarkdf, filename) {
  afile <- file( filename, open='w')
  for(i in 1:length(sexmarkdf$LocName)){
  cat(paste('>', sexmarkdf$Species[i], "|",sexmarkdf$Population[i], '|',sexmarkdf$method[i],"|", sexmarkdf$LocName[i], '\n',sexmarkdf$Sequence[i],'\n', sep = ''), sep='', file=afile)
  }
  close( afile)
}


# Input format for DEMEtics
inputformat.FDD <- function (filename, object) 
{
  cat("The table format is transformed...")
  if (object == TRUE) {
    y <- get(filename)
  }
  else {
    y <- read.table(filename, header = TRUE, sep = "")
  }
  
  locs <- y[,3:ncol(y)]
  locvect <- unlist(locs)
  
  inds <- as.character(y[,1])
  inds2 <- c(inds,inds)
  indrep <- rep(inds2, ncol(locs)/2)
  
  pop <- as.character(y[,2])
  pop2 <- c(pop,pop)
  poprep <- rep(pop2, ncol(locs)/2)

  locnum <- c()
  for(i in 1:(ncol(locs)/2)){
    locnum <- c(locnum, paste0("L",i, sep = ""))
  }
  locnum2 <- mapply(rep, times = (nrow(locs)*2), x = locnum )
  locnum3 <- as.vector(t(locnum2))
  
  
  df <- data.frame(individual = indrep, population = poprep, fragment.length = locvect, locus = locnum3)
  
  return(df)
    
}


gl.read.silico.count.FDD <- function (datafile = f.in, topskip = 5, nmetavar = 15, nas = "-", ind.metafile = fm) 
{
  cat("Reading data from file:", datafile, "\n")
  cat("  This may take some time, please wait!\n")
  x <- read.csv(datafile, na.strings = nas, skip = topskip, 
                check.names = FALSE)
  cat("The following metadata for loci was identified: ", names(x[1:nmetavar]), 
      "\n")
  if (any(names(x) == "CloneID")) {
    cat("  includes key variable CloneID\n")
  }
  else {
    cat("Fatal Error: Dataset does not include key variable CloneID!\n")
    stop()
  }
  if (any(names(x) == "Reproducibility")) {
    cat("  includes key variable Reproducibility\n")
  }
  else {
    cat("Warning: Dataset does not include variable Reproducibility which may limit your filtering options in later analyses!\n")
  }
  ind.names <- colnames(x)[(nmetavar + 1):ncol(x)]
  cat("Data identified for ", ncol(x) - nmetavar, "individuals, ", 
      nrow(x), "loci")
  if (length(ind.names) != length(unique(ind.names))) {
    cat("Warning: Specimen names are not unique!\n")
    cat("         Duplicated names:\n")
    noccur <- table(ind.names)
    cat(paste("              ", names(noccur[noccur > 1])), 
        "\n")
    cat("         Rendering specimen names unique with sequential suffix _1, _2 etc\n")
    ind.names <- make.unique(ind.names, sep = "_")
  }
  snpdata <- x[, (nmetavar + 1):ncol(x)]
  locus.metadata <- x[, 1:nmetavar]
  # if (max(snpdata, na.rm = TRUE) != 1 || min(snpdata, na.rm = TRUE) != 
  #     0) {
  #   cat("Fatal Error: SNP data must be 0 or 1!\n")
  #   stop()
  # }
  nind <- ncol(snpdata)
  nloci <- nrow(locus.metadata)
  cat("\nStarting conversion to genInd object ....\n")
  cat("Please note conversion of bigger data sets will take some time!\n")
  locname <- x$CloneID
  if (length(locname) != length(unique(locname))) {
    cat("Warning: Locus names are not unique!\n")
    cat("         Duplicated names:\n")
    noccur <- table(locname)
    cat(paste("              ", names(noccur[noccur > 1])), 
        "\n")
    cat("         Rendering locus names unique with sequential suffix _1, _2 etc\n")
    locname <- make.unique(as.character(x$CloneID), sep = "_")
  }
  glind <- new("genind", tab = t(snpdata), ploidy = 2, ind.names = colnames(snpdata), 
               loc.names = locname, parallel = F)
  colnames(glind@tab) <- locname
  glind@other$loc.metrics <- locus.metadata
  if (!is.null(ind.metafile)) {
    cat("Adding population assignments and other additional individual metadata from file :", 
        ind.metafile, "\n")
    ind.metadata <- read.csv(ind.metafile, header = T, stringsAsFactors = T)
    ind.metadata$id <- gsub("^\\s+|\\s+$", "", ind.metadata$id)
    id.col = match("id", names(ind.metadata))
    if (is.na(id.col)) {
      cat("Fatal Error: No id column present!\n")
      stop()
    }
    else {
      if (length(ind.metadata[, id.col]) != length(names(snpdata))) {
        cat("Ids of covariate file does not match the number of ids in the genetic file. Maybe this is fine if a subset matches.\n")
      }
      ord <- match(names(snpdata), ind.metadata[, id.col])
      ord <- ord[!is.na(ord)]
      if (length(ord) > 1 & length(ord) <= nind) {
        cat(paste("Ids of covariate file (at least a subset of) are matching!\nFound ", 
                  length(ord == nind), "matching ids out of", 
                  nrow(ind.metadata), "ids provided in the covariate file. Subsetting snps now!.\n "))
        ord2 <- match(ind.metadata[ord, id.col], indNames(glind))
        glind <- glind[ord2, ]
      }
      else {
        cat("Fatal Error: Ids in files ", datafile, "and ", 
            ind.metafile, " do not match!\n\n")
        stop()
      }
    }
    pop.col = match("pop", names(ind.metadata))
    if (is.na(pop.col)) {
      cat("Warning: No pop column present\n")
    }
    else {
      pop(glind) <- as.factor(ind.metadata[ord, pop.col])
      cat("Populations assigned to individuals\n")
    }
    lat.col = match("lat", names(ind.metadata))
    lon.col = match("lon", names(ind.metadata))
    if (is.na(lat.col)) {
      cat("Warning: No lat column present\n")
    }
    if (is.na(lon.col)) {
      cat("Warning: No lon column present\n")
    }
    if (!is.na(lat.col) & !is.na(lon.col)) {
      glind@other$latlong <- ind.metadata[ord, c(lat.col, lon.col)]
      rownames(glind@other$latlong) <- ind.metadata[ord, id.col]
      cat("Added latlon data.\n")
    }
    known.col <- names(ind.metadata) %in% c("id", "pop", 
                                            "lat", "lon")
    other.col <- names(ind.metadata)[!known.col]
    if (length(other.col > 0)) {
      glind@other$ind.metrics <- ind.metadata[ord, other.col, 
                                        drop = FALSE]
      rownames(glind@other$ind.metrics) <- ind.metadata[ord, 
                                                  id.col]
      cat(paste("Added ", other.col, " to the other$ind.metrics slot.\n"))
    }
  }
  cat("GenInd object created")
  return <- glind
}




# 
# gl.filter.callrate.FD<-
# function (gl, method = "loc", threshold = 1) 
# {
#   x <- glneoFS
#   if (class(x) == "genlight") {
#     cat("Reporting for a genlight object\\n")
#   }
#   else if (class(x) == "genind") {
#     cat("Reporting for a genind object\\n")
#   }
#   else {
#     cat("Fatal Error: Specify either a genlight or a genind object\\n")
#     stop()
#   }
#   cat("Note: Missing values most commonly arise from restriction site mutation.\\n\\n")
#   if (method == "loc") {
#     n0 <- nLoc(x)
#     cat("Initial no. of loci =", n0, "\\n")
#     if (class(x) == "genlight") {
#       x2 <- x[, glNA(x, alleleAsUnit = FALSE) <= ((1 - 
#                                                      threshold) * nInd(x))]
#       x2@other$loc.metrics <- x@other$loc.metrics[glNA(x, 
#                                                        alleleAsUnit = FALSE) <= ((1 - threshold) * nInd(x)), 
#                                                   ]
#       cat("  No. of loci deleted =", (n0 - nLoc(x2)), "\\n")
#     }
#     else if (class(x) == "genind") {
#       x2 <- x[, (colSums(is.na(tab((x))))/nInd(x)) <= (1 - 
#                                                          threshold)]
#       idx <- which((colSums(is.na(tab((x))))/nInd(x)) <= 
#                      (1 - threshold))
#       x2@other$loc.metrics <- x@other$loc.metrics[c(idx), 
#                                                   ]
#       cat("  No. of loci deleted =", (n0 - nLoc(x2)), "\\n")
#     }
#     else {
#       cat("Fatal Error: genlight or genind objects required for call rate filtering!\\n")
#       stop()
#     }
#   }
#   else if (method == "ind") {
#     n0 <- nInd(x)
#     cat("Initial no. of individuals =", n0, "\\n")
#     ind.call.rate <- 1 - rowSums(is.na(as.matrix(x)))/nLoc(x)
#     if (sum(ind.call.rate >= threshold) == 0) 
#       stop(paste("Maximum individual call rate =", max(ind.call.rate), 
#                  ". Nominated threshold of", threshold, "too stringent.\\n No individuals remain.\\n"))
#     x2 <- x[ind.call.rate >= threshold, ]
#     if (class(x) == "genlight") {
#       cat("Filtering a genlight object\\n  no. of individuals deleted =", 
#           (n0 - nInd(x2)), "\\nIndividuals retained =", 
#           nInd(x2), "\\n")
#     }
#     if (class(x) == "genind") {
#       cat("Filtering a genind object\\n  No. of individuals deleted =", 
#           (n0 - nInd(x2)), "\\n  Individuals retained =", 
#           nInd(x2), "\\n")
#     }
#     if (any(ind.call.rate <= threshold)) {
#       x3 <- x[ind.call.rate <= threshold, ]
#       if (length(x3) > 0) {
#         cat("List of individuals deleted because of low call rate\\n", 
#             indNames(x3), "\\n")
#         cat("   from populations\\n", as.character(pop(x3)), 
#             "\\n")
#       }
#     }
#   }
#   else {
#     cat("Fatal Error: the method parameter must be specified as either loc or ind\\n")
#     stop()
#   }
#   cat("Summary of filtered dataset\\n")
#   cat(paste("  Call Rate >", threshold, "\\n"))
#   cat(paste("  No. of loci:", nLoc(x2), "\\n"))
#   cat(paste("  No. of individuals:", nInd(x2), "\\n"))
#   cat(paste("  No. of populations: ", length(levels(factor(pop(x2)))), 
#             "\\n"))
#   return(x2)
# }
