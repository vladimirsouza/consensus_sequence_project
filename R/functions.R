


#' load data into R
#'
#' Load data from bam and fasta files at once
#'
#' @param bam_file #### to add ####
#' @param fasta_file #### to add ####
#' @param what #### to add ####
#' @param tag #### to add ####
#'
#' @return #### to add ####
#' @export
#'
#' @examples #### to add ####
load_pacbio_data <- function(bam_file, fasta_file, what="qname", tag=c("zm", "pw")){
  require(Rsamtools)
  
  dat <- scanBam( bam_file, param=ScanBamParam(what=what, tag=tag) ) [[1]]
  dat <- c(
    dat[names(dat) != "tag"],
    dat[["tag"]],
    seq=readDNAStringSet(fasta_file)
  )
  stopifnot( setequal(names(dat$seq), dat$qname) )
  dat$seq <- dat$seq[dat$qname]
  dat
}












#' Get pulse widths sorted when same mucleotides in a row
#'
#' @param pws #### to add ####
#' @param seqs #### to add ####
#'
#' @return #### to add ####
#' @export
#'
#' @examples #### to add ####
sort_pws <- function(pws, seqs, is_reverse=FALSE, gap_to="right"){
  require(dplyr)
  require(Biostrings)
  
  stopifnot( gap_to %in% c("right", "left") )
  
  if( !is.character(seqs) )
    seqs <- as.character(seqs) %>%
      unname
  x <- strsplit(seqs, NULL) %>%
    RleList %>%
    runLength
  y <- vector("list", length(pws))
  
  if(gap_to == "right"){
    op_noRev <- `-`
    op_rev <- `+`
  }else{
    op_noRev <- `+`
    op_rev <- `-`
  }
  
  k <- unlist(pws) %>%
    { max(.) * 10 }
  
  y[!is_reverse] <- mapply( function(x, pw){
    rep(seq_along(x), x) %>%
      op_noRev( pw/k ) %>%
      order %>%
      pw[.]
  },
  x[!is_reverse],
  pws[!is_reverse],
  SIMPLIFY=FALSE )
  
  y[is_reverse] <- mapply( function(x, pw){
    rep(seq_along(x), x) %>%
      op_rev( pw/k ) %>%
      order %>%
      pw[.]
  },
  x[is_reverse],
  pws[is_reverse],
  SIMPLIFY=FALSE )
  
  y
}



# This is a function to create the colors matrix needed by DECIPHER, based on pulse widths
#' Title
#'
#' @param seqs #### to add ####
#' @param pws #### to add ####
#' @param pw_levels #### to add ####
#'
#' @return #### to add ####
#' @export
#'
#' @examples #### to add ####
decipher_colours <- function(seqs, pws, pw_levels=NULL){
  require(dplyr)
  if( !is.character(seqs) ) seqs <- as.character(seqs)
  stopifnot( !is.null(pw_levels) )
  stopifnot( length(pw_levels)==4 )
  te <- list(c("A","G"), NULL, c("A","C"))
  mapply(function(s,p){
    param <- c( list(p), setNames( list(.9,.75,.5,0), pw_levels ) )
    y <- do.call( dplyr::recode, param)
    t( sapply(te, function(x){
      replace(y, strsplit(s, NULL)[[1]] %in% x, 1)
    }) )
  },seqs, pws, SIMPLIFY=FALSE)
}





#' Identify weird holes
#'
#' @param aln 
#' @param fa_hls_i 
#'
#' @return
#' @export
#'
#' @examples
identify_weird_holes <- function(aln, fa_hls_i){
  k <- setNames( as.vector(strand(aln)), mcols(aln)$qname )
  k <- k[ names(fa_hls_i) ]
  k <- split( k, rep(1:2, len=length(k)) )
  k <- lapply( k, function(u) unique(na.omit(u)) )
  if( any( lengths(k) == 0 ) ){
    k <- lengths(k)
    return(
      ! identical( unname(k[k != 0]), 1L )
    )
  }else{
    return(
      (k[[1]] == k[[2]]) | any(lengths(k) != 1)
    )
  }
}

#' Transform pulse widths
#'
#' @param pws 
#' @param max_pw 
#'
#' @return
#' @export
#'
#' @examples
transform_pws <- function(pws, max_pw){
  p <- unlist(pws, use.names=FALSE)
  p[p > max_pw] <- max_pw
  relist(p, pws)
}

#' Insert2
#'
#' @param string 
#' @param positions 
#' @param insertions 
#'
#' @return
#' @export
#'
#' @examples
insert2 <- function(string, positions, insertions){
  paste(
    insert( strsplit(string, NULL)[[1]], positions, insertions ),
    collapse=""
  )
}

#' Sort pulse widths without transformation
#' 
#' Sort the pulse widths of same nucleotides in a row
#'
#' @param pws 
#' @param seqs 
#' @param gap_to 
#'
#' @return
#' @export
#'
#' @examples
only_sort_pws <- function(pws, seqs, gap_to="left"){
  stopifnot( gap_to %in% c("right", "left") )
  x <- runLength( RleList(seqs) )
  if(gap_to == "left"){
    operation <- `+`
  }else{
    operation <- `-`
  }
  mapply(
    function(xi, pw){
      pw[order(
        operation(
          rep(seq_along(xi), xi),
          pw / (max(pw)*10)
        )
      )]
    },
    x, pws, SIMPLIFY=FALSE
  )
}










#' Find subread of reference
#'
#' @param temp_rootDirs_i 
#' @param fa_hls_i 
#' @param fa_hls_files_i 
#' @param cmds_ovl_i 
#' @param ovl_files_i 
#'
#' @return
#' @export
#'
#' @examples
find_subread_of_reference <- function(temp_rootDirs_i, fa_hls_i, fa_hls_files_i, cmds_ovl_i, ovl_files_i){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  # create temporary directory
  dir.create(temp_rootDirs_i, showWarnings=FALSE)
  
  # write one fasta per hole
  writeXStringSet(fa_hls_i, fa_hls_files_i)
  
  # overlap subreads
  system(cmds_ovl_i)
  
  # find the best subread of reference
  suppressWarnings(
    ovl <- tryCatch(
      scanBam(ovl_files_i, param=sbp_ovl) [[1]],
      error=function(e){NULL}
    )
  )
  if( is.null(ovl) ){
    return(NULL)
  }
  ovl <- with(
    ovl,
    data.frame(
      rname=as.vector(rname),
      qname=qname,
      AS=tag$AS,
      stringsAsFactors=FALSE
    )
  )
  ovl <- ovl[ ovl$rname != ovl$qname, ]
  if( nrow(ovl) == 0 ){
    return(NULL)
  }
  
  k <- format( seq( nrow(ovl) ) )
  k <- matrix(
    sub(".*?_", "",
        sort( c(paste(k, ovl$rname, sep="_"),
                paste(k, ovl$qname, sep="_")) )
    ),
    ncol=2, byrow=TRUE
  )
  k <- paste(k[,1], k[,2])
  ovl <- sapply(
    unname( split(ovl, k) ),
    function(x){
      x[which.max(x$AS),]
    }
  )
  ovl <- data.frame(
    sr=unlist( t(ovl[1:2,]), use.names=FALSE ),
    AS=unlist( ovl[3,], use.names=FALSE ),
    stringsAsFactors=FALSE
  )
  ovl <- sapply(
    split(ovl, ovl$sr),
    function(x) sum(x$AS)
  )
  ovl <- names( ovl[which.max(ovl)] )
  
  return( fa_hls_i[ovl] )
}







#' Get truth and training matrix for insertions
#'
#' @param Seqs 
#' @param PWs 
#' @param Gaps 
#' @param temp_rootDirs_i 
#' @param nCover 
#' @param indx 
#'
#' @return
#' @export
#'
#' @examples
truth_and_training_matrix_for_insertions <- function(Seqs, PWs, Gaps, temp_rootDirs_i, nCover, indx){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  # when there is only a single insertion in a position, there is no msa, just
  # compare the nucleotide weights with the sum of all gap weights of that position
  k <- lengths(PWs) == 1
  PWs <- transform_pws(PWs, max_pw)
  nts <- do.call(
    "rbind",
    mapply(function(s, i){
      k <- matrix("-", width(s), nCover)
      k[,i] <- strsplit(as.character(s), NULL)[[1]]
      k
    }, Seqs[k], indx[k], SIMPLIFY=FALSE)
  )
  aln_pws <- matrix(gap_weight_ins, nrow(nts), ncol(nts))  ##### can i use gap_weight_ins here, or should it be 1??
  aln_pws[ nts != "-" ] <- unlist(PWs[k])
  truth <- rep("-", nrow(aln_pws))
  
  Seqs <- Seqs[!k]
  PWs <- PWs[!k]
  Gaps <- Gaps[!k]
  indx <- indx[!k]
  
  if(length(PWs) > 0){
    # when there are more than one insertion in a position, we compute the consensus
    # from a msa
    temp_dir <- tempfile(tmpdir=temp_rootDirs_i)
    dir.create(temp_dir)
    fls <- tempfile(
      as.character( seq_along(Seqs) ),
      temp_dir,
      ".fasta"
    )
    
    nts1 <- mapply(function(s, f, j){
      writeXStringSet(s, f)
      k <- as.matrix(
        msaClustalW(f, gapOpening=0, gapExtension=0, type="dna", order="input")
      )
      nts <- matrix("-", nCover, ncol(k))
      nts[j,] <- k
      t(nts)
    }, Seqs, fls, indx, SIMPLIFY=FALSE)
    aln_pws1 <- do.call(
      "rbind",
      mapply(function(u, p){
        replace(
          matrix(gap_weight_ins, nrow(u), ncol(u)),
          u != "-",
          unlist(p)
        )
      }, nts1, PWs, SIMPLIFY=FALSE)
    )
    nts1 <- do.call("rbind", nts1)
    
    k <- c("A", "C", "G", "T", "-")
    cons <- sapply( k, function(u){
      rowSums( ifelse(nts1 == u, aln_pws1, 0) )
    }, USE.NAMES=FALSE)
    truth <- c(
      truth,
      k[ apply(cons, 1, which.max) ]
    )
    
    # using_pws <- FALSE
    # if(using_pws){
      aln_pws[nts == "-"] <- 1
      aln_pws1[nts1 == "-"] <- 1
    # }else{
    #   aln_pws[] <- 1
    #   aln_pws1[] <- 1
    # }
    
    nts <- rbind(nts, nts1)
    aln_pws <- rbind(aln_pws, aln_pws1)
  }
  
  return(
    list(nts=nts,
         pws=aln_pws,
         truth=truth,
         isIns=rep( 1, length(truth)) )
  )
}









#' Get truth and training matrix for NON insertions
#'
#' @param aln_no_ins 
#' @param no_ins_pws 
#' @param clipping_sign 
#'
#' @return
#' @export
#'
#' @examples
truth_and_training_matrix_for_NON_insertions <- function(aln_no_ins, no_ins_pws, clipping_sign="+"){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  aln_no_ins <-t(
    as.matrix( unname(aln_no_ins) )
  )
  aln_pws <- matrix(gap_weight_del, nrow(aln_no_ins), ncol(aln_no_ins))
  k <- !(aln_no_ins %in% c(clipping_sign, "-"))
  aln_pws[k] <- unlist( transform_pws(no_ins_pws, max_pw) )
  aln_pws[ aln_no_ins == clipping_sign ] <- NA
  
  k <- c("A","C","G","T","-")
  cons <- sapply(k, function(u){
    rowSums( aln_pws * (aln_no_ins == u), na.rm=TRUE )
  }, USE.NAMES=FALSE)
  
  truth <- k[ apply(cons, 1, which.max) ]
  
  # using_pws <- FALSE
  # if(using_pws){
    aln_pws[aln_no_ins == "-"] <- 1
  # }else{
  #   aln_pws[] <- 1
  # }
  
  return(
    list(nts=aln_no_ins,
         pws=aln_pws,
         truth=truth,
         isIns=rep( 0, length(truth)) )
  )
}









#' Find the truth
#'
#' @param aln 
#' @param nCover 
#' @param fa_ref 
#' @param pws_i 
#' @param aln_files_i 
#' @param temp_rootDirs_i 
#'
#' @return
#' @export
#'
#' @examples
find_the_truth <- function(aln, nCover, fa_ref, pws_i, aln_files_i, temp_rootDirs_i){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  ### for the training, use only regions that the coverage is `n`
  cover <- coverage(aln)[[1]]
  
  # find `n`-covered regions
  n_cover <- cover == nCover
  n_cover_val <- runValue(n_cover)
  
  # it might doesn't contain any `n`-coverage region at all
  if( any(n_cover_val) ){
    # if there are more than one `n`-covered region, take the longest one
    if( length(n_cover_val[n_cover_val]) > 1 ){
      n_cover_len <- runLength(n_cover)
      max_len <- max(n_cover_len[n_cover_val])
      runValue(n_cover) [n_cover_val] [n_cover_len[n_cover_val] != max_len] <- FALSE
    }
    # the region where it will be computed the consensus sequence
    n_cover <- IRanges(n_cover)
    k <- width(n_cover) >= min_obs_training
  }else{
    k <- FALSE
  }
  
  # if there is a `n`-coverage region and it's long enough
  if(k){
    # if a subread is outside the `n`-coverege region, eliminate it
    qry_aln_rgs <- ranges(aln)
    sr_ovl_n_cover <- queryHits( findOverlaps( qry_aln_rgs, n_cover ) )
    
    stopifnot( length(sr_ovl_n_cover) == length(aln) ) #### <<<<<<<<----------===
    
    aln <- aln[sr_ovl_n_cover]
    qry_aln_rgs <- qry_aln_rgs[sr_ovl_n_cover]
    
    # all insertions
    rrs <- cigarRangesAlongReferenceSpace( cigar(aln), ops="I", pos=start(aln) )
    
    # test whether there are insertions inside the high coverage region
    ins_pos <- start( unlist(rrs) )
    is_in_n_cover <- (ins_pos > start(n_cover)) & (ins_pos <= end(n_cover))
    
    # correct pulse width sequences of subreads in reverse strand
    pws_i <- pws_i[ mcols(aln)$qname ]
    k <- as.vector(strand(aln) == "-")
    pws_i[k] <- lapply( pws_i[k], rev )
    
    # sort pulse widths of consecutive repeated nucleotides (for minimap2, insertions to the left)
    pws_i <- RleList(
      only_sort_pws(pws_i, mcols(aln)$seq)
    )
    
    if( any(is_in_n_cover) ){
      # keep only the insertion positions inside the high coverage region
      ins_pos <- ins_pos[is_in_n_cover]
      
      # number of insertions per position
      num_ins <- as.integer(table(ins_pos))
      
      # number of gaps per position
      num_gap <- as.integer(cover[IRanges(sort(unique(ins_pos)), width=1)]) - num_ins
      
      # insertion sequences per position
      qrs <- cigarRangesAlongQuerySpace( cigar(aln), ops="I" )
      qrs_unlist <- unlist(qrs) [is_in_n_cover]
      indx <- rep( seq_along(qrs), lengths(qrs) ) [is_in_n_cover]
      ins_seq <- compact(
        subseq(
          mcols(aln)$seq[indx],
          start(qrs_unlist),
          end(qrs_unlist)
        )
      )
      ins_seq <- split(
        # setNames( ins_seq, seq_along(ins_seq) ),
        setNames( ins_seq, indx ),
        ins_pos
      )
      
      # get pulse widths of insertions
      k <- width(qrs_unlist)
      ins_pws <- split(
        unname(
          split(
            as.vector(
              unlist(
                pws_i[unname(
                  split(
                    qrs_unlist,
                    factor( indx, levels=seq_along(qrs) )
                  )
                )]
              )
            ),
            rep(seq_along(k), k)
          )
        ),
        ins_pos
      )
      
      indx <- split(indx, ins_pos)
      
      # Consensus of insertions and input matrix for svm training using all subreads:
      svm_args <- truth_and_training_matrix_for_insertions(ins_seq, ins_pws, num_gap, temp_rootDirs_i, nCover, indx)
      
    }else{
      svm_args <- NULL
    }
    
    # to take only the pulse widths of those nucleotedes that are inside the high-
    # coveraged region, we need to convert the ranges of `n_cover` (that are based
    # on the reference) to be based on the queries, getting different ranges for each
    # subread. But before that, we need to fix the ranges when the position where a
    # subread start/stop aligning to the reference is inside the high-coveraged region.
    start(qry_aln_rgs)[ start(qry_aln_rgs) < start(n_cover) ] <- start(n_cover)
    end(qry_aln_rgs)[ end(qry_aln_rgs) > end(n_cover) ] <- end(n_cover)
    qry_aln_rgs <- pmapToAlignments(
      qry_aln_rgs,
      setNames(
        aln,
        rep( "x", length(aln) )
      )
    )
    
    # get pulse widths of non-insertion nucleotides and in the high-coveraged region
    no_ins_pws <- pws_i [gaps(
      cigarRangesAlongQuerySpace( cigar(aln), ops=c("S","H","I") ),
      start(qry_aln_rgs), end(qry_aln_rgs)
    )]
    
    # consensus, but not about insertions
    invisible(
      indexBam(
        sortBam(
          aln_files_i,
          sub(".bam", "", aln_files_i)
        )
      )
    )
    aln_no_ins <- stackStringsFromBam(
      BamFile(aln_files_i),
      param=GRanges( names(fa_ref), n_cover),
      use.names=TRUE
    ) [mcols(aln)$qname]
    
    svm_args <- mapply(
      function(a, b){
        if( is.matrix(a) ){
          rbind(a, b)
        }else{
          c(a, b)
        }
      },
      svm_args,
      truth_and_training_matrix_for_NON_insertions(aln_no_ins, no_ins_pws)
    )
    
    
    
  }else{
    unlink(temp_rootDirs_i, recursive=TRUE)
    return(NULL)
  }
}



#' Get variables for the training
#'
#' @param fa_hls_i 
#' @param pws_i 
#' @param fa_hls_files_i 
#' @param ovl_files_i 
#' @param cmds_ovl_i 
#' @param fa_ref_files_i 
#' @param cmds_aln_i 
#' @param aln_files_i 
#' @param temp_rootDirs_i 
#' @param transform_pws 
#' @param only_sort_pws 
#' @param find_subread_of_reference 
#' @param truth_and_training_matrix_for_insertions 
#' @param truth_and_training_matrix_for_NON_insertions 
#' @param find_the_truth 
#'
#' @return
#' @export
#'
#' @examples
variables_for_the_training <- function(fa_hls_i, pws_i, fa_hls_files_i,
                                       ovl_files_i, cmds_ovl_i, fa_ref_files_i,
                                       cmds_aln_i, aln_files_i, temp_rootDirs_i,
                                       transform_pws, only_sort_pws,
                                       find_subread_of_reference,
                                       truth_and_training_matrix_for_insertions,
                                       truth_and_training_matrix_for_NON_insertions,
                                       find_the_truth){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  fa_hls_i <- compact(fa_hls_i) ########### <<<<<<<<-------------====
  
  # find subread of reference
  fa_ref <- find_subread_of_reference(temp_rootDirs_i, fa_hls_i, fa_hls_files_i,
                                      cmds_ovl_i, ovl_files_i)
  
  
  if( is.null(fa_ref) ){
    unlink(temp_rootDirs_i, recursive=TRUE)
    return(NULL)
  }
  
  
  
  # write the subread of reference as a fasta file
  writeXStringSet(fa_ref, fa_ref_files_i)
  
  
  
  # align all subreads to the subread of reference
  system(cmds_aln_i)
  aln <- readGAlignments( aln_files_i, param=ScanBamParam(what=c("qname", "seq")) )
  
  
  # If before the iteration only one sequence could be aligned, that means that it is the reference itself. So, there are no consensus to do.
  if( length(aln) == 1 ){
    unlink(temp_rootDirs_i, recursive=TRUE)
    return(NULL)
  }
  
  
  # if there are two consecutive subreads with the same direction, it's a weird
  # hole and we discarde it.
  if( identify_weird_holes(aln, fa_hls_i) ){
    unlink(temp_rootDirs_i, recursive=TRUE)
    return(NULL)
  }
  
  
  nCover <- length(fa_hls_i)
  
  return(
    find_the_truth(aln, nCover, fa_ref, pws_i, aln_files_i, temp_rootDirs_i)
  )
  
}





#' Train the models
#'
#' @param N 
#' @param x_training 
#' @param truth 
#' @param isIns 
#' @param THREADS 
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
train_the_models <- function(N, x_training, truth, isIns, THREADS, model){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  nts <- x_training$nts
  pws <- x_training$pws
  truth <- factor(truth)
  
  model <- switch(
    model,
    "svm"=svm,
    "rf"=randomForest,
    stop("model must be 'svm' or 'rf'")
  )
  
  registerDoParallel(cores=THREADS)
  f <- function(n, nts, pws, truth, isIns, model){
    k <- sample( ncol(pws), n)
    nam <- c("A","C","G","T","-")
    x <- sapply(nam, function(u){
      rowSums( pws[,k] * (nts[,k] == u), na.rm=TRUE )
    }, USE.NAMES=FALSE)
    x <- cbind(x, isIns)
    colnames(x) <- NULL
    
    model(x, truth)
  }
  
  foreach(n=N) %dopar%
    f(n, nts, pws, truth, isIns, model)
}







#' Get consensus of insertions
#'
#' @param Seqs 
#' @param PWs 
#' @param Gaps 
#' @param n 
#' @param temp_rootDirs_i 
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
consensus_insertions <- function(Seqs, PWs, Gaps, n, temp_rootDirs_i, model){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  # if there is only a single insertion in a position, and i always use 3 subreads
  # at the least, we can say that the consensus is a gap at once
  k <- lengths(PWs) != 1
  Seqs <- Seqs[k]
  PWs <- PWs[k]
  Gaps <- Gaps[k]
  
  PWs <- transform_pws(PWs, max_pw)
  
  # when there are more than one insertion in a position, we compute the consensus
  # from a msa
  if(length(PWs) > 0){
    temp_dir <- tempfile(tmpdir=temp_rootDirs_i)
    dir.create(temp_dir, showWarnings=FALSE)
    fls <- tempfile(
      as.character( seq_along(Seqs) ),
      temp_dir,
      ".fasta"
    )
    
    nts1 <- mapply(function(s, f){
      writeXStringSet(s, f)
      k <- as.matrix(
        msaClustalW(f, gapOpening=0, gapExtension=0, type="dna", order="input")
      )
      nts <- matrix("-", n, ncol(k))
      nts[ seq(nrow(k)), ] <- k
      t(nts)
    }, Seqs, fls, SIMPLIFY=FALSE)
    x_svm_pos <- rep( names(PWs), sapply(nts1, nrow) )
    aln_pws1 <- do.call(
      "rbind",
      mapply(function(u, p){
        replace(
          matrix(1, nrow(u), ncol(u)),
          u != "-",
          unlist(p)
        )
      }, nts1, PWs, SIMPLIFY=FALSE)
    )
    nts1 <- do.call("rbind", nts1)
    
    k <- c("A", "C", "G", "T", "-")
    cons <- sapply( k, function(u){
      rowSums( ifelse(nts1 == u, aln_pws1, 0) )
    }, USE.NAMES=FALSE)
    if( !is.matrix(cons) ) cons <- matrix(cons, nrow=1)
    cons <- cbind(cons, 1)
    
    cons <- sapply(
      split(
        as.character( predict(model, cons) ),
        x_svm_pos
      ),
      paste, collapse=""
    )
    cons <- gsub("-", "", cons)
    return( cons[cons != ""] )
  }else{
    return(NULL)
  }
}






#' Get insertion consensus per coverage area
#'
#' @param cover_start_i 
#' @param cover_end_i 
#' @param n 
#' @param model 
#' @param ins_pos 
#' @param cover 
#' @param aln 
#' @param pws_i 
#' @param temp_rootDirs_i 
#'
#' @return
#' @export
#'
#' @examples
insertion_consensus_per_coverage <- function(cover_start_i, cover_end_i, n, model, ins_pos, cover, aln, pws_i, temp_rootDirs_i){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  # test whether there are insertions inside the high coverage region
  is_in_n_cover <- (ins_pos > cover_start_i) & (ins_pos <= cover_end_i)
  
  if( any(is_in_n_cover) ){
    # keep only the insertion positions inside the high coverage region
    ins_pos <- ins_pos[is_in_n_cover]
    
    # number of insertions per position
    num_ins <- as.integer( table(ins_pos) )
    
    # number of gaps per position
    num_gap <- as.integer(cover[IRanges(sort(unique(ins_pos)), width=1)]) - num_ins
    
    # insertion sequences per position
    qrs <- cigarRangesAlongQuerySpace( cigar(aln), ops="I" )
    qrs_unlist <- unlist(qrs) [is_in_n_cover]
    indx <- rep( seq_along(qrs), lengths(qrs) ) [is_in_n_cover]
    ins_seq <- compact(
      subseq(
        mcols(aln)$seq[indx],
        start(qrs_unlist),
        end(qrs_unlist)
      )
    )
    ins_seq <- split(
      # setNames( ins_seq, seq_along(ins_seq) ),
      setNames( ins_seq, indx ),
      ins_pos
    )
    
    
    # get pulse widths of insertions
    k <- width(qrs_unlist)
    ins_pws <- split(
      unname(
        split(
          as.vector(
            unlist(
              pws_i[unname(
                split(
                  qrs_unlist,
                  factor( indx, levels=seq_along(qrs) )
                )
              )]
            )
          ),
          rep(seq_along(k), k)
        )
      ),
      ins_pos
    )
    
    
    # Consensus of insertions:
    # To save time, don't compute the consensus of positions that it's obvious.
    # If in a position the number of gaps multiplied by the gap_weight is higher than
    # the number of any insertion multiplied by nt_weight, let the consensus be a gap.
    # k <- (num_gap * gap_weight_ins) < (num_ins * max_pw)
    
    return(
      consensus_insertions(ins_seq, ins_pws, num_gap, n, temp_rootDirs_i, model)
    )
    
  }else{
    return(NULL)
  }
}







#' Get consensus of NON insertions
#' 
#' Get consensus of those alingment positions where here are deletions or mismatches.
#'
#' @param aln_no_ins 
#' @param no_ins_pws 
#' @param model 
#' @param clipping_sign 
#'
#' @return
#' @export
#'
#' @examples
consensus_NON_insertions <- function(aln_no_ins, no_ins_pws, model, clipping_sign="+"){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  aln_no_ins <-t(
    as.matrix( unname(aln_no_ins) )
  )
  aln_pws <- matrix(1, nrow(aln_no_ins), ncol(aln_no_ins))
  k <- !(aln_no_ins %in% c(clipping_sign, "-"))
  aln_pws[k] <- unlist( transform_pws(no_ins_pws, max_pw) )
  aln_pws[ aln_no_ins == clipping_sign ] <- NA
  
  k <- c("A","C","G","T","-")
  cons <- sapply(k, function(u){
    rowSums( aln_pws * (aln_no_ins == u), na.rm=TRUE )
  }, USE.NAMES=FALSE)
  if( !is.matrix(cons) ) cons <- matrix(cons, nrow=1)
  cons <- cbind(cons, 0)
  
  return(
    paste(
      as.character( predict(model,cons) ),
      collapse=""
    )
  )
}





#' Get non insertion consensus per coverage region
#'
#' @param aln_no_ins 
#' @param pws_i 
#' @param aln 
#' @param qry_aln_rgs 
#' @param cover_start_i 
#' @param cover_end_i 
#' @param model 
#'
#' @return
#' @export
#'
#' @examples
non_insertion_consensus_per_coverage <- function(aln_no_ins, pws_i, aln, qry_aln_rgs, cover_start_i, cover_end_i, model){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  # sequences region of a coverage
  aln_no_ins <- subseq(aln_no_ins, cover_start_i, cover_end_i)
  
  # eliminate subreads that don't overlap the coverage region
  k <- subjectHits(
    findOverlaps( IRanges(cover_start_i, cover_end_i), qry_aln_rgs )
  )
  qry_aln_rgs <- qry_aln_rgs[k]
  pws_i <- pws_i[k]
  aln <- aln[k]
  aln_no_ins <- aln_no_ins[k]
  
  # pulse widths of a coverage
  # start(qry_aln_rgs)[ start(qry_aln_rgs) < cover_start_i & end(qry_aln_rgs) > cover_start_i ] <- cover_start_i
  # end(qry_aln_rgs)[ end(qry_aln_rgs) > cover_end_i & start(qry_aln_rgs) < cover_end_i ] <- cover_end_i
  start(qry_aln_rgs)[ start(qry_aln_rgs) < cover_start_i ] <- cover_start_i
  end(qry_aln_rgs)[ end(qry_aln_rgs) > cover_end_i ] <- cover_end_i
  qry_aln_rgs <- pmapToAlignments(
    qry_aln_rgs,
    setNames(
      aln,
      rep( "x", length(aln) )
    )
  )
  no_ins_pws <- pws_i [gaps(
    cigarRangesAlongQuerySpace( cigar(aln), ops=c("S","H","I") ),
    start(qry_aln_rgs), end(qry_aln_rgs)
  )]
  
  # consensus
  return(
    consensus_NON_insertions(aln_no_ins, no_ins_pws, model)
  )
}





#' Conpute consensus sequences using machine learning
#'
#' @param fa_hls_i 
#' @param pws_i 
#' @param fa_hls_files_i 
#' @param ovl_files_i 
#' @param cmds_ovl_i 
#' @param fa_ref_files_i 
#' @param cmds_aln_i 
#' @param aln_files_i 
#' @param temp_rootDirs_i 
#' @param min_coverage 
#' @param insert2 
#' @param consensus_insertions 
#' @param insertion_consensus_per_coverage 
#' @param svm_consensus_NON_insertions 
#' @param non_insertion_consensus_per_coverage 
#'
#' @return
#' @export
#'
#' @examples
consensus_using_ml <- function(fa_hls_i,
                               pws_i,
                               fa_hls_files_i,
                               ovl_files_i,
                               cmds_ovl_i,
                               fa_ref_files_i,
                               cmds_aln_i,
                               aln_files_i,
                               temp_rootDirs_i,
                               min_coverage,
                               insert2,
                               consensus_insertions,
                               insertion_consensus_per_coverage,
                               svm_consensus_NON_insertions,
                               non_insertion_consensus_per_coverage){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  
  fa_hls_i <- compact(fa_hls_i) ########### <<<<<<<<-------------====
  
  # find subread of reference
  fa_ref <- find_subread_of_reference(temp_rootDirs_i, fa_hls_i, fa_hls_files_i,
                                      cmds_ovl_i, ovl_files_i)
  
  if( is.null(fa_ref) ){
    unlink(temp_rootDirs_i, recursive=TRUE)
    return( sub(".+/([0-9]+)/.+", "no_ref__\\1", names(fa_hls_i[1])) )
  }
  
  # write the subread of reference as a fasta file
  writeXStringSet(fa_ref, fa_ref_files_i)
  
  # align all subreads to the subread of reference
  system(cmds_aln_i)
  aln <- readGAlignments( aln_files_i, param=ScanBamParam(what=c("qname", "seq")) )
  
  
  
  # If before the iteration only one sequence could be aligned, that means that it is the reference itself. So, there are no consensus to do.
  if( length(aln) == 1 ){
    
    unlink(temp_rootDirs_i, recursive=TRUE)
    return( sub(".+/([0-9]+)/.+", "no_aln__\\1", names(fa_ref)) )
  }
  
  
  # to find whether a hole is weird. we check whether the direction of the
  # subreads always alternate, considering that some subreads might not align
  # to the reference
  if( identify_weird_holes(aln, fa_hls_i) ){
    unlink(temp_rootDirs_i, recursive=TRUE)
    return( sub( ".+/([0-9]+)/.+", "no_only_normal_srs__\\1", names(fa_ref) ) )
  }
  
  
  
  ### for the training, use only regions that the coverage is `n`
  cover <- coverage(aln)[[1]]
  
  # find `n`-covered regions
  n_cover <- cover >= min_coverage
  n_cover_val <- runValue(n_cover)
  
  # it might doesn't contain any `n`-coverage region at all
  if( any(n_cover_val) ){
    # if there are more than one `n`-covered region, take the longest one
    if( length(n_cover_val[n_cover_val]) > 1 ){
      n_cover_len <- runLength(n_cover)
      max_len <- max(n_cover_len[n_cover_val])
      runValue(n_cover) [n_cover_val] [n_cover_len[n_cover_val] != max_len] <- FALSE
    }
    # the region where it will be computed the consensus sequence
    n_cover <- IRanges(n_cover)
    k <- (width(n_cover) / width(fa_ref)) >= min_percentage_high_coverage_area
  }else{
    k <- FALSE
  }
  
  # if there is a `n`-coverage region and it's long enough
  if(k){
    # if a subread is outside the `n`-coverege region, eliminate it
    qry_aln_rgs <- ranges(aln)
    sr_ovl_n_cover <- queryHits( findOverlaps( qry_aln_rgs, n_cover ) )
    
    # stopifnot( length(sr_ovl_n_cover) == length(aln) ) #### <<<<<<<<----------===
    
    aln <- aln[sr_ovl_n_cover]
    qry_aln_rgs <- qry_aln_rgs[sr_ovl_n_cover]
    
    # all insertions
    rrs <- cigarRangesAlongReferenceSpace( cigar(aln), ops="I", pos=start(aln) )
    
    # test whether there are insertions inside the high coverage region
    ins_pos <- start( unlist(rrs) )
    is_in_n_cover <- (ins_pos > start(n_cover)) & (ins_pos <= end(n_cover))
    
    # correct pulse width sequences of subreads in reverse strand
    pws_i <- pws_i[ mcols(aln)$qname ]
    k <- as.vector(strand(aln) == "-")
    pws_i[k] <- lapply( pws_i[k], rev )
    
    # sort pulse widths of consecutive repeated nucleotides (for minimap2, insertions to the left)
    pws_i <- RleList(
      only_sort_pws(pws_i, mcols(aln)$seq)
    )
    
    # divide the alignment according the coverage
    n_cover2 <- cover[n_cover]
    cover_end <- cumsum( runLength(n_cover2) ) + start(n_cover) -1
    cover_start <- cover_end - runLength(n_cover2) +1
    
    # consensus, but not about insertions
    invisible(
      indexBam(
        sortBam(
          aln_files_i,
          sub(".bam", "", aln_files_i)
        )
      )
    )
    aln_no_ins <- stackStringsFromBam(
      BamFile(aln_files_i),
      param=GRanges( names(fa_ref), IRanges(1, width(fa_ref)) ),
      use.names=TRUE
    ) [mcols(aln)$qname]
    
    # consensus of insertions separeted per coverage
    ins_tbl <- do.call(
      c,
      mapply(
        insertion_consensus_per_coverage,
        cover_start,
        cover_end,
        runValue(n_cover2),
        models[ runValue(n_cover2) -2 ],
        MoreArgs= list(ins_pos=ins_pos,
                       cover=cover,
                       aln=aln,
                       pws_i=pws_i,
                       temp_rootDirs_i=temp_rootDirs_i),
        SIMPLIFY=FALSE
      )
    )
    
    # consensus of non insertions separeted per coverage
    cons_seq <- paste(
      mapply(
        non_insertion_consensus_per_coverage,
        cover_start_i= cover_start,
        cover_end_i= cover_end,
        model= models[ runValue(n_cover2) -2 ],
        MoreArgs=list(aln_no_ins= aln_no_ins,
                      pws_i= pws_i,
                      aln= aln,
                      qry_aln_rgs= qry_aln_rgs)
      ),
      collapse=""
    )
    
    
    if( length(ins_tbl) > 0 ){
      cons_seq <- insert2(
        cons_seq,
        as.integer(names(ins_tbl)) - start(n_cover) +1,
        ins_tbl
      )
    }
    cons_seq <- setNames(
      gsub("[-+]", "", cons_seq),
      sub( "[^/]+$", "ccs", names(fa_ref) )
    )
  }else{
    unlink(temp_rootDirs_i, recursive=TRUE)
    return( sub( ".+/([0-9]+)/.+", "no_cover__\\1", names(fa_ref) ) )
  }
  unlink(temp_rootDirs_i, recursive=TRUE)
  return(cons_seq)
}










#' Consensus using machine learning
#' 
#' This is the main function to compute consensus. It calls many of the other functions.
#'
#' @param raw_srs_bam_file 
#' @param raw_srs_fa_file 
#' @param subread_lengths_median 
#' @param min_subread_number_per_hole 
#' @param min_obs_training 
#' @param min_percentage_high_coverage_area 
#' @param mm2 
#' @param num_of_top_holes 
#' @param THREADS 
#' @param min_coverage 
#' @param model 
#' @param max_pw 
#' @param gap_weight_ins 
#' @param gap_weight_del 
#'
#' @return
#' @export
#'
#' @examples
consensus_using_machine_learning <- function(raw_srs_bam_file,
                                             raw_srs_fa_file,
                                             subread_lengths_median,
                                             min_subread_number_per_hole,
                                             min_obs_training,
                                             min_percentage_high_coverage_area,
                                             mm2,
                                             num_of_top_holes,
                                             THREADS,
                                             min_coverage,
                                             model,
                                             max_pw,
                                             gap_weight_ins,
                                             gap_weight_del){
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  stopifnot( model %in% c("svm", "rf") )
  
  #### LOAD PACKAGES ####
  
  require(R.utils)
  require(msa)
  require(Rsamtools)
  require(GenomicAlignments)
  require(e1071)
  require(randomForest)
  require(foreach)
  require(iterators)
  require(doParallel)
  
  
  
  
  
  
  # register parallel backend
  registerDoParallel(cores=THREADS)
  
  
  
  
  
  #### SELECT HOLES FOR THE TRAINNING STEP ####
  
  ### load data
  dat <- load_pacbio_data(raw_srs_bam_file, raw_srs_fa_file)
  fa_hls <- dat$seq
  fa_hls <- split(
    fa_hls,
    sub(".*/([0-9]+)/.*", "\\1", names(fa_hls))
  )
  
  
  
  
  
  # ### only holes that the number of subreads is equal to or between subread_number_per_hole values
  # fa_hls <- fa_hls[lengths(fa_hls) >= subread_number_per_hole[1] & lengths(fa_hls) <= subread_number_per_hole[2]]
  ### only holes that contain `subread_number_per_hole` or more subreads
  fa_hls <- fa_hls[lengths(fa_hls) >= min_subread_number_per_hole]
  
  
  ### Only holes that the median of the subread lengths is between `subread_lengths_median[1]`
  ### and `subread_lengths_median[2]`.
  k <- sapply( fa_hls, function(x) median(width(x)) )
  fa_hls <- fa_hls[ (k >= subread_lengths_median[1]) & (k <= subread_lengths_median[2]) ]
  
  
  ### take to top holes that contain more subreads
  if( length(fa_hls) < num_of_top_holes ) num_of_top_holes <- length(fa_hls)
  fa_hls <- compact( fa_hls[ sample(length(fa_hls), num_of_top_holes) ] )
  k <- dat$zm %in% names(fa_hls)
  pws <- setNames(dat$pw, dat$qname) [k]
  pws <- split(pws, dat$z[k]) [names(fa_hls)]
  stopifnot( identical(names(fa_hls), names(pws)) )
  stopifnot( identical(lengths(fa_hls), lengths(pws)) )
  
  
  
  fa_hls <- compact(fa_hls)
  
  
  
  ### temp dir
  temp_rootDir <- tempdir()
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### eliminate the first and the last subread as they potentially don't cover the whole truth
  fa_hls <- compact(
    lapply( fa_hls, function(u) u[ -c(1, length(u)) ] )
  )
  fa_hls <- compact(fa_hls)
  pws <- lapply( pws, function(u) u[ -c(1, length(u)) ] )
  
  ### temporary directory
  temp_rootDirs <- tempfile( pattern=as.character(seq_along(pws)) )
  
  ### to save fasta files of holes
  fa_hls_files <- file.path(temp_rootDirs, "hl.fasta")
  
  ### to save fasta files of subreads of reference
  fa_ref_files <- file.path(temp_rootDirs, "ref.fasta")
  
  ### to overlap subreads
  ovl_files <- file.path(temp_rootDirs, "ovl.bam")
  cmds_ovl <- gettextf("%s -a -x ava-pb -n 0 -m 1 -s 1 -O 5 -E 6 -A 8 -B 8 %s %s | samtools view -bSh -F 2052 - > %s",
                       mm2,
                       fa_hls_files,
                       fa_hls_files,
                       ovl_files)
  sbp_ovl <- ScanBamParam( what=c("rname", "qname"), tag="AS" )
  
  ### to align to the subread of reference
  aln_files <- file.path(temp_rootDirs, "aln.bam")
  cmds_aln <- gettextf("%s -a -x map-pb -w 5 --no-long-join -O 5 -E 6 -A 8 -B 8 --secondary=no %s %s | samtools view -bSh -F 2308 - > %s",
                       mm2,
                       fa_ref_files,
                       fa_hls_files,
                       aln_files)
  
  
  
  
  
  
  # iterators?????? #################### <<<<<<<<<------===
  ### for each hole, get the svm input variables to train when using 3:length(fa_hls)
  ### subreads
  x_training <- foreach(fa_hls_i= fa_hls,
                        pws_i= pws,
                        fa_hls_files_i= fa_hls_files,
                        ovl_files_i= ovl_files,
                        cmds_ovl_i= cmds_ovl,
                        fa_ref_files_i= fa_ref_files,
                        cmds_aln_i= cmds_aln,
                        aln_files_i= aln_files,
                        temp_rootDirs_i= temp_rootDirs) %dopar%
    variables_for_the_training(fa_hls_i,
                               pws_i,
                               fa_hls_files_i,
                               ovl_files_i,
                               cmds_ovl_i,
                               fa_ref_files_i,
                               cmds_aln_i,
                               aln_files_i,
                               temp_rootDirs_i,
                               transform_pws,
                               only_sort_pws,
                               find_subread_of_reference,
                               truth_and_training_matrix_for_insertions,
                               truth_and_training_matrix_for_NON_insertions,
                               find_the_truth)
  
  
  
  
  
  # for each hole the consensus sequences is computed using total_subread_number - 2
  # subreads, but the training is done using, at the most, minimum subreads aligned
  # to reference per holes. Also, join holes data.
  x_training <- unlist(x_training, recursive=FALSE)
  
  k <- names(x_training) == "isIns"
  isIns <- x_training[k]
  x_training <- x_training[!k]
  isIns <- do.call("c", unname(isIns))
  
  k <- names(x_training) == "truth"
  truth <- x_training[k]
  x_training <- x_training[!k]
  truth <- do.call("c", unname(truth))
  
  x_training <- split(x_training, names(x_training))
  k <- sapply(x_training$nts, ncol)
  N <- min(k)
  k <- lapply( k, sample, size=N )
  
  x_training$nts <- do.call(
    "rbind",
    mapply(function(x,j){
      x[,j]
    },x_training$nts, k, SIMPLIFY=FALSE)
  )
  
  x_training$pws <- do.call(
    "rbind",
    x_training$pws <- mapply(function(x,j){
      x[,j]
    },x_training$pws, k, SIMPLIFY=FALSE)
  )
  
  
  
  # number of samples used for the training step
  nrow( x_training[[1]] )
  
  
  
  # get the models
  models <- train_the_models(3:N, x_training, truth, isIns, THREADS, model)
  
  
  
  
  
  #### COMPUTE CONSENSUS SEQUENCE ####
  
  ### load data
  fa_hls <- dat$seq
  fa_hls <- split(
    fa_hls,
    sub(".*/([0-9]+)/.*", "\\1", names(fa_hls))
  )
  
  ### take holes that contain `N` subreads at the most
  fa_hls <- compact( fa_hls[lengths(fa_hls) <= N] )
  
  k <- dat$zm %in% names(fa_hls)
  pws <- setNames(dat$pw, dat$qname) [k]
  pws <- split(pws, dat$zm[k])
  pws <- pws[ names(fa_hls) ]
  stopifnot( identical(names(fa_hls), names(pws)) )
  stopifnot( identical(length(fa_hls), length(pws)) )
  
  
  
  
  
  
  ### temporary directory
  temp_rootDirs <- tempfile( pattern=as.character(seq_along(pws)) )
  
  ### to save fasta files of holes
  fa_hls_files <- file.path(temp_rootDirs, "hl.fasta")
  
  ### to save fasta files of subreads of reference
  fa_ref_files <- file.path(temp_rootDirs, "ref.fasta")
  
  ### to overlap subreads
  ovl_files <- file.path(temp_rootDirs, "ovl.bam")
  cmds_ovl <- gettextf("%s -a -x ava-pb -n 0 -m 1 -s 1 -O 5 -E 4 -A 9 -B 5 %s %s | samtools view -bSh -F 2052 - > %s",
                       mm2,
                       fa_hls_files,
                       fa_hls_files,
                       ovl_files)
  sbp_ovl <- ScanBamParam( what=c("rname", "qname"), tag="AS" )
  
  ### to align to the subread of reference
  aln_files <- file.path(temp_rootDirs, "aln.bam")
  cmds_aln <- gettextf("%s -a -x map-pb -w 5 --no-long-join -O 5 -E 4 -A 9 -B 5 --secondary=no %s %s | samtools view -bSh -F 2308 - > %s",
                       mm2,
                       fa_ref_files,
                       fa_hls_files,
                       aln_files)
  
  
  
  
  
  
  ### compute consensus using machine learning
  cons_seqs <- foreach(fa_hls_i= fa_hls,
                       pws_i= pws,
                       fa_hls_files_i= fa_hls_files,
                       ovl_files_i= ovl_files,
                       cmds_ovl_i= cmds_ovl,
                       fa_ref_files_i= fa_ref_files,
                       cmds_aln_i= cmds_aln,
                       aln_files_i= aln_files,
                       temp_rootDirs_i= temp_rootDirs,
                       .combine="c") %dopar%
    consensus_using_ml(fa_hls_i,
                       pws_i,
                       fa_hls_files_i,
                       ovl_files_i,
                       cmds_ovl_i,
                       fa_ref_files_i,
                       cmds_aln_i,
                       aln_files_i,
                       temp_rootDirs_i,
                       min_coverage,
                       insert2,
                       consensus_insertions,
                       insertion_consensus_per_coverage,
                       svm_consensus_NON_insertions,
                       non_insertion_consensus_per_coverage)
  
  return(cons_seqs)
  
}






