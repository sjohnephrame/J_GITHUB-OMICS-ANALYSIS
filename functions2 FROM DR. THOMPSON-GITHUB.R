# library(gageData)
library(org.Mm.eg.db)     # contains symbols
library(annotate)
library(org.Hs.eg.db)     # contains symbols
#library(AnnotationFuncs)  # for translation of symbols

#' Translate between different identifiers
#'
#' Function for translating from one annotation to another, eg. from RefSeq to  
#' Ensemble.  This function takes a vector of annotation values and translates     
#' first to the primary annotation in the Biocore Data Team package (ie. entrez gene identifier for org.Bt.eg.db)
#' and then to the desired product, while removing non-translated annotations   
#' and optionally reducing the result so there is only a one-to-one relation. 
#'
#' If you want to do some further mapping on the result, you will have to use
#' either \code{unlist} or \code{lapply}, where the first returns all the end-products
#' of the first mapping, returning a new list, and the latter produces a list-within-list.
#'
#' If \code{from} returns GO identifiers (e.g. \code{from = org.Bt.egGO}), then the 
#' returned resultset is more complex and consists of several layers of lists instead of 
#' the usual list of character vectors. If \code{to} has also been specified, the GO IDs 
#' must be extracted (internally) and you have the option of filtering for evidence and category at this point.
#' See \code{\link{pickGO}}.
#'
#' @note Requires user to deliver the annotation packages such as org.Bt.egREFSEQ.   
#' @param values Vector of annotations that needs translation. Coerced to character vector. 
#' @param from Type of annotation \code{values} are given in. NB! take care in the       
#'         orientation of the package, ie. if you have RefSeq annotations, use  
#'         \code{org.Bt.egREFSEQ2EG} or (in some cases) \code{revmap(org.Bt.egREFSEQ)}.
#' @param to Desired goal, eg. \code{org.Bt.egENSEMBLPROT}. If \code{NULL} (default), goal 
#'         if the packages primary annotation (eg. entrez gene for org.Bt.eg.db).
#'         Throws a warning if the organisms in \code{from} and \code{to} are not the same.
#' @param reduce Reducing method, either return all annotations (one-to-many relation)
#'         or the first or last found annotation. The reducing step is applied  
#'         after translating to the goal:                                       
#'         \code{all}: returns all annotations                                       
#'         \code{first} or \code{last}: choose first or last of arbitrarily ordered list.   
#' @param return.list Logical, when \code{TRUE}, returns the translation as a list where names  
#         are the given values that could be translated. The function works    
#         on lists, so this is just as easy to return.
#         When \code{FALSE}, a table.  Convenience for calling on the result.             
#' @param remove.missing Logical, whether to remove non-translated values, defaults \code{TRUE}.
#' @param simplify Logical, unlists the result. Defaults to FALSE. Usefull when using \code{translate} in 
#'                 a \code{lapply} or \code{sapply}.
#' @param ... Additional arguments sent to \code{\link{pickGO}} if \code{from} returns GO set.
#' @return List; names of elements are \code{values} and the elements are the translated elements,
#'        or \code{NULL} if not translatable with \code{remove.missing = TRUE}.
#' @author Stefan McKinnon Edwards \email{stefan.hoj-edwards@@agrsci.dk}
#' @seealso \code{\link{pickRefSeq}}, \code{\link{pickGO}}
#' @export
#' @examples
#' library(org.Bt.eg.db)
#' genes <- c(280705, 280706, 100327208)
#' translate(genes, org.Bt.egSYMBOL)
#'
#' symbols <- c("SERPINA1","KERA","CD5")
#' refseq <- translate(symbols, from=org.Bt.egSYMBOL2EG, to=org.Bt.egREFSEQ)
#' # Pick the proteins:
#' pickRefSeq(refseq, priorities=c('NP','XP'), reduce='all')
#'
#' # If you wanted do do some further mapping on the result from 
#' # translate, simply use lapply.
#' 
#' library(GO.db)
#' GO <- translate(genes, org.Bt.egGO)
#' # Get all biological processes:
#' pickGO(GO, category='BP')
#' # Get all ontologies with experimental evidence:
#' pickGO(GO, evidence=c('IMP','IGI','IPI','ISS','IDA','IEP','IEA'))
translate <- function(values, from, to=NULL, 
                      reduce=c('all','first','last'), 
                      return.list = TRUE,
                      remove.missing=TRUE, 
                      simplify=FALSE,
                      ...) {
  # Roadmap:
  # Check validity of attributes. 
  # Check that the organisms of from and to match.
  # Translate to primary annotation.
  # Optionally translate to end.
  # Optionally reduce result.
  # Optionally unstack.
  
  ###  
  # Check validity of attributes
  ###
  values <- unique(as.vector(values, mode='character'))
  values <- values[!is.na(values)]
  values <- values[which(sapply(values, nchar)>0)]  # Removes zero-length entries.
  reduce <- match.arg(reduce)
  return.list <- as.logical(return.list)
  
  #l <- NULL # R compiler complains about l not being binded? Half a year later and I am still not sure what this does.
  
  ###  
  # Check organisms - throw warning if not.
  ###
  if (!is.null(to)) {
    # Compare taxid:
    org.from <- NULL
    org.to <- NULL
    try( org.from <- dbmeta(dbconn(from), "ORGANISM"), silent=TRUE)
    try( org.to <- dbmeta(dbconn(to), "ORGANISM"), silent=TRUE)
    if (!is.null(org.from) & !is.null(org.to)) {
      if (org.from != org.to) {
        warning(sprintf("TAXID for 'to' and 'from' does not match! We found:
\t From:\t%s
\tTo:\t%s
\tFor cross species (ie. homologes) look at e.g. hom.Hs.inp.db.", org.from, org.to))
        # End of warning message.
      }
    }
  }
  
  ###
  # Translate to primary id:
  ###
  primary <- AnnotationDbi::mget(values, from, ifnotfound=NA)
  # primary is now a list.
  # Remove nulls and NAs
  primary <- primary[!is.na(primary)]
  # primary is now a list with non-empty, multiple length entries
  
  # Provide a method for reducing.
  pickOne <- switch(reduce,
                    all = function(x) x, 
                    first = function(x) x[1], 
                    last = function(x) x[length(x)])
  
  ###
  # Translate all primary ids to goal annotation 
  ###
  if (is.null(to) == FALSE) {
    # Check for special case where the from-package is into GO.
    # If it is, then the primary must be reduced to only the GO identifiers and not the set of GO ID, evidence and category.
    isGO <- FALSE
    try(isGO <- from@objName == 'GO', silent=TRUE)
    if (isGO) {
      primary <- pickGO(primary, ...)
    }  
    # Remember to unlist, so we get every entry of primary, which is every 
    # possible translation of the starting values.
    goal <- AnnotationDbi::mget(unlist(primary, use.names=F), to, ifnotfound=NA)
    # Map all goal-annotations to the starting annotations:
    # lapply over primary, so we get primary's names.
    goal <- lapply(primary, function(x) {
      x <- unlist(goal[x], use.names=FALSE) # x is vector, so we get all entries in goal for all x.
      if (length(x) == 0 || is.na(x[1])) 
        return(NA)
      x <- pickOne(x)
      return(x)
    })
  } else {
    # If goal is primary annotation, reduce accordingly. 
    goal <- primary
    if (reduce != 'all')
      goal <- lapply(goal, pickOne)     
  }
  
  ###
  # Remove nulls and NAs
  ###
  if (remove.missing) {
    goal <- goal[!is.na(goal)]
  } else {
    missing <- values[!(values %in% names(goal))]
    #missing.list <- as.list(rep(NA, length(missing)))
    #names(missing.list) <- missing
    #goal <- c(goal, missing.list)
    goal[missing] <- NA
  }
  
  ###
  # Return result if a list is desired, 
  # else unstack to a table.
  ###
  if (simplify) return(unlist(goal, use.names = FALSE))
  if (return.list) {
    return(goal)
  } else {
    if (length(goal) == 0)
      return(data.frame(from = factor(), to = factor()))
    goal <- stack(goal)
    new.c <- c('from', 'to') # Matching 'ind' and 'values'.
    colnames(goal) <- new.c[match(c('ind', 'values'), colnames(goal))]
    return(goal)
  }                      
}


#' Load a library while suppressing all messages
#' 
#' @param lib (character): name of the library
#' 
#' @return list of attached packages
load_lib <- function(lib) {
	suppressWarnings(
		suppressMessages(
			suppressPackageStartupMessages(library(lib, character.only = T))))
}


#' Create a data.frame containing significant PTMs and whether they have sig
#' differential expression in at least one dataset
#' 
#' @param sig_ptms (data.frame): Accession column and pvalue column
#' @param d1tt (data.frame): should have an Accession column 
#' @param d2tt (data.frame): should have an Accession column 
#' @param pval_col (character): name of column containing p-values
#' 
#' @return data.frame with PTM accessions, there significance, and whether
#' they are differentially expressed in at least one dataset.
#' 
annotate_ptms <- function(sig_ptms, d1tt, d2tt, pval_col) {
  sig_ptms <- sig_ptms[sig_ptms$pvalue < 0.05,]
  sig_status <- check_for_one_sig(sig_ptms, d1tt, d2tt, pval_col)
	one_sigs <- data.frame(Accession = names(sig_status), Sig = sig_status)
	sig_ptms <- sig_ptms %>% left_join(one_sigs)
	return(na.omit(sig_ptms))
}

#' Check if at least one of the results is significant for two data.frames
#' given a specific accession.
#' 
#' @param d1tt (data.frame): should have an Accession column 
#' @param d2tt (data.frame): should have an Accession column 
#' @param acc (character): accession of the protein to search for
#' @param pval_col (character): name of column containing p-values
#' @param threshold (numeric): threshold for significance
#' 
#' @return Boolean value
#' 
at_least_one_sig <- function(d1tt, d2tt, acc, pval_col, threshold = 0.05) {
	return(min(d1tt[d1tt$Accession == acc, pval_col], 
						 d2tt[d2tt$Accession == acc, pval_col], na.rm = T) < threshold)
}

#' Given a data.frame of PTMs and their significance, check for each if it is 
#' differentially expressed in at least one dataset. For example, this could
#' be the treatment vs control for the protein and the treatment vs. control
#' for the PTM. This is because a PTM could be significantly different compared
#' to the protein as a whole, without either being individually significant. 
#' 
#' @param sig_ptms (data.frame): Accession column and pvalue column
#' @param d1tt (data.frame): should have an Accession col and a P.Value column
#' @param d2tt (data.frame): should have an Accession col and a P.Value column
#' @param pval_col (character): name of column containing p-values
#' @param threshold (numeric): threshold for significance
#' 
#' @return vector of Boolean values
#' 
check_for_one_sig <- function(sig_ptms, d1tt, d2tt, pval_col, threshold = .05) {
	res <- sapply(sig_ptms[!is.na(sig_ptms$pvalue) & sig_ptms$pvalue < threshold,]$Accession,
								function(x) at_least_one_sig(d1tt, d2tt, x, pval_col))
	return(res)
}


#' Perform overrepresentation analysis
#' 
#' @param react (list): named list of gene sets for pathways
#' @param ptms_sig (character): vector of ptm symbols as shown in gene sets
#' @param background (character): vector of ptm's that could have been found
#' @param min_size (integer): min number of proteins in a pathway
#' 
#' @return table with results
#' 
get_enrichment <- function(react, ptms_sig, background, min_size = 5) { 
	require(fgsea)
	
	res <- fora(react, ptms_sig, background, minSize = min_size)
	
	return(res)
}

#' Determine which PTMs have changed significantly compared to the overall
#' change for that protein compared to the control
#' 
#' @param full_dat (data.frame): protein expression data, with proteins as rows
#' @param ptm_dat (data.frame): ptm expression data, with proteins as rows
#' @param n1 (integer): number of samples in full_dat
#' @param n2 (integer): number of samples in ptm_dat
#' @param method (character): type of t-test to perform, either "welch" or "student"
#' @param logFC_col (character): name of column in topTable results containing logFC values
#' 
#' @return data.frame of protein accessions, gene symbols, and pvalues
#' 
sig_ptms <- function(full_dat, ptm_dat, n1, n2, method, logFC_col) {
	common <- intersect(full_dat$tt$Accession, ptm_dat$tt$Accession)
	accs <- ptm_dat$tt$Accession[ptm_dat$tt$Accession %in% common]
	pvalues <- lapply(accs, 
										function(x) lfc_ttest(ptm_dat$mod, 
																					full_dat$mod, 
																					ptm_dat$tt, 
																					full_dat$tt, 
																					x,
																					n1,
																					n2,
																					method,
																					logFC_col))
	pvalues <- do.call(rbind, pvalues)
	colnames(pvalues) <- c('Accession', 'pvalue', 'tstat') # colnames(pvalues) <- c('Accession', 'pvalue')
	temp <- ptm_dat$tt[,c('Accession', 'Symbol')]
	pvalues <- pvalues %>% left_join(temp)
	return(pvalues)
}

#' Perform Student's t-test or Welch's t-test of two log-fold changes.
#' 
#' @param m1 (limma): limma eBayes model
#' @param m2 (limma): limma eBayes model
#' @param d1tt (data.frame): topTable results 
#' @param d2tt (data.frame): topTable results
#' @param acc (character): Accession number to test
#' @param n1 (integer): number of samples in first dataset
#' @param n2 (integer): number of samples in second dataset
#' @param method (character): type of t-test to perform, either "welch" or "student"
#' @param logFC_col (character): name of column in topTable results containing logFC values
#' 
#' @return named vector of log fold changes
#'
lfc_ttest = function(m1, m2, d1tt, d2tt, acc, n1, n2 , method = c( "welch", "student"), logFC_col) {
  method = match.arg(method)
  
	d1_sd = m1[acc, ]$sigma
	d2_sd = m2[acc, ]$sigma
	
	d1_lfc = d1tt[d1tt$Accession == acc, logFC_col] # $logFC
	d2_lfc = d2tt[d2tt$Accession == acc, logFC_col] # $logFC
	
	# Calculate the standard error
	if(method == "student"){
	  numer_s2 =  (n1 - 1) * (d1_sd)^2 + (n2 - 1) * (d2_sd)^2
	  denom_s2 = (n1 + n2 - 2)
	  s2 = numer_s2/denom_s2
	  s = sqrt(s2)
	  sem = s * sqrt(1/n1 + 1/n2)
	  dof = n1 + n2 - 2
	}else if(method == "welch"){
	  sem = sqrt((d1_sd)^2/n1 + (d2_sd)^2/n2)
	  dof = (d1_sd^2/n1 + d2_sd^2/n2)^2 / ((d1_sd^2/n1)^2/(n1 - 1) + (d2_sd^2/n2)^2/(n2 - 1))
	}else stop("argument \'method\' requires either \'student\' or \'welch\'")
  
	# tstar = (d1_lfc - d2_lfc) / sem
	tstar = abs(d1_lfc - d2_lfc) / sem
	
	# return(data.frame(Accession = acc, pvalue = 2*pt(abs(tstar), dof, lower = FALSE)))
	return(data.frame(Accession = acc, 
	                  pvalue = 2*pt(tstar, dof, lower = FALSE),
	                  tstat = tstar))
} # end lfc_ttest


#' Adds the gene symbols from one dataframe onto another
#' 
#' @param tt (data.frame): topTable results from limma 
#' @param dat (data.frame): a data.frame with an Accession and a Symbol column
#' 
#' @return data.frame containing topTable results with a gene symbol column
#'
add_annotation <- function(tt, dat) {
	require(tidyverse)
	
	tt$Accession <- rownames(tt)
	tt <- tt %>% left_join(dat[,c('Accession', 'Symbol')],
												 by = c('Accession' = 'Accession'))
	return(tt)
} # end add_annotation


#' Builds a limma model and gets the topTable results, with gene symbols added
#' 
#' @param dat (data.frame): a data.frame with an Accession and a Symbol column
#' @param dat_mat (matrix): a data.frame with an genes as rows
#' @param design (matrix): a design matrix for limma
#' @param cnt_mat (matrix): a contrast matrix for limma
#' @param annotate (logical): flag indicating if data should be annotated
#' @param trend (logical): flag indicating if mean-variance trend should be accounted for in eBayes (limma-trend method)
#' 
#' @return list with the topTable results and the model
#'
top_table_results <- function(dat, dat_mat, design, cnt_mat, annotate = F, robust = F, trend = F) {
  ###
  if(0){ # For debugging
    dat = mapping; dat_mat = dat[[1]]; design = DM[[1]]; cnt_mat = cnt_mat[[1]]; annotate = FALSE; trend = TRUE
  }
  ###
  require(limma)
  
  mod <- get_model(dat_mat, design, cnt_mat, robust = robust, trend = trend)
  tt <- topTable(mod, number = 100000, adjust = 'BH') # topTable - Extract a table of the top-ranked genes from a linear model fit.
  tt = tt[,!colnames(tt) %in% c("F", "P.Value", "adj.P.Val")]
  
  cnt_mat_names = colnames(cnt_mat)
  p.vals = vector(mode = "list", length = length(cnt_mat_names))
  for(i in 1:ncol(cnt_mat)){
    tt.temp = topTable(mod, coef = cnt_mat_names[i], number = 100000, adjust = 'BH')
    p.vals[[i]] = tt.temp[match(row.names(tt), row.names(tt.temp)), c("P.Value", "adj.P.Val")]
    colnames(p.vals[[i]]) = str_replace(paste(colnames(p.vals[[i]]), cnt_mat_names[i], sep = "."), "-", ".")
  }
  tt = cbind(tt, do.call(cbind, p.vals))
  
  if(annotate) {
    tt <- add_annotation(tt, dat)	
  } else {
    tt$Symbol = rownames(tt)
    tt$Accession = rownames(tt)
  }
  return(list(tt = tt, mod = mod))
} # end top_table_results


#' Perform PCA on data and visualized with ggplot2
#' 
#' @param dat (data.frame): with genes as rows, observations as columns
#' @param labels (character): vector of group labels
#' @param legend_title (character): title for legend
#' @param cnt_mat (matrix): a contrast matrix for limma
#' @param pep (Boolean): flag indicating if data contains peptides
#' 
#' @return list with the topTable results and the model
#'
get_pca <- function(dat, labels, legend_title = 'Treatment') {
	require(ggplot2)
	
	pca_res <- prcomp(t(dat), center = T, scale. = F)
	
	U <- data.frame(pca_res$x)
	
	p <- ggplot(data = U, aes(x = PC1, y = PC2, color = as.factor(labels))) +
		geom_point(size = 3, alpha = 0.5) + 
		theme_bw() +
		labs(color = legend_title)
	
	return(list(p = p , summ = summary(pca_res)))
} # end get_pca

#' Add HUGO gene symbols to a data.frame with genes and RefSeq gene ids.
#' 
#' Assumes that data.frame passed has an Acc column with RefSeq gene ids. 
#' 
#' @author Jeffrey A. Thompson
#' 
#' @param df (data.frame): data.frame with genes as rows
#' @param acc_col_name (character): name of column with accession number
#' 
#' @return data.frame with new Symbol column
#'
add_symbols <- function(df, acc_col_name = 'Accession') {
	accs <- as.matrix(df[,acc_col_name])[,1]
	
	# Many of the accessions will have a version, indicated by a dot and
	# a number at the end. We'll strip these off.
	dot_locs <- unlist(gregexpr(pattern = '\\.', accs))
	df[,acc_col_name] <- ifelse(dot_locs == -1, accs, substr(accs, 1, dot_locs - 1))
	
	accs <- as.matrix(df[,acc_col_name])[,1]
	
	# translate RefSeq accessions to gene symbols
	symbols <- unlist(translate(accs, org.Mm.egREFSEQ2EG, org.Mm.egSYMBOL))
	
	# Make symbols unique by numbering any duplicates, be aware that this
	# will break matching these symbols with anything else though. 
	dup_symbols <- list()
	for(i in which(duplicated(symbols))) {
		if(is.null(dup_symbols[[symbols[i]]])) {
			dup_symbols[[symbols[i]]] <- 1
			symbols[i] <- paste0(symbols[i], '.', dup_symbols[[symbols[i]]] + 1)
		} else {
			dup_symbols[[symbols[i]]] = dup_symbols[[symbols[i]]] + 1
			symbols[i] <- paste0(symbols[i], '.', dup_symbols[[symbols[i]]] + 1)
		}
	}
	
	# Create data.frame with symbols and accession and then join back to
	# original data.frame.
	symbols <- data.frame(Acc = names(symbols), Symbol = symbols)
	colnames(symbols)[1] <- acc_col_name
	df <- symbols %>% left_join(df)
	
	return(df)
}

#' Get geometric means.
#' 
#' Note that if there are missing values for all technical replicates but one
#' then just that one will be returned.
#' 
#' @author Jeffrey A. Thompson
#' 
#' @param dat (data.frame): data.frame with genes as rows
#' @param sample_ids (character): which sample ids to find the mean for
#' @param cols (integer): which columns have the numeric data
#' @param label_col (character): name of column with gene ids
#' @param rowname_length (integer): max length of a gene name, used to make unique row names
#' 
#' @return data.frame with geometric means
#'
get_geo_means <- function(dat, sample_ids, cols = 3:ncol(dat), label_col = 'Accession', rowname_length = 15) {
	unique_ids <- unique(sample_ids)
	
	new_dat <- matrix(0L, nrow = nrow(dat), ncol = length(unique_ids))
	
	for(i in 1:nrow(dat)) {
		for(j in 1:length(unique_ids)) {
			new_dat[i,j] <- geoMean(as.numeric(dat[i,cols][which(sample_ids == unique_ids[j])]), na.rm=T)
		}
		#message(i, appendLF = T)
	}
	
	new_dat <- data.frame(new_dat)
	colnames(new_dat) <- unique_ids
	rownames(new_dat) <- make.unique(substr(as.character(as.matrix((dat[,label_col]))), 1, rowname_length))
	
	return(new_dat)
}

#' This functions does all the steps to use limma to fit an eBayes model.
#' 
#' @param mat (matrix): matrix of gene expression w/ genes as rows
#' @param design: design matrix from model.matrix
#' @param cnt_mat: contrast matrix from MakeContrasts
#' @param trend (logical): flag indicating if mean-variance trend should be accounted for in eBayes (limma-trend method)
#' 
#' @return differential expression data.frame
#'
get_model = function(mat, design, cnt_mat, robust = F, trend = F) {
  require(limma)
  
  fit = lmFit(mat, design)
  fit_contrast = contrasts.fit(fit, cnt_mat)
  if(robust) {
    fit_ebayes = eBayes(fit_contrast, trend = T, robust = T)
  } else {
    fit_ebayes = eBayes(fit_contrast, trend = trend, robust = F)	
  }
  
  return(fit_ebayes)
}

#' Creates a named vector of log fold changes for fgsea.
#' 
#' @param df (data.frame): data.frame with diff expression data.
#' @param lfc_col (character): name of column with log fold changes
#' @param name_col (characeter): name of column with protein/gene names
#' 
#' @return named vector of log fold changes
#'
get_gsea_stat = function(df, lfc_col = 'logFC', name_col = 'Symbol') {
	stat = df[,lfc_col]
	stat = setNames(stat, df[,name_col])
	return(stat)
}

#' Creates a named vector of p-values for ECEA.
#' 
#' @param df (data.frame): data.frame with diff expression data.
#' @param pval_col (character): name of column with p-values
#' @param name_col (characeter): name of column with protein/gene names
#' 
#' @return named vector of p-values
#'
get_p_values = function(df, pval_col = 'P.Value', name_col = 'Symbol') {
	stat = df[,pval_col]
	stat = setNames(stat, df[,name_col])
	return(stat)
}

#' Nicely format the results from FGSEA.
#' 
#' 
#' @param df (data.frame): FGSEA results
#' @param cap (character): table caption
#' @param digits (integer): the number of digits of significance
#' @param rem_cols (vector): vector of columns to remove from table
#' 
#' @return named vector of log fold changes
#'
format_enrichment = function(df, cap = '', digits = 3, rem_cols = c()) {
	tbl = flextable(df)
	tbl = colformat_num(tbl, 
											j = setdiff(c('pval', 'padj', 'log2err', 'ES', 'NES'), 
																	rem_cols),
											digits = digits)
	tbl = set_caption(tbl, cap)
	tbl = fontsize(tbl, size = 9)
	tbl = width(tbl, j = c('pathway'), width = 4)
	return(tbl)
}

#' Nicely format the results from eBayes topTable.
#' 
#' 
#' @param df (data.frame): topTable results
#' @param cap (character): table caption
#' @param digits (integer): the number of digits of significance
#' @param rem_cols (vector): vector of columns to remove from table
#' 
#' @return named vector of log fold changes
#'
format_de = function(df, cap = '', digits = 3, rem_cols = c()) {
	tbl = flextable(df)
	tbl = colformat_num(tbl, 
											j = setdiff(c('logFC', 'AveExpr', 't', 'P.Value', 
																		'adj.P.Val', 'B'), 
																	rem_cols),
											digits = digits)
	tbl = set_caption(tbl, cap)
	tbl = fontsize(tbl, size = 9)
	return(tbl)
}

get_KEGG = function() {
  data("kegg.sets.mm")
  
  kegg.sets = list()
  for(i in 1:length(kegg.sets.mm)) {
    kegg.sets[[i]] = as.vector(na.omit(getSYMBOL(kegg.sets.mm[[i]], data='org.Mm.eg')))
    message(".", appendLF=FALSE)
  }
  
  names(kegg.sets) = names(kegg.sets.mm)
  
  return(kegg.sets)
}

get_GO = function() {
	data("go.sets.mm")
	data("go.subs.mm")
	
	pathways = go.sets.mm
	pathways = pathways[go.subs.mm$BP]
	
	pb = txtProgressBar(min = 0, max = length(pathways), style = 3)
	
	go.sets = list()
	for(i in 1:length(pathways)) {
		go.sets[[i]] = as.vector(na.omit(getSYMBOL(pathways[[i]], data='org.Mm.eg')))
		setTxtProgressBar(pb, i)
	}
	close(pb)
	names(go.sets) = names(pathways)
	
	return(go.sets)
}

getReactome <- function(species = 'human') {
	#' Gets gene sets with Hugo gene symbols to use with FGSEA analysis.
	#'
	#' @param species currently accepts 'human' or 'mouse'
	#'
	require(reactome.db)
	require(annotate)
	
	db = ''
	
	if(species == 'human') {
		require(org.Hs.eg.db)
		db <- 'org.Hs.eg'
	} else if(species == 'mouse') {
		require(org.Mm.eg.db)
		db <- 'org.Mm.eg'
	} else {
		stop(paste0('Species ', species, ' not supported.'))
	}
	
	reactome_sets_full <- as.list(reactomePATHID2EXTID)
	
	pb <- txtProgressBar(min = 0, max = length(reactome_sets_full), style = 3)
	
	reactome_sets <- list()
	
	for(i in 1:length(reactome_sets_full)) {
		reactome_sets[[i]] <- as.vector(na.omit(getSYMBOL(reactome_sets_full[[i]], data=db)))
		setTxtProgressBar(pb, i)
	}
	close(pb)
	names(reactome_sets) <- names(reactome_sets_full)
	
	xx = as.list(reactomePATHID2NAME)
	
	names(reactome_sets) = xx[names(reactome_sets)]
	
	reactome_sets <- reactome_sets[lapply(reactome_sets, length)>0]
	
	return(reactome_sets)
}
