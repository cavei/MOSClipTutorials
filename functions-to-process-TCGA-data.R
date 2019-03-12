#### Function 
imputeKNN <- function(...) {
  sink("/dev/null")
  data <- impute.knn(...)
  sink()
  data
}

selectMutation <- function(smallMaf, type){
  mutationSelection <- smallMaf$type %in% type
  # hugo <- smallMaf[[1]][mutationSelection]
  # entrez <- smallMaf[[2]][mutationSelection]
  # samples <- smallMaf[[16]][mutationSelection]
  # patients <- substr(samples, 1,12)
  # type <- smallMaf[[9]][mutationSelection]
  # data.frame(hugo=hugo, entrez=entrez, patient=patients, type=type, stringsAsFactors = FALSE)
  smallMaf[mutationSelection,]
}

selectImpact <- function(smallMaf, type){
  impactSelection <- smallMaf$impact %in% type
  # hugo <- smallMaf$hugo[impactSelection]
  # entrez <- smallMaf$entrez[[2]][impactSelection]
  # samples <- smallMaf$1[[16]][impactSelection]
  # patients <- substr(samples, 1,12)
  # type <- smallMaf[[9]][impactSelection]
  # impact <- smallMaf[[94]][impactSelection]
  # data.frame(hugo=hugo, entrez=entrez, patient=patients, type=type, impact=impact, stringsAsFactors = FALSE)
  smallMaf[impactSelection,]
}

summarizeMutation <- function(reducedMaf) {
  patient2gene2mutation <- tapply(1:NROW(reducedMaf), paste(reducedMaf$entrez, reducedMaf$patient, sep="_"),
                                  function(idx) {
                                    mat <- reducedMaf[idx, , drop=F]
                                    type <- paste(reducedMaf$type[idx], collapse=";")
                                    impact <- paste(reducedMaf$impact[idx], collapse=";")
                                    cbind(mat[1, c("hugo", "entrez")], impact, type, mat[1, c("patient")])
                                    })
  do.call(rbind,patient2gene2mutation)
}

set_expand <- function(df, cols, sep=";", pairs=TRUE){
  assert_class(df,classes = "DataFrame")
  
  if (is.numeric(cols)) {
    stop("Colnames needed")
    # if (max(cols) > NCOL(df))
    #   stop("Column index out of range")
  } else {
    unk <- setdiff(cols, colnames(df))
    if (length(unk))
      stop(paste0("Some columns were not found: ", paste(unk, collapse=", ")))
  }
  
  expandedCols <- lapply(cols, function(c) {
    strsplit(df[[c]], sep)
  })
  
  expColsList <- lapply(expandedCols, function(x) {
    do.call(c, x)
  })
  
  names(expColsList) <- cols
  
  col1list <- lapply(expandedCols[[1]], length)
  
  if (length(cols) > 1 & pairs) {
    a <- lapply(expandedCols[-c(1)], function(x) {
      collist <- lapply(x, length)
      if (!identical(collist, col1list))
        stop("Pair mode needs the same length when columns are expanded.")
      collist
    })
  }
  
  mCols <- setdiff(colnames(df), cols)
  mCOlList <- lapply(mCols, function(m) {
    rep(df[[m]], times=col1list)
  })
  names(mCOlList) <- mCols
  
  DataFrame(mCOlList, expColsList)[,colnames(df)]
}

assignCGToGenePromoter <- function(data, promoterUP=-2000, promoterDW=500) {
  discretePos <- rep("geneBody", length(data$Position_to_TSS))
  discretePos[data$Position_to_TSS >= promoterUP & data$Position_to_TSS <= promoterDW] <- "Promoter"
  discretePos
}

prepareMutations <- function(maf, impact=NULL, mutType=NULL, filterByThisEntrez=NULL, patients=NULL,
                             onlyPrimary=TRUE, select_columns=c(hugo=1,entrez=2, impact=94, type=9, patient=16)) {
  
  if (is.null(names(select_columns)))
    stop("you need names associated to the columns to select")
  
  if (!is.numeric(select_columns))
    stop("select_columns must be culumn index")
  
  if (NCOL(maf) < max(select_columns))
    stop("your imput does not appear to be a valid maf file")
  
  if (onlyPrimary) {
    bcode <- maf$Tumor_Sample_Barcode
    select <- substr(bcode, 14,15) == "01"
    maf <- maf[select, , drop=F]
  }
  
  allPatientsMeasured <- unique(extractTCGAPatientsName(maf[[16]]))
  smallMaf <- maf[, select_columns]
  colnames(smallMaf) <- names(select_columns)
  smallMaf$patient <- extractTCGAPatientsName(smallMaf$patient)
  
  if (!is.null(impact))
    smallMaf <- selectImpact(smallMaf, impact)
  
  if (!(is.null(mutType)))
    smallMaf <- selectMutation(maf, mutType)
  
  mut <- summarizeMutation(smallMaf)
  mut$entrez <- as.character(mut$entrez)
  
  if (!is.null(filterByThisEntrez))
    mut <- mut[mut$entrez %in% filterByThisEntrez, , drop=F]
  
  if (!is.null(patients)) {
    patients <- intersect(unique(mut$patient), patients)
    mut <- mut[mut$patient %in% patients, , drop=F]
  }
  entrez <- unique(mut$entrez)
  patients <- unique(mut$patient)
  mutations <- matrix(0, nrow=length(entrez), ncol=length(patients), dimnames = list(entrez, patients))
  for (i in seq_along(mut$entrez)) {
    x = as.character(mut$entrez[i])
    y = mut$patient[i]
    mutations[x,y] <- 1
  }
  if (length(setdiff(patients, allPatientsMeasured)) > 0)
    warning("some patients were excluded beacuse they do not have any of your favourite mutations. 
            Consider adding an 0 column to prevent patients loss.")
  list(data=mutations, allPatients=allPatientsMeasured, mutationsTypes=mut)
}

createEntrezMethylationMatrix <- function(methylationAssay) {
  require(org.Hs.eg.db)
  
  meta = rowData(methylationAssay)
  
  expandedMeta <- set_expand(meta, cols=c("Gene_Symbol", "Position_to_TSS"))
  expandedMeta$Position_to_TSS <- as.numeric(expandedMeta$Position_to_TSS)
  
  expandedMeta <- na.omit(expandedMeta)
  expandedMeta$discretePos <- assignCGToGenePromoter(expandedMeta)
  
  promoterMeta <- unique(expandedMeta[expandedMeta$discretePos=="Promoter",c("Composite.Element.REF", "Gene_Symbol")])
  
  # symbol2entrez <- mapIds(org.Hs.eg.db, keys=promoterMeta$Gene_Symbol, column="ENTREZID", keytype="SYMBOL", multiVals="list")
  symbol2entrez <- mapIds(org.Hs.eg.db, keys=promoterMeta$Gene_Symbol, column="ENTREZID", keytype="SYMBOL")
  
  promoterMeta$entrez <- symbol2entrez
  promoterMeta <- na.omit(promoterMeta)
  
  methylationPromoter <- methylationAssay[promoterMeta$Composite.Element.REF] ## filter with valid CG
  assayPromoter <- methylationPromoter@assays[[1]]
  
  summarized <- tapply(row.names(assayPromoter), promoterMeta$entrez, function(cg){
    colMeans(assayPromoter[cg, ,drop=F], na.rm = T)
  })
  
  summarized <- do.call(rbind, summarized)
  colnames(summarized) <- substr(colnames(summarized), 1, 12)
  keep <- apply(summarized, 1, function(x) { !all(is.na(x) )})
  list(assay=summarized[keep, , drop=F], annotation=promoterMeta)
}
createDiscreteVersionOfMethylation <- function(data, breaks= c(0, 0.2, 0.8, 1), labels=c(0,0.5,1)) {
  if (!(length(breaks)-1 == length(labels)))
    stop("Labels length must be the same of length minus 1.")
  
  data[is.na(data)] <- NA

  binary <- t(apply(data, 1, function(x) {
    as.character(cut(x, breaks=breaks, labels=labels))
  }))
  
  colnames(binary) <- colnames(data)
  binary
}

# Clinical Function
resolveDuplicatedAndMakeDataFrame<- function(m) {
  if(NCOL(m)!=3)
    stop("Incorrect format. 'm' must have 3 columns: barcode, status, days")
  lu <- tapply(seq_len(NROW(m)), m[,1], function(idx) {
    r <- unique(m[idx,, drop=F])
    if (NROW(r)==1)
      return(r)
    
    if (all(r[,2]=="0")){
      j <- which.max(as.numeric(r[,3]))
      return(r[j, ,drop=F])
    }
    
    r <- r[which(r[,2]=="1"),,drop=F]
    j <- which.min(as.numeric(r[,3]))
    r[j, , drop=F]
  })
  lu <- do.call(rbind, lu)
  data.frame(status=as.numeric(lu[,2]),
             days = as.numeric(lu[,3]),
             row.names=lu[,1], stringsAsFactors=F)
}

createFollowUp <- function(followup, newtumor) {
  dropped <- setdiff(followup$bcr_patient_barcode,
                     newtumor$bcr_patient_barcode)
  recDaysNA <- rep(NA, length(dropped))
  
  names(recDaysNA) <- dropped
  
  recDays <- newtumor$days_to_new_tumor_event_after_initial_treatment
  names(recDays) <- newtumor$bcr_patient_barcode
  recDays <- c(recDays, recDaysNA)[followup$bcr_patient_barcode]
  
  barcode = followup$bcr_patient_barcode
  status  = followup$vital_status
  fup     = followup$days_to_last_followup
  death   = followup$days_to_death
  rec     = followup$new_tumor_events
  recDays = recDays
  remission = followup$followup_treatment_success
  remission1 = followup$primary_therapy_outcome_success
  cancerStatus = followup$person_neoplasm_cancer_status
  
  os <- t(sapply(seq_len(length(barcode)), function(idx) {
    if (status[idx]=="Dead") {
      os_binary = 1
      days = death[idx]
    } else if (status[idx]=="Alive"){
      os_binary = 0
      days = fup[idx]
    } else {
      stop("unrecognized stauts")
    }
    c(barcode[idx],os_binary, days)
  }))
  
  pfs <- t(sapply(seq_len(length(barcode)), function(idx) {
    if (rec[idx]!="" & rec[idx]!="NA" & rec[idx]!="NO"){
      pfs_binary = 1
      days = recDays[idx]
    } else if (status[idx]=="Dead") {
      # old condition for former data was only remission[idx] == "Progressive Disease"
      if (remission[idx] == "Progressive Disease" |
          remission1[idx] == "Progressive Disease" |
          cancerStatus[idx] == "WITH TUMOR") {
        pfs_binary=1
      } else {
        pfs_binary = 0
      }
      days = death[idx]
      
    } else if (status[idx]=="Alive"){
      pfs_binary = 0
      days = fup[idx]
    } else {
      stop("unrecognized stauts")
    }
    c(barcode[idx], pfs_binary, days)
  }))
  
  survAnnot.os  <- resolveDuplicatedAndMakeDataFrame(os)
  survAnnot.pfs <- resolveDuplicatedAndMakeDataFrame(pfs)
  return(list(os=survAnnot.os, pfs=survAnnot.pfs))
}

extractTCGAPatientsName <- function(nms){
  substr(nms, 1, 12)
}
