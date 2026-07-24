UVMAttributes <-
  c("SYMBOL","CHROM","POS","STRAND","Protein_position","Amino_acids","Codons","Consequence","VARIANT_CLASS","REF","ALT","BIOTYPE","SYMBOL","UNIPROT_ISOFORM","Feature","ENSP","SWISSPROT", "TREMBL", "UNIPARC","Genotype","NORMAL_DEPTH","TUMOR_DEPTH","NORMAL_AD_REF","NORMAL_AD_ALT","TUMOR_AD_REF","TUMOR_AD_ALT","AF","gnomADg_AF","gnomADe_AF","SOMATIC","Existing_variation")
UVMNewAttributes <- c("Hugo_Symbol","Chromosome","Start_Position","STRAND","Protein_position","Amino_acids","Codons","Variant_Classification","Variant_Type","REF","ALT","BIOTYPE","SYMBOL","UNIPROT_ISOFORM","Transcript_ID","ENSP","SWISSPROT", "TREMBL", "UNIPARC","Genotype","NORMAL_DEPTH","TUMOR_DEPTH","NORMAL_AD_REF","NORMAL_AD_ALT","TUMOR_AD_REF","TUMOR_AD_ALT","AF","gnomADg_AF","gnomADe_AF", "SOMATIC","Existing_variation","CALLER")
inocrasAttributes <- c("SYMBOL","CHROM","POS","Codons","Consequence","VARIANT_CLASS","REF","ALT","BIOTYPE","Feature","MANE_SELECT","SWISSPROT")
inocrasNewAttributes <- c("Hugo_Symbol","Chromosome","Start_Position","Codons","Variant_Classification","Variant_Type","REF","ALT","BIOTYPE","Transcript_ID","Transcript_ID_ENST","SWISSPROT","CALLER","Amino_acids","Protein_position")


getUniprotFromDB <- function(ids=NULL,dbFile) {
  db <- read.csv(dbFile,header=TRUE)
  df <- db[db$From %in% ids,c("From","Entry")]
  rownames(df) <- 1:nrow(df)
  df
}

getUniprotFromDB_old <- function(ids=NULL,dbFile) {
  db <- read.csv(dbFile)
  if (!is.null(ids)) {
    db <- db[db$From %in% ids,]
  }
  db
}

purify <- function(file) {
  tempDir <- tempdir(check=TRUE)
  tempFile <- tempfile(tmpdir=tempDir)
  cmd <- paste0("grep -E \"PASS\" ",file," > ",tempFile)
  system(cmd)
  if (file.info(tempFile)$size==0) {
    return(NULL)
  }
  return(tempFile)
}


processDataset <- function(file,attributes,createTabixFile,variantTypes,type,caller,normalName=NULL,
                           tumorName=NULL) {
  #get VCF, combine if necessary
  print("FILTERING ONLY VARIANTS THAT PASS")
  isPASSPrefilter <- function(x) grepl("PASS", x, fixed=TRUE) & grepl("missense_variant", x, fixed=TRUE)
  vcf <- readVEPFromVCF(file=file,
                        onlyPassing = TRUE,
                        createTabixFile = createTabixFile,
                        variantTypes = variantTypes,
                        preFilterFunction = isPASSPrefilter,
                        genome="grc38",
                        filterFunction=NA,
                        normalName=normalName,
                        tumorName=tumorName)
  if (is.null(vcf) || nrow(vcf)==0) return(NULL)
  if (nrow(vcf)>0) {
    vcfTmp <- vcf[,attributes]
    vcfTmp$Caller <- caller
    if (type=="inocras") {
      protein <- vcf$HGVSp
      proteinAA <- gsub("\\d+","/",gsub(".*p\\.","",protein))
      coordinate <- gsub(".*?(\\d+).*","\\1",gsub(".*p\\.","",protein))
      vcfTmp$Amino_Acids <- sapply(strsplit(proteinAA,"/"),function(x) {paste0(a(x[1]),"/",a(x[2]))})
      vcfTmp$Protein_position <- paste0(coordinate,"/",coordinate)
    }

  }
  return(vcfTmp)
}

transformToDataFrame <- function(results) {
  lens <- sapply(results,length)
  idx <- which(lens != length(results[[1]]))
  r <- results[[1]]
  if (length(idx)>0) {
    for (i in idx) {
      row <- results[[i]]
      extraColumns <- setdiff(colnames(r),colnames(row))
      for (col in extraColumns) {
        row[col] <- NA
      }
      row <- row[colnames(r)]
      results[[i]] <- row
    }
  }
  r <- do.call(rbind,results)

  resultsWithMutationInformation <- addTopologicalInformationAndMutation(r)
  rowOrder <- c( "subcellularLocationComplete", "topologicalDomain", c("intramembraneComplete",  "transmembraneComplete", "intramembrane", "transmembrane","cytoplasm", "extracellular", "extracellularMutationQ", "secretedProteinQ", "membraneKeywords","secretedButNotMutationInECDQ","membraneButNotMutationInECDorSecretedKeywords") )
  x<-c(setdiff(colnames(resultsWithMutationInformation),rowOrder),
    rowOrder)
  resultsWithMutationInformation <- resultsWithMutationInformation[,x]
  resultsWithMutationInformation
}


processSarekResults <- function(datasetLocation,resultsFile,redoTabix,uniprotDBFile=NULL,geneSequenceTableFile=NULL) {
  #gather all the variant information
  callers <- c("strelka","mutect2")
  print("Collecting variants")
  completeVCF <- do.call(rbind,lapply(1:2,function(i) {
    print(callers[i])
    caller=callers[i]
    location <- paste0(datasetLocation,"/",caller)
    if (caller %in% c("mutect2")) {
      file <- list.files(path=location,recursive=TRUE,pattern=".*vs.*.vcf.gz$",full.names = TRUE)
    } else if (caller %in% c("strelka")) {
      file <- list.files(path=location,recursive=TRUE,pattern=".*vs.*snvs.*.vcf.gz$",full.names = TRUE)
      fileIndels <- list.files(path=location,recursive=TRUE,pattern=".*indels.*vcf.gz$",full.names = TRUE)
    }
    VEPAttributes <- attributes
    m1<-processDataset(file=file,
                    attributes=VEPAttributes,
                    format=VEPFormat,
                    createTabixFile=redoTabix,
                    variantTypes=c("missense_variant"),
                    caller=caller)
    m1
  }))
  colnames(completeVCF) <- newAttributes
  resultsWithCellularLocation <- addSubcellularLocation(completeVCF,uniprotDBFile,geneSequenceTableFile)
  resultsWithCellularLocation
}



readBroadMafFile <- function(file,attributes) {
      firstLine <- grep('Hugo_Symbol',readLines(file))
      df = read.csv(file, skip = firstLine-1, header = T, sep="\t")
      df <- df[df$Variant_Classification=="Missense_Mutation",]
      tmp<-df[,c("Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Hugo_Symbol","Protein_Change","Codon_Change","SwissProt_acc_Id","Annotation_Transcript","t_ref_count","t_alt_count","n_alt_count", "n_ref_count","HGNC_UniProt_ID.supplied_by_UniProt.", "HGNC_Ensembl_ID.supplied_by_Ensembl.")]
      tmp$SYMBOL <- tmp$Hugo_Symbol
      tmp$Amino_acids <- gsub("p.","",tmp$Protein_Change)
      aaPosition <- gsub(".*?([0-9]+).+","\\1",tmp$Amino_acids)
      tmp$Protein_position <- paste0(aaPosition,"/",aaPosition)
      tmp$Amino_acids <- gsub("\\d+","/",tmp$Amino_acids)
      tmp$Codons <- gsub(".*)","",tmp$Codon_Change)
      tmp$Codons <- gsub(">","/",tmp$Codon_Change)
      tmp$BIOTYPE <- "protein_coding"
      tmp$ENSP <- NA
      tmp$SWISSPROT <- NA
      tmp$TREMBL <- NA
      tmp$UNIPARC <- NA
      tmp$CALLER <- "mutect"
      tmp$REF <- tmp$Reference_Allele
      tmp$ALT <- tmp$Tumor_Seq_Allele2
      tmp$Variant_Type <- "SNV"
      tmp$UNIPROT_ISOFORM <- tmp$HGNC_Ensembl_ID.supplied_by_Ensembl.
      tmp$Transcript_ID <- tmp$Annotation_Transcript
      tmp$REF_DP <- tmp$t_ref_count
      tmp$ALT_DP <- tmp$t_alt_count
      tmp$NORMAL_DEPTH <- NA
      tmp$TUMOR_DEPTH  <- NA
      tmp$NORMAL_AD_REF <- tmp$n_ref_count
      tmp$NORMAL_AD_ALT <- tmp$n_alt_count
      tmp$TUMOR_AD_REF <- tmp$t_ref_count
      tmp$TUMOR_AD_ALT <- tmp$t_alt_count
      tmp$AF <- NA
      tmp$Genotype <- NA
      tmp <- tmp[,attributes]
      tmp
}

processVCFResults <- function(vcfGZFile,caller,redoTabix,uniprotDBFile=NULL,geneSequenceTableFile=NULL,attributes,newAttributes,type, IDFilter, normalName=NULL,tumorName=NULL) {
                                        #gather all the variant information
  if (type=="Broad") {
    completeVCF<-readBroadMafFile(vcfGZFile,newAttributes)
  } else {
    VEPAttributes <- attributes
    completeVCF<-processDataset(file=vcfGZFile,
                     attributes=VEPAttributes,
                     createTabixFile=redoTabix,
                     variantTypes=c("missense_variant"),
                     type=type,
                                caller=caller,
                                normalName=normalName,
                                tumorName=tumorName)
    if (is.null(completeVCF)) return(NULL)
    colnames(completeVCF) <- newAttributes
  }
  resultsWithCellularLocation <- addSubcellularLocation(completeVCF,
                                                        uniprotDBFile,
                                                        geneSequenceTableFile,
                                                        IDFilter=IDFilter)
  resultsWithCellularLocation
}


addSubcellularLocation <- function(completeVCF,uniprotDBFile,geneSequenceTableFile=NULL,IDFilter="ensembl_transcript_id") {

  print("Fetching Uniprot ID")

  uniprotDB <- readRDS(uniprotDBFile)
  #subcellularInformation <- unique(uniprotDB[,c("Entry","Subcellular.location..CC.","Intramembrane", "Topological.domain","Transmembrane")])
  #rownames(subcellularInformation) <- subcellularInformation$Entry
  subcellularInformation <- uniprotDB

  if (!is.null(geneSequenceTableFile)) {
    geneSequenceTable <- read.csv(geneSequenceTableFile)
  } else {
    geneSequenceTable <- NULL
  }

  vcf <- completeVCF
  dataWithUniprotID <- assignUniprotIDToGeneSNPCoordinate(
      data=vcf,
      uniprot=uniprotDB,
      genomeVersion="grch38",
      IDAttribute="Transcript_ID",
      IDFilter=IDFilter,
      type="TCGA", #This is not TCGA but it has the same format
      skipErrors=FALSE,
                         geneSequenceTable=geneSequenceTable,
                         skipGenome=TRUE)

  print(paste0("IDs: ",length(unique(vcf$Transcript_ID)),"; assigned: ",length(unique(dataWithUniprotID$Transcript_ID))))

#  subcellularInformation <- GetSubcellular_location(unique(uniprotTable$Entry))
  resultsWithCellularLocationList <- addUniprotCellularLocationAndTCGAPatientInformationFromManifest(
    allData=dataWithUniprotID,
    uniprotTable=uniprotDB,
    patientInformationTable=NA,
    attributesToAdd=NA,
    allSubcellularInformation=uniprotDB,
    uniprotAttribute="Entry")
  tmp<-resultsWithCellularLocationList[which(sapply(resultsWithCellularLocationList,length)>0)]

  resultsWithCellularLocation <- addTopologicalInformation(tmp)
  return(resultsWithCellularLocation)
}

processRNASeqResults <- function(resultsFile) {
  countdata <- read.table(resultsFile,header = TRUE) #Load with read.csv
  #Creates a gene_name df with names and ids
  cols <- setdiff(colnames(countdata), c("gene_id", "gene_name"))
  geneNames <- countdata$gene_id
  countsMatrix <- countdata[,cols]
  rownames(countsMatrix) <- geneNames

  #filter out low count genes
  keep <- rowSums(countsMatrix) > 10
  countsMatrix <- countsMatrix[keep,]
  design <- data.frame(sample=colnames(countsMatrix), condition=c("Tumor","Normal"))
  rownames(design) <- design$sample
  design$condition <- factor(design$condition, levels=c("Normal","Tumor"))
  dds <- DESeqDataSetFromMatrix(countData = round(countsMatrix), colData = design, design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds <- estimateSizeFactors(dds)
  normalizedCounts <- data.frame(counts(dds,normalized=TRUE))
  normalizedCounts[normalizedCounts$KF53==0,"KF53"] <- 0.00001
  normalizedCounts[normalizedCounts$HF51==0,"HF51"] <- 0.00001
  normalizedCounts$ratio <- normalizedCounts$HF51/normalizedCounts$KF53
  normalizedCounts
}

createDifferenceMatrix <- function(variantFiles,variantListNames,variantListCallers,normalNames,tumorNames) {

  variantList <- lapply(1:length(variantFiles), function(i) {
    print(i)
    file <- variantFiles[i]
    caller <- variantListCallers[i]
    if (caller=="broad") {
      r<-readBroadMafFileForDifferenceMatrix(file)
    } else {
      tmpDir <- tempdir(check = TRUE)
      tmpFile <- file.path(tmpDir, basename(file))
      file.copy(file, tmpFile, overwrite = TRUE)
      is_bgzipped <- function(f) identical(readBin(f, raw(), 4L), as.raw(c(0x1f, 0x8b, 0x08, 0x04)))
      if (is_bgzipped(tmpFile)) {
        file <- tmpFile
      } else {
        system(paste0("gunzip --verbose -q -f ", tmpFile))
        filename <- tools::file_path_sans_ext(tmpFile)
        if (!file.exists(filename)) filename <- tmpFile
        gzipFile <- paste0(tools::file_path_sans_ext(filename), ".gz")
        Rsamtools::bgzip(filename, gzipFile, overwrite = TRUE)
        file <- gzipFile
      }
      if (caller %in% c("mutect","strelka","inocras","both","azenta","tnscope")) {
        isPASSPrefilter <- function(x) grepl("PASS", x, fixed=TRUE) & grepl("missense_variant", x, fixed=TRUE)
        postParseFilter <- NA
      } else if (caller=="muse2") {
        isPASSPrefilter <- function(x) grepl("PASS", x, fixed=TRUE)
        postParseFilter <- NA
      } else if (caller=="freebayes") {
        isPASSPrefilter <- function(x) grepl("missense_variant", x, fixed=TRUE)
        postParseFilter <- function(x) { fixed(x)$QUAL >= 1000 }
      }
      r<-readVEPFromVCF(file=file,
                        createTabixFile=TRUE,
                        onlyPassing=TRUE,
                        variantTypes="missense_variant",
                        genes="ALL",
                        genome="hg38",
                        preFilterFunction=isPASSPrefilter,
                        filterFunction=postParseFilter,
                        normalName=normalNames[i],
                        tumorName=tumorNames[i])
      if (!is.null(r)) {
        if (caller=="muse2") {
          r$CHROM <- r$seqnames
          r$POS <- r$start
        }
        r$change <- paste0(r$CHROM,"_",r$POS,"_",r$REF,"_",sapply(1:nrow(r),function(l) { toString(r$ALT[[l]])}))
      }
    }
    r
  })

  variantListNamesTmp <- sapply(1:length(variantListNames) ,function(i) { paste0(variantListNames[i],"(",length(unique(variantList[[i]]$change)), ")")})
  n <- length(variantListNames)
  differenceMatrix <- matrix(data=rep(0,n*n), nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in 1:n) {
      r1 <- variantList[[i]]
      r2 <- variantList[[j]]
      inCommon <- intersect(r1$change,r2$change)
      differenceMatrix[i,j] <- paste0(length(inCommon)," (",round(100*length(inCommon)/length(unique(r1$change)),1),"%)")
    }
  }
  colnames(differenceMatrix) <- variantListNamesTmp
  rownames(differenceMatrix) <- variantListNamesTmp
  differenceMatrix
}

readBroadMafFileForDifferenceMatrix <- function(file) {
      firstLine <- grep('Hugo_Symbol',readLines(file))
      df = read.csv(file, skip = firstLine-1, header = T, sep="\t")
      df <- df[df$Variant_Classification=="Missense_Mutation",]
      tmp<-df[,c("Chromosome","Start_Position","Strand","Reference_Allele","Tumor_Seq_Allele2","Variant_Classification","Hugo_Symbol")]
      colnames(tmp) <- c("CHROM","POS", "strand","REF", "ALT", "Consequence", "SYMBOL")
      tmp$width<-1
      tmp$VCFRowID <- 1:nrow(tmp)
      tmp$IMPACT <- NA
      tmp$change <- paste0(tmp$CHROM,"_",tmp$POS,"_",tmp$REF,"_",tmp$ALT)
      tmp
}


processVariantFiles <- function(variantFiles,variantListCallers,sequenceDB,uniprotDBFile,normalNames,tumorNames) {
  allVariants <- c()
  for (i in 1:length(variantFiles)) {
    file <- variantFiles[i]
    paste0("Working on ",file)
    caller <- variantListCallers[i]
    variantListName <- variantListNames[i]
    if (grepl("UVM",variantListNames[i])) {
      attributes <- UVMAttributes
      newAttributes <- UVMNewAttributes
      type <- "UVM"
      IDFilter <- "ensembl_transcript_id"
    } else if (grepl("inocras",variantListNames[i])) {
      attributes <- inocrasAttributes
      newAttributes <- inocrasNewAttributes
      type <- "inocras"
      IDFilter <- "refseq_mrna"
    } else if (grepl("Azenta",variantListNames[i])) {
      attributes <- UVMAttributes
      newAttributes <- UVMNewAttributes
      type <- "Azenta"
      IDFilter <- "ensembl_transcript_id"
    } else if (grepl("Broad",variantListNames[i])) {
      attributes <- UVMAttributes
      newAttributes <- UVMNewAttributes
      type <- "Broad"
      IDFilter <- "Transcript_ID"
    }

    r<-processVCFResults(vcfGZFile=file,
	caller=caller,
	redoTabix=TRUE,
	uniprotDBFile=uniprotDBFile,
  	geneSequenceTableFile=sequenceDB,
        attributes=attributes,
        newAttributes=newAttributes,
        type=type,
        IDFilter=IDFilter,
        normalName=normalNames[i],
        tumorName=tumorNames[i]
 	)
    if (is.null(r)) next
    r$ALT <- sapply(r$ALT,function(x) {as.character(x)})
    r$change <- paste0(r$Chromosome,"_",r$Start_Position,"_",r$Amino_acids)
    r$Caller <- caller
    r$nCaller <- 1

    if (is.null(allVariants) || length(allVariants)==0) {
      allVariants <- r
    } else {
      for (j in 3:nrow(r)) {
        idx <- which(allVariants$change %in% r[j,"change"])
        if (length(idx)==0) {
          allVariants <- rbind(allVariants,r[j,colnames(allVariants)])
        } else {
          allVariants[idx,"Caller"] <- paste0(allVariants[idx,"Caller"],";",caller)
          allVariants[idx,"nCaller"] <- allVariants[idx,"nCaller"] + 1
        }
      }
    }
  }
  allVariants
}

calculateAndAddDifferentialExpression <- function(allVariants, sampleSheetFile, countsMatrixFile ) {

#,quantFiles,tx2geneFile,samplesheetFile) {
    samples <- read.csv(sampleSheetFile,header=TRUE)
    samples$condition <- ifelse(grepl("Blood",samples$sample),"Normal","Tumor")
    samples$condition <- factor(samples$condition)

    counts <- read.csv(countsMatrixFile,sep="\t")
    allVariants<-merge(allVariants,counts,by.x="Hugo_Symbol",by.y="gene_name")
    return(allVariants)
}

createSubcellularLocationTable <- function(allVariants) {
  callers <- unique(allVariants$Caller)
  subcellularLocationTable <- data.frame()
  for (i in 1:length(callers)) {
    caller <- callers[i]
    r <- allVariants[allVariants$Caller==caller,]
    r$change <- paste0(r$Chromosome,"_",r$Start_Position,"_",r$Amino_acids)
    df <- data.frame(Name=caller,Total=length(unique(r$change)),ECD=sum(r$extracellularMutationQ==TRUE),Secreted=sum(r$secretedButNotMutationInECDQ),NonDefined=0.5*sum(grepl("MEMBRANE",r$membraneButNotMutationInECDorSecretedKeywords)))
    df$TotalMSPs=df$ECD+df$Secreted+df$NonDefined
    df$TotalMSPs=paste0(df$TotalMSPs," (",round(100*df$TotalMSPs/df$Total),"%)")
    subcellularLocationTable <- rbind(subcellularLocationTable,df)
  }
  subcellularLocationTable
}


addVariantCountInformation <- function(variants,countsFile, columnName) {
  counts <- NA
  for (i in 1:length(countsFile)) {
    c<- read.csv(countsFile[i],sep="\t")
    if (i==1) {
      counts <- c
    } else {
      counts <- unique(rbind(counts,c))
    }
  }

  columnName <- paste0(columnName,"_Ref_Alt")

  for (i in 1:nrow(variants)) {
    v <- variants[i,]
    if ( (!columnName %in% colnames(variants)) ||
        (is.na(variants[i,columnName])) )  {
      idx <- counts$contig==v$Chromosome & counts$position==v$Start_Position & counts$refAllele==v$REF & counts$altAllele==v$ALT
      if (sum(idx)==0) {
        variants[i,columnName]=NA
        } else {
          count <- counts[idx,c("refCount","altCount")]
          variants[i,columnName]=paste(count,collapse="_")
        }
    }
  }
  variants
}


addVariantCountInformation_old <- function(variants,tumorCountsFile, normalCountsFile) {
    tumorCounts <- read.csv(tumorCountsFile,sep="\t")
    normalCounts <- read.csv(normalCountsFile,sep="\t")
    for (i in 1:nrow(variants)) {
      v <- variants[i,]
      if ( (!"NormalAlleleCount" %in% colnames(variants)) ||
           (is.na(variants[i,"NormalAlleleCount"])) )  {
	normalIdx <- normalCounts$contig==v$Chromosome & normalCounts$position==v$Start_Position & normalCounts$refAllele==v$REF & normalCounts$altAllele==v$ALT
        if (sum(normalIdx)==0) {
           variants[i,"NormalAlleleCount"]=NA
        } else {
          normalCount <- normalCounts[normalIdx,c("refCount","altCount")]
          variants[i,"NormalAlleleCount"]=paste(normalCount,collapse="_")
        }
      }
      if ( (!"TumorAlleleCount" %in% colnames(variants)) ||
           (is.na(variants[i,"TumorAlleleCount"]) ))  {
        tumorIdx <- tumorCounts$contig==v$Chromosome & tumorCounts$position==v$Start_Position & tumorCounts$refAllele==v$REF & tumorCounts$altAllele==v$ALT
        if (sum(tumorIdx)==0) {
          variants[i,"TumorAlleleCount"]=NA
        } else {
          tumorCount <- tumorCounts[tumorIdx,c("refCount","altCount")]
          variants[i,"TumorAlleleCount"]=paste(tumorCount,collapse="_")
        }
      }
    }
    variants
}

parseTumorPurity <- function(file_path) {
  lines <- readLines(file_path)

  extractValue <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)
    if (length(line) == 0) return(NA_character_)
    val <- trimws(sub(pattern, "", line[1]))
    if (val == "NA") return(NA_character_)
    val
  }

  header_idx <- grep("^Sample\t", lines)
  if (length(header_idx) > 0) {
    vals  <- strsplit(lines[header_idx + 1], "\t")[[1]]
    sample_id      <- vals[1]
    tumor_fraction <- as.numeric(vals[2])
    ploidy         <- as.numeric(vals[3])
  } else {
    sample_id      <- extractValue("^Sample:\t*")
    tumor_fraction <- as.numeric(extractValue("Tumor Fraction:\t*"))
    ploidy         <- as.numeric(extractValue("Ploidy:\t*"))
  }

  data.frame(Sample = sample_id, Tumor_Fraction = tumor_fraction, Ploidy = ploidy)
}

compareAndSave <- function(var, file, replaceMaster = FALSE, verbose = TRUE) {
  extension <- tools::file_ext(file)
  masterFile <- paste0(tools::file_path_sans_ext(file), "_Master.", extension)

  read_file <- function(f) {
    if (toupper(extension) == "RDS") readRDS(f)
    else if (toupper(extension) == "CSV") read.csv(f)
    else stop(paste0("Unsupported extension: ", extension))
  }

  save_file <- function(x, f) {
    if (toupper(extension) == "RDS") saveRDS(x, file = f)
    else if (toupper(extension) == "CSV") write.csv(x, file = f, row.names = FALSE)
  }

  if (file.exists(masterFile)) {
    masterData <- read_file(masterFile)
    if (!identical(var, masterData)) {
      stop(paste0("Not the same: ", file, " and ", masterFile))
    }
    #cols <- intersect(colnames(var),colnames(masterData))
    #result <- all.equal(var[, sort(cols)], masterData[, sort(cols)])
    #if (!isTRUE(result)) {
    #  stop(paste0("Differences found:\n", paste(result, collapse = "\n")))
    #

    if (verbose) message("Files are identical: ", masterFile)
  } else {
    if (verbose) message("Master file not found. Saving current data as master: ", masterFile)
  }

  if (!file.exists(masterFile) || replaceMaster) {
    save_file(var, masterFile)
  }
  save_file(var, file)
}

