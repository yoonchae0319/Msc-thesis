## -----------------------------------------------------------------------------
library(eSMC2)
library(e1071) ##hamming distance calculation
#?hamming.distance() ##calculate the number of different elements between the chosen two vectors
#hamming.distance(os[1:(n - 2), 1],os[1:(n - 2), 2]  )

## -----------------------------------------------------------------------------
#' Reads a vcf file
#'
#' @param path : location of the vcf file
#' @return A list containing all lines of the vcf file
Get_vcf_data<-function(path){
  DNAseqfile=list()
  count_DNA=0
  con = file(path, "r")
  while ( TRUE ) {
    count_DNA= count_DNA+1
    line = readLines(con, n = 1)
    if ( length(line) == 0 ){
      break
    }
    DNAseqfile[[count_DNA]]=as.vector(unlist(strsplit(as.vector(line),"\t")))
  }
  close(con)
  return(DNAseqfile)
}


## -----------------------------------------------------------------------------
#' Builds a segregatic matrix from a vcf file
#'
#' @param path : location of the vcf file
#' @export
#' @return The segregating matrix of the vcf file
Process_multiplemarkers_vcf<-function(path){
  data=Get_vcf_data(path)
  
  i=0
  
  for(search_id in 1:length(data)){
    if(data[[search_id]][1]=="#CHROM"){
      start_position=search_id+1
      break()
    }
  }
  for(search_col in 1:length(data[[start_position]])){
    if(data[[(start_position-1)]][search_col]=="FORMAT"){
      start_col=search_col+1
      break()
    }
  }
  
  M=length(data[[start_position]][c(start_col:length(data[[start_position]]))])*2
  
  if(length(which(as.vector(strsplit(data[[start_position]][(start_col-1)],"")[[1]])==":"))>1){
    not_good=T
  }else{
    not_good=F
  }
  
  #checking if all rows belong to the same chromosome
  if(length(unique(matrix(as.vector(unlist(data[start_position:length(data)])),nrow=(length(data)-(start_position-1)),ncol=length(data[[start_position]]),byrow = T)[,1]))==1){
    
    multiplemarkers_mat=matrix(0,ncol = length(data),nrow = M+2)
    #moving on to the position
    for(ii in start_position:length(data)){
      if(data[[ii]][4] != "."){ #where there is no reference allele
        i=i+1
        multiplemarkers_mat[(M+2),i]=as.numeric(data[[ii]][2])
        if(i==1){
          multiplemarkers_mat[(M+1),i]=as.numeric(data[[ii]][2])
        }else{
          multiplemarkers_mat[(M+1),i]=as.numeric(data[[ii]][2])- as.numeric(multiplemarkers_mat[(M+2),(i-1)])
        }
        
        for(ccc in start_col:length(data[[ii]])){
          pos_temp=which(strsplit(data[[ii]][ccc],"")[[1]]=="|")
          if(length(pos_temp)>1){
            data[[ii]][ccc]=substr(data[[ii]][ccc],1,min(pos_temp))
          }
        }
        
        temp_seq=strsplit(paste(data[[ii]][c(start_col:length(data[[ii]]))],collapse = ""),"")[[1]]
        #print(temp_seq)
        if(not_good){
          pos=which(temp_seq=="|")
          pos=sort(c((pos+1),(pos-1)))
          multiplemarkers_mat[1:M,i]=temp_seq[pos]
        }else{
          multiplemarkers_mat[1:M,i]=temp_seq[which(temp_seq!="|")]
        }
        for (allele in 0:7) {
          pos = which(as.numeric(multiplemarkers_mat[1:M, i]) == allele)
          if (length(pos) > 0) {
            if (allele == 0) {
              multiplemarkers_mat[pos, i] = data[[ii]][4]  # Reference allele
            } else {
              alt_alleles = strsplit(data[[ii]][5], ",")[[1]]
              if (allele <= length(alt_alleles)) {
                multiplemarkers_mat[pos, i] = alt_alleles[allele]
              }
            }
          }
        }
      }
    }
    
    pos_0=which(as.numeric(multiplemarkers_mat[dim(multiplemarkers_mat)[1],])==0)
    if(length(pos_0)>1){
      multiplemarkers_mat=multiplemarkers_mat[,-pos_0] #here they are restracting everything so need to fix this
    }
  }else{
    multiplemarkers_mat=list()
    #print(dim(multiplemarkers_mat))
    count_chr=0
    for(chr in unique(matrix(as.vector(unlist(data[start_position:length(data)])),nrow=(length(data)-(start_position-1)),ncol=length(data[[start_position]]),byrow = T)[,1])){
      count_chr=count_chr+1
      i=0
      pos_chr=which(matrix(as.vector(unlist(data[start_position:length(data)])),nrow=(length(data)-(start_position-1)),ncol=length(data[[start_position]]),byrow = T)[,1]==chr)
      multiplemarkers_mat[[count_chr]]=matrix(0,ncol = length(pos_chr),nrow = (M+2))
      for(ii in sort(pos_chr)){
        ii=start_position+ii-1
        if(nchar(data[[ii]][4])==1&nchar(data[[ii]][5])==1){
          
          i=i+1
          
          multiplemarkers_mat[[count_chr]][(M+2),i]=as.numeric(data[[ii]][2])
          if(i==1){
            multiplemarkers_mat[[count_chr]][(M+1),i]=as.numeric(data[[ii]][2])
          }else{
            multiplemarkers_mat[[count_chr]][(M+1),i]=as.numeric(data[[ii]][2])- as.numeric(multiplemarkers_mat[[count_chr]][(M+2),(i-1)])
          }
          for(ccc in start_col:length(data[[ii]])){
            pos_temp=which(strsplit(data[[ii]][ccc],"")[[1]]=="|")
            if(length(pos_temp)>1){
              data[[ii]][ccc]=substr(data[[ii]][ccc],1,min(pos_temp))
            }
          }
          temp_seq=strsplit(paste(data[[ii]][c(start_col:length(data[[ii]]))],collapse = ""),"")[[1]]
          if(not_good){
            pos=which(temp_seq=="|")
            pos=sort(c((pos+1),(pos-1)))
            multiplemarkers_mat[[count_chr]][1:M,i]=temp_seq[pos]
          }else{
            multiplemarkers_mat[[count_chr]][1:M,i]=temp_seq[which(temp_seq!="|")]
          }
          for (allele in 0:7) {
            pos = which(as.numeric(multiplemarkers_mat[1:M, i]) == allele)
            if (length(pos) > 0) {
              if (allele == 0) {
                multiplemarkers_mat[pos, i] = data[[ii]][4]  # Reference allele
              } else {
                alt_alleles = strsplit(data[[ii]][5], ",")[[1]]
                if (allele <= length(alt_alleles)) {
                  multiplemarkers_mat[pos, i] = alt_alleles[allele]
                }
              }
            }
          }
        }
      }
      pos_0=which(as.numeric(multiplemarkers_mat[[count_chr]][dim(multiplemarkers_mat[[count_chr]])[1],])==0)
      if(length(pos_0)>1){
        multiplemarkers_mat[[count_chr]]=multiplemarkers_mat[[count_chr]][,-pos_0]
      }
    }
  }
}


## -----------------------------------------------------------------------------
##obtaining the distance matrix along the sequence (window by window)
##we got vcf file and segregating could be considered if that's easier
##this is to calculate the distribution of the first betty number and the mean length of the zero homology group per window along the sequence

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
sim_number <- args[2]

os <- Process_multiplemarkers_vcf(vcf_file)
print(dim(os))
n <- nrow(os)
position1 <- 0
window <- 5e4
positions <- as.numeric(os[n, ])  # Last row contains positions
a <- 0
allele_map <- c(
  "A" = "1", "G" = "2", "C" = "3", "T" = "4",
  "100" = "5", "101" = "6", "102" = "7", "103" = "8", "104" = "9",
  "10" = "10", "11" = "11", "12" = "12", "13" = "13", "14" = "14",
  "15" = "15", "16" = "16", "17" = "17"
)
segSites <- c()

while (position1 < max(positions)) {
  subset <- c()
  position2 <- position1 + window
  a <- a + 1
  
  for (i in 1:length(positions)) {
    if (positions[i] >= position1 & positions[i] < position2) {
      subset <- c(subset, i)
    }
  }

  sites <- 0

  if (length(subset) > 0) {
    segregating_matrix <- os[1:(n - 2), subset, drop = FALSE]

    # Hamming distance calculation
    iii <- c()
    for (i in 1:(n - 2)) {
      for (j in 1:(n - 2)) {
        iii <- c(iii, hamming.distance(segregating_matrix[i, ], segregating_matrix[j, ]))
      }
    }
    distance_matrix <- matrix(unlist(iii), nrow = n - 2, ncol = n - 2, byrow = TRUE)
    write.table(distance_matrix, file = paste0("distanceMatrix_", sim_number, "_", a, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

    # Allele mapping
    segregating_matrix_mapped <- apply(segregating_matrix, c(1, 2), function(x) {
      if (x %in% names(allele_map)) allele_map[[x]] else x
    })

    # Check row pairing
    if (nrow(segregating_matrix_mapped) %% 2 != 0) {
      warning("Uneven number of rows in window ", a, " – skipping this window.")
      position1 <- position1 + window
      next
    }

    # EMIBD input matrix
    final <- t(sapply(seq(1, nrow(segregating_matrix_mapped), by = 2), function(i) {
    pair <- segregating_matrix_mapped[c(i, i + 1), , drop = FALSE]
    if (ncol(pair) == 1) {
    # When only one column, return a named vector
      setNames(paste(pair[1, 1], pair[2, 1]), colnames(pair))
    } else {
      apply(pair, 2, function(x) paste(x, collapse = " "))
    }
  }))

  if (is.matrix(final) && nrow(final) == 1) {
    final <- t(final)
  }
    input_file <- paste0("EMIBD_input_", sim_number, "_", a, ".dat")
    write.table(final, file = input_file, row.names = TRUE, col.names = FALSE, quote = FALSE)

    # Write parameter file
    n_indiv <- (nrow(os) - 2) / 2
    n_loci <- ncol(segregating_matrix_mapped)
    output_file <- paste0("EMIBD_output_", sim_number, "_", a, ".txt")
    par_lines <- c(
      sprintf("%d                    !Integer, NumIndiv, #Individuals", n_indiv),
      sprintf("%d                !Integer, NumLoci, #Loci", n_loci),
      "0                    !Integer, DataForm, 0/1/2",
      "1                    !Boolean, Inbreed, Inbreeding: 1/0=Y/N",
      paste0("\"", input_file, "\"      !String, GtypeFile, Gtype file path & name"),
      paste0("\"", output_file, "\"            !String, OutFileName, output file path & name"),
      "111                  !Integer, ISeed, Random number seed",
      "1                    !Boolean, RndDelta0, 1/0=Y/N",
      "1                    !EM_Method, 0/1=Update Delta only/Update delta & allele freq jointly",
      "0                    !Boolean, OutAlleleFre, 1/0=Y/N"
    )
    writeLines(par_lines, con = paste0("EMIBD_", sim_number, "_", a, ".par"))

    for (j in 1:length(subset)) {
      if (length(unique(segregating_matrix[, j])) != 1) {
        sites <- sites + 1
      }
    }

  } else {
    # Handle windows with no data
    segregating_matrix <- matrix(0, nrow = n - 2, ncol = 1)
    iii <- c()
    for (i in 1:(n - 2)) {
      for (j in 1:(n - 2)) {
        iii <- c(iii, hamming.distance(segregating_matrix[i, ], segregating_matrix[j, ]))
      }
    }
    distance_matrix <- matrix(unlist(iii), nrow = n - 2, ncol = n - 2, byrow = TRUE)
    write.table(distance_matrix, file = paste0("distanceMatrix_", sim_number, "_", a, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  

    final <- t(sapply(seq(1, nrow(segregating_matrix), by = 2), function(i) {
    pair <- segregating_matrix[c(i, i + 1), , drop = FALSE]
    if (ncol(pair) == 1) {
    # When only one column, return a named vector
      setNames(paste(pair[1, 1], pair[2, 1]), colnames(pair))
    } else {
      apply(pair, 2, function(x) paste(x, collapse = " "))
    }
  }))

if (is.matrix(final) && nrow(final) == 1) {
  final <- t(final)
}

    input_file <- paste0("EMIBD_input_", sim_number, "_", a, ".dat")
    write.table(final, file = input_file, row.names = TRUE, col.names = FALSE, quote = FALSE)

  # Parameter file (same structure)
    n_indiv <- (nrow(os) - 2) / 2
    n_loci <- ncol(segregating_matrix)
    output_file <- paste0("EMIBD_output_", sim_number, "_", a, ".txt")
    par_lines <- c(
      sprintf("%d                    !Integer, NumIndiv, #Individuals", n_indiv),
      sprintf("%d                !Integer, NumLoci, #Loci", n_loci),
      "0                    !Integer, DataForm, 0/1/2",
      "1                    !Boolean, Inbreed, Inbreeding: 1/0=Y/N",
      paste0("\"", input_file, "\"      !String, GtypeFile, Gtype file path & name"),
      paste0("\"", output_file, "\"            !String, OutFileName, output file path & name"),
      "111                  !Integer, ISeed, Random number seed",
      "1                    !Boolean, RndDelta0, 1/0=Y/N",
      "1                    !EM_Method, 0/1=Update Delta only/Update delta & allele freq jointly",
      "0                    !Boolean, OutAlleleFre, 1/0=Y/N"
    )
    writeLines(par_lines, con = paste0("EMIBD_", sim_number, "_", a, ".par"))
  }

  segSites <- c(segSites, sites)
  segSites_matrix <- matrix(segSites, ncol = 1)
  write.table(segSites_matrix, file = paste0("segSites_", sim_number, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

  position1 <- position1 + window
}

##nohup Rscript distancematrixgen.R > distancematrixgen_output.log 2>&1 &