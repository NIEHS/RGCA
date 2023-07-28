#extract_mixturedata.R
#Nov. 15, 2022
#Extract mixture and component data from files downloaded from
#NCATS web site.
#Uses code from read_tox21_assay_forR.R
#----------------
library(tidyverse)
#--------------------
setwd("/Users/zilberds/Desktop/tox_mixtures/code/Input/")
#Input CAS numbers of components & mixtures
file_name<-"AllCAS_ARER.txt"

caslist <- read_table(file_name,col_names=FALSE)
caslist <- unlist(caslist)
#-------------
infiles <- vector('list',5)
infiles[[1]] <- c("tox21-er-bla-agonist-p2_ratio.txt")
infiles[[2]] <- c("tox21-er-luc-bg1-4e2-agonist-p2_set1.txt",
                  "tox21-er-luc-bg1-4e2-agonist-p2_set2.txt",
                  "tox21-er-luc-bg1-4e2-agonist-p2_set3.txt")
infiles[[3]] <- c("tox21-ar-bla-agonist-p1_ratio_1.txt",
                  "tox21-ar-bla-agonist-p1_ratio_2.txt",
                  "tox21-ar-bla-agonist-p1_ratio_3.txt")
infiles[[4]] <- c("tox21-ar-mda-kb2-luc-agonist-p1_set1.txt",
                  "tox21-ar-mda-kb2-luc-agonist-p1_set2.txt",
                  "tox21-ar-mda-kb2-luc-agonist-p1_set3.txt") 
infiles[[5]] <- c("tox21-ar-mda-kb2-luc-antagonist-p1_antagonist1.txt",
                  "tox21-ar-mda-kb2-luc-antagonist-p1_antagonist2.txt",
                  "tox21-ar-mda-kb2-luc-antagonist-p1_antagonist3.txt")


assaynames <- c("tox21-er-bla-agonist-p2",
                "tox21-er-luc-bg1-4e2-agonist-p2",
                "tox21-ar-bla-agonist-p1",
                "tox21-ar-mda-kb2-luc-agonist-p1",
                "tox21-ar-mda-kb2-luc-antagonist-p1")
assayreceptors <- c('ER','ER','AR','AR', 'AR')
assaytech <- c('bla','luc','bla','luc','luc_ant')
#-----------
#Read in and process data for each assay
#Use read_delim(), from tidyverse package, not read.delim(), which
#misreads the files (not sure why).
count <- 0
librarynames <- vector('list',0)
data_allassays <- vector('list',5)

for (i in 1:5) {
  files_i <- infiles[[i]]
  nfiles_i <- length(files_i)
  all_replicate_data <- list()
  for (j in 1:nfiles_i) {
    file_ij <- files_i[j]
    
    data_ij <- read_delim(file_ij,delim='\t')
    
    #Keep only the data for the ER and AR mixtures & components
    mix_only_data <- data_ij %>%
      filter(CAS %in% caslist)
    
    #Also, the data are on the NTP-C plates; use only those 
    #plate to eliminate replicates of the components on
    #other plates.
    mix_only_data <- mix_only_data %>%
      filter(str_detect(Library,'NTP-C'))
    
    #Drop replicates of components that were not used to create the mixtures
    #Replicates of AR and ER components
    drops<-c("Tox21_400088","Tox21_400057","Tox21_400004")
    
    mix_only_data <- mix_only_data %>%
      filter(!(`Tox21 ID` %in% drops))
    
    #Add replicate (experiment) number
    if (i > 1 ) {
      #These are the assays with multiple files/assay
      mix_only_data <- 
        mix_only_data %>%
        mutate(ReplicateSet=j)
    } else {
      #Replicate number is last characer of library name
      mix_only_data <-
        mix_only_data %>%
        mutate(ReplicateSet=str_sub(Library,
                                    start=-1L,
                                    end=-1L))
    }
    
    #Keep only needed columns
    #CONC and RESP variable go from 0 to 14 for these
    CONCnames <- paste('CONC',0:14,sep='')
    RESPnames <- paste('RESP',0:14,sep='')
    
    mix_only_data <- mix_only_data %>%
      select(Assay,CAS,`Sample Name`,
             ReplicateSet,
             `Curve Class`,`Data Points`,
             `Mask No`,`Mask Flags`,
             Library,any_of(CONCnames),
             any_of(RESPnames))
    #If there are multiple input files, put everything in 
    #one big table
    all_replicate_data <- rbind(all_replicate_data,
                                mix_only_data)
    # print(nrow(all_replicate_data))
  } #end loop over j
  data_allassays[[i]] <- all_replicate_data
  
}
#Only filtering on NTP-C and the CAS list gives 294 curves.
#Some of these were from samples not included in the mixtures;
#dropping those component samples, as above, gives 256
#curves (rows in the table).
#--------------------
#Convert mask vector to separate columns of mask indicator
#variables.

#First add the variables
for (assay in 1:5) {
  for (i in 0:14) {
    newname <- str_c('mask',i)
    data_allassays[[assay]] <- 
      data_allassays[[assay]] %>%
      mutate("{newname}" := 0)
  }
}

#Then go down each table one row at a time
for (assay in 1:5) {
  for (r in 1:nrow(data_allassays[[assay]])) {
    masknumber_r <- data_allassays[[assay]]$`Mask No`[r]
    if (!is.na(masknumber_r) & masknumber_r>0) {
      #Get vector of 0s and 1s as a character string
      maskvec_r <- 
        data_allassays[[assay]]$`Mask Flags`[r]
      #Convert mask to a vector of 0s and 1s
      maskasvector <- unlist(strsplit(maskvec_r,split=" "))
      #This has 15 elements; but concentration numbers are
      #0-14 so need to offset by 1
      for (i in 0:14) {
        varname <- str_c('mask',i)
        data_allassays[[assay]][r,varname] <-
          as.numeric(maskasvector[i+1])
      } #end loop over elements of mask vector
    }
  } #end loop over rows
} #end loop over assays
#-------------------
#Eliminate columns no longer needed
for (i in 1:5) {
  data_allassays[[i]] <-
    data_allassays[[i]] %>%
    select(-`Data Points`,
           -`Mask No`,
           -`Mask Flags`,
           -Library)
  filename <- str_c(assayreceptors[i],'-',
                    assaytech[i],'.txt')
  write_tsv(data_allassays[[i]],
            file=filename)
}

