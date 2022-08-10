# useful blog 
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# load topGO
library("topGO")

#if(!file.exists("topgo")){
#  dir.create("topgo")
#}

tg_func <- function(universeFile, interestingGenesFile, ontology, output_file) {
  
  # test
  # universeFile = "annotations_background.txt"
  # interestingGenesFile = "genes_of_interest.txt"
  # ontology = "CC"
  # output_file = "topgo/output"
  
  # read in the 'gene universe' file
  geneID2GO <- readMappings(file = universeFile)
  geneUniverse <- names(geneID2GO)
  
  # read in the genes of interest
  genesOfInterest <- read.table(interestingGenesFile,header=FALSE)
  genesOfInterest <- as.character(genesOfInterest$V1)
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  
  # build the GOdata object in topGO
  myGOdata <- new("topGOdata", description="My project", ontology=ontology, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
  myGOdata
  
  # run the Fisher's exact tests
  resultClassic     <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  resultElim        <- runTest(myGOdata, algorithm="elim", statistic="fisher")
  resultTopgo       <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
  
  # see how many results we get where weight01 gives a P-value <= 0.001:
  mysummary <- summary(attributes(resultTopgo)$score <= 0.001)
  if(length(mysummary) > 2) {
    numsignif <- as.integer(mysummary[[3]]) # how many terms is it true that P <= 0.001
  } else {
    numsignif <- 0 
  }
  
  # print out the top 'numsignif' results:
  allRes <- GenTable(myGOdata, 
                     classicFisher = resultClassic, 
                     elimFisher = resultElim, 
                     topgoFisher = resultTopgo, 
                     parentchildFisher = resultParentchild, 
                     orderBy = "topgoFisher", 
                     ranksOf = "classicFisher", 
                     topNodes = 10) # used to be numsignif, change to 10 to avoid error
  
  allRes
  
  # write results table
  write.table(allRes, paste(output_file, ontology, "Topgo_results.txt", sep="_"), row.names = F, sep = "\t", quote = F)
  
  # print a graph (to a pdf file) with the top 'numsignif' results:
  output_file2 = paste(output_file, ontology, "Topgo", sep="_")
  printGraph(myGOdata, resultTopgo, firstSigNodes = numsignif, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)
  
  # print out the genes that are annotated with the GO terms:
  myterms <- allRes$GO.ID
  mygenes <- genesInTerm(myGOdata, myterms)
  for (i in 1:length(myterms))
  {
    myterm <- myterms[i]
    mygenesforterm <- mygenes[myterm][[1]]
    myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
    mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
    print(paste("Term",myterm,"genes:",mygenesforterm2))
  }
  
}

tg_func("top_go/annotations_background.txt", "top_go/genes_of_interest.txt", "BP", "top_go/output")
tg_func("top_go/annotations_background.txt", "top_go/genes_of_interest.txt", "MF", "top_go/output")
tg_func("top_go/annotations_background.txt", "top_go/genes_of_interest.txt", "CC", "top_go/output")

