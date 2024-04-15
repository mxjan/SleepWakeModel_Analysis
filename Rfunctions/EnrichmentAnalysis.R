ScoreEnrichment<-function(geneList,Score="decreasing",data="Mouse"){
  require(topGO)
  require(ggrepel)
  
  if (data=="Mouse"){
    require(org.Mm.eg.db)
    
    # Generate topGO data
    GOdata <- new("topGOdata",
                  description = "GO analysis",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSel = function(x)x,
                  annot = annFUN.org,
                  mapping = "org.Mm.eg.db",
                  ID = "symbol",nodeSize = 10)#
  }else if (data == "Human"){
    
    require(org.Hs.eg.db)
    
    # Generate topGO data
    GOdata <- new("topGOdata",
                  description = "GO analysis",
                  ontology = "BP",
                  allGenes = geneList,
                  geneSel = function(x)x,
                  annot = annFUN.org,
                  mapping = "org.Hs.eg.db",
                  ID = "symbol",nodeSize = 10)#
  }

  
  
  # Run KS test weight01 algorithm
  resultks <- runTest(GOdata, algorithm = "weight01", statistic = "ks", scoreOrder = Score)
  
  # Keep significant only
  allRes <- GenTable(GOdata, classicFisher = resultks,
                     orderBy = "resultks", ranksOf = "resultks", topNodes = length(resultks@score))
  allResf<-allRes[allRes$classicFisher<0.05,]
  
  # Get all gene in terms
  GeneInTermL<-genesInTerm(GOdata)
  
  # # Rank GeneList
  # if (Score=="decreasing"){
  #   RankGeneList<-rank(geneList)
  # }else{
  #   RankGeneList<-rank(geneList*-1)
  # }

  
  # Get mean score per term
  allResf$MeanRank<-sapply(allResf$GO.ID,function(x){
    genes <- GeneInTermL[[x]] # get genes associated with term
    scores.in <- geneList[names(geneList) %in% genes]
    return(sum(scores.in)/length(scores.in))
  })
  
  allResf$classicFisher<-as.numeric(allResf$classicFisher)
  allResf$MeanRank<-as.numeric(allResf$MeanRank)
  allResf$Significant<-as.numeric(allResf$Significant)
  
  if (Score == "decreasing"){
    lowcol<-"red"
    highcol<-"purple"
    FiltScore<-allResf$MeanRank>=sort(allResf$MeanRank,decreasing=T)[5]#quantile(allResf$MeanRank,prob=.95)
  }else{
    lowcol<-"orange"
    highcol<-"red"
    FiltScore<-allResf$MeanRank<=sort(allResf$MeanRank,decreasing=F)[5]#quantile(allResf$MeanRank,prob=.05)
  }
  gg<-ggplot(aes(x=MeanRank,y=-log10(classicFisher),label=Term),data=allResf)+ #,color=MeanRank
    geom_point() + #aes(size=Significant)
    geom_text_repel(aes(label=ifelse(FiltScore | classicFisher <= sort(classicFisher)[5],strtrim(as.character(Term),25),'')),size=2,min.segment.length=.1,max.time = 5)+
    scale_color_gradient(low=lowcol, high=highcol) #classicFisher < quantile(classicFisher,prob=0.05)
  gg<-gg+theme_classic()+ylab("-log10 p-value")+xlab("Mean Contribution [%]")
  return(list(plot=gg,df=allResf[order(allResf$classicFisher,decreasing = F),],dfAll=allResf,GeneInTermL=GeneInTermL))
  
}



enrichment<-function(Genes,GL.universe){
  require(topGO)
  require(org.Mm.eg.db)
  
  geneList <- factor(as.integer(GL.universe %in% Genes))
  names(geneList) <- as.factor(GL.universe)
  
  # Generate topGO data
  GOdata <- new("topGOdata",
                description = "GO analysis",
                ontology = "BP",
                allGenes = geneList,
                geneSel = function(x)x,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "symbol",nodeSize = 10)#
  
  
  resultclassic <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultclassic,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(resultclassic@score))
  allRes$FC<-allRes$Significant/allRes$Expected
  allRes$Term <-as.character(topGO:::.getTermsDefinition(allRes$GO.ID,ontology(GOdata),numChar = 10000))
  allRes$classicFisheradj<-p.adjust(allRes$classicFisher,method = "fdr")
  return(allRes)
}

enrichmentHS<-function(Genes,GL.universe,onto="BP"){
  require(topGO)
  require(org.Hs.eg.db)
  
  geneList <- factor(as.integer(GL.universe %in% Genes))
  names(geneList) <- as.factor(GL.universe)
  
  # Generate topGO data
  GOdata <- new("topGOdata",
                description = "GO analysis",
                ontology = onto,
                allGenes = geneList,
                geneSel = function(x)x,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "symbol",nodeSize = 10)#
  
  
  resultclassic <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultclassic,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(resultclassic@score))
  allRes$FC<-allRes$Significant/allRes$Expected
  allRes$Term <-as.character(topGO:::.getTermsDefinition(allRes$GO.ID,ontology(GOdata),numChar = 10000))
  allRes$classicFisheradj<-p.adjust(allRes$classicFisher,method = "fdr")
  return(allRes)
}


enrichmentMF<-function(Genes,GL.universe){
  geneList <- factor(as.integer(names(GL.universe) %in% Genes))
  names(geneList) <- as.factor(names(GL.universe))
  # Generate topGO data
  GOdata <- new("topGOdata",
                description = "GO analysis",
                ontology = "MF",
                allGenes = geneList,
                geneSel = function(x)x,
                annot = annFUN.org,
                mapping = "org.Mm.eg.db",
                ID = "symbol",nodeSize = 10)#
  resultclassic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = resultclassic,
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(resultclassic@score))
  allRes$FC<-allRes$Significant/allRes$Expected
  allRes$Term <-as.character(topGO:::.getTermsDefinition(allRes$GO.ID,ontology(GOdata),numChar = 10000))
  allRes$classicFisheradj<-p.adjust(allRes$classicFisher,method = "fdr")
  return(allRes)
}