library(bibliometrix)
library(tidyverse)
library(Matrix)

#import data
citations <- convert2df('/Users/au620441/Desktop/Meta_Analysis_Documents/citations_24.bib', dbsource = "wos", format = "bibtex")

#summary results
results <- biblioAnalysis(citations, sep = ";")
S=summary(object = results, k = 24, pause = FALSE)

# Traditional network analysis of the papers:
set.seed(23)
NetMatrix <- biblioNetwork(citations, analysis = "coupling", network = "references", sep = ". ")
a <- b <- c("PEJOVIC, 2020", 
            "DORN, 2018", 
            "RICHOZ, 2017", 
            "GUELLAI, 2016", 
            "STRERI, 2016", 
            "ALTVATER-MACKENSEN, 2016", 
            "DE BOISFERON, 2015", 
            "LEWKOWICZ, 2015",
            "KUBICEK, 2014-1",
            "KUBICEK, 2014-2" ,
            "LEWKOWICZ, 2013",
            "KUBICEK, 2013", 
            "PONS, 2009",
            "TREHUB, 2009",
            "BAHRICK, 2005",
            "PATTERSON, 2002",
            "ALDRIDGE, 1999",
            "PATTERSON, 1999",
            "BAHRICK, 1998",
            "PICKENS, 1994",
            "POULINDUBOIS, 1994",
            "WALKER-ANDREWS, 1991",
            "KUHL, 1984",
            "MACKAIN, 1983")
dimnames(NetMatrix) <- list(a, b)
net2 <- networkPlot(NetMatrix, 
                    normalize = "association", 
                    n = 24, 
                    Title = "Citation Network Analysis of Studies on Infants' Ability to Perceive Audio-Visual Congruence", 
                    size=5, 
                    labelsize = 0.9,
                    label.color = "black",
                    cluster = "louvain",
                    type = "fruchterman", 
                    label.cex = F,
                    remove.multiple=FALSE, 
                    halo = F,
                    edgesize = 3,
                    curved = F,
                    alpha = 0.9,
                    size.cex = TRUE)

#Historical citation network:
x <- localCitations(TEST, fast.search = FALSE, sep = ";")
H=histNetwork(TEST,min.citations = 1, sep=";", network=FALSE)
histResults$NetMatrix
histResults <- histNetwork(TEST, min.citations = 1, sep = ";", network = TRUE, verbose = TRUE)
histPlot(histResults, n = 24, size = 12, labelsize = 4, title_as_label = FALSE, verbose = FALSE)
