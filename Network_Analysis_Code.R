library(bibliometrix)
library(tidyverse)

#import data
TEST <- convert2df('savedrecs(2).bib', dbsource = "wos", format = "bibtex")

#summary results
results <- biblioAnalysis(TEST, sep = ";")
S=summary(object = results, k = 24, pause = FALSE)

# Traditional network analysis of the papers:
NetMatrix <- biblioNetwork(TEST, analysis = "coupling", network = "references", sep = ". ")
net2 <- networkPlot(NetMatrix, 
                 normalize = "association", 
                 n = 24, 
                 Title = "Citation Network for Studies on Infants' Ability to Perceive Audio-Visual Congruence", 
                 size=7, 
                 labelsize = 1,
                 label.color = "black",
                 cluster = "louvain",
                 type = "fruchterman", 
                 remove.multiple=FALSE, 
                 halo = F,
                 edgesize = 2,
                 curved = F,
                 alpha = 0.9,
                 size.cex = TRUE)

#Historical citation network:
histResults <- histNetwork(TEST, min.citations = 0, sep = ";", network = TRUE, verbose = TRUE)
histPlot(histResults, n = 24, size = 12, labelsize = 4, title_as_label = FALSE, verbose = FALSE)
