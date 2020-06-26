
library(tidyverse)
library(ggplot2)
data <- read.delim("test_Ecl_chr2/gene_distribution_tip_v_anc.txt", header = TRUE, stringsAsFactors = TRUE)

head(data)

anc1 <- data[data$AncestralChr==1, ]
anc1 <- anc1[anc1$TipChr!="all",]
anc1$TipChr <- as.integer(anc1$TipChr)
anc1$TipChrGenes <- as.integer(anc1$TipChrGenes)
anc1 <- as_tibble(anc1)
ggplot(data = anc1, aes(x=TipChr, y=Overlap)) +
    facet_wrap(.~TipGenome) +
    geom_line()+
    ggtitle("MGLO Ancestral genome 6 Chromosome 1")
ggsave("MGLO_ancestral_genome6_chromosome_1.pdf")

anc2 <- data[data$AncestralChr==2, ]
anc2 <- anc2[anc2$TipChr!="all",]
anc2$TipChr <- as.integer(anc2$TipChr)
anc2$TipChrGenes <- as.integer(anc2$TipChrGenes)
anc2 <- as_tibble(anc2)
ggplot(data = anc2, aes(x=TipChr, y=Overlap)) +    
  facet_wrap(.~TipGenome) +
  geom_line()+
  ggtitle("MGLO Ancestral genome 6 Chromosome 2")
ggsave("MGLO_ancestral_genome6_chromosome_2.pdf")

anc3 <- data[data$AncestralChr==3, ]
anc3 <- anc3[anc3$TipChr!="all",]
anc3$TipChr <- as.integer(anc3$TipChr)
anc3$TipChrGenes <- as.integer(anc3$TipChrGenes)
anc3 <- as_tibble(anc3)
ggplot(data = anc3, aes(x=TipChr, y=Overlap)) +
  facet_wrap(TipGenome~.) +
  geom_line()+
  ggtitle("MGLO Ancestral genome 6 Chromosome 3")
ggsave("MGLO_ancestral_genome6_chromosome_3.pdf")


anc4 <- data[data$AncestralChr==4, ]
anc4 <- anc4[anc4$TipChr!="all",]
anc4$TipChr <- as.integer(anc4$TipChr)
anc4$TipChrGenes <- as.integer(anc4$TipChrGenes)
anc4 <- as_tibble(anc4)
ggplot(data = anc4, aes(x=TipChr, y=Overlap)) +
  facet_wrap(TipGenome~.) +
  geom_line()+
  ggtitle("MGLO Ancestral genome 6 Chromosome 4")
ggsave("MGLO_ancestral_genome6_chromosome_4.pdf")
  
anc5 <- data[data$AncestralChr==5, ]
anc5 <- anc5[anc5$TipChr!="all",]
anc5$TipChr <- as.integer(anc5$TipChr)
anc5$TipChrGenes <- as.integer(anc5$TipChrGenes)
anc5 <- as_tibble(anc5)
ggplot(data = anc5, aes(x=TipChr, y=Overlap)) +
  facet_wrap(TipGenome~.) +
  geom_line()+
  ggtitle("MGLO Ancestral genome 6 Chromosome 5")
ggsave("MGLO_ancestral_genome6_chromosome_5.pdf")

anc6 <- data[data$AncestralChr==6, ]
anc6 <- anc6[anc6$TipChr!="all",]
anc6$TipChr <- as.integer(anc6$TipChr)
anc6$TipChrGenes <- as.integer(anc6$TipChrGenes)
anc6 <- as_tibble(anc6)

  ggtitle("MGLO Ancestral genome 6 Chromosome 6")
ggsave("MGLO_ancestral_genome6_chromosome_6.pdf")

anc <- data[data$AncestralChr != 2, ]
anc <- anc[anc$AncestralChr != 3,]
anc <- anc[anc$TipChr != "all",]
anc$AncestralChr <- as.factor(anc$AncestralChr)
ggplot(data = anc, aes(x=TipChr, y=Overlap)) +
  facet_wrap(TipGenome~.) +
  geom_point(aes(color = AncestralChr))+
  ggtitle("MGLO Ancestral genome 6")
ggsave("MGLO_ancestral_genome6_distribution_on_all_tips.pdf")
