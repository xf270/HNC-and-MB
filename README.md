# HNC-and-MB

#LOADING R PACKAGES#
pac <- c("survival","plyr","ggplot2","reshape2","phyloseq","elasticnet","doBy","Kendall", "vegan",
         "ade4", "lme4", "ggtree", "compositions", "scales", "DESeq2", "PMCMR", "gplots", "Hmisc")
lapply(pac, library, character.only = T)

#SET UP WORKING DIRECTORY#
setwd("H:/Dissertation/HNC/R")
load("map.Rdata")
load("map.sq.Rdata")

##TABLE 2 DESEQ2##
#ALL CaseS##
load("phy.phylum.abd.Rdata")
row.names(map) <- as.character(map$SampleID)
sample_data(phy.phylum.abd) <- map

phy.subs <- subset_taxa(phy.phylum.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 50,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:28], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:6,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)
#save(res_all, file = "res_phylum_all0804.Rdata")

#TABLE2 & TABLE 3 SCHNC#
#SQ CELL#
load("phy.phylum.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.phylum.abd) <- map.sq

phy.subs <- subset_taxa(phy.phylum.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 50,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:28], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:6,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(phylum.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(phylum.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.s))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj
#save(res_all, file = "res_phylum_sq0804.Rdata")



##CLASS LEVEL##
#SQ CELL#
load("phy.class.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.class.abd) <- map.sq

phy.subs <- subset_taxa(phy.class.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

class.deseq=phyloseq_to_deseq2(phy.subs, ~1)
class.deseq <- class.deseq[rowSums( counts(class.deseq) > 2 ) >= 50,]
class.deseq

design(class.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(class.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:34], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:12,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(class.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(class.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(class.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(class.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.s))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj
#save(res_all, file = "res_class_sq0804.Rdata")


##ORDER LEVEL##
#SQ CELL#
load("phy.order.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.order.abd) <- map.sq

phy.subs <- subset_taxa(phy.order.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                     "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

order.deseq=phyloseq_to_deseq2(phy.subs, ~1)
order.deseq <- order.deseq[rowSums( counts(order.deseq) > 2 ) >= 50,]
order.deseq

design(order.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(order.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:40], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:18,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(order.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(order.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(order.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(order.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.a))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj

#save(res_all, file = "res_order_sq0804.Rdata")


##FAMILY LEVEL##
#SQ CELL#
load("phy.family.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.family.abd) <- map.sq

phy.subs <- subset_taxa(phy.family.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                     "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

family.deseq=phyloseq_to_deseq2(phy.subs, ~1)
family.deseq <- family.deseq[rowSums( counts(family.deseq) > 2 ) >= 50,]
family.deseq

design(family.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(family.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T,independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:51], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:29,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(family.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(family.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(family.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(family.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.s))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj

save(res_all, file = "res_fam_sq0804.Rdata")


##GENUS LEVEL##
#SQ CELL#
load("phy.genus.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.genus.abd) <- map.sq

phy.subs <- subset_taxa(phy.genus.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

genus.deseq=phyloseq_to_deseq2(phy.subs, ~1)
genus.deseq <- genus.deseq[rowSums( counts(genus.deseq) > 2 ) >= 50,]
genus.deseq

design(genus.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(genus.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T,independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:72], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:50,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(genus.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(genus.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(genus.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                               ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(genus.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.s))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj
#save(res_all, file = "res_gen_sq0804.Rdata")


##SPECIES LEVEL##
#SQ CELL#
load("phy.species.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.species.abd) <- map.sq

phy.subs <- subset_taxa(phy.species.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                     "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

species.deseq=phyloseq_to_deseq2(phy.subs, ~1)
species.deseq <- species.deseq[rowSums( counts(species.deseq) > 2 ) >= 50,]
species.deseq

design(species.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case)
dds_phy2<-DESeq(species.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = T,independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res <- data.frame(res)
row.names(res) <- res$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:184], by=list(map.all[,21]), mean)), check.names = F)
phy.mean <- y.mean[2:162,]
colnames(phy.mean) <- c("Control", "Case")
res_all <- merge(res, phy.mean, by=0)


#TEST FOR INTERACTION#
#SMOKING#
design(species.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                 ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(species.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res.s <- data.frame(results(dds))
res.s <-merge(res.s, taxa, by =0)
res.s$padj <- p.adjust(res.s$pvalue, method="fdr")
res.s

#ALCOHOL#
design(species.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                 ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(species.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res.a <- data.frame(results(dds))
res.a$padj <- p.adjust(res.a$pvalue, method="fdr")
res.a<-merge(res.a, taxa, by =0)
res.a

identical(row.names(res_all), row.names(res.a))
res_all$s.q <- res.s$padj
res_all$a.q <- res.a$padj
#save(res_all, file = "res_spe_sq0804.Rdata")


##SELECTED TAXA CORRELATION MATRIX##
load("res_spe_sq0804.Rdata")
sub.spe <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
spe <- as.character(sub.spe$taxa)
spe.names <- paste(substring(as.character(sub.spe$gen), 4), substring(as.character(sub.spe$spe), 4), "(S)", sep =" ")

load("res_gen_sq0804.Rdata")
sub.gen <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
gen <- as.character(sub.gen$taxa)

load("res_fam_sq0804.Rdata")
sub.fam <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
fam <- as.character(sub.fam$taxa)

load("res_order_sq0804.Rdata")
sub.ord <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
ord <- as.character(sub.ord$taxa)
ord.names <-paste(substring(as.character(sub.ord$ord),4), "(O)", sep = " ")

load("res_class_sq0804.Rdata")
sub.cla <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
cla <- as.character(sub.cla$taxa)
cla.names <-paste(substring(as.character(sub.cla$cla),4), "(C)", sep = " ")

load("res_phylum_sq0804.Rdata")
sub.phy <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
phy <- as.character(sub.phy$taxa)
phy.names <- paste(substring(as.character(sub.phy$phy),4), "(P)", sep = " ")



load("count_spe.Rdata")
count.spe <- count_spe[spe]
colnames(count.spe) <- paste(colnames(count.spe), "_s", sep="")
load("count_ord.Rdata")
count.ord <- count_ord[ord]
colnames(count.ord) <- paste(colnames(count.ord), "_o", sep="")
load("count_cla.Rdata")
count.cla <- count_cla[cla]
colnames(count.cla) <- paste(colnames(count.cla), "_c", sep="")
load("count_phy.Rdata")
count.phy <- count_phy[phy]
colnames(count.phy) <- paste(colnames(count.phy), "_p", sep="")

identical(row.names(count.phy), row.names(count.cla))
count.all <- as.matrix(cbind(count.phy, count.cla, count.ord, count.spe))
taxa.all <- c(phy.names, cla.names, ord.names, spe.names)
colnames(count.all) <- taxa.all
colnames(count.all)
seq <- c("Proteobacteria (P)", "Betaproteobacteria (C)", "Neisseriales (O)", "Neisseria sicca (S)", 
         "Actinobacteria (P)", "Corynebacteriales (O)",  "Corynebacterium durum (S)", "Actinomyces lingnae_[NVP] (S)", "Actinomyces sp._oral_taxon_170 (S)", "Rothia dentocariosa (S)", 
         "Firmicutes (P)", "Streptococcus gordonii (S)",  "Streptococcus sanguinis (S)", "Lachnoanaerobaculum saburreum (S)", "Selenomonas sputigena (S)",
         "Bacteroidetes (P)", "Prevotella histicola (S)", "Prevotella nanceiensis (S)", "Capnocytophaga leadbetteri (S)")

count.all <- count.all[,seq]

cont.rs=rcorr(count.all, type="spearman")
cor <- cont.rs$r

corm.r <- cor[nrow(cor):1,ncol(cor):1]
melted_cormat <- melt(corm.r)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(corm.r){
  corm.r[upper.tri(corm.r)] <- NA
  return(corm.r)
}

lower_tri <- get_lower_tri(corm.r)

melted_corm.rat <- melt(lower_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_corm.rat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "darkslategray", high = "red4", mid = "white", midpoint = 0, limit = c(-1,1), 
                       space = "Lab", name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 13, hjust = 1),
        axis.text.y = element_text(size = 13, face = "italic"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) +
  coord_fixed()




##SELECTED TAXA WITH SMOKING##
load("phy.phylum.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.phylum.abd) <- map.sq

phy.subs <- subset_taxa(phy.phylum.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 50,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

#CURRENT#
res_phy2=results(dds_phy2, contrast = c("CIG_STAT", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.c = round(2^(subs$log2FoldChange), digits = 2), 
             CI.c = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.c = round(subs$pvalue, digits = 4), qadj.c  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

#FORMER#
res_phy3=results(dds_phy2, contrast = c("CIG_STAT", "2","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy3))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy3$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy3$lci <-res_phy3$log2FoldChange - 1.96*abs(res_phy3$lfcSE)
res_phy3$uci <-res_phy3$log2FoldChange + 1.96*abs(res_phy3$lfcSE)

res_phy3=data.frame(res_phy3)
res_phy3<-merge(res_phy3, taxa, by=0)
subs <- res_phy3

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.f = round(2^(subs$log2FoldChange), digits = 2), 
             CI.f = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.f = round(subs$pvalue, digits = 4), qadj.f  = round(subs$padj,digits = 4))
res2 <- data.frame(res)
row.names(res2) <- res2$taxa

identical(row.names(res1), row.names(res2))
res.all<- cbind(res1, res2[,8:11])

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:28], by=list(map.all[,6]), mean)), check.names = F)
phy.mean <- y.mean[2:6,]
colnames(phy.mean) <- c("Non", "Current", "Former")
row.names(phy.mean)
res_all <- merge(res.all, phy.mean, by=0)
save(res_all, file = "res_phylum_smk.Rdata")



#CLASS#
load("phy.class.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.class.abd) <- map.sq

phy.subs <- subset_taxa(phy.class.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

class.deseq=phyloseq_to_deseq2(phy.subs, ~1)
class.deseq <- class.deseq[rowSums( counts(class.deseq) > 2 ) >= 50,]
class.deseq

design(class.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT)
dds_phy2<-DESeq(class.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

#CURRENT#
res_phy2=results(dds_phy2, contrast = c("CIG_STAT", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.c = round(2^(subs$log2FoldChange), digits = 2), 
             CI.c = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.c = round(subs$pvalue, digits = 4), qadj.c  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

#FORMER#
res_phy3=results(dds_phy2, contrast = c("CIG_STAT", "2","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy3))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy3$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy3$lci <-res_phy3$log2FoldChange - 1.96*abs(res_phy3$lfcSE)
res_phy3$uci <-res_phy3$log2FoldChange + 1.96*abs(res_phy3$lfcSE)

res_phy3=data.frame(res_phy3)
res_phy3<-merge(res_phy3, taxa, by=0)
subs <- res_phy3

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.f = round(2^(subs$log2FoldChange), digits = 2), 
             CI.f = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.f = round(subs$pvalue, digits = 4), qadj.f  = round(subs$padj,digits = 4))
res2 <- data.frame(res)
row.names(res2) <- res2$taxa

identical(row.names(res1), row.names(res2))
res.all<- cbind(res1, res2[,8:11])

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:34], by=list(map.all[,6]), mean)), check.names = F)
phy.mean <- y.mean[2:12,]
colnames(phy.mean) <- c("Non", "Current", "Former")
row.names(phy.mean)
res_all <- merge(res.all, phy.mean, by=0)
save(res_all, file = "res_class_smk.Rdata")


#ORDER#
load("phy.order.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.order.abd) <- map.sq

phy.subs <- subset_taxa(phy.order.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

order.deseq=phyloseq_to_deseq2(phy.subs, ~1)
order.deseq <- order.deseq[rowSums( counts(order.deseq) > 2 ) >= 50,]
order.deseq

design(order.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT)
dds_phy2<-DESeq(order.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

#CURRENT#
res_phy2=results(dds_phy2, contrast = c("CIG_STAT", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.c = round(2^(subs$log2FoldChange), digits = 2), 
             CI.c = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.c = round(subs$pvalue, digits = 4), qadj.c  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

#FORMER#
res_phy3=results(dds_phy2, contrast = c("CIG_STAT", "2","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy3))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy3$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy3$lci <-res_phy3$log2FoldChange - 1.96*abs(res_phy3$lfcSE)
res_phy3$uci <-res_phy3$log2FoldChange + 1.96*abs(res_phy3$lfcSE)

res_phy3=data.frame(res_phy3)
res_phy3<-merge(res_phy3, taxa, by=0)
subs <- res_phy3

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.f = round(2^(subs$log2FoldChange), digits = 2), 
             CI.f = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.f = round(subs$pvalue, digits = 4), qadj.f  = round(subs$padj,digits = 4))
res2 <- data.frame(res)
row.names(res2) <- res2$taxa

identical(row.names(res1), row.names(res2))
res.all<- cbind(res1, res2[,8:11])

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:40], by=list(map.all[,6]), mean)), check.names = F)
phy.mean <- y.mean[2:18,]
colnames(phy.mean) <- c("Non", "Current", "Former")
row.names(phy.mean)
res_all <- merge(res.all, phy.mean, by=0)
save(res_all, file = "res_order_smk.Rdata")


#SPECIES#
load("phy.species.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.species.abd) <- map.sq

phy.subs <- subset_taxa(phy.species.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                     "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

species.deseq=phyloseq_to_deseq2(phy.subs, ~1)
species.deseq <- species.deseq[rowSums( counts(species.deseq) > 2 ) >= 50,]
species.deseq

design(species.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT)
dds_phy2<-DESeq(species.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

#CURRENT#
res_phy2=results(dds_phy2, contrast = c("CIG_STAT", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.c = round(2^(subs$log2FoldChange), digits = 2), 
             CI.c = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.c = round(subs$pvalue, digits = 4), qadj.c  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

#FORMER#
res_phy3=results(dds_phy2, contrast = c("CIG_STAT", "2","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy3))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy3$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy3$lci <-res_phy3$log2FoldChange - 1.96*abs(res_phy3$lfcSE)
res_phy3$uci <-res_phy3$log2FoldChange + 1.96*abs(res_phy3$lfcSE)

res_phy3=data.frame(res_phy3)
res_phy3<-merge(res_phy3, taxa, by=0)
subs <- res_phy3

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est.f = round(2^(subs$log2FoldChange), digits = 2), 
             CI.f = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p.f = round(subs$pvalue, digits = 4), qadj.f  = round(subs$padj,digits = 4))
res2 <- data.frame(res)
row.names(res2) <- res2$taxa

identical(row.names(res1), row.names(res2))
res.all<- cbind(res1, res2[,8:11])

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:184], by=list(map.all[,6]), mean)), check.names = F)
phy.mean <- y.mean[2:162,]
colnames(phy.mean) <- c("Non", "Current", "Former")
row.names(phy.mean)
res_all <- merge(res.all, phy.mean, by=0)
save(res_all, file = "res_species_smk.Rdata")




##SELECTED TAXA WITH ALCOHOL DRINKING##
#PHYLUM#
load("phy.phylum.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.phylum.abd) <- map.sq

phy.subs <- subset_taxa(phy.phylum.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 50,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day +  ETHA_GRAMS_PER_DAY + Case + CIG_STAT + alc_yn)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("alc_yn", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:28], by=list(map.all[,22]), mean)), check.names = F)
phy.mean <- y.mean[2:6,]
colnames(phy.mean) <- c("Non", "Ever","Missing")
row.names(phy.mean)
res_all <- merge(res1, phy.mean, by=0)
save(res_all, file = "res_phylum_alc.Rdata")



#class#
load("phy.class.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.class.abd) <- map.sq

phy.subs <- subset_taxa(phy.class.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

class.deseq=phyloseq_to_deseq2(phy.subs, ~1)
class.deseq <- class.deseq[rowSums( counts(class.deseq) > 2 ) >= 50,]
class.deseq

design(class.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day +  ETHA_GRAMS_PER_DAY + Case + CIG_STAT + alc_yn)
dds_phy2<-DESeq(class.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

#CURRENT#
res_phy2=results(dds_phy2, contrast = c("alc_yn", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:34], by=list(map.all[,22]), mean)), check.names = F)
phy.mean <- y.mean[2:12,]
colnames(phy.mean) <- c("Non", "Ever","Missing")
row.names(phy.mean)
res_all <- merge(res1, phy.mean, by=0)
save(res_all, file = "res_class_alc.Rdata")


#ORDER#
load("phy.order.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.order.abd) <- map.sq

phy.subs <- subset_taxa(phy.order.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

order.deseq=phyloseq_to_deseq2(phy.subs, ~1)
order.deseq <- order.deseq[rowSums( counts(order.deseq) > 2 ) >= 50,]
order.deseq

design(order.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day +  ETHA_GRAMS_PER_DAY + Case + CIG_STAT + alc_yn)
dds_phy2<-DESeq(order.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("alc_yn", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                          UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:40], by=list(map.all[,22]), mean)), check.names = F)
phy.mean <- y.mean[2:18,]
colnames(phy.mean) <- c("Non", "Ever","Missing")
row.names(phy.mean)
res_all <- merge(res1, phy.mean, by=0)
save(res_all, file = "res_order_alc.Rdata")




#SPECIES#
load("phy.species.abd.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.species.abd) <- map.sq

phy.subs <- subset_taxa(phy.species.abd, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                     "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.subs))

species.deseq=phyloseq_to_deseq2(phy.subs, ~1)
species.deseq <- species.deseq[rowSums( counts(species.deseq) > 2 ) >= 50,]
species.deseq

design(species.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day +  ETHA_GRAMS_PER_DAY + Case + CIG_STAT + alc_yn)
dds_phy2<-DESeq(species.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("alc_yn", "1","0"), cooksCutoff = T, independentFiltering = T)

identical(rownames(dds_phy2),rownames(res_phy2))
mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by=0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, ord = subs$Order, fam = subs$Family,
             gen = subs$Genus, spe = subs$Species,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             CI = paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                        UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p = round(subs$pvalue, digits = 4), qadj  = round(subs$padj,digits = 4))
res1 <- data.frame(res)
row.names(res1) <- res1$taxa

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.sq, count_phy, by.x = "SampleID", by.y = "row.names")
y.mean <- data.frame(t(aggregate(map.all[,24:184], by=list(map.all[,22]), mean)), check.names = F)
phy.mean <- y.mean[2:162,]
colnames(phy.mean) <- c("Non", "Ever","Missing")
row.names(phy.mean)
res_all <- merge(res1, phy.mean, by=0)
save(res_all, file = "res_species_alc.Rdata")



##EXTRACT OR##
load("res_spe_sq0804.Rdata")
sub.spe <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
spe <- as.character(sub.spe$taxa)
spe.names <- paste(substring(as.character(sub.spe$gen), 4), substring(as.character(sub.spe$spe), 4), "(S)", sep =" ")

load("res_order_sq0804.Rdata")
sub.ord <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
ord <- as.character(sub.ord$taxa)
ord.names <-paste(substring(as.character(sub.ord$ord),4), "(O)", sep = " ")

load("res_class_sq0804.Rdata")
sub.cla <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
cla <- as.character(sub.cla$taxa)
cla.names <-paste(substring(as.character(sub.cla$cla),4), "(C)", sep = " ")

load("res_phylum_sq0804.Rdata")
sub.phy <- subset(res_all, as.numeric(as.character(res_all$qadj)) <=0.1)
phy <- as.character(sub.phy$taxa)
phy.names <- paste(substring(as.character(sub.phy$phy),4), "(P)", sep = " ")


load("res_phylum_smk.Rdata")
OR.smk1<-res_all[,c("taxa", "phy", "Est.c", "qadj.c", "Est.f", "qadj.f")]
OR.smk1<-OR.smk1[sub.phy$taxa,]

load("res_class_smk.Rdata")
OR.smk2<-res_all[,c("taxa", "cla", "Est.c", "qadj.c", "Est.f", "qadj.f")]
OR.smk2<-OR.smk2[sub.cla$taxa,]

load("res_order_smk.Rdata")
OR.smk3<-res_all[,c("taxa", "ord", "Est.c", "qadj.c", "Est.f", "qadj.f")]
OR.smk3<-OR.smk3[sub.ord$taxa,]

load("res_species_smk.Rdata")
OR.smk4<-res_all[,c("taxa", "spe", "Est.c", "qadj.c", "Est.f", "qadj.f")]
OR.smk4<-OR.smk4[sub.spe$taxa,]

colnames(OR.smk1) <- c("taxa", "name",  "Est.c", "qadj.c", "Est.f", "qadj.f")
colnames(OR.smk2) <- c("taxa", "name",  "Est.c", "qadj.c", "Est.f", "qadj.f")
colnames(OR.smk3) <- c("taxa", "name",  "Est.c", "qadj.c", "Est.f", "qadj.f")
colnames(OR.smk4) <- c("taxa", "name",  "Est.c", "qadj.c", "Est.f", "qadj.f")
OR.smk <- rbind(OR.smk1, OR.smk2, OR.smk3, OR.smk4)
#save(OR.smk, file = "OR.smk.Rdata")


load("res_phylum_alc.Rdata")
OR.alc1<-res_all[,c("taxa", "phy", "Est", "qadj")]
OR.alc1<-OR.alc1[sub.phy$taxa,]

load("res_class_alc.Rdata")
OR.alc2<-res_all[,c("taxa", "cla", "Est", "qadj")]
OR.alc2<-OR.alc2[sub.cla$taxa,]

load("res_order_alc.Rdata")
OR.alc3<-res_all[,c("taxa", "ord", "Est", "qadj")]
OR.alc3<-OR.alc3[sub.ord$taxa,]

load("res_species_alc.Rdata")
OR.alc4<-res_all[,c("taxa", "spe", "Est", "qadj")]
OR.alc4<-OR.alc4[sub.spe$taxa,]

colnames(OR.alc1) <- c("taxa", "name",  "Est", "qadj")
colnames(OR.alc2) <- c("taxa", "name",  "Est", "qadj")
colnames(OR.alc3) <- c("taxa", "name",  "Est", "qadj")
colnames(OR.alc4) <- c("taxa", "name",  "Est", "qadj")
OR.alc <- rbind(OR.alc1, OR.alc2, OR.alc3, OR.alc4)
#save(OR.alc, file = "OR.alc.Rdata")

seq <- c("p__Proteobacteria", "c__Betaproteobacteria", "o__Neisseriales", "s__sicca", 
         "p__Actinobacteria", "o__Corynebacteriales",  "s__durum", "s__lingnae_[NVP]", "s__sp._oral_taxon_170", "s__dentocariosa", 
         "p__Firmicutes", "s__gordonii",  "s__sanguinis", "s__saburreum", "s__sputigena",
         "p__Bacteroidetes", "s__histicola", "s__nanceiensis", "s__leadbetteri")

#OR HEATMAP#
load("OR.smk.Rdata")
OR.s <- OR.smk[, c("name", "Est.c", "Est.f")]
load("OR.alc.Rdata")
OR.a <- OR.alc[, c("name", "Est")]
identical(OR.s$name, OR.a$name)
OR <- cbind(OR.s, OR.a)
OR <- OR[order(seq),]
row.names(OR) <- OR$name
OR.h<-OR[,c("Est.c", "Est.f", "Est")]
#ADD SCHNC ORs#
OR.h$hnc <- c(0.90, 0.76, 0.70, 0.47, 1.37, 0.61, 0.61, 1.53, 2.09, 1.75,
              1.14, 1.48, 1.57, 0.64, 0.53, 1.21, 1.87, 0.49, 0.60)
OR.h<-apply(OR.h, 2, as.numeric)
row.names(OR.h) <- seq

lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.3,4)
lhei = c(0.3,5,1)
breaks = c(seq(0.3, 0.999, length=100),seq(1, 1.01, length=2), seq(1.02, 3, length = 100))
my_palette <- colorRampPalette(c("navyblue", "white", "orange4"))

#PLOT#
group <- c("Current Smokers", "Former Smokers", "Alcohol Drinkers", "SCHNC case vs. control")
or <-     heatmap.2(OR.h, col=my_palette, Colv=F, Rowv=F, breaks=breaks,
                    dendrogram="none", symkey=F, labRow=F, srtCol=90,
                    symm=F, symbreaks=T, scale="none", trace="none", labCol=group, key = F,
                    margins=c(8,10), cexCol = 1.2, cexRow = 2.0,
                    colsep=1:ncol(as.matrix(OR.h)), rowsep=1:nrow(as.matrix(OR.h)),
                    lmat = lmat, lwid = lwid, lhei = lhei)

cor <-    heatmap.2(OR.h, col= my_palette, breaks = breaks, trace="none",  density.info="none", dendrogram=c("row"), 
                    key.title = "Fold Changes", keysize = 1.2, key.ylab = NULL,
                    symm=F, symkey=F,symbreaks=T, scale="none",  densadj = 0.25)

