
library(yaml)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
masterDir <- "/Volumes/SNPFISH_O/DNMT1 Concentration Ratio Experiments/"
outdir <- "/Volumes/SNPFISH_O/TestOutput/"


AnalyzeSNPDirectory(masterDir, outdir, SNPTableFlag = TRUE)


test <- GetAllYamlGroupsinDir("/Volumes/SNPFISH_O/DNMT1 Concentration Ratio Experiments/")
groups <- ExtractFilesFromYamlGroup(test[[1]], masterDir)
outdir <- "/Volumes/SNPFISH_O/TestOutput/"



#Load data files
exp.data <- groups[[1]]
guide.table <- groups[[2]]
cy.table <- groups[[3]]
tmr.table <- groups[[4]]


#Make sure 3-color and undetec are at the end
guide.table <- ReorderSNPLevels(guide.table)
cy.table <- ReorderSNPLevels(cy.table)
tmr.table <- ReorderSNPLevels(tmr.table)

MakeGuideTablePlots(exp.data, guide.table, outdir)

MakeStackedCellCountPlot(exp.data, guide.table, outdir)
MakeTotalCellCountPie(exp.data, guide.table, outdir)
MakeOverallImbalanceBar(exp.data, guide.table, outdir)
MakeAlleleScatterPlot(exp.data, guide.table, outdir)



MakeGuideTablePlots(exp.data, guide.table, outdir)

##
by.cell.counts.guide.filt <- guide.table %>% group_by(cellID, labels) %>%
  summarize(rnaCount = n()) %>% spread(labels, rnaCount) %>%
  gather(labels, rnaCount, -cellID) %>%
  mutate(rnaCount = replace(rnaCount, is.na(rnaCount), 0)) %>%
  filter((!labels %in% c("undetec", "3-color"))) %>% spread(labels, rnaCount)

by.cell.counts.cy.tot <- cy.table %>% group_by(cellID) %>%
  summarize(cySpots = n())

by.cell.counts.tmr.tot <- tmr.table %>% group_by(cellID) %>%
  summarize(tmrSpots = n())

by.cell.counts.SNPSpots <- left_join(by.cell.counts.tmr.tot, by.cell.counts.cy.tot)
by.cell.counts.forPlot <- left_join(by.cell.counts.SNPSpots, by.cell.counts.guide.filt)


cy.probe <- exp.data$channels$cy$probe


detected.vs.total.bar.plot <- ggplot() +
  geom_bar(data = by.cell.counts.forPlot, aes(x = cellID, y= ggname(exp.data$channels$cy$probe))
                                              , stat = "identity", position = "identity") +
  xlab("Individual Cells") +
  ylab("RNA Count") +
  theme_bw(base_size = 20, base_family = "Helvetica") +
  theme(axis.ticks = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("B6" ="deepskyblue1", "C7" = "orange")) +
  ggtitle(exp.data$name)
plot(detected.vs.total.bar.plot)


## Makes stacked cell count plot
frac.detec <- by.cell.counts.guide %>% group_by(cellID) %>%
  summarize(fracDetec = 1 - rnaCount[labels == 'undetec']/sum(rnaCount),
                          totalRNACount = sum(rnaCount))

stacked.cell.count.plot <- ggplot() +
  geom_bar(data = by.cell.counts.cy, aes(x = cellID, y = rnaCount, fill = labels), stat = "identity") +
  geom_text(data = frac.detec, aes(label = format(100*fracDetec, digits=1, drop0trailing=TRUE),
                 x = cellID, y= totalRNACount ), stat= "identity", vjust = -1) +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle(exp.data$name)
plot(stacked.cell.count.plot)

currdir <- getwd()
setwd(outdir)
ggsave(plot =  stacked.cell.count.plot,  width = 9.89, height = 9,
       filename =  paste(exp.data$name, "_Stacked_Cell_Count_Plot.pdf", sep = ""))
setwd(currdir)


##Total Cell Count Ratio
totals <- guide.table %>% group_by(labels) %>% summarize(rnaCount = n()) %>% mutate(totalFrac = rnaCount/sum(rnaCount))



## Plot Pie Graph
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

overall.pie.chart <- ggplot(data = totals) +
  geom_bar(aes(x="", y = totalFrac, fill = labels), stat = "identity", width = 1) +
  coord_polar(theta="y") + blank_theme +
  theme(axis.text.x=element_blank()) +
  geom_text(aes(x = 1.8,  y = totalFrac/3 + c(0,cumsum(totalFrac[-length(totalFrac)])),
                label = percent(totalFrac)), size=5) +
  scale_fill_brewer(palette="Set1") +
  ggtitle(exp.data$name)
plot(overall.pie.chart)

currdir <- getwd()
setwd(outdir)
ggsave(plot = overall.pie.chart,  width = 9.89, height = 9,
       filename =  paste(exp.data$name, "_Overall_Pie_Chart.pdf", sep = ""))
setwd(currdir)


## SNPA vs SNPB
total.SNP <- guide.table %>% group_by(labels) %>% summarize(rnaCount = n()) %>%
  filter((!labels %in% c("undetec", "3-color"))) %>% mutate(rnaFrac = rnaCount/sum(rnaCount))

stacked.allelic.ratio <- ggplot() +
  geom_bar(data = total.SNP, aes(x = labels, y = rnaFrac), stat = "identity") +
  ggtitle(exp.data$name) +
  ylim(0,1) +
  geom_text(data = total.SNP, aes(x = labels, y = rnaFrac + 0.1, label = percent(rnaFrac)), size = 5) +
  theme_minimal() +
  scale_fill_brewer(palette="Set1")
plot(stacked.allelic.ratio)

setwd(outdir)
ggsave(plot = StackedAllelicRatio,  width = 9.89, height = 9,
       filename =  paste(exp.data$name, "_Stacked_Alleilic_Ratio.pdf", sep = ""))
setwd(currdir)

## Scatter plot of SNPA vs SNPB

by.cell.counts.allele <- guide.table %>% group_by(cellID, labels) %>%
  summarize(rnaCount = n()) %>% filter(!labels %in% c('undetec', '3-color') ) %>%
  spread(labels, rnaCount)

by.cell.counts.allele[is.na(by.cell.counts.allele)] <- 0

labels <- levels(guide.table$labels)

correlation <- cor(by.cell.counts.allele[labels[1]], by.cell.counts.allele[labels[2]])

allele.scatterplot <- ggplot() +
  geom_point(data = by.cell.counts.allele, aes_string(x = ggname(labels[1]), y = ggname(labels[2]))) +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle(paste(exp.data$name, " Correlation = ", correlation, sep = ""))
plot(allele.scatterplot)

setwd(outdir)
ggsave(plot = StackedAllelicRatio,  width = 9.89, height = 9,
       filename =  paste(exp.data$name, "_Alleilic_Counts_Scatter.pdf", sep = ""))
setwd(currdir)


## Scatterplot in D

#Make sure 3-color and undetec are at the end
initial.levels <- levels(cy.table$labels)

undetec.index <-  which(initial.levels == "undetec")
three.color.index <-  which(initial.levels == "3-color")
initial.levels <- initial.levels[-undetec.index]
initial.levels <- initial.levels[-three.color.index]
final.levels <- append(initial.levels, c("undetec", "3-color"),
                       after = length(initial.levels))

guide.table$labels <- factor(cy.table$labels, levels = final.levels)


by.cell.counts <- guide.table %>% group_by(cellID, labels) %>% summarize(rnaCount = n())



## Makes stacked cell count plot
stacked.cell.count.plot <- ggplot() +
  geom_bar(data = by.cell.counts, aes(x = cellID, y = rnaCount, fill = labels), stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  ggtitle(exp.data$name)
plot(stacked.cell.count.plot)



