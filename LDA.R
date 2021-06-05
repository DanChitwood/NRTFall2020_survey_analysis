library(ggplot2) # for visualization
library(MASS) # read in for LDA

data <- read.table("./survey_results.txt", header=TRUE) # read in data

data$update <- factor(data$update, levels=c("bio1", "bio2", "bio3", "quant3","quant2", "quant1")) # reorder factor levels

r <- lda(formula =type~coding+statistics+modeling+bioinformatics+comp+mol_bio+gen_breed+plant_devo+phylogen+team+management+interdisc+sci_comm+sci_writing, data=data ) # LDA model

# LDA model results

r$prior
r$counts
r$means
r$scaling
r$svd

# check proportion of variance explained is 1

prop = r$svd^2/sum(r$svd^2)
prop

plda = predict(object=r, newdata=data) # prediction to calculate LD scores

plda$posterior #posterior probabilities

plda$x #projections/scores

lda_scores <- cbind(data, plda$x) # merge LD scores with original data

# box plot of student groups at 3 survey points

p <- ggplot(lda_scores, aes(y=LD1, fill=update))
p + geom_boxplot(color="gray30") + coord_flip() + scale_fill_manual(values=c("#4d9221","#7fbc41","#b8e186","#f1b6da","#de77ae","#c51b7d")) + theme_bw()

ggsave("boxplot.jpg")

# plot of LD1 scalings

scalings <- as.data.frame(r$scaling)

p <- ggplot(scalings, aes(x=reorder(rownames(scalings), LD1),y=LD1, fill=LD1))
p + geom_bar(stat="identity", alpha=0.8, width=0.6)+ theme_bw() + theme(axis.text.x=element_text(angle=90)) + scale_fill_gradient2(low="#4d9221",high="#c51b7d", space="Lab")

ggsave("scalings.jpg")


