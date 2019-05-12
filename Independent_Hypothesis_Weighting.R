## Expression
library(readr)
Expression <- read_delim("Documents/BMI776/Project/051009_Expression_Match.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
# Transformation
Expression_t <- as.data.frame(t(Expression))
Expression_t <- Expression_t[-1,]
Names_exp <- as.character(unlist(Expression_t[1,]))
colnames(Expression_t) <- Names_exp
Expression_t <- Expression_t[-1,]

# Sorting
Expression_t1 <- Expression_t[order(Expression_t$BRAF_Status),]
#Expression_t2 <- sapply(Expression_t1[,4:20533], function(x)as.numeric(x))



# T-test my myself

Expression_t3 = Expression_t1[,-1:-2]
Expression_t4 = as.data.frame(t(Expression_t3))
write.csv(Expression_t4, "/Users/user/Desktop/Expression_ready.csv")

library(readr)
Expression_t5 <- read_csv("Desktop/Expression_ready.csv", col_names = FALSE)
Expression_t6 = Expression_t5[,-1]
Expression_t6
pvalue_exp = apply(Expression_t6, 1, function(x) {t.test(x[1:274], x[275:470])$p.value})

Expression_t5[,1]
# adjust
adj_pvalue_exp<-data.frame(p.adjust(pvalue_exp, method = "BH"))

pvalue_non_adj_exp <- data.frame(Gene = row.names(Expression_t4), pvalue = pvalue_exp)
pvalue_adj_exp <- data.frame(Gene = row.names(Expression_t4), pvalue = adj_pvalue_exp)
hist(pvalue_adj_exp$p.adjust.pvalue_exp..method....BH..)
hist(pvalue_exp)

## Methylation

Methylation <- read_delim("Documents/BMI776/Project/051009_Methylation_Match.txt", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
Methylation_t <- as.data.frame(t(Methylation))
Methylation_t <- Methylation_t[-1,]
Names_meth <- as.character(unlist(Methylation_t[1,]))
colnames(Methylation_t) <- Names_meth
Methylation_t <- Methylation_t[-1,]

# Sort
Methylation_t1 <- Methylation_t[order(Methylation_t$BRAF_Status),]
Methylation_t3 = Methylation_t1[,-1:-2]
Methylation_t4 = as.data.frame(t(Methylation_t3))
write.csv(Methylation_t4, "/Users/user/Desktop/Methylation_ready.csv")
Methylation_t5 <- read_csv("Desktop/Methylation_ready.csv", col_names = FALSE)
Methylation_t6 = Methylation_t5[,-1]

# Differential Expression
pvalue_meth = apply(Methylation_t6, 1, function(x) {t.test(x[1:274], x[275:470])$p.value})
adj_pvalue_meth<-data.frame(p.adjust(pvalue_meth, method = "BH"))

pvalue_non_adj_methy <- data.frame(Gene = row.names(Methylation_t4), pvalue = pvalue_meth)
pvalue_adj_methy <- data.frame(Gene = row.names(Methylation_t4), pvalue = adj_pvalue_meth)

hist(pvalue_adj_methy$p.adjust.pvalue_meth..method....BH.., xlab = "p_BH", main = "Differential Methylation")
write.csv(pvalue_adj_methy, "/Users/user/Desktop/pvalue_adj_methy.csv")
write.csv(pvalue_non_adj_exp, "/Users/user/Desktop/pvalue_non_adj_exp.csv")
write.csv(pvalue_adj_exp, "/Users/user/Desktop/pvalue_adj_exp.csv")


# Plot
hist(pvalue_non_adj_methy$pvalue, xlab = "p-value", main = "Differential Methylation")
hist(pvalue_non_adj_exp$pvalue, xlab = "p-value", main = "Differential Expression")

hist(pvalue_final$p.adjust.pvalue_exp..method....BH.., xlab = "q-value", main = "Differential Expression")
write.csv(pvalue_final, "/Users/user/Desktop/mydata_exp_qvalue.csv")

?hist

# IHW

# Check if Methylation is good variate
p_expression <- read_csv("Desktop/pvalue_non_adj_exp.csv")
p_expression_g <- split(p_expression,p_expression$DM) 
hist(p_expression_g$N$pvalue, breaks = 20)
hist(p_expression_g$Y$pvalue, breaks = 20)
hist(pvalue_non_adj_exp$pvalue, breaks = 20)
?hist

# Apply IHW
IHW = data.frame(p_exp = pvalue_exp, p_meth_BH = adj_pvalue_meth)
write.csv(IHW, "/Users/user/Desktop/IHW.csv")
IHW <- read_csv("Desktop/pvalue_non_adj_exp.csv")
library("IHW")
ihwRes <- ihw(pvalue ~ p_meth_BH, data = IHW, alpha = 0.1)
rejections(ihwRes)

sum(adj_pvalue_exp <= 0.1, na.rm = TRUE)
head(weights(ihwRes))
weights(ihwRes, levels_only = TRUE)
plot(ihwRes)
plot(ihwRes, what = "decisionboundary") 

# IHW using base mean
ihwRes_2 <- ihw(pvalue ~ baseMean, data = IHW, alpha = 0.1)
rejections(ihwRes_2)
head(weights(ihwRes_2))
weights(ihwRes, levels_only = TRUE)
plot(ihwRes_2)
plot(ihwRes_2, what = "decisionboundary") 




ihwRes_3 <- ihw(pvalue ~ baseMean + p_meth_BH, data = IHW, alpha = 0.1)
rejections(ihwRes_3)

sum(adj_pvalue_exp <= 0.1, na.rm = TRUE)
head(weights(ihwRes_3))
weights(ihwRes, levels_only = TRUE)
plot(ihwRes_3)
plot(ihwRes_3, what = "decisionboundary") 
