library(dplyr)
library(ggplot2)
library(reshape2)
filenames = c("AR-bla.txt", "AR-luc.txt", "ER-bla.txt", "ER-luc.txt")

df = read.delim("/Users/zilberds/Desktop/tox_mixtures/ER-AR-data/AR-luc.txt")
df_alt =read.delim("/Users/zilberds/Desktop/tox_mixtures/ER-AR-data/AR-bla.txt")
df = data.frame(df)
head(df)

mk_vec = function(arow) unlist(array(arow))
mix_df <- df %>% filter(substr(CAS, 1,5)=="NOCAS")
pure_df <- df %>% filter(substr(CAS, 1,5)!="NOCAS")
conc_idx = 6:20
resp_idx = 21:35
mask = 36:50
plot(mk_vec(df[4, conc_idx]), mk_vec(df[4, resp_idx]), log = "x")


sample_names = unique(pure_df$Sample.Name)

flu_idx = which(df$Sample.Name == "Progesterone")
flu_idx_alt = which(df$Sample.Name == "Fluoxymestrone")
x_vals = cbind(t(df[flu_idx, conc_idx]), t(df[flu_idx_alt, conc_idx]))
y_vals = cbind(t(df[flu_idx, resp_idx]), t(df[flu_idx_alt, resp_idx]))

cbind(t(df[flu_idx, conc_idx]), 1)

#*(1-df[flu_idx, mask])
matplot(x=x_vals, y=y_vals, log = "x", type ="l", col = c(rep(1, length(flu_idx)), rep(2, length(flu_idx_alt))), 
        main = "Hydroxyflutamide", xlab = "Conc", ylab = "Resp")
legend("topleft", legend= c("BLA", "LUC"), lty =c(1,1), col=c(1,2))
df$Curve.Class[flu_idx]
df$ReplicateSet[flu_idx]
df_alt$Curve.Class[flu_idx_alt]
df_alt$ReplicateSet[flu_idx_alt]




### GGPLOT version
xpt1 = cbind(melt(t(df[flu_idx, conc_idx])),"chem"= "Progesterone (57-83-0)")
xpt2 = cbind(melt(t(df[flu_idx_alt, conc_idx])), "chem" = "Fluoxymestrone (76-43-7)")
x_vals = rbind(xpt1, xpt2)
ypt1 = cbind(melt(t(df[flu_idx, resp_idx])))
ypt2 = cbind(melt(t(df[flu_idx_alt, resp_idx])))
y_vals = rbind(ypt1, ypt2)

#pdf(file = "replicate_ex.pdf", onefile= T, width=6, height = 4)
gg_df = cbind(x_vals, y_vals)[,c(4,2,3,7)]
names(gg_df) = c("chem","idx","Concentration", "Response")
ggplot(data = gg_df, aes(Concentration, Response, col = as.factor(chem), group = idx ))+
  geom_point()+geom_line()+scale_x_continuous(trans = "log10")+
  theme(legend.position = "top")+
  labs(colour='Tox21 AR-luc Replicates')
#dev.off()

# Questions  --------------
# Difference between assays?  ex: progesterone
# high variability between replicates ex: Fluoxymestrone
# max dose of mixture has higher concentration than max dose of indivd?

# Answers ------ 
# BLA has better specificity (correctly rejecting if not agonist), 
# LUC has better sensitivity (correctly identifying if an agonist)

# Create data matrix with relevant chemicals to test with DPMCMC ####
colnames(pure_df)
Cx_T21 = pure_df[,6:20]
y_i_T21 = pure_df[,21:35]



