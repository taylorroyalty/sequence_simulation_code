#load in fit object and predict sequencing effort (Figure 5)
genome_size<-1e6 #in base pairs
target_fraction<-0.5 #desired fraction of the target genome (values: 0 to 1)
metagenome_prop<-0.5 #proportion of metagenome of the target genome (values: 0 to 1)

load(file = "GRASE_FIT_OBJECT.rda")

#note probability is the legacy term for proportion of metagenome variable used above
df_tmp<-data.frame(genome_size=genome_size,target=target_fraction,probability=metagenome_prop)
pred<-predict.gam(fit_object, newdata = df_tmp)
print(pred*100) #multiply by 100 due to simulation being ran with 100 bp read length; converts to number of bases needed to be sequeced
