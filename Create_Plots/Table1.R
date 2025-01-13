
meta_data <- read.csv("Pipeline_All/results_TP1/metadata.filtered.csv")
length(unique(meta_data$ADRC.study.ID))

p_val <- lapply(unique(meta_data$Diagnosis), function(d) {
  print(d)
  res <- t.test(meta_data$Age[meta_data$Diagnosis == d & meta_data$Sex == "female"], meta_data$Age[meta_data$Diagnosis == d & meta_data$Sex == "male"])
  res$p.value

})
names(p_val) <- unique(meta_data$Diagnosis)
p.adjust(p_val)


p_val <- lapply(c("White", "Asian"), function(d) {
   print(d)
   res <- t.test(meta_data$Age[meta_data$Race == d & meta_data$Sex == "female"], meta_data$Age[meta_data$Race == d & meta_data$Sex == "male"])
   res$p.value
} )

res <- t.test(meta_data$Age[(!meta_data$Race %in% c("White", "Asian"))& meta_data$Sex == "female"], meta_data$Age[(!meta_data$Race %in% c("White", "Asian")) & meta_data$Sex == "male"])
res$p.value

p_val <- c(p_val, res$p.value)

names(p_val) <- c("White", "Asian", "Other")
p.adjust(p_val)



p_val_w <- lapply(unique(meta_data$Diagnosis), function(d) {
  print(d)
  res <- t.test(meta_data$Age[meta_data$Race == "White" & meta_data$Diagnosis == d & meta_data$Sex == "female"], meta_data$Age[meta_data$Race == "White" & meta_data$Diagnosis == d & meta_data$Sex == "male"])
  res$p.value
  
})
names(p_val_w) <- unique(meta_data$Diagnosis)

p_val_a <- lapply(unique(meta_data$Diagnosis), function(d) {
  print(d)
  if (d == "Parkinson's Disease with MCI" | d == "Alzheimer's disease" | d == "Parkinson's Disease only") return(2)
  res <- t.test(meta_data$Age[meta_data$Race == "Asian" & meta_data$Diagnosis == d & meta_data$Sex == "female"], meta_data$Age[meta_data$Race == "Asian" & meta_data$Diagnosis == d & meta_data$Sex == "male"])
  res$p.value
  
})
names(p_val_a) <- unique(meta_data$Diagnosis)



p_val_o <- lapply(unique(meta_data$Diagnosis), function(d) {
  print(d)
  res <- t.test(meta_data$Age[(!meta_data$Race %in% c("White", "Asian")) & meta_data$Diagnosis == d & meta_data$Sex == "female"], meta_data$Age[(!meta_data$Race %in% c("White", "Asian"))& meta_data$Diagnosis == d & meta_data$Sex == "male"])
  res$p.value
  
})
names(p_val_o) <- unique(meta_data$Diagnosis)

p <- c(p_val_w[p_val_w != 2], p_val_a[p_val_a != 2])
p.adjust(p, n = 15)
