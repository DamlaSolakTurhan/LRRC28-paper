# seahorse experiment.
library(readxl)
library(data.table)
library(scales)


# Calculate the Unpaired two-sided Student’s t-test at each time point for Replicate 1
Huh_3.4_rawdata= read_xlsx("Huh7.5_3426_raw.xlsx") 
head(Huh_3.4_rawdata)
Huh_3.4_rawdata=as.data.table(Huh_3.4_rawdata)
Huh_3.4_rawdata_28= Huh_3.4_rawdata[grep("^LRRC28$",Group)]
Huh_3.4_rawdata_42= Huh_3.4_rawdata[grep("^LRRC42$", Group)]
Huh_3.4_rawdata_Cont= Huh_3.4_rawdata[grep("Huh 7.5",Group)]
Huh_3.4_rawdata=rbind(Huh_3.4_rawdata_28, Huh_3.4_rawdata_42,Huh_3.4_rawdata_Cont)
dim(Huh_3.4_rawdata)


#compare 28 and 42 baseline 1

Huh28and42_3.4test_base1_msmt4=t.test(Huh_3.4_rawdata[19:24,]$`ECAR Base1`, Huh_3.4_rawdata[73:78,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_3.4test_base1_msnt5=t.test(Huh_3.4_rawdata[25:30,]$`ECAR Base1`, Huh_3.4_rawdata[79:84,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_3.4test_base1_msnt6=t.test(Huh_3.4_rawdata[31:36,]$`ECAR Base1`, Huh_3.4_rawdata[85:90,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_3.4test_base1_msnt7=t.test(Huh_3.4_rawdata[37:42,]$`ECAR Base1`, Huh_3.4_rawdata[91:96,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
#compare 28 and control baseline 1

Huh28cont_3.4test_base1_msnt4=t.test(Huh_3.4_rawdata[19:24,]$`ECAR Base1`, Huh_3.4_rawdata[127:132,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_3.4test_base1_msnt5=t.test(Huh_3.4_rawdata[25:30,]$`ECAR Base1`, Huh_3.4_rawdata[133:138,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_3.4test_base1_msnt6=t.test(Huh_3.4_rawdata[31:36,]$`ECAR Base1`, Huh_3.4_rawdata[139:144,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_3.4test_base1_msnt7=t.test(Huh_3.4_rawdata[37:42,]$`ECAR Base1`, Huh_3.4_rawdata[145:150,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)


#compare 42 and control baseline 1
Huh42cont_3.4test_base1_msnt4=t.test(Huh_3.4_rawdata[73:78,]$`ECAR Base1`, Huh_3.4_rawdata[127:132,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_3.4test_base1_msnt5=t.test(Huh_3.4_rawdata[79:84,]$`ECAR Base1`, Huh_3.4_rawdata[133:138,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_3.4test_base1_msnt6=t.test(Huh_3.4_rawdata[85:90,]$`ECAR Base1`, Huh_3.4_rawdata[139:144,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_3.4test_base1_msnt7=t.test(Huh_3.4_rawdata[91:96,]$`ECAR Base1`, Huh_3.4_rawdata[145:150,]$`ECAR Base1`, alternative = "two.sided", var.equal = TRUE)


#Huh 3.4 base1 norm t.test table

all_huh_Base1_3.4 <- ls(pattern = "base1")# Retrieve the objects and put them in a list
huh_3.4_Base1_list <- mget(all_huh_Base1_3.4) # Create a container for the names and p-values
test_names <- names(huh_3.4_Base1_list)
Base1_p_values <- numeric(length(huh_3.4_Base1_list))# Run a loop to extract the p-values
for (i in 1:length(huh_3.4_Base1_list)) {
  Base1_p_values[i] <- huh_3.4_Base1_list[[i]]$p.value
}

# Combine them into a clean table
final_Base1_p_table <- data.frame(
  Test = test_names,
  P_Value = Base1_p_values,
  stringsAsFactors = FALSE
)

# Add a column that says "Yes" if p < 0.05 and "No" otherwise
final_Base1_p_table$Significant <- ifelse(final_Base1_p_table$P_Value < 0.05, "Yes", "No")


Huh_34_final_p_table$P_Value_Base1 <- pvalue(Huh_34_final_p_table$P_Value_Base1)


#Huh-7.5 Replicate 2 - Unpaired two-sided Student’s t-test

Huh_313_rawdata= read_xlsx("Huh7.5_31326_raw.xlsx")
head(Huh_313_rawdata)
Huh_313_rawdata=as.data.table(Huh_313_rawdata)
Huh_313_rawdata_28= Huh_313_rawdata[grep("^LRRC28$",Group)]
Huh_313_rawdata_42= Huh_313_rawdata[grep("^LRRC42$", Group)]
Huh_313_rawdata_Cont= Huh_313_rawdata[grep("Huh 7.5 Control",Group)]
Huh_313_rawdata=rbind(Huh_313_rawdata_28, Huh_313_rawdata_42,Huh_313_rawdata_Cont)
dim(Huh_313_rawdata)


# to compare between 28 and 42 base2
Huh28and42_313test_base2_msnt4=t.test(Huh_313_rawdata[19:24,]$`ECAR Base2`, Huh_313_rawdata[73:78,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_313test_base2_msnt5=t.test(Huh_313_rawdata[25:30,]$`ECAR Base2`, Huh_313_rawdata[79:84,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_313test_base2_msnt6=t.test(Huh_313_rawdata[31:36,]$`ECAR Base2`, Huh_313_rawdata[85:90,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28and42_313test_base2_msnt7=t.test(Huh_313_rawdata[37:42,]$`ECAR Base2`, Huh_313_rawdata[91:96,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)

#compare 28 and control base2

Huh28cont_313test_base2_msnt4=t.test(Huh_313_rawdata[19:24,]$`ECAR Base2`, Huh_313_rawdata[127:132,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_313test_base2_msnt5=t.test(Huh_313_rawdata[25:30,]$`ECAR Base2`, Huh_313_rawdata[133:138,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_313test_base2_msnt6=t.test(Huh_313_rawdata[31:36,]$`ECAR Base2`, Huh_313_rawdata[139:144,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh28cont_313test_base2_msnt7=t.test(Huh_313_rawdata[37:42,]$`ECAR Base2`, Huh_313_rawdata[145:150,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)


#compare 42 and control base2
Huh42cont_313test_base2_msnt4=t.test(Huh_313_rawdata[73:78,]$`ECAR Base2`, Huh_313_rawdata[127:132,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_313test_base2_msnt5=t.test(Huh_313_rawdata[79:84,]$`ECAR Base2`, Huh_313_rawdata[133:138,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_313test_base2_msnt6=t.test(Huh_313_rawdata[85:90,]$`ECAR Base2`, Huh_313_rawdata[139:144,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)
Huh42cont_313test_base2_msnt7=t.test(Huh_313_rawdata[91:96,]$`ECAR Base2`, Huh_313_rawdata[145:150,]$`ECAR Base2`, alternative = "two.sided", var.equal = TRUE)



all_huh_test_3.13_base2 <- ls(pattern = "313test_base2")# Retrieve the objects and put them in a list
huh_3.13_list_base2 <- mget(all_huh_test_3.13_base2)# Create a container for the names and p-values
test_names_base2 <- names(huh_3.13_list_base2)
p_values_base2 <- numeric(length(huh_3.13_list_base2))# Run a loop to extract the p-values
for (i in 1:length(huh_3.13_list_base2)) {
  p_values_base2[i] <- huh_3.13_list_base2[[i]]$p.value
}

Huh_313_final_p_table <- data.frame(
  Test = test_names_raw,
   P_Value_Base2 = p_values_base2,
   stringsAsFactors = FALSE
)

Huh_313_final_p_table$P_Value_Base2 <- pvalue(Huh_313_final_p_table$P_Value_Base2)

Huh_34_final_p_table
Huh_313_final_p_table



