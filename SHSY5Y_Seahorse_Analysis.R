
library(readxl)
library(data.table)
library(scales)


SHSY5Y_4.10_analysis= read_xlsx("SHSY5Y_41026_raw&base2_data.xlsx") #I used SHSY5Y Replicate 1(4.10) base2 correction
head(SHSY5Y_4.10_analysis)
SHSY5Y_4.10_analysis=as.data.table(SHSY5Y_4.10_analysis)
unique(SHSY5Y_4.10_analysis$Group)
SHSY5Y_4.10_analysis_28_MOI1= SHSY5Y_4.10_analysis[grep("^LRRC28 MOI<1$",Group)]
SHSY5Y_4.10_analysis_28_MOI10= SHSY5Y_4.10_analysis[grep("^LRRC28 MOI10$",Group)]
SHSY5Y_4.10_analysis_42_MOI1= SHSY5Y_4.10_analysis[grep("^LRRC42 MOI<1$", Group)]
SHSY5Y_4.10_analysis_42_MOI10= SHSY5Y_4.10_analysis[grep("^LRRC42 MOI10$", Group)]
SHSY5Y_4.10_analysis_NTD_MOI1= SHSY5Y_4.10_analysis[grep("^LRRCNTD MOI<1$", Group)]
SHSY5Y_4.10_analysis_NTD_MOI10= SHSY5Y_4.10_analysis[grep("^LRRCNTD MOI10$", Group)]
SHSY5Y_4.10_analysis_Cont= SHSY5Y_4.10_analysis[grep("Cont-SHSY5Y",Group)]
SHSY5Y_4.10_analysis_bind=rbind(SHSY5Y_4.10_analysis_28_MOI1,SHSY5Y_4.10_analysis_28_MOI10,SHSY5Y_4.10_analysis_42_MOI1,SHSY5Y_4.10_analysis_42_MOI10,SHSY5Y_4.10_analysis_NTD_MOI1,SHSY5Y_4.10_analysis_NTD_MOI10,SHSY5Y_4.10_analysis_Cont)
dim(SHSY5Y_4.10_analysis_bind)
head(SHSY5Y_4.10_analysis_bind)

#compare 28 and 42 base2
SHSY5Y_MOI1_28and42_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[16:20,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[106:110,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28and42_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[21:25,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[111:115,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28and42_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[26:30,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[116:120,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)

#compare 28 and control base2

SHSY5Y_MOI1_28andcont_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[16:20,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andcont_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[21:25,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andcont_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[26:30,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)


#compare 28 and dNTD base2
SHSY5Y_MOI1_28andNTD_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[16:20,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[208:213,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andNTD_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[21:25,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[214:219,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andNTD_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[26:30,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[220:225,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)


#compare cont and dNTD base2
SHSY5Y_MOI1_contandNTD_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[208:213,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contandNTD_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[214:219,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contandNTD_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[220:225,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)

#compare cont and 42 base2
SHSY5Y_MOI1_contand42_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[106:110,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contand42_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[111:115,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contand42_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[116:120,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)

#compare 28 and 42 base2
SHSY5Y_MOI10_28and42_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[61:65,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[154:159,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28and42_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[66:70,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[160:165,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28and42_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[71:75,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[166:171,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)

#compare 28 and control base2

SHSY5Y_MOI10_28andcont_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[61:65,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andcont_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[66:70,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andcont_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[71:75,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)


#compare 28 and dNTD base2
SHSY5Y_MOI10_28andNTD_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[61:65,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[262:267,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andNTD_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[66:70,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[268:273,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andNTD_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[71:75,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[274:279,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)


#compare cont and dNTD base2
SHSY5Y_MOI10_contandNTD_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[262:267,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contandNTD_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[268:273,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contandNTD_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[274:279,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)

#compare cont and 42 base2
SHSY5Y_MOI10_contand42_4.10testbase2_msmt4=t.test(SHSY5Y_4.10_analysis_bind[313:317,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[154:159,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contand42_4.10testbase2_msnt5=t.test(SHSY5Y_4.10_analysis_bind[318:322,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[160:165,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contand42_4.10testbase2_msnt6=t.test(SHSY5Y_4.10_analysis_bind[323:327,]$ECAR_baselined2, SHSY5Y_4.10_analysis_bind[166:171,]$ECAR_baselined2, alternative = "two.sided",var.equal = TRUE)


#SHSY5Y 4.10 base2 normalizaation t.test
all_shsy5y_base2_4.10 <- ls(pattern = "4.10testbase2")# Retrieve the objects and put them in a list
shsy5y_4.10_base2ttest_list <- mget(all_shsy5y_base2_4.10)#  Create a container for the names and p-values
shsy5y_4.10base2_test_names <- names(shsy5y_4.10_base2ttest_list)
shsy5y_4.10base2_p_values <- numeric(length(shsy5y_4.10_base2ttest_list))#  Run a loop to extract the p-values
for (i in 1:length(shsy5y_4.10_base2ttest_list)) {
  shsy5y_4.10base2_p_values[i] <- shsy5y_4.10_base2ttest_list[[i]]$p.value
}

#  Combine them into a clean table
shsy5y_final_4.10base2_p_table <- data.frame(
  Test = shsy5y_4.10base2_test_names,
  P_Value = shsy5y_4.10base2_p_values,
  stringsAsFactors = FALSE
)


shsy5y_final_4.10base2_p_table$Significant <- ifelse(shsy5y_final_4.10base2_p_table$P_Value < 0.05, "Yes", "No")


SHSY5Y_4.10_final_p_table$P_Value_Base2 <- pvalue(SHSY5Y_4.10_final_p_table$P_Value_Base2)


SHSY5Y_4.17_analysis= read_xlsx("SHSY5Y_41726_raw&base1_data.xlsx") 
head(SHSY5Y_4.17_analysis)
SHSY5Y_4.17_analysis=as.data.table(SHSY5Y_4.17_analysis)
unique(SHSY5Y_4.17_analysis$Group)
SHSY5Y_4.17_analysis_28_MOI1= SHSY5Y_4.17_analysis[grep("^LRRC28 MOI<1$",Group)]
SHSY5Y_4.17_analysis_28_4326_MOI10= SHSY5Y_4.17_analysis[grep("^LRRC28 MOI10 4.3.26$",Group)]
SHSY5Y_4.17_analysis_28_3926_MOI10= SHSY5Y_4.17_analysis[grep("^LRRC28 MOI10 3.9.26$",Group)]
SHSY5Y_4.17_analysis_42_MOI1= SHSY5Y_4.17_analysis[grep("^LRRC42 MOI<1$", Group)]
SHSY5Y_4.17_analysis_42_MOI10= SHSY5Y_4.17_analysis[grep("^LRRC42 MOI10$", Group)]
SHSY5Y_4.17_analysis_NTD_MOI1= SHSY5Y_4.17_analysis[grep("^LRRC dNTD MOI<1$", Group)]
SHSY5Y_4.17_analysis_NTD_MOI10= SHSY5Y_4.17_analysis[grep("^LRRC dNTD MOI10$", Group)]
SHSY5Y_4.17_analysis_Cont= SHSY5Y_4.17_analysis[grep("SH-SY5Y Control",Group)]
SHSY5Y_4.17_analysis_bind=rbind(SHSY5Y_4.17_analysis_28_MOI1,SHSY5Y_4.17_analysis_28_4326_MOI10,SHSY5Y_4.17_analysis_28_3926_MOI10,SHSY5Y_4.17_analysis_42_MOI1,SHSY5Y_4.17_analysis_42_MOI10,SHSY5Y_4.17_analysis_NTD_MOI1,SHSY5Y_4.17_analysis_NTD_MOI10,SHSY5Y_4.17_analysis_Cont)
dim(SHSY5Y_4.17_analysis_bind)
head(SHSY5Y_4.17_analysis_bind)


#compare 28 and 42 MOI1 4.17.26 Base1
SHSY5Y_MOI1_28and42_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[19:24,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[127:132,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28and42_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[25:30,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[133:138,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28and42_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[31:36,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[139:144,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)

#compare 28 and control MOI1 4.17.26 Base1

SHSY5Y_MOI1_28andcont_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[19:24,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andcont_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[25:30,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andcont_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[31:36,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)


#compare 28 and dNTD MOI1 4.17.26 Base1
SHSY5Y_MOI1_28andNTD_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[19:24,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[235:240,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andNTD_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[25:30,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[241:246,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_28andNTD_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[31:36,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[247:252,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)


#compare cont and dNTD MOI1 4.17.26 Base1
SHSY5Y_MOI1_contandNTD_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[235:240,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contandNTD_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[241:246,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contandNTD_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[247:252,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)

#compare cont and 42 MOI1 4.17.26 Base1
SHSY5Y_MOI1_contand42_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[127:132,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contand42_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[133:138,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI1_contand42_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[139:144,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)

#compare 28 and 42 MOI10 4.17.26 Base1
SHSY5Y_MOI10_28and42_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[64:66,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[181:186,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28and42_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[67:69,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[187:192,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28and42_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[70:72,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[193:198,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)

#compare 28 and control MOI10 4.17.26 Base1

SHSY5Y_MOI10_28andcont_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[64:66,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andcont_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[67:69,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andcont_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[70:72,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)


#compare 28 and dNTD MOI10 4.17.26 Base1
SHSY5Y_MOI10_28andNTD_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[64:66,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[289:294,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andNTD_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[67:69,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[295:300,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_28andNTD_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[70:72,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[301:306,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)


#compare cont and dNTD MOI10 4.17.26 Base1
SHSY5Y_MOI10_contandNTD_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[289:294,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contandNTD_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[295:300,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contandNTD_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[301:306,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)

#compare cont and 42 MOI10 4.17.26 Base1
SHSY5Y_MOI10_contand42_4.17testbase1_msmt4=t.test(SHSY5Y_4.17_analysis_bind[343:348,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[181:186,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contand42_4.17testbase1_msnt5=t.test(SHSY5Y_4.17_analysis_bind[349:354,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[187:192,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)
SHSY5Y_MOI10_contand42_4.17testbase1_msnt6=t.test(SHSY5Y_4.17_analysis_bind[355:360,]$ECAR_base1, SHSY5Y_4.17_analysis_bind[193:198,]$ECAR_base1, alternative = "two.sided",var.equal = TRUE)


all_shsy5y_test_4.17_base1 <- ls(pattern = "4.17testbase1")  # Retrieve the objects and put them in a list
shsy5y_4.17_list_base1 <- mget(all_shsy5y_test_4.17_base1)  # Create a container for the names and p-values
shsy5y_4.17_test_names_base1 <- names(shsy5y_4.17_list_base1)
shsy5y_4.17_p_values_base1 <- numeric(length(shsy5y_4.17_list_base1)) #Run a loop to extract the p-values
for (i in 1:length(shsy5y_4.17_list_base1)) {
  shsy5y_4.17_p_values_base1[i] <- shsy5y_4.17_list_base1[[i]]$p.value
}


final_shsy5y_4.17_p_table_base1 <- data.frame(
  Test = shsy5y_4.17_test_names_base1,
  P_Value = shsy5y_4.17_p_values_base1,
  stringsAsFactors = FALSE
)


final_shsy5y_4.17_p_table_base1$Significant <- ifelse(final_shsy5y_4.17_p_table_base1$P_Value < 0.05, "Yes", "No")

