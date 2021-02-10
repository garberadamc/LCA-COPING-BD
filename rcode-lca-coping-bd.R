## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

## TITLE: "Appendix - Implement LCA Analyses with MplusAutomation"

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("rhdf5")


## -----------------------------------------------------------------------------------

## library(extrafont)
## loadfonts()

library(MplusAutomation)
library(tidyverse)
library(rhdf5)
library(here)
library(glue)
library(gt)
library(reshape2)
library(cowplot)
library(patchwork)


## -----------------------------------------------------------------------------------

ies_data <- read_csv(here("data","ies_bd_data.csv"),
                     na=c("",".", "999"))


## -----------------------------------------------------------------------------------

##  1.1 Enumeration: 

## -----------------------------------------------------------------------------------

lca_k1_6  <- lapply(1:6, function(k) {
  lca_enum  <- mplusObject(

  TITLE = glue("Class-{k} LCA Enumeration - Youth Coping Strategies"),

  VARIABLE = glue(
    "categorical = do1 do2 do3 do5 do6; !!! Coping strategy items for measurement model !!!
     usevar = do1 do2 do3 do5 do6;

     classes = c({k});"),

  ANALYSIS =
    "estimator = mlr;
     type = mixture;
     processors=10;",

  PLOT =
    "type = plot3;
     series = do1 do2 do3 do5 do6(*);",

  OUTPUT = "sampstat tech11 tech14;",

  usevariables = colnames(ies_data),
  rdata = ies_data)

lca_enum_fit <- mplusModeler(lca_enum,
                  dataout=glue(here("enum_mplus", "c_lca_enum_bd.dat")),
                  modelout=glue(here("enum_mplus", "c{k}_lca_enum_bd.inp")) ,
                  check=TRUE, run = TRUE, hashfilename = FALSE)
})



## -----------------------------------------------------------------------------------

##  1.2 Generate Model Fit Summary Table 

## -----------------------------------------------------------------------------------

output_enum <- readModels(here("enum_mplus"), quiet = TRUE)

enum_extract <- LatexSummaryTable(output_enum,                                 
                keepCols=c("Title", "Parameters", "LL", "BIC", "aBIC",
                           "BLRT_PValue", "T11_VLMR_PValue","Observations"))


## -----------------------------------------------------------------------------------

allFit <- enum_extract %>% 
  mutate(aBIC = -2*LL+Parameters*log((Observations+2)/24)) %>% 
  mutate(CIAC = -2*LL+Parameters*(log(Observations)+1)) %>% 
  mutate(AWE = -2*LL+2*Parameters*(log(Observations)+1.5)) %>%
  mutate(SIC = -.5*BIC) %>% 
  mutate(expSIC = exp(SIC - max(SIC))) %>% 
  mutate(BF = exp(SIC-lead(SIC))) %>% 
  mutate(cmPk = expSIC/sum(expSIC)) %>% 
  select(1:5,9:10,6:7,13,14) %>% 
  arrange(Parameters)


## -----------------------------------------------------------------------------------

allFit %>% 
  mutate(Title = str_remove(Title, " LCA Enumeration - Youth Coping Strategies")) %>% 
  gt() %>%
  tab_header(
    title = md("**Model Fit Summary Table**"), subtitle = md("&nbsp;")) %>% 
  cols_label(
    Title = "Classes",
    Parameters = md("Par"),
    LL = md("*LL*"),
    T11_VLMR_PValue = "VLMR",
    BLRT_PValue = "BLRT",
    BF = md("BF"),
    cmPk = md("*cmP_k*")) %>%
  tab_footnote(
    footnote = md(
    "*Note.* Par = parameters; *LL* = log likelihood;
      BIC = bayesian information criterion;
      aBIC = sample size adjusted BIC; CAIC = consistent Akaike information criterion;
      AWE = approximate weight of evidence criterion;
      BLRT = bootstrapped likelihood ratio test p-value;
      VLMR = Vuong-Lo-Mendell-Rubin adjusted likelihood ratio test p-value;
      cmPk = approximate correct model probability."), 
    locations = cells_title()) %>% 
  tab_options(column_labels.font.weight = "bold") %>% 
  fmt_number(10,decimals = 2,
             drop_trailing_zeros=TRUE,
             suffixing = TRUE) %>% 
  fmt_number(c(3:9,11), decimals = 0) %>% 
  fmt_missing(1:11, missing_text = "--") %>% 
  fmt(c(8:9,11), fns = function(x) 
    ifelse(x<0.001, "<.001",
           scales::number(x, accuracy = 0.01))) %>%
  fmt(10, fns = function(x) 
    ifelse(x>100, ">100",
           scales::number(x, accuracy = .1))) 

# --------------------------------------------------------------------------------------

##  1.3 Plot Information Criteria 

# --------------------------------------------------------------------------------------

allFit %>% select(2:7) %>% 
  rowid_to_column() %>% 
  pivot_longer(`BIC`:`AWE`, 
    names_to = "Index", 
    values_to = "ic_value") %>% 
  mutate(Index = factor(Index,
    levels = c("AWE","CIAC","BIC","aBIC"))) %>%
  ggplot(aes(x = rowid, y = ic_value,
    color = Index, shape = Index,
    group = Index, lty = Index)) + 
  geom_point(size = 2.0) + geom_line(size = .8) +
  scale_x_continuous(breaks = 1:6) +
  scale_colour_grey(end = .5) +
  theme_cowplot() + 
  labs(x = "Number of Classes", y = "Information Criteria Value") +
  theme(legend.title = element_blank(), 
        legend.position = "top")


ggsave(here("figures", "Fig2_IC_plot.png"), dpi=300, height=5, width=7, units="in")


# --------------------------------------------------------------------------------------

##  1.4 Compare Conditional Item Probability Plots

# --------------------------------------------------------------------------------------

model_results <- data.frame()

for (i in 1:length(output_enum)) {
  temp <- output_enum[[i]]$parameters$unstandardized
  temp <- data.frame(unclass(temp)) %>%
    mutate(model = paste0(i, "-Class Model")) 
  model_results <- rbind(model_results, temp)
}

model_results <- model_results %>% 
  filter(paramHeader == "Thresholds") %>% 
  select(est, model, LatentClass, param) %>% 
  mutate(prob = (1 / (1 + exp(est)))) 

ggplot(model_results,
       aes(x = param, y = prob,
           color = LatentClass,
           shape = LatentClass,
           group = LatentClass)) + 
  geom_point() + geom_line() +
  facet_wrap(~ model, ncol = 2) +
  labs(title = "LCA Posterior Probability Plot",
       x= "Coping Strategy Items", y = "Probability") +
  theme_minimal()


## -----------------------------------------------------------------------------------

ggsave(here("figures","FigA1_compare_kclass_LCAs.png"),
       dpi=300, height=4, width=6, units="in")

# --------------------------------------------------------------------------------------

##  1.5 Plot Final Model - Conditional Item Probability Plot 

# --------------------------------------------------------------------------------------

model_step1 <- readModels(here("3step_mplus", "Step1_3step_BD.out"), quiet = TRUE)
                           
## -----------------------------------------------------------------------------------

plot_lca_function <- function(model_name,item_num,class_num,item_labels,
                              class_labels,class_legend_order,plot_title){

mplus_model <- as.data.frame(model_name$gh5$means_and_variances_data$estimated_probs$values)
plot_data <- mplus_model[seq(2, 2*item_num, 2),]

c_size <- as.data.frame(model_name$class_counts$modelEstimated$proportion)
colnames(c_size) <- paste0("cs")
c_size <- c_size %>% mutate(cs = round(cs*100, 2))
colnames(plot_data) <- paste0(class_labels, glue(" ({c_size[1:class_num,]}%)"))
plot_data <- plot_data %>% relocate(class_legend_order)

plot_data <- cbind(Var = paste0("U", 1:item_num), plot_data)
plot_data$Var <- factor(plot_data$Var,
               labels = item_labels)
plot_data$Var <- fct_inorder(plot_data$Var)

pd_long_data <- melt(plot_data, id.vars = "Var") 

# This syntax uses the date-frame created above to produce the plot with `ggplot()`

p <- pd_long_data %>%
  ggplot(aes(x = as.integer(Var), y = value,
  shape = variable, colour = variable, lty = variable)) +
  geom_point(size = 4) + geom_line() + 
  scale_x_continuous("", breaks = 1:5, labels = plot_data$Var) + 
  scale_colour_grey() + 
  labs(title = plot_title, y = "Probability") +
  theme_cowplot() +
  theme(legend.title = element_blank(), 
        legend.position = "top")

p
return(p)
}

## -----------------------------------------------------------------------------------

plot_lca_function(
  model_name = model_step1, 
  item_num = 5,
  class_num = 3,
  item_labels = c("Talked to Someone","Got Mad", "Avoidance","Problem Solving", "Worried/Cried"),
  class_labels = c("Adaptive Coping","Externalizing Behavior","No Coping"),
  class_legend_order = c(2,1,3),
  plot_title = "LCA Posterior Probability Plot"
  )

## -----------------------------------------------------------------------------------

ggsave(here("figures", "Fig1_LCA_3-Class.png"), dpi=300, height=5, width=7, units="in")

# --------------------------------------------------------------------------------------

##  1.6 Manual "3-Step" ML Auxiliary Variable Integration Method

# --------------------------------------------------------------------------------------

### Step 1 - Estimate the unconditional model with all covariate & distal outcome variables mentioned in the `auxiliary` statement. 

m_step1  <- mplusObject(
  TITLE = "Step1_3step_automation Behavioral Disorder",
  VARIABLE =
   "categorical = do1 do2 do3 do5 do6;

    usevar = do1 do2 do3 do5 do6;

    classes = c(3);

    !!! NOTE: All auxiliary variables to be considered in the final model should be listed here !!!
    auxiliary =
    FEMALE ETHN_CMP SOC_STRS
    BOTHR_U negmood1 posmood1;",

  ANALYSIS =
   "estimator = mlr;
    type = mixture;
    starts = 500 100;",

  SAVEDATA =
   "!!! NOTE: This saved dataset will contain class probabilities and modal assignment columns !!!
    File=3step_BD_savedata_012020.dat;
    Save=cprob;
    Missflag= 999;",

  MODEL = "",
  OUTPUT = "",

  PLOT =
    "type = plot3;
    series = do1 do2 do3 do5 do6(*);",

  usevariables = colnames(ies_data),
  rdata = ies_data)

m_step1_fit <- mplusModeler(m_step1,
                 dataout=here("3step_mplus", "Step1_3step_BD.dat"),
                 modelout=here("3step_mplus", "Step1_3step_BD.inp") ,
                 check=TRUE, run = TRUE, hashfilename = FALSE)


## -----------------------------------------------------------------------------------

### Step 2 - Extract logits & saved data from the step 1 unconditional model.

logit_cprobs <- as.data.frame(m_step1_fit[["results"]]
                                         [["class_counts"]]
                                         [["logitProbs.mostLikely"]])

## -----------------------------------------------------------------------------------

savedata <- as.data.frame(m_step1_fit[["results"]]
                                     [["savedata"]])

## -----------------------------------------------------------------------------------

colnames(savedata)[colnames(savedata)=="C"] <- "N"

## -----------------------------------------------------------------------------------

### Step 3 (part 1) - Estimate the unconditional model with logits from step 2. 

m_step2  <- mplusObject(
  TITLE = "Step2_3step_automation Behavioral Disorder",

  VARIABLE =
 "nominal=N;
  USEVAR = n;
  missing are all (999);
  classes = c(3); ",

  ANALYSIS =
 "estimator = mlr;
  type = mixture;
  starts = 0;",

  MODEL =
    glue(
 "%C#1%
  [n#1@{logit_cprobs[1,1]}];
  [n#2@{logit_cprobs[1,2]}];

  %C#2%
  [n#1@{logit_cprobs[2,1]}];
  [n#2@{logit_cprobs[2,2]}];

  %C#3%
  [n#1@{logit_cprobs[3,1]}];
  [n#2@{logit_cprobs[3,2]}];"),

  OUTPUT = "!tech11  tech14 res;",

  PLOT =
 "!type = plot3;
  !series = do1 do2 do3 do5 do6(*);",

  usevariables = colnames(savedata),
  rdata = savedata)

m_step2_fit <- mplusModeler(m_step2,
                 dataout=here("3step_mplus", "Step2_3step_BD.dat"),
                 modelout=here("3step_mplus", "Step2_3step_BD.inp"),
                 check=TRUE, run = TRUE, hashfilename = FALSE)


# --------------------------------------------------------------------------------------

### Step 3 (part 2) - Add covariates & distal outcomes to the model. 

##  1.7 Moderation - Estimate the final SEM Model

# --------------------------------------------------------------------------------------

m_step3  <- mplusObject(
  TITLE = "Step3_3step_automation Behavioral Disorder",

  VARIABLE =
 "nominal = N;
  usevar = n;
  missing are all (999);

  usevar = SOC_STRS POSMOOD1 NEGMOOD1;
  classes = c(3); ",

  DEFINE =
 "Center SOC_STRS (Grandmean);",

  ANALYSIS =
 "estimator = mlr;
  type = mixture;
  starts = 0;",

  MODEL =
  glue(
 "!DISTAL = POSMOOD1 NEGMOOD1
  !MODERATOR = SOC_STRS

  %OVERALL%
  POSMOOD1 on SOC_STRS;
  POSMOOD1;

  NEGMOOD1 on SOC_STRS;
  NEGMOOD1;

  %C#1%
  [n#1@{logit_cprobs[1,1]}];
  [n#2@{logit_cprobs[1,2]}];

  [NEGMOOD1](m01);
  NEGMOOD1;                      !!! estimate conditional intercept !!!
  NEGMOOD1 on SOC_STRS (s01);    !!! estimate conditional regression !!!

  [POSMOOD1] (m1);
  POSMOOD1;
  POSMOOD1 on SOC_STRS (s1);

  %C#2%
  [n#1@{logit_cprobs[2,1]}];
  [n#2@{logit_cprobs[2,2]}];

  [NEGMOOD1](m02);
  NEGMOOD1;
  NEGMOOD1 on SOC_STRS (s02);

  [POSMOOD1] (m2);
  POSMOOD1;
  POSMOOD1 on SOC_STRS (s2);

  %C#3%
  [n#1@{logit_cprobs[3,1]}];
  [n#2@{logit_cprobs[3,2]}];

  [NEGMOOD1](m03);
  NEGMOOD1;
  NEGMOOD1 on SOC_STRS (s03);

  [POSMOOD1] (m3);
  POSMOOD1;
  POSMOOD1 on SOC_STRS (s3);"),

  MODELCONSTRAINT =
 "New (diff12 diff13
  diff23 slope12 slope13
  slope23 ndiff12 ndiff13
  ndiff23 nslope12 nslope13
  nslope23);

  diff12 = m1-m2;   ndiff12 = m01-m02;
  diff13 = m1-m3;   ndiff13 = m01-m03;
  diff23 = m2-m3;   ndiff23 = m02-m03;
  slope12 = s1-s2;  nslope12 = s01-s02;
  slope13 = s1-s3;  nslope13 = s01-s03;
  slope23 = s2-s3;  nslope23 = s02-s03;",

  MODELTEST =
  ## NOTE: Only a single Wald test can be conducted per model run. Therefore,
  ## this example requires running separate models for each omnibus test (e.g.,
  ## 4 models; 2 outcomes and 2 slope coefficients). This can be done by
  ## commenting out all but one test and then making multiple input/output files.

 "m1=m2;       !!! Distal outcome omnibus Wald test for `POSMOOD1` !!!
  m2=m3;

  !s1=s2;      !!! Slope difference omnibus Wald test `POSMOOD1 on SOC_STRS` !!!
  !s2=s3;

  !m01=m02;    !!! Distal outcome omnibus Wald test for `NEGMOOD1` !!!
  !m02=m03;

  !s01=s02;   !!! Slope difference omnibus Wald test for `POSMOOD1 on SOC_STRS` !!!
  !s02=s03;",

  usevariables = colnames(savedata),
  rdata = savedata)

m_step3_fit <- mplusModeler(m_step3,
                 dataout=here("3step_mplus", "Step3_3step_BD.dat"),
                 modelout=here("3step_mplus", "Step3_3step_BD.inp"),
                 check=TRUE, run = TRUE, hashfilename = FALSE)


# --------------------------------------------------------------------------------------

## Estimate step 3 model with covariate un-centered for simple-slopes plots. 

m_uncen <- update(m_step3,
  DEFINE = ~" ") # This update removes the centering syntax from the model object `m_step3`

m_uncen_fit <- mplusModeler(m_uncen,
   dataout=here("3step_mplus", "Step3_3step_BD.dat"),
   modelout=here("3step_mplus", "Step3_uncentered_BD.inp"),
   check=TRUE, run = TRUE, hashfilename = FALSE)UCSB_Navy_mark.png


# --------------------------------------------------------------------------------------

## 1.8 Distal Outcome Plot

# --------------------------------------------------------------------------------------

model_step3 <- readModels(here("3step_mplus", "Step3_3step_BD.out"), quiet = TRUE)
  
model_step3 <- data.frame(model_step3$parameters$unstandardized)

## -----------------------------------------------------------------------------------

distal_data <- model_step3 %>% 
  filter(paramHeader == "Intercepts") %>% 
  mutate(param = case_when(
      param == "POSMOOD1" ~ "Positive Mood",
      param == "NEGMOOD1" ~ "Negative Mood")) %>% 
  mutate(LatentClass = factor(LatentClass,
      labels = c("Adaptive Coping", "Externalizing Behavior","No Coping"))) %>% 
  mutate(value_labels = c("3.16 a", "1.82 d", "2.67 b", "2.31 c", "2.72 b", "1.67 d"))
  
## -----------------------------------------------------------------------------------

ggplot(distal_data, 
       aes(fill=LatentClass, y=est, x=fct_rev(param))) + 
  geom_bar(position="dodge", stat="identity", color="black", size=.3) + 
  geom_errorbar( aes(ymin=est-se, ymax=est+se),
                 width=.2, position=position_dodge(.9)) +
  geom_text(aes(y = est -.5, label = value_labels),
            family="Times New Roman", size=4,
            position=position_dodge(.9)) +
  scale_fill_grey(start = 0.6, end = 1.0) +
  theme_cowplot() +
  xlab("") + ylab("Mean Value") +
  theme(text=element_text(family="Times New Roman", size=12),
        axis.text.x=element_text(size=10)) +
        coord_cartesian(expand = FALSE)


ggsave(here("figures","Fig3_distal_barplot.png"), dpi=300, height=4, width=6, units="in")

# --------------------------------------------------------------------------------------

## 1.9 Simple Slope Plots

# --------------------------------------------------------------------------------------


model_uncen <- readModels(here("3step_mplus", "Step3_uncentered_BD.out"), quiet = TRUE)
  
model_uncen <- data.frame(model_uncen$parameters$unstandardized)

slope_data <- model_uncen %>% 
  filter(str_detect(paramHeader, 'ON|Inter')) %>% 
  unite("param", paramHeader:param, remove = TRUE) %>% 
  mutate(param = str_replace(param, "MOOD1.ON_SOC_STRS", "_COEF")) %>% 
  mutate(param = str_remove_all(param, "Intercepts_|MOOD1")) %>% 
  mutate(LatentClass = factor(LatentClass,
    labels = 
    c("Adaptive Coping (63.6%)", "Externalizing Behavior (5.3%)","No Coping (31.1%)")))

## -----------------------------------------------------------------------------------

pos_data <- slope_data %>% 
  filter(str_detect(param, 'POS')) 

pos_wide <- pos_data %>%
  select(param,est, LatentClass) %>% 
  pivot_wider(names_from = param, values_from = est) %>% 
  rename("No.Social.Stress" = 'POS') %>% 
  mutate(Social.Stress = No.Social.Stress + POS_COEF) %>% # calc. condit. means `SOC_STRS = 1`
  select(-POS_COEF)

plot_pos <- melt(pos_wide, id.vars = "LatentClass") %>% 
  mutate(variable = factor(variable,
                      levels = c("No.Social.Stress","Social.Stress"),
                      labels = c("No Social Stress","Social Stress")))

## -----------------------------------------------------------------------------------

p_plot <- ggplot(plot_pos,
            aes(y=value, x=variable,
                color=LatentClass,
                group=LatentClass,
                shape=LatentClass,
                lty=LatentClass)) + 
  geom_point(size = 4) + geom_line() + 
  xlab("") + ylab("Positive Mood")  + ylim(2.5,3.5)+
  scale_colour_grey() +
  theme_classic() + 
  theme(text=element_text(family="Times New Roman", size=12), 
        axis.text.x=element_text(size=12),
        legend.text = element_text(family="Times New Roman", size=10),
        legend.position = "top", legend.title = element_blank()) +
      annotate(geom = "text",
           x = 1.8, y = 2.77,
           label = "N.S.", color = "black") + 
      annotate(geom = "text",
           x = 1.8, y = 2.60,
           label = "N.S.", color = "black")

## -----------------------------------------------------------------------------------

neg_data <- slope_data %>% 
  filter(str_detect(param, 'NEG')) 

neg_wide <- neg_data %>%
  select(param,est, LatentClass) %>% 
  pivot_wider(names_from = param, values_from = est) %>% 
  rename("No.Social.Stress" = 'NEG') %>% 
  mutate(Social.Stress = No.Social.Stress + NEG_COEF) %>%   # calculate means for `SOC_STRS = 1`
  select(-NEG_COEF)

plot_neg <- melt(neg_wide, id.vars = "LatentClass") %>% 
  mutate(variable = factor(variable,
                      levels = c("No.Social.Stress","Social.Stress"),
                      labels = c("No Social Stress","Social Stress")))



## -----------------------------------------------------------------------------------

n_plot <- ggplot(plot_neg,
            aes(y=value, x=variable,
            color=LatentClass,
            group=LatentClass,
            shape=LatentClass,
            lty=LatentClass)) + 
  geom_point(size=4) + geom_line() + 
  xlab("") + ylab("Negative Mood") + ylim(1.5,3)+
  scale_colour_grey() + theme_classic() +
  theme(text=element_text(family="Times New Roman", color = "black", size=12),
        axis.text.x=element_text(size=12),
        legend.text = element_text(family="Times New Roman", size=10),
        legend.position = "top", legend.title = element_blank()) +
    annotate(geom = "text",
       x = 1.8, y = 1.6,
       label = "N.S.", color = "black")



## -----------------------------------------------------------------------------------

p_plot / n_plot


ggsave(here("figures", "Fig5_simple_slopes.png"), dpi=300, height=8.5, width=6.5, units="in")

# --------------------------------------------------------------------------------------

## 2.0 Response Pattern Table 

# To replicate the response pattern table produced in this study (*Table 3*) see the 
# associated `R/MplusAutomation` code example here:
  
### https://garberadamc.github.io/project-site/Lab10.1-response-patterns
  
## Additional examples using `R/MplusAutomation` for a range of analyses can be found here:
  
### https://garberadamc.github.io/project-site
  
# --------------------------------------------------------------------------------------

## 2.1 Plot of the 2-Class Conditional Item Plot for Comparison

# --------------------------------------------------------------------------------------

plot_lca_function(
  model_name = output_enum[[2]],
  item_num = 5,
  class_num = 2,
  item_labels = c("Talked to Someone","Got Mad", "Avoidance","Problem Solving", "Worried/Cried"),
  class_labels = c("Class 1","Class 2"),
  class_legend_order = c(1,2),
  plot_title = ""
  )

## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------



