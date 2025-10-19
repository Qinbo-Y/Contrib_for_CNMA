library(MASS)
library(netmeta)
library(tidyverse)
library(devtools)

###1.Data from Table 1-----------
treat1a <- c("A","B","A","A+B", "A","C","A+C", "A+C")
treat2a <- c("Placebo","Placebo","B","Placebo","C","Placebo","Placebo","C")
studlab <- paste(treat1a, treat2a, sep = ":")
TEa <- c(-0.43, -0.39, -0.12, -0.82, 0.26, -0.51, -0.22, -0.19)
seTEa <- c(0.34, 0.32, 0.33, 0.18, 0.23, 0.46, 0.29, 0.45)
## NMA model
net_a <- netmeta(TEa, seTEa, treat1a, treat2a, studlab, ref = "Placebo", sm='OR', random = F)
# Construct NMA network plot
netgraph(net_a)

######## additive model
nc0 <- netcomb(net_a, inactive = "Placebo")
######## The component estimates for B from CNMA additive model
OR_B <- nc0$Comp.common[2]

# Xa matrix and comparisons in CNMA
reference <- "Placebo"
X_matrix <- nc0$X.matrix
Xa<- as.data.frame(X_matrix)
Xa <- Xa %>%
  mutate(!!reference := apply(., 1, function(row) {  
    if (all(row == 1 | row == 0)) {
      return(-1)
    } else if (all(row == -1 | row == 0)) {
      return(1)
    } else if (all(row %in% c(0, 1, -1))) {
      return(0)
    } else {
      return(NA)
    }
  }))

titles <- colnames(Xa)


Xa$treat1 <- apply(Xa, 1, function(row) {
  indices <- which(row == 1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$treat2 <- apply(Xa, 1, function(row) {
  indices <- which(row == -1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$comparison <- paste(Xa$treat1, Xa$treat2, sep = ":")

# Ha matrix at the component level
TE <- nc0$TE
seTE <- nc0$seTE.adj.common
weight <- 1/(seTE^2)
W <- diag(weight)
C <- nc0$C.matrix
studylab <- nc0$studlab
X_matrix <- nc0$X.matrix
pseudo_inverse <- ginv(t(X_matrix) %*% W %*% X_matrix)
Ha <- pseudo_inverse %*% t(X_matrix) %*% W
rownames(Ha) <- nc0$comps
colnames(Ha) <- Xa$comparison
Ha[,1] <- Ha[,1] + Ha[,8] 
Ha_merged <- Ha[,-8]  

## Construct CNMA network plot
merged_data <- data.frame(TE, net_a$seTE.adj.common, treat1 = Xa$treat1, treat2 = Xa$treat2, studylab)
names(merged_data)[2] <- "seTE"
cnet_a <- netmeta(merged_data$TE, merged_data$seTE, merged_data$treat1, merged_data$treat2, 
                  merged_data$studlab, ref = "Placebo", sm='OR', random = F)
# CNMA network
netgraph(cnet_a)

### Identify pseudo paths by functions
contrib_B <- cnma_contrib (nc0, component = "B", 
                           inactive = "Placebo", random = FALSE)
paths_B <- contrib_B$Paths

### Estimating component B from each stream
##Stream1
S1_dt <- merged_data[2, ]
#log(OR) from Stream1
OR_B1 <- S1_dt$TE

##Stream2
S2_dt <- merged_data[c(1,4,8), ]
S2_Indirct <- discomb(S2_dt$TE,S2_dt$seTE,S2_dt$treat1,S2_dt$treat2,
                      S2_dt$studylab,sm= "OR",inactive = "Placebo")
#log(OR) from Stream2
OR_B2 <- S2_Indirct$Comp.common[2]

##Stream3
S3_dt <- merged_data[c(3,4), ]
S3_Indirct <- discomb(S3_dt$TE,S3_dt$seTE,S3_dt$treat1,S3_dt$treat2,
                      S3_dt$studylab,sm= "OR",inactive = "Placebo")
#log(OR) from Stream3
OR_B3 <- S3_Indirct$Comp.common[2]

##Stream4
S4_dt <- merged_data[c(4:6), ]
S4_Indirct <- discomb(S4_dt$TE,S4_dt$seTE,S4_dt$treat1,S4_dt$treat2,
                      S4_dt$studylab,sm= "OR",inactive = "Placebo")
#log(OR) from Stream4
OR_B4 <- S4_Indirct$Comp.common[2]

##Stream5
S5_dt <- merged_data[c(4,5,7), ]
S5_Indirct <- discomb(S5_dt$TE,S5_dt$seTE,S5_dt$treat1,S5_dt$treat2,
                      S5_dt$studylab,sm= "OR",inactive = "Placebo")
#log(OR) from Stream5
OR_B5 <- S5_Indirct$Comp.common[2]


# all pseudo paths' estimates weighted by flow equal to the CNMA estimate
paths_B$Estimate <- c(OR_B1, OR_B2, OR_B3, OR_B4, OR_B5)
round(sum(paths_B$Estimate * paths_B$Contrib_for_path),3)
round(OR_B,3)


###2.Data for Heart failure -----------
#HF_dt is the data in HF_dt.csv
HF_wide <- pairwise(treat = HF_dt$treat, event = HF_dt$acm_r, n = HF_dt$acm_n, 
        data = HF_dt, studlab = HF_dt$studyid, sm='RR', reference.group = 'PCO')
### Disconnected networks
HF_net0 <- netconnection(data = HF_wide,
              treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              TE = TE,
              seTE = seTE)

## Subnetwork1
HF_net1 <- netmeta(HF_wide$TE, HF_wide$seTE,  HF_wide$treat1, HF_wide$treat2, 
                   HF_wide$studlab, subset = HF_net0$subnet == 1)
netgraph(
  HF_net1,
  points = TRUE,
  number.of.studies = TRUE,
  col.points = "#4682b4",
  cex.points = 6   
)

## Subnetwork2
HF_net2 <- netmeta(HF_wide$TE, HF_wide$seTE,  HF_wide$treat1, HF_wide$treat2, 
                   HF_wide$studlab, subset = HF_net0$subnet == 2)
netgraph(
  HF_net2,
  points = TRUE,
  number.of.studies = TRUE,
  col.points = "#4682b4",
  cex.points = 6   
)

#### CNMA with interaction between MRA and RASi
HF_nc <- discomb(HF_wide$TE, HF_wide$seTE,  HF_wide$treat1, HF_wide$treat2, 
                  HF_wide$studlab, inactive = "PCO", sm='OR', random = F)
C <- HF_nc$C.matrix
C.int <- cbind(C, int_MS = C[, "MRA"] * C[, "RASi"])
HF_int <- discomb(C.matrix = C.int, HF_wide$TE, HF_wide$seTE,  HF_wide$treat1, HF_wide$treat2, 
                  HF_wide$studlab, inactive = "PCO", sm='OR', random = F)
HF_int
HF_int$comps
######The component estimate for MRA
RR_MRA <- HF_int$Comp.common[4]

# Xa matrix and comparisons in CNMA
reference <- "PCO"
X_matrix <- HF_int$X.matrix
Xa<- as.data.frame(X_matrix)
Xa <- Xa %>%
  mutate(!!reference := apply(., 1, function(row) {  
    if (all(row == 1 | row == 0)) {
      return(-1)
    } else if (all(row == -1 | row == 0)) {
      return(1)
    } else if (all(row %in% c(0, 1, -1))) {
      return(0)
    } else {
      return(NA)
    }
  }))

titles <- colnames(Xa)


Xa$treat1 <- apply(Xa, 1, function(row) {
  indices <- which(row == 1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$treat2 <- apply(Xa, 1, function(row) {
  indices <- which(row == -1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$comparison <- paste(Xa$treat1, Xa$treat2, sep = ":")

# Ha matrix at the component level
TE <- HF_int$TE
seTE <- HF_int$seTE.adj.common
weight <- 1/(seTE^2)
W <- diag(weight)
C <- HF_int$C.matrix
studylab <- HF_int$studlab
X_matrix <- HF_int$X.matrix
pseudo_inverse <- ginv(t(X_matrix) %*% W %*% X_matrix)
Ha <- pseudo_inverse %*% t(X_matrix) %*% W
rownames(Ha) <- HF_int$comps
colnames(Ha) <- Xa$comparison
 

## Construct CNMA network plot
merged_data <- data.frame(TE = HF_int$TE, seTE = HF_int$seTE.adj.common, 
                          treat1 = Xa$treat1, treat2 = Xa$treat2, 
                          studlab = HF_int$studlab, comparison = Xa$comparison)


cnet_HF <- netmeta(merged_data$TE, merged_data$seTE, merged_data$treat1, merged_data$treat2, 
                  merged_data$studlab, ref = "PCO", sm='RR', random = F)
# CNMA network
netgraph(
  cnet_HF,
  points = TRUE,
  number.of.studies = TRUE,
  col.points = "#4682b4",
  cex.points = 6   
)

### Identify pseudo paths by functions
contrib_MRA <- cnma_contrib (HF_int, component = "MRA", 
                           inactive = "PCO", random = FALSE)
path_MRA <- contrib_MRA$Paths


###Estimating component MRA from each pseudo path
##stream 1
S1_MRA <- merged_data[merged_data$comparison == "MRA:PCO", ]
#log(RR) from Stream1
RR_MRA1 <- S1_MRA$TE
#stream 2
S2_MRA <- merged_data[merged_data$comparison %in% 
                        c("MRA+int_MS:PCO", "PCO:RASi", "PCO:RASi+int_MS"),]
S2_MRA_Indirct <- discomb(S2_MRA$TE,S2_MRA$seTE,S2_MRA$treat1,S2_MRA$treat2,
                          S2_MRA$studylab,sm= "RR",inactive = "PCO")
RR_MRA2 <- S2_MRA_Indirct$Comp.common[2]
#stream 3
S3_MRA <- merged_data[merged_data$comparison %in% 
                        c("MRA+int_MS:PCO", "ARNI:RASi+int_MS", "ARNI:RASi"),]
S3_MRA_Indirct <- discomb(S3_MRA$TE,S3_MRA$seTE,S3_MRA$treat1,S3_MRA$treat2,
                          S3_MRA$studylab,sm= "RR",inactive = "PCO")
RR_MRA3 <- S3_MRA_Indirct$Comp.common[3]
#stream 4
S4_MRA <- merged_data[merged_data$comparison %in% 
                        c("MRA+int_MS:PCO", "BB:PCO", "BB:RASi","PCO:RASi+int_MS"),]
S4_MRA_Indirct <- discomb(S4_MRA$TE,S4_MRA$seTE,S4_MRA$treat1,S4_MRA$treat2,
                          S4_MRA$studylab,sm= "RR",inactive = "PCO")
RR_MRA4 <- S4_MRA_Indirct$Comp.common[3]

# all pseudo paths' estimates weighted by flow equal to the CNMA estimate
path_MRA$Estimate <- c(RR_MRA1, RR_MRA2, RR_MRA3, RR_MRA4)
round(sum(path_MRA$Estimate * path_MRA$Contrib_for_path),3)
round(RR_MRA,3)


###3.Data for Adverse event-----------------------------
#AEdata is the data from the study (doi: 10.1186/s12874-023-01959-9)
AEdata <- AEdata[!is.na(AEdata$lnRR), ]
## NMA model
net_AE <- netmeta(AEdata$lnRR, AEdata$selnRR, AEdata$t1, AEdata$t2, 
                  AEdata$id, ref = "plac", sm='RR', random = F)
# Construct NMA network plot
netgraph(
  net_AE,
  points = TRUE,
  number.of.studies = TRUE,
  col.points = "#4682b4",
  cex.points = 6,
  scale = 0.8
)

######## additive model
nc0_AE <- netcomb(net_AE, inactive = "plac")
######## The component estimate for trop 
RR_trop <- nc0_AE$Comp.common[15]

# Xa matrix and comparisons in CNMA
reference <- "plac"
X_matrix <- nc0_AE$X.matrix
Xa<- as.data.frame(X_matrix)
Xa <- Xa %>%
  mutate(!!reference := apply(., 1, function(row) {  # 动态设置列名
    if (all(row == 1 | row == 0)) {
      return(-1)
    } else if (all(row == -1 | row == 0)) {
      return(1)
    } else if (all(row %in% c(0, 1, -1))) {
      return(0)
    } else {
      return(NA)
    }
  }))

titles <- colnames(Xa)


Xa$treat1 <- apply(Xa, 1, function(row) {
  indices <- which(row == 1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$treat2 <- apply(Xa, 1, function(row) {
  indices <- which(row == -1)
  if (length(indices) > 0) {
    return(paste(titles[indices], collapse = "+"))
  } else {
    return(NA)
  }
})
Xa$comparison <- paste(Xa$treat1, Xa$treat2, sep = ":")

# Ha matrix at the component level
TE <- nc0_AE$TE
seTE <- nc0_AE$seTE.adj.common
weight <- 1/(seTE^2)
W <- diag(weight)
C <- nc0_AE$C.matrix
studylab <- nc0_AE$studlab
X_matrix <- nc0_AE$X.matrix
pseudo_inverse <- ginv(t(X_matrix) %*% W %*% X_matrix)
Ha <- pseudo_inverse %*% t(X_matrix) %*% W
rownames(Ha) <- nc0_AE$comps
colnames(Ha) <- nc0_AE$comparison


## Construct CNMA network plot
merged_data <- data.frame(TE = nc0_AE$TE, seTE = nc0_AE$seTE.adj.common, 
                          treat1 = Xa$treat1, treat2 = Xa$treat2, 
                          comparison = Xa$comparison, studlab= nc0_AE$studlab)
merged_data$studlab <- c(1:89)
AE_cnet <- netmeta(merged_data$TE, merged_data$seTE, merged_data$treat1, merged_data$treat2, 
                  merged_data$studlab, sm='RR', random = F)
# CNMA network
netgraph(AE_cnet,
         points = TRUE,
         number.of.studies = TRUE,
         col.points = "#4682b4",
         cex.points = 6,
         scale = 0.8)
### Identify pseudo paths by functions
contrib_trop <- cnma_contrib (nc0_AE, component = "trop", 
                             inactive = "plac", random = FALSE)
path_trop <- contrib_trop$Paths
path_trop$Flow_for_path <- round(path_trop$Flow_for_path, 3)
path_trop <- path_trop[path_trop$Flow_for_path != 0, ]

##Estimating component MRA from each network-flow stream
#stream 1
S1_trop <- merged_data[merged_data$comparison == "plac:trop", ]
S1_trop_dirct <- discomb(S1_trop$TE,S1_trop$seTE,S1_trop$treat1,S1_trop$treat2,
                         S1_trop$studylab,sm= "RR",inactive = "plac")
S1_trop_dirct$Comp.common[1]
#stream 2
S2_trop <- merged_data[merged_data$comparison %in% c("onda:plac", "onda:trop"),]
S2_trop_indirct <- discomb(S2_trop$TE,S2_trop$seTE,S2_trop$treat1,S2_trop$treat2,
                           S2_trop$studylab,sm= "RR",inactive = "plac")
S2_trop_indirct$Comp.common[2]
#stream 3
S3_trop <- merged_data[merged_data$comparison %in% c("dexa:plac", "dexa+trop:plac"),]
S3_trop_indirct <- discomb(S3_trop$TE,S3_trop$seTE,S3_trop$treat1,S3_trop$treat2,
                           S3_trop$studylab,sm= "RR",inactive = "plac")
S3_trop_indirct$Comp.common[2]
#stream 4
S4_trop <- merged_data[merged_data$comparison %in% c("dexa:trop", "dexa+trop:plac"),]
S4_trop_indirct <- discomb(S4_trop$TE,S4_trop$seTE,S4_trop$treat1,S4_trop$treat2,
                           S4_trop$studylab,sm= "RR",inactive = "plac")
S4_trop_indirct$Comp.common[2]
#stream 5
S5_trop <- merged_data[merged_data$comparison %in% c("plac:onda","onda:trop"),]
S5_trop_indirct <- discomb(S5_trop$TE,S5_trop$seTE,S5_trop$treat1,S5_trop$treat2,
                           S5_trop$studylab,sm= "RR",inactive = "plac")
S5_trop_indirct$Comp.common[2]
#stream 6
S6_trop <- merged_data[merged_data$comparison %in% c("caso:plac",
                                                     "onda:trop", "caso:onda"),]
S6_trop_indirct <- discomb(S6_trop$TE,S6_trop$seTE,S6_trop$treat1,S6_trop$treat2,
                           S6_trop$studylab,sm= "RR",inactive = "plac")
S6_trop_indirct$Comp.common[3]

#stream 7
S7_trop <- merged_data[merged_data$comparison %in% c("onda:ramo", "plac:ramo",
                                                     "onda:trop"),]
S7_trop_indirct <- discomb(S7_trop$TE,S7_trop$seTE,S7_trop$treat1,S7_trop$treat2,
                           S7_trop$studylab,sm= "RR",inactive = "plac")
S7_trop_indirct$Comp.common[3]

#stream 8
S8_trop <- merged_data[merged_data$comparison %in% c("gran:plac", "gran:onda",
                                                     "onda:trop"),]
S8_trop_indirct <- discomb(S8_trop$TE,S8_trop$seTE,S8_trop$treat1,S8_trop$treat2,
                           S8_trop$studylab,sm= "RR",inactive = "plac")
S8_trop_indirct$Comp.common[3]

#stream 9
S9_trop <- merged_data[merged_data$comparison %in% c("palo:plac", "onda:palo",
                                                     "onda:trop"),]
S9_trop_indirct <- discomb(S9_trop$TE, S9_trop$seTE, S9_trop$treat1, S9_trop$treat2,
                           S9_trop$studylab, sm = "RR", inactive = "plac")
S9_trop_indirct$Comp.common[3]

#stream 10
S10_trop <- merged_data[merged_data$comparison %in% c("drop:plac", "drop:onda",
                                                      "onda:trop"),]
S10_trop_indirct <- discomb(S10_trop$TE, S10_trop$seTE, S10_trop$treat1, S10_trop$treat2,
                            S10_trop$studylab, sm = "RR", inactive = "plac")
S10_trop_indirct$Comp.common[3]

#stream 11
S11_trop <- merged_data[merged_data$comparison %in% c("meto:onda", "meto:plac",
                                                      "onda:trop"),]
S11_trop_indirct <- discomb(S11_trop$TE, S11_trop$seTE, S11_trop$treat1, S11_trop$treat2,
                            S11_trop$studylab, sm = "RR", inactive = "plac")
S11_trop_indirct$Comp.common[3]

#stream 12
S12_trop <- merged_data[merged_data$comparison %in% c("onda:trop", "dexa+trop:plac",
                                                      "dexa:drop","drop:onda"),]

S12_trop_indirct <- discomb(S12_trop$TE, S12_trop$seTE, S12_trop$treat1, S12_trop$treat2,
                            S12_trop$studylab, sm = "RR", inactive = "plac")

S12_trop_indirct$Comp.common[4]

path_estimates <- c(
  S1_trop_dirct$Comp.common[1],
  S2_trop_indirct$Comp.common[2],
  S3_trop_indirct$Comp.common[2],
  S4_trop_indirct$Comp.common[2],
  S5_trop_indirct$Comp.common[2],
  S6_trop_indirct$Comp.common[3],
  S7_trop_indirct$Comp.common[3],
  S8_trop_indirct$Comp.common[3],
  S9_trop_indirct$Comp.common[3],
  S10_trop_indirct$Comp.common[3],
  S11_trop_indirct$Comp.common[3],
  S12_trop_indirct$Comp.common[4]
)
path_trop$Estimate <- path_estimates
# all pseudo paths' estimates weighted by flow equal to the CNMA estimate
round(sum(path_trop$Estimate * path_trop$Contrib_for_path),3)
round(RR_trop,3)

