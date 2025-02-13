#### Data preparation ####

# Create a logical OR column
OR <- matchRTComp$radiation_treatment_adjuvant | matchRTComp$targeted_molecular_therapy

# Check if OR_AB matches C
match <- OR == matchRTComp$newTRT

# Summarize results
if (all(match)) {
  print("Variable C is the logical OR of A and B for all observations")
} else {
  print("Variable C is NOT the logical OR of A and B for some observations")
}

rm(match, OR)
# newTRT = radiation_treatment_adjuvant OR targeted_molecular_therapy !!!

sum(matchRTComp$newTRT)
sum(matchRTComp$radiation_treatment_adjuvant)
sum(matchRTComp$targeted_molecular_therapy)
# 79 people got newTRT (at least 1 of the 2 treatments)

matchRTComp[matchRTComp[,c("as.factor(myRsp)")]==1,c("treatment_outcome_first_course","as.factor(myRsp)")]
# treatment_outcome_first_course IS THE SAME AS as.factor(myRsp)
# Progressive Disease = 0
# Stable Disease AND Partial Remission/Response = 1
# Complete Remission/Response = 2


matchRTComp$treatment_outcome_first_course=factor(
  matchRTComp$treatment_outcome_first_course,
  levels = c("Progressive Disease","Stable Disease", "Partial Remission/Response", 
             "Complete Remission/Response")
)

# Plot heatmap
table_df=table(matchRTComp$newTRT, matchRTComp$`as.factor(myRsp)`, matchRTComp$patients.tumor_grade)
table_df

table_df2=table(matchRTComp$newTRT, matchRTComp$treatment_outcome_first_course, matchRTComp$patients.tumor_grade)
table_df2

table_df3=table(matchRTComp$newTRT, matchRTComp$`as.factor(myRsp)`)
table_df3

# What is variable "patients.tumor_grade"?

#library(BART)
#data=matchRTComp
#data=data[,-c(1,2,3,4,9)]
#y=matchRTComp[,"as.factor(myRsp)"]

#WE'RE GOING TO TUNE: BASE & POWER
#WE JUST FIND VALUES FOR ALPHA, BETA AND TAU

#USE THE SAME DATA FOR PREDICTION TO LOOK AT YHAT CHAINS
#WHEN CHAINS CONVERGE PASS AS TEST SET THE SAME THING WITH NEWTREAT SWAPPED
#USE OUTPUTS TO COMPUTE DELTAJ (treat-eff-4-ord-outc pg.22)

#### useless - plot inv-gamma ####

library(extraDistr)

# Values for x
x <- seq(0.01, 3, length.out = 250)

# Compute density
a=7
b=2
y <- dinvgamma(x, alpha = a, beta = b)
yy <- pinvgamma(x, alpha = a, beta = b)
# Plot the density
plot(x, y, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Density",
     main = paste("Inverse Gamma Density (alpha =", a, ", beta =", b, ")"),
     ylim=c(0,2.5))
grid()
lines(x, yy, type = "l", col = "red", lwd = 2)
print(b/(a-1))


#### RUN THE MCMC ####

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
# y has 1 2 3 values

#var(y)

#x=mbart(data, y,data,
#        power=2.5, base=0.95,
#        ndpost=5000, nskip=1000, printevery=5000,
#        ntree=50, k=2.0, a=0.5, b=2, sparse=FALSE)

x=mbart(data, y, data, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
        base=0.95, power=2, sparse = FALSE)
#x=mbart(data, y, data, ndpost=5000, nskip=1000, printevery=5000, ntree=50,
#        base=0.95, power=3, sparse = TRUE)

#sensitivity analysis in report

#names(x)

prob_test_mean=x$prob.test.mean
prob_test_mean=matrix(prob_test_mean, nrow=158, byrow=TRUE)

pred_class <- apply(prob_test_mean, 1, which.max)

sum(y==pred_class)
table(pred_class, y)

table(y)
table(pred_class)

#### CI on failed obs (2 classes) ####

# Identify misclassified observations
misclassified_idx <- which(y != pred_class)

# Extract the posterior probabilities for the misclassified observations
# x$prob.test contains the probabilities for each class (rows correspond to MCMC iterations)
prob_test <- x$prob.test

# Reorganize the matrix to align probabilities with patients and classes
# Rows: Iterations, Columns: Classes (in blocks of 158 for each class)
#prob_test_matrix <- array(prob_test, dim = c(5000, 158, 3)) # Shape: [Iterations, Patients, Classes]

# For misclassified patients, extract the probabilities of the correct and predicted classes
credible_intervals <- lapply(misclassified_idx, function(idx) {
  correct_class <- y[idx]
  predicted_class <- pred_class[idx]
  
  # Probabilities for the correct and predicted classes
  prob_correct <- prob_test[, 3*(idx-1)+correct_class]
  prob_predicted <- prob_test[, 3*(idx-1)+predicted_class]
  
  # Compute 2.5% and 97.5% quantiles for both classes
  correct_ci <- quantile(prob_correct, probs = c(0.025, 0.975))
  predicted_ci <- quantile(prob_predicted, probs = c(0.025, 0.975))
  
  list(
    patient_id = idx,
    correct_class = correct_class,
    predicted_class = predicted_class,
    correct_ci = correct_ci,
    predicted_ci = predicted_ci
  )
})

# Print credible intervals for each misclassified patient
# for (interval in credible_intervals) {
#   cat(sprintf(
#     "Patient %d:\n  Correct Class: %d, Predicted Class: %d\n  CI (Correct): [%.3f, %.3f]\n  CI (Predicted): [%.3f, %.3f]\n\n",
#     interval$patient_id,
#     interval$correct_class,
#     interval$predicted_class,
#     interval$correct_ci[1], interval$correct_ci[2],
#     interval$predicted_ci[1], interval$predicted_ci[2]
#   ))
# }

# Load required library
library(ggplot2)

# Prepare data for plotting
plot_data <- do.call(rbind, lapply(credible_intervals, function(interval) {
  data.frame(
    patient_id = interval$patient_id,
    class = c("Correct", "Predicted"),
    lower = c(interval$correct_ci[1], interval$predicted_ci[1]),
    upper = c(interval$correct_ci[2], interval$predicted_ci[2]),
    mean_value = c(mean(prob_test[, 3*(interval$patient_id-1) + interval$correct_class]),
                   mean(prob_test[, 3*(interval$patient_id-1) + interval$predicted_class]))
  )
}))
# Plot credible intervals with mean values
ggplot(plot_data, aes(x = factor(patient_id), y = mean_value, color = class)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Credible Intervals for Correct and Predicted Classes",
    x = "Patient ID",
    y = "Probability",
    color = "Class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("Correct" = "blue", "Predicted" = "red"))

#### CI on failed obs (3 classes) ####


# Prepare data for plotting (all 3 classes for misclassified patients)
# Prepare data for plotting (all 3 classes with custom labels)
plot_data_custom_labels <- do.call(rbind, lapply(misclassified_idx, function(idx) {
  correct_class <- y[idx]
  predicted_class <- pred_class[idx]
  other_class <- setdiff(1:3, c(correct_class, predicted_class))
  
  # Custom labels for the classes
  labels <- c(
    "Correct" = correct_class,
    "Predicted" = predicted_class,
    "Other" = other_class
  )
  
  data.frame(
    patient_id = idx,
    class = factor(c("Correct", "Predicted", "Third"), levels = c("Correct", "Predicted", "Third")),
    lower = sapply(labels, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.025)),
    upper = sapply(labels, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.975)),
    mean_value = sapply(labels, function(cls) mean(prob_test[,3*(idx-1)+cls]))
  )
}))

# Plot credible intervals with mean values for all 3 custom-labeled classes
ggplot(plot_data_custom_labels, aes(x = factor(patient_id), y = mean_value, color = class)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Credible Intervals for Misclassified Observations (Custom Labels)",
    x = "Patient ID",
    y = "Probability",
    color = "Class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  
  scale_color_manual(values = c("Correct" = "blue", "Predicted" = "red", "Third" = "green"))

#### CI on good obs (3 classes) ####

# Identify classified observations
classified_idx <- which(y == pred_class)

halfway <- ceiling(length(classified_idx) / 2)
classified_idx_1 <- classified_idx[1:halfway]
classified_idx_2 <- classified_idx[(halfway + 1):length(classified_idx)]

# Prepare data for plotting (all 3 classes for classified patients)
# Prepare data for plotting (all 3 classes with custom labels)
plot_data_all_classes_1 <- do.call(rbind, lapply(classified_idx_1, function(idx) {
  data.frame(
    patient_id = idx,
    class = factor(1:3, labels = c("Class 1", "Class 2", "Class 3")),
    lower = sapply(1:3, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.025)),
    upper = sapply(1:3, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.975)),
    mean_value = sapply(1:3, function(cls) mean(prob_test[,3*(idx-1)+cls]))
  )
}))


# Plot credible intervals with mean values for all 3 custom-labeled classes
ggplot(plot_data_all_classes_1, aes(x = factor(patient_id), y = mean_value, color = class)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Credible Intervals for Correcly Classified Observations",
    x = "Patient ID",
    y = "Probability",
    color = "Class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  
  scale_color_manual(values = c("Class 1" = "blue", "Class 2" = "red", "Class 3" = "green"))

plot_data_all_classes_2 <- do.call(rbind, lapply(classified_idx_2, function(idx) {
  data.frame(
    patient_id = idx,
    class = factor(1:3, labels = c("Class 1", "Class 2", "Class 3")),
    lower = sapply(1:3, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.025)),
    upper = sapply(1:3, function(cls) quantile(prob_test[, 3*(idx-1)+cls], probs = 0.975)),
    mean_value = sapply(1:3, function(cls) mean(prob_test[,3*(idx-1)+cls]))
  )
}))


# Plot credible intervals with mean values for all 3 custom-labeled classes
ggplot(plot_data_all_classes_2, aes(x = factor(patient_id), y = mean_value, color = class)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  labs(
    title = "Credible Intervals for Correcly Classified Observations",
    x = "Patient ID",
    y = "Probability",
    color = "Class"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  
  scale_color_manual(values = c("Class 1" = "blue", "Class 2" = "red", "Class 3" = "green"))


#### ex of mcmc ####

#EX_1, obs 1 belongs to class 2

# Define observation index and class
observation <- 1  # Observation index (1 to 158)
class <- 2        # Class index (1, 2, or 3)
# Compute the starting column index for the observation
start_col <- (observation - 1) * 3 + (class)
prob_test=x$prob.test
trace=prob_test[, start_col]
# Plot the trace plot
plot(trace, type = "l", col = "green",
     xlab = "MCMC Iteration", ylab = "Predicted Probability",
     main = paste("Trace Plot for Observation", observation, "Class", class))

#EX_2, obs 158 belongs to class 2

# Define observation index and class
observation <- 158  # Observation index (1 to 158)
class <- 3      # Class index (1, 2, or 3)
# Compute the starting column index for the observation
start_col <- (observation - 1) * 3 + (class)
#prob_test=x$prob.test
trace2=prob_test[, start_col]
# Plot the trace plot
plot(trace2, type = "l", col = "green",
     xlab = "MCMC Iteration", ylab = "Predicted Probability",
     main = paste("Trace Plot for Observation", observation, "Class", class))

#### Var importance ####

var_imp1=x$varcount.mean[1,]
var_imp2=x$varcount.mean[2,]

var_imp=(var_imp1+var_imp2)/2

#var_importance = (var_importance1+var_importance2)/2
#names(var_importance)=colnames(x$varcount[[1]])
  
#print(var_importance)

par(mfrow=c(2,1))
barplot(var_imp1, main="Variable importance from BART", las=2)
barplot(var_imp2, main="Variable importance from BART", las=2)

par(mfrow=c(1,1))
barplot(var_imp, main="Variable importance from BART", las=2)

par(mfrow=c(1,1))
barplot(sort(var_imp, TRUE), main="Variable importance from BART", las=2)

plot(sort(var_imp, TRUE))
abline(h=2.65)
abline(v=5.5)
abline(h=3.9)
abline(v=21.5)

#### Computation of tau_i and eta_i ####

#for each i and for each iteration we compute pr{Yi(1)>=Yi(0)}
#we obtain a matrix 5000x158 of tau_i
#we compute 158 credible intervals of tau_i
#we can also compute eta, it's the same machinery, just pr{Yi(1)>Yi(0)}

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)

data_test1=data
data_test1$newTRT=1
data_test0=data
data_test0$newTRT=0

data_new=rbind(data_test0, data_test1)

x2=mbart(data, y, data_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
         base=0.95, power=2, sparse = FALSE)

#x2=mbart(data, y, data_new, ndpost=5000, nskip=1000, printevery=5000, ntree=50,
#        base=0.95, power=3, sparse = TRUE)

str(x2$prob.test)

mat_for_Y_0=x2$prob.test[,1:474]
mat_for_Y_1=x2$prob.test[,475:948]

tau_i=matrix(0,5000,158)

probabilities_tau<-function(x1,x0){
  #compute probability that X>=Y according to the distributions: x,y
  p=x0[1]*(x1[1]+x1[2]+x1[3])+x0[2]*(x1[2]+x1[3])+x0[3]*x1[3]
  return(p)
}

probabilities_eta<-function(x1,x0){
  #compute probability that X>Y according to the distributions: x,y
  p=x0[1]*(x1[2]+x1[3])+x0[2]*(x1[3])
  return(p)
}

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

credible_intervals_eta[1,]
credible_intervals_tau[1,]

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
abline(h=0.09)

plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])
abline(h=0.8)


#### Computation of deltaj_i ####

delta_j_i=matrix(0,5000,474)

compute_delta=function(x1,x0,j){
  #compute probability that X>=j - probability that Y>=j
  #according to the distributions: x,y
  p1=sum(x1[j:3])
  p0=sum(x0[j:3])
  return(p1-p0)
}

for(ii in 1:5000){
  for(j in 1:3){
    for(i in 1:158){
      delta_j_i[ii,3*(i-1)+j]=compute_delta(mat_for_Y_1[ii,(3*(i-1)+1):(3*(i-1)+3)],mat_for_Y_0[ii,(3*(i-1)+1):(3*(i-1)+3)],j)
    }
  }
}
ind0=seq(from = 1, to = 472, by = 3)
ind1=seq(from = 2, to = 473, by = 3)
ind2=seq(from = 3, to = 474, by = 3)

delta_0=delta_j_i[,ind0]
delta_1=delta_j_i[,ind1]
delta_2=delta_j_i[,ind2]

credible_intervals_delta_0 <- apply(delta_0, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_0 <- t(credible_intervals_delta_0)
credible_intervals_delta_0

credible_intervals_delta_1 <- apply(delta_1, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_1 <- t(credible_intervals_delta_1)
credible_intervals_delta_1

credible_intervals_delta_2 <- apply(delta_2, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_2 <- t(credible_intervals_delta_2)
credible_intervals_delta_2


plot(credible_intervals_delta_1[,1], ylim=c(-0.5,0.5))
points(credible_intervals_delta_1[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_1[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_1[, 2], col = "gray")
abline(h=-0.015)

plot(credible_intervals_delta_2[,1], ylim=c(-0.5,0.5))
points(credible_intervals_delta_2[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_2[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_2[, 2], col = "gray")
abline(h=-0.1)

#### Computing general causal effects ####

str(eta_i)
mean_eta=rowMeans(eta_i)
quantiles_general_eta=quantile(eta_i, probs=c(0.025, 0.975))
quantiles_general_eta

mean_tau=rowMeans(tau_i)
quantiles_general_tau=quantile(tau_i, probs=c(0.025, 0.975))
quantiles_general_tau

mean_delta_1=rowMeans(delta_1)
quantiles_general_delta_1=quantile(delta_1, probs=c(0.025, 0.975))
quantiles_general_delta_1

mean_delta_2=rowMeans(delta_2)
quantiles_general_delta_2=quantile(delta_2, probs=c(0.025, 0.975))
quantiles_general_delta_2

#### Propensity score ####

data_for_pscore = data[,-5]
y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
pscore = glm(data$newTRT ~ ., data=data_for_pscore, family = binomial)$fitted.values
colors <- ifelse(data$newTRT == 1, "blue", "red")
# First plot
plot(pscore, col = colors, pch = 16, 
     main = "Pscore by newTRT", 
     xlab = "Index", 
     ylab = "Pscore")
legend("topright", legend = c("newTRT = 1", "newTRT = 0"), 
       col = c("blue", "red"), pch = 16, 
       title = "Group")
# Second plot (sorted pscore)
plot(sort(pscore), col = colors[order(pscore)], pch = 16, 
     main = "Sorted Pscore by newTRT", 
     xlab = "Index", 
     ylab = "Sorted Pscore")
legend("topleft", legend = c("newTRT = 1", "newTRT = 0"), 
       col = c("blue", "red"), pch = 16, 
       title = "Group")

Neyman_SRE=function(z,y,x){
  xlevels = unique(x)
  K = length(xlevels)
  PiK = rep(0,K)
  TauK = rep(0,K)
  varK = rep(0,K)
  for(k in 1:K){
    xk = xlevels[k]
    zk = z[x==xk]
    yk = y[x==xk]
    PiK = length(zk)/length(z)
    TauK = mean(yk[zk==1])-mean(yk[zk==0])
    varK = var(yk[zk==1])/sum(zk) + var(yk[zk==0])/sum(1-zk)
  }
  return(c(sum(PiK*TauK), sum(PiK^2*varK)))
}

n.strata = c(3, 4, 5, 10, 20)
z=data[,5]
strat.res = sapply(n.strata, FUN = function(nn){
  q.pscore = quantile(pscore, (1:(nn-1))/nn)
  ps.strata = cut(pscore, breaks = c(0,q.pscore,1),
                  labels = 1:nn) 
  Neyman_SRE(z, y, ps.strata)})

rownames(strat.res) = c("est", "se") 
colnames(strat.res) = n.strata

round(strat.res, 3)

#Increasing K from 3 to 20 reduces the standard error. 
#However, we cannot go as extreme as K = 25 because 
#the standard error is not well-defined in some strata 
#with only a few control units. 
#The above estimators show a negative but insignificant 
#effect of the treatment program on y.

# Here we try to run mbart and group tau_i from the same group
# together, and see if we get something out of it

n_groups=3 #4 did now work

quantiles <- quantile(pscore, probs = seq(0, 1, length.out = n_groups + 1))

# Plot the data
plot(sort(pscore), col = colors[order(pscore)], pch = 16, 
     main = paste("Data with", n_groups, "Quantile Groups"), 
     xlab = "Index", 
     ylab = "Sorted Pscore")
legend("topleft", legend = c("newTRT = 1", "newTRT = 0"), 
       col = c("blue", "red"), pch = 16, 
       title = "Group")

# Add quantile lines
abline(h = quantiles[-c(1, length(quantiles))], col = "red", lty = 2)

groups <- cut(pscore, breaks = quantiles, include.lowest = TRUE, labels = FALSE)

# Find the indexes for each group
gr_1=which(groups==1)
gr_2=which(groups==2)
gr_3=which(groups==3)
gr_4=which(groups==4)

data_pscore=data
data_pscore[,"pscore"]=pscore

x3=mbart(data_pscore, y, data_pscore, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
        base=0.95, power=2, sparse = FALSE)

prob_test_mean=x3$prob.test.mean
prob_test_mean=matrix(prob_test_mean, nrow=158, byrow=TRUE)

pred_class <- apply(prob_test_mean, 1, which.max)

sum(y==pred_class)
table(pred_class, y)

table(y)
table(pred_class)

data_test1=data_pscore
data_test1$newTRT=1
data_test0=data_pscore
data_test0$newTRT=0

data_pscore_new=rbind(data_test0, data_test1)

x4=mbart(data_pscore, y, data_pscore_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
         base=0.95, power=2, sparse = FALSE)

mat_for_Y_0=x4$prob.test[,1:474]
mat_for_Y_1=x4$prob.test[,475:948]

tau_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
abline(h=0.2)

plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])
abline(h=0.7)

# Now average the rows wrt our groups

# Create a list of index vectors for simplicity
index_list <- list(gr_1, gr_2, gr_3)#, gr_4)

# Initialize the result matrix
eta_i_groups<- matrix(NA, nrow = nrow(eta_i), ncol = length(index_list))
tau_i_groups<- matrix(NA, nrow = nrow(tau_i), ncol = length(index_list))

# Compute the averages
for (i in seq_along(index_list)) {
  eta_i_groups[, i] <- rowMeans(eta_i[, index_list[[i]], drop = FALSE])
  tau_i_groups[, i] <- rowMeans(tau_i[, index_list[[i]], drop = FALSE])
}

credible_intervals_eta_groups <- apply(eta_i_groups, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta_groups <- t(credible_intervals_eta_groups)

credible_intervals_tau_groups <- apply(tau_i_groups, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau_groups <- t(credible_intervals_tau_groups)

plot(credible_intervals_eta_groups[,1], ylim=c(0,1))
points(credible_intervals_eta_groups[,2])
abline(h=0.2)

plot(credible_intervals_tau_groups[,1], ylim=c(0,1))
points(credible_intervals_tau_groups[,2])
abline(h=0.7)

# Propensity Score v2 #####
# Here we try to run mbart and group tau_i from the same group
# together, and see if we get something out of it

n_groups=3 #4 did now work

quantiles <- quantile(pscore, probs = seq(0, 1, length.out = n_groups + 1))

# Plot the data
plot(sort(pscore), col = colors[order(pscore)], pch = 16, 
     main = paste("Data with", n_groups, "Quantile Groups"), 
     xlab = "Index", 
     ylab = "Sorted Pscore")
legend("topleft", legend = c("newTRT = 1", "newTRT = 0"), 
       col = c("blue", "red"), pch = 16, 
       title = "Group")

# Add quantile lines
abline(h = quantiles[-c(1, length(quantiles))], col = "red", lty = 2)

groups <- cut(pscore, breaks = quantiles, include.lowest = TRUE, labels = FALSE)

# Find the indexes for each group
gr_1=which(groups==1)
gr_2=which(groups==2)
gr_3=which(groups==3)
gr_4=which(groups==4)

data_pscore=data
data_pscore[,"pscore"]=pscore

x3=mbart(data_pscore, y, data_pscore, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
         base=0.95, power=2, sparse = FALSE)

prob_test_mean=x3$prob.test.mean
prob_test_mean=matrix(prob_test_mean, nrow=158, byrow=TRUE)

pred_class <- apply(prob_test_mean, 1, which.max)

sum(y==pred_class)
table(pred_class, y)

table(y)
table(pred_class)

data_test1=data_pscore
data_test1$newTRT=1
data_test0=data_pscore
data_test0$newTRT=0

data_pscore_new=rbind(data_test0, data_test1)

x4=mbart(data_pscore, y, data_pscore_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
         base=0.95, power=2, sparse = FALSE)

mat_for_Y_0=x4$prob.test[,1:474]
mat_for_Y_1=x4$prob.test[,475:948]

tau_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
abline(h=0.2)

plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])
abline(h=0.7)

# Now average the rows wrt our groups

# Create a list of index vectors for simplicity
index_list <- list(gr_1, gr_2, gr_3)#, gr_4)

# Initialize the result matrix
eta_i_groups<- matrix(NA, nrow = nrow(eta_i), ncol = length(index_list))
tau_i_groups<- matrix(NA, nrow = nrow(tau_i), ncol = length(index_list))

# Compute the averages
for (i in seq_along(index_list)) {
  eta_i_groups[, i] <- rowMeans(eta_i[, index_list[[i]], drop = FALSE])
  tau_i_groups[, i] <- rowMeans(tau_i[, index_list[[i]], drop = FALSE])
}

credible_intervals_eta_groups <- apply(eta_i_groups, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta_groups <- t(credible_intervals_eta_groups)

credible_intervals_tau_groups <- apply(tau_i_groups, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau_groups <- t(credible_intervals_tau_groups)

plot(credible_intervals_eta_groups[,1], ylim=c(0,1))
points(credible_intervals_eta_groups[,2])
abline(h=0.2)

plot(credible_intervals_tau_groups[,1], ylim=c(0,1))
points(credible_intervals_tau_groups[,2])
abline(h=0.7)


#### simulated dataset #####

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
data_sim=data
for (i in 1:158){
  if (y[i]==3)
    data_sim$newTRT[i]=1
  else
    data_sim$newTRT[i]=0
}

data_test_sim1=data_sim
data_test_sim1$newTRT=1
data_test_sim0=data_sim
data_test_sim0$newTRT=0

data_new=rbind(data_test_sim0, data_test_sim1)

x_sim=mbart(data_sim, y, data_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
            base=0.95, power=2, sparse = FALSE)

mat_for_Y_0=x_sim$prob.test[,1:474]
mat_for_Y_1=x_sim$prob.test[,475:948]

delta_j_i=matrix(0,5000,474)

for(ii in 1:5000){
  for(j in 1:3){
    for(i in 1:158){
      delta_j_i[ii,3*(i-1)+j]=compute_delta(mat_for_Y_1[ii,(3*(i-1)+1):(3*(i-1)+3)],mat_for_Y_0[ii,(3*(i-1)+1):(3*(i-1)+3)],j)
    }
  }
}
ind0=seq(from = 1, to = 472, by = 3)
ind1=seq(from = 2, to = 473, by = 3)
ind2=seq(from = 3, to = 474, by = 3)

delta_0=delta_j_i[,ind0]
delta_1=delta_j_i[,ind1]
delta_2=delta_j_i[,ind2]

credible_intervals_delta_0 <- apply(delta_0, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_0 <- t(credible_intervals_delta_0)
credible_intervals_delta_0

credible_intervals_delta_1 <- apply(delta_1, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_1 <- t(credible_intervals_delta_1)
credible_intervals_delta_1

credible_intervals_delta_2 <- apply(delta_2, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_2 <- t(credible_intervals_delta_2)
credible_intervals_delta_2


plot(credible_intervals_delta_1[,1], ylim=c(0,1))
points(credible_intervals_delta_1[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_1[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_1[, 2], col = "gray")

plot(credible_intervals_delta_2[,1], ylim=c(0,1))
points(credible_intervals_delta_2[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_2[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_2[, 2], col = "gray")

tau_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

credible_intervals_eta[1,]
credible_intervals_tau[1,]

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])


plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])


####simulated dataset v2 ####

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
data_sim=data
for (i in 1:158){
  if (data_sim$patients.tumor_grade[i]==1){
    if (y[i]>1)
      data_sim$newTRT[i]=1
    else
      data_sim$newTRT[i]=0
  }
}

data_test_sim1=data_sim
data_test_sim1$newTRT=1
data_test_sim0=data_sim
data_test_sim0$newTRT=0

data_new=rbind(data_test_sim0, data_test_sim1)

x_sim=mbart(data_sim, y, data_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
            base=0.95, power=2, sparse = FALSE)

mat_for_Y_0=x_sim$prob.test[,1:474]
mat_for_Y_1=x_sim$prob.test[,475:948]

delta_j_i=matrix(0,5000,474)

for(ii in 1:5000){
  for(j in 1:3){
    for(i in 1:158){
      delta_j_i[ii,3*(i-1)+j]=compute_delta(mat_for_Y_1[ii,(3*(i-1)+1):(3*(i-1)+3)],mat_for_Y_0[ii,(3*(i-1)+1):(3*(i-1)+3)],j)
    }
  }
}
ind0=seq(from = 1, to = 472, by = 3)
ind1=seq(from = 2, to = 473, by = 3)
ind2=seq(from = 3, to = 474, by = 3)

delta_0=delta_j_i[,ind0]
delta_1=delta_j_i[,ind1]
delta_2=delta_j_i[,ind2]

credible_intervals_delta_0 <- apply(delta_0, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_0 <- t(credible_intervals_delta_0)
credible_intervals_delta_0

credible_intervals_delta_1 <- apply(delta_1, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_1 <- t(credible_intervals_delta_1)
credible_intervals_delta_1

credible_intervals_delta_2 <- apply(delta_2, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_2 <- t(credible_intervals_delta_2)
credible_intervals_delta_2

par(mfrow=c(1,2))
plot(credible_intervals_delta_1[which(data_sim$patients.tumor_grade==1),1], ylim=c(-0.5,0.5),col='red')
points(credible_intervals_delta_1[which(data_sim$patients.tumor_grade==1),2], col='red')
plot(credible_intervals_delta_1[which(data_sim$patients.tumor_grade==0),1],ylim=c(-0.5,0.5))
points(credible_intervals_delta_1[which(data_sim$patients.tumor_grade==0),2])
segments(x0 = 1:158, y0 = credible_intervals_delta_1[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_1[, 2], col = "gray")
t.test(credible_intervals_delta_1[which(data_sim$patients.tumor_grade==1),2],credible_intervals_delta_1[which(data_sim$patients.tumor_grade==0),2], alternative='greater')
plot(credible_intervals_delta_2[which(data_sim$patients.tumor_grade==1),1], ylim=c(-0.5,0.5),col='red')
points(credible_intervals_delta_2[which(data_sim$patients.tumor_grade==1),2], col='red')
plot(credible_intervals_delta_2[which(data_sim$patients.tumor_grade==0),1],ylim=c(-0.5,0.5))
points(credible_intervals_delta_2[which(data_sim$patients.tumor_grade==0),2])
t.test(credible_intervals_delta_2[which(data_sim$patients.tumor_grade==1),2],credible_intervals_delta_2[which(data_sim$patients.tumor_grade==0),2], alternative='greater')

segments(x0 = 1:158, y0 = credible_intervals_delta_2[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_2[, 2], col = "gray")

tau_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

credible_intervals_eta[1,]
credible_intervals_tau[1,]

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
abline(h=0.09)

plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])


library(bayesreg)

#standardize
data_sim[,-c(1,2,5)] <- scale(data_sim[,-c(1,2,5)])

y_bin1 = as.factor(ifelse(y==1,0,1))
fit1 = bayesreg(y_bin1~.,data_sim,n.samples=1500,model="logistic",prior="horseshoe")

y2 = y[y_bin1==1]
y_bin2 = as.factor(ifelse(y2==2,0,1))
fit2 = bayesreg(y_bin2~.,data_sim[y_bin1==1,],n.samples=1500,model="logistic",prior="horseshoe")

#### SIM DATASET V3 ####

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
data_sim=data
plot(data_sim[,12])
for (i in 1:79){
  if (y[i]==3)
    data_sim$`FOXO3a_pS318_S321-R-C`[i]=
      1+data_sim$`FOXO3a_pS318_S321-R-C`[i]
}
for (i in 80:158){
  if (y[i]==3){
    data_sim$newTRT[i]=1
    data_sim$`FOXO3a_pS318_S321-R-C`[i]=
      1+data_sim$`FOXO3a_pS318_S321-R-C`[i]
  }
}
plot(data_sim[,12])


data_test_sim1=data_sim
data_test_sim1$newTRT=1
data_test_sim0=data_sim
data_test_sim0$newTRT=0

data_new=rbind(data_test_sim0, data_test_sim1)

x_sim=mbart(data_sim, y, data_new, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
            base=0.95, power=2, sparse = FALSE)

mat_for_Y_0=x_sim$prob.test[,1:474]
mat_for_Y_1=x_sim$prob.test[,475:948]

delta_j_i=matrix(0,5000,474)

for(ii in 1:5000){
  for(j in 1:3){
    for(i in 1:158){
      delta_j_i[ii,3*(i-1)+j]=compute_delta(mat_for_Y_1[ii,(3*(i-1)+1):(3*(i-1)+3)],mat_for_Y_0[ii,(3*(i-1)+1):(3*(i-1)+3)],j)
    }
  }
}
ind0=seq(from = 1, to = 472, by = 3)
ind1=seq(from = 2, to = 473, by = 3)
ind2=seq(from = 3, to = 474, by = 3)

delta_0=delta_j_i[,ind0]
delta_1=delta_j_i[,ind1]
delta_2=delta_j_i[,ind2]

credible_intervals_delta_0 <- apply(delta_0, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_0 <- t(credible_intervals_delta_0)
credible_intervals_delta_0

credible_intervals_delta_1 <- apply(delta_1, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_1 <- t(credible_intervals_delta_1)
credible_intervals_delta_1

credible_intervals_delta_2 <- apply(delta_2, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_2 <- t(credible_intervals_delta_2)
credible_intervals_delta_2


plot(credible_intervals_delta_1[,1], ylim=c(-0.2,0.5))
points(credible_intervals_delta_1[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_1[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_1[, 2], col =1+as.integer(data_sim[,12]>0.5))

plot(credible_intervals_delta_2[,1], ylim=c(-0.2,0.5))
points(credible_intervals_delta_2[,2])
segments(x0 = 1:158, y0 = credible_intervals_delta_2[, 1], 
         x1 = 1:158, y1 = credible_intervals_delta_2[, 2], col =1+as.integer(data_sim[,12]>0.5))

tau_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    tau_i[i,j]=probabilities_tau(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
#credible_intervals_tau

eta_i=matrix(0,5000,158)

for(i in 1:5000){
  for(j in 1:158){
    eta_i[i,j]=probabilities_eta(mat_for_Y_1[i,(3*(j-1)+1):(3*(j-1)+3)],mat_for_Y_0[i,(3*(j-1)+1):(3*(j-1)+3)])
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

credible_intervals_eta[1,]
credible_intervals_tau[1,]

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
segments(x0 = 1:158, y0 = credible_intervals_eta[, 1], 
         x1 = 1:158, y1 = credible_intervals_eta[, 2], col =1+as.integer(data_sim[,12]>0.5))


plot(credible_intervals_tau[,1], ylim=c(0,1))
points(credible_intervals_tau[,2])
segments(x0 = 1:158, y0 = credible_intervals_tau[, 1], 
         x1 = 1:158, y1 = credible_intervals_tau[, 2], col =1+as.integer(data_sim[,12]>0.5))

#### New model: bayesreg v1 ####

library(bayesreg)

data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y = as.factor(as.integer(as.character(matchRTComp[,9]))+1)

y_bin1 = as.factor(ifelse(y==1,1,0))
fit1 = bayesreg(y_bin1~.,data,n.samples=1500,model="logistic",prior="horseshoe")

y_bin2 = as.factor(ifelse(y==2,1,0))
fit2 = bayesreg(y_bin2~.,data,n.samples=1500,model="logistic",prior="horseshoe")

y_bin3 = as.factor(ifelse(y==3,1,0))
fit3 = bayesreg(y_bin3~.,data,n.samples=1500,model="logistic",prior="horseshoe")

p_i1 = predict (fit1, newdata=data, type="response")
p_i2 = predict (fit2, newdata=data, type="response")
#p_i3 = 1 - p_i1 - p_i2
p_i3 = predict (fit3, newdata=data, type="response")

tot_prob = p_i1 + p_i2 + p_i3
p_i1=p_i1/tot_prob
p_i2=p_i2/tot_prob
p_i3=p_i3/tot_prob

predicted_probs = data.frame(p_i1=p_i1,p_i2=p_i2,p_i3=p_i3)
predicted_probs = apply(predicted_probs, 2, function(x) pmax(0,pmin(1,x)))

#head(predicted_probs)

chosen_class <- apply(predicted_probs, 1, which.max)

table(y)
table(chosen_class)
table(chosen_class,y)

sum(y==chosen_class)

#### New model: bayesreg v2 ####

library(bayesreg)

data=matchRTComp
data=data[,-c(1,2,3,4,9)]

y = as.factor(as.integer(as.character(matchRTComp[,9]))+1)

#standardize
data[,-c(1,2,5)] <- scale(data[,-c(1,2,5)])

y_bin1 = as.factor(ifelse(y==1,0,1))
fit1 = bayesreg(y_bin1~.,data,n.samples=1500,model="logistic",prior="horseshoe")

y2 = y[y_bin1==1]
y_bin2 = as.factor(ifelse(y2==2,0,1))
fit2 = bayesreg(y_bin2~.,data[y_bin1==1,],n.samples=1500,model="logistic",prior="horseshoe")

#rho
rho1 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho1[i,] = 1/ ( 1 + exp(-fit1$beta0[i] - as.matrix(data)%*%as.matrix(fit1$beta[,i])))
}

rho2 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho2[i,] = 1/ ( 1 + exp(-fit2$beta0[i] - as.matrix(data)%*%as.matrix(fit2$beta[,i])))
}

#p{i,j}
pi1 <- 1-rho1
pi2 <- rho1*(1-rho2)
pi3 <- rho1*rho2

predicted_class <- matrix(0, 1500, 158)

for (i in 1:1500) {
  for (j in 1:158) {
    max_value <- max(pi1[i, j], pi2[i, j], pi3[i, j])
    
    if (max_value == pi1[i, j]) {
      predicted_class[i, j] <- 1
    } else if (max_value == pi2[i, j]) {
      predicted_class[i, j] <- 2
    } else {
      predicted_class[i, j] <- 3
    }
  }
}

pred_class <- apply(predicted_class, 2, function(col) {
  as.numeric(names(which.max(table(col))))
})

sum(y==pred_class)
table(pred_class, y)

table(y)
table(pred_class)


#causal effects
data_test1=data
data_test1$newTRT=1
data_test0=data
data_test0$newTRT=0

data_new=rbind(data_test0, data_test1)

rho1 <- matrix (NA, 1500, 316)
for (i in 1:1500){
  rho1[i,] = 1/ ( 1 + exp(-fit1$beta0[i] - as.matrix(data_new)%*%as.matrix(fit1$beta[,i])))
}

rho2 <- matrix (NA, 1500, 316)
for (i in 1:1500){
  rho2[i,] = 1/ ( 1 + exp(-fit2$beta0[i] - as.matrix(data_new)%*%as.matrix(fit2$beta[,i])))
}

pi1 <- 1-rho1
pi2 <- rho1*(1-rho2)
pi3 <- rho1*rho2

mat_for_Y_0 <- cbind(pi1[,1:158],pi2[,1:158],pi3[,1:158])
mat_for_Y_1 <- cbind(pi1[,159:316],pi2[,159:316],pi3[,159:316])

eta_i=matrix(0,1500,158)

for(i in 1:1500){
  for(j in 1:158){
    eta_i[i,j] = mat_for_Y_0[i,j]*(mat_for_Y_1[i,158+j]+mat_for_Y_1[i,158*2+j]) + mat_for_Y_0[i,158+j]*mat_for_Y_1[i,158*2+j]
  }
}

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
#credible_intervals_eta

plot(credible_intervals_eta[,1], ylim=c(0,1))
points(credible_intervals_eta[,2])
abline(h=0.28)

#### Code to check if the code is overall correct ####

#### New model: bayesreg v3 ####

# Here we study 2 different models for newTrT = 0 and newTrT = 1

library(bayesreg)

data=matchRTComp
data=data[,-c(1,2,3,4,9)]

data_for_pscore = data[,-5]
pscore = glm(data$newTRT ~ ., data=data_for_pscore, family = binomial)$fitted.values
data$pscore <- pscore

y = as.factor(as.integer(as.character(matchRTComp[,9]))+1)

#standardize
data[,-c(1,2,5)] <- scale(data[,-c(1,2,5)])

data_newTRT_1 = data[1:79,-5]
data_newTRT_0 = data[80:158,-5]
data_fit = data[,-5]

# innovative treatment

y_1 = y[1:79]

y_bin1 = as.factor(ifelse(y_1==1,0,1))
fit1 = bayesreg(y_bin1~.,data_newTRT_1,n.samples=1500,model="logistic",prior="horseshoe")
for (i in 1:10){
  matplot(fit1$beta[i,],type='l')
}

y2 = y_1[y_bin1==1]
y_bin2 = as.factor(ifelse(y2==2,0,1))
fit2 = bayesreg(y_bin2~.,data_newTRT_1[y_bin1==1,],n.samples=1500,model="logistic",prior="horseshoe")
for (i in 1:10){
  matplot(fit2$beta[i,],type='l')
}

#rho
rho1 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho1[i,] = 1/ ( 1 + exp(-fit1$beta0[i] - as.matrix(data_fit)%*%as.matrix(fit1$beta[,i])))
}
for (i in 1:10){
  matplot(rho1[,i],type='l')
}

rho2 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho2[i,] = 1/ ( 1 + exp(-fit2$beta0[i] - as.matrix(data_fit)%*%as.matrix(fit2$beta[,i])))
}
for (i in 1:10){
  matplot(rho2[,i],type='l')
}

#p{i,j}
pi1 <- 1-rho1
pi2 <- rho1*(1-rho2)
pi3 <- rho1*rho2

mat_for_Y_1 <- cbind(pi1,pi2,pi3)

# traditional treatment

y_0 = y[80:158]

y_bin1 = as.factor(ifelse(y_0==1,0,1))
fit1 = bayesreg(y_bin1~.,data_newTRT_0,n.samples=1500,model="logistic",prior="horseshoe")
for (i in 1:10){
  matplot(fit1$beta[i,],type='l')
}

y2 = y_0[y_bin1==1]
y_bin2 = as.factor(ifelse(y2==2,0,1))
fit2 = bayesreg(y_bin2~.,data_newTRT_0[y_bin1==1,],n.samples=1500,model="logistic",prior="horseshoe")

#rho
rho1 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho1[i,] = 1/ ( 1 + exp(-fit1$beta0[i] - as.matrix(data_fit)%*%as.matrix(fit1$beta[,i])))
}
for (i in 1:10){
  matplot(rho1[,i],type='l')
}

rho2 <- matrix (NA, 1500, 158)
for (i in 1:1500){
  rho2[i,] = 1/ ( 1 + exp(-fit2$beta0[i] - as.matrix(data_fit)%*%as.matrix(fit2$beta[,i])))
}
for (i in 1:10){
  matplot(rho2[,i],type='l')
}

#p{i,j}
pi1 <- 1-rho1
pi2 <- rho1*(1-rho2)
pi3 <- rho1*rho2

mat_for_Y_0 <- cbind(pi1,pi2,pi3)

#propensity score
n_groups=4
quantiles <- quantile(pscore, probs = seq(0, 1, length.out = n_groups + 1))
groups <- cut(pscore, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
colours <- rainbow(n_groups)
col=colours[groups]
col=col[order(pscore)]

# Find the indexes for each group
gr_1=which(groups==1)
gr_2=which(groups==2)
gr_3=which(groups==3)
gr_4=which(groups==4)


eta_i=matrix(0,1500,158)

for(i in 1:1500){
  for(j in 1:158){
    eta_i[i,j] = mat_for_Y_0[i,j]*(mat_for_Y_1[i,158+j]+mat_for_Y_1[i,158*2+j]) + mat_for_Y_0[i,158+j]*mat_for_Y_1[i,158*2+j]
  }
}

tau_i=matrix(0,1500,158)

for(i in 1:1500){
  for(j in 1:158){
    tau_i[i,j] = mat_for_Y_0[i,j]*(mat_for_Y_1[i,j]+mat_for_Y_1[i,158+j]+mat_for_Y_1[i,158*2+j]) + mat_for_Y_0[i,158+j]*(mat_for_Y_1[i,158+j]+mat_for_Y_1[i,158*2+j]) + mat_for_Y_0[i,158*2+j]*mat_for_Y_1[i,158*2+j]
  }
}

credible_intervals_tau <- apply(tau_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_tau <- t(credible_intervals_tau)
plot(credible_intervals_tau[order(pscore),1], ylim=c(0,1))
points(credible_intervals_tau[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_tau[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_tau[order(pscore), 2], col = col)

credible_intervals_eta <- apply(eta_i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_eta <- t(credible_intervals_eta)
plot(credible_intervals_eta[order(pscore),1], ylim=c(0,1))
points(credible_intervals_eta[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_eta[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_eta[order(pscore), 2], col = col)

delta_1i <- matrix (0, 1500, 158)
for(i in 1:1500){
  for(j in 1:158){
    p0 = mat_for_Y_0[i,158+j] + mat_for_Y_0[i, 158*2+j]
    p1 = mat_for_Y_1[i,158+j] + mat_for_Y_1[i, 158*2+j]
    delta_1i[i,j] = p1-p0
  }
}

delta_2i <- matrix (0, 1500, 158)
for(i in 1:1500){
  for(j in 1:158){
    p0 = mat_for_Y_0[i, 158*2+j]
    p1 = mat_for_Y_1[i, 158*2+j]
    delta_2i[i,j] = p1-p0
  }
}


credible_intervals_delta_1 <- apply(delta_1i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_1 <- t(credible_intervals_delta_1)
credible_intervals_delta_1

credible_intervals_delta_2 <- apply(delta_2i, 2, quantile, probs = c(0.025, 0.975))
credible_intervals_delta_2 <- t(credible_intervals_delta_2)
credible_intervals_delta_2


plot(credible_intervals_delta_1[order(pscore),1], ylim=c(-1,1))
points(credible_intervals_delta_1[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_delta_1[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_delta_1[order(pscore), 2], col = col)

plot(credible_intervals_delta_2[order(pscore),1], ylim=c(-1,1))
points(credible_intervals_delta_2[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_delta_2[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_delta_2[order(pscore), 2], col = col)

#mean by groups
credible_intervals_tau_groups <- aggregate(credible_intervals_tau, by = list(groups), FUN = mean)
credible_intervals_eta_groups <- aggregate(credible_intervals_eta, by = list(groups), FUN = mean)
credible_intervals_delta_1_groups <- aggregate(credible_intervals_delta_1, by = list(groups), FUN = mean)
credible_intervals_delta_2_groups <- aggregate(credible_intervals_delta_2, by = list(groups), FUN = mean)

plot(credible_intervals_tau_groups[,2], ylim=c(0,1))
points(credible_intervals_tau_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_tau_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_tau_groups[, 3], col = colours)

plot(credible_intervals_eta_groups[,2], ylim=c(0,1))
points(credible_intervals_eta_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_eta_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_eta_groups[, 3], col = colours)

plot(credible_intervals_delta_1_groups[,2], ylim=c(-0.5,0.5))
points(credible_intervals_delta_1_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_delta_1_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_delta_1_groups[, 3], col = colours)

plot(credible_intervals_delta_2_groups[,2], ylim=c(-0.5,0.5))
points(credible_intervals_delta_2_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_delta_2_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_delta_2_groups[, 3], col = colours)

credible_intervals_eta_groups
credible_intervals_tau_groups

#Useful plots

plot(credible_intervals_tau[,1], ylim=c(0,1), ylab="",xlab="", main="Credible intervals for tau / T-learner")
points(credible_intervals_tau[,2])
segments(x0 = 1:158, y0 = credible_intervals_tau[, 1], 
         x1 = 1:158, y1 = credible_intervals_tau[, 2])

plot(credible_intervals_tau[order(pscore),1], ylim=c(0,1), ylab="",xlab="", 
     main="Credible intervals for tau / T-learner \n grouped by pscore")
points(credible_intervals_tau[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_tau[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_tau[order(pscore), 2], col = col)

plot(credible_intervals_eta[,1], ylim=c(0,1), ylab="",xlab="", main="Credible intervals for eta / T-learner")
points(credible_intervals_eta[,2])
segments(x0 = 1:158, y0 = credible_intervals_eta[, 1], 
         x1 = 1:158, y1 = credible_intervals_eta[, 2])

plot(credible_intervals_eta[order(pscore),1], ylim=c(0,1), ylab="",xlab="", 
     main="Credible intervals for eta / T-learner \n grouped by pscore")
points(credible_intervals_eta[order(pscore),2])
segments(x0 = 1:158, y0 = credible_intervals_eta[order(pscore), 1], 
         x1 = 1:158, y1 = credible_intervals_eta[order(pscore), 2], col = col)

plot(credible_intervals_tau_groups[,2], ylim=c(0,1))
points(credible_intervals_tau_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_tau_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_tau_groups[, 3], col = colours)

plot(credible_intervals_eta_groups[,2], ylim=c(0,1))
points(credible_intervals_eta_groups[,3])
segments(x0 = 1:4, y0 = credible_intervals_eta_groups[, 2], 
         x1 = 1:4, y1 = credible_intervals_eta_groups[, 3], col = colours)


#### Characterization of the 4 pscore groups ####

library(BART)
data=matchRTComp
data=data[,-c(1,2,3,4,9)]

data_for_pscore = data[,-5]
y=as.integer(as.integer(as.character(matchRTComp[,9]))+1)
mod1=glm(data$newTRT ~ ., data=data_for_pscore, family = binomial)
pscore = mod1$fitted.values

n_groups=4
quantiles <- quantile(pscore, probs = seq(0, 1, length.out = n_groups + 1))
groups <- cut(pscore, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
colours <- rainbow(n_groups)
col=colours[groups]
col=col[order(pscore)]

# Find the indexes for each group
gr_1=which(groups==1)
gr_2=which(groups==2)
gr_3=which(groups==3)
gr_4=which(groups==4)

summary(mod1)

mod2=glm(data$newTRT ~ . - `ER-alpha-R-V` - `HER2_pY1248-R-C` - patients.tumor_grade
         - `Bad_pS112-R-V` - `PAI-1-M-E` - `Src_pY416-R-C` - `RBM15-R-V` - `Bcl-2-M-V`
         - `Paxillin-R-C` - `CD31-M-V` - `14-3-3_epsilon-M-C` - `HSP70-R-C` - patients.age_at_initial_pathologic_diagnosis
         - `YAP_pS127-R-E` - `Ku80-R-C` - `Caspase-7_cleavedD198-R-C` - `Rab25-R-V` - `Lck-R-V` 
         - `PRAS40_pT246-R-V` - `Smad1-R-V` - `C-Raf-R-V` - `ACVRL1-R-C` - `Akt_pS473-R-V` - patients.gender
         - `FOXO3a_pS318_S321-R-C` - inYear - `SF2-M-V` - `Cyclin_E2-R-C` - `GSK3-alpha-beta-M-V`
           ,data=data_for_pscore, family = binomial)
summary(mod2)

# I keep `MYH11-R-V`, `Caveolin-1-R-V` and `Claudin-7-R-V`

# Combine the groups
ordered_indices <- c(gr_1, gr_2, gr_3, gr_4)

# Choose a single row

##### MYH11-R-V #####
variable_index <- "MYH11-R-V"
values <- data[ordered_indices, variable_index]

mean_1 <- mean(values[1:length(gr_1)])
mean_2 <- mean(values[(length(gr_1) + 1):(length(gr_1) + length(gr_2))])
mean_3 <- mean(values[(length(gr_1) + length(gr_2) + 1):(length(gr_1) + length(gr_2) + length(gr_3))])
mean_4 <- mean(values[(length(gr_1) + length(gr_2) + length(gr_3) + 1):158])

# Plot points
plot(values, type = "p", col = rep(1:4, times = c(length(gr_1), length(gr_2), length(gr_3), length(gr_4))),
     xlab = "Individuals (Ordered by Group)", ylab = "Variable Value",
     main = paste("Variable", variable_index, "Across Groups"))
#legend("topright", legend = paste("Group", 1:4), col = 1:4, pch = 1)

# Define group boundaries
boundaries <- c(1, length(gr_1), length(gr_1) + length(gr_2), 
                length(gr_1) + length(gr_2) + length(gr_3), 158)

# Add horizontal mean lines
segments(boundaries[1], mean_1, boundaries[2], mean_1, col = 1, lwd = 2, lty = 1)
segments(boundaries[2] + 1, mean_2, boundaries[3], mean_2, col = 2, lwd = 2, lty = 1)
segments(boundaries[3] + 1, mean_3, boundaries[4], mean_3, col = 3, lwd = 2, lty = 1)
segments(boundaries[4] + 1, mean_4, boundaries[5], mean_4, col = 4, lwd = 2, lty = 1)

##### Caveolin-1-R-V #####
variable_index <- "Caveolin-1-R-V"
values <- data[ordered_indices, variable_index]

mean_1 <- mean(values[1:length(gr_1)])
mean_2 <- mean(values[(length(gr_1) + 1):(length(gr_1) + length(gr_2))])
mean_3 <- mean(values[(length(gr_1) + length(gr_2) + 1):(length(gr_1) + length(gr_2) + length(gr_3))])
mean_4 <- mean(values[(length(gr_1) + length(gr_2) + length(gr_3) + 1):158])

# Plot points
plot(values, type = "p", col = rep(1:4, times = c(length(gr_1), length(gr_2), length(gr_3), length(gr_4))),
     xlab = "Individuals (Ordered by Group)", ylab = "Variable Value",
     main = paste("Variable", variable_index, "Across Groups"))
#legend("topright", legend = paste("Group", 1:4), col = 1:4, pch = 1)

# Define group boundaries
boundaries <- c(1, length(gr_1), length(gr_1) + length(gr_2), 
                length(gr_1) + length(gr_2) + length(gr_3), 158)

# Add horizontal mean lines
segments(boundaries[1], mean_1, boundaries[2], mean_1, col = 1, lwd = 2, lty = 1)
segments(boundaries[2] + 1, mean_2, boundaries[3], mean_2, col = 2, lwd = 2, lty = 1)
segments(boundaries[3] + 1, mean_3, boundaries[4], mean_3, col = 3, lwd = 2, lty = 1)
segments(boundaries[4] + 1, mean_4, boundaries[5], mean_4, col = 4, lwd = 2, lty = 1)

##### Claudin-7-R-V #####
variable_index <- "Claudin-7-R-V"
values <- data[ordered_indices, variable_index]

mean_1 <- mean(values[1:length(gr_1)])
mean_2 <- mean(values[(length(gr_1) + 1):(length(gr_1) + length(gr_2))])
mean_3 <- mean(values[(length(gr_1) + length(gr_2) + 1):(length(gr_1) + length(gr_2) + length(gr_3))])
mean_4 <- mean(values[(length(gr_1) + length(gr_2) + length(gr_3) + 1):158])

# Plot points
plot(values, type = "p", col = rep(1:4, times = c(length(gr_1), length(gr_2), length(gr_3), length(gr_4))),
     xlab = "Individuals (Ordered by Group)", ylab = "Variable Value",
     main = paste("Variable", variable_index, "Across Groups"))
#legend("topright", legend = paste("Group", 1:4), col = 1:4, pch = 1)

# Define group boundaries
boundaries <- c(1, length(gr_1), length(gr_1) + length(gr_2), 
                length(gr_1) + length(gr_2) + length(gr_3), 158)

# Add horizontal mean lines
segments(boundaries[1], mean_1, boundaries[2], mean_1, col = 1, lwd = 2, lty = 1)
segments(boundaries[2] + 1, mean_2, boundaries[3], mean_2, col = 2, lwd = 2, lty = 1)
segments(boundaries[3] + 1, mean_3, boundaries[4], mean_3, col = 3, lwd = 2, lty = 1)
segments(boundaries[4] + 1, mean_4, boundaries[5], mean_4, col = 4, lwd = 2, lty = 1)

##### Extra #####
variable_index <- 33
values <- data[ordered_indices, variable_index]

mean_1 <- mean(values[1:length(gr_1)])
mean_2 <- mean(values[(length(gr_1) + 1):(length(gr_1) + length(gr_2))])
mean_3 <- mean(values[(length(gr_1) + length(gr_2) + 1):(length(gr_1) + length(gr_2) + length(gr_3))])
mean_4 <- mean(values[(length(gr_1) + length(gr_2) + length(gr_3) + 1):158])

# Plot points
plot(values, type = "p", col = rep(1:4, times = c(length(gr_1), length(gr_2), length(gr_3), length(gr_4))),
     xlab = "Individuals (Ordered by Group)", ylab = "Variable Value",
     main = paste("Variable", variable_index, "Across Groups"))
#legend("topright", legend = paste("Group", 1:4), col = 1:4, pch = 1)

# Define group boundaries
boundaries <- c(1, length(gr_1), length(gr_1) + length(gr_2), 
                length(gr_1) + length(gr_2) + length(gr_3), 158)

# Add horizontal mean lines
segments(boundaries[1], mean_1, boundaries[2], mean_1, col = 1, lwd = 2, lty = 1)
segments(boundaries[2] + 1, mean_2, boundaries[3], mean_2, col = 2, lwd = 2, lty = 1)
segments(boundaries[3] + 1, mean_3, boundaries[4], mean_3, col = 3, lwd = 2, lty = 1)
segments(boundaries[4] + 1, mean_4, boundaries[5], mean_4, col = 4, lwd = 2, lty = 1)

# Relevant ones: 6_X, 12_OK, 20!_OK, 23!_OK, 26_X, 27!_OK, 28_X, 30_X, 33_X

# ANOVA on groups
variable_index <- 33
values <- data[ordered_indices,variable_index]
# Creazione del data frame con i gruppi
group_labels <- rep(1:4, times = c(length(gr_1), length(gr_2), length(gr_3), length(gr_4)))
df_anova <- data.frame(Values = values, Group = as.factor(group_labels))
# Esegui l'ANOVA
anova_result <- aov(Values ~ Group, data = df_anova)
# Riassunto dei risultati
summary(anova_result)
TukeyHSD(anova_result)

# Relevant ones: 12, 20!, 23!, 27!, 30?, 33? OKOK

# fare boxplot

##### Radarchart plot #####

library(fmsb)

Data=as.data.frame(scale(data))

# Define covariate indices
cov_indices <- c(12, 20, 23, 27, 30, 33)

# Compute overall mean for selected covariates
overall_mean <- colMeans(Data[, cov_indices])

# Function to compute mean for a given group
compute_group_mean <- function(group_indices) {
  colMeans(Data[group_indices, cov_indices])
}

# Compute mean for each group
group_means <- list(
  "Group 1" = compute_group_mean(gr_1),
  "Group 2" = compute_group_mean(gr_2),
  "Group 3" = compute_group_mean(gr_3),
  "Group 4" = compute_group_mean(gr_4)
)

# Prepare data for radar plot
max_vals <- apply(Data[, cov_indices], 2, max)  # Maximum values for scaling
min_vals <- apply(Data[, cov_indices], 2, min)  # Minimum values for scaling

# Set up a 2x2 grid layout for the plots
#par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
par(mfrow = c(2, 2), mar = c(1, 1, 2, 1), oma = c(2, 2, 2, 2))  # Reduce margins, add outer margins

# Define colors (Overall in black, Groups in different colors)
group_colors <- c(1,2,3,4)

# Plot each group vs. overall mean
for (i in seq_along(group_means)) {
  plot_data <- as.data.frame(rbind(max_vals, min_vals, overall_mean, group_means[[i]]))
  rownames(plot_data) <- c("Max", "Min", "Overall", names(group_means)[i])
  
  radarchart(plot_data, 
             axistype = 0, 
             pcol = c("gray", group_colors[i]),   # Black for overall, colored for group
             plwd = 3, 
             plty = 1,
             cglcol = "gray",  # Makes grid less distracting
             cglwd = 0.5,  # Thinner grid lines
             vlcex = 0.7,  # Makes labels (variable names) bigger
             title = paste("Comparison:", names(group_means)[i]) # Title for each plot
             #centerzero = TRUE, # Expands hexagons to fill more space)  
             #caxislabels = c(0, 0.5, 1, 1.5) # Increase axis range to stretch out the hexagons) # Increase radius to make hexagons fill more space
             )  
  
  # Add a small legend inside each plot
  #legend("topright", legend = c("Overall Mean", names(group_means)[i]), 
  #       col = c("black", group_colors[i]), lty = 1, lwd = 2, cex = 0.8)
}

# Reset plot layout
par(mfrow = c(1,1)) 

###

group_means <- list(
  "Group 2" = compute_group_mean(gr_2),
  "Group 3" = compute_group_mean(gr_3)
)

# Prepare data for radar plot
max_vals <- apply(Data[, cov_indices], 2, max)  # Maximum values for scaling
min_vals <- apply(Data[, cov_indices], 2, min)  # Minimum values for scaling

# Set up a 1x2 grid layout for the plots
par(mfrow = c(1, 2))

# Define colors (Overall in black, Groups in different colors)
group_colors <- c(2,3)

# Plot each group vs. overall mean
for (i in seq_along(group_means)) {
  plot_data <- as.data.frame(rbind(max_vals, min_vals, overall_mean, group_means[[i]]))
  rownames(plot_data) <- c("Max", "Min", "Overall", names(group_means)[i])
  
  radarchart(plot_data, 
             axistype = 0, 
             pcol = c("gray", group_colors[i]),   # Black for overall, colored for group
             plwd = 3, 
             plty = 1,
             cglcol = "gray",  # Makes grid less distracting
             cglwd = 1.5,  # Thinner grid lines
             vlcex = 0.9,  # Makes labels (variable names) bigger
             title = paste("Comparison:", names(group_means)[i]) # Title for each plot
             )  
}

# Reset plot layout
par(mfrow = c(1,1)) 


# Add a small legend inside the plot
legend("toplet", legend = c("Overall Mean", names(group_means)[1], names(group_means)[2]), 
       col = c("black", group_colors[1],group_colors[2] ), lty = 1, lwd = 2, cex = 0.8)


#### old stuff ####



prob_test_mean2=x$prob.test.mean
prob_test_mean2=matrix(prob_test_mean2, nrow=316, byrow=TRUE)

pred_class2 <- apply(prob_test_mean2, 1, which.max)
Y_0=pred_class2[1:158]
Y_1=pred_class2[159:316]
table(Y_0)
table(Y_1)

delta0=sum(Y_1>=1)/158-sum(Y_0>=1)/158
delta0
delta1=sum(Y_1>=2)/158-sum(Y_0>=2)/158
delta1
delta2=sum(Y_1>=3)/158-sum(Y_0>=3)/158
delta2

sort(Y_1-Y_0)



#### old wrong stuff #####

data2=data
trt_to_swap=data2$newTRT
trt_to_swap=1-trt_to_swap
data2$newTRT=trt_to_swap

x2=mbart(data, y, data2, ndpost=5000, nskip=1000, printevery=5000, ntree=100,
           base=0.95, power=2, sparse = FALSE)

prob_test2=x2$prob.test

tau_i=numeric(158)

for (i in 1:79){
  Y_1=y[i]
  pred=prob_test2[,(3*(i-1)+1):(3*(i-1)+3)]
  pred <- apply(pred, 1, which.max)
  tot=Y_1>=pred
  tau_i[i]=sum(tot)/5000
}
for (i in 80:158){
  Y_0=y[i]
  pred=prob_test2[,(3*(i-1)+1):(3*(i-1)+3)]
  pred <- apply(pred, 1, which.max)
  tot=pred>=Y_0
  tau_i[i]=sum(tot)/5000
}
tau_i

eta_i=numeric(158)

for (i in 1:79){
  Y_1=y[i]
  pred=prob_test2[,(3*(i-1)+1):(3*(i-1)+3)]
  pred <- apply(pred, 1, which.max)
  tot=Y_1>pred
  eta_i[i]=sum(tot)/5000
}
for (i in 80:158){
  Y_0=y[i]
  pred=prob_test2[,(3*(i-1)+1):(3*(i-1)+3)]
  pred <- apply(pred, 1, which.max)
  tot=pred>Y_0
  eta_i[i]=sum(tot)/5000
}
eta_i

class_colors <- c("red", "blue", "green")

plot(eta_i, col = class_colors[y])
legend("topright", legend = c("Class 1", "Class 2", "Class 3"), 
       col = class_colors, pch = 19)
abline(h=0.5)
plot(sort(eta_i), col = class_colors[y])
legend("topleft", legend = c("Class 1", "Class 2", "Class 3"), 
       col = class_colors, pch = 19)
abline(h=0.8)
abline(h=0.5)

patients_to_trt_80=which(eta_i>0.8)
y[patients_to_trt_80]
x2$prob.test.mean[(3*(104-1)+1):(3*(104-1)+3)]
x2$prob.test.mean[(3*(129-1)+1):(3*(129-1)+3)]
x2$prob.test.mean[(3*(131-1)+1):(3*(131-1)+3)]

patients_to_trt_50=which(eta_i>0.5)
y[patients_to_trt_50]




#### don't run these ones ####

head(x$varprob.mean)

library(coda)
#mcmc_varcount1=as.mcmc(x$varcount[[1]])
#par(mar=c(4,4,2,1))
#plot(mcmc_varcount1, main="MCMC chains for variable counts (contrast 1)")

#mcmc_varcount2=as.mcmc(x$varcount[[2]])
#plot(mcmc_varcount2, main="MCMC chains for variable counts (contrast 2)")

mcmc_treedraws=as.mcmc(x$treedraws)
plot(mcmc_treedraws, main="MCMC chains for the tree draws (depths/counts)")

print(x$varprob.mean)

hist(x$varprob.mean[,1], main="Predicted probability for class 0")
hist(x$varprob.mean[,2], main="Predicted probability for class 1")
hist(x$varprob.mean[,3], main="Predicted probability for class 2")

#____
prob_train=x$prob.test
pred_class_train=apply(prob_train, 1, function(x) which.max(x))
pred_class_majority=apply(pred_class_train,2,function(x) names(sort(table(x), decreasing=TRUE)[1]))

accuracy=mean(pred_class_majority==y)
print(paste("Training accuracy:", accuracy))

help(mbart)