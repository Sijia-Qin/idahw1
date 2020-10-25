
require(MASS)
n <- 500
Sigma <- matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), nrow = 3, byrow = T)
#simulating the data
set.seed(1)
Z <- mvrnorm(n = n, mu = c(0,0,0), Sigma = Sigma)
Z1 <- Z[,1]; Z2 <- Z[,2]; Z3 <- Z[,3]
Y1= 1+Z1
Y2= 5 + 2*Z1 + Z2
a=2
b=0
r_ = a*(Y1-1) + b * (Y2-5) + Z3
ind_ <- which(r_ >= 1)
Y2_obs=Y2[ind_]
Y2_mis=Y2[-ind_]
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y2),
     main = "MAR", ylim = c(0, 0.7))
lines(density(Y2_obs), lwd = 2, col = "red")
lines(density(Y2_mis), lwd = 2, col = "darkgreen")
legend(1, 0.7, legend = c("Complete data", "Observed data", "Missing data"),
       col = c("blue", "red", "darkgreen"), lty = c(1,1,1), lwd = c(2,2,2), bty ="n")
#replace the missing value by NA
Y2_new=Y2
for (i in (1:500)){
  if((2*(Y1[i]-1) + Z3[i]) <0){
    Y2_new[i]=NA
  }
}
data=data.frame(Y1,Y2_new)
fit = lm(Y2_new ~ Y1, data = data)
predicted_ = predict(fit, newdata = data) + rnorm(nrow(data), 0, sigma(fit))
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y2),
     main=" Stochastic" ,ylim = c(0, 0.7))
lines(density(predicted_), lwd = 2, col = "red")
legend(1, 0.7, legend = c("Original data", "Complete data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")



a=0
b=2
r_ = a*(Y1-1) + b * (Y2-5) + Z3
ind_ <- which(r_ >= 1)
Y2_obs=Y2[ind_]
Y2_mis=Y2[-ind_]
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y2),
     main = "MNAR", ylim = c(0, 0.7))
lines(density(Y2_obs), lwd = 2, col = "red")
lines(density(Y2_mis), lwd = 2, col = "darkgreen")
legend(1, 0.7, legend = c("Complete data", "Observed data", "Missing data"),
       col = c("blue", "red", "darkgreen"), lty = c(1,1,1), lwd = c(2,2,2), bty ="n")
#replace the missing value by NA
Y2_new=Y2
for (i in (1:500)){
  if((2*(Y1[i]-1) + Z3[i]) <0){
    Y2_new[i]=NA
  }
}
data=data.frame(Y1,Y2_new)
fit = lm(Y2_new ~ Y1, data = data)
predicted_ = predict(fit, newdata = data) + rnorm(nrow(data), 0, sigma(fit))
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y2),
     main=" Stochastic" ,ylim = c(0, 0.7))
lines(density(predicted_), lwd = 2, col = "red")
legend(1, 0.7, legend = c("Original data", "Complete data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")









databp = read.table("databp.txt", header = TRUE)
#compelete value
ind <- which(is.na(databp$recovtime) == FALSE) 
mccoverall <- mean(databp$recovtime, na.rm = TRUE)
seccoverall <- sd(databp$recovtime, na.rm = TRUE)/sqrt(length(ind))
mccoverall; seccoverall
cor(databp$recovtime, databp$logdose, use = "complete")
cor(databp$recovtime, databp$bloodp, use = "complete")

#mean imputation
mean(databp$recovtime, na.rm = TRUE)
recovtime_mi <- ifelse(is.na(databp$recovtime), 19.27273, databp$recovtime)
sd(recovtime_mi)
cor(recovtime_mi,databp$logdose)
cor(recovtime_mi,databp$bloodp)

#mean regression
fit <- lm(recovtime ~ logdose+bloodp, data = databp)
predicted_ri <- predict(fit, newdata = databp)
recovtime_ri <- ifelse(is.na(databp$recovtime), predicted_ri, databp$recovtime)
sd(recovtime_ri)
cor(recovtime_ri,databp$logdose)
cor(recovtime_ri,databp$bloodp)

#stochastic regression imputation
predicted_sri <- predict(fit, newdata = databp) + rnorm(nrow(databp), 0, sigma(fit))
recovtime_sri <- ifelse(is.na(databp$recovtime), predicted_sri, databp$recovtime)
sd(recovtime_sri)
cor(recovtime_sri,databp$logdose)
cor(recovtime_sri,databp$bloodp)

#predictive mean matching
predic=c(recovtime_ri[4],recovtime_ri[10],recovtime_ri[22])
ind_donors_1 <- predic
recovtime_=databp$recovtime[-4]
recovtime_=recovtime_[-9]
recovtime_=recovtime_[-20]
min1=sqrt(abs(ind_donors_1[1]^2-recovtime_[1]^2))
for (i in recovtime_){
    a=sqrt(abs(ind_donors_1[1]^2-i^2))
    if(a<min1){
      min1=a
      donor1=i
    }
}
min2=sqrt(abs(ind_donors_1[2]^2-recovtime_[1]^2))
for (i in recovtime_){
  a=sqrt(abs(ind_donors_1[2]^2-i^2))
  if(a<min2){
    min2=a
    donor2=i
  }
}
min3=sqrt(abs(ind_donors_1[3]^2-recovtime_[1]^2))
for (i in recovtime_){
  a=sqrt(abs(ind_donors_1[3]^2-i^2))
  if(a<min3){
    min3=a
    donor3=i
  }
}
recovtime_hotdeck <- c(recovtime_[1:3],donor1,recovtime_[4:8],donor2,recovtime_[9:19],donor3,recovtime_[20:22])
sd(recovtime_hotdeck)
cor(recovtime_hotdeck,databp$logdose)
cor(recovtime_hotdeck,databp$bloodp)

