# Import boston.txt dataset as 'data'

data=read.csv("C:/Users/Sayan/Downloads/boston.txt", sep="")
# Loading data and describing data

attach(data)

# y and X

y = MEDV
detach(data)
data = as.matrix(data)
n = length(data[, 1])
p = length(data[1, ]) - 3
X = cbind(rep(1,n), data[,c(-2,-4,-9,-14)])



# Standardizing the regressors ( X ---> Z )

Z = X
for(i in 2:p)
{
  Z[,i] = (X[, i] - mean(X[, i]))/sd(X[, i])
}


# Splitting the dataset ( X.train,y.train and X.test,y.test )

set.seed(99999)
all = 1:n
test_i = sort(sample(all , 100))

X.test = Z[test_i,]
y.test = y[test_i]
X.train = Z[-test_i,]
y.train = y[-test_i]


# OLS estimate of y.test ( y.test.ht )

beta.train.ht = solve(crossprod(X.train))%*%crossprod(X.train,y.train) 
y.test.ht = as.vector(X.test%*%beta.train.ht)

# RMSE

RMSE = as.vector(sqrt(crossprod(y.test - y.test.ht)/length(y.test)))

# Projection matrix (P_x.train) and Leverage pts ( Lev.pt )

P_x.train = X.train%*%solve(crossprod(X.train))%*%t(X.train)
hii = as.vector(diag(P_x.train))
Lev.pt = hii > 0.1
Lev.pt1 <- as.factor(Lev.pt)

## Leverage plot

par(mfrow=c(1,1))
plot(1:406, hii,xlab = "sample_points", pch = 16, cex = 0.8, col = c("black","blue")[Lev.pt1])

## New regressors and response without Leverage pts ( Z.new,y.new )

X.new = X.train[!Lev.pt,]
y.new = y.train[!Lev.pt]
n.new= length(y.new)

## Cook's Distance ( cooks.dist )

beta.new.ht = as.vector(solve(crossprod(X.new))%*%crossprod(X.new,y.new))
y.new.ht = as.vector(X.new%*%beta.new.ht)
ei.new = y.new - y.new.ht

RSS = as.vector(crossprod(ei.new))
MSres = RSS/(n.new-p)


beta.ht.i = matrix(rep(0,n.new*p),nrow = n.new)
for(i in 1:n.new)
{
  beta.ht.i[i,] = solve(crossprod(X.new[-i,]))%*%crossprod(X.new[-i,],y.new[-i])
}

C_i = y.new
for(i in 1:n.new)
{
  C_i[i] = t(beta.ht.i[i,] - beta.new.ht)%*%crossprod(X.new)%*%(beta.ht.i[i,] - beta.new.ht)/(p*MSres)  
}

cooks.dist = C_i > 1
ind.1 = as.numeric(as.factor( C_i > 1 ))

## DFBETA

RSS.i = rep(0,n.new)
for( i in 1:n.new )
{
  RSS.i[i] = RSS - ei.new[i]^2/(1-hii[i])
}
S_i.sq = RSS.i/(n.new-1-p)

Cjj = as.vector(diag(solve(crossprod(X.new))))

dfbeta.ji = matrix(rep(0,n.new*p), nrow = n.new)
for(i in 1:n.new)
{
  for(j in 1:p)
  {
    dfbeta.ji[i,j] = (beta.new.ht[j] - beta.ht.i[i,j])/(sqrt(S_i.sq[i]*Cjj[j]))
  }
}


## DFFITS


y.ht.i = rep(0,n.new)
for(i in 1:n.new)
{
  y.ht.i[i] = as.vector(X.new[i,]%*%beta.ht.i[i,])
}

dffit.i = rep(0,n.new)
for(i in 1:n.new)
{
  dffit.i[i] = (y.new.ht[i] - y.ht.i[i])/sqrt(S_i.sq[i]*hii[i])
}
which.max(dffit.i)

ind.2 = as.numeric(as.factor(abs(dffit.i) > 0.5 ))

## COVRATIO

covratio.i = rep(0,n.new)
for(i in 1:n.new)
{
  covratio.i[i] = ((S_i.sq[i]/MSres)^p)/(1-hii[i])
}

ind.3 = as.numeric(as.factor( covratio.i > 1.1 | covratio.i < 0.9 ))
which.max(covratio.i)

##plots

colors = c( "black", "blue" )
colors.1 = colors[ind.1]
colors.2 = colors[ind.2]
colors.3 = colors[ind.3]

par(mfrow=c(1,3))
plot(1:n.new , C_i, ylab = "Cook's distance", pch = 16,cex = 0.8, col = colors.1 )
plot(1:n.new , dffit.i, ylab = "DFFITS", pch = 16,cex = 0.8, col = colors.2 )
plot(1:n.new , covratio.i,ylab = "COVRATIO", pch = 16,cex = 0.8, col = colors.3 )


## Determining the influential observations

index = which((abs(dffit.i) > 0.5) & (covratio.i > 1.1 | covratio.i < 0.9))
suspects = dfbeta.ji[index,]
M = abs(suspects) > 0.1

inf <- rep(0,nrow(M))

for(i in 1:nrow(M)){
  inf[i] <- sum(M[i,])
}

influ_pts = index[inf > 6]

sample_points <- 1:n.new

par(mfrow = c(1,3))
plot(1:n.new, C_i , ylab = "Cooks distance" , pch = 16 , cex =0.8, ylim = c(0,0.3), col = ifelse(sample_points == 293 | sample_points == 295 | sample_points == 296 |sample_points == 297 |sample_points == 299|sample_points == 302 |sample_points == 332 ,"blue","black"))
plot(1:n.new, dffit.i,ylab = "DFFITS" , pch = 16 ,cex = 0.8, col = ifelse(sample_points == 293 | sample_points == 295 | sample_points == 296 |sample_points == 297 |sample_points == 299|sample_points == 302 |sample_points == 332 ,"blue","black"))
plot(1:n.new, covratio.i,ylab = "COVRATIO",pch = 16 ,cex = 0.8, col = ifelse(sample_points == 293 | sample_points == 295 | sample_points == 296 |sample_points == 297 |sample_points == 299|sample_points == 302 |sample_points == 332 ,"blue","black"))


y.new1 = y.new[-influ_pts  ]
X.new1 = X.new[-influ_pts ,]



# OLS estimate of y.test1 ( y.test1.ht )

beta.new1.ht = as.vector(solve(crossprod(X.new1))%*%crossprod(X.new1,y.new1))
y.test1.ht = as.vector(X.test%*%beta.new1.ht)

# new RMSE

RMSE1 = as.vector(sqrt(crossprod(y.test - y.test1.ht)/length(y.test)))


paste0("RMSE = ", RMSE)
paste0("RMSE_new = ", RMSE1)

## Data to be used in question 2

data.train1 <- cbind(X.new1,y.new1)
data.test <- cbind(X.test,y.test)

write.table(data.train1, file = "data.wo.influential.csv" , row.names = F, sep = "," )
write.table(data.test, file = "data.test.csv" , row.names = F, sep = "," )

##########****************Q2***************################################

data_train=read.csv("C:\\Users\\Sayan\\Downloads\\data.wo.influential.csv",header=T)
data_test=read.csv("C:\\Users\\Sayan\\Downloads\\data.test.csv",header=T)

X.train = as.matrix(data_train[, -ncol(data_train)])
y.train = as.vector(data_train[, ncol(data_train)])

X.test = as.matrix(data_test[, -ncol(data_test)])
y.test = as.vector(data_test[, ncol(data_test)])

#__________________________________________________________________________________________________________
###########################################################################################################
#__________________________________________________________________________________________________________

#---------------------------------------Dealing with the Curvature-----------------------------------------

#Running the OLS rgression and checking the RMSE
n = length(y.train)
p = length(X.train[1,])
H = X.train%*%solve(crossprod(X.train))%*%t(X.train)
In = diag(rep(1, n))
e = as.vector((In - H)%*%y.train)
beta = as.vector(solve(crossprod(X.train))%*%crossprod(X.train, y.train))

y.test_pred = X.test%*%beta
RMSE = sqrt(crossprod(y.test - y.test_pred)/length(y.test))
paste0("RMSE = ", RMSE)

#Checking the correlation between y and regressors
cor_init = cbind(cor(X.train[, -1], y.train))
cor_init

par(mfrow = c(2,5))
for(i in 2:p){
  plot(X.train[, i], e, col = 'blue', pch = 16, cex = 0.8, xlab = colnames(X.train)[i], ylab = "residuals")
}
par(mfrow = c(2,5))
for(i in 2:p){
  plot(X.train[, i], y.train, col = 'blue', pch = 16, cex = 0.8, xlab = colnames(X.train)[i], ylab = "MEDV")
}
suspects = c(2, 5, 6, 7, 10, 11)  # identifying suspects for non-linearity

#####################################  CPR & APR plots
par(mfrow = c(1,2))
cpr = e + X.train%*%diag(beta)

for(i in suspects){
  X1 = cbind(X.train, X.train[, i]^2)
  b = solve(crossprod(X1))%*%crossprod(X1, y.train)
  apr = y.train - X1[, c(-i, -(p+1))]%*%b[c(-i, -(p+1))]
  plot(X.train[, i], cpr[, i], col = 'blue', pch = 16, cex = 0.8, xlab = colnames(X.train)[i], ylab = "cpr")
  plot(X.train[, i], apr, col = 'blue', pch = 16, cex = 0.8, xlab = colnames(X.train)[i], ylab = "apr")
}

suspects = c(5, 7, 11) 

##################################### Partial residual plots
par(mfrow = c(1, 3))
for(i in suspects){
  X = X.train[, -i]
  Hi = X%*%solve(crossprod(X))%*%t(X)
  ey = (In - Hi)%*%y.train
  ei = (In - Hi)%*%X.train[, i]
  plot(ei, ey, col = 'blue', pch = 16, cex = 0.8, ylab = "Partial y", xlab = colnames(X.train)[i])
}

##################################### Transformations (determined empirically)
X = X.train
X[, 7] = exp(-3*X[, 7])
X[, 11] = exp(-1*X[, 11])
cor_final = cbind(cor(X[, -1], y.train))
cor_final

##################################### Calculating the RMSE
beta_trans = as.vector(solve(crossprod(X))%*%crossprod(X, y.train))
y.test_pred = X.test%*%beta_trans
RMSE = sqrt(crossprod(y.test - y.test_pred)/length(y.test))
paste0("RMSE = ", RMSE)

#Since RMSE increases on making the transformation, we drop them and continue with the original matrix.

############*****************Q3*********************#####################
x=read.csv("C:\\Users\\Sayan\\Downloads\\data.wo.influential.csv",header=T)
X=data.matrix(x) #Converting the data into matrix
p=ncol(X)  #Number of columns of X
y=X[,p]  #Vector of respones after removing influential points
X=X[,-p]  #deleting the column of the responses and creating the design matrix 
p=ncol(X)  #number of columns of X
n=nrow(X)  #number of rows of X
beta_hat=solve(crossprod(X),t(X)%*%y)#Estimate of beta_hat
y_hat=X%*%beta_hat     #Vector of fitted values
e=y-y_hat    #vector of residuals
H=X%*%solve(t(X)%*%X)%*%t(X)
h=diag(H)
RSS= sum((y-y_hat)^2)
r=e/sqrt((1-h)*RSS)
#........3(a).............

d=(e*e)/mean(e*e)   #calculating di's as said in part(b)

lab=c("CRIM","INDUS","NOX","RM","AGE","DIS","TAX","PTRATIO","B","LSTAT")
par(mfrow=c(1,1))
#.......Plotting Fitted values vs Residuals
plot(y_hat,e,xlab="Fitted Values",ylab="Residuals",col="green")

#.......Plotting Fitted Values v/s Studentized Residuals
plot(y_hat,r,xlab="Fitted Values",ylab="Studentized Residuals",col="red")


#Plotting the regressors vs residuals
par(mfrow=c(2,(p-1)/2))
for(i in 2:p)
{
  plot(X[,i],e,xlab=lab[i-1],ylab="Residuals",main="Scatter Plot",col="red")
}

#......Plotting the regressors vs Estimate of variances
par(mfrow=c(2,(p-1)/2))
for(i in 2:p)
{
  plot(X[,i],e*e,xlab=lab[i-1],ylab="sqrd Residuals",main="Scatter PLot",col="green")
}

#........Plotting the regressors vs di's
par(mfrow=c(2,(p-1)/2))
for(i in 2:p)
{
  plot(X[,i],d,xlab=lab[i-1],ylab="d",main="Scatter PLot",col="blue")
}

#...............3(b).............
Z=matrix(c(rep(1,times=n),X[,2],X[,5],X[,6],X[,7],X[,10],X[,11]),nrow=n) #Forming Z matrix after dropping some of the regressors

d_hat=Z%*%(solve(crossprod(Z),t(Z)%*%d)) #fitted values of di's
r=(cov(d,d_hat))/(sqrt(var(d))*sqrt(var(d_hat))) 
obs_chisq=n*r*r #Required Test Statistic
critical= qchisq(0.95,7) 
if(obs_chisq>critical) print("(Ho: model is homoscedastic) may be rejected") else
  print("(Ho: model is homascedastic) is accepted")

#The hypothesis of homoscedasticty is rejected

#..............3(c)..............

#Iterative Algorithm
count=0
beta_hatold=beta_hat
e_new=e
alpha_hat=solve(crossprod(Z))%*%t(Z)%*%log(e_new^2)
alpha_new=NA
repeat{
  beta_hatnew=beta_hatold
  alpha_new=alpha_hat
  sig=exp(Z%*%alpha_new)
  sigma=diag(as.vector(sig))
  beta_hatold=solve(t(X)%*%solve(sigma)%*%X)%*%t(X)%*%solve(sigma)%*%y
  e_new=y-X%*%beta_hatold
  alpha_hat=solve(crossprod(Z))%*%t(Z)%*%log(e_new^2)
  count=count+1
  if(((norm(beta_hatnew-beta_hatold,"f"))^2<=0.1 & (norm(alpha_new-alpha_hat,"f"))^2<=0.1)|count>200)
  {break
  }
}
beta_hatupdated=beta_hatold
#End of the Algorithm

#...........3(d)..................
z=read.csv("C:\\Users\\Sayan\\Downloads\\data.test.csv",header=T)
X_test=data.matrix(z) #COnverting the test data into matrix
p1=ncol(X_test)
n1=nrow(X_test)
y_test=X_test[,p1]
X_test=X_test[,-p1]
y_test_hat=X_test%*%beta_hatupdated  #fitted y values of the test data based on the updated estimate of beta_hat

RMSE_updated=norm(as.matrix(y_test-y_test_hat),"f")/sqrt(n1) #Updated RMSE
RMSE_updated


###******************Q4***********************###
data=read.csv("C:\\Users\\Sayan\\Downloads\\data.wo.influential.csv",header = T)
Z=data.matrix(data) #Converting the data into matrix
p=ncol(Z)  #Number of columns of X
y=Z[,p]  #Vector of respones after removing influential points
Z=Z[,-p]  #deleting the column of the responses and creating the design matrix
R=solve(t(Z)%*%Z)
P=Z%*%R%*%t(Z) #projection matrix
n=length(y)
#########################################################################
#Q(4.a)

yh=P%*%y #predicted values
e=y-yh #residuals

RSS=sum(e*e) #RSS
#Constructing RSS(i) vector
RSSi=rep(NULL,times=dim(Z)[1])
for(j in 1:n)
{
  RSSi[j]=RSS-((e[j]^2)/(1-P[j,j]))
}

#Constructing R-Student residual vector
ri=rep(NULL,times=n)
for(j in 1:n)
{
  ri[j]=e[j]/(sqrt(RSSi[j]*(1-P[j,j])/495))
}

#qqplot

##Here the population distribution of R-Student residual is t distribution with df n-p-1
##As n=506,p=10 so the distribution of R-student is t distribution with df 495
par(mfrow=c(1,1))
set.seed(1)
rs=rt(1000,385)
qqplot(rs,ri,xlab='Theoretical quantiles',ylab='R-Student residuals',main='Fig:4.1 - Q-Q plot of R-Student Residuals');abline(0,1,col='red',lty=2)

##The distribution is positively skewed, i.e. has longer right tail##

##############################################################################

#Q(4.b)

#defining transformation function
ulambda=function(u,l)
{
  if(l!=0) return ((u^l-1)/l)
  else return (log(u))
}

#RSS(lambda)
u=rep(NULL,times=n)
RSSl=function(l)
{
  for(i in 1:n)
    u[i]=ulambda(y[i],l)
  
  uh=P%*%u
  eu=u-uh
  return (sum(eu*eu))
}

#profile log-liklihood of lambda
boxcox=function(l)
{
  c=(-n/2)*log(2*pi)-(n/2) #constant part
  d=(l-1)*sum(log(y))+(-n/2)*log(RSSl(l)/n) #variable part
  return (c+d)
}

#plotting profile log-liklihood of lambda in the interval (-3,3)
lambdain=seq(-3,3,length=100)
lambda=c(NULL,length=n)
for(i in 1:length(lambdain))
{lambda[i]=boxcox(lambdain[i])}

plot(lambdain,lambda,type='l',lwd=2,xlab = "Value of lambda",ylab = "log-liklihood",main="Fig:4.2 - Profile log-liklihood of lambda")
#From the plot we can see that, at approximately lambda=0, profile log-liklihood of lambda is maximized
#So, the optimal value of lambd is 0 in Box-Cox transformation

################################################################################

#Q(4.c)

#As lambda=0, we transform y to log(y)
yt=(y^(0.25)-1)/(0.25)

yth=P%*%yt #predicted values of transformed y
et=yt-yth #residuals

RSSt=sum(et*et) #RSS of transformed variable

#Constructing RSS(i) vector of transformed variable
RSSti=rep(NULL,times=n)
for(j in 1:n)
{
  RSSti[j]=RSSt-((et[j]^2)/(1-P[j,j]))
}

#Constructing R-Student residual vector of transformed variable
rti=rep(NULL,times=n)
for(j in 1:n)
{
  rti[j]=et[j]/(sqrt(RSSti[j]*(1-P[j,j])/495))
}

#qqplot of the transformed variable#

##Here the population distribution of R-Student residual is t distribution with df n-p-1##
##As n=506,p=10 so the distribution of R-student is t distribution with df 495##

set.seed(1)
rs=rt(1000,385)
qqplot(rs,rti,xlab='Theoretical Quantiles',ylab='R-Student residuals',main='Fig:4.3 - Q-Q plot of R-Student residuals \n of the transformed variables');abline(0,1,col='red',lty=2)

##After applying Box-Cox transformation we can still observe that normality assumption on transformed variable is not entirely correct

############################################################################
############################################################################

