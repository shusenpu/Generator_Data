rm(list=ls())
## Add libraries
library(tidyverse)
library(ggfortify)
library(survival)
library(patchwork)
library(ggsurvfit)
library(ggsci)
library(AdequacyModel)
source(file = "functions.R")

x <- sort(c(19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, 32, 32, 32, 32, 33, 34, 35, 35, 35, 35, 36, 36, 37, 37, 39, 42, 44))

t2gwe_result <- goodness.fit(pdf=dt2gwe,cdf= pt2gwe,starts=c(0.11265398,5.01217839,0.02231709), 
                             data=x,method = "BFGS", domain = c(0,100),
                             lim_inf = c(0,0,0),lim_sup = c(100,100,100))

wge_result <- goodness.fit(pdf=wge_pdf,cdf= wge_cdf,starts=c(19.941702395,0.866902664,0.009036835), 
                           data=x,method = "BFGS", domain = c(0,100), 
                           lim_inf = c(0,0,0),lim_sup = c(100,100,100))

lg2_result <- goodness.fit(pdf=lg2_pdf,cdf= lg2_cdf,starts=c(2.737394,6.577782e-03,1.612461e+02,1.019259),
                           data=x,method = "BFGS", domain = c(0,100),
                           lim_inf = c(0,0,0,0),lim_sup = c(100,100,100,100))

t2g_result <- goodness.fit(pdf=t2g_pdf,cdf= t2g_cdf,starts=c(1.684370e+03,6.431474), data=x,
                           method = "BFGS", lim_inf = c(0,0),
                           lim_sup = c(100,100))

eg2_result <- goodness.fit(pdf=eg2_pdf,cdf= eg2_cdf,
                           starts=c(1.553229,4.803123,0.169095e+04), data=x,
                           method = "BFGS", domain = c(0,100), lim_inf = c(0,0,0),
                           lim_sup = c(100,100,100))

ewl_result <- goodness.fit(pdf=ewl_pdf,cdf=ewl_cdf,
                           starts=c(0.8,0.02,0.8,0.2), data=x,
                           method = "BFGS", domain = c(0,100), lim_inf = c(0,0,0,0),
                           lim_sup = c(100,100,100,100))

ll_result <- goodness.fit(pdf=dloglogis,cdf=ploglogis,
                          starts=c(0.8, 0.8), data=x, method = "BFGS", domain = c(0,1000000),
                          lim_inf = c(0,0), lim_sup = c(100,100))

# t2gwe_result$W;eg2_result$W;wge_result$W;lg2_result$W;t2g_result$W;ewl_result$W
# t2gwe_result$A;eg2_result$A;wge_result$A;lg2_result$A;t2g_result$A;ewl_result$A
# t2gwe_result$AIC;eg2_result$AIC;wge_result$AIC;lg2_result$AIC;t2g_result$AIC;ewl_result$AIC
# t2gwe_result$BIC;eg2_result$BIC;wge_result$BIC;lg2_result$BIC;t2g_result$BIC;ewl_result$BIC
# t2gwe_result$`CAIC `;eg2_result$`CAIC `;wge_result$`CAIC `;lg2_result$`CAIC `;t2g_result$`CAIC `;ewl_result$`CAIC `
# t2gwe_result$HQIC;eg2_result$HQIC;wge_result$HQIC;lg2_result$HQIC;t2g_result$HQIC;ewl_result$HQIC
# t2gwe_result$Value;eg2_result$Value;wge_result$Value;lg2_result$Value;t2g_result$Value;ewl_result$Value
# t2gwe_result$KS$p.value;eg2_result$KS$p.value;wge_result$KS$p.value;lg2_result$KS$p.value;t2g_result$KS$p.value;ewl_result$KS$p.value
# t2gwe_result$KS$statistic;eg2_result$KS$statistic;wge_result$KS$statistic;lg2_result$KS$statistic;t2g_result$KS$statistic;ewl_result$KS$statistic



# GT-II

## Sorted Data
df = data.frame(x=x,evt=1)
km = survfit(Surv(x,evt)~1,data=df)

## Survival data
y_surv= km$surv


############################################
## PLOT 1: Survival Function of the Model ##
############################################
t2gwe <- function(alpha,beta,gamma,x){
  1 - exp(-alpha*(exp(gamma*x)-1)^(-beta))
}


## survival model #############
y_main=t2gwe(x=x, alpha=t2gwe_result$mle[1], 
             beta=t2gwe_result$mle[2],
             gamma=t2gwe_result$mle[3])

## group everything
df$y_main = y_main
df.surv = data.frame(x=unique(x),y_surv=y_surv)

df.grp = left_join(df,df.surv)

gg = df.grp %>%
  ggplot(aes(x=x))+
  geom_line(aes(y=y_main,col=I("black")),size=1) +
  geom_step(aes(y=y_surv,col="red"),size=1)+
  scale_colour_manual(name=NULL, #legend name
                      values = c("black","red"), #specify and change colors of all lines
                      labels = c("T2GWE","Kaplan-Meier estimate"))+
  labs(y="Survival Function", x = "x")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(0.95,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)

ggsave("~/Desktop/paper/code/plots/survival_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution


##################
## PLOT 2: ECDF ##
##################
gg = df.grp %>%
  ggplot(aes(x=x))+
  geom_line(aes(y=1-y_main,col=I("black")),size=1) +
  geom_step(aes(y=1-y_surv,col="red"),size=1)+
  scale_colour_manual(name=NULL, #legend name
                      values = c("black","red"), #specify and change colors of all lines
                      labels = c("T2GWE","ECDF"))+
  
  labs(y="F(x)", x = "x")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(0.95,0.3), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
ggsave("~/Desktop/paper/code/plots/ecdf_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution


##################################
## PLOT 3: Scaled TTT Transform ##
##################################
xs=x
df = tibble(x=(1:length(xs))/length(xs), 
            f = (cumsum(xs) + (length(xs)-(1:length(x)))*xs)/ sum(xs))
gg = df %>%
  ggplot(aes(x=x,y=f))+
  geom_line(linewidth=1)+
  geom_point(col="green",size=2)+
  geom_abline(linetype=2)+
  labs(y="Scaled TTT-Transform")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(1,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
ggsave("~/Desktop/paper/code/plots/TTT_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution




##############################
## PLOT 4: HRF of the model ##
##############################

x=sort(x)
hrf.t2gwe = function(x,alpha,beta,gamma){
  alpha*beta*gamma*(exp(gamma*x)-1)^(-beta-1)*exp(gamma*x-alpha*(exp(gamma*x)-1)^(-beta)
  )/(1-exp(-alpha*(exp(gamma*x)-1)^(-beta)))
}

hrf.df=tibble(x=x , y = hrf.t2gwe(x,
                                  alpha =t2gwe_result$mle[1], 
                                  beta = t2gwe_result$mle[2], 
                                  gamma = t2gwe_result$mle[3]))

gg= hrf.df %>%
  ggplot(aes(x=x,y=y))+
  geom_line(linewidth=1,col="black")+
  labs(y="h(x) of T2GWE")+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        #plot.margin = margin(1, 1, 1, 1, "cm"), # t, r, b, l
        legend.position = c(1,1), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1,1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
ggsave("~/Desktop/paper/code/plots/hrf_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution



## PLOT 5: Expected Probability Plots ##
########################################
F_observed = ((1:length(x))-0.375)/(length(x)+0.25)
F.t2gwe = 1- t2gwe(x, 
                   alpha =t2gwe_result$mle[1], 
                   beta = t2gwe_result$mle[2], 
                   gamma = t2gwe_result$mle[3])
SS.t2gwe = sum((F.t2gwe - F_observed)^2)

# weibull generalized exponential
wge_cdf<-function(x, alpha, theta,gamma){
  1 - exp(-alpha*(exp(gamma*x)-1)^theta)
}

F.wge = wge_cdf(x,alpha=wge_result$mle[1],
                theta=wge_result$mle[2],
                gamma=wge_result$mle[3])
SS.wge = sum((F.wge - F_observed)^2)


# Lomax Gumbel Type -Two
lg2_cdf <- function(alpha, beta, theta, k, x){
  1 - beta^alpha*(beta - log(1-exp(-theta*x^(-k))))^(-alpha)
}

F.lg2 = lg2_cdf(x, alpha=lg2_result$mle[1], 
                beta=lg2_result$mle[2], 
                theta=lg2_result$mle[3], 
                k=lg2_result$mle[4])
SS.lg2 = sum((F.lg2 - F_observed)^2)


# weibull generalized exponential
t2g_cdf <- function(alpha, nu,x){
  exp(-alpha*x^(-nu))
}

F.t2g = t2g_cdf(x,alpha=t2g_result$mle[1], nu=t2g_result$mle[2])
SS.t2g = sum((F.t2g - F_observed)^2)



# Exponentiated Gumbel Type-2

eg2_cdf <- function(alpha,phi,theta,x){
  1 - (1-exp(-theta*x^(-phi)))^alpha
}

F.eg2 = eg2_cdf(x,alpha=eg2_result$mle[1], 
                phi=eg2_result$mle[2], 
                theta=eg2_result$mle[3])
SS.eg2  = sum((F.eg2 - F_observed)^2)   



### EWL
ewl_cdf <- function(alpha, beta, lambda, theta, x){
  (1-exp(-alpha*exp(lambda*beta*x)))^theta
}

F.ewl = ewl_cdf(x, alpha=ewl_result$mle[1],
                beta=ewl_result$mle[2],
                lambda=ewl_result$mle[3],
                theta=ewl_result$mle[4])
SS.ewl = sum((F.ewl - F_observed)^2)   



# log-logistic
ll_cdf<-function(x, alpha, beta){
  (x^beta)/(x^beta+alpha^beta)
}

F.ll = ll_cdf(x, alpha=ll_result$mle[1], beta=ll_result$mle[2])
SS.ll = sum((F.ll - F_observed)^2)





## group data
df.cdf = tibble(F_observed,F.t2gwe,F.wge,F.ewl,F.eg2,F.lg2,F.t2g, F.ll) %>%
  pivot_longer(cols=-F_observed, names_to = "dist", values_to = "y")

gg= df.cdf %>%
  ggplot(aes(x=F_observed, y=y, col=dist))+
  geom_line()+
  scale_colour_manual(name="", #legend name with a title "Parameters" if you want
                      labels = c(paste0("(SS=",round(SS.eg2,4),") EGT"),
                                 paste0("(SS=",round(SS.ewl,4),") EWL"),
                                 paste0("(SS=",round(SS.lg2,4),") LGT"),
                                 paste0("(SS=",round(SS.ll,4),") LL"),
                                 paste0("(SS=",round(SS.t2g,4),") T2G"),
                                 paste0("(SS=",round(SS.t2gwe,4),") T2GWE"),
                                 paste0("(SS=",round(SS.wge,4),") WGE")
                      ),
                      values = c("orange","magenta","green","cyan","purple","red","blue"))+
  geom_line(aes(x=F_observed,y=F_observed), color="black")+
  theme_bw()+
  labs(x="Observed Probability", y="Expected Probability")+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 18),
        legend.position = c(0.95,0.3), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
ggsave("~/Desktop/paper/code/plots/pp_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution





#######################
## PLOT 6: Histogram ##
#######################

domain = seq(0, 45, 0.05)
# domain = seq(min(x)*0.95, max(x)*1.05, 0.01)

####  type2 - gumble weibull exponential
dt2gwe <- function(alpha,beta,gamma,x){
  alpha*beta*gamma*(exp(gamma*x)-1)^(-beta-1)*exp(gamma*x-alpha*(exp(gamma*x)-1)^(-beta))
}

y.t2gwe = dt2gwe(domain, alpha=t2gwe_result$mle[1], 
                 beta = t2gwe_result$mle[2], 
                 gamma = t2gwe_result$mle[3])



#### weibull generalized exponential
wge_pdf<-function(alpha, theta, gamma, x){
  alpha*theta*gamma*exp(gamma*x)*(exp(gamma*x)-1)^(theta-1)*exp(-alpha*(exp(gamma*x)-1)^theta)
}

y.wge = wge_pdf(domain, alpha=wge_result$mle[1],
                theta=wge_result$mle[2],
                gamma=wge_result$mle[3])



#### Lomax Gumbel Type -Two
lg2_pdf <- function(alpha, beta, theta, k, x){
  alpha*beta^alpha*theta*k*x^(-k-1)*exp(-theta*x^(-k))*(beta-log(1-exp(-theta*x^(-k))))^(-alpha-1)/(1-exp(-theta*x^(-k)))
}

y.lg2 = lg2_pdf(domain, 
                alpha=lg2_result$mle[1], 
                beta=lg2_result$mle[2], 
                theta=lg2_result$mle[3], 
                k=lg2_result$mle[4])



### Exponentiated Gumbel Type-2
t2g_pdf <- function(alpha,nu,x){
  alpha*nu*x^(-nu-1)*exp(-alpha*x^(-nu))
}

y.t2g = t2g_pdf(domain,alpha=t2g_result$mle[1],
                nu=t2g_result$mle[2])



### Exponentiated Gumbel Type-2
eg2_pdf <- function(alpha,phi,theta,x){
  alpha*phi*theta*x^(-phi-1)*exp(-theta*x^(-phi))*(1-exp(-theta*x^(-phi)))^(alpha-1)
}

y.eg2 = eg2_pdf(domain,alpha = eg2_result$mle[1], 
                phi = eg2_result$mle[2], 
                theta = eg2_result$mle[3])


### EWL
ewl_pdf <- function(alpha, beta, lambda, theta, x){
  theta*(1-exp(-alpha*exp(lambda*beta*x)))^(theta-1)*(lambda*alpha*beta*exp(lambda*beta*x-alpha*exp(lambda*beta*x)))
}

y.ewl = ewl_pdf(domain,
                alpha=ewl_result$mle[1],
                beta=ewl_result$mle[2],
                lambda=ewl_result$mle[3],
                theta=ewl_result$mle[4])


### log-logistic
ll_pdf <- function(alpha, beta, x){
  ((beta/alpha)*(x/alpha)^(beta-1))/((1+(x/alpha)^beta)^2)
}

y.ll = ll_pdf(domain, alpha=ewl_result$mle[1], beta=ewl_result$mle[2])

max_y = max(y.t2gwe,y.wge,y.ewl,y.eg2,y.lg2, y.t2g, y.ll)*1.05

## group models
new_densities = data.frame(domain,y.t2gwe,y.wge,y.ewl,y.eg2,y.lg2, y.t2g, y.ll) %>%
  pivot_longer(cols = -domain,names_to = "Densities", values_to = "y")

new_densities$x=c(x,rep(NA,nrow(new_densities)-length(x)))

gg<-ggplot(new_densities) +   
  geom_histogram(aes(x=x,y=after_stat(density)),
                 bins=15, 
                 colour=1,fill="grey")+ # breaks=seq(0,5.6,0.7)
  geom_line(aes(x=domain,y=y,col=Densities),lwd=0.8) +
  scale_colour_manual(name="", #legend name with a title "Parameters" if you want
                      labels = c("EGT","EWL","LGT","LL","T2G","T2GWE","WGE"),
                      values = c("purple","magenta","green","cyan","blue","red","orange"))+
  scale_linetype_manual("",
                        values = c(1:6),
                        labels = c("EGT","EWL","LGT","LL","T2G","T2GWE","WGE"))+
  labs(y="Density", x = "x")+
  # xlim(0.1,4)+
  ylim(0, max_y)+
  theme_bw()+
  theme(axis.text = element_text(face="bold",size = 20,color="black"), #axis font 
        axis.title = element_text(size=20,face="bold"),
        axis.line.x.bottom = element_line(color = "black",size = 0.8), #control axis and its color  
        axis.line.y.left = element_line(color = "black", size = 0.8), #control axis
        # legend.text = element_text(face="bold",size = 10), #uncomment if you want to bold
        legend.text = element_text(size = 15),  #legend font
        legend.title = element_text(face="bold",size = 19),
        legend.text.align = 0,
        legend.position = c(0.95,1.025), #position of the legend
        legend.justification = c("right", "top"), # location of the legend within the specified location
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1), #change this will change the layout of the legend in the box created
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), #remove all grids, delete if you want grid mesh on
        legend.background = element_rect(fill = "transparent")) #make transparent background in legend
print(gg)
ggsave("~/Desktop/paper/code/plots/hist_depressive_new.eps", gg, width=18, height=18, units="cm", dpi=1080) #save with a specific dimension and resolution






