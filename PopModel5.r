###################################################################################
#					White-tailed deer and coyote population model				  #
###################################################################################



##########			          Set up parameter values 					 ##########

it <- 10000						#Number of iterations
yrs <- 15						#Number of years
coyote.removal.goal <- 0.0		#Set coyote removal goal to be used in the 		simulation. Proportion of each class removed is then random and sums to 1/2 the goal because half removed are male/female
doe.harvest <- 0.0				#Set doe harvest percent to be used for the simulation. Proportion of yearling and adult females harvested is then random and sums to the % harvest. 
buck.harvest <- 0.0				#Set buck harvest percent to be used for the simulaiton. Proportion of yearling and adult males harvested is then random and sums to the % harvest. 
f.fawn.harvest <- 0				#Proportion of female fawns harvested
m.fawn.harvest <- 0				#Proportion of male fawns harvested

initial.deer.density <- 6		#Intial deer density to determine starting deer population size
area <- 250						#Area of interest to determine starting population size
buck.doe.ratio <- 0.67			#Buck:doe ratio expressed as proportion, i.e., 0-1, to determine starting deer population
fawn.doe.ratio <- 0.4			#Fawn:doe ratio expressed as proportion to determine starting deer population
yearling.adult.ratio <- 1.25	#Yearling:adult ratio expressed as proportion to determine starting deer population
resident.density <- 0.92*0.5
transient.density <- 0.7*0.5




##########			          Deer population parameters				 ##########

#Determine inital population size for deer
total.pop <- initial.deer.density*area
females <- initial.deer.density*area/(buck.doe.ratio+1)
Nffi<-Nfmi <- 0.5*(fawn.doe.ratio*females/(0.5*fawn.doe.ratio+1))	#Initial number of fawn females and fawn males age 0.5
Nyfi <- yearling.adult.ratio*(females-Nffi)/(1+yearling.adult.ratio) #Intial number of yearling females age 1.5
Nymi <- yearling.adult.ratio*(total.pop-females-Nfmi)/(1+yearling.adult.ratio) #Intial number of yearling males age 1.5
Nafi <- females-Nyfi-Nffi		#Initial number of adult females age 2.5+
Nami <- total.pop-females-Nymi-Nfmi		#Initial number of adult males age 2.5+

#Empty matrices for population size
Nff = matrix(,it,yrs)			#Fawn females
Nfm = matrix(,it,yrs)			#Fawn males
Nyf = matrix(,it,yrs)			#Yearling females
Nym = matrix(,it,yrs)			#Yearling males
Naf = matrix(,it,yrs)			#Adult females
Nam = matrix(,it,yrs)			#Adult males
Nnf = matrix(,it,yrs)			#Newborn females, < 0.5
Nnm = matrix(,it,yrs)			#Newborn males, < 0.5
Nd = matrix(,it,yrs)			#Total deer population size 

prefd = matrix(,it,yrs)			#Fawn:doe ratios pre harvest
postfd = matrix(,it,yrs)		#Fawn:doe ratios post harvest
bd = matrix(,it,yrs)			#Buck:doe ratios
DeerDens = matrix(,it,yrs)		#Empty matrix for deer density
FawnSurv = matrix(,it,yrs)		#Empty matrix for total fawn survival 
YrlDoeSurv = matrix(,it,yrs)		#Total yearling doe survival
YrlBuckSurv = matrix(,it,yrs)		#Totaly yearling buck survival
AdDoeSurv = matrix(,it,yrs)		#Total adult doe survival
AdBuckSurv = matrix(,it,yrs)		#Total adult buck survival


###Survival###
#Empty matrices for survival
Sff = matrix(,it,yrs)			
Sfm = matrix(,it,yrs)			
Syf = matrix(,it,yrs)			
Sym = matrix(,it,yrs)			
Saf = matrix(,it,yrs)			
Sam = matrix(,it,yrs)			
So = matrix(,it,yrs)			#Survival of newborns till 6 months of age
Safmax = matrix(,it,yrs)		#Maximum survival of adult females			
Sammax = matrix(,it,yrs)		#Maximum survival of adult males

#Density-dependent parameters for adult survival
sdalpha = runif(1:it, 19000, 21000)
sdbeta = runif(1:it, 2, 4)

#Beta distribution for mean survival and the variance to draw from for each iteration
Sffa = rbeta(1:it, 24.8, 6.2)	#Mean	
SffVar = 0.001				#Variance
Sfma = rbeta(1:it, 24.8, 6.2)		
SfmVar = 0.001				
Syyfa = rbeta(1:it, 44.175, 2.325)		
SyfVar = 0.001				
Syyma = rbeta(1:it, 27.257, 0.843)		
SymVar = 0.001				
Safa = rbeta(1:it, 59.613, 4.487)	
SafVar = 0.001				
Sama = rbeta(1:it, 32.936, 2.864)	
SamVar = 0.001				
Soa = rbeta(1:it, 11.67475, 5.75025)	
SoVar = 0.005				

#This prevents NaNs from being produced when the mean value for survival is too close to 1
Syfa = matrix(,it)
Syma = matrix(,it)
for(i in 1:it){
  if(Syyfa[i] > 0.990) Syfa[i] = 0.990 else Syfa[i]=Syyfa[i]
  if(Syyma[i] > 0.990) Syma[i] = 0.990 else Syma[i]=Syyma[i]
}

#Empty matrices for Alpha and Beta shape parameters to create survival distributions for each iteration
aSff = matrix(,it,1)			
aSfm = matrix(,it,1)			
aSyf = matrix(,it,1)			
aSym = matrix(,it,1)			
aSaf = matrix(,it,1)			
aSam = matrix(,it,1)			
aSo = matrix(,it,1)			
bSff = matrix(,it,1)			
bSfm = matrix(,it,1)			
bSyf = matrix(,it,1)			
bSym = matrix(,it,1)			
bSaf = matrix(,it,1)			
bSam = matrix(,it,1)			
bSo = matrix(,it,1)			



###Harvest###
#Empty matrix for harvest numbers
H = matrix(,it,yrs)			#Total harvest array

#Empty matrices for the porportion of deer to be harvested each year
hff = matrix(,it,yrs)			
hfm = matrix(,it,yrs)			
hyf = matrix(,it,yrs)			
hym = matrix(,it,yrs)			
haf = matrix(,it,yrs)			
ham = matrix(,it,yrs)			



###Fecundity###
#Empty matrices for female fecundity 
Ff = matrix(,it,yrs)			#Fecundity of fawn females
Fy = matrix(,it,yrs)			#Fecundity of yearling females
Fa = matrix(,it,yrs)			#Fecundity of adult females

#Empty matrices for productivity
mf = matrix(,it,yrs)			#Number of newborns produced by fawns (bred as fawns)
my = matrix(,it,yrs)			#Number of newborns produced by yearlings (bred as yearlings)
ma = matrix(,it,yrs)			#Number of newborns produced by adults
mymax = matrix(,it,yrs)			#Max productivity of yearlings to calculate productivity in presence of density-dependence
mamax = matrix(,it,yrs)			#Max productivity of adults

mdalpha = runif(1:it, 12000, 16000)				#Half saturation constant for adult and yearling productivity density dependence 
mdbeta = runif(1:it, 4, 6)				#Shape parameter for adult and yearling productivity density dependence

ff = matrix(,it,yrs)			#fertility function values for yearlings and adults
fff = matrix(,it,yrs)			#fertility function values for fawns

#Proportion of fawns bred, density-dependent
pff = matrix(,it,yrs)			#Empty matrix for % fawns bred
pffmax = matrix(,it,yrs)		#Max % fawns bred
pf = matrix(,it,yrs)			#Empty matrix for % yearlings and adults bred

pfaa = rbeta(1:it, 52.076, 3.324)	
pfVar = 0.001		
pffa = rbeta(1:it, 8.9, 80.1)
pffVar = 0.002			

pfa = matrix(,it)
for(i in 1:it){
  if(pfaa[i] > 0.990) pfa[i] = 0.990 else pfa[i]=pfaa[i]
}



#Empty matrices for alpha and beta shape parameters to creat % females bred for each iteration
apf = matrix(,it,1)
apff = matrix(,it,1)
bpf = matrix(,it,1)
bpff = matrix(,it,1)

pffalpha = runif(1:it, 1000, 3000)				#Half saturation constant for % fawns bred
pffbeta = runif(1:it, 2, 4)				#shape parameter

#Log-normal distribution for mean productivity and the variance to draw from for each iteration
mfa = rlnorm(1:it,-0.0004997502,0.03161488)	#Mean	
mfVar = 0.001				#Variance
mya = rlnorm(1:it,0.2594144,0.07680965)		
myVar = 0.223				
maa = rlnorm(1:it,0.5862458,0.05551276)		
maVar = 0.278				

#Empty matrices for Alpha and Beta shape parameters to create productivity distributions for each iteration
amf = matrix(,it,1)			
amy = matrix(,it,1)			
ama = matrix(,it,1)			
bmf = matrix(,it,1)			
bmy = matrix(,it,1)			
bma = matrix(,it,1)			





##########			          Coyote population parameters				 ##########

###Population###
Ntyi =66 				#Initial number of transient yearling coyotes
Nryi =39				#Initial number of resident yearling coyotes
Ntai = 22				#Initial number of transient adult coyotes
Nrai = 76				#Initial number of resident adult coyotes

#Empty matrices for population size
Nty = matrix(,it,yrs)			#Transient yearling
Nry = matrix(,it,yrs)			#Resident yearling
Nta = matrix(,it,yrs)			#Transient adult
Nra = matrix(,it,yrs)			#Resident adults
Nc = matrix(,it,yrs)			#Total coyote population size 
ResAdSurv = matrix(,it,yrs)
ResYrlSurv = matrix(,it,yrs)
TrnAdSurv = matrix(,it,yrs)
TrnYrlSurv = matrix(,it,yrs)

TTr=77					#Number of territories

###Survival###
#Empty matrices for survival
Sty = matrix(,it,yrs)			
Sry = matrix(,it,yrs)
Sta = matrix(,it,yrs)
Sra = matrix(,it,yrs)
Sj = matrix(,it,yrs)			#Empty matrix for calculated juvenile survival
Sjmax = matrix(,it,yrs)			#Empty matrix for max juvenile survival from beta distribution used to calculate juvenile survival

#Beta distribution for survival and the variance to draw from for each iteration
Stya = rbeta(1:it, 8.8881, 13.9019)	#Mean	
StyVar = 0.005				#Variance 
Sryya = rbeta(1:it, 12.8736, 10.9664)		
SryVar = 0.005	
Staa = rbeta(1:it, 8.736, 9.464)		
StaVar = 0.005	
Sraaa = rbeta(1:it, 6.407211, 1.1913842)		
SraVar = 0.005	
Sja = rbeta(1:it, 8.327273, 12.49091)		
SjVar = 0.01

Sraa = matrix(,it)
Srya = matrix(,it)
for(i in 1:it){
  if(Sraaa[i] > 0.990) Sraa[i] = 0.990 else Sraa[i]=Sraaa[i]
  if(Sryya[i] > 0.990) Srya[i] = 0.990 else Srya[i]=Sryya[i]
}

#Empty matrices for Alpha and Beta shape parameters to create survival distributions for each iteration
aSty = matrix(,it,1)			
aSry = matrix(,it,1)
aSta = matrix(,it,1)			
aSra = matrix(,it,1)	
aSj = matrix(,it,1)		
bSty = matrix(,it,1)			
bSry = matrix(,it,1)			
bSta = matrix(,it,1)			
bSra = matrix(,it,1)
bSj = matrix(,it,1)

#Shape parameters for the effect of kill rate on juvenile survival
salpha=matrix(runif(1:it, 1, 4), it, 1)	#Half saturation constant for juvenile survival calculation
sbeta=matrix(runif(1:it, 1, 4), it, 1)	#Shape parameter for juvenile survival calculation
E = runif(1:it, 0.25, 0.75)		#Min kill rate before the the deficit function kicks in and decreases juvenile survival with decreasing kill rate
#salpha = 2.5
#sbeta = 2.5
#E = 0.50
L = matrix(,it,yrs)			#Empty matrix for the deficit between kill rate and the min required (E)



###Kill rate###
k=matrix(,it,yrs)			#Empty matrix for kill rate per coyote
kalpha = runif(1:it, 0.2, 0.5)		#Half saturation constant for coyote kill rate calculation
kbeta = runif(1:it, 3, 4)		#Shape parameter for kill rate calculation
#kbeta = 3.5
kmax = runif(1:it, 4, 7)		#Max kill rate per coyote
w = rbeta(1:it, 1.7625, 9.9875)		#Random number for transient kill rate adjustment


###Coyote removal###
R = matrix(,it,yrs)			#Empty matrix for number of coyotes removed

#Empty matrices for proportion of each stage removed
rty = matrix(,it,yrs)			
rry = matrix(,it,yrs)			
rta = matrix(,it,yrs)			
rra = matrix(,it,yrs)



###Fecundity###
#Empty matrices for fucundity
Fry = matrix(,it,yrs)			
Fra = matrix(,it,yrs)	

#Empty matrices for productivity					
mry = matrix(,it,yrs)			#Calculated productivity
mra = matrix(,it,yrs)
mrymax = matrix(,it,yrs)		#Max productivity from log-normal distribution used to calclute productivity			
mramax = matrix(,it,yrs)

#Proportion of yearlings and adults bred
pry = matrix(,it,yrs)			#Empty matrix for % of yearlings bred
prymax = matrix(,it,yrs)		#Empty matrix for max % of yearlings bred
pra = matrix(,it,yrs)			#Empty matrix for % of adults bred

apry = matrix(,it,1)			#Empty matrces for alpha and beta shape parameters
bpry = matrix(,it,1)
apra = matrix(,it,1)		
bpra = matrix(,it,1)

pryaa = rbeta(1:it, 1.756851, 2.426128)	#Distribution to drawn mean from for each iteration	for % yearlings bred
pryVar = 0.005				#Variance to use within each iteration to create distribution for time steps
praa = rbeta(1:it, 7.5652778, 4.07361111)	#Mean	for % adults bred
praVar = 0.005				#Variance

prya = matrix(,it)
for(i in 1:it){
  if(pryaa[i] < 0.010) prya[i] = 0.010 else prya[i]=pryaa[i]
}


#palpha = 4.5
#pbeta = 4

palpha =runif(1:it, 3, 6)		#Half saturation constant for % yearlings bred, used for sensitivity analysis but a constant was set for simulations
pbeta =runif(1:it, 3, 5)		#Shape parameter for % yearlings bred, used for sensitivity analysis but a constant was set for simulations

#Log-normal distribution for mean max productivity and the variance to draw from for each iteration
mrya = rlnorm(1:it,1.468855,0.1596827)	#Mean	
mryVar = 0.05				#Variance	
mraa = rlnorm(1:it,1.86592,0.1084659)		
mraVar = 0.05

#Empty matrices for Alpha and Beta shape parameters to create mean max productivity distributions for each iteration
amry = matrix(,it,1)			
amra = matrix(,it,1)						
bmry = matrix(,it,1)						
bmra = matrix(,it,1)

#Density-dependent shape parameters for productivity
malpha=7.5				#Half saturation constant for adult and yearling productivity
mbeta=3					#Shape parameter for adult and yearling productivity



###Transition rates from transient to resident and juvenile dispersal###
#Empty matrices
Dj = matrix(,it,yrs)			#Juvenile dispersal
Ty = matrix(,it,yrs)			#Transition from transient yearling to resident adult
Ta = matrix(,it,yrs)			#Transition from transient adult to resident adult

talpha = runif(1:it, 4, 6)				#Half saturation constant for transition rate
Dcalpha = 0.02				#Juvenile dispersal constant





##########			          Run simulations						 ##########


#####Iteration loop#####
for(i in 1:it){
  
  
  
  ###Deer Population###
  #Generate Alpha and Beta shape parameters for survival to create survival distributions for each iteration
  aSff[i] = Sffa[i]*((Sffa[i]*(1-Sffa[i])/SffVar)-1)		#Alpha shape parameters for each deer class
  aSfm[i] = Sfma[i]*((Sfma[i]*(1-Sfma[i])/SfmVar)-1)
  aSyf[i] = Syfa[i]*((Syfa[i]*(1-Syfa[i])/SyfVar)-1)
  aSym[i] = Syma[i]*((Syma[i]*(1-Syma[i])/SymVar)-1)
  aSaf[i] = Safa[i]*((Safa[i]*(1-Safa[i])/SafVar)-1)
  aSam[i] = Sama[i]*((Sama[i]*(1-Sama[i])/SamVar)-1)
  aSo[i] = Soa[i]*((Soa[i]*(1-Soa[i])/SoVar)-1)
  bSff[i] = (1-Sffa[i])*((Sffa[i]*(1-Sffa[i])/SffVar)-1)		#Beta shape parameters for each deer class	
  bSfm[i] = (1-Sfma[i])*((Sfma[i]*(1-Sfma[i])/SfmVar)-1)
  bSyf[i] = (1-Syfa[i])*((Syfa[i]*(1-Syfa[i])/SyfVar)-1)
  bSym[i] = (1-Syma[i])*((Syma[i]*(1-Syma[i])/SymVar)-1)
  bSaf[i] = (1-Safa[i])*((Safa[i]*(1-Safa[i])/SafVar)-1)
  bSam[i] = (1-Sama[i])*((Sama[i]*(1-Sama[i])/SamVar)-1)
  bSo[i] = (1-Soa[i])*((Soa[i]*(1-Soa[i])/SoVar)-1)
  
  #Generate Alpha and Beta shape parameters for productivity to create productivity distributions for each iteration
  amf[i] = log((mfa[i]^2)/(sqrt(mfVar+(mfa[i]^2))))		#Alpha shape parameters
  amy[i] = log((mya[i]^2)/(sqrt(myVar+(mya[i]^2))))
  ama[i] = log((maa[i]^2)/(sqrt(maVar+(maa[i]^2))))
  bmf[i] = sqrt(log(1+(mfVar/(mfa[i]^2))))			#Beta shape parameters
  bmy[i] = sqrt(log(1+(myVar/(mya[i]^2))))
  bma[i] = sqrt(log(1+(maVar/(maa[i]^2))))
  
  #Generate Alpha and Beta shape parameters for % bred distributions to be used in each iteration
  apf[i] = pfa[i]*((pfa[i]*(1-pfa[i])/pfVar)-1)
  apff[i] = pffa[i]*((pffa[i]*(1-pffa[i])/pffVar)-1)
  bpf[i] = (1-pfa[i])*((pfa[i]*(1-pfa[i])/pfVar)-1)
  bpff[i] = (1-pffa[i])*((pffa[i]*(1-pffa[i])/pffVar)-1)
  
  
  ###Coyote Population###
  #Generate Alpha and Beta shape parameters for survival distribtions to be used in each iteration
  aSry[i] = Srya[i]*((Srya[i]*(1-Srya[i])/SryVar)-1)		#Alpha shape parameters
  aSty[i] = Stya[i]*((Stya[i]*(1-Stya[i])/StyVar)-1)
  aSra[i] = Sraa[i]*((Sraa[i]*(1-Sraa[i])/SraVar)-1)
  aSta[i] = Staa[i]*((Staa[i]*(1-Staa[i])/StaVar)-1)
  aSj[i] = Sja[i]*((Sja[i]*(1-Sja[i])/SjVar)-1)
  bSry[i] = (1-Srya[i])*((Srya[i]*(1-Srya[i])/SryVar)-1)		#Beta shape parameters
  bSty[i] = (1-Stya[i])*((Stya[i]*(1-Stya[i])/StyVar)-1)
  bSra[i] = (1-Sraa[i])*((Sraa[i]*(1-Sraa[i])/SraVar)-1)
  bSta[i] = (1-Staa[i])*((Staa[i]*(1-Staa[i])/StaVar)-1)
  bSj[i] = (1-Sja[i])*((Sja[i]*(1-Sja[i])/SjVar)-1)
  
  #Generate Alpha and Beta shape parameters for productivity distributions to be used in each iteration
  amry[i] = log((mrya[i]^2)/(sqrt(mryVar+(mrya[i]^2))))
  amra[i] = log((mraa[i]^2)/(sqrt(mraVar+(mraa[i]^2))))
  bmry[i] = sqrt(log(1+(mryVar/(mrya[i]^2))))
  bmra[i] = sqrt(log(1+(mraVar/(mraa[i]^2))))
  
  #Generate Alpha and Beta shape parameters for % bred distributions to be used in each iteration
  apry[i] = prya[i]*((prya[i]*(1-prya[i])/pryVar)-1)
  apra[i] = praa[i]*((praa[i]*(1-praa[i])/praVar)-1)
  bpry[i] = (1-prya[i])*((prya[i]*(1-prya[i])/pryVar)-1)
  bpra[i] = (1-praa[i])*((praa[i]*(1-praa[i])/praVar)-1)
  
  
  #####Year loop#####
  for(j in 1:yrs){
    
    
    
    ###Initial population size, and determining survival and fecundity values for the year###
    
    ##Deer Population##
    #Specifcy population size for first 3 years
    Nff[i,1:3] = Nffi
    Nfm[i,1:3] = Nfmi
    Nyf[i,1:3] = Nyfi
    Nym[i,1:3] = Nymi
    Naf[i,1:3] = Nafi
    Nam[i,1:3] = Nami

    #Determine survival for each year by drawing from a beta distribution of survival from the iteration loop
    Sff[i,j] = rbeta(1,aSff[i], bSff[i])
    Sfm[i,j] = rbeta(1,aSfm[i], bSfm[i])
    Syf[i,j] = rbeta(1,aSyf[i], bSyf[i])
    Sym[i,j] = rbeta(1,aSym[i], bSym[i])
    Safmax[i,j] = rbeta(1,aSaf[i], bSaf[i])
    Sammax[i,j] = rbeta(1,aSam[i], bSam[i])
    #Syf[i,j] = rbeta(1,35.904, 1.496)
    #Sym[i,j] = rbeta(1,27.257, 0.843)
    #Safmax[i,j] = rbeta(1,59.613, 4.487)
    #Sammax[i,j] = rbeta(1,32.936, 2.864)
	  So[i,j] = rbeta(1,aSo[i], bSo[i])
    
    #Determine productivity for each year by drawing from a log-normal distribution of productivity from the iteration loop
    mf[i,j] = rlnorm(1,amf[i],bmf[i])
    mymax[i,j] = rlnorm(1,amy[i],bmy[i])
    mamax[i,j] = rlnorm(1,ama[i],bma[i])
    
    #Calculate productivity in presence of density-dependence
    if(j>1){ 
      my[i,j]=mymax[i,j]*(mdalpha[i]^mdbeta[i])/((mdalpha[i]^mdbeta[i])+((Nff[i,j-1]+Nfm[i,j-1]+Nyf[i,j-1]+Nym[i,j-1]+Naf[i,j-1]+Nam[i,j-1])^mdbeta[i]))
      ma[i,j]=mamax[i,j]*(mdalpha[i]^mdbeta[i])/((mdalpha[i]^mdbeta[i])+((Nff[i,j-1]+Nfm[i,j-1]+Nyf[i,j-1]+Nym[i,j-1]+Naf[i,j-1]+Nam[i,j-1])^mdbeta[i]))
      
      #Calculate % fawns, yearlings, and adults bred
      pf[i,j] = rbeta (1, apf[i], bpf[i])
      pffmax[i,j] = rbeta(1, apff[i], bpff[i])
      
      pff[i,j]=pffmax[i,j]*(pffalpha[i]^pffbeta[i])/((pffalpha[i]^pffbeta[i])+((Nff[i,j-1]+Nfm[i,j-1]+Nyf[i,j-1]+Nym[i,j-1]+Naf[i,j-1]+Nam[i,j-1])^pffbeta[i]))
      
      #Fertility function calculation for yearlings and adults
      if((Nam[i,j-1]+Nym[i,j-1]+Nfm[i,j-1])/(Nff[i,j-1]+Nyf[i,j-1]+Naf[i,j-1])<0.02) ff[i,j]=(pf[i,j]-.02)*((Nam[i,j-1]+Nym[i,j-1]+Nfm[i,j-1])/((Naf[i,j-1]+Nyf[i,j-1]+Nff[i,j-1])*.02)) else ff[i,j]=pf[i,j]
      
      #Fertility function calculation for fawns
      if((Nam[i,j-1]+Nym[i,j-1]+Nfm[i,j-1])/(Nff[i,j-1]+Nyf[i,j-1]+Naf[i,j-1])<0.02) fff[i,j]=(pff[i,j]-.02)*((Nam[i,j-1]+Nym[i,j-1]+Nfm[i,j-1])/((Naf[i,j-1]+Nyf[i,j-1]+Nff[i,j-1])*.02)) else fff[i,j]=pff[i,j]
      
      #Density-dependent adult survival
      Saf[i,j]=Safmax[i,j]*(sdalpha[i]^sdbeta[i])/((sdalpha[i]^sdbeta[i])+((Nff[i,j-1]+Nyf[i,j-1]+Naf[i,j-1])^sdbeta[i]))
      Sam[i,j]=Sammax[i,j]*(sdalpha[i]^sdbeta[i])/((sdalpha[i]^sdbeta[i])+((Nfm[i,j-1]+Nym[i,j-1]+Nam[i,j-1])^sdbeta[i]))
      
      
      ##Coyote Population##
      #Specify population size for first 3 years
      Nry[i,1:3] = Nryi
      Nty[i,1:3] = Ntyi
      Nra[i,1:3] = Nrai
      Nta[i,1:3] = Ntai
      
      #Determine survival for each year by drawing from a beta distribution of survival from the iteration loop
      Sry[i,j] = rbeta(1,aSry[i], bSry[i])
      Sty[i,j] = rbeta(1,aSty[i], bSty[i])
      Sra[i,j] = rbeta(1,aSra[i], bSra[i])
      Sta[i,j] = rbeta(1,aSta[i], bSta[i])
      Sjmax[i,j] = rbeta(1,aSj[i], bSj[i])
      
      #Determine max productivity for each year by drawing from a log-normal distribution of productivity from the iteration loop
      mrymax[i,j] = rlnorm(1,amry[i],bmry[i])
      mramax[i,j] = rlnorm(1,amra[i],bmra[i])
      
      #Density-dependent feedback for coyote population; calculate productivity based on number of resident coyotes
      mry[i,j] = mrymax[i,j]*(malpha^mbeta)/((malpha^mbeta)+(((Nry[i,j-1]+Nra[i,j-1])/TTr*2)^mbeta))
      mra[i,j] = mramax[i,j]*(malpha^mbeta)/((malpha^mbeta)+(((Nry[i,j-1]+Nra[i,j-1])/TTr*2)^mbeta))
      
      
      
      ###Deer Reproduction###
      #Calculte fecundity for deer
      if(0.5*fff[i,j]*mf[i,j]*sqrt(Sff[i,j])<0) Ff[i,j]=0 else Ff[i,j]=0.5*fff[i,j]*mf[i,j]*sqrt(Sff[i,j])	#Prevents negative numbers
      Fy[i,j]=0.5*ff[i,j]*my[i,j]*sqrt(Syf[i,j])
      Fa[i,j]=0.5*ff[i,j]*ma[i,j]*sqrt(Saf[i,j])
    }
    
    # Calculte number of newborn fawns
    if(j>2){
      Nnf[i,j]= Nff[i,j-1]*Ff[i,j-1]+Nyf[i,j-1]*Fy[i,j-1]+Naf[i,j-1]*Fa[i,j-1]
      Nnm[i,j]= Nff[i,j-1]*Ff[i,j-1]+Nyf[i,j-1]*Fy[i,j-1]+Naf[i,j-1]*Fa[i,j-1] 
      
      
      
      ###Predator Removal###
      rgoal=coyote.removal.goal			#Percent of the coyote population to remove set at the begining of the code
      
      #Random numbers to determine how many of each segment are removed which sum to the goal
      a=runif(1,0,rgoal)
      b=runif(1,0,rgoal-a)
      c=runif(1,0,rgoal-a-b)
      d=rgoal-a-b-c
      rrand=c(a,b,c,d)
      rrands=sample(rrand)
      
      #Proportion of each segment removed
      rry[i,j] = rrands[1]
      rty[i,j] = rrands[2]
      rra[i,j] = rrands[3]
      rta[i,j] = rrands[4]
      
      
      
      ###Predator kill rates and juvenile survival###
      #Calculate kill rate per coyote
      k[i,j] = kmax[i]*(kalpha[i]^kbeta[i])/((kalpha[i]^kbeta[i])+((((Nry[i,j-1]*2)+(Nra[i,j-1]*2)+(w[i]*((Nty[i,j-1]*2)+(Nta[i,j-1]*2))))/((Nnf[i,j]+Nnm[i,j])*So[i,j-1]))^kbeta[i]))
      
      #Juvenile coyote survival
      if(E[i]-k[i,j]<=0) L[i,j]=0 else L[i,j]=E[i]-k[i,j]					#Determines deficit, i.e. not enough fawns killed to result in max juvenile survival
      Sj[i,j] = Sjmax[i,j]*(1-((L[i,j]^sbeta[i])/((salpha[i]^sbeta[i])+(L[i,j]^sbeta[i]))))	#Calculte juvenile survival
      
      
      
      ###Calculate fecundity for coyotes###
      #Proportion of coyotes bred
      prymax[i,j] = rbeta(1,apry[i], bpry[i])	#Max % yearlings bred
      #% yearlings bred, prevents zeros
      if(prymax[i,j]*(palpha[i]^pbeta[i])/((palpha[i]^pbeta[i])+((((Nry[i,j-1])+(Nra[i,j-1]))/TTr*2)^pbeta[i]))<0) pry[i,j] = 0 else pry[i,j] = prymax[i,j]*(palpha[i]^pbeta[i])/((palpha[i]^pbeta[i])+((((Nry[i,j-1])+(Nra[i,j-1]))/TTr*2)^pbeta[i]))
      pra[i,j] = rbeta(1,apra[i], bpra[i])		#% adults bred		
      
      #Calculate fecundity
      Fry[i,j] = 0.5*pry[i,j]*mry[i,j]*Sj[i,j]*(Sry[i,j]^(1/4))
      Fra[i,j] = 0.5*pra[i,j]*mra[i,j]*Sj[i,j]*(Sra[i,j]^(1/4))
      
      
      
      ###Deer Harvest###
      
      dgoal=doe.harvest			#% females harvested set at begining of code
      
      #Random numbers to determine proportion of yearling and adult females harvested
      adoe=runif(1,0,dgoal)
      bdoe=dgoal-adoe
      drand=c(adoe,bdoe)
      drands=sample(drand)
      
      bgoal=buck.harvest			#% males harvested set at begining of code
      
      #Random numbers to determine proportion of yearling and adult males harvested
      abuck=runif(1,0,bgoal)
      bbuck=bgoal-abuck
      brand=c(abuck,bbuck)
      brands=sample(brand)
      
      #Proportion of each class harvested
      hff[i,j]= f.fawn.harvest
      hyf[i,j]= drands[1]
      haf[i,j]= drands[2]
      hfm[i,j]= m.fawn.harvest
      hym[i,j]= brands[1]
      ham[i,j]= brands[2]
      
      
      
      ###Transition and dispersal for coyotes###
      #Juvenile dispersal
      if(Dcalpha*((((Nry[i,j-1]-rry[i,j]*Nry[i,j-1])+(Nra[i,j-1]-rra[i,j]*Nra[i,j-1]))/TTr*2)^2)>1) Dj[i,j] = sum(rbinom(((Fry[i,j]*Nry[i,j-1]*(1-rry[i,j]))+(Fra[i,j]*Nra[i,j-1]*(1-rra[i,j]))),1,0.99))/((Fry[i,j]*Nry[i,j-1]*(1-rry[i,j]))+(Fra[i,j]*Nra[i,j-1]*(1-rra[i,j]))) else Dj[i,j] = sum(rbinom(((Fry[i,j]*Nry[i,j-1]*(1-rry[i,j]))+(Fra[i,j]*Nra[i,j-1]*(1-rra[i,j]))),1,Dcalpha*((((Nry[i,j-1]-rry[i,j]*Nry[i,j-1])+(Nra[i,j-1]-rra[i,j]*Nra[i,j-1]))/TTr*2)^2)))/((Fry[i,j]*Nry[i,j-1]*(1-rry[i,j]))+(Fra[i,j]*Nra[i,j-1]*(1-rra[i,j])))
      
      #Transition rate from transient to resident
      if(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2)<0) Ty[i,j] = 0 else if(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2)>1) Ty[i,j] = .5 else Ty[i,j] = .5*(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2))
      if(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2)<0) Ta[i,j] = 0 else if(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2)>1) Ta[i,j] = .5 else Ta[i,j] = .5*(((((rry[i,j]*Nry[i,j-1])+(rra[i,j]*Nra[i,j-1])))/((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1])))-((((Nry[i,j-1]+Nra[i,j-1])/TTr*2)/talpha[i])^2))
    }
    
    
    
    ###Population projection###
    #Calculte population size of each age/sex class for deer
    if(j>3){
      if(((Nnf[i,j]*So[i,j-1])-(0.5*k[i,j-1]*((Nry[i,j-1]-Nry[i,j-1]*rry[i,j-1])+(Nra[i,j-1]-Nra[i,j-1]*rra[i,j-1])+(w[i]*((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1]))))))*(1-hff[i,j-1])<0) Nff[i,j] = 0 else Nff[i,j]= ((Nnf[i,j]*So[i,j-1])-(0.5*k[i,j-1]*((Nry[i,j-1]-Nry[i,j-1]*rry[i,j-1])+(Nra[i,j-1]-Nra[i,j-1]*rra[i,j-1])+(w[i]*((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1]))))))*(1-hff[i,j-1])
      if(((Nnm[i,j]*So[i,j-1])-(0.5*k[i,j-1]*((Nry[i,j-1]-Nry[i,j-1]*rry[i,j-1])+(Nra[i,j-1]-Nra[i,j-1]*rra[i,j-1])+(w[i]*((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1]))))))*(1-hfm[i,j-1])<0) Nfm[i,j] = 0 else Nfm[i,j]= ((Nnm[i,j]*So[i,j-1])-(0.5*k[i,j-1]*((Nry[i,j-1]-Nry[i,j-1]*rry[i,j-1])+(Nra[i,j-1]-Nra[i,j-1]*rra[i,j-1])+(w[i]*((Nty[i,j-1]-rty[i,j]*Nty[i,j-1])+(Nta[i,j-1]-rta[i,j]*Nta[i,j-1]))))))*(1-hfm[i,j-1])
      Nyf[i,j]= Nff[i,j-1]*Sff[i,j-1]*(1-hyf[i,j-1])
      Nym[i,j]= Nfm[i,j-1]*Sfm[i,j-1]*(1-hym[i,j-1])
      Naf[i,j]= Nyf[i,j-1]*Syf[i,j-1]*(1-haf[i,j-1])+Naf[i,j-1]*Saf[i,j-1]*(1-haf[i,j-1]) 
      Nam[i,j]= Nym[i,j-1]*Sym[i,j-1]*(1-ham[i,j-1])+Nam[i,j-1]*Sam[i,j-1]*(1-ham[i,j-1])
      
      #Calculate population size of each age/stage class for coyotes
      Nty[i,j] = ((Fry[i,j-1]*Nry[i,j-1]*(1-rry[i,j-1]))+(Fra[i,j-1]*Nra[i,j-1]*(1-rra[i,j-1])))*Dj[i,j-1]
      Nry[i,j] = ((Fry[i,j-1]*Nry[i,j-1]*(1-rry[i,j-1]))+(Fra[i,j-1]*Nra[i,j-1]*(1-rra[i,j-1])))*(1-Dj[i,j-1])
      Nta[i,j] = Nty[i,j-1]*Sty[i,j-1]*(1-Ty[i,j-1])*(1-rty[i,j-1])+Nta[i,j-1]*Sta[i,j-1]*(1-Ta[i,j-1])*(1-rta[i,j-1])
      Nra[i,j] = Nty[i,j-1]*Sty[i,j-1]*Ty[i,j-1]*(1-rty[i,j-1])+Nry[i,j-1]*Sry[i,j-1]*(1-rry[i,j-1])+Nta[i,j-1]*Sta[i,j-1]*Ta[i,j-1]*(1-rta[i,j-1])+Nra[i,j-1]*Sra[i,j-1]*(1-rra[i,j-1])
    }
    
    #Calculte number of deer harvested and coyotes removed
    if(j>2){
      H[i,j]=Nff[i,j]*hff[i,j]+Nfm[i,j]*hfm[i,j]+Nyf[i,j]*hyf[i,j]+Nym[i,j]*hym[i,j]+Naf[i,j]*haf[i,j]+Nam[i,j]*ham[i,j]		#Number of deer harvested
      R[i,j] = Nry[i,j]*rry[i,j]+Nty[i,j]*rty[i,j]+Nra[i,j]*rra[i,j]+Nta[i,j]*rta[i,j]						#Number of coyotes removed
      
      #Calculate parameters of interest
      Nd[i,j]=Nff[i,j]+Nfm[i,j]+Nyf[i,j]+Nym[i,j]+Naf[i,j]+Nam[i,j]									#Total number of deer
      DeerDens[i,j] = Nd[i,j]/250													#Deer density
      Nc[i,j]=Nry[i,j]+Nty[i,j]+Nra[i,j]+Nta[i,j]											#Total number of coyotes
      prefd[i,j]=(Nff[i,j]+Nfm[i,j])/(Nyf[i,j]+(Nyf[i,j]*hyf[i,j-1])+Naf[i,j]+(Naf[i,j]*haf[i,j-1]))					#Pre-harvest fawn:doe ratio
      postfd[i,j]=(Nff[i,j]+Nfm[i,j])/(Nyf[i,j]+Naf[i,j])										#Post-harvest fawn:doe ratio
      bd[i,j]=(Nym[i,j]+Nam[i,j])/(Nyf[i,j]+Naf[i,j])											#Buck:doe ratio
      FawnSurv[i,j] = (Nff[i,j]+Nfm[i,j])/(Nnf[i,j-1]+Nnm[i,j-1])									#Total fawn survival to 6 months of age
      YrlDoeSurv[i,j] = (Nyf[i,j]*Syf[i,j]*(1-hyf[i,j]))/(Nff[i,j-1]*Sff[i,j-1]*(1-hyf[i,j-1]))					#Yearling doe survival with harvest
      YrlBuckSurv[i,j] = (Nym[i,j]*Sym[i,j]*(1-hym[i,j]))/(Nfm[i,j-1]*Sfm[i,j-1]*(1-hym[i,j-1]))					#Yearling buck survival with harvest
      AdDoeSurv[i,j] = (Naf[i,j]*Saf[i,j]*(1-haf[i,j]))/(Nyf[i,j-1]*Syf[i,j-1]*(1-haf[i,j-1])+Naf[i,j-1]*Saf[i,j-1]*(1-haf[i,j-1]))	#Adult doe survival with harvest
      AdBuckSurv[i,j] = (Nam[i,j]*Sam[i,j]*(1-ham[i,j]))/(Nym[i,j-1]*Sym[i,j-1]*(1-ham[i,j-1])+Nam[i,j-1]*Sam[i,j-1]*(1-ham[i,j-1]))	#Adult buck survival with harvest
      ResAdSurv[i,j] = (Nra[i,j]*Sra[i,j]*(1-rra[i,j]))/(Nty[i,j-1]*Sty[i,j-1]*Ty[i,j-1]*(1-rty[i,j-1])+Nry[i,j-1]*Sry[i,j-1]*(1-rry[i,j-1])+Nta[i,j-1]*Sta[i,j-1]*Ta[i,j-1]*(1-rta[i,j-1])+Nra[i,j-1]*Sra[i,j-1]*(1-rra[i,j-1]))
      ResYrlSurv[i,j] = (Nry[i,j]*Sry[i,j]*(1-rry[i,j]))/(((Fry[i,j-1]*Nry[i,j-1]*(1-rry[i,j-1]))+(Fra[i,j-1]*Nra[i,j-1]*(1-rra[i,j-1])))*(1-Dj[i,j-1]))
      TrnAdSurv[i,j] = (Nta[i,j]*Sta[i,j]*(1-rta[i,j])*(1-Ty[i,j]))/(Nty[i,j-1]*Sty[i,j-1]*(1-Ty[i,j-1])*(1-rty[i,j-1])+Nta[i,j-1]*Sta[i,j-1]*(1-Ta[i,j-1])*(1-rta[i,j-1]))
      TrnYrlSurv[i,j] = (Nty[i,j]*Sty[i,j]*(1-rty[i,j])*(Dj[i,j]))/(((Fry[i,j-1]*Nry[i,j-1]*(1-rry[i,j-1]))+(Fra[i,j-1]*Nra[i,j-1]*(1-rra[i,j-1])))*Dj[i,j-1])
    }
  }
  }






###Plot parameters of interest

#Plot deer density
se.deerdens <- apply(DeerDens, 2, function(x) sqrt(var(x)/length(x)))
plot(1:yrs, colMeans(DeerDens, na.rm=TRUE), pch=18, cex=1.5, xlab="Simulation years", ylab="Deer density", ylim=c(0,20))
arrows(1:yrs, colMeans(DeerDens, na.rm=TRUE), 1:yrs, colMeans(DeerDens, na.rm=TRUE)-(1.96*se.deerdens), length=.05, angle=90)
arrows(1:yrs, colMeans(DeerDens, na.rm=TRUE), 1:yrs, colMeans(DeerDens, na.rm=TRUE)+(1.96*se.deerdens), length=.05, angle=90)


#Plot number of coyotes
se.Nc <- apply(Nc, 2, function(x) sqrt(var(x)/length(x)))
plot(1:yrs, colMeans(Nc, na.rm=TRUE), pch=18, cex=1.5, xlab="Simulation years", ylab="Number of coyotes", ylim=c(0,250))
arrows(1:yrs, colMeans(Nc, na.rm=TRUE), 1:yrs, colMeans(Nc, na.rm=TRUE)-(1.96*se.Nc), length=.05, angle=90)
arrows(1:yrs, colMeans(Nc, na.rm=TRUE), 1:yrs, colMeans(Nc, na.rm=TRUE)+(1.96*se.Nc), length=.05, angle=90)


#Plot pre f:d ratios
se.prefd <- apply(prefd, 2, function(x) sqrt(var(x)/length(x)))
plot(1:yrs, colMeans(prefd, na.rm=TRUE), pch=18, cex=1.5, xlab="Simulation years", ylab="F:d ratios", ylim=c(0,1))
arrows(1:yrs, colMeans(prefd, na.rm=TRUE), 1:yrs, colMeans(prefd, na.rm=TRUE)-(1.96*se.prefd), length=.05, angle=90)
arrows(1:yrs, colMeans(prefd, na.rm=TRUE), 1:yrs, colMeans(prefd, na.rm=TRUE)+(1.96*se.prefd), length=.05, angle=90)





#Prepare output from for sensitivity analysis
SensOutput=data.frame(
  Sra=Sraa, 
  Sry=Srya, 
  Sta=Staa, 
  Sty=Stya, 
  Sj=Sja, 
  SjAplpha=salpha, 
  SjBeta=sbeta, 
  pra=praa, 
  pry=prya, 
  pCalpha=palpha, 
  pCbeta=pbeta, 
  mra=mraa, 
  mry=mrya,   
  E=E, 
  kmax=kmax, 
  kalpha=kalpha, 
  kbeta=kbeta, 
  kw=w, 
  I=talpha,
  Saf=Safa, 
  Sam=Sama, 
  Syf=Syfa, 
  Sym=Syma, 
  Sff=Sffa, 
  Sfm=Sfma, 
  So=Soa,
  sdalpha=sdalpha,
  sdbeta=sdbeta,
  pf=pfa, 
  pff=pffa, 
  mf=mfa, 
  my=mya, 
  ma=maa,
  mdalpha=mdalpha,
  mdbeta=mdbeta,
  pffalpha=pffalpha, 
  pffbeta=pffbeta,
  kill=apply(k, 1, mean, na.rm=TRUE),
  Prefd=apply(prefd, 1, mean, na.rm=TRUE), 
  DeerHarvest=apply(H, 1, mean, na.rm=TRUE), 
  DeerDensity=apply(DeerDens, 1, mean, na.rm=TRUE), 
  FawnSurvival=apply(FawnSurv, 1, mean, na.rm=TRUE),
  NCoyotes=apply(Nc, 1, mean, na.rm=TRUE), 
  YrlDoeSurv=apply(YrlDoeSurv, 1, mean, na.rm=TRUE), 
  YrlBuckSurv=apply(YrlBuckSurv, 1, mean, na.rm=TRUE),
  AdDoeSurv=apply(AdDoeSurv, 1, mean, na.rm=TRUE),
  AdBuckSurv=apply(AdBuckSurv, 1, mean, na.rm=TRUE),
  ResAdSurv=apply(ResAdSurv, 1, mean, na.rm=TRUE),
  ResYrlSurv=apply(ResYrlSurv, 1, mean, na.rm=TRUE),
  TrnAdSurv=apply(TrnAdSurv, 1, mean, na.rm=TRUE),
  TrnYrlSurv=apply(TrnYrlSurv, 1, mean, na.rm=TRUE))

#Write output file
write.csv(SensOutput, "SensitivityOutput.csv")


















