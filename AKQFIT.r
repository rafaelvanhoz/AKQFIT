#AKFIT R Edition
#Rafel Vanhoz Ribeiro

#INICIO
#escolhendo a biblioteca.

library(tidyverse)  #Tidtverse e uma coletanea de biblioteca para analise de dados
library(minpack.lm) #minpack e uma biblioteca para realizacao do levenberg - marquardt
library(matlib) #matlib e uma biblioteca que realiza calculos matriciais e inversoes de matrizes
library(scales) #Scales e uma bibioteca que configura, de modo facil, as notacoes cientificas
library(corrplot)
library(Hmisc)
library(broom)

options(width = 1000)
nls

#Criando as variaveis a aprtir do dado de entrada.
#O dado de entrada e um arquivo CVS delimitado por virgula e sem espaco.
#A variavel data pega o arquivo de entrada por meio da funcao read.csv(). O file.choose() é uma funcao que abre uma janela para selecionar o arquivo desejado, nao precisando escrever o nome do input direto no codigo. O header indica se o arquivo possui cabecalho.

data = read.csv(file.choose(), header = TRUE)
    Isotopo = as.character(data$Isotopo)
    Na = data$Na
    sNa = data$sNa
    fz = data$Fz
    sfz = data$sFz
    fa = data$Fa
    sfa = data$sFa
    D = data$D
    sD = data$sD
    C = data$C
    sC = data$sC
    S = data$S
    sS = data$sS
    w = data$W
    sw = data$sW
    NaCd = data$NaCd
    sNaCd = data$sNaCd
    fzCd = data$FZCd
    sfzCd = data$sFZCd
    faCd = data$FaCd
    sfaCd = data$sFaCd
    DCd = data$DCd
    sDCd = data$sDCd
    CCd = data$CCd
    sCCd = data$sCCd
    SCd = data$SCd
    sSCd = data$sSCd
    wCd = data$WCd
    swCd = data$sWCd
    Egama =  data$Egama
    IN =  data$IN
    Efic =  data$Efic
    sEfic =  data$sEfic
    Eres =  data$Eres
    sEres =  data$sEres
    k0 =  data$k0
    sk0 =  data$sk0
    FCd =  data$FCd
    sFCd =  data$sFCd
    Gepi =  data$Gepi
    sGepi =  data$sGepi
    Gth =  data$Gth
    sGth =  data$sGth
    Alfa = data$Alfa
    sAlfa = data$sAlfa
    FI = data$FI
    sFI = data$sFI
    idouro = data$idouro

#O dataau e identico ao data, porem e para os dados do ouro.
    
dataau = read.csv(file.choose(), header = TRUE)
    Isotopoau = dataau$Isotopoau
    Naau = dataau$Naau
    sNaau = dataau$sNaau
    fzau = dataau$Fzau
    sfzau = dataau$sFzau
    faau = dataau$Faau
    sfaau = dataau$sFaau
    Dau = dataau$Dau
    sDau = dataau$sDau
    Cau = dataau$Cau
    sCau = dataau$sCau
    Sau = dataau$Sau
    sSau = dataau$sSau
    wau = dataau$Wau
    swau = dataau$sWau
    NaCdau = dataau$NaCdau
    sNaCdau = dataau$sNaCdau
    fzCdau = dataau$FzCdau
    sfzCdau = dataau$sFzCdau
    faCdau = dataau$FaCdau
    sfaCdau = dataau$sFaCdau
    DCdau = dataau$DCdau
    sDCdau = dataau$sDCdau
    CCdau = dataau$CCdau
    sCCdau = dataau$sCCdau
    SCdau = dataau$SCdau
    sSCdau = dataau$sSCdau
    wCdau = dataau$WCdau
    swCdau = dataau$sWCdau
    Egamaau =  dataau$Egamaau
    INau =  dataau$INau
    Eficau =  dataau$Eficau
    sEficau =  dataau$sEficau
    Eresau =  dataau$Eresau
    sEresau =  dataau$sEresau
    k0au =  dataau$k0au
    sk0au =  dataau$sk0au
    FCdau =  dataau$Fcdau
    sFCdau =  dataau$sFcdau
    Gepiau =  dataau$Gepiau
    sGepiau =  dataau$sGepiau
    Gthau =  dataau$Gthau
    sGthau =  dataau$sGthau
    Q0au = dataau$Q0au
    sQ0au = dataau$sQ0au
    anoau = dataau$anoau

#Calcula o Asp, AspCd, Aspau, AspCdau, RCd e RCdau
Asp = (Na*fz*fa)/(D*C*S*w)
AspCd = (NaCd*fzCd*faCd)/(DCd*CCd*SCd*wCd)
Aspau = (Naau*fzau*faau)/(Dau*Cau*Sau*wau)
AspCdau = (NaCdau*fzCdau*faCdau)/(DCdau*CCdau*SCdau*wCdau)
RCd = Asp/AspCd
RCdau = Aspau/AspCdau

#calculando os erros pela propagacao de erros
Aspexp = expression((Na*fz*fa)/(D*C*S*w))
    dAspNa = eval(D(Aspexp, 'Na'))
    dAspfz = eval(D(Aspexp, 'fz'))
    dAspfa = eval(D(Aspexp, 'fa'))
    dAspD = eval(D(Aspexp, 'D'))
    dAspC = eval(D(Aspexp, 'C'))
    dAspS = eval(D(Aspexp, 'S'))
    dAspw = eval(D(Aspexp, 'w'))

AspCdexp = expression((NaCd*fzCd*faCd)/(DCd*CCd*SCd*wCd))
    dAspNaCd = eval(D(AspCdexp, 'NaCd'))
    dAspfzCd = eval(D(AspCdexp, 'fzCd'))
    dAspfaCd = eval(D(AspCdexp, 'faCd'))
    dAspDCd = eval(D(AspCdexp, 'DCd'))
    dAspCCd = eval(D(AspCdexp, 'CCd'))
    dAspSCd = eval(D(AspCdexp, 'SCd'))
    dAspwCd = eval(D(AspCdexp, 'wCd'))

Aspexpau = expression((Naau*fzau*faau)/(Dau*Cau*Sau*wau))
    dAspNaau = eval(D(Aspexpau, 'Naau'))
    dAspfzau = eval(D(Aspexpau, 'fzau'))
    dAspfaau = eval(D(Aspexpau, 'faau'))
    dAspDau = eval(D(Aspexpau, 'Dau'))
    dAspCau = eval(D(Aspexpau, 'Cau'))
    dAspSau = eval(D(Aspexpau, 'Sau'))
    dAspwau = eval(D(Aspexpau, 'wau'))

AspCdexpau = expression((NaCdau*fzCdau*faCdau)/(DCdau*CCdau*SCdau*wCdau))
    dAspNaCdau = eval(D(AspCdexpau, 'NaCdau'))
    dAspfzCdau = eval(D(AspCdexpau, 'fzCdau'))
    dAspfaCdau = eval(D(AspCdexpau, 'faCdau'))
    dAspDCdau = eval(D(AspCdexpau, 'DCdau'))
    dAspCCdau = eval(D(AspCdexpau, 'CCdau'))
    dAspSCdau = eval(D(AspCdexpau, 'SCdau'))
    dAspwCdau = eval(D(AspCdexpau, 'wCdau'))

sAsp = sqrt(dAspNa^2*sNa^2 + dAspfz^2*sfz^2 + dAspfa^2*sfa^2 + dAspD^2*sD^2 + dAspC^2*sC^2 + dAspS^2*sS^2 + dAspw^2*sw^2)
sAspCd = sqrt(dAspNaCd^2*sNaCd^2 + dAspfzCd^2*sfzCd^2 + dAspfaCd^2*sfaCd^2 + dAspDCd^2*sDCd^2 + dAspCCd^2*sCCd^2 + dAspSCd^2*sSCd^2 + dAspwCd^2*swCd^2)
sAspau = sqrt(dAspNaau^2*sNaau^2 + dAspfzau^2*sfzau^2 + dAspfaau^2*sfaau^2 + dAspDau^2*sDau^2 + dAspCau^2*sCau^2 + dAspSau^2*sSau^2 + dAspwau^2*swau^2)
sAspCdau = sqrt(dAspNaCdau^2*sNaCdau^2 + dAspfzCdau^2*sfzCdau^2 + dAspfaCdau^2*sfaCdau^2 + dAspDCdau^2*sDCdau^2 + dAspCCdau^2*sCCdau^2 + dAspSCdau^2*sSCdau^2 + dAspwCdau^2*swCdau^2)

#teste para verificar se a incerteza esta batendo com o outro metodo de calculo:
#sAsp2 = sqrt(((sNa/Na)^2 + (sfz/fz)^2 + (sfa/fa)^2 + (sD/D)^2 + (sC/C)^2 + (sS/S)^2 + (sw/w)^2)) * Asp 

#Define uma variavel que conta o total de linhas do arquivo de dados
#Essa variavel sera utilizada para determinar o n dos loops, for e if.
N = nrow(data)

#Determina os valores de K0 e Q0 usando as mesmas expressoes utilizadas no COVAR
#Estes valores serao utilizados para os valores iniciais dos porametros do levenberg-marquardt
Aa = 25.538
Alfa = Alfa[1]
Q0alfac = ((Q0au - 0.429)/(Eresau**Alfa)+0.429/((2*Alfa+1)*(0.55**Alfa)))
k0b=0
Q0Ialfa=0
Q0b = 0
for (i in 1:N) {
  k0b[i] = (Asp[i]-AspCd[i]/FCd[i])/(Aspau[idouro[i]]-AspCdau[idouro[i]]/FCdau[idouro[i]])*(Gthau[idouro[i]]*Eficau[idouro[i]])/(Gth[i]*Efic[i])
  Q0Ialfa[i] = (( FCdau[idouro[i]]*RCdau[idouro[i]] - 1 )/(FCd[i]*RCd[i] - 1 )* (Gth[i]*Gepiau[idouro[i]] /( Gthau[idouro[i]]*Gepi[i]))*Q0alfac[idouro[i]])
  Q0b[i] = ( Q0Ialfa[i]-0.429/((2*Alfa+1)*(0.55**Alfa)))*Eres[i]**Alfa + 0.429
}

#funcao para fazer matriz de q0
monta_matriz_por_isotopo <- function(m) {
  d = diag(m)
  j = 2
  for (i in 2:nrow(d)) {
    if (isTRUE(Isotopo[i-1] == Isotopo[i])) {
      j = j-1
    }
    d[i,j] = d[i,i]
    if (isFALSE(i == j)) {
      d[i,i] = 0
    }
    j = j+1
  }
  d = d[,1:j-1]
  result <- d
}



Q0b = monta_matriz_por_isotopo(Q0b)

nQ0 = ncol(Q0b)
nrQ0 = nrow(Q0b)
for(i in 1:nrQ0) {
  for(j in 1:nQ0) {
    if(Q0b[i,j]==0) {
      Q0b[i,j] = NA
    }
  }
}
Q0b = colMeans(Q0b, na.rm = TRUE)

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#Calcula o Y Superior da matrix Y
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

#Ys = expression(AspCd / Efic / FCd / Gepi)
Ys = expression(log(AspCd / Efic / FCd / Gepi))
Ysuperior = eval(Ys)

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#Calcula o Y Medio da Matrix Y
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

Ym = expression(log(((FCdau * (Aspau / AspCdau) - 1) * (Gth * Gepiau)) / ((FCd * (Asp / AspCd) - 1) * Gthau * Gepi)))
#Ym = expression(((FCdau * (Aspau / AspCdau) - 1) * (Gth * Gepiau)) / ((FCd * (Asp / AspCd) - 1) * Gthau * Gepi))
Ymedio = 0
for(i in 1:N){
   Ymedio[i] = log(((FCdau[idouro[i]] * (Aspau[idouro[i]] / AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) / ((FCd[i] * (Asp[i] / AspCd[i]) - 1) * (Gthau[idouro[i]] * Gepi[i])))
  #Ymedio[i] = ((FCdau[idouro[i]] * (Aspau[idouro[i]] / AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) / ((FCd[i] * (Asp[i] / AspCd[i]) - 1) * (Gthau[idouro[i]] * Gepi[i]))
}

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#Calcula o Y Inferior da matrix Y
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

Yi = expression(log(((Asp - AspCd / FCd) / (Aspau - AspCdau / FCdau)) * (Gthau / Gth * Eficau / Efic)))
Yinferior = 0
for(i in 1:N){
  Yinferior[i] = log(((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]))
  #Yinferior[i] = ((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i])
}


Y = c(Ysuperior, Ymedio, Yinferior)

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#derivadas parciais de Ysuperior, Ymedio e Yinferior
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

#Ysuperior
dYsAspCd = eval(D(Ys,'AspCd'))
dYsEfic = eval(D(Ys,'Efic'))
dYsFCd = eval(D(Ys,'FCd'))
dYsGepi = eval(D(Ys,'Gepi'))


#Ymedio
dYmAsp = 0
dYmAspCd = 0
dYmFCd = 0
dYmGth = 0
dYmGepi = 0
#dYmQ0au = 0
dYmAspau = 0
dYmAspCdau = 0
dYmFCdau = 0
dYmGthau = 0
dYmGepiau = 0
for (i in 1 : N) {
  dYmAsp[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * (FCd[i] * (1/AspCd[i]) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])))
  dYmAspCd[i] = ((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * (FCd[i] * (Asp[i]/AspCd[i]^2) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  dYmFCd[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((Asp[i]/AspCd[i]) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])))
  dYmGth[i] = (FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * Gepiau[idouro[i]]/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  dYmGepi[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])))
  dYmAspau[i] = FCdau[idouro[i]] * (1/AspCdau[idouro[i]]) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  dYmAspCdau[i] = -(FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]^2) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])))
  dYmFCdau[i] = (Aspau[idouro[i]]/AspCdau[idouro[i]]) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  dYmGthau[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])))
  dYmGepiau[i] = (FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * Gth[i]/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])/(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]]))/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  #dYmAsp[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * (FCd[i] * (1/AspCd[i]) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2)
  #dYmAspCd[i] = ((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * (FCd[i] * (Asp[i]/AspCd[i]^2) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2
  #dYmFCd[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((Asp[i]/AspCd[i]) * Gthau[idouro[i]] * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2)
  #dYmGth[i] = (FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * Gepiau[idouro[i]]/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])
  #dYmFCdau[i] = (Aspau[idouro[i]]/AspCdau[idouro[i]]) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])
  #dYmGepi[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2)
  #dYmAspau[i] = FCdau[idouro[i]] * (1/AspCdau[idouro[i]]) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])
  #dYmAspCdau[i] = -(FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]^2) * (Gth[i] * Gepiau[idouro[i]])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i]))
  #dYmGthau[i] = -(((FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * (Gth[i] * Gepiau[idouro[i]])) * ((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gepi[i])/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])^2)
  #dYmGepiau[i] = (FCdau[idouro[i]] * (Aspau[idouro[i]]/AspCdau[idouro[i]]) - 1) * Gth[i]/((FCd[i] * (Asp[i]/AspCd[i]) - 1) * Gthau[idouro[i]] * Gepi[i])
}


#Yinferior
dYiAsp = 0
dYiAspCd = 0
dYiFCd = 0
dYiGth = 0
dYiEfic = 0
dYiAspau = 0
dYiAspCdau = 0
dYiFCdau = 0
dYiGthau = 0
dYiEficau = 0
for (i in 1 : N) {
  dYiAsp[i] = 1/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])/(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]))
  dYiAspCd[i] = -(1/FCd[i]/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])/(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])))
  dYiFCd[i] = AspCd[i]/FCd[i]^2/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])/(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]))
  dYiGth[i] = -(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i]^2 * Eficau[idouro[i]]/Efic[i])/(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])))
  dYiEfic[i] = -(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]^2)/(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])))
  dYiAspau[i] = - ((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])^2 * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]) / (((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] /Gth[i] * Eficau[idouro[i]] / Efic[i])))
  dYiAspCdau[i] = (Asp[i] - AspCd[i] / FCd[i]) * (1 / FCdau[idouro[i]]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])^2 * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]) / (((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] /FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]))
  dYiFCdau[i] = - ((Asp[i] - AspCd[i] / FCd[i]) * (AspCdau[idouro[i]] / FCdau[idouro[i]]^2) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])^2 * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]) / (((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i])))
  dYiGthau[i] = ((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (1 / Gth[i] * Eficau[idouro[i]] / Efic[i]) / (((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]))
  dYiEficau[i] = ((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] / Efic[i]) / (((Asp[i] - AspCd[i] / FCd[i]) / (Aspau[idouro[i]] - AspCdau[idouro[i]] / FCdau[idouro[i]])) * (Gthau[idouro[i]] / Gth[i] * Eficau[idouro[i]] / Efic[i]))
  #dYiAsp[i] = 1/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])
  #dYiAspCd[i] = -(1/FCd[i]/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]))
  #dYiFCd[i] = AspCd[i]/FCd[i]^2/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]]) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])
  #dYiGth[i] = -(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i]^2 * Eficau[idouro[i]]/Efic[i]))
  #dYiEfic[i] = -(((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]^2))
  #dYiAspau[i] = -((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])^2 * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]))
  #dYiAspCdau[i] = (Asp[i] - AspCd[i]/FCd[i]) * (1/FCdau[idouro[i]])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])^2 * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i])
  #dYiFCdau[i] = -((Asp[i] - AspCd[i]/FCd[i]) * (AspCdau[idouro[i]]/FCdau[idouro[i]]^2)/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])^2 * (Gthau[idouro[i]]/Gth[i] * Eficau[idouro[i]]/Efic[i]))
  #dYiGthau[i] = ((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (1/Gth[i] * Eficau[idouro[i]]/Efic[i])
  #dYiEficau[i] = ((Asp[i] - AspCd[i]/FCd[i])/(Aspau[idouro[i]] - AspCdau[idouro[i]]/FCdau[idouro[i]])) * (Gthau[idouro[i]]/Gth[i]/Efic[i])
}

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#Calculo da Variancia
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

vY = matrix(0, nrow = dim(matrix(Y))[1], ncol = dim(matrix(Y))[1])

#covariancias entre YSuperior e YSuperior - Parte Diagonal
for(i in 1:N) {
  for(j in i:N) {
    if(isTRUE(Isotopo[i] == Isotopo[j]) & isTRUE(Egama[i] == Egama[j]) & isTRUE(idouro[i] == idouro[j])) {
      vY[i,j] = dYsAspCd[i] * dYsAspCd[j] * sAspCd[i]^2 + 
                dYsFCd[i] * dYsFCd[j] * sFCd[i]^2 +
                #dYsEfic[i] * dYsEfic[j] * sEfic[i]^2 +
                dYsGepi[i] * dYsGepi[j] * sGepi[i]^2
      vY[j,i] = vY[i,j]
    }
  }
}


#covariancias entre YSuperior e YSuperior - Parte nao diagonal, mesmo isotopo, gama diferente
for(i in 1:N) {
  for(j in i + 1:N) {
    if(isTRUE(Isotopo[i] == Isotopo[j]) & isTRUE(Egama[i] != Egama[j]) & isTRUE(idouro[i] == idouro[j])) {
      vY[i,j] = dYsFCd[i] * dYsFCd[j] * sFCd[i]^2 + 
                #dYsEfic[i] * dYsEfic[j] * sEfic[i]^2 +
                dYsGepi[i] * dYsGepi[j] * sGepi[i]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YMedio e YMedio - Parte diagonal
for(i in (N + 1):(N * 2)) {
  for(j in (i):(N * 2)) {
    if(isTRUE(Isotopo[i - N] == Isotopo[j - N]) & isTRUE(Egama[i - N] == Egama[j - N]) & isTRUE(idouro[i - N] == idouro[j - N])) {
      k = i - N
      l = j - N
      vY[i,j] = dYmAsp[k] * dYmAsp[l] * sAsp[k]^2 +
                dYmAspCd[k] * dYmAspCd[l] * sAspCd[k]^2 + 
                dYmFCd[k] * dYmFCd[l] * sFCd[k]^2 + 
                dYmGth[k] * dYmGth[l] * sGth[k]^2 + 
                dYmGepi[k] * dYmGepi[l] * sGepi[k]^2 +
                dYmAspau[idouro[k]] * dYmAspau[idouro[l]] * sAspau[idouro[k]]^2 +
                dYmAspCdau[idouro[k]] * dYmAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYmFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYmGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYmGepiau[idouro[k]] * dYmGepiau[idouro[l]] * sGepiau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YMedio e YMedio - Parte nao diagonal, mesmo isotopo, gama diferente
for(i in (N + 1):(N * 2)) {
  for(j in (i + 1):(N * 2)) {
    if(isTRUE(Isotopo[i - N] == Isotopo[j - N]) & isTRUE(Egama[i - N] != Egama[j - N]) & isTRUE(idouro[i - N] == idouro[j - N])) {
      k = i - N
      l = j - N
      vY[i,j] = dYmFCd[k] * dYmFCd[l] * sFCd[k]^2 + 
                dYmGth[k] * dYmGth[l] * sGth[k]^2 + 
                dYmGepi[k] * dYmGepi[l] * sGepi[k]^2 + 
                dYmAspau[idouro[k]] * dYmAspau[idouro[l]] * sAspau[idouro[k]]^2 +
                dYmAspCdau[idouro[k]] * dYmAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYmFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYmGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYmGepiau[idouro[k]] * dYmGepiau[idouro[l]] * sGepiau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}



#covariancias entre YMedio e YMedio - Parte nao diagonal, com Isotopo diferente
for(i in (N + 1):(N * 2)) {
  for(j in (i + 1):(N * 2)) {
    if(isTRUE(Isotopo[i - N] != Isotopo[j - N]) & isTRUE(Egama[i - N] != Egama[j - N]) & isTRUE(idouro[i - N] == idouro[j - N])) {
      k = i - N
      l = j - N
      vY[i,j] = dYmAspau[idouro[k]] * dYmAspau[idouro[l]] * sAspau[idouro[k]]^2 + 
                dYmAspCdau[idouro[k]] * dYmAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYmFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYmGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYmGepiau[idouro[k]] * dYmGepiau[idouro[l]] * sGepiau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YInferior e YInferior - Parte Diagonal
for(i in (N * 2 + 1):(N * 3)) {
  for(j in (i):(N * 3)) {
    if(isTRUE(Isotopo[i - N * 2] == Isotopo[j - N * 2]) & isTRUE(Egama[i - N * 2] == Egama[j - N * 2]) & isTRUE(idouro[i - N * 2]==idouro[j - N * 2])) {
      k = i - N * 2
      l = j - N * 2
      vY[i,j] = dYiAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYiAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 + 
                dYiFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 + 
                dYiGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYiEficau[idouro[k]] * dYiEficau[idouro[l]] * sEficau[idouro[k]]^2 +
                dYiAsp[k] * dYiAsp[l] * sAsp[k]^2 +
                #dYiEfic[k] * dYiEfic[l] * sEfic[k]^2 +
                dYiAspCd[k] * dYiAspCd[l] * sAspCd[k]^2 
                dYiGth[k] * dYiGth[l] * sGth[k]^2 + 
                dYiFCd[k] * dYiFCd[l] * sFCd[k]^2
      vY[j,i] = vY[i,j]
    }
  }
}


#covariancias entre YInferior e YInferior - Parte nao Diagonal, mesmo isotopo, gama diferente
for(i in (N * 2 + 1):(N * 3)) {
  for(j in (i + 1):(N * 3)) {
    if(isTRUE(Isotopo[i - N * 2] == Isotopo[j - N * 2]) & isTRUE(Egama[i - N * 2] != Egama[j - N * 2]) & isTRUE(idouro[i - N * 2]==idouro[j - N * 2])) {
      k = i - N * 2
      l = j - N * 2
      vY[i,j] = dYiAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYiAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 + 
                dYiFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 + 
                dYiGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYiEficau[idouro[k]] * dYiEficau[idouro[l]] * sEficau[idouro[k]]^2 +
                dYiGth[k] * dYiGth[l] * sGth[k]^2 +
                #dYiEfic[k] * dYiEfic[l] * sEfic[k]^2 +
                dYiFCd[k] * dYiFCd[l] * sFCd[k]^2
      vY[j,i] = vY[i,j]
    }
  }
} 

#covariancias entre YInferior e YInferior - Parte nao Diagonal, isotopo diferente
for(i in (N * 2 + 1):(N * 3)) {
  for(j in (i + 1):(N * 3)) {
    if(isTRUE(Isotopo[i - N * 2] != Isotopo[j - N * 2]) & isTRUE(Egama[i - N * 2] != Egama[j - N * 2]) & isTRUE(idouro[i - N * 2]==idouro[j - N * 2])) {
      k = i - N * 2
      l = j - N * 2
      vY[i,j] = dYiAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYiAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 + 
                dYiFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 + 
                dYiGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2 + 
                dYiEficau[idouro[k]] * dYiEficau[idouro[l]] * sEficau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YSuperior e YMedio - Parte diagonal - com Egama igual
for(i in 1:N) {
  for(j in (N + 1):(N * 2)) {
    if(isTRUE(Isotopo[i]==Isotopo[j - N]) & isTRUE(Egama[i]==Egama[j - N]) & isTRUE(idouro[i]==idouro[j - N])) {
      k = j - N
      vY[i,j] = dYsAspCd[i] * dYmAspCd[k] * sAspCd[i]^2 +
                dYsFCd[i] * dYmFCd[k] * sFCd[i]^2 +
                dYsGepi[i] * dYmGepi[k] * sGepi[i]^2 
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YSuperior e YMedio - Parte nao diagonal - Isotopo Igual, com Egama diferente
for(i in 1:N) {
  for(j in (N + 1):(N * 2)) {
    if(isTRUE(Isotopo[i]==Isotopo[j - N]) & isTRUE(Egama[i]!=Egama[j - N]) & isTRUE(idouro[i]==idouro[j - N])) {
      k = j - N
      vY[i,j] = dYsFCd[i] * dYmFCd[k] * sFCd[i]^2 + 
                dYsGepi[i] * dYmGepi[k] * sGepi[i]^2 
      vY[j,i] = vY[i,j]
    }
  }
}


#covariancias entre YSuperior e YInferior - Parte diagonal - Isotopo Igual, com Egama Igual
for(i in 1:N) {
  for(j in (N * 2 + 1):(N * 3)) {
    if(isTRUE(Isotopo[i] == Isotopo[j - N * 2]) & isTRUE(Egama[i] == Egama[j - N * 2]) & isTRUE(idouro[i] == idouro[j - N * 2])) {
      k = j - N * 2
      vY[i,j] = dYsAspCd[i] * dYmAspCd[k] * sAspCd[i]^2 +
                dYsFCd[i] * dYmFCd[k] * sFCd[i]^2 +
                #dYsEfic[i] * dYiEfic[k] * sEfic[i]^2 +
                dYsGepi[i] * dYmGepi[k] * sGepi[i]^2 
                
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YSuperior e YInferior - Parte nao diagonal - Isotopo Igual, com Egama diferente
for(i in 1:N) {
  for(j in (N * 2 + 1):(N * 3)) {
    if(isTRUE(Isotopo[i] == Isotopo[j - N * 2]) & isTRUE(Egama[i] != Egama[j - N * 2]) & isTRUE(idouro[i] == idouro[j - N * 2])) {
      k = j - N * 2
      vY[i,j] = #dYsEfic[i] * dYiEfic[k] * sEfic[i]^2
                dYsFCd[i] * dYiFCd[k] * sFCd[i]^2
      vY[j,i] = vY[i,j]
    }
  }
}


#covariancias entre YMedio e YInferior - Parte Diagonal - Isotopo Igual, com Egama igual
for(i in (N + 1):(N * 2)) {
  for(j in (N * 2):(N * 3)) {
    if(isTRUE(Isotopo[i - N] == Isotopo[j - N * 2]) & isTRUE(Egama[i - N] == Egama[j - N * 2]) & isTRUE(idouro[i - N] == idouro[j - N * 2])) {
      k = i - N
      l = j - N * 2
      vY[i,j] = dYmAsp[k] * dYiAsp[l] * sAsp[k]^2 + 
                dYmAspCd[k] * dYiAspCd[l] * sAspCd[k]^2 + 
                dYmFCd[k] * dYiFCd[l] * sFCd[k]^2 + 
                dYmGth[k] * dYiGth[l] * sGth[k]^2 + 
                dYmAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 +
                dYmAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YMedio e YInferior - Parte Diagonal - Isotopo Igual, com Egama diferente
for(i in (N + 1):(N * 2)) {
  for(j in (N * 2):(N * 3)) {
    if(isTRUE(Isotopo[i - N] == Isotopo[j - N * 2]) & isTRUE(Egama[i - N] != Egama[j - N * 2]) & isTRUE(idouro[i - N] == idouro[j - N * 2])) {
      k = i - N
      l = j - N * 2
      vY[i,j] = dYmFCd[k] * dYiFCd[l] * sFCd[k]^2 + 
                dYmGth[k] * dYiGth[l] * sGth[k]^2 + 
                dYmAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 +
                dYmAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

#covariancias entre YMedio e YInferior - Parte Diagonal - Isotopo diferente, com Egama diferente
for(i in (N + 1):(N * 2)) {
  for(j in (N * 2):(N * 3)) {
    if(isTRUE(Isotopo[i - N] != Isotopo[j - N * 2]) & isTRUE(Egama[i - N] != Egama[j - N * 2]) & isTRUE(idouro[i - N] == idouro[j - N * 2])) {
      k = i - N
      l = j - N * 2
      vY[i,j] = dYmAspau[idouro[k]] * dYiAspau[idouro[l]] * sAspau[idouro[k]]^2 +
                dYmAspCdau[idouro[k]] * dYiAspCdau[idouro[l]] * sAspCdau[idouro[k]]^2 + 
                dYmFCdau[idouro[k]] * dYiFCdau[idouro[l]] * sFCdau[idouro[k]]^2 +
                dYmGthau[idouro[k]] * dYiGthau[idouro[l]] * sGthau[idouro[k]]^2
      vY[j,i] = vY[i,j]
    }
  }
}

varsup = (dYsAspCd * sAspCd)^2 + (dYsEfic * sEfic)^2 + (dYsFCd * sFCd)^2 + (dYsGepi * sGepi)^2

varmed = 0
for (i in 1 : N) {
  varmed[i] = (dYmAsp[i] * sAsp[i])^2 +
              (dYmAspCd[i] * sAspCd[i])^2 +
              (dYmFCd[i] * sFCd[i])^2 + 
              (dYmGth[i] * sGth[i])^2 +
              (dYmGepi[i] * sGepi[i])^2 +
              (dYmAspau[idouro[i]] * sAspau[idouro[i]])^2 + 
              (dYmAspCdau[idouro[i]] * sAspCdau[idouro[i]])^2 + 
              (dYmFCdau[idouro[i]] * sFCdau[idouro[i]])^2 +
              (dYmGepiau[idouro[i]] * sGepiau[idouro[i]])^2 +
              (dYmGthau[idouro[i]] * sGthau[idouro[i]])^2
}

varinf = 0
for (i in 1 : N) {
  varinf[i] = (dYiAsp[i] * sAsp[i])^2 + 
              (dYiAspCd[i] * sAspCd[i])^2 + 
              (dYiFCd[i] * sFCd[i])^2 + 
              (dYiGth[i] * sGth[i])^2 + 
              (dYiEfic[i] * sEfic[i])^2 +
              (dYiAspau[idouro[i]] * sAspau[idouro[i]])^2 + 
              (dYiAspCdau[idouro[i]] * sAspCdau[idouro[i]])^2 + 
              (dYiFCdau[idouro[i]] * sFCdau[idouro[i]])^2 +
              (dYiGthau[idouro[i]] * sGthau[idouro[i]])^2 + 
              (dYiEficau[idouro[i]] * sEficau[idouro[i]])^2
}

varY = c(varsup,varmed, varinf)
# vY = diag(varY)
invVy = inv(vY)

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#construindo a matrix de planejamento
#Definindo as variaveis para o chute inicial
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

A2 = matrix(c(Aa, Alfa, Q0b, k0b))
Nb = nrow(A2)

#valores de lambda e parametros para o loop na hora de fazer a analise
chidif = 1
lambda = 0.0001
chi2 = -1
recalcularY = 1
chi2novo = 0

while(1) {
  if(recalcularY == 1) {
    
    Aa = A2[1]
    Alfa1 = A2[2]
    Q0b = 0
    for (i in 1:nQ0) {
      Q0b[i] = A2[i+2]
    }
    for (i in 1 : N) {
       k0b[i] = A2[(i + nQ0 + 2)]
  }

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
# DERIVADAS PARCIAIS DO YSUPERIOR
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

#Yajusupexp = expression((exp(Aa) * Eres^(2*Alfa1) * k0b) * (((Q0b-0.429)/Eres^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))))
Yajusupexp = expression(Aa + (2 * Alfa1 * log(Eres)) + log(k0b) + log((((Q0b-0.429)/Eres^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1)))))
            

colsa = 0
for (i in 1:N) {
  #colsa[i] = exp(Aa) * Eres[i]^(2 * Alfa1) * k0b[i] * (((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
  colsa[i] = 1
}
colsa = matrix(colsa)

colsalfa = 0
for (i in 1:N) {
  #colsalfa[i] = exp(Aa) * (Eres[i]^(2 * Alfa1) * (log(Eres[i]) * 2)) * k0b[i] * (((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1))) - (exp(Aa) * Eres[i]^(2 * Alfa1) * k0b[i]) * (0.429 * (2 * 0.55^Alfa1 + (2 * Alfa1 + 1) * (0.55^Alfa1 * log(0.55)))/((2 * Alfa1 + 1) * 0.55^Alfa1)^2 + (Q0b[IN[i]] - 0.429) * (Eres[i]^Alfa1 * log(Eres[i]))/(Eres[i]^Alfa1)^2)
  colsalfa[i] = 2 * log(Eres[i]) - (0.429 * (2 * 0.55^Alfa1 + (2 * Alfa1 + 1) * (0.55^Alfa1 * log(0.55)))/((2 * Alfa1 + 1) * 0.55^Alfa1)^2 + (Q0b[IN[i]] - 0.429) * (Eres[i]^Alfa1 * log(Eres[i]))/(Eres[i]^Alfa1)^2)/(((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
  #colsalfa[i] = (-(log(Eres[i])*(Q0b[IN[i]]-0.429))/Eres[i]^Alfa1+0.2564720733241612/((2*Alfa1+1)*0.55^Alfa1)-0.858/((2*Alfa1+1)^2*0.55^Alfa1))/((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1))+2*log(Eres[i])
}
colsalfa = matrix(colsalfa)

colsQ0exp = D(Yajusupexp, 'Q0b')
colsQ0=0
for (i in 1:N) {
  #colsQ0[i] = (exp(Aa) * Eres[i]^(2 * Alfa1) * k0b[i]) * (1/Eres[i]^Alfa1)
  colsQ0[i] = 1/Eres[i]^Alfa1/(((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
}
colsq0Diag = monta_matriz_por_isotopo(colsQ0)
colsq0Diag2 = diag(colsQ0)

colsK0 = 0
for (i in 1:N) {
  #colsK0[i] = exp(Aa) * Eres[i]^(2 * Alfa1) * (((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
  colsK0[i] = 1/k0b[i]
}
colsK0 = diag(colsK0)
# colsK0 = eval(D(Yajusupexp, 'k0b'))


#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
# DERIVADAS PARCIAIS DO YMEDIO
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

Yajumedexp = expression(log((((Q0b-0.429)/Eres^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))) / (((Q0au-0.429)/Eresau^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1)))))

colma = matrix(0, N)

D(Yajumedexp, 'Alfa1')
colmalfa = 0
for (i in 1 : N) {
    colmalfa[i] = -(((0.429 * (2 * 0.55^Alfa1 + (2 * Alfa1 + 1) * (0.55^Alfa1 * log(0.55)))/((2 * Alfa1 + 1) * 0.55^Alfa1)^2 + (Q0b[IN[i]] - 0.429) * (Eres[i]^Alfa1 * log(Eres[i]))/(Eres[i]^Alfa1)^2)/(((Q0au[idouro[i]] - 0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1))) - (((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1))) * (0.429 * (2 * 0.55^Alfa1 + (2 * Alfa1 + 1) * (0.55^Alfa1 * log(0.55)))/((2 * Alfa1 + 1) * 0.55^Alfa1)^2 + (Q0au[idouro[i]] - 0.429) * (Eresau[idouro[i]]^Alfa1 * log(Eresau[idouro[i]]))/(Eresau[idouro[i]]^Alfa1)^2)/(((Q0au[idouro[i]] - 0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))^2)/((((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))/(((Q0au[idouro[i]] - 0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))))
    #colmalfa[i] = (((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1))*((-(log(Eres[i])*(Q0b[IN[i]]-0.429))/Eres[i]^Alfa1+0.2564720733241612/((2*Alfa1+1)*0.55^Alfa1)-0.858/((2*Alfa1+1)^2*0.55^Alfa1))/((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1))-((-(log(Eresau[idouro[i]])*(Q0au[idouro[i]]-0.429))/Eresau[idouro[i]]^Alfa1+0.2564720733241612/((2*Alfa1+1)*0.55^Alfa1)-0.858/((2*Alfa1+1)^2*0.55^Alfa1))*((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1)))/((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1))^2))/((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1+0.429/((2*Alfa1+1)*0.55^Alfa1))
}
colmalfa = matrix(colmalfa)


 # colmQ0 = 0
 # for (i in 1 : N) {
 #   colmQ0[i] = 1/Eres[i]^Alfa1/(((Q0[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
 # }
colmQ0 = 0
for (i in 1 : N) {
   colmQ0[i] = 1/Eres[i]^Alfa1/(((Q0au[idouro[i]] - 0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))/((((Q0b[IN[i]] - 0.429)/Eres[i]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))/(((Q0au[idouro[i]]  - 0.429)/Eresau[idouro[i]] ^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1))))
  #colmQ0[i] = 1/Eres[i]^Alfa1/(((Q0au[idouro[i]] - 0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2 * Alfa1 + 1) * 0.55^Alfa1)))
}
colmQ0Diag = monta_matriz_por_isotopo(colmQ0)
colmQ0Diag2 = diag(colmQ0)


# for (i in 1 : N) {
# colmQ0[i] = 1/(Eres[i]^Alfa1*(((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1)+(0.429/((2*Alfa1+1)*0.55^Alfa1))))
# }

colmK0 = diag(rep(0, nrow(colmQ0Diag)))

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
# DERIVADAS PARCIAIS DO INFERIOR
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

colialfa = matrix(0, nrow(colmQ0Diag))

colia = matrix(0, nrow(colmQ0Diag))

coliQ0 = matrix(0, N, ncol(colsq0Diag))

coliK0 = 0
for (i in 1:N) {
  #coliK0[i] = 1
  coliK0[i] = 1/k0b[i]
}
coliK0 = diag(coliK0)
#coliK0 = diag(rep(1, N))

col1 = rbind(colsa, colma, colia)
col2 = rbind(colsalfa, colmalfa, colialfa)
col3 = rbind(colsq0Diag, colmQ0Diag, coliQ0)
col4 = rbind(colsK0, colmK0, coliK0)

X = cbind(col1, col2, col3, col4)

vpar = Ginv(t(X) %*% Ginv(vY, tol = 0.00000000005*sqrt(.Machine$double.eps)) %*% X)
par = vpar %*% t(X) %*% Ginv(vY, tol = 0.00000000005*sqrt(.Machine$double.eps)) %*% Y

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#COMECANDO O LEVENBERG
#VETOR A COM A VARIAVEIS
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

Yexp = matrix(Y)
Yajusup = 0
Yajumed = 0
Yajuinf = 0

R = t(X) %*% invVy %*% X

    for (i in 1 : N) {
      Yajusup[i] = Aa + (2 * Alfa1 * log(Eres[i])) + log(k0b[i]) + log((((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))))
      #Yajusup[i] = exp(Aa) * Eres[i]^(2*Alfa1) * k0b[i] * (((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1)))
      #Yajumed[i] = (((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))) / (((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1)))
      Yajumed[i] = log((((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))) / (((Q0au[idouro[i]]-0.429)/Eresau[idouro[i]]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))))
      #Yajuinf[i] = k0b[i]
      Yajuinf[i] = log(k0b[i])
    }
    Yaju2 = matrix(c(Yajusup, Yajumed, Yajuinf))
    datb = Yexp - Yaju2
    Yy = Yexp - Yaju2
    ss = cbind(Yexp, Yaju2)
    if (chi2 == -1) {
      chi2 = t(datb) %*% invVy %*% datb
    }
  }
  Rlambda = matrix(0, nrow = nrow(R), ncol = ncol(R))
  for (i in 1: nrow(Rlambda)) {
    for (j in 1: ncol(Rlambda)) {
      if (isTRUE(i == j)) {
        Rlambda[i,j] = (1+lambda) * R[i,j]
      } else {
          Rlambda[i,j] = R[i,j]
        }
    }
  }
  invRlambda = Ginv(Rlambda)
  DA = invRlambda %*% t(X) %*% invVy %*% Yy
  Ynovo = X %*% DA
  DD = Yy - Ynovo
  chi2novo = t(DD) %*% invVy %*% DD
  chidif = chi2novo - chi2
  if (chidif >= 0) {
    lambda = lambda * 6
    recalcularY = 1
  } else {
      lambda = lambda / 4
      Anovo = A2 + DA
      A2 = Anovo
      chi2 = chi2novo
      recalcularY = 1
  }
  if (isTRUE(abs(chidif) < 2.0) & isTRUE(abs(chidif) > 0.0000001)) {
    break
  }
}

Afinal = round(A2, digits = 9)

numpar = ncol(X)
Ntotal = N * 3
gl = Ntotal - numpar
quirednovo = abs(chi2novo/gl)
qchisq(0.05, df=gl, lower.tail=FALSE) #valor critico do xhi-quadrado com nivel de significancia de 0.05
pchisq(chi2novo, df=gl, lower.tail=FALSE) # p-valor 


a = Afinal[1]
alfa = Afinal[2]
Q0 = 0
for (i in 1:nQ0) {
  Q0[i] = matrix(Afinal[i+2])
}
k0=0
for (i in 1 : N) {
  k0[i] = matrix(Afinal[(i + nQ0 + 2)])
}

Q0 = 0
for (i in 1:N) {
  Q0[i] = Q0b[IN[i]]
}


#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#
#VERIFICACAO DE ALFA E F UTILIZANDO FUNCOES DO R
#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

param = c(unique(Isotopo), Isotopo)
parQ0 = 0
for (i in 1 : nQ0) {
  parQ0[i] = 'Q0'
}

park0 = 0
for (i in 1 : N) {
  park0[i] = 'k0'
}


#UTILIZANDO FUNCAO PARA VERIFICAR O INTERCEPTO DAS RETAS DO METODOS CD-COVERED MULTI-MONITOR (CDCMM) E CD-RATIO MULTI-MONITOR (CDRMM)

ycdcmm =log(Eres^-Alfa1*AspCd/(k0*Efic*FCd*Q0*Gepi))
xcdcmm = log(Eres)
cdcmm = cbind(ycdcmm, xcdcmm)
cdcmm = as.data.frame(cdcmm)
lm_cdcmm = tidy(summary(lm(cdcmm)))

ycdrmm = log(Eres^alfa*Gth/((FCd*(Asp/AspCd)-1)*Q0*Gepi))
xcdrmm = log(Eres)
cdrmm = cbind(ycdrmm, xcdrmm)
cdrmm = as.data.frame(cdrmm)
summary(lm(cdrmm))

sY = sqrt(diag(vY))

sAfinal = abs(matrix(sqrt(abs(diag(Ginv(Rlambda, tol = 0.00000000005*sqrt(.Machine$double.eps)))))))

alfa = lm_cdcmm$estimate[2]
sa = sAfinal[1]
salfa = lm_cdcmm$std.error[2]
sQ0 = 0
for (i in 1:nQ0) {
  sQ0[i] = matrix(sAfinal[i+2])
}
sk0=0
for (i in 1 : N) {
  sk0[i] = matrix(sAfinal[(i + nQ0 + 2)])
}


param2 = cbind(matrix(param), c(parQ0, park0), matrix(scientific(Afinal[c(-1,-2)], digits = 5)), matrix(scientific(sAfinal[c(-1,-2)], digits = 5)))
param3 = c('a', 'alfa')
param4 = cbind(param3, c(a, alfa), scientific(c(sa, salfa), digits = 4), c('-','-'))

paramfinal = rbind(param4, param2, colnames(c('Alvos','Valores', 'Valores', 'Incerteza')))
colnames(paramfinal) = c('Alvos','Valores', 'Valores', 'Incerteza')

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

#gerando graficos

#------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------------------#

vetorx=0
for (i in 1 : Ntotal) {
  vetorx[i] = 0.95+i*1.005
}

xplot=0
for (i in 1:Ntotal) {
  xplot[i] = i
}

pointname = c(Isotopo, Isotopo)
pointname2 = c(Egama, Egama)

plot(xplot, (DD/sY), main = 'RESIDUOS PONDERADOS PELO DESVIO PADRAO - k0 e Q0', ylab = 'Y-Yaju-sigmaY', xlab = 'elementos', ylim = c(-7, 7))
text(xplot, (DD/sY), labels = pointname, cex = 0.7, pos = 3)
text(xplot, (DD/sY), labels = pointname2, cex = 0.7, pos = 1)
abline(3, 0, lty=3)
abline(-3, 0, lty=3)

# ggplot(cdcmm, aes(x=xcdcmm, y=ycdcmm)) + geom_point() + geom_smooth(method="lm")
# ggplot(cdrmm, aes(x=xcdrmm, y=ycdrmm)) + geom_point() + geom_smooth(method="lm")

paramfinal

covcor = matrix(data = 0, nrow = nrow(abs(Ginv(Rlambda))), ncol = ncol(abs(Ginv(Rlambda))))
for (i in 1:numpar) {
  for (j in 1:numpar) {
    covcor[i,j] = Ginv(Rlambda, tol = 0.00000000005*sqrt(.Machine$double.eps))[i,j]/sqrt(Ginv(Rlambda, tol = 0.00000000005*sqrt(.Machine$double.eps))[i,i]*Ginv(Rlambda, tol = 0.00000000005*sqrt(.Machine$double.eps))[j,j])
  }
}

#, tol = 0.00000000005*sqrt(.Machine$double.eps)

string = 0
for (i in 1:numpar-2) {
  string[i] = '-'
}

rownames(covcor) = rbind(matrix(param3), matrix(str_c(param, string, c(parQ0, park0))))
colnames(covcor) = t(rbind(matrix(param3), matrix(str_c(param, string, c(parQ0, park0)))))

corrplot(covcor, method = 'color', type = 'lower')
covcor = round(covcor, digits = 4)

sink('Resultados.txt', append = FALSE, split = FALSE)
'VALORES DOS PARAMETROS'
'----------------------------------------------------------------'
'----------------------------------------------------------------'
paramfinal
'----------------------------------------------------------------'
'----------------------------------------------------------------'
covcor
sink()

print(paramfinal)

for (i in 1:N) {
  aa = sum(Yy[i]^2/Yexp[i])
}


# ass = 0
# for (i in 1:N) {
#   ass[i] = log((((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1))))  
# }
# 
# ass2 = 0
# for (i in 1:N) {
#   ass2[i] = log((((Q0b[IN[i]]-0.429)/Eres[i]^Alfa1) + (0.429/((2*Alfa1+1)*0.55^Alfa1)))) 
# }