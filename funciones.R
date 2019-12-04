library(BacArena)
library(sybilSBML)
library(animation)
library(parallel)

emptyHoodMod <- function(object, pos, n, m, x, y){
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  #xp = xp[xp>0 & xp<=n]
  #xp = xp[yp>0 & yp<=m]
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
}

emptyHoodMod2 <- function(object, pos, n, m, x, y){
  xp = c(x-1,x,x+1)
  yp = c(y-2,y-1,y)
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  #xp = xp[xp>0 & xp<=n]
  #xp = xp[yp>0 & yp<=m]
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
}


chemotaxisMod <- function(object, population, j, chemo){
	popvec <- population@orgdat[j,]
	attract <- population@media[[chemo]]@diffmat
	freenb <- emptyHoodMod(object, population@orgdat[,c('x','y')],population@n, population@m, popvec$x, popvec$y)	
	if(length(freenb) != 0)
	{
		conc <- sapply(freenb, function(x, attract)
		{
			npos = as.numeric(unlist(strsplit(x,'_')))
			return(attract[npos[2],npos[1]])#Aqui es donde se equivocaron
		}, attract=attract)
		abs <- freenb[which(conc==max(conc))]
		if(length(abs)!=1)
		{
			abs <- abs[sample(1:length(abs),1)]
		}
		npos = as.numeric(unlist(strsplit(abs,'_')))
    #eval.parent(substitute(population@orgdat[j,]$x <- npos[1]))
    #eval.parent(substitute(population@orgdat[j,]$y <- npos[2]))
		return(npos)
	}
}

estimate_lrw <- function(grid_n, grid_m){
  x=c(10*10, 25*25, 51*51, 61*61, 71*71, 81*81, 91*91, 101*101)
  y=c(3901, 29911, 160000, 230000, 330000, 430000, 580000, 710000)
  lm <- lm(y~x)
  #summary(lm)
  #plot(x,y)
  #abline(coef(lm))
  #abline(coef=c(0, lm$coefficients[2]))
  lrw <- as.numeric(lm$coefficients[2]*grid_n*grid_m + grid_n*grid_m*100)
  #lrw <- ((grid_n*grid_m)*18.5 + 20)*10 -> alternative function
  return(lrw)
}

simEnvMod<-function(object, time, lrw=NULL, continue=FALSE, reduce=FALSE, diffusion=TRUE, diff_par=FALSE, cl_size=2, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE, verbose=TRUE,potspace,guardar,res,vecb,E0,bf,bdes){
celdas<-0*potspace
celdas[1,]=10
corriente<-0
#La matriz que maneja el potencial en el medio

	if(length(object@media)==0) stop("No media present in Arena!")
	#Se especifica si la simulación se continuará o se iniciará desde 0.
 	switch(class(object),
 	        "Arena"={arena <- object; evaluation <- Eval(arena)},
 	        "Eval"={arena <- getArena(object); evaluation <- object},
 	        stop("Please supply an object of class Arena or Eval."))
 	if( sec_obj=="none" & any(sapply(object@specs, function(s){s@limit_growth})))
 	   warning("If growth is limitted by maximum weight for an organism (max_weight=TRUE) it is recommended to use minimize total flux (sec_obj='mtf').")
  
 	if(is.null(lrw)){lrw=estimate_lrw(arena@n,arena@m)}
	for(i in names(arena@specs)){
   		phensel <- arena@phenotypes[which(names(arena@phenotypes)==i)]
    		if(length(phensel)==0){
			test = getPhenotype(arena@specs[[i]], cutoff=pcut, fbasol=arena@specs[[i]]@fbasol)
			pvec = rep(0,length(arena@mediac))
	      		names(pvec) = arena@mediac
	     		pvec[names(test)] = test
	      		pvec <- paste(pvec,collapse='')
	      		names(pvec) = i
	      		arena@phenotypes <- c(arena@phenotypes,pvec)
    		}
  	}
  	if(class(object)!="Eval"){addEval(evaluation, arena)}
	#sublb corresponde a las concentraciones de todas las sustancias en las celdas donde hay microorganismos
	arena@sublb <- getSublb(arena)
	arena@exchanges <- data.frame() # remember exchanges
	diff_t=0
	#Se contempla la posibilidad de un reactor agitado
	if(arena@stir){ #create all possible positions on arena
	    allxy = expand.grid(1:arena@n,1:arena@m)
	    colnames(allxy) = c("x","y")
	  }
 	if(length(arena@specs) > 0) biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
	arena0=arena 
	sumelec=0
	Farad<-96500 #Constante de Faraday
	sumelec<-0
  pos <- arena@orgdat[,c('x','y')]
		        for(j in 1:nrow(arena@orgdat)){
m<-arena@m
n<-arena@n
 x=pos[j,1]
 y=pos[j,2]
  xp = c(x-2,x-1,x,x+1,x+2)
  yp = c(y-2,y-1,y,y+1,y+2)
  xp=ifelse(xp<=0,x,xp)
  xp=na.omit(ifelse(xp>n,x,xp))
  yp=ifelse(yp<=0,y,yp)
  yp=na.omit(ifelse(yp>m,y,yp))

if(sum(celdas[yp[1:2],xp])>0)
	{
		celdas[yp,xp]<-celdas[yp,xp]+1
	}
#La escencia del biofilm. Si el potencial en dirección al electrodo cerca de una celula es alto, se establece una red y el biofilm crece, porque la celula tiene la capacidad de conducir

}
		#Se crean los directorios para guardar la información de la simulación 
		#Para guardar la información de cada sustancia en el medio
		rutaS<-paste(guardar,'/sustrato',sep='')
		#Para guardar la información de cada microorganismo
		rutaMO<-paste(guardar,'/MO',sep='')
		#Un resumen de la simulación en cada momento, contiene la misma información del anterior item pero escrita de otra forma. 
		rutasim<-paste(guardar,'/simlist',sep='')
		rutapoMO<-paste(guardar,'/posMO',sep='')
		rutaAC<-paste(rutaS,'/ac',sep='')
		rutaPS<-paste(guardar,'/potencial',sep='')
		dir.create(rutaPS)
		dir.create(rutaS)
		dir.create(rutaMO)
		dir.create(rutasim)
		dir.create(rutapoMO)
		dir.create(rutaAC)
i=0
	pos <- arena@orgdat
	rutasa<-paste(rutapoMO,'/',toString(i),'.csv',sep='')
	write.table(pos,file=rutasa,row.names=FALSE,col.names=FALSE,sep=',')
	rutaACd<-paste(rutaAC,'/',toString(i),'.csv',sep='')
	matAC<-matrix(nrow=m,ncol=n,arena@media[['EX_ac__40__e__41__']]@diffmat)
	write.table(matAC,file=rutaACd,row.names=FALSE,col.names=FALSE,sep=',')
	rutpsT<-paste(rutaPS,'/',toString(i),'.csv',sep='')
	write.table(potspace,file=rutpsT,row.names=FALSE,col.names=FALSE,sep=',')
	#corriente<-Farad*arena@n*arena@m*(mean(arena@media[['EX_e(e)']]@diffmat)-sumelec)/arena@tstep/3600*1e-6
	file<-file(paste(guardar,'/Corriente.csv',sep=''),'a')
	writeLines(paste(toString(arena@tstep*(i)),toString(corriente),sep=' '),file)
	close(file)
			m='\t biomass, type, phenotype, x, y '
			#Se crea un archivo para almacenar toda la información de todos los microorganismos
  pos <- arena@orgdat[,c('x','y')]


			for (x in 1:length(arena@orgdat[,1])){
 xb=pos[x,1]
 yb=pos[x,2]
				file2=file(paste(rutaMO,'/',toString(rownames(arena@orgdat[toString(x),])),'.csv',sep=''),'a')
				writeLines(paste(toString(arena@tstep*(i)),',',toString(arena@orgdat[toString(x),]),',',toString(potspace[yb,xb]),sep=' '),file2)
				close(file2)
				#Se van conectando, para guardarlo en un solo archivo simlist
				m=paste(m,paste(toString(rownames(arena@orgdat[toString(x),])),',',toString(arena@orgdat[toString(x),]),',',toString(potspace[yb,xb]),sep=''),sep='\n')
			}
			#El archivo lleva el conteo de los microorganismos presentes hasta ese momento
			file=file(paste(rutasim,'/',toString(i),'-MO',toString(length(arena@orgdat[,1])),'.csv',sep=''),'a')
			writeLines(m,file)	
			close(file)
#Importante revisar	#Se guarda la suma de toda la biomasa, pero se deberia diferenciar por tipo de microorganismo
			file3=file(paste(rutaMO,'/BO','.csv',sep=''),'a')
			writeLines(paste(toString(arena@tstep*(i)),toString(sum(arena@orgdat[,'biomass'])) ,toString(length(arena@orgdat[,'biomass'])),str=' '),file3)	
			close(file3)
			for (x in arena@media)	{
				file=file(paste(rutaS,'/',x@name,'.csv',sep=''),'a')	
				writeLines(paste(toString(arena@tstep*(i)),toString(mean(x@diffmat)),sep=' '),file)
#Importante revisar		#Se almacena el dato de unicamente la concentración promedio, pero seria mejor crear una carpeta para cada sustancia y tener toda la información en todo el espacio
				close(file)
			}
	sublb <- arena@sublb
	upbnd<-arena@specs[[arena@orgdat[1,'type']]]@ubnd
	upbnd[['EX_e(e)']]<-1000
	lb=-sublb[1,arena@specs[[arena@orgdat[1,'type']]]@medium]
	lobnd <- arena@specs[[arena@orgdat[1,'type']]]@lbnd
	reacts<-unique(arena@specs[[arena@orgdat[1,'type']]]@medium)
	lobnd[reacts]<- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb)
	#print(lobnd[reacts])
	#print(upbnd[reacts])
	#print(upbnd)
	#print(lobnd)
	upbnd[['Biomass']]<-1e-6
	rem<-optimizeLP(arena@specs[[arena@orgdat[1,'type']]],sec_obj='mtf', 	ub=upbnd,lb=lobnd,j=1)[[1]]$fluxes
	#print(rem[reacts])
	#print(rem['Biomass'])
	flu<-rem
	flu2<-flu


	for ( x in 1:1000){
	flu<-cbind(flu,flu2)}

	contE<-matrix(nrow=1000,ncol=1,0)


	for(i in 1:time){

matSusF<-0*potspace
corriente<-0
celdas[1,]=10
#Debo revisar detenidamente que hace esta parte del codigo
    		init_t <- proc.time()[3]

		sublb <- arena@sublb
#print('content')
#print(sublb[,'EX_ac__40__e__41__'])
    		if(nrow(arena@orgdat) > 1){
			new_ind = sample(1:nrow(arena@orgdat),nrow(arena@orgdat)) #shuffle through all bacteria to increase randomness
      			arena@orgdat = arena@orgdat[new_ind,]
			sublb = sublb[new_ind,] #apply shuffeling also to sublb to ensure same index as orgdat
		}
		#if(verbose) cat("\niteration:", i, "\t organisms:",nrow(arena@orgdat), "\t biomass:", sum(arena@orgdat$biomass), "pg \n")
		org_stat <- sapply(seq_along(arena@specs), function(x){dim(arena@orgdat[which(arena@orgdat$type==x),])[1]})
		if(length(arena@specs) > 0){
			old_biomass<-biomass_stat
			biomass_stat <- sapply(seq_along(arena@specs), function(x){sum(arena@orgdat$biomass[which(arena@orgdat$type==x)])})
      			org_stat <- cbind(org_stat, biomass_stat, 100*(biomass_stat-old_biomass)/old_biomass); rownames(org_stat) <- names(arena@specs); colnames(org_stat) <- c("count", "biomass", "%")
      			#if(verbose) print(org_stat)
}
    		arena@mflux <- lapply(arena@mflux, function(x){numeric(length(x))}) # empty mflux pool
    		arena@shadow <-lapply(arena@shadow, function(x){numeric(length(x))}) # empty shadow pool
    		if(nrow(arena@orgdat) > 0){ # if there are organisms left
      			org.count <- nrow(arena@orgdat)
			#Aquí se analiza cada microorganismos, uno por uno
print(i*arena@tstep)
		        for(j in 1:org.count){ # for each organism in arena
        			if(!verbose) cat("\rOrganims",j,"/",org.count) #Para no mostrar el conteo de todos los microorganismos
			        org <- arena@specs[[arena@orgdat[j,'type']]]
			        xp<-arena@orgdat[j,'x']
			        yp<-arena@orgdat[j,'y']
        			potencial<-potspace[yp,xp]
			        bacnum <-1 #round((arena@scale/(org@cellarea*10^(-8)))) #calculate the number of bacteria individuals per gridcell
L<-length(arena@media)
m<-arena@m
n<-arena@n
				#matSusF=array(dim=c(m,n,L),0)

bioM<-arena@orgdat[j,'biomass']
if (bioM>0){
#print(Sys.time())
#print(arena@orgdat[j,])
        			switch(class(org),
					#La variable sublb se modifica dentro de la función, el potencial que siente la celula se define por la
          			     "Bac"= {arena = simBacMod(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow,potencial=potencial,corriente=corriente,potspace=potspace,E0=E0,bf=bf,matSusF=matSusF,celdas=celdas,flu=flu,i=i,contE=contE,rem=rem)}, #the sublb matrix will be modified within this function
			             "Human"= {arena = simHum(org, arena, j, sublb, bacnum, sec_obj=sec_obj, cutoff=cutoff, pcut=pcut, with_shadow=with_shadow)}, #the sublb matrix will be modified within this function
              			stop("Simulation function for Organism object not defined yet."))

#print('alfinal')
#print(sublb[,'EX_ac__40__e__41__'])
}
			      }
			 test <- is.na(arena@orgdat$biomass)
		         if(sum(test)!=0) arena@orgdat <- arena@orgdat[-which(test),]
      rm("test")
    }
    if(verbose) cat("\r")
maxbf<-matrix(ncol=length(potspace[1,]),nrow=1,0)
for (cf in 1:length(potspace[1,]))
	{
		maxbf[cf]<-m+sum(potspace[,cf])/1000
	}
#print(maxbf)
#print(max(maxbf))
Lbf<-round(max(maxbf))
    if(diffusion && !arena@stir){
      if(diff_par){
        diff_t <- system.time(arena <- diffuse_parMod(arena, cluster_size=cl_size, lrw=lrw, sublb=sublb) )[3]
      }else {#diff_t <- system.time(arena <- diffuseMod(arena, lrw=lrw, sublb=sublb, verbose=verbose,matSusF) )[3] #Usar función modificada
#  for(s in seq_along(arena@media)){	
#	dx<-arena@Lx/arena@n
#Cmat<-arena@media[[s]]@diffmat
#	mat<-estadoestable(C=Cmat,J=matSusF[,,s],D=arena@media[[s]]@difspeed,dx=dx)
      	submat <- matrix(arena@media[['EX_ac__40__e__41__']]@diffmat, nrow=object@m, ncol=object@n)
        submat[sublb[,c("y","x")]] <- sublb[,'EX_ac__40__e__41__']
	dx<-arena@Lx/arena@n
	Cmat<-Matrix::Matrix(submat, sparse=TRUE)
	Cmat[1:(Lbf+1),]<-estadoestable(C=Cmat[1:(Lbf+1),],J=matSusF[,,'EX_ac__40__e__41__'],D=arena@media[['EX_ac__40__e__41__']]@difspeed,dx=dx)
	solucion<-Cmat[Lbf:m,]
        sumc <- sum(solucion) #sum of all concentrations
	largo<-length(solucion[,1])
	ancho<-length(solucion[1,])
        meanc <- sumc/(largo*ancho)
	Cmat[Lbf:m,]<-meanc
	arena@media[['EX_ac__40__e__41__']]@diffmat<-Cmat
	arena@sublb <- getSublb(arena)
#}
}

    }
    if(!diffusion){
      if(nrow(sublb)>0){
        for(k in 1:length(arena@media)){
          for(l in 1:nrow(sublb)){
            arena@media[[k]]@diffmat[sublb[l,"y"],sublb[l,"x"]] = sublb[l,k+2] # first two columns are coordinates
          }
        }
      }
      
    }
if (i%%1 ==FALSE)
{
    if(arena@stir){ #stir environment -> random movement of bacteria + perfect diffusion #Mezcla perfecta
      sublb_tmp = arena@orgdat[,c("x","y")]

      for(sub in names(arena@media)){ #go through each metabolite in medium

if(arena@media[[sub]]@id=='EX_ac__40__e__41__'){
      	submat <- matrix(arena@media[['EX_ac__40__e__41__']]@diffmat, nrow=object@m, ncol=object@n)
        submat[sublb[,c("y","x")]] <- sublb[,'EX_ac__40__e__41__']
	dx<-arena@Lx/arena@n
	Cmat<-Matrix::Matrix(submat, sparse=TRUE)
	Cmat[1:(Lbf+1),]<-estadoestable(C=Cmat[1:(Lbf+1),],J=matSusF[1:(Lbf+1),],D=arena@media[['EX_ac__40__e__41__']]@difspeed,dx=dx)
va=31059-m
	solucion<-Cmat[(Lbf+1):m,]
        sumc <- sum(solucion) #sum of all concentrations
	largo<-length(solucion[,1])
	ancho<-length(solucion[1,])
        meanc <- sumc/(largo*ancho)
	sumt<-sumc+Cmat[m,n]*va*ancho+sum(Cmat[Lbf,])
	meant<-sumt/((largo+va+1)*ancho)
	Cmat[Lbf:m,]<-meanc
	arena@media[['EX_ac__40__e__41__']]@diffmat<-Cmat
	arena@sublb <- getSublb(arena)
}
#else{
      #  sumc = sum(arena@media[[sub]]@diffmat) #sum of all concentrations
      #  meanc = sumc/(arena@n*arena@m) #mean per grid cell
      #  conc = ((sumc-(meanc*nrow(sublb)))+sum(sublb[,sub]))/(arena@n*arena@m) #remove concentrations where bacteria are sitting + add the current concentration in their position
      #  arena@media[[sub]]@diffmat = Matrix::Matrix(conc,nrow=arena@m,ncol=arena@n,sparse=TRUE) #create matrix with homogen concentration
      #  sublb_tmp[,sub] = conc #create a new sublb matrix
      #}
}
 for(j in 1:length(arena@orgdat[,1])){
org <- arena@specs[[arena@orgdat[j,'type']]]
if (org@type != 'geo'){
bioM<-arena@orgdat[j,'biomass']
if (bioM>0){
      arena@orgdat[j,c('x','y')] = allxy[sample(1:nrow(allxy),1),]
      #arena@orgdat[,c('x','y')] = newpos los microorganismos no se mueven igual
     # sublb_tmp[,c('x','y')] = newpos
      #arena@sublb = as.matrix(sublb_tmp)
}
}
}
    }
}
else{
arena@media[['EX_ac__40__e__41__']]@diffmat[sublb[,c("y","x")]]<-sublb[,'EX_ac__40__e__41__']
arena@sublb <- getSublb(arena)
}

numerodebichos<-length(arena@orgdat[arena@orgdat[,'biomass']>0,'biomass'])
#print(arena@orgdat[arena@orgdat[,'biomass']>0,])
#print(numerodebichos)
celdes<-round(bdes*numerodebichos*arena@tstep)
maxy<-max(arena@orgdat[,'y'])
desata<-sample(rownames(arena@orgdat[arena@orgdat[,'y']>(maxy-5),]),celdes)
arena@orgdat[desata,'biomass']<-0
arena@orgdat[desata,c('x','y')]<-1
#print('eliminar')
#print(celdes)
#print(bdes*sum(arena@orgdat[,'biomass']))
#print('biomasa')
#print(sum(arena@orgdat[,'biomass']))
#print('elegidos')
#print(desata)
#arena@sublb <- getSublb(arena)
#Para evitar que se pierda el gradiente y la chemotaxis funciones

#for (co in 1:m){
#	arena@media[['EX_elec']]@diffmat[co,]=matrix(nrow=1,ncol=n,(m-co))
#}

#Bacte <- matrix(nrow=m, ncol=n, 0)

#for ( cba in 1:length(arena@orgdat[,'x']))
#{
#	xba<-arena@orgdat[cba,'x']
#	yba<-arena@orgdat[cba,'y']
#	Bacte[yba,xba]<-Bacte[yba,xba]+1
#}

#print(arena@media[['EX_h(e)']]@diffmat)
#Calculos para estimar la eficiencia coulumbimetrica y la corriente
if (i%%res ==FALSE)
{
#print(arena@media[['EX_ac__40__e__41__']]@diffmat)
#print(i*arena@tstep)	
#print(Sys.time())
#print(i)
	pos <- arena@orgdat
	rutasa<-paste(rutapoMO,'/',toString(i),'.csv',sep='')
	write.table(pos,file=rutasa,row.names=FALSE,col.names=FALSE,sep=',')
	rutaACd<-paste(rutaAC,'/',toString(i),'.csv',sep='')
	matAC<-matrix(nrow=m,ncol=n,arena@media[['EX_ac__40__e__41__']]@diffmat)
	write.table(matAC,file=rutaACd,row.names=FALSE,col.names=FALSE,sep=',')
	rutpsT<-paste(rutaPS,'/',toString(i),'.csv',sep='')
	write.table(potspace,file=rutpsT,row.names=FALSE,col.names=FALSE,sep=',')
	#corriente<-Farad*arena@n*arena@m*(mean(arena@media[['EX_e(e)']]@diffmat)-sumelec)/arena@tstep/3600*1e-6
	file<-file(paste(guardar,'/Corriente.csv',sep=''),'a')
	writeLines(paste(toString(arena@tstep*(i)),toString(corriente),toString(Lbf),toString(celdes),toString(length(arena@orgdat[arena@orgdat[,'biomass']>0,'biomass'])),sep=','),file)
	close(file)
			m='\t biomass, type, phenotype, x, y '
			#Se crea un archivo para almacenar toda la información de todos los microorganismos
  pos <- arena@orgdat[,c('x','y')]


			for (x in rownames(arena@orgdat[arena@orgdat[,'biomass']>0,])){

 xb=pos[x,1]
 yb=pos[x,2]

	rutaMet<-paste(rutaMO,'/',toString(rownames(arena@orgdat[x,])),sep='')
dir.create(rutaMet)	
rutaMetp<-paste(rutaMet,'/',toString(i),'.csv',sep='')
write.table(cbind(names(flu[,as.integer(rownames(arena@orgdat[x,]))]),flu[,as.integer(rownames(arena@orgdat[x,]))]),file=rutaMetp,row.names=FALSE,col.names=FALSE,sep=',')



				file2=file(paste(rutaMO,'/',toString(rownames(arena@orgdat[x,])),'.csv',sep=''),'a')
				writeLines(paste(toString(arena@tstep*(i)),',',toString(arena@orgdat[x,]),',',toString(potspace[yb,xb]),sep=''),file2)
				close(file2)
				#Se van conectando, para guardarlo en un solo archivo simlist
				m=paste(m,paste(toString(rownames(arena@orgdat[x,])),',',toString(arena@orgdat[x,]),',',toString(potspace[yb,xb]),sep=''),sep='\n')
			}
			#El archivo lleva el conteo de los microorganismos presentes hasta ese momento
			file=file(paste(rutasim,'/',toString(arena@tstep*i),'-MO',toString(length(arena@orgdat[arena@orgdat[,'biomass']>0,'biomass'])),'.csv',sep=''),'a')
			writeLines(m,file)	
			close(file)
#Importante revisar	#Se guarda la suma de toda la biomasa, pero se deberia diferenciar por tipo de microorganismo
			file3=file(paste(rutaMO,'/BO.csv',sep=''),'a')
			writeLines(paste(toString(arena@tstep*(i)),',',toString(sum(arena@orgdat[,'biomass'])),',',toString(length(arena@orgdat[arena@orgdat[,'biomass']>0,'biomass'])),str=','),file3)	
			close(file3)
			for (x in arena@media)	{
				file=file(paste(rutaS,'/',x@name,'.csv',sep=''),'a')	
				writeLines(paste(toString(arena@tstep*(i)),toString(mean(x@diffmat)),sep=' '),file)
#Importante revisar		#Se almacena el dato de unicamente la concentración promedio, pero seria mejor crear una carpeta para cada sustancia y tener toda la información en todo el espacio
				close(file)
			}

	#Corriente en µA
#He omitido la sección de escritura de archivos para acelerar el calculo, dado que solo me interesa la posición de las bacterias
	#sumelec<-mean(arena@media[['EX_e(e)']]@diffmat)
	#sumcomp<-0
	#for(j in seq_along(arena@media)){
	#	if(vecb[j]>0)
	#	{
	#		sumcomp=sumcomp+vecb[j]*(mean(arena0@media[[j]]@diffmat)-mean(arena@media[[j]]@diffmat))
	#	}
	#}
	#Ec<-sumelec/sumcomp*100
	#file<-file(paste(guardar,'/EficienciaC.csv',sep=''),'a')
	#writeLines(paste(toString(arena@tstep*(i)),toString(Ec),sep=' '),file)
	#close(file)

	#m='\t biomass, type, phenotype, x, y '
	#for (x in 1:length(arena@orgdat[,1]))
	#{
	#	file2=file(paste(rutaMO,'/',toString(x),'.csv',sep=''),'a')
	#	writeLines(paste(toString(arena@tstep*(i)),',',toString(arena@orgdat[toString(x),]),sep=' '),file2)
	#	close(file2)
	#	m=paste(m,paste(toString(rownames(arena@orgdat[toString(x),])),',',toString(arena@orgdat[toString(x),]),sep=' '),sep='\n')	
	#}
	#file=file(paste(rutasim,'/',toString(i),'-MO',toString(length(arena@orgdat[,1])),'.csv',sep=''),'a')
	#writeLines(m,file)	
	#close(file)
	#file3=file(paste(rutaMO,'/BO','.csv',sep=''),'a')
	#writeLines(paste(toString(arena@tstep*(i)),toString(sum(arena@orgdat[,'biomass'])) ,toString(length(arena@orgdat[,'biomass'])),str=' '),file3)	
	#close(file3)

	#for (x in arena@media)
	#{
	#file=file(paste(rutaS,'/',x@name,'.csv',sep=''),'a')	
	#writeLines(paste(toString(arena@tstep*(i)),toString(mean(x@diffmat)),sep=' '),file)
	#close(file)
	#}

}
    #addEval(evaluation, arena) #esta cosa destruye la ram
#print(arena@media[['EX_ac__40__e__41__']]@diffmat)
#print(arena@media[['EX_h(e)']]@diffmat)
#print(potspace)
#print(celdas)
#print(Bacte)
#print(arena@orgdat)
	evaluation=Eval(arena)
    if(reduce && i<time){evaluation = redEval(evaluation)}
    if(nrow(arena@orgdat)==0 && !continue){
      if(verbose) print("All organisms died!")
      break
    }
    step_t <- proc.time()[3] - init_t
    #if(verbose) cat("\ttime total: ", round(step_t,3), "\tdiffusion: ", round(diff_t,3), " (", 100*round(diff_t/step_t,3),"%)\n" )
  }
  return(evaluation)
}
estadoestable <- function(C,J,D,dx)
	{	
		m<-length(C[,1])
		n<-length(C[1,])
		mat<-matrix(nrow=m+1,ncol=n+2,0)
		mat[2:(m+1),2:(n+1)]<-matrix(C,nrow=m,ncol=n)
		er<-1
#print(mat)
#print(J)
		#while(er>1e-10)
		for (te in 1:20)
			{	
				mat[1,]<-mat[2,]
				#C[m,]<-C[m-1,]
				mat[,1]<-mat[,2]
				mat[,n+2]<-mat[,n+1]
				Cold=mat
#print(mat)
				for (i in 2:m)
					{
						for (j in 2:(n+1))
							{
								mat[i,j]<-1/4*(Cold[i-1,j]+Cold[i+1,j]+Cold[i,j+1]+Cold[i,j-1]+J[i-1,j-1])
							}	
					}

				#print(C)
				er<-sum(((Cold-mat)/mat)^2)
			}
#print(er)
#print(C)
#print(J)
		C<-mat[2:(m+1),2:(n+1)]
		return(C)
	}

consumeMod<-function(object, sublb, cutoff=1e-6, bacnum, fbasol){
  #if(fbasol$obj>=cutoff && !is.na(fbasol$obj)){
    flux = fbasol[object@medium[object@medium !='EX_ac__40__e__41__']] #scale flux to whole population size
    #flux = na.omit(ifelse(abs(flux)<=cutoff,NA,flux)) ?
#print(sublb['EX_ac__40__e__41__'])
    sublb[names(flux)] = round(sublb[names(flux)]+flux, round(-log10(cutoff))) # use cutoff also in this case
 # }

#print(sublb['EX_ac__40__e__41__'])
  return(sublb)
}

simBacMod <-function(object, arena, j, sublb, bacnum, sec_obj="none", cutoff=1e-6, pcut=1e-6, with_shadow=FALSE,potencial=0, corriente=0,potspace,E0,bf,matSusF,celdas,flu,i,contE,rem){

nflu<-flu
potspa<-potspace
cells<-celdas
  pos <- arena@orgdat[,c('x','y')]
#print(object@type)
if (object@type == 'geo'){

 x=pos[j,1]
 y=pos[j,2]
n=arena@n
m=arena@m
  xp = c(x-2,x-1,x,x+1,x+2)
  yp = c(y-2,y-1,y,y+1,y+2)
  xp=ifelse(xp<=0,x,xp)
  xp=na.omit(ifelse(xp>n,x,xp))
  yp=ifelse(yp<=0,y,yp)
  yp=na.omit(ifelse(yp>m,y,yp))


#print('inicial')
#print(cells)
for (cx in xp){
 for (cy in yp)
{
	if(cells[cy,cx]>0){
	potspa[cy,cx]<-E0}
else{potspa[cy,cx]<--1000}	
}
}
potencial<-potspa[y,x]

}
const <- constrainMod(object, object@medium, lb=-sublb[j,object@medium], #scale to population size
                     dryweight=arena@orgdat[j,"biomass"], tstep=arena@tstep, scale=arena@scale, j)
    #La función constrain toca revisarla mucho mas detenidamente
    tstep=arena@tstep
    lobnd <- const[[1]]
    upbnd <- const[[2]]
    #print(lobnd)
    #print(object@lbnd)
    corriente0<-1000 #Calibrar
    const<-10 #Calibrar
if (object@type == 'geo'){ #Este error no lo recordaba, no me habia pasado antes
    upbnd[['EX_e(e)']]<-corriente0/(exp(-38.9*potencial)+1)
	
	#print(upbnd[['EX_e(e)']])
}
#print(arena@orgdat[j,])
#print(sublb[j,])
Dc<-0.00024 #cm
Dif<-arena@media[['EX_ac__40__e__41__']]@difspeed
#print(lobnd[['EX_ac__40__e__41__']])
   # upbnd[['EX_e(e)']]<-corriente0/(exp(-38.9*potencial)+1)
#Para analizar los perfiles metabolicos debo usar la función de sybil directamente!
popvec <- arena@orgdat[j,]
#print('pop')
#print(popvec)
#print(j)
jp<-as.integer(rownames(arena@orgdat[j,]))
c<-contE[jp]
    freenb <- emptyHoodMod(object, arena@orgdat[,c('x','y')],
              arena@n, arena@m, popvec$x, popvec$y)
#print(freenb)
 # popvec$biomass >= object@maxweight
#print(jp)
#print(freenb)
#print(arena@orgdat)
    if(length(freenb) == 0 & popvec$biomass >= object@maxweight/1.05){# & popvec$biomass > object@maxweight){

upbnd[['Biomass']]<-1e-6
c<-c+1
nflu[,jp]<-rem
#print('in')
}else{
c<-0}
if ((i==1) |(i%%10 ==FALSE &upbnd[['Biomass']]>1e-4)){
#print(upbnd[['Biomass']])
#print('opt')
  optimization <- optimizeLP(object, lb=lobnd, 	ub=upbnd, j=j, sec_obj=sec_obj, cutoff=cutoff, with_shadow=with_shadow)[[1]]$fluxes
  #eval.parent(substitute(flu[,jp]<-optimization))
nflu[,jp]<-optimization

}
else{
optimization<-nflu[,jp]
}

eval.parent(substitute(contE[jp]<-c))
#print(optimization['EX_ac__40__e__41__'])
#print(optimization['EX_ac__40__e__41__'])
#print(Dc)
#print(Dif)
#print(arena@orgdat[j,'biomass'])
eval.parent(substitute(matSusF[y,x] <-optimization['EX_ac__40__e__41__']/Dc/1e9/Dif*arena@orgdat[j,"biomass"]))

	miu<-optimization[object@rbiomass]
#print(optimization['EX_ac__40__e__41__'])
    optimization<-optimization*arena@orgdat[j,"biomass"]*tstep/(Dc**3/1000)/1e12

#Esta parte la cambie, sin embargo, toca revisar mas a fondo              #*arena@orgdat[j,"biomass"]/object@cellweight_mean
#También se debe tener en cuenta el volumen de la celda
  fbasol <- optimization
#print(c(i,j))
#print('hay')
#print(sublb[,'EX_ac__40__e__41__'])
#print(sublb[j,'EX_ac__40__e__41__'])
#print('come')
#print(optimization['EX_ac__40__e__41__'])
#print('biomass')
#print(miu)
#print(arena@orgdat[j,])
#sub <- consumeMod(object, sublb[j,], bacnum=bacnum, fbasol=fbasol, cutoff)
#print(sub['EX_ac__40__e__41__'])
#  eval.parent(substitute(sublb[j,] <- sub )) #scale consumption to the number of cells?
#sublb[j,] <- sub
#print(sublb[,'EX_ac__40__e__41__'])
#print('comio')
#print(sublb[j,'EX_ac__40__e__41__'])
    if(length(freenb) == 0 & popvec$biomass >= object@maxweight)
{
dead<-FALSE 
}
else
{
dead <- growthMod(object, arena, j, arena@occupyM, fbasol=fbasol, 
tstep=arena@tstep,celdas=cells,E0=E0,potspace=potspa,miu=miu,flu=nflu)
}


eval.parent(substitute(flu<-nflu))

# for(s in seq_along(arena@media)){	

#}

  #arena@orgdat[j,'phenotype'] <- as.integer(checkPhen(arena, org=object, fbasol=fbasol, cutoff=pcut))

  type <- object@type
  arena@mflux[[type]]  <- arena@mflux[[type]] + fbasol # remember active fluxes
#  arena@shadow[[type]] <- arena@shadow[[type]]+ optimization[[2]]
  idx <- match(arena@mediac, names(fbasol))
  exchanges <- data.frame(type,t(fbasol[idx]))
  colnames(exchanges) <- c("species", unname(arena@mediac))
  arena@exchanges <- rbind(arena@exchanges, exchanges) # remember exchanges
  
  if(dead && object@lyse){
    eval.parent(substitute(sublb[j,] <- lysis(object, sublb[j,])))
  }
  pos <- arena@orgdat[,c('x','y')]

 if(!dead && object@speed != 0){
    if(object@chem[1] == ''){
      mov_pos <- move(object, pos, arena@n, arena@m, j, arena@occupyM)
      arena@orgdat[,c('x','y')] <- mov_pos
    }else{
      for (v in seq_along(object@chem)){
      chemo <- object@chem[[v]]
      chemo_pos <- chemotaxisMod(object, arena, j, chemo)
      if(!is.null(chemo_pos)){arena@orgdat[j,c('x','y')] <- chemo_pos}
      }
    }
  }
if (object@type == 'geo'){


  eval.parent(substitute(corriente <- corriente+fbasol[['EX_e(e)']]*(Dc**3/1000)/tstep))
if(!dead){
if(sum(cells[yp,xp])>0)
	{
		cells[yp,xp]<-cells[yp,xp]-1
	}

#print('movimiento')
#print(cells)

for (cx in xp){
 for (cy in yp)
{
	if(cells[cy,cx]>0){
	potspa[cy,cx]<-E0}
else{potspa[cy,cx]<--1000}	
}
}

  pos <- arena@orgdat[,c('x','y')]

 x=pos[j,1]
 y=pos[j,2]
  xp = c(x-2,x-1,x,x+1,x+2)
  yp = c(y-2,y-1,y,y+1,y+2)
  xp=ifelse(xp<=0,x,xp)
  xp=na.omit(ifelse(xp>n,x,xp))
  yp=ifelse(yp<=0,y,yp)
  yp=na.omit(ifelse(yp>m,y,yp))

if(sum(cells[yp[1:2],xp])>0)
	{
		cells[yp,xp]<-cells[yp,xp]+1
	}

#print('nuevo')
#print(cells)


for (cx in xp){
 for (cy in yp)
{
	if(cells[cy,cx]>0){
	potspa[cy,cx]<-E0}
else{potspa[cy,cx]<--1000}	
}
}

    freenb2 <- emptyHoodMod2(object, arena@orgdat[,c('x','y')],
              arena@n, arena@m, popvec$x, popvec$y)

if(length(freenb2) != 0){
      npos = freenb2[sample(length(freenb2),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
   #   if(occupyM[npos[2], npos[1]] == 0){ # check if there is no obstacle #Por el momento no hay obtaculos en mi ambiente
        arena@orgdat[j,]$x <- npos[1]
        arena@orgdat[j,]$y <- npos[2]
}


}
nflu<-flu
eval.parent(substitute(celdas<-cells))
eval.parent(substitute(potspace<-potspa))
}

  return(arena)
}
growExpMod<-function(object, biomass, fbasol, tstep,miu){
  growth <-exp(miu*tstep)
vd<-matrix(nrow=length(fbasol),ncol=1,0)
s<-sum(fbasol==vd)/length(vd)

  if(s == 1){
	grow_accum <- biomass - object@deathrate*biomass*tstep
#print(growth)
  } else {
	grow_accum <- growth*biomass
}

  #cat("\t", growth, biomass, grow_accum, "\n")
  return(grow_accum)
}

emptyHoodMod <- function(object, pos, n, m, x, y){
  xp = c(x-1,x,x+1)
  yp = c(y-1,y,y+1)
  xp=ifelse(xp<=0,NA,xp)
  xp=na.omit(ifelse(xp>n,NA,xp))
  yp=ifelse(yp<=0,NA,yp)
  yp=na.omit(ifelse(yp>m,NA,yp))
  #xp = xp[xp>0 & xp<=n]
  #xp = xp[yp>0 & yp<=m]
  nb=sapply(xp,function(x,y){return(paste(x,y,sep='_'))},y=yp)
  pos = pos[which(pos$x %in% xp),]
  pos = pos[which(pos$y %in% yp),]
  freenb=setdiff(nb,paste(pos$x,pos$y,sep='_'))
  if(length(freenb)==0){return(NULL)}else{return(freenb)}
}

growthMod <-function(object, population, j, occupyM, fbasol, tstep,celdas,potspace,E0,miu,flu){
n=population@n
m=population@m
cells<-celdas
potspa<-potspace
  neworgdat <- population@orgdat
  popvec <- neworgdat[j	,]
  switch(object@growtype,
         "linear"= {popvec$biomass <- growLin(object, popvec$biomass, fbasol, tstep)},
         "exponential"= {popvec$biomass <- growExpMod(object, popvec$biomass, fbasol, tstep,miu)},
         stop("Growth type must be either linear or exponential"))
  dead <- F
  neworgdat[j,'biomass'] <- popvec$biomass
  while( popvec$biomass > object@maxweight ){
  #if(popvec$biomass > object@maxweight){ 
#print('aqui')
    freenb <- emptyHoodMod(object, neworgdat[,c('x','y')],
              population@n, population@m, popvec$x, popvec$y)
#print(freenb)
    if(length(freenb) != 0){

if (length(rownames(neworgdat[neworgdat[,'biomass']==0,]))>0){
	kO<-sample(rownames(neworgdat[neworgdat[,'biomass']==0,]),1)
}else{	kO<-nrow(neworgdat)+1}

      npos = freenb[sample(length(freenb),1)]
      npos = as.numeric(unlist(strsplit(npos,'_')))
   #   if(occupyM[npos[2], npos[1]] == 0){ # check if there is no obstacle #Por el momento no hay obtaculos en mi ambiente
        popvec$biomass <- popvec$biomass/2
        daughter <- popvec
        daughter$biomass <- popvec$biomass
        daughter$x <- npos[1]
        daughter$y <- npos[2]
#print('aqui')
        neworgdat[kO,] <- daughter
	jp<-as.integer(rownames(population@orgdat[j,]))
#print(c(kO,jp))
	eval.parent(substitute(flu[,as.integer(kO)]<-flu[,jp]))
        neworgdat[j,'biomass'] <- popvec$biomass

#print(neworgdat)
if (object@type == 'geo'){
x <- npos[1]
y <- npos[2]
  xp = c(x-2,x-1,x,x+1,x+2)
  yp = c(y-2,y-1,y,y+1,y+2)
  xp=ifelse(xp<=0,x,xp)
  xp=na.omit(ifelse(xp>n,x,xp))
  yp=ifelse(yp<=0,y,yp)
  yp=na.omit(ifelse(yp>m,y,yp))

if(sum(cells[yp[1:2],xp])>0)
	{
		cells[yp,xp]<-cells[yp,xp]+1
	}
}

     # }
    }else
      break
  }
  if(popvec$biomass < object@minweight){
pos=neworgdat[,c('x','y')]
if (object@type == 'geo'){
 x=pos[j,1]
 y=pos[j,2]
  xp = c(x-2,x-1,x,x+1,x+2)
  yp = c(y-2,y-1,y,y+1,y+2)
  xp=ifelse(xp<=0,x,xp)
  xp=na.omit(ifelse(xp>n,x,xp))
  yp=ifelse(yp<=0,y,yp)
  yp=na.omit(ifelse(yp>m,y,yp))


if(sum(cells[yp,xp])>0)
	{
		cells[yp,xp]<-cells[yp,xp]-1
	}

for (cx in xp){
 for (cy in yp)
{
if(cells[cy,cx]>0){
	potspa[cy,cx]<-E0}
else{potspa[cy,cx]<--1000}	
}
}



}
print('muerto')
    neworgdat[j,'biomass'] <- 0
    neworgdat[j,'x'] <- 1
    neworgdat[j,'y'] <- 1


    dead <- T 
  }

if (object@type == 'geo'){
eval.parent(substitute(celdas<-cells))
eval.parent(substitute(potspace<-potspa))
}
  eval.parent(substitute(population@orgdat <- neworgdat))
  return(dead)
}



diffuseMod <- function(object, lrw, sublb, verbose=TRUE,matSusF){
  arena <- object
  diff_init_t <- proc.time()[3]	
  sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
  #diff_pre_t <- system.time({
  if(!all(is.na(sublb)) & dim(sublb)[1] > 0){ # if there are organisms
    testdiff <- t(sublb[,-c(1,2)]) == unlist(lapply(arena@media,function(x,n,m){return(mean(x@diffmat))})) #check which mets in sublb have been changed by the microbes
    changed_mets <- which(apply(testdiff,1,sum)/nrow(sublb) < 1) #find the metabolites which are changed by at least one microbe
  } else changed_mets <- list()#})[3]
  #diff_pde_t=0; diff_sublb_t=0
  #diff_loop_t <- system.time({for(j in seq_along(arena@media)){
  for(j in seq_along(arena@media)){
    if(!verbose) cat("\rSubstances",j,"/",length(arena@media)) #Para no mostrar todas las sustancias
    #skip diffusion if already homogenous (attention in case of boundary/source influx in pde!)
    if(length(changed_mets)>0) homogenous = !(j %in% changed_mets) else homogenous = FALSE
    diffspeed  = arena@media[[j]]@difspeed>0
    #eval.parent(substitute(object@media[[j]]@pde<-"MFCDiff"))
    diff2d     = TRUE
    if(diff2d&&!homogenous || !diff2d){
      submat <- matrix(arena@media[[j]]@diffmat, nrow=object@m, ncol=object@n)
      if(!all(is.na(sublb)) && dim(sublb)[1] > 0 && (nrow(sublb) != sum(sublb[,j+2]==mean(submat)))){
        submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id]
      }
      #diff_pde_t <- diff_pde_t + system.time(switch(arena@media[[j]]@difunc,

      if(diffspeed || !diff2d){
k<-'EX_ac__40__e__41__'
	dx<-arena@Lx/arena@n
Cmat<-arena@media[[j]]@diffmat
	arena@media[[j]]@diffmat<-estadoestable(C=Cmat,J=matSusF[,,k],D=arena@media[[k]]@difspeed,dx=dx)
        #switch(arena@media[[j]]@difunc,
        #       "pde"  = {submat <- diffusePDE(arena@media[[j]], submat, gridgeometry=arena@gridgeometry, lrw, tstep=object@tstep)},
        #       "pde2" = {diffuseSteveCpp(submat, D=arena@media[[j]]@difspeed, h=1, tstep=arena@tstep)},
        #       "naive"= {diffuseNaiveCpp(submat, donut=FALSE)},
        #       "r"    = {for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
        #       stop("Diffusion function not defined yet.")
        #)
      }
      arena@media[[j]]@diffmat <- Matrix::Matrix(submat, sparse=TRUE)
    }else submat <- arena@media[[j]]@diffmat
    tryCatch({
      sublb_tmp[,j] <- submat[cbind(arena@orgdat$y,arena@orgdat$x)]
    }, error=function(cond){
      print(cond)
      browser()
    })
  }#})[3]
  if(verbose) cat("\r")
  sublb <- cbind(as.matrix(arena@orgdat[,c("y","x")]),sublb_tmp)
  colnames(sublb) <- c('y','x',arena@mediac)
  arena@sublb <- sublb
  
  #diff_t <- proc.time()[3] - diff_init_t
  #print(paste("diffusion time total", round(diff_t,3), "pre", round(diff_pre_t,3), "loop", round(diff_loop_t,3), "pde", round(diff_pde_t,3), "sublb", round(diff_sublb_t,3), "post", round(diff_post_t,3) ))

  
  #return(list(arena, sublb))
  return(arena)
}

diffusePDEMod <- function(object, init_mat, gridgeometry, lrw=NULL, tstep){
  if(is.null(lrw)){
    lrw=estimate_lrw(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N)}
  solution <- deSolve::ode.2D(y = init_mat, func = get(object@pde), times=c(1,1+tstep), parms = c(gridgeometry=gridgeometry, diffgeometry=object@diffgeometry, boundS=object@boundS),
                     dimens = c(gridgeometry$grid2D$x.N, gridgeometry$grid2D$y.N), method="lsodes", lrw=lrw)#160000
  diff_mat <- matrix(data=solution[2,][-1], ncol=ncol(init_mat), nrow=nrow(init_mat))
  return(diff_mat)
}

diffuse_parMod<-function(object, lrw, cluster_size, sublb){
  diff_init_t <- proc.time()[3]
  cl <- parallel::makeCluster(cluster_size, type="FORK")
  #parallel::clusterExport(cl, c(""))
  arena <- object
  sublb_tmp <- matrix(0,nrow=nrow(arena@orgdat),ncol=(length(arena@mediac)))
  #diff_pre_t <- system.time({ if(dim(sublb)[1] > 0){
  if(!all(is.na(sublb)) && dim(sublb)[1] > 0){ # if there are organisms
      testdiff <- t(sublb[,-c(1,2)]) == unlist(lapply(arena@media,function(x,n,m){return(mean(x@diffmat))})) #check which mets in sublb have been changed by the microbes
      changed_mets <- which(apply(testdiff,1,sum)/nrow(sublb) < 1) #find the metabolites which are changed by at least one microbe
  } else changed_mets <- list()#})[3]
  #diff_pde_t=0; diff_sublb_t=0
  #diff_loop_t <- system.time(parallel_diff <- parallel::mclapply(seq_along(arena@media), function(j){
  #parallel_diff <- parallel::mclapply(seq_along(arena@media), function(j){
  parallel_diff  <- parallel::parLapply(cl, seq_along(arena@media), function(j){
  #parallel_diff <- lapply(seq_along(arena@media), function(j){
  #diff_loop_t <- system.time(parallel_diff <-  foreach(j=seq_along(arena@media)) %dopar% {
    #skip diffusion if already homogenous (attention in case of boundary/source influx in pde!)
    if(length(changed_mets)>0) homogenous = !(j %in% changed_mets) else homogenous = FALSE
    diffspeed  = arena@media[[j]]@difspeed>0
    #diff2d     = arena@media[[j]]@pde=="Diff2d"
	diff2d=TRUE

    if(diff2d&&!homogenous || !diff2d){
      submat <- matrix(arena@media[[j]]@diffmat, nrow=arena@m, ncol=arena@n)
      if(!all(is.na(sublb)) && dim(sublb)[1] > 0 && (nrow(sublb) != sum(sublb[,j+2]==mean(submat)))){
        #diff_sublb_t <<- diff_sublb_t + system.time(submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id])[3]}
        submat[sublb[,c("y","x")]] <- sublb[,arena@media[[j]]@id]}
      #browser()
      #diff_pde_t <<- diff_pde_t + system.time(switch(arena@media[[j]]@difunc,
      if(diffspeed || !diff2d){
        switch(arena@media[[j]]@difunc,
               "pde"  = {submat <- diffusePDE(arena@media[[j]], submat, gridgeometry=arena@gridgeometry, lrw, tstep=object@tstep)},
               "pde2" = {diffuseSteveCpp(submat, D=arena@media[[j]]@difspeed, h=1, tstep=arena@tstep)},
               "naive"= {diffuseNaiveCpp(submat, donut=FALSE)},
               "r"    = {for(k in 1:arena@media[[j]]@difspeed){diffuseR(arena@media[[j]])}},
               stop("Diffusion function not defined yet."))#)[3]
      }
        diffmat_tmp <- Matrix::Matrix(submat, sparse=TRUE)
    }else{
      diffmat_tmp <- arena@media[[j]]@diffmat
      submat <- matrix(arena@media[[j]]@diffmat, nrow=arena@m, ncol=arena@n)
    }
    sublb_tmp  <- submat[cbind(arena@orgdat$y,arena@orgdat$x)]
    list("diffmat"=diffmat_tmp, "sublb"=sublb_tmp)
  #})#)[3]
  #}, mc.cores=cluster_size)#)[3]
  })
  parallel::stopCluster(cl)
  
  #diff_post_t <- system.time({ for(j in seq_along(arena@media)){
  for(j in seq_along(arena@media)){
      arena@media[[j]]@diffmat <- parallel_diff[[j]][[1]]
      sublb_tmp[,j] <- parallel_diff[[j]][[2]]
  }#})[3]
  sublb <- cbind(as.matrix(arena@orgdat[,c(4,5)]),sublb_tmp)
  colnames(sublb) <- c('x','y',arena@mediac)
  arena@sublb <- sublb
  #diff_t <- proc.time()[3] - diff_init_t
  #print(paste("diffusion time total", round(diff_t,3), "pre", round(diff_pre_t,3), "loop", round(diff_loop_t,3), "pde", round(diff_pde_t,3), "sublb", round(diff_sublb_t,3), "post", round(diff_post_t,3) ))
  return(arena)
}


MFCDiff <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    nrow <- gridgeometry.grid2D$x.N
    ncol <- gridgeometry.grid2D$y.N
    CONC  <- matrix(nrow , ncol, data = y)
    fxup=rep(0,ncol)
    fxdown=rep(0,ncol)
    fyup=rep(0,nrow)
    fydown=rep(0,nrow)
    dCONC <- ReacTran::tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid,
                               flux.x.up = fxup, flux.x.down = fxdown,  flux.y.up = fyup, flux.y.down = fydown)$dC
    return (list(dCONC))
  })
}
#Se modifican las condiciones de frontera para crear un ambiente aislado

MFCDiff2 <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    nrow <- gridgeometry.grid2D$x.N
    ncol <- gridgeometry.grid2D$y.N
    CONC  <- matrix(nrow , ncol, data = y)
    fxup=rep(0,ncol)
    dCONC <- ReacTran::tran.2D(CONC, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid,flux.x.up = fxup)$dC
    return (list(dCONC))
  })
}
#Se modifican las condiciones de frontera para crear un ambiente aislado

MFCDiffelec <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    nrow <- gridgeometry.grid2D$x.N
    ncol <- gridgeometry.grid2D$y.N
    CONC  <- matrix(nrow , ncol, data = y)
    Cxup=rep(0,ncol)
    fxdown=rep(0,ncol)
    fyup=rep(0,nrow)
    fydown=rep(0,nrow)
    ablxup=10
    dCONC <- ReacTran::tran.2D(CONC,C.x.up = Cxup, a.bl.x.up=ablxup,grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid,
                                flux.x.down = fxdown,  flux.y.up = fyup, flux.y.down = fydown)$dC
    return (list(dCONC))
  })
}

#Para el caso de los electrones, la condicion de frontera define el electrodo al poner 0 a su concentración

MFCDiffelec2 <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    nrow <- gridgeometry.grid2D$x.N
    ncol <- gridgeometry.grid2D$y.N
    CONC  <- matrix(nrow , ncol, data = y)
    Cxup=rep(0,ncol)
    #ablxup=10 , a.bl.x.up=ablxup
    dCONC <-ReacTran::tran.2D(CONC,C.x.up = Cxup, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC

    return (list(dCONC))
  })
}
#Para el caso de los electrones, la condicion de frontera define el electrodo al poner 0 a su concentración

MFCelec <- function (t, y, parms){
  # geometry values are in parms
  with (as.list(parms), {
    nrow <- gridgeometry.grid2D$x.N
    ncol <- gridgeometry.grid2D$y.N
    CONC  <- matrix(nrow , ncol, data = y)
    Cxup=rep(10,ncol)
    Cxdown=rep(0,ncol)
    #ablxup=10 , a.bl.x.up=ablxup
    dCONC <- ReacTran::tran.2D(CONC,C.x.up = Cxup,C.x.down=Cxdown, grid = gridgeometry.grid2D, D.grid = diffgeometry.Dgrid)$dC
    return (list(dCONC))
  })
}
Movimiento<-function(){
return(x,y)}

constrainMod<-function(object, reacts, lb, dryweight, tstep, scale, j, cutoff=1e-6){
  reacts = unique(reacts)
#print(reacts)
 Sac<--lb['EX_ac__40__e__41__']
Km=0.21206000322316076
qm=7.611210408114023
#print(Sac)
Dc<-0.00024 #cm
  lb = 0.8*lb[reacts]/(dryweight*tstep)*Dc**3/1000*1e12
#lb['EX_ac__40__e__41__']<--qm*Sac/(Sac+Km)
#print(lb['EX_ac__40__e__41__'])
  growth_limit <- (object@maxweight*1.5) - dryweight
  lobnd <- object@lbnd
lb['EX_ac__40__e__41__']<--qm*Sac/(Sac+Km)
#print('bicho')
#print(lobnd['EX_ac__40__e__41__'])
#print('ambiente')
#print(lb['EX_ac__40__e__41__'])
  upbnd <- object@ubnd
 #Comparar las cantidades que puede ingerir---hay que mejorarlo, porque solo funciona para cuando hay muy poquito sustrato
    #Comparemos mas bien flujos
#  if(dryweight<Inf){lobnd[reacts] <- object@lbnd[reacts]}#*dryweight*(dryweight/object@cellweight_mean)*tstep} #costrain according to flux definition: mmol/(gDW*hr)
  #lobnd[reacts] <- ifelse(lb<=lobnd[reacts], ifelse(lobnd[reacts]==0, lb, lobnd[reacts]), lb) #check if lower bounds in biological relevant range
  lobnd[reacts] <- ifelse(lb<=lobnd[reacts], lobnd[reacts], lb)

 #check if lower bounds in biological relevant range
  lobnd <-lobnd#*growth_limit/1000
  upbnd <-upbnd#*growth_limit/1000
  #if(j==1) browser()
 # if(length(object@kinetics) != 0){
 #   lobnd[names(object@kinetics)] <- unlist(lapply(names(object@kinetics), function(name){
 #     Km  <- (object@kinetics[[name]][["Km"]]*0.01*scale)*10^12 #scale mM to fmol/gridcell
 #     vmax <- object@kinetics[[name]][["vmax"]]*(dryweight/object@cellweight_mean)*tstep #scale to fmol/pgDW*hr
 #     s   <- -lb[name] # change sign to get concentrations
 #     lnew = -(vmax*s/(Km + s))*tstep
 #     if(abs(lnew)>s){if(-s>=lobnd[name]){lnew=-s}else{lnew=lobnd[name]}}
 #     return(lnew)
 #   }))
 # }
  if( object@limit_growth ){ # set upper bound for growth
    growth_limit <- (object@maxweight*1.5) - dryweight # 1.5 is factor to allow some room for the biomass to accumulate
    # if growth is too fast, then growth_limit could become negative (hight dryweight). Use low cutoff value so that dryweight can be decreased by cell division.
    upbnd[which(object@model@react_id == object@rbiomass)] <- ifelse( growth_limit>0, growth_limit, cutoff )
  }
  return(list(lobnd, upbnd))
}

