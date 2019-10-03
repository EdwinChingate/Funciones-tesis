closeAllConnections()
rm(list=ls())
source('funciones.R')

simM<- function (I='1',n=5, m=40, dt=1e-3, t=2e3, bf=1 ,res=1000,nMO=5)
{ #Todo el codigo se define como una función para facilitar el desarrollo de experimentos
#Descripción de parámetros: I es para facilitar el nombre del experimento realizado; n y m son el ancho y largo en número de celdas del espacio; dt es el tamaño de paso en la discretización del tiempo (h); t es el número de pasos a desarrollar; bf es el ancho del biofilm en número de celdas; nMO es el número de microorganismos que se especifican inicialmente.
	#se toma el tiempo en el que inicia la simulación para verificar cuanto tardó
	hora0<-Sys.time()
	horas<-gsub(' ','_',toString(hora0))
	ruta<-getwd()
	#setwd(ruta)
	#Con la fecha, se crea un directorio para guardar los resultados de cada experimento para poder diferenciarlos y analizarlos uno por uno   
print(n)
print(m)
print(dt)
print(t)
print(res)
print(nMO)

	guardar<-paste(ruta,'/resultadosM',sep='')
	dir.create(guardar)
	guardar<-paste(guardar,'/',horas,'-',I,sep='')
	dir.create(guardar)
	#Se carga el modelo del microorganismo
	GEO<-readSBMLmod('geo.xml')
GEO@react_id[533]<-'EX_ac__40__e__41__'
#load("GEO.RData")
	#Pf<-readSBMLmod('PF.xml')
	#Bq<-readSBMLmod('BQ.xml')
	#Cs<-readSBMLmod('CS.xml')
	#Ea<-readSBMLmod('EA.xml')
	#Se carga el medio en el que se desarrollaran los microorganismos / las concentraciones, difusividades y funciones asociadas a cada sustancia
	#data('Ec_core')
#Se carga el medio en el que se desarrollará el experimento
	medium <- read.csv(paste('MFCmed',I,'.csv',sep=''))
	geo<-Bac(GEO)#,chem='EX_elec')
	geo@speed<-0
#Geobacter tiene la capacidad de desplazarse en dirección al electrodo por el gradiente de elec.
	#eco<-Bac(Ec_core)
	#chem='EX_elec' le otorga a geo la capacidad de desplazarce preferiblemente en dirección al electrodo
	#geo@limit_growth=FALSE
	#geo@maxweight=536

	#pf<-Bac(Pf)
	#bq<-Bac(Bq)
	#cs<-Bac(Cs)
	#ea<-Bac(Ea)
	
	#Se crea el espacio donde interactuaran los microorganismos
	arena<-Arena(n=n,m=m,Lx=0.00024*n,Ly=0.00024*m,tstep=dt,stir=TRUE)

	arena<-addOrg(arena,geo,amount=nMO)
	#arena<-addOrg(arena,pf,amount=nMO)
	#arena<-addOrg(arena,bq,amount=nMO)
	#arena<-addOrg(arena,cs,amount=nMO)
	#arena<-addOrg(arena,ea,amount=nMO)
	arena@orgdat$y[1]=1
	arena@orgdat$y[2]=1
	arena@orgdat$y[3]=1
	arena@orgdat$y[4]=1
	arena@orgdat$y[5]=1
	arena@orgdat$x[1]=1
	arena@orgdat$x[2]=2
	arena@orgdat$x[3]=3
	arena@orgdat$x[4]=4
	arena@orgdat$x[5]=5

#	arena@orgdat$y[1]=m
	#Se carga el medio al espacio recien creado
	arena <- addSubs(arena, smax=medium$Concentracion, mediac=medium$Sustancia, difspeed=medium$Difusividad,unit='fmol/cell')
	vecb=matrix(nrow=length(arena@mediac),ncol=1,1) #La lleno con unos para sacarle el cuerpo a muchas cosas
	for(j in seq_along(arena@media)) 
	{#Se especifica la función que define el espacio para cada una de las sustancias
		arena@media[[j]]@pde<-toString(medium[j,4])
		vecb[j]<-medium[j,5]
	}
	arena<-addSubs(object=arena,smax=1,mediac='EX_elec',add=FALSE,difspeed=10,unit='fmol/cell')
	#Para evitar errores, se perturba la concentración de la "sustancia" EX_elec
	for (co in 1:m){
	arena@media[['EX_elec']]@diffmat[co,]=matrix(nrow=1,ncol=n,(m-co))
}


	#Se limita la movilidad de los electrones unicamente al biofilm, dado que es el unico con capacidad conductora
	Dbf<-arena@media[['EX_e(e)']]@difspeed
	#Aunque en realidad no importa la concentración de los electrones
	x.mid<-matrix(nrow=m,ncol=n,0)
	x.int<-matrix(nrow=m+1,ncol=n,0)
	y.mid<-matrix(nrow=m,ncol=n,0)
	y.int<-matrix(nrow=m,ncol=n+1,0)
	Vad<-list(x.mid = x.mid, y.mid = y.mid, x.int = x.int, y.int = y.int)
	x.int[1:bf,1:(n)]<-Dbf
	y.int[1:bf,1:(n+1)]<-Dbf
	x.mid[1:bf,1:(n)]<-Dbf
	y.mid[1:bf,1:(n)]<-Dbf
	Difelec<-list(x.mid = x.mid, y.mid = y.mid, x.int = x.int, y.int = y.int)
	Difgelec<-list(Dgrid=Difelec, Vgrid=Vad)
	arena@media[['EX_e(e)']]@diffgeometry<-Difgelec

	#Potencial en el espacio, que regulará el metabolismo de geobacter, de tal manera que solo hay en el biofilm, los microorganismos solo pueden transferir electrones en el biofilm
	E0=20
	potspace<-matrix(nrow=m,ncol=n,-1000)
	potspace[1:bf,1:n]<-matrix(nrow=bf,ncol=n,E0)

	#Aqui se corre la simulación

	print(arena@orgdat)
	eva<-simEnvMod(arena,time=t,potspace=potspace,cutoff=1e-20,guardar=guardar,res=res,vecb=vecb,E0=E0,bf=bf,sec_obj='mtf')#,diff_par=TRUE)

	hora1=Sys.time()

	h=hora1-hora0
	print(h)
	#En esta sección se crean copias de los archivos que se usaron para guardar todo en una misma carpeta con la fecha y hora, y poder llevar registros de los cambios en el proyecto
	write(paste(horas,'\n',toString(hora1),'\n',toString(h),sep=''),paste(guardar,'/tiempo.txt',sep=''))
	write(ruta,paste(guardar,'/ruta.txt',sep=''))
	simulacion<-'/simM.R'
	medio<-paste('/MFCmedM',I,'.csv',sep='')
	MO<-'/GeoM.xml'
	funciones<-'/funciones.R'
	result<-'/resultadom.R'
	figuras<-'/figuras.R'

	copiar<-c(simulacion,medio,funciones,result,figuras)#,MO) copiar el microorganismo solo si es necesario
#Guardar una copia de los archivos que se usaron para construir la simulación 

	for (x in copiar)
	{
		from=paste(ruta,x,sep='')
		to=paste(guardar,x,sep='')
		if (x==result){
		file.rename(from,to)}else{
		file.copy(from,to)}
	}

	##Crear el archivo para mostrar los perfiles temporales
#	new<-paste('cd', guardar,'\n','R CMD BATCH figuras.R result.R')
#	write(new,paste(ruta,'/fig.sh',sep=''))
	}
