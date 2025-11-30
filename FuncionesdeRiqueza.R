
################################################################
### Curvas de Acumulación de Especies Observadas y Estimadas ###
################################################################


### Estimador de la Riqueza Total Chao 1 ###
# Necesaria para calcular la varianza no condicional de rarefacción

chao1 <- function(data) {
	
	data <- data[data > 0]

	data <- as.numeric(data)
	
	Sobs <- length(data)
	
	singlet <- length(data[which(data==1)])
	doublet <- length(data[which(data==2)])
	
	if(doublet == 0) {
	
	Sest <- Sobs + singlet*(singlet-1)/2
		
	} else
	
	Sest <- Sobs + singlet^2/(2*doublet)
	
	return(Sest)
}




### Riqueza Esperada en Muestras Enrarecidas (Rarefaction) ###

# Calcula el numero esperado de especies en una muestra de menos individuos. (Hurlbert 1971, Brewer & Williamson 1994). La formula de Brewer & Williamson (1994) es más robusta a los problemas numéricos que pueden surgir en el calculo de combinaciones requerida por la formula de Hurlbert (1971) y es entonces preferida. El método de Coleman (1981) es casi equivalente y es aun más sencillo de calcular. La varianze es la "varianza no condicional" estimada con la ecuación de Colwell et al. (2012).

# Colwell, R. K., Chao, A., Gotelli, N. J., Lin, S. Y., Mao, C. X., Chazdon, R. L., & Longino, J. T. (2012). Models and estimators linking individual-based and sample-based rarefaction, extrapolation and comparison of assemblages. Journal of Plant Ecology 5(1):3-21.

Rarefac <- function(data, n=sum(data), method="Brewer & Williamson 1994", varianza=FALSE, Stot=NULL) {
	
	# Eliminar ceros, si están presentes
	
	data <- data[data > 0]

	# Abundancia total
	N <- sum(data)
	
	if(n > N) stop("El número de individuos evaluados (n = ",n, ") debe ser menor que el número total de individuos (N = ", N, ").")
	
	# Numero de especies
	Sobs <- length(data)

	if(method=="Hurlbert 1971") {

	B <- choose(N, n)

	comb <- numeric()
	for(i in 1:length(data)) { comb[i] <- choose(N-data[i], n) }
	
	Comb <- sum(comb, na.rm =TRUE)

	# Calculo final
	Rare <- Sobs - Comb/B
	
	}
	
	if(method=="Brewer & Williamson 1994") {
		
		SUM <- numeric()
		
		for(i in 1:Sobs) {
			
			PROD <- numeric()

			for(j in 0:data[i]) { PROD[j] <- 1-n/(N-j) }
			
			SUM[i] <- prod(PROD, na.rm =TRUE)
		}
		
		Rare <- Sobs - sum(SUM)
	}
	
	if(method=="Coleman 1981") { Rare <- Sobs - sum((1-n/N)^data) }

	if(varianza==TRUE) {
		#  unconditional variance of Colwell et al. (2012)
		# compute in logarithmic scale to avoid numerical overflow
		B <- lchoose(N, n)

		comb2 <- numeric()
		
		for(i in 1:length(data)) {	
			comb2[i] <- (1 - exp(lchoose(N-data[i], n) - B) )^2 # substract logarithms
		}
		
		Comb2 <- sum(comb2, na.rm =TRUE)
		
		if(is.null(Stot))	{ Stot <- chao1(data) }
		
		var <- Comb2 - (Rare^2)/Stot
		
		Rare <- cbind(Rare, sqrt(var))
	}

	return(Rare)	
}




# Genera el numero esperado de especies de n=1 hasta n=N, para graficar una "curva de enrarecimiento (rarefacción)".
# Opcionalmente la riqueza esperada se calcula sólo para algunos valores de número de individuos, sequencia que se puede definir por ejemplo con seq(1,N, length.out=100)
# Notar que la riqueza total necesaria para calcular la varianza del número esperado de especie se calcula una vez para la muestra completa, en lugar de hacer estimaciones con muestras menores.

RareCurve <- function(data, x=NULL, varianza=FALSE, Stot=NULL, ...) {
	
	if(is.null(x)) { x <- 1:sum(data) }

	if(varianza==FALSE) {
		Sexp <- numeric()
		for(i in 1:length(x)) { Sexp[i] <- Rarefac(data, n=x[i], varianza=FALSE)}
		res <- cbind(x, Sexp)
		res <- rbind(c(0,0), res)
		colnames(res) <- c("n", "Sexp")
	}
	if(varianza==TRUE) {
		
		if(is.null(Stot))	{ Stot <- chao1(data) }
		
		Sexp <- matrix(nrow=length(x), ncol=2)
		for(i in 1:length(x)) { Sexp[i,] <- Rarefac(data, n=x[i], varianza=TRUE, Stot=Stot)}
		res <- cbind(x, Sexp)
		res <- rbind(c(0,0,0), res)
		colnames(res) <- c("n", "Sexp", "Ssd")
		}
	return(as.data.frame(res))
}



plot.RareCI <- function(rare, type=1, ...) {
	
	rare <- as.data.frame(rare)
	
	final <- nrow(rare)
	
	if(type==1) {
		polygon(x=c(rare$n, rev(rare$n)), y=c(rare$Sexp+2*rare$Ssd, rev(rare$Sexp-2*rare$Ssd)),...)
	}

	if(type==2) { 
		
		arrows(rare$n[final], rare$Sexp[final]+2*rare$Ssd[final], rare$n[final], rare$Sexp[final]-2*rare$Ssd[final], code=3, length=0.05, angle=90, col=color, ...)
	
		points(rare$n[final], rare$Sexp[final], pch=16, cex=0.75, col=color, ...)
	}
}



### Acumulación de Especies Estimadas ###

# Esta función toma submuestras de individuos y estima la riqueza total utilizando la función SpecAbunAce del paquete SpadeR.

# Uso

# data		vector con abundancias de cada especie (como en ChaoSpecies)
# valores	valores de número total de individuos para los que se quiere evaluar las estimaciones.
# replicas	cuantas replicas de las estimaciones para cada valor de abundancia total

Srarefac <- function(data, valores=seq(from=max(c(12, round(sum(data)/20))), to=sum(data), length.out=20), replicas=10) {
		
	# Eliminar ceros, si están presentes
	
	data <- data[data > 0]	
	
	# Abundancia total
	N <- sum(data)
		
	# Expand data to a factor of individuals using numbers to identiy species
	
	individ <- factor()
	for(i in 1:length(data)) {
	individ <- append(individ, rep.int(i, times=data[i]))	
	}
	
	RES <- matrix(ncol=4, nrow=length(valores))
	
	for(i in 1:length(valores)) {

		resj <- matrix(ncol=4, nrow=replicas)
		
		for(j in 1:replicas) {
		#genera una pseudoréplica de i individuos y lista la abundancia de cada especie
		PR <- table(sample(individ, size=valores[i]))
		
		#Estimar Parametros
		resj[j,] <- SpadeR:::SpecAbunAce(PR, k=10, conf=0.95)		
		}		
		
		# Resumir resultados de réplicas y guardar
		
		RES[i,] <- apply(resj, 2, median, na.rm=TRUE)
	}
		
	RES <- cbind(valores, RES)

	colnames(RES) <- c("Individuos", "Riqueza", "Dev.Est.", "95%I", "95%S")
	
	return(RES)

}



### Gráfico de Evaluación de Suficiencia de Muestreo para Estimación de Riqueza Total de Especies ###

# data		vector con abundancias de cada especie (como en ChaoSpecies)
# valores	valores de número total de individuos para los que se quiere evaluar las estimaciones.
# replicas	cuantas replicas de las estimaciones para cada valor de abundancia total

Acumulacion.plot <- function(data, Sest=NULL, valores=seq(from=max(c(15, round(sum(data)/10))), to=sum(data), length.out=10), replicas=10, anotaciones=TRUE, pos=4, lineas=TRUE, ...) {
	
	if(is.null(Sest)) { Sest <- Srarefac(data, valores, replicas) }
	Srare <- RareCurve(data, varianza=FALSE)
	
	N  <- sum(data)
	SobsN <- Srare[dim(Srare)[1],"Sexp"]
	SestN <- Sest[dim(Sest)[1], "Riqueza"]
	
	plot(c(Srare[,"n"],Sest[,"Individuos"]), c(Srare[,"Sexp"], Sest[,"95%S"]), type='n', lwd=2, col="red", yaxs="i", xaxs='i', ylim=c(0, max(Sest[,"95%S"])*1.03), xlab="Individuos Muestreados", ylab="Número de Especies", las=1)
	
	lines(Srare, lwd=1.5)
	polygon(c(Sest[,"Individuos"], rev(Sest[,"Individuos"])), c(Sest[,"95%I"], rev(Sest[,"95%S"])), border=NA, col=rgb(1,0,0,0.3))

	lines(Sest[,"Individuos"],Sest[,"95%I"], col="red")
	lines(Sest[,"Individuos"], Sest[,"95%S"], col="red")
	lines(Sest[,"Individuos"],Sest[,"Riqueza"], col="red", lwd=2)

	if(anotaciones) {
	text(N, SestN, labels=round(SestN), col="red", pos=pos, xpd=TRUE)
	text(N, SobsN, labels=round(SobsN), pos=pos, xpd=TRUE)
	}

	if(lineas) {
	abline(h=SestN, lty="dashed", col="red")
	abline(h=SobsN, lty="dashed")
	}

}
