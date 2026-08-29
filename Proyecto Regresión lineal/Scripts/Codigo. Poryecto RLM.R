## ProyectoRLM 
library("car")
library("corrplot")
library(lmtest)
library(corrplot)
library(readr)
library(leaps)
library(segmented)


getwd()
list.files()

## Importar datos 
Turismo <- read.csv(file="Turismo.csv", header = TRUE)
Turismo

TUR <- Turismo$TUR
IPC <- Turismo$IPC
PH <- Turismo$PH
GH <- Turismo$GH
PC <- Turismo$PC
IGAE <- Turismo$IGAE
PROD <- Turismo$PROD
CL <- Turismo$CL
INF <- Turismo$INF

regresion_original<- lm(formula = TUR ~ IPC + PH + GH + PC + IGAE + PROD + CL + INF, data=Turismo)


###------- Modelo 1 -------###

### Criterio selección hacia delante y eliminación hacia atras
step <- step(regresion_original, direction="both")
regresion_1 <- lm(formula=TUR ~ IPC + PH + GH + PC + IGAE + CL +INF, data=Turismo)
summary(regresion_1)
anova(regresion_1)
correlacion <- cor(Turismo[, c("IPC", "PH", "GH", "PC", "IGAE", "INF")]) 
correlacion
corrplot(correlacion, method="color", col = colorRampPalette(c("skyblue", "white", "pink"))(200))
vif(regresion_1) ## No ya que solo hay tres variables predictoras que cumplen 1< VIF < 5
##1 <VIF < 5 correlacion modeara, >5 niveles criticos de multicolinealidad
##Varianza contanste 
bptest(regresion_1)

## box pirce independencia (no autocorrelacion)
residuales_regresion_1 <- regresion_1$residuals
Box.test(residuales_regresion_1, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_1, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.1)

## Normalidad
shapiro.test(residuales_regresion_1)



### Criterio R^2 ajustada
modelo0 <- regsubsets(TUR ~ IPC + PH + GH + PC + IGAE + PROD + CL + INF, data=Turismo, nvmax=8)
resultado <- summary(modelo0)
which.max(resultado$adjr2)
vif(regresion_original) ## No ya que solo hay dos variables predictoras que cumplen 1< VIF < 5
correlacion1 <- cor(Turismo[, c("IPC", "PH", "GH", "PC", "IGAE", "PROD", "CL", "INF")]) 
correlacion1
corrplot(correlacion1, method="color", col = colorRampPalette(c("skyblue", "white", "pink"))(200))
##Varianza contanste 
bptest(regresion_original)

## box pierce independencia (no autocorrelacion)
residuales_regresion_original <- regresion_original$residuals
Box.test(residuales_regresion_original, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_original, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_original)



### Por prueba y error 
regresionpe <- lm(formula=TUR~PH + PC + IGAE + PROD + CL + INF, data=Turismo)
step <- step(regresionpe, direction="both")

regresion_5 <- lm(formula = TUR ~ PH + IGAE + INF, data=Turismo)
summary(regresion_5)
anova(regresion_5)
vif(regresion_5) #Sí pasan la prueba del VIF 

correlacion2 <- cor(Turismo[, c("PH", "IGAE", "INF")]) 
correlacion2
corrplot(correlacion2, method="color", col = colorRampPalette(c("skyblue", "white", "pink"))(200))
# VARIANZA CONSTANTE
bptest(regresion_5)
#RESIDUALES
residuos <- regresion_5$residuals
# INDEPENDENCIA
Box.test(residuos, lag=1, type="Box-Pierce")
# NORMALIDAD
shapiro.test(residuos)
# MEDIA CERO
t.test(residuos, alternative="two.sided", mu=0, paired=FALSE, conf.level=0.95)

# Valores influyentes y atípicos
regresion_atip <- lm(formula = TUR ~ PH + IGAE + INF, data=Turismo)
regresion_atip
# Crear las matrices para modelo
X1 <- as.matrix(cbind(1, Turismo$PH, Turismo$IGAE, Turismo$INF))
X1
y1 <- as.matrix(Turismo$TUR, ncol=1)
y1
n1 <- nrow(X1)
n1
p1 <- ncol(X1)
p1

# Matriz sombrero H
H1 <- X1 %*% solve(t(X1) %*% X1) %*% t(X1)
H1
h_ii1 <- diag(H1)  # apalancamiento
h_ii1


# Usando el criterio hi>3*((p)/n), identifique las observaciones 
# influyentes.
criterio <- 3*((p)/n) # 0.1690141
criterio

# Valores de influencia
h_ii1

# Identificar observaciones influyentes
obs_influ <- which(h_ii1 > criterio)

# Mostrar resultados
tabla_obs_influ <- data.frame(
  Observación = 1:n,
  Apalancamiento = round(h_ii1, 4),
  Criterio = criterio,
  Influyente = ifelse(h_ii1 > criterio, "SÍ", "NO")
)

tabla_obs_influ

if(length(obs_influ) > 0) {
  Obs_sin_influ <- Turismo[-obs_influ, ] # Datos sin las observaciones 
  #influyentes
  
  regresion_final<- lm(formula = TUR ~ PH + IGAE + INF, data=Obs_sin_influ)
  
  cat("=== MODELO ORIGINAL ===\n")
  print(summary(regresion_5))
  
  cat("=== MODELO SIN OBSERVACIONES INFLUYENTES ===\n")
  print(summary(regresion_final))
  
} else {
  print(summary(regresion_5))
  cat("No hay observaciones influyentes según el criterio, por lo que los 
modelos son idénticos.\n")
}

# VARIANZA CONSTANTE
bptest(regresion_final)
#RESIDUALES
residuos_modfinal <- regresion_final$residuals
# INDEPENDENCIA
Box.test(residuos_modfinal, lag=1, type="Box-Pierce")
# NORMALIDAD
shapiro.test(residuos_modfinal)
# MEDIA CERO
t.test(residuos_modfinal, alternative="two.sided", mu=0, paired=FALSE, conf.level=0.95)

# Crear las matrices para modelo final
X0 <- as.matrix(cbind(1, Obs_sin_influ$PH, Obs_sin_influ$IGAE, Obs_sin_influ$INF))
X0
y0 <- as.matrix(Obs_sin_influ$TUR, ncol=1)
y0
n0 <- nrow(X0)
n0
p0 <- ncol(X0)
p0






# Nos quedamos con el de prueba y error, el que únicamente tiene 3 variables
# I. ¿Cómo cambia la variable respuesta al incrementar en una unidad cada variable regresora manteniendo fijas las demás? (En las unidades originales)
regresion_oficial <- lm(formula = TUR ~ PH + IGAE + INF, data=Obs_sin_influ)
regresion_oficial
summary(regresion_oficial)
anova(regresion_oficial)
# Intercepto = -8.98713. Cuando PH, IGAE E INF = 0, la actividad turística esperada sería de -8.98713
# Coeficiente de PH = 0.08326. Si se mantienen constantes IGAE e INF, y si el PH aumenta en 1 unidad, entonces la actividad turística aumenta 0.08326 unidades en promedio
# Coeficiente de IGAE = 1.01924. Si se mantienen constantes PH e INF, y si el IGAE aumenta en una unidad, entonces la actividad turística aumenta en 1.01924 unidades en promedio
# Coeficiente de INF = 0.84398. Si se mantienen constantes PH e IGAE, y si INF aumenta en una unidad, entonces la actividad turística aumenta en 0.84398 unidades en promedio 
# El coeficiente de determinación es de R^2 = 0.984, que significa que las variabls IGAE, INF y PH explican 98.4% de la variabilidad de la actividad turística



# II. ¿Tiene sentido utilizar una transformación Log-Log? De ser afirmativo, ¿Cómo se interpreta en este caso cada coeficiente de la RLM?

# En este caso, no sirve hacer transformación Log-Log, ya que la variable INF contiene valores negativos, lo cual nos arroja NaN's, que significa que no se puede aplicar este modelo sin arruinar el modelo
# Y como únicamente son 11 observaciones <0, entonces podemos simplemente quitarlas
regresion_log <- lm(log(TUR) ~ log(PH) + log(IGAE) + log(INF), data = Obs_sin_influ)
regresion_log
vif(regresion_log)
# Varianza constante
bptest(regresion_log)
# Residuales
residuos2 <- regresion_log$residuals
# Independencia
Box.test(residuos2, lag=1, type="Box-Pierce")
# Normalidad
shapiro.test(residuos2)
# Media cero
t.test(residuos2, alternative="two.sided", mu=0, paired=FALSE, conf.level=0.90)
# Cumple con 3/4 supuestos de un modelo de RLM 
summary(regresion_log)
# Como una regresión log-log todos los coeficientes se interpretan como elasticidades (cambio porcentual), no es necesario aplicar exponencial
# Intercepto = -0.487740. Cuando todas las variables independientes son iguales a 1 (Ln(1)=0), mi variable respuesta es -0.487740
# log(PH) = 0.097248. Un aumento del 1% de PH se asocia con un aumento del 0.097248% en la variable dependiente (TUR). Ceteris Paribus.
# log(IGAE) = 1.013011. Un aumento del 1% de IGAE está asociado con un aumento del 1.013011% en la variable TUR. Ceteris paribus.
# log(INF) = 0.001625. Un aumento del 1% de INF es asociado con un aumento del 0.001625% en la variable respuesta. Ceteris Paribus. 
### PREGUNTAR SI ESTO ES POSIBLE ###

# Comentario Shaiel: Dado que la variable de inflación presenta valores negativos, no es posible aplicar la 
# transformación logarítmica sobre dicha variable, ya que el logaritmo natural sólo está definido para valores 
# positivos. Por esta razón, no es posible estimar un modelo log-log completo. En su lugar, puede utilizarse un modelo 
# mixto donde únicamente se aplican logaritmos a las variables estrictamente positivas.





# III. ¿Tiene sentido utilizar una transformación Log-Lin? De ser afirmativo, ¿Cómo se interpreta en este caso cada coeficiente de la RLM?
regresion_loglin <- lm(log(TUR) ~ PH + IGAE + INF, data = Obs_sin_influ)
regresion_loglin
vif(regresion_loglin)
# Varianza constante
bptest(regresion_loglin)
# Residuos
residuos3 <- regresion_loglin$residuals
# Independencia
Box.test(residuos3, lag=1, type="Box-Pierce")
# Normalidad
shapiro.test(residuos3)
# Media cero
t.test(residuos3, alternative="two.sided", mu=0, paired=FALSE, conf.level=0.90)
# Cumple con 2/4 supuestos de RLM
summary(regresion_loglin)



# Intercepto = 3.4523780. Cuando las demás variables son 0, el logaritmo de TUR es 3.4523780. Es igual a decir que el valor de Y = exp{3.4523780} = 31.57539
# PH = 0.0007401. Si las demás variables se mantienen constantes y PH aumenta en una unidad, entonces log(T(UR) aumenta aproximadamente 0.07401%, o es equivalente a decir que TUR aumenta exp{0.0007401}-1=0.0007743096
# IGAE = 0.0109002. Si las demás variables se mantienen constantes e IGAE aumenta en una unidad, entonces log(TUR) aumenta aproximadamente 1.09002%, o es equivalente a decir que TUR aumenta exp{0.0109002}-1=0.01095982
# INF = 0.0093570. Si las demás variables se mantienen constantes e INF aumenta en una unidad, entonces log(TUR) aumenta aproximadamente 0.93570%, o es equivalente a decir que TUR aumenta exp{0.0093570}-1=0.009400914

## IV. Intervalo de predicción y respuesta media
## IVa. Intervalo de respuesta media

## Regresión original
# Crear el vector x0
x0 <- matrix(c(1,125,95,0.3), ncol=1)
x0
Turismo
# Obtener y ajustada
beta_gorro0 <- solve(t(X0) %*% X0) %*% t(X0) %*% y0
y_ajust0 <- as.numeric(t(x0) %*% beta_gorro0)
y_ajust0

# Residuales
e <- regresion_final$residuals
SSres <- as.numeric(t(e) %*% e)
SSres

varMCO <- SSres / (n0-p0)
varMCO

# Cuantil
alpha <- 0.05
t0 <- qt(alpha/2, df=n0-p0, lower.tail=FALSE)
t0

## Matriz x0'(X'X)x0
M0 <- as.numeric(t(x0) %*% solve(t(X0) %*% X0) %*% x0)
M0

# Intervalo de confianza para la respuesta media al 95% de confianza

lim_inf0 <- y_ajust0 - (t0 * sqrt(varMCO*M0))
lim_inf0

lim_sup0 <- y_ajust0 + (t0 * sqrt(varMCO*M0))
lim_sup0

# Con comando
x0.1 <- data.frame(
  PH = 125,
  IGAE = 95,
  INF = 0.3
)

resp_med1 <- predict(regresion_oficial, newdata = x0.1,
                     interval = "confidence", level = 0.95)
resp_med1


## Log-Lin

x0.2 <- data.frame(
  PH = 125,
  IGAE = 95,
  INF = 0.3
)

resp_med2 <- predict(regresion_loglin, newdata = x0.2, interval = "confidence", level = 0.95)
resp_med2

exp(resp_med2)

## IVb. Intervalo de predicción
# Intervalo de confianza para la predicción al 95% de confianza

## Regresión original
y_ajust0

lim_inf0.1 <- y_ajust0 - (t0 * sqrt(varMCO*(1+M0)))
lim_inf0.1

lim_sup0.1 <- y_ajust0 + (t0 * sqrt(varMCO*(1+M0)))
lim_sup0.1

# Con comando
x0.3 <- data.frame(
  PH = 125,
  IGAE = 95,
  INF = 0.3
)

predic1 <- predict(regresion_oficial, newdata = x0.3,
                   interval = "prediction", level = 0.95)
predic1

## Log-Lin

x0.4 <- data.frame(
  PH = 125,
  IGAE = 95,
  INF = 0.3
)

predic2 <- predict(regresion_loglin, newdata = x0.4, interval = "prediction", level = 0.95)
predic2

exp(predic2)


# V. ¿Cómo interpretamos los resultados de la RLM cuando las variables son estandarizadas?
TURstd <- (TUR-mean(TUR))/sd(TUR)
PHstd <- (PH-mean(PH))/sd(PH)
IGAEstd <- (IGAE-mean(IGAE))/sd(IGAE)
INFstd <- (INF-mean(INF))/sd(INF)

regresion_std <- lm(TURstd ~ 0 + PHstd + IGAEstd + INFstd, data = Obs_sin_influ)
regresion_std
summary(regresion_std)
# Varianza constante
bptest(regresion_std)
# Residuos
residuos4 <- regresion_std$residuals
# Independencia
Box.test(residuos4, lag=1, type="Box-Pierce")
# Normalidad
shapiro.test(residuos4)
# Media cero
t.test(residuos4, alternative="two.sided", mu=0, paired=FALSE, conf.level=0.90)
# Cumple con 2/4 supuestos de RLM
summary(regresion_std)