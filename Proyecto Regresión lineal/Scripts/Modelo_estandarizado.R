## ProyectoRLM 
library(car)
library(corrplot)
library(lmtest)
library(corrplot)
library(readr)
library(leaps)
library(segmented)

## Importar datos 
Turismo <- read.csv(file="Turismo1.csv", header = TRUE)
Turismo

TUR <- Turismo$TUR
IPC <- Turismo$IPC
PH <- Turismo$PH
GH <- Turismo$GH
PC <- Turismo$PC
PROD <- Turismo$PROD
CL <- Turismo$CL
INF <- Turismo$INF

#----
regresion_original<- lm(formula = TUR ~ IPC + PH + GH + PC + PROD + CL + INF, data=Turismo)
summary(regresion_original)
anova(regresion_original)
correlacion <- cor(Turismo[, c("IPC", "PH", "GH", "PC", "PROD", "CL","INF")]) 
correlacion
corrplot(correlacion)

## Multicolinealidad
vif(regresion_original)

##Varianza contante 
bptest(regresion_original)

## box pierce independencia (no autocorrelacion)
residuales_regresion_original <- regresion_original$residuals
Box.test(residuales_regresion_original, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_original, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_original)

## Prueba de Hipótesis Individual
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL con α =", alpha, "===\n")
for(i in 1:length(coef(regresion_original))) {
  # Obtener el coeficiente y su error estándar
  beta1 <- coef(regresion_original)[i]
  error_std1 <- summary(regresion_original)$coefficients[i, 2]
  
  t_calculado1 <- beta1 / error_std1 # Estadístico de prueba
  
  df1 <- df.residual(regresion_original) # Grados de libertad
  
  p_valor1 <- 2 * pt(abs(t_calculado1), df1, lower.tail = FALSE) # p-value
  
  cat(names(coef(regresion_original))[i], "con p =", format(p_valor1, scientific = TRUE, digits = 4), "es", 
      ifelse(p_valor1 < alpha, "Significativo", "No significativo"), "\n")
}

## Prueba de Hipótesis conjunta
# H0: B1 = ... = Bi = 0
# HA: Existe al menos un i tal que Bi != 0, con i=1,...,k
# Crear las matrices para modelo
X1 <- as.matrix(cbind(1,Turismo$IPC, Turismo$PH, Turismo$GH, Turismo$PC, Turismo$PROD, Turismo$CL, Turismo$INF))
y1 <- as.matrix(Turismo$TUR, ncol=1)
n3 <- nrow(X1)
k3 <- ncol(X1)-1

# Matriz sombrero H
H1 <- X1 %*% solve(t(X1) %*% X1) %*% t(X1)
H1

beta_m1 <- solve(t(X1) %*% X1) %*% t(X1) %*% y1
beta_m1

# Obtener matriz y
y_matriz1 <- X1 %*% beta_m1
y_matriz1

# Vector de residuos
e1 <- y1 - H1 %*% y1
e1

# Calcular Sumas de Cuadrados
SSreg1 <- as.numeric((t(beta_m1) %*% t(X1) %*% y_matriz1) - ((sum(y1)^2)/n3))
SSreg1
SSres1 <- as.numeric(t(e1) %*% e1)
SSres1

# Calcular estadístico de prueba F
F0_1 <- (SSreg1/(k3))/(SSres1/(n3-k3-1))
F0_1 # 113.3084

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_1, df1=k3, df2=n3-k3-1, lower.tail=FALSE) 
# 3.207159e-33 < 0.05 => Se rechaza H0


## Intervalos de respuesta media y de predicción
# Intervalos de respuesta media
# Crear el vector x0
x0 <- matrix(c(1,75,100,105,95,106,120,0.2), ncol=1)
x0
Turismo

# Obtener y ajustada
y_ajust1 <- as.numeric(t(x0) %*% beta_m1)
y_ajust1

# Calcular varianza por MCO
SSres1
varMCO1 <- SSres1/(n3-k3-1)
varMCO1

# Cuantil
t1 <- qt(alpha/2, df=n3-k3-1, lower.tail=FALSE)
t1

## Matriz x0'(X'X)x0
M1 <- as.numeric(t(x0) %*% solve(t(X1) %*% X1) %*% x0)
M1

# Intervalo de confianza para la respuesta media al 95% de confianza

lim_inf1 <- y_ajust1 - (t1 * sqrt(varMCO1*M1))
lim_inf1

lim_sup1 <- y_ajust1 + (t1 * sqrt(varMCO1*M1))
lim_sup1


# Con comando
x0.1 <- data.frame(
  IPC = 75,
  PH = 100,
  GH = 105,
  PC = 95,
  PROD = 106,
  CL = 120,
  INF = 0.2
)

resp_med1 <- predict(regresion_original, newdata = x0.1,
                     interval = "confidence", level = 0.95)
resp_med1

# Intervalo de predicción con comando
int_pred1 <- predict(regresion_original, newdata = x0.1,
                     interval = "predict", level = 0.95)
int_pred1


#Akaike----
step(regresion_original, direction="both")
regresion_step <- lm(formula=TUR ~ IPC + PH + GH + PC + CL, data=Turismo)
summary(regresion_step)
anova(regresion_step)
correlacion <- cor(Turismo[, c("IPC", "PH", "GH", "PC", "CL")]) 
correlacion
corrplot(correlacion, method="color", col = colorRampPalette(c("skyblue", "white", "pink"))(200))
## Multicolinealidad
vif(regresion_step)

##Varianza contante 
bptest(regresion_step)

## box pierce independencia (no autocorrelacion)
residuales_regresion_step <- regresion_step$residuals
Box.test(residuales_regresion_step, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_step, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_step)

# valores atípicos
n <- length(IPC)
k <- 5
studentizado_step <- rstudent(regresion_step)
t_alpha <- qt(0.05/2, n-k-1, lower.tail = FALSE)
t_alpha

atipicos <- which(abs(studentizado_step)>t_alpha)
atipicos

# valores influyentes
Dist_cook <- cooks.distance(regresion_step)
which(Dist_cook>1)

hat_values <- hatvalues(regresion_step)
influyentes0 <- which(hat_values>3*(k+1)/n)
influyentes0

atipicos_influyentes <- union(atipicos,influyentes0)
atipicos_influyentes
Turismo_limpio_step <- Turismo[-atipicos_influyentes,]
Turismo_limpio_step
modelo_limpio_AIC <- lm(TUR ~ IPC + PH + GH + PC + CL, data=Turismo_limpio_step)
summary(modelo_limpio_AIC)
residuales_regresion_limpia_AIC <- modelo_limpio_AIC$residuals

#Multicolinealidad sin valores influyentes
vif(modelo_limpio_AIC)

##Varianza contante 
bptest(modelo_limpio_AIC)

## box pierce independencia (no autocorrelacion)
Box.test(residuales_regresion_limpia_AIC, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_limpia_AIC, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_limpia_AIC)

#MPE
y_Akaike <- Turismo_limpio_step$TUR
PM_Akaike <- 100/length(y_Akaike) * sum(residuales_regresion_limpia_AIC / y_Akaike)
PM_Akaike

##############################################
#BIC----
modelo_regsubsets <- regsubsets(TUR ~ IPC + PH + GH + PC + PROD + CL + INF, data=Turismo, nvmax=7)
resultado_regsubsets <- summary(modelo_regsubsets)
resultado_regsubsets

#Criterio de información Bayesiana
which.min(resultado_regsubsets$bic)
resultado_regsubsets$which[3,]

modelo_BIC <- lm(TUR ~ IPC + GH + PC, data=Turismo)
summary(modelo_BIC)

## Prueba de Hipótesis Individual
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL con α =", alpha, "===\n")
for(i in 1:length(coef(modelo_BIC))) {
  # Obtener el coeficiente y su error estándar
  beta2 <- coef(modelo_BIC)[i]
  error_std2 <- summary(modelo_BIC)$coefficients[i, 2]
  
  t_calculado2 <- beta2 / error_std2 # Estadístico de prueba
  
  df2 <- df.residual(modelo_BIC) # Grados de libertad
  
  p_valor2 <- 2 * pt(abs(t_calculado2), df2, lower.tail = FALSE) # p-value
  
  cat(names(coef(modelo_BIC))[i], "con p =", format(p_valor2, scientific = TRUE, digits = 4), "es", 
      ifelse(p_valor2 < alpha, "Significativo", "No significativo"), "\n")
}

## Prueba de Hipótesis conjunta
# H0: B1 = ... = Bi = 0
# HA: Existe al menos un i tal que Bi != 0, con i=1,...,k
# Crear las matrices para modelo
X2 <- as.matrix(cbind(1,Turismo$IPC, Turismo$GH, Turismo$PC))
y2 <- as.matrix(Turismo$TUR, ncol=1)
n4 <- nrow(X2)
k4 <- ncol(X2)-1

# Matriz sombrero H
H2 <- X2 %*% solve(t(X2) %*% X2) %*% t(X2)
H2

beta_m2 <- solve(t(X2) %*% X2) %*% t(X2) %*% y2
beta_m2

# Obtener matriz y
y_matriz2 <- X2 %*% beta_m2
y_matriz2

# Vector de residuos
e2 <- y2 - H2 %*% y2
e2

# Calcular Sumas de Cuadrados
SSreg2 <- as.numeric((t(beta_m2) %*% t(X2) %*% y_matriz2) - ((sum(y2)^2)/n4))
SSreg2
SSres2 <- as.numeric(t(e2) %*% e2)
SSres2

# Calcular estadístico de prueba F
F0_2 <- (SSreg2/(k4))/(SSres2/(n4-k4-1))
F0_2 # 240.2578

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_2, df1=k4, df2=n4-k4-1, lower.tail=FALSE) 
# 8.813258e-36 < 0.05 => Se rechaza H0


## Intervalos de respuesta media y de predicción
# Intervalos de respuesta media con comando
x0.2 <- data.frame(
  IPC = 75,
  GH = 105,
  PC = 95
)

resp_med2 <- predict(modelo_BIC, newdata = x0.2,
                     interval = "confidence", level = 0.95)
resp_med2

# Intervalo de predicción con comando
int_pred2 <- predict(modelo_BIC, newdata = x0.2,
                     interval = "predict", level = 0.95)
int_pred2

# valores atípicos
n2 <- length(IPC)
k2 <- 3
studentizado_BIC <- rstudent(modelo_BIC)
t_alpha2 <- qt(0.05/2, n2-k2-1, lower.tail = FALSE)
t_alpha2

atipicos2 <- which(abs(studentizado_BIC)>t_alpha2)
atipicos2

# valores influyentes
Dist_cook_BIC <- cooks.distance(modelo_BIC)
which(Dist_cook_BIC>1)

hat_values_1 <- hatvalues(modelo_BIC)
influyentes1 <- which(hat_values_1>3*(k2+1)/n2)
influyentes1

atipicos_influyentes_1 <- union(atipicos2,influyentes1)
atipicos_influyentes_1
Turismo_limpio_BIC <- Turismo[-atipicos_influyentes_1,]
Turismo_limpio_BIC

# Verificación de supuestos
modelo_limpio_BIC <- lm(TUR ~ IPC + GH + PC, data=Turismo_limpio_BIC)
summary(modelo_limpio_BIC)
residuales_regresion_limpia_BIC <- modelo_limpio_BIC$residuals

#Multicolinealidad sin valores influyentes
vif(modelo_limpio_BIC)

##Varianza constante 
bptest(modelo_limpio_BIC)

## box pierce independencia (no autocorrelacion)
Box.test(residuales_regresion_limpia_BIC, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_limpia_BIC, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_limpia_BIC)

## Prueba de Hipótesis Individual
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL con α =", alpha, "===\n")
for(i in 1:length(coef(modelo_limpio_BIC))) {
  # Obtener el coeficiente y su error estándar
  beta3 <- coef(modelo_limpio_BIC)[i]
  error_std3 <- summary(modelo_limpio_BIC)$coefficients[i, 2]
  
  t_calculado3 <- beta3 / error_std3 # Estadístico de prueba
  
  df3 <- df.residual(modelo_limpio_BIC) # Grados de libertad
  
  p_valor3 <- 2 * pt(abs(t_calculado3), df3, lower.tail = FALSE) # p-value
  
  cat(names(coef(modelo_limpio_BIC))[i], "con p =", format(p_valor3, scientific = TRUE, digits = 4), "es", 
      ifelse(p_valor3 < alpha, "Significativo", "No significativo"), "\n")
}

## Prueba de Hipótesis conjunta
# H0: B1 = ... = Bi = 0
# HA: Existe al menos un i tal que Bi != 0, con i=1,...,k
# Crear las matrices para modelo
X3 <- as.matrix(cbind(1,Turismo_limpio_BIC$IPC, Turismo_limpio_BIC$GH, Turismo_limpio_BIC$PC))
y3 <- as.matrix(Turismo_limpio_BIC$TUR, ncol=1)
n5 <- nrow(X3)
k5 <- ncol(X3)-1

# Matriz sombrero H
H3 <- X3 %*% solve(t(X3) %*% X3) %*% t(X3)
H3

beta_m3 <- solve(t(X3) %*% X3) %*% t(X3) %*% y3
beta_m3

# Obtener matriz y
y_matriz3 <- X3 %*% beta_m3
y_matriz3

# Vector de residuos
e3 <- y3 - H3 %*% y3
e3

# Calcular Sumas de Cuadrados
SSreg3 <- as.numeric((t(beta_m3) %*% t(X3) %*% y_matriz3) - ((sum(y3)^2)/n5))
SSreg3
SSres3 <- as.numeric(t(e3) %*% e3)
SSres3

# Calcular estadístico de prueba F
F0_3 <- (SSreg3/(k5))/(SSres3/(n5-k5-1))
F0_3 # 302.0185

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_3, df1=k5, df2=n5-k5-1, lower.tail=FALSE) 
# 6.1753e-37 < 0.05 => Se rechaza H0


## Intervalos de respuesta media y de predicción
# Intervalos de respuesta media con comando
x0.3 <- data.frame(
  IPC = 75,
  GH = 105,
  PC = 95
)

resp_med3 <- predict(modelo_limpio_BIC, newdata = x0.3,
                     interval = "confidence", level = 0.95)
resp_med3

# Intervalo de predicción con comando
int_pred3 <- predict(modelo_limpio_BIC, newdata = x0.3,
                     interval = "predict", level = 0.95)
int_pred3


#MPE
y_cib <- Turismo_limpio_BIC$TUR
PM_cib <- 100/length(y_cib) * sum(residuales_regresion_limpia_BIC / y_cib)
PM_cib

############### MEJOR MPE
#MPE Akaike = -0.04326615
#MPE BIC = -0.04390477
#Akaike está más cercano a 0

#MODELO GANADOR AIC