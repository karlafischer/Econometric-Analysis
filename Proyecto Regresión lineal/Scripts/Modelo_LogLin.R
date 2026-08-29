## ProyectoRLM 
library(car)
library(corrplot)
library(lmtest)
library(corrplot)
library(readr)
library(leaps)
library(segmented)
library(glmnet)

## Importar datos 
Turismo <- read.csv(file="Turismo1.csv", header = TRUE)
Turismo

TUR_log <- log(Turismo$TUR)
IPC <- Turismo$IPC
PH <- Turismo$PH
GH <- Turismo$GH
PC <- Turismo$PC
PROD <- Turismo$PROD
CL <- Turismo$CL
INF <- Turismo$INF

#----
regresion_loglin <- lm(TUR_log ~ IPC + PH + GH + PC + PROD + CL + INF, data=Turismo)
summary(regresion_loglin)
#Verificación de supuestos
#Multicolinealidad
vif(regresion_loglin)

##Varianza contante 
bptest(regresion_loglin)

## box pierce independencia (no autocorrelacion)
residuales_regresion_loglin <- regresion_loglin$residuals
Box.test(residuales_regresion_loglin, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_loglin, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.05)

## Normalidad
shapiro.test(residuales_regresion_loglin)

## Prueba de Hipótesis Individual
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL con α =", alpha, "===\n")
for(i in 1:length(coef(regresion_loglin))) {
  # Obtener el coeficiente y su error estándar
  beta1 <- coef(regresion_loglin)[i]
  error_std1 <- summary(regresion_loglin)$coefficients[i, 2]
  
  t_calculado1 <- beta1 / error_std1 # Estadístico de prueba
  
  df1 <- df.residual(regresion_loglin) # Grados de libertad
  
  p_valor1 <- 2 * pt(abs(t_calculado1), df1, lower.tail = FALSE) # p-value
  
  cat(names(coef(regresion_loglin))[i], "con p =", format(p_valor1, scientific = TRUE, digits = 4), "es", 
      ifelse(p_valor1 < alpha, "Significativo", "No significativo"), "\n")
}

## Prueba de Hipótesis conjunta
# H0: B1 = ... = Bi = 0
# HA: Existe al menos un i tal que Bi != 0, con i=1,...,k
# Crear las matrices para modelo
X1 <- as.matrix(cbind(1,IPC, PH, GH, PC, PROD, CL,INF))
y1 <- as.matrix(TUR_log, ncol=1)
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

anova(regresion_loglin)
# Calcular estadístico de prueba F
F0_1 <- (SSreg1/(k3))/(SSres1/(n3-k3-1))
F0_1 # 99.06794

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_1, df1=k3, df2=n3-k3-1, lower.tail=FALSE) 
# 1.543485e-31 < 0.05 => Se rechaza H0


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

resp_med1 <- predict(regresion_loglin, newdata = x0.1,
                     interval = "confidence", level = 0.95)
resp_med1

# Intervalo de predicción con comando
int_pred1 <- predict(regresion_loglin, newdata = x0.1,
                     interval = "predict", level = 0.95)
int_pred1

#Akaike----
step(regresion_loglin, direction="both")
regresion_step_loglin <- lm(log(TUR) ~ IPC + PH + GH + PC + INF, data=Turismo)
summary(regresion_step_loglin)
anova(regresion_step_loglin)
## Multicolinealidad
vif(regresion_step_loglin)

##Varianza contante 
bptest(regresion_step_loglin)

## box pierce independencia (no autocorrelacion)
residuales_regresion_step_loglin <- regresion_step_loglin$residuals
Box.test(residuales_regresion_step_loglin, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_step_loglin, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.05)

## Normalidad
shapiro.test(residuales_regresion_step_loglin)

# valores atípicos
n <- length(IPC)
k <- 5
studentizado_step_loglin <- rstudent(regresion_step_loglin)
t_alpha <- qt(0.05/2, n-k-1, lower.tail = FALSE)

atipicos <- which(abs(studentizado_step_loglin)>t_alpha)
atipicos

# valores influyentes
Dist_cook_step_loglin <- cooks.distance(regresion_step_loglin)
which(Dist_cook_step_loglin>1)

hat_values <- hatvalues(regresion_step_loglin)
influyentes <- which(hat_values>3*(k+1)/n)
influyentes

atipicos_influyentes <- union(atipicos,influyentes)
atipicos_influyentes
Turismo_limpio_step_loglin <- Turismo[-atipicos_influyentes,]

# Verificación de supuestos
modelo_limpio_step <- lm(log(TUR) ~ IPC + PH + GH + PC + INF, data=Turismo_limpio_step_loglin)
summary(modelo_limpio_step)
residuales_regresion_limpia_step <- modelo_limpio_step$residuals

#Multicolinealidad sin valores influyentes
vif(modelo_limpio_step)

##Varianza constante 
bptest(modelo_limpio_step)

## box pierce independencia (no autocorrelacion)
Box.test(residuales_regresion_limpia_step, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_regresion_limpia_step, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_regresion_limpia_step)

#MPE
y_Akaike <- Turismo_limpio_step$TUR
PM_Akaike <- 100/length(y_Akaike) * sum(residuales_regresion_limpia_step / y_Akaike)
PM_Akaike

#BIC----
modelo_regsubsets <- regsubsets(TUR_log ~ IPC + PH + GH + PC + PROD + CL + INF, data=Turismo, nvmax=7)
resultado_regsubsets <- summary(modelo_regsubsets)
resultado_regsubsets

#Criterio de información Bayesiana
which.min(resultado_regsubsets$bic)
resultado_regsubsets$which[4,]

modelo_BIC <- lm(TUR_log ~ IPC + PH + GH + PC, data=Turismo)
summary(modelo_BIC)
anova(modelo_BIC)
residuales_BIC <- modelo_BIC$residuals

#Multicolinealidad
vif(modelo_BIC)

##Varianza constante 
bptest(modelo_BIC)

## box pierce independencia (no autocorrelacion)
Box.test(residuales_BIC, lag=5, type="Box-Pierce")

##Media cero 
t.test(residuales_BIC, alternative= "two.sided", mu=0, paired= FALSE, conf.level=0.95)

## Normalidad
shapiro.test(residuales_BIC)

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
X2 <- as.matrix(cbind(1,IPC, PH, GH, PC))
y2 <- as.matrix(TUR_log, ncol=1)
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
summary(modelo_BIC)
# Calcular estadístico de prueba F
F0_2 <- (SSreg2/(k4))/(SSres2/(n4-k4-1))
F0_2 # 167.6521

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_2, df1=k4, df2=n4-k4-1, lower.tail=FALSE) 
# 8.281416e-34 < 0.05 => Se rechaza H0


## Intervalos de respuesta media y de predicción
# Intervalos de respuesta media con comando
x0.2 <- data.frame(
  IPC = 75,
  PH = 100,
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
k2 <- 4
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

modelo_limpio_BIC <- lm(log(TUR) ~ IPC + PH + GH + PC, data=Turismo_limpio_BIC)
summary(modelo_limpio_BIC)
anova(modelo_limpio_BIC)
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
TUR_log_limpio <- log(Turismo_limpio_BIC$TUR)
X3 <- as.matrix(cbind(1,Turismo_limpio_BIC$IPC, Turismo_limpio_BIC$PH,Turismo_limpio_BIC$GH, Turismo_limpio_BIC$PC))
y3 <- as.matrix(TUR_log_limpio, ncol=1)
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
F0_3 # 194.1152

# p value para H0: B0 = ... = Bi = 0
pf(q=F0_3, df1=k5, df2=n5-k5-1, lower.tail=FALSE) 
# 5.870681e-34 < 0.05 => Se rechaza H0


## Intervalos de respuesta media y de predicción
# Intervalos de respuesta media con comando
x0.3 <- data.frame(
  IPC = 75,
  PH = 100,
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
#MPE Akaike = -0.0005170334
#MPE BIC = -0.0004994916
#BIC está más cercano a 0

#MEJOR MODELO BIC