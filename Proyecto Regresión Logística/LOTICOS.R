#Paquetes 
library(car)
library(corrplot)
library(lmtest)
library(corrplot)
library(readr)
library(leaps)
library(segmented)
library(glmnet)
library(rsample)
library(DescTools)
library(caret)
library(splines)

#####---------------------------------------
#LOTICOS
loticos_semaforo <- read.csv("LOTICO-SEMAFORO.csv")

# División de la base de datos en 80% / 20%
set.seed(2508)

split <- initial_split(loticos_semaforo, prop = 0.8, strata = "SEMAFORO")

entrenamiento <- training(split)
prueba  <- testing(split)
prop.table(table(loticos_semaforo$SEMAFORO))
prop.table(table(entrenamiento$SEMAFORO))
prop.table(table(prueba$SEMAFORO))
# La estratificación funcionó correctamente (las proporciones de 0 y 1 para
# training y testing son idénticas a la proporción de la base original)

# Declaración de variable respuesta
SEMAFORO_glm <- entrenamiento$SEMAFORO

# Declaración de variables regresoras
TOX_V_15 <- entrenamiento$TOX_V_15_UT..NORMAL.
E_COLI <- entrenamiento$E_COLI
COLI_FEC <- entrenamiento$COLI_FEC
pH_CAMPO <- entrenamiento$pH_CAMPO
TEMP_AGUA <- entrenamiento$TEMP_AGUA
COT <- entrenamiento$COT
DBO_TOT <- entrenamiento$DBO_TOT
DQO_TOT <- entrenamiento$DQO_TOT
N_TOT  <- entrenamiento$N_TOT
P_TOT <- entrenamiento$P_TOT
COLOR_VER <- entrenamiento$COLOR_VER
ABS_UV <- entrenamiento$ABS_UV
SST <- entrenamiento$SST
TURBIEDAD <- entrenamiento$TURBIEDAD
AS_TOT <- entrenamiento$AS_TOT 
HG_TOT <- entrenamiento$HG_TOT 
PB_TOT <- entrenamiento$PB_TOT
DUR_TOT <- entrenamiento$DUR_TOT
CONDUC_CAMPO <- entrenamiento$CONDUC_CAMPO 
OD_PORC <- entrenamiento$OD_.

# Modelo lótico completo
modelo_lotico_completo <- glm(SEMAFORO_glm ~ TOX_V_15 + E_COLI + COLI_FEC + pH_CAMPO + TEMP_AGUA + COT + DBO_TOT + DQO_TOT + N_TOT + P_TOT + COLOR_VER + ABS_UV + SST + TURBIEDAD + AS_TOT + HG_TOT + PB_TOT + DUR_TOT + CONDUC_CAMPO + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_lotico_completo)
deviance(modelo_lotico_completo)

# Modelo lótico depurado
modelo_lotico <- glm(SEMAFORO_glm ~ TOX_V_15 + COLI_FEC + pH_CAMPO + TEMP_AGUA + COT + N_TOT + P_TOT + COLOR_VER + ABS_UV + TURBIEDAD + AS_TOT + HG_TOT + PB_TOT + DUR_TOT + CONDUC_CAMPO + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_lotico)
anova(modelo_lotico)
deviance(modelo_lotico)

## --- AIC --- ## 
modelo_aic_lo <- step(modelo_lotico, direction="both")
modelo_aic_lo

modelo_lotico_AIC <- glm(SEMAFORO_glm ~ COLI_FEC + TEMP_AGUA + COT + N_TOT + COLOR_VER 
                         + ABS_UV + HG_TOT + PB_TOT, 
                         data = entrenamiento, family="binomial")
summary(modelo_lotico_AIC)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lotico_AIC,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lotico_AIC)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + COLOR_VER + ABS_UV , data = entrenamiento)

## --- Modelo AIC con transformaciones --- ##
# Aplicar transformación en COLOR_VER con lambda = 0.25695 y en ABS_UV con 
# lambda = 0.74388
entrenamiento$COLOR_VER_AIC_trans <- entrenamiento$COLOR_VER^0.25695
entrenamiento$ABS_UV_AIC_trans <- entrenamiento$ABS_UV^0.74388

modelo_lotico_AIC_trans <- glm(SEMAFORO_glm ~ COLI_FEC + TEMP_AGUA + COT + N_TOT + 
                                 COLOR_VER_AIC_trans + ABS_UV_AIC_trans + HG_TOT + PB_TOT,
                               data = entrenamiento, 
                               family = "binomial") # F-Parcial para COLI_FEC

summary(modelo_lotico_AIC_trans)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lotico_AIC_trans,which=c("CoxSnell","Nagelkerke","McFadden"))

## Prueba de hipotesis sobre subconjuntos de parámetros (COLI_FEC no es significativa)
alpha <- 0.05
r <- 1
modelo_AICtrans_sincolifec <- glm(SEMAFORO_glm ~ TEMP_AGUA + COT + N_TOT + 
                                    COLOR_VER_AIC_trans + ABS_UV_AIC_trans + HG_TOT + PB_TOT,
                                  data = entrenamiento, 
                                  family = "binomial")
summary(modelo_AICtrans_sincolifec)
lambda_colifec <- deviance(modelo_AICtrans_sincolifec) - deviance(modelo_lotico_AIC_trans)
lambda_colifec < qchisq(alpha, r, lower.tail = FALSE) # TRUE
# Se concluye que se puede quitar COLI_FEC del modelo con transformaciones

## --- Modelo AIC con transformaciones sin COLI_FEC --- ##
modelo_lot_AIC_trans <- glm(SEMAFORO_glm ~  TEMP_AGUA + COT + N_TOT + 
                              COLOR_VER_AIC_trans + ABS_UV_AIC_trans + HG_TOT + PB_TOT,
                            data = entrenamiento, 
                            family = "binomial")
summary(modelo_lot_AIC_trans)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lot_AIC_trans,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lot_AIC_trans)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + COLOR_VER_AIC_trans + ABS_UV_AIC_trans, 
           data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA, data = entrenamiento) # Hay problemas con esta variable
boxTidwell(SEMAFORO_glm ~ COLOR_VER_AIC_trans, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ ABS_UV_AIC_trans, data = entrenamiento)

# Probar la desviación del modelo con transformaciones
deviance(modelo_lot_AIC_trans)

n <- length(SEMAFORO_glm)
n
k <- 7

# Desviación del modelo 
deviance(modelo_lot_AIC_trans) < qchisq(alpha,n-k,lower.tail = FALSE) # TRUE
# Se concluye que el modelo ajustado es adecuado

## Predicción sobre la PRUEBA
# Aplicar las mismas transformaciones a prueba
prueba$COLOR_VER_AIC_trans <- prueba$COLOR_VER^0.25695
prueba$ABS_UV_AIC_trans <- prueba$ABS_UV^0.74388

proba_predic_AIC <- predict(modelo_lot_AIC_trans, newdata = prueba, type = "response")

# Convertir probabilidades a variable dicotómica 
pred_AIC <- ifelse(proba_predic_AIC > 0.75, 1, 0)

## Matriz de confusión para Modelo AIC con transformaciones
MC_AIC <- confusionMatrix(as.factor(pred_AIC), as.factor(prueba$SEMAFORO), positive="1")
MC_AIC

## Interpretación de la Matriz de confusión
MC_AIC_prop <- prop.table(MC_AIC$table)
MC_AIC_prop
cat("--- EVALUACIÓN DE LA MATRIZ DE CONFUSIÓN AIC (%) ---\n")
cat("Accuracy:", round(MC_AIC$overall["Accuracy"] * 100, 2), "%\n")
cat("Falsos Positivos (FP):", round(MC_AIC_prop[2,1]*100, 2), "%\n")
cat("Falsos Negativos (FN):", round(MC_AIC_prop[1,2]*100, 2), "%\n")
cat("Predicciones Correctas (TN + TP):", round((MC_AIC_prop[1,1] + MC_AIC_prop[2,2])*100, 2), "%\n")


## --- BIC --- ##
## BIC con variable pH_CAMPO
n <- nrow(entrenamiento)
modelo_bic_lo <- step(modelo_lotico, direction = "both", k=log(n))
modelo_bic_lo

modelo_lotico_BIC <- glm(SEMAFORO_glm ~ pH_CAMPO + N_TOT + COLOR_VER + ABS_UV + AS_TOT + HG_TOT, data = entrenamiento, family = "binomial")
summary(modelo_lotico_BIC)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lotico_BIC,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lotico_BIC)

# Linealidad
boxTidwell(SEMAFORO_glm ~ pH_CAMPO + COLOR_VER + ABS_UV , data = entrenamiento) # NO FUNCIONA
boxTidwell(SEMAFORO_glm ~ pH_CAMPO, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ COLOR_VER, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ ABS_UV, data = entrenamiento)


## --- BIC sin variable pH_CAMPO --- ##
modelo_lot_BIC_sinph <- glm(SEMAFORO_glm ~ N_TOT + COLOR_VER + ABS_UV + AS_TOT + HG_TOT, data = entrenamiento, family = "binomial")
summary(modelo_lot_BIC_sinph)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lot_BIC_sinph,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lot_BIC_sinph)

# Linealidad
boxTidwell(SEMAFORO_glm ~  COLOR_VER + ABS_UV , data = entrenamiento)
boxTidwell(SEMAFORO_glm ~  COLOR_VER, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~  ABS_UV , data = entrenamiento)

## --- BIC sin pH con transformaciones --- ##
# Aplicar transformación en COLOR_VER con lambda = 0.23686 y en ABS_UV con 
# lambda = 0.70773
entrenamiento$COLOR_VER_BIC_trans <- entrenamiento$COLOR_VER^0.23686
entrenamiento$ABS_UV_BIC_trans <- entrenamiento$ABS_UV^0.70773

modelo_lotico_BIC_trans <- glm(SEMAFORO_glm ~ N_TOT + COLOR_VER_BIC_trans + ABS_UV_BIC_trans + AS_TOT + HG_TOT,
                               data = entrenamiento, family = "binomial")

summary(modelo_lotico_BIC_trans)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lotico_BIC_trans,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lotico_BIC_trans)

# Linealidad
boxTidwell(SEMAFORO_glm ~ COLOR_VER_BIC_trans + ABS_UV_BIC_trans, 
           data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ COLOR_VER_BIC_trans, 
           data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ ABS_UV_BIC_trans, 
           data = entrenamiento)

# Probar la desviación del modelo con transformaciones
deviance(modelo_lotico_BIC_trans)

n <- length(SEMAFORO_glm)
n
k2 <- 5

# Desviación del modelo 
deviance(modelo_lotico_BIC_trans) < qchisq(alpha,n-k2,lower.tail = FALSE) # TRUE
# Se concluye que el modelo ajustado es adecuado


## Predicción sobre la PRUEBA
# Aplicar las mismas transformaciones a prueba
prueba$COLOR_VER_BIC_trans <- prueba$COLOR_VER^0.23686
prueba$ABS_UV_BIC_trans <- prueba$ABS_UV^0.70773

proba_predic_BIC <- predict(modelo_lotico_BIC_trans, newdata = prueba, type = "response")

# Convertir probabilidades a variable dicotómica 
pred_BIC <- ifelse(proba_predic_BIC > 0.75, 1, 0)

## Matriz de confusión para Modelo BIC con transformaciones
MC_BIC <- confusionMatrix(as.factor(pred_BIC), as.factor(prueba$SEMAFORO), positive="1")
MC_BIC

## Interpretación de la Matriz de confusión
MC_BIC_prop <- prop.table(MC_BIC$table)
MC_BIC_prop
cat("--- EVALUACIÓN DE LA MATRIZ DE CONFUSIÓN BIC (%) ---\n")
cat("Accuracy:", round(MC_BIC$overall["Accuracy"] * 100, 2), "%\n")
cat("Falsos Positivos (FP):", round(MC_BIC_prop[2,1]*100, 2), "%\n")
cat("Falsos Negativos (FN):", round(MC_BIC_prop[1,2]*100, 2), "%\n")
cat("Predicciones Correctas (TN + TP):", round((MC_BIC_prop[1,1] + MC_BIC_prop[2,2])*100, 2), "%\n")


AIC(modelo_lot_AIC_trans, modelo_lotico_BIC_trans)
BIC(modelo_lot_AIC_trans, modelo_lotico_BIC_trans)

## Se decide que el modelo que mejor se ajusta a los datos muestrales es el 
## modelo AIC transformado.

## --- MEJOR MODELO --- ##
modelo_lot_AIC_trans

# Prueba de Hipótesis Individual para modelo logístico (Wald test)
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL (Wald test) con α =", alpha, "===\n")
coeficientes <- coef(modelo_lot_AIC_trans)
summary_mod <- summary(modelo_lot_AIC_trans)$coefficients

for(i in 1:length(coeficientes)) {
  beta <- coeficientes[i]
  error_std <- summary_mod[i, 2]
  z_calc <- beta / error_std
  p_valor <- 2 * pnorm(abs(z_calc), lower.tail = FALSE)
  
  cat(names(coeficientes)[i], "→ z =", round(z_calc, 4), 
      ", p =", format(p_valor, scientific = TRUE, digits = 4), 
      "→", ifelse(p_valor < alpha, "Significativo", "No significativo"), "\n")
}

## --- Cálculo de Probabilidades --- ##
beta <- as.matrix(coef(modelo_lot_AIC_trans))
beta
# Matriz de diseño
x_i <- cbind(1, 
             entrenamiento$TEMP_AGUA, 
             entrenamiento$COT, 
             entrenamiento$N_TOT, 
             entrenamiento$COLOR_VER_AIC_trans, 
             entrenamiento$ABS_UV_AIC_trans, 
             entrenamiento$HG_TOT, 
             entrenamiento$PB_TOT)
x_i

# Log-odds
nu <- x_i %*% beta

# Probabilidad estimada
prob_estimada <- exp(nu)/(1+exp(nu))
prob_estimada

# Cambios en el cociente de posibilidades (odds ratio) para incrementos y
# de disminución 1 unidad en TEMP_AGUA
beta_temp <- coef(modelo_lot_AIC_trans)["TEMP_AGUA"]
OR_temp_mas <- exp(beta_temp)
OR_temp_menos <- exp(-beta_temp)
cat("Cociente de Posibilidades para TEMP_AGUA (+1°C):", round(OR_temp_mas, 4))
cat("Cociente de Posibilidades para TEMP_AGUA (-1°C):", round(OR_temp_menos, 4))

# Cambios en el cociente de posibilidades (odds ratio) para incrementos de 
# y disminuciones de 1 unidad en COLOR_VER
beta_color <- coef(modelo_lot_AIC_trans)["COLOR_VER_AIC_trans"]
# Percentiles de COLOR_VER en entrenamiento
c_vals <- quantile(entrenamiento$COLOR_VER, probs = c(0.25, 0.5, 0.75))
cat("\n=== COLOR_VER (transformación ^0.25695) ===\n")
for(c in c_vals) {
  delta_mas <- (c+1)^0.25695 - c^0.25695 # Aumento +1
  OR_mas <- exp(beta_color * delta_mas)
  if(c - 1 > 0) {
    delta_menos <- (c-1)^0.25695 - c^0.25695 # Disminución -1 (solo si c-1 > 0)
    OR_menos <- exp(beta_color * delta_menos)
  } else {
    OR_menos <- NA
  }
  cat(sprintf("  COLOR_VER = %.1f: +1 → OR = %.4f", c, OR_mas))
  if(!is.na(OR_menos)) {
    cat(sprintf(",   -1 → OR = %.4f\n", OR_menos))
  } else {
    cat("   -1 → no factible (valor inicial muy pequeño)\n")
  }
}

# Cambios en el cociente de posibilidades (odds ratio) para incrementos de 
# y disminuciones de 0.1 unidades en ABS_UV
beta_abs <- coef(modelo_lot_AIC_trans)["ABS_UV_AIC_trans"]
# Percentiles de ABS_UV en entrenamiento
a_vals <- quantile(entrenamiento$ABS_UV, probs = c(0.25, 0.5, 0.75))
cat("\n=== ABS_UV (transformación ^0.74388) ===\n")
for(a in a_vals) {
  delta_mas <- (a+0.1)^0.74388 - a^0.74388 # Aumento +0.1
  OR_mas <- exp(beta_abs * delta_mas)
  if(a - 0.1 > 0) {
    delta_menos <- (a-0.1)^0.74388 - a^0.74388 # Disminución de -0.1 (si a - 0.1 > 0)
    OR_menos <- exp(beta_abs * delta_menos)
  } else {
    OR_menos <- NA
  }
  cat(sprintf("  ABS_UV = %.3f: +0.1 → OR = %.4f", a, OR_mas))
  if(!is.na(OR_menos)) {
    cat(sprintf(",   -0.1 → OR = %.4f\n", OR_menos))
  } else {
    cat("   -0.1 → no factible (valor muy pequeño)\n")
  }
}
