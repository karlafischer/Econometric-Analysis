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

#LENTICOS
lenticos_semaforo <- read.csv("LENTICO-SEMAFORO.csv")

# División de la base de datos en 80% / 20%
set.seed(2508)

split <- initial_split(lenticos_semaforo, prop = 0.8, strata = "SEMAFORO")

entrenamiento <- training(split)
prueba  <- testing(split)
prop.table(table(lenticos_semaforo$SEMAFORO))
prop.table(table(entrenamiento$SEMAFORO))
prop.table(table(prueba$SEMAFORO))
# La estratificación funcionó correctamente (las proporciones de 0 y 1 para
# training y testing son casi idénticas a la proporción de la base original)

# Declaración de variable respuesta
SEMAFORO_glm <- entrenamiento$SEMAFORO

# Declaración de variables regresoras
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
TRANSPARENCIA <- entrenamiento$TRANSPARENCIA
ABS_UV <- entrenamiento$ABS_UV
SST <- entrenamiento$SST
TURBIEDAD <- entrenamiento$TURBIEDAD
AS_TOT <- entrenamiento$AS_TOT 
HG_TOT <- entrenamiento$HG_TOT 
PB_TOT <- entrenamiento$PB_TOT
DUR_TOT <- entrenamiento$DUR_TOT
CONDUC_CAMPO <- entrenamiento$CONDUC_CAMPO 
OD_PORC <- entrenamiento$OD_.

# Modelo léntico completo
modelo_lentico_completo <- glm(SEMAFORO_glm ~ E_COLI + COLI_FEC + pH_CAMPO + TEMP_AGUA + COT + DBO_TOT + DQO_TOT + N_TOT + P_TOT + COLOR_VER + TRANSPARENCIA + ABS_UV + SST + TURBIEDAD + AS_TOT + HG_TOT + PB_TOT + DUR_TOT + CONDUC_CAMPO + OD_PORC, data = entrenamiento, family = "binomial")
summary(modelo_lentico_completo)
deviance(modelo_lentico_completo)
# Se decide quitar P_TOT porque todos sus valores son 1

# Modelo léntico depurado 
modelo_lentico <- glm(SEMAFORO_glm ~ COLI_FEC + pH_CAMPO + TEMP_AGUA + COT + DBO_TOT + N_TOT + COLOR_VER + TRANSPARENCIA + ABS_UV + TURBIEDAD + AS_TOT + HG_TOT + PB_TOT + DUR_TOT + CONDUC_CAMPO + OD_PORC, data = entrenamiento, family = "binomial")
summary(modelo_lentico)
anova(modelo_lentico)
deviance(modelo_lentico)

## --- AIC --- ##
modelo_aic_len <- step(modelo_lentico, direction="both")
modelo_aic_len

modelo_lentico_AIC <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + ABS_UV + PB_TOT + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_lentico_AIC)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lentico_AIC,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lentico_AIC)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + ABS_UV + DUR_TOT + OD_PORC, data = entrenamiento)

## Prueba de hipotesis sobre subconjuntos de parámetros (ABS_UV no es significativa)
alpha <- 0.05
r <- 1
modelo_AIC_sinuv <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + PB_TOT + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_AIC_sinuv)
lambda_uv <- deviance(modelo_AIC_sinuv) - deviance(modelo_lentico_AIC)
lambda_uv < qchisq(alpha, r, lower.tail = FALSE) # TRUE
# Se concluye que se puede quitar ABS_UV del modelo

## Prueba de hipotesis sobre subconjuntos de parámetros (PB_TOT no es significativa)
modelo_AIC_sinpb <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + ABS_UV + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_AIC_sinpb)
lambda_pb <- deviance(modelo_AIC_sinpb) - deviance(modelo_lentico_AIC)
lambda_pb < qchisq(alpha, r, lower.tail = FALSE) # FALSE
# Se concluye que no se puede quitar PB_TOT del modelo

## Prueba de hipotesis sobre subconjuntos de parámetros (DBO_TOT no es significativa)
modelo_AIC_sindbo <- glm(SEMAFORO_glm ~ TEMP_AGUA + PB_TOT + ABS_UV + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_AIC_sindbo)
lambda_dbo <- deviance(modelo_AIC_sindbo) - deviance(modelo_lentico_AIC)
lambda_dbo < qchisq(alpha, r, lower.tail = FALSE) # FALSE
# Se concluye que no se puede quitar DBO_TOT del modelo

## Prueba de hipotesis sobre subconjuntos de parámetros (DBO_TOT y PB_TOT no son significativas)
r2 <- 2
modelo_AIC_sindos <- glm(SEMAFORO_glm ~ TEMP_AGUA + ABS_UV + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_AIC_sindos)
lambda_dos <- deviance(modelo_AIC_sindos) - deviance(modelo_lentico_AIC)
lambda_dos < qchisq(alpha, r2, lower.tail = FALSE) # FALSE
# Se concluye que no se puede quitar DBO_TOT y PB_TOT del modelo simultáneamente

## --- Modelo AIC sin ABS_UV --- ##
modelo_lent_AIC <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + PB_TOT + DUR_TOT + OD_PORC, data = entrenamiento, family="binomial")
summary(modelo_lent_AIC)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lent_AIC,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lent_AIC)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + DUR_TOT + OD_PORC, data = entrenamiento)

## --- Modelo AIC sin ABS_UV con transformaciones --- ##
# Aplicar transformación en DUR_TOT con lambda = 0.133289
entrenamiento$DUR_TOT_AIC_trans <- entrenamiento$DUR_TOT^0.133289

modelo_lentico_AIC_trans <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + PB_TOT + DUR_TOT_AIC_trans + OD_., data = entrenamiento, family = "binomial")
summary(modelo_lentico_AIC_trans)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lentico_AIC_trans,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lentico_AIC_trans)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + DUR_TOT_AIC_trans + OD_PORC, 
           data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA, data = entrenamiento) 
boxTidwell(SEMAFORO_glm ~ DUR_TOT_AIC_trans, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ OD_PORC, data = entrenamiento)

# Probar la desviación del modelo con transformaciones
deviance(modelo_lentico_AIC_trans)

n <- length(SEMAFORO_glm)
n
k <- 5

# Desviación del modelo 
deviance(modelo_lentico_AIC_trans) < qchisq(alpha,n-k,lower.tail = FALSE) # TRUE
# Se concluye que el modelo ajustado es adecuado

## Predicción sobre la PRUEBA
# Aplicar las mismas transformaciones a prueba
prueba$DUR_TOT_AIC_trans <- prueba$DUR_TOT^0.133289

proba_predic_AIC <- predict(modelo_lentico_AIC_trans, newdata = prueba, type = "response")

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
n <- nrow(entrenamiento)
n
modelo_bic_len <- step(modelo_lentico, direction = "both", k=log(n))
modelo_bic_len
modelo_lentico_BIC <- glm(formula = SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + DUR_TOT, 
                          data = entrenamiento, family = "binomial")
summary(modelo_lentico_BIC)


# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lentico_BIC,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lentico_BIC)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + DUR_TOT, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ DUR_TOT, data = entrenamiento)

## Prueba de hipotesis sobre subconjuntos de parámetros (DBO_TOT no es significativa)
alpha <- 0.05
r <- 1
modelo_BIC_sindbo <- glm(SEMAFORO_glm ~ TEMP_AGUA + DUR_TOT, data = entrenamiento, family="binomial")
summary(modelo_BIC_sindbo)
lambda_dbo2 <- deviance(modelo_BIC_sindbo) - deviance(modelo_lentico_BIC)
lambda_dbo2 < qchisq(alpha, r, lower.tail = FALSE) # FALSE
# Se concluye que no se puede quitar DBO_TOT del modelo

## --- Modelo BIC con transformación --- ###
# Aplicar transformación en DUR_TOT con lambda = 0.11728
entrenamiento$DUR_TOT_BIC_trans <- entrenamiento$DUR_TOT^0.11728

modelo_lentico_BIC_trans <- glm(SEMAFORO_glm ~ TEMP_AGUA + DBO_TOT + DUR_TOT_BIC_trans,
                                data = entrenamiento, 
                                family = "binomial")
summary(modelo_lentico_BIC_trans)

# Coeficientes de determinación (Pseudo R2)
PseudoR2(modelo_lentico_BIC_trans,which=c("CoxSnell","Nagelkerke","McFadden"))

## Supuestos
# Multicolinealidad
vif(modelo_lentico_BIC_trans)

# Linealidad
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA + DUR_TOT_BIC_trans, 
           data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ TEMP_AGUA, data = entrenamiento)
boxTidwell(SEMAFORO_glm ~ DUR_TOT_BIC_trans, data = entrenamiento)

# Probar la desviación del modelo con transformaciones
deviance(modelo_lentico_BIC_trans)

n <- length(SEMAFORO_glm)
n
k2 <- 3

# Desviación del modelo 
deviance(modelo_lentico_BIC_trans) < qchisq(alpha,n-k2,lower.tail = FALSE) # TRUE
# Se concluye que el modelo ajustado es adecuado

## Predicción sobre la PRUEBA
# Aplicar las mismas transformaciones a prueba
prueba$DUR_TOT_BIC_trans <- prueba$DUR_TOT^0.11728

proba_predic_BIC <- predict(modelo_lentico_BIC_trans, newdata = prueba, type = "response")

# Convertir probabilidades a variable dicotómica 
pred_BIC <- ifelse(proba_predic_BIC > 0.75, 1, 0)

## Matriz de confusión para Modelo AIC con transformaciones
MC_BIC <- confusionMatrix(as.factor(pred_BIC), as.factor(prueba$SEMAFORO), positive="1")
MC_BIC

## Interpretación de la Matriz de confusión
MC_BIC_prop <- prop.table(MC_BIC$table)
MC_BIC_prop
cat("--- EVALUACIÓN DE LA MATRIZ DE CONFUSIÓN AIC (%) ---\n")
cat("Accuracy:", round(MC_BIC$overall["Accuracy"] * 100, 2), "%\n")
cat("Falsos Positivos (FP):", round(MC_BIC_prop[2,1]*100, 2), "%\n")
cat("Falsos Negativos (FN):", round(MC_BIC_prop[1,2]*100, 2), "%\n")
cat("Predicciones Correctas (TN + TP):", round((MC_BIC_prop[1,1] + MC_BIC_prop[2,2])*100, 2), "%\n")

AIC(modelo_lentico_AIC_trans, modelo_lentico_BIC_trans)
BIC(modelo_lentico_AIC_trans, modelo_lentico_BIC_trans)

## Se decide que el modelo que mejor se ajusta a los datos muestrales es el 
## modelo AIC transformado.

## ---- MEJOR MODELO --- ###
modelo_lentico_AIC_trans
summary(modelo_lentico_AIC_trans)

# Prueba de Hipótesis Individual para modelo logístico (Wald test)
alpha <- 0.05
cat("=== SIGNIFICANCIA INDIVIDUAL (Wald test) con α =", alpha, "===\n")
coeficientes <- coef(modelo_lentico_AIC_trans)
summary_mod <- summary(modelo_lentico_AIC_trans)$coefficients

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
beta <- as.matrix(coef(modelo_lentico_AIC_trans))
beta
# Matriz de diseño
x_i <- cbind(1, 
             entrenamiento$TEMP_AGUA, 
             entrenamiento$DBO_TOT, 
             entrenamiento$PB_TOT, 
             entrenamiento$DUR_TOT_AIC_trans, 
             entrenamiento$OD_.)
x_i

# Log-odds
nu <- x_i %*% beta

# Probabilidad estimada
prob_estimada <- exp(nu)/(1+exp(nu))
prob_estimada

# Cambios en el cociente de posibilidades (odds ratio) para incrementos y
# de disminución 1 unidad en TEMP_AGUA
beta_temp <- coef(modelo_lentico_AIC_trans)["TEMP_AGUA"]
OR_temp_mas <- exp(beta_temp)
OR_temp_menos <- exp(-beta_temp)
cat("Cociente de Posibilidades para TEMP_AGUA (+1°C):", round(OR_temp_mas, 4))
cat("Cociente de Posibilidades para TEMP_AGUA (-1°C):", round(OR_temp_menos, 4))

# Cambios en el cociente de posibilidades (odds ratio) para incrementos de 
# y disminuciones de 1 unidad en DUR_TOT
beta_dur <- coef(modelo_lentico_AIC_trans)["DUR_TOT_AIC_trans"]
# Percentiles de DUR_TOT en entrenamiento
dur_vals <- quantile(entrenamiento$DUR_TOT, probs = c(0.25, 0.5, 0.75))
cat("\n=== DUR_TOT (transformación ^0.133289) ===\n")
for(d in dur_vals) {
  # Aumento de +1 unidad
  delta_mas <- (d+1)^0.133289 - d^0.133289
  OR_mas <- exp(beta_dur * delta_mas)
  # Disminución de -1 unidad (solo si d-1 > 0)
  if(d - 1 > 0) {
    delta_menos <- (d-1)^0.133289 - d^0.133289
    OR_menos <- exp(beta_dur * delta_menos)
    cat(sprintf("  DUR_TOT = %.1f: +1 → OR = %.4f, -1 → OR = %.4f\n", d, OR_mas, OR_menos))
  } else {
    cat(sprintf("  DUR_TOT = %.1f: +1 → OR = %.4f, -1 → no factible (valor negativo)\n", d, OR_mas))
  }
}

# Cambios en el cociente de posibilidades (odds ratio) para incrementos de 
# y disminuciones de 0.1 unidades en OD_PORC
beta_od <- coef(modelo_lentico_AIC_trans)["OD_."]
OR_od_mas <- exp(beta_od)       
OR_od_menos <- exp(-beta_od)     
cat("\n=== OD_PORC (OD_.) ===\n")
cat(sprintf("  +1 unidad → OR = %.4f\n", OR_od_mas))
cat(sprintf("  -1 unidad → OR = %.4f\n", OR_od_menos))
