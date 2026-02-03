---
title: Visualización resultados del Mk2
layout: home
nav_order: 11
index: true
redirect: false
parent: Tutorials
math: katex
---

Creado por Rosana Zenil-Ferguson (Enero 2025)
Actualizado Nicolás Castillo-Rodríguez (Febrero 2026)

## Organizando tus librerias

Para este tutorial necesitaremos los siguientes paquetes. **Importante**: El paquete RevGadgets esta disponible directamente con ``install.packages("RevGadgets")`` pero vamos a revisar una versión especial que es la que grafica mapas estocásticos. Por esa razón vamos a bajar la versión que se indica aquí, así que. sigue las instrucciones como se indican abajo. 

```
# install.packages("devtools")
# install.packages("coda")
# install.packages("phytools")
# install.packages ("ggplot2")

## Instalar treeio
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("treeio")

##Instalar ggtree
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# Importante: no instales install.packages("RevGadgets") porque esto instala desde CRAN. Queremos una version especifica. Por favor sigue la instruccion de abajo. Si ya habias instalado revgadgets por favor reinicia y reinstala asi

# devtools::install_github("revbayes/RevGadgets@stochastic_map",force=TRUE)
library(coda)
library(RevGadgets)
library(phytools)
library(ggplot2)

setwd("./mk2/output/")
```

## Convergencia del MCMC

Revisemos la convergencia básica

Baja estos archivos
+ [Primera corrida del MCMC](downloads/mk2_polinizador_run_1.log)
+ [Segunda corrida del MCMC](downloads/mk2_polinizador_run_2.log)

```
# Agrega to directorio de trabajo
# Por ejemplo el mio esta en 
# setwd("./mk2/output/")

# Con el paquete coda
mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
mcmc_trace <- as.mcmc(mcmc_run1[[1]])
traceplot(mcmc_trace)
effectiveSize(mcmc_trace)
```

Ahora hagamos un mejor trabajo con los gráficos de convergencia

```
# Esta parte se hace con ggplot2 y con RevGadgets
mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
mcmc_run1<-data.frame(mcmc_run1)
mcmc_run1<- cbind(mcmc_run1,run=rep("run 1",length(mcmc_run1$Iteration)))

mcmc_run2 <- readTrace(path ="mk2_polinizador_run_2.log", burnin = 0.1)
mcmc_run2<-data.frame(mcmc_run2)
mcmc_run2<- cbind(mcmc_run2,run=rep("run 2",length(mcmc_run2$Iteration)))

mcmc_table<-rbind(mcmc_run1,mcmc_run2)

trace_plot<- ggplot(mcmc_table, aes(x=Iteration,y=Posterior,group=run))+
  geom_line(aes(color=run))+
  theme_classic()
trace_plot
```

## Estimadores de las tasas

Queremos graficar y obtener estadísticas resumen de los parametros de interés. Recordemos que para el modelo mk2 los parámetros que estimamos son

+ $$q_{01}$$
+ $$q_{10}$$
+ Frecuencias de la raíz

```
# Esta parte se hace con RevGadgets

mcmc_run1 <- readTrace(path ="mk2_polinizador_run_1.log", burnin = 0.1)
# Estadisticas resumen
## MAP= maximo a posteriori
## Mean= media de la distribucion
## Cuantiles 2.5 y 97.5 define el intervalo de credibilidad
summarizeTrace(trace = mcmc_run1, vars =c("q_01","q_10","root_frequencies[1]","root_frequencies[2]"))

## Graficando las distribuciones posteriores de las transiciones
plotTrace(trace = mcmc_run1, vars = c("q_01","q_10"))[[1]]

plotTrace(trace = mcmc_run1, vars = c("root_frequencies[1]","root_frequencies[2]"))[[1]]
```

En lo personal prefiero graficar violines, y tener más control en los colores especialmente cuando hay muchos parámetros que comparar (veremos en el ejercicio de HiSSE que esto es super útil).

```
# Esta parte se hace con ggplot2

traitcols<-c("#7678ED","#F35B04")

mk2 <- read.table("mk2_polinizador_run_1.log", header=TRUE)
mk2 <- mk2[-seq(1,5000,1),] # nunca olvidemos quitar el burn-in con esta instruccion tengo mas contro de cuantas iteraciones remuevo.


transition_rates <- data.frame(dens=c(mk2$q_01,mk2$q_10),
                               rate=rep(c("q_01","q_10"),
                               each=length(mk2$q_01)))

violin_transitions<- ggplot(transition_rates,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE, quantiles = c(0.025, 0.975), quantile.linetype = 1) +
  labs(title="Transition rates")+
  scale_fill_manual( values = traitcols)+
  xlab("Transition type")+
  ylab("Rate")+
  scale_x_discrete(labels = c("q_01" = expression(italic(q)["0,1"]),
                              "q_10" = expression(italic(q)["1,0"]))) + 
  theme_classic() +
  theme(legend.position = "none")

violin_transitions
```

## Prueba de hipótesis: Las tasas son iguales?

Una de las ventajas más importante de las distribuciones posteriores es que rápidamente podemos concluir si dos parámetros son iguales y distintos. Pero en general, si queremos publicar estos resultados debemos dar una estadística formal (como un p-valor) para que acepten nuestro trabajo para publicación. ¿Cómo formalizamos esta idea con Bayesiana?

Construimos una estadística de resumen llamémosla $$D= q_{01}-q_{10}$$ y grafiquemos.

```
# Esta parte se hace con ggplot2

D <- data.frame(dens=(mk2$q_01-mk2$q_10),
                rate=rep(c("Difference"),
                         length(mk2$q_01)))

violin_difference<- ggplot(D,aes(x=rate,y=dens, fill=rate))+
  geom_violin(trim=FALSE, quantiles = c(0.025, 0.975), quantile.linetype = 1) +
  labs(title="Hypothesis testing")+
  scale_fill_manual( values = "hotpink")+
  xlab("Test Statistic")+
  ylab("Difference")+
  geom_hline(yintercept = 0,linetype="dashed",lwd=1)+
  theme_classic() +
  theme(legend.position = "none")

# Observemos que 0 cruza la diferencia entre las tasas lo que quiere decir que pueden ser iguales
violin_difference
```

Formalizando la prueba de hipótesis
$$ H_0: D=0$$

```
cuantiles <- quantile(D$dens, probs=c(0.025,0.975))
cuantiles
```

Observamos que el intervalo de credibilidad del 95% de la distribución de $$D$$ son (-0.19, 0.01). Cómo 0 pertenece a este intervalo entonces $$P(H_0\lvert Datos)>0.05$$ lo que significa que las transiciones son iguales con probabilidad mayor al 5%. **Importante**: Noten que en las pruebas de hipótesis bayesiana no utilizo, p-valores, ni significancia, ni rechazo, solo probabilidad e intervalos de credibilidad. Esto sucede porque la manera en que interpretamos probabilidad en bayesiana es distinto. Tengan cuidado cuando interpreten. 

## Reconstrucciones ancestrales

Hay dos tipos de reconstrucciones ancestrales

1. Reconstruccion ancestral de estados marginal: En estadística bayesiana condicionamos las probabilidades a lo que pudo pasar en cada nodo interno y promediamos lo que sucedió en el resto, hacemos esto con cada nodo y obtenemos las distribuciones de probabilidad marginal para cada nodo. Lo que graficamos es el MAP (máximo *a posteriori*) de cada uno de los nodos con la probabilidad que se le asigna. 

    + [Árbol con reconstruccion ancestral](downloads/asr_mk2.tree)

    ```
       # Esta parte se hace con RevGadgets

       anc_states <- processAncStates(path ="asr_mk2.tree",state_labels=c("0"="insect","1"="wind"))
       plotAncStatesMAP(t = anc_states, tree_layout="rectangular",
                 state_transparency = 0.5,
                 node_size = c(0.1, 5),
                 tip_labels_size = 2,
                 tip_states_size=2,
                 node_color=c("red","blue"))
       # produce the plot object, showing MAP states at nodes.
       # color corresponds to state, size to the state's posterior probability
    ```

2. Mapas estocásticos: Esta reconstrucción se enfoca en describir que es lo que pudo suceder en las ramas del árbol. Bajo una serie de simulaciones que va de la raíz a la punta bajo nuestro modelo (en este caso Mk2), se simula lo que sucedio en las ramas y se grafica el MAP de las simulaciones en pedacitos de las ramas que son divididos de manera igual. 

    + [Árbol filogenético en Nexus](downloads/poliniza_arbol.nex)
    + [Mapas estocásticos](downloads/stochmap_mk2_polinizador_run_1.log)

```
      #Todo esto funciona con la rama en desarrollo de RevGadgets
      #Lee el arbol
      file <- "../data/poliniza_arbol.nex" #Check path, remember that it has to be the nexus file as given by RevBayes
      pol.tree <- readTrees(paths = file)[[1]][[1]]

      # Lee los mapas estocasticos
      mapsfile <- "stochmap_mk2_polinizador_run_1.log" 

      #Resume los mapas estocasticos 
      stoch_map_df <- processStochMaps(pol.tree,
                                 mapsfile, 
                                 state_labels=c("Insecto" = "0",
                                                "Viento" = "1"),
                                 burnin = 0.1)

     mycolors= setNames(c("red","blue"),c("Insecto","Viento"))

      # Grafica los mapas estocasticos
      plotStochMaps(tree=pol.tree,
              maps=stoch_map_df, 
              color_by="MAP",
              colors = mycolors,
              tip_labels_size=1, 
              tree_layout = "rectangular",
              line_width = 0.8)
              
```
+ [Código completo de R](downloads/mk2_visualization.R)
