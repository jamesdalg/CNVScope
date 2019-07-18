## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(CNVScope)
library(pwr)
library(magrittr)

## ---- knit='asis'--------------------------------------------------------
large.effect.size<-pwr::cohen.ES(test="f2",size="large")$effect.size
large.effect.size
f2.res<-pwr::pwr.f2.test(u = 1, f2 = large.effect.size/(1 - large.effect.size), sig.level = 0.05, power = 0.8)
f2.res
n<-ceiling(f2.res$v+f2.res$u+1)
n

## ------------------------------------------------------------------------
pwr::pwr.r.test(r = 0.3, sig.level = 0.01, power = 0.8, alternative = "greater")$n %>% ceiling()

