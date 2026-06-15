# =====================================================================
#  A Six-Week Heat-Accumulation Window for Weekly EHEC Notifications
#  in Korea, 2005-2024 — REPRODUCIBLE ANALYSIS CODE (v2, 2026-06-15)
#  Seongdae Kim, Byung Chul Chun
# ---------------------------------------------------------------------
#  v2 supersedes v1: residual autocorrelation is modelled CORRECTLY with
#  gamm()+corAR1 (+ a quarterly harmonic). NOTE: bam(family=nb(), rho=)
#  silently IGNORES AR(1) for generalized models — the corrected
#  single-lag and cumulative estimates below come from gamm()+corAR1 and
#  match the EHPM manuscript.
#
#  R 4.5.x | dplyr, mgcv, nlme, MASS, dlnm, splines | License: MIT
#  Data (no PII): KDCA weekly EHEC counts; KMA ASOS climate; AirKorea PM10.
#
#  Reproduces (autocorrelation-adjusted, corAR1 + quarterly):
#    single-lag temp lag3 IRR 1.05 (p~0.010); cumulative 5wk +9.2% (sig),
#    6wk +8.6% (p~0.055, borderline peak).
#
#  OUTPUT MAP:
#    Part 2 Table 1 (descriptives+Spearman) | Part 3 Table 2 (single-lag)
#    Part 4 Table 3 + Table S2 (temp profile) | Part 5 S1/S3/S4 + family + concurvity + rolling-CV
#    Part 6 autocorrelation diagnostic | Figures via companion figure script
# =====================================================================
suppressMessages({library(dplyr); library(mgcv); library(nlme); library(MASS); library(dlnm); library(splines)})

## 0. Config
BASE_IV <- "PATH/TO/FBD_DATA_ZIP"      # <-- set to data directory (see README)
DISEASE <- "장출혈성대장균감염증"; Y0 <- 2005; Y1 <- 2024; MAXL <- 8
ctrl <- lmeControl(opt="optim", maxIter=200, msMaxIter=200, returnObject=TRUE)
adj_terms <- "s(time_idx,k=20,fx=TRUE)+sin52+cos52+sin26+cos26+sin13+cos13"   # time spline + Fourier(52,26)+quarterly(13)
varbase <- list(temp="avg_temp", humidity="humidity", pm10="pm10", precipitation="precipitation")
others  <- list(
  temp          = c("s(humidity,k=6)","s(pm10,k=6)","s(precipitation_lag4,k=6)"),
  humidity      = c("s(avg_temp_lag3,k=6)","s(pm10,k=6)","s(precipitation_lag4,k=6)"),
  pm10          = c("s(avg_temp_lag3,k=6)","s(humidity,k=6)","s(precipitation_lag4,k=6)"),
  precipitation = c("s(avg_temp_lag3,k=6)","s(humidity,k=6)","s(pm10,k=6)"))

## 1. Load & preprocess (weekly, 2005-2024)
dv <- read.csv(file.path(BASE_IV,"식중독최종.csv"), fileEncoding="UTF-8-BOM", check.names=FALSE, stringsAsFactors=FALSE)
food <- dv %>% filter(disease==DISEASE, year>=Y0, year<=Y1) %>% group_by(year,week) %>%
  summarise(cases=sum(cases,na.rm=TRUE), pop=sum(population,na.rm=TRUE)/n_distinct(region), .groups="drop") %>%
  filter(!is.na(week), week>=1, week<=52) %>% arrange(year,week)
wx <- read.csv(file.path(BASE_IV,"기상데이터_전처리.csv"), fileEncoding="UTF-8", check.names=FALSE, stringsAsFactors=FALSE)
wxw <- wx %>% filter(year>=Y0,year<=Y1) %>% group_by(year,week) %>%
  summarise(avg_temp=mean(`평균기온(°C)`,na.rm=TRUE), humidity=mean(`평균 상대습도(%)`,na.rm=TRUE),
            wind_speed=mean(`평균 풍속(m/s)`,na.rm=TRUE), precipitation=mean(`일강수량(mm)`,na.rm=TRUE), .groups="drop") %>%
  filter(!is.na(week), week>=1, week<=52)
hs <- read.csv(file.path(BASE_IV,"황사_전처리.csv"), fileEncoding="UTF-8", check.names=FALSE, stringsAsFactors=FALSE)
hsw <- hs %>% filter(year>=Y0,year<=Y1) %>% group_by(year,week) %>%
  summarise(pm10=mean(`일 미세먼지 농도(㎍/㎥)`,na.rm=TRUE), .groups="drop") %>% filter(!is.na(week), week>=1, week<=52)
ts <- food %>% left_join(wxw,c("year","week")) %>% left_join(hsw,c("year","week")) %>% arrange(year,week) %>%
  mutate(time_idx=row_number(), sin52=sin(2*pi*week/52),cos52=cos(2*pi*week/52),
         sin26=sin(2*pi*week/26),cos26=cos(2*pi*week/26), sin13=sin(2*pi*week/13),cos13=cos(2*pi*week/13),
         period=factor(ifelse(year>=2020,"Post","Pre"),levels=c("Pre","Post")))
for(L in 0:MAXL) for(v in c("avg_temp","humidity","pm10","precipitation","wind_speed"))
  ts[[paste0(v,"_lag",L)]] <- dplyr::lag(ts[[v]], L)
base_ok <- ts %>% filter(!is.na(cases),!is.na(avg_temp_lag3),!is.na(humidity),!is.na(pm10),!is.na(precipitation_lag4)) %>% arrange(time_idx)
cat(sprintf("[Data] full weeks=%d | analytic(base_ok)=%d\n", nrow(ts), nrow(base_ok)))

## 2. Table 1 + Spearman  (computed on the analytic set, n = base_ok, to match the manuscript)
cat("\n[Table 1]\n")
for(v in c("cases","avg_temp","humidity","precipitation","wind_speed","pm10")){
  m<-mean(base_ok[[v]],na.rm=T); s<-sd(base_ok[[v]],na.rm=T)
  rho<-if(v=="cases") NA else suppressWarnings(cor(base_ok$cases,base_ok[[v]],method="spearman",use="complete.obs"))
  cat(sprintf("  %-14s %.2f ± %.2f | Spearman %s\n", v, m, s, ifelse(is.na(rho),"-",sprintf("%+.3f",rho))))}

## NB theta (fixed in gamm for stability)
TH <- bam(cases~s(avg_temp_lag3,k=6)+s(humidity,k=6)+s(pm10,k=6)+s(precipitation_lag4,k=6)+s(time_idx,k=29)+sin52+cos52+sin26+cos26,
          family=nb(), data=base_ok, method="fREML")$family$getTheta(TRUE)
cat(sprintf("\n[theta] %.3f\n", TH))

## 3. Table 2 — single-lag IRR (gamm + corAR1 + quarterly); focal linear, others smooth
single_lag <- function(vkey,L){
  vb<-varbase[[vkey]]; fv<-if(L==0) vb else paste0(vb,"_lag",L)
  dd<-base_ok %>% filter(!is.na(.data[[fv]])) %>% arrange(time_idx)
  fm<-as.formula(paste0("cases ~ ",fv,"+",paste(others[[vkey]],collapse="+"),"+",adj_terms))
  g<-tryCatch(gamm(fm,family=negbin(TH),data=dd,correlation=corAR1(form=~time_idx),control=ctrl),error=function(e)NULL)
  if(is.null(g)) return(NULL); co<-summary(g$gam)$p.table; if(!fv%in%rownames(co)) return(NULL)
  b<-co[fv,1]; se<-co[fv,2]
  data.frame(var=vkey,lag=L,IRR=round(exp(b),3),lo=round(exp(b-1.96*se),3),hi=round(exp(b+1.96*se),3),p=signif(co[fv,4],3))
}
cat("\n[Table 2] single-lag IRR (corAR1+quarterly)\n")
T2 <- do.call(rbind, lapply(c("temp","humidity","pm10","precipitation"), function(v) do.call(rbind, lapply(0:8, function(L) single_lag(v,L)))))
print(T2, row.names=FALSE)

## 4. Table 3 + S2 — cumulative via distributed-lag-LINEAR in gamm (sum of lag coefs; delta-method CI)
cum_gamm <- function(vkey, maxlag=8, winsor=FALSE, offset=FALSE){
  vb<-varbase[[vkey]]; tl<-paste0(vb,"_lag",0:maxlag); tl[1]<-vb
  needcol<-paste0(vb,"_lag",maxlag)
  dd<-base_ok %>% filter(!is.na(.data[[needcol]]),!is.na(pop)) %>% arrange(time_idx)
  if(winsor) dd$cases<-pmin(dd$cases, as.numeric(quantile(dd$cases,0.99,na.rm=TRUE)))
  rhs<-paste(c(tl, others[[vkey]]),collapse="+"); if(offset) rhs<-paste0(rhs,"+offset(log(pop))")
  g<-gamm(as.formula(paste0("cases ~ ",rhs,"+",adj_terms)),family=negbin(TH),data=dd,correlation=corAR1(form=~time_idx),control=ctrl)
  co<-summary(g$gam)$p.table; V<-vcov(g$gam)
  do.call(rbind, lapply(0:maxlag, function(K){
    idx<-tl[1:(K+1)]; idx<-idx[idx%in%rownames(co)]; b<-sum(co[idx,1]); se<-sqrt(sum(V[idx,idx]))
    data.frame(endlag=K, pct=round((exp(b)-1)*100,2), lo=round((exp(b-1.96*se)-1)*100,2),
               hi=round((exp(b+1.96*se)-1)*100,2), p=signif(2*pnorm(-abs(b/se)),3))}))
}
cat("\n[Table S2] temperature single-lag & cumulative profile (corAR1+quarterly)\n")
Ttemp <- cum_gamm("temp"); print(Ttemp, row.names=FALSE)
cat("\n[Table 3] best cumulative-lag per variable\n")
for(v in c("temp","humidity","pm10","precipitation")){
  cp<-cum_gamm(v); sig<-cp[cp$lo>0|cp$hi<0,]; best<-if(nrow(sig)) sig[which.max(abs(sig$pct)),] else cp[which.max(abs(cp$pct)),]
  cat(sprintf("  %-13s best lag %d: %+.2f%% (%.2f,%.2f) p=%.3f\n", v, best$endlag,best$pct,best$lo,best$hi,best$p))}

## 5. Sensitivity / robustness
cat("\n[Table S1] sensitivity sets (cumulative temp, corAR1+quarterly)\n")
for(win in c("LONG","SHORT")) for(tr in c("Original","De-spiked")){
  K<-if(win=="LONG") 8 else 4; cp<-cum_gamm("temp", maxlag=K, winsor=(tr=="De-spiked"))
  sig<-cp[cp$lo>0|cp$hi<0,]; best<-if(nrow(sig)) sig[which.max(abs(sig$pct)),] else cp[which.max(abs(cp$pct)),]
  cat(sprintf("  %-5s/%-9s best lag %d: %+.2f%% (%.2f,%.2f) %s\n", win,tr,best$endlag,best$pct,best$lo,best$hi, ifelse(best$lo>0|best$hi<0,"SIG","ns")))}

cat("\n[Table S4] population-offset sensitivity (6-wk cumulative temp)\n")
o0<-cum_gamm("temp",maxlag=6,offset=FALSE); o1<-cum_gamm("temp",maxlag=6,offset=TRUE)
cat(sprintf("  no offset %.2f%% | +offset %.2f%% (delta %.3f pp)\n", o0$pct[7], o1$pct[7], o1$pct[7]-o0$pct[7]))

cat("\n[Table S3] Pre vs Post single-lag-3 temp\n")
for(pp in c("Pre","Post")){
  dd<-base_ok %>% filter(period==pp) %>% arrange(time_idx)
  g<-gamm(as.formula(paste0("cases ~ avg_temp_lag3+",paste(others$temp,collapse="+"),"+",adj_terms)),
          family=negbin(TH),data=dd,correlation=corAR1(form=~time_idx),control=ctrl)
  co<-summary(g$gam)$p.table["avg_temp_lag3",]
  cat(sprintf("  %-4s N=%d IRR=%.3f (%.3f,%.3f) p=%.3f\n", pp,nrow(dd),exp(co[1]),exp(co[1]-1.96*co[2]),exp(co[1]+1.96*co[2]),co[4]))}

cat("\n[Family] Poisson vs NB\n")
mp<-glm(cases~avg_temp_lag3+humidity+pm10+precipitation_lag4+sin52+cos52, family=poisson, data=base_ok)
mn<-glm.nb(cases~avg_temp_lag3+humidity+pm10+precipitation_lag4+sin52+cos52, data=base_ok)
cat(sprintf("  Poisson AIC=%.1f | NB AIC=%.1f (dAIC %.1f, NB preferred)\n", AIC(mp),AIC(mn),AIC(mp)-AIC(mn)))

cat("\n[Table S5/S6] worst-case concurvity (full-smooth diagnostic spec)\n")
dc<-base_ok %>% filter(!is.na(wind_speed_lag3))
mc<-bam(cases~s(time_idx,k=29)+s(week,bs="cc",k=24)+s(avg_temp_lag3,k=6)+s(wind_speed_lag3,k=6), family=nb(), data=dc, method="fREML")
print(round(mgcv::concurvity(mc, full=FALSE)$worst,3))

## 6. AUTOCORRELATION DIAGNOSTIC — proper NB+AR(1) via gamm+corAR1 (+quarterly)
## (bam(nb,rho=) would silently ignore AR(1); we verify residual whitening here)
cat("\n[Autocorrelation] gamm+corAR1+quarterly, normalized-residual Ljung-Box\n")
g_ar<-gamm(cases~s(avg_temp_lag3,k=6)+s(humidity,k=6)+s(pm10,k=6)+s(precipitation_lag4,k=6)+s(time_idx,k=20,fx=TRUE)+sin52+cos52+sin26+cos26+sin13+cos13,
           family=negbin(TH),data=base_ok,correlation=corAR1(form=~time_idx),control=ctrl)
r<-residuals(g_ar$lme,type="normalized"); phi<-as.numeric(coef(g_ar$lme$modelStruct$corStruct,unconstrained=FALSE))
cat(sprintf("  phi_hat=%.3f | Ljung p lag12=%.3f lag26=%.3f lag52=%.3f\n", phi,
            Box.test(r,12,"Ljung-Box")$p.value, Box.test(r,26,"Ljung-Box")$p.value, Box.test(r,52,"Ljung-Box")$p.value))

cat("\n[DONE] v2 reproducible pipeline complete.\n")
# Figures (Fig1 time series; Fig2 single-lag/cumulative/forest; S1 forest; S2 smooths; S3 heatmap)
# are produced by the companion figure script; Fig S4 (data-flow) by a separate diagram script. See README.
