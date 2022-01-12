# Carregando Pacotes
library(BatchGetSymbols)
library(plyr)
library(dplyr)
library(xts)
library(PerformanceAnalytics)
library(Benchmarking)
library(nloptr)
library(ggplot2)

# Funções
sdp = function(x){
  return(as.vector(sqrt(t(as.matrix(x)) %*% covMat %*% as.matrix(x))*sqrt(252)))
}
rtp = function(x){
  return(1/as.vector((1+(x %*% mean_return))^252-1))
}
ig = function(x){
  return(sum(x)-1)
}

# Parâmetros
di = as.Date("01/01/2010",format = "%d/%m/%Y")
df = Sys.Date()
neg = 10000
lim = 0.85 
simulac = 5000

# Coletando Dados
ibov = c('MGLU3.SA','VIIA3.SA','B3SA3.SA','ITUB4.SA','CPLE6.SA','VALE3.SA','PETR4.SA','WEGE3.SA','ABEV3.SA','RADL3.SA')
ibov = unique(ibov)
ibov_i = c('^BVSP')
colect = c(ibov_i,ibov)
benchmark = '^BVSP'

dados_ibov = BatchGetSymbols(
  tickers = colect,
  first.date = di,
  last.date = df,
  thresh.bad.data = lim,
  bench.ticker = benchmark,)

dados_ibov = dados_ibov$df.tickers
dados_ibov = dlply(dados_ibov,.(ticker), function(x){rownames(x)=x$row;x$row = NULL;x})

data = dados_ibov[[1]][,c(7,4)]
vol = dados_ibov[[1]][,c(7,5)]
colnames(data) = c("Data", dados_ibov[[1]][1,8])
colnames(vol) = c("Data", dados_ibov[[1]][1,8])


for(i in 2:length(dados_ibov)){
  acao = dados_ibov[[i]][,c(7,4,5)]
  colnames(acao) = c("Data", dados_ibov[[i]][1,8],paste0(dados_ibov[[i]][1,8],'v'))
  data = merge(data, acao[,1:2], by = "Data")
  vol = merge(vol, acao[,c(1,3)], by = 'Data')
}
colnames(data) = gsub(x = colnames(data), pattern = '\\^', replacement = '')
colnames(vol) = gsub(x = colnames(vol), pattern = '\\^', replacement = '')

# Selecionando datas aleatórias
ref1 = data
datin = min(data$Data)
datmax = max(data$Data)-365
datini = datin + 365
datmaxi = datmax - 365


# Rodando simulações
results = data.frame()
for (t in 1:simulac) {
  print(paste0('Inciando simulação ', t))
  # Criando insumo para sorteio
  interfin = data.frame(data = seq.Date(from = datini, to = datmax, by = 'days'), 
                        sorteio = runif(n = length(seq.Date(from = datini, to = datmax, by = 'days'))))
  interini = data.frame(data = seq.Date(from = datin, to = datmaxi, by = 'days'), 
                        sorteio = runif(n = length(seq.Date(from = datin, to = datmaxi, by = 'days'))))
  
  # Sorteio
  print(paste0('Sorteando datas'))
  dini = interini[interini$sorteio == max(interini$sorteio),1]
  sel = interfin %>% filter(data > dini+60)
  dfin = sel[sel$sorteio == max(sel$sorteio), 1]
  dprev = dfin + 365
  
  # Seleção de dados no intervalo
  data = ref1 %>% filter(Data >= dini,
                         Data <= dfin)
  dados = xts(data[,-1],order.by = data$Data)
  
  # Calculando Estatísticas
  print(paste0('Calculando estatísticas'))
  BETA = CalculateReturns(dados, method = 'log')
  BETA = BETA[-1,]
  
  return = as.vector(Return.annualized(BETA, scale = 252))
  risk = as.vector(sd.annualized(BETA, scale = 252))
  REL = data.frame(acaon = colnames(BETA), return, risk)

  colnames(REL) = c("Ativo","Retorno","Risco")
  REL = REL[-1,]
  Dados = REL
    
  # Selecionando Pesos
  print(paste0('Calculando carteira'))
  ativos = as.vector(Dados$Ativo)
  returns = xts(data.frame(BETA) %>% select(ativos), order.by = data$Data[2:length(data$Data)])
    
  mean_return = colMeans(returns)
  covMat = cov(returns)
    
  xo = rep(1,dim(returns)[2])
  lw = rep(0,dim(returns)[2])
  up = rep(1,dim(returns)[2])
  
  # Selecionando pesos ideais, aleatórios e mercado
  shp = function(x){
    return(1/((1/sdp(x))/rtp(x)))
  }
  mshp = auglag(xo,shp,gr=NULL,lower = lw,upper = up,hin=NULL,heq = ig,localsolver = "MMA")
  ale = runif(10)
  ale = ale/sum(ale)
  market = vol %>% filter(Data >= dfin,
                         Data <= dprev)
  market = market[,-c(1,2)]
  market = colSums(market)/sum(colSums(market))
    
  # Calculando métricas futuras
  print(paste0('Calculando métricas'))
  prev = ref1 %>% filter(Data >= dfin-1,
                         Data <= dprev) %>% select(Data, ativos)
  prev_r = xts(prev[,-1],order.by = prev$Data)
  prev_r = CalculateReturns(prev_r)
  prev_r = prev_r[-1,]
    
  # Risco e Retorno anterior
  risco_a = as.vector(sd.annualized(as.matrix(returns) %*% mshp$par, scale = 252))
  retorno_a = as.vector(Return.annualized(as.matrix(returns) %*% mshp$par, scale = 252))
  retorno_a_a = as.vector(Return.annualized(as.matrix(returns) %*% ale, scale = 252))
  risco_a_a = as.vector(sd.annualized(as.matrix(returns) %*% ale, scale = 252))
  retorno_m_a = as.vector(Return.annualized(as.matrix(returns) %*% market, scale = 252))
  risco_m_a = as.vector(sd.annualized(as.matrix(returns) %*% market, scale = 252))
    
  data_inicial = min(data$Data)
  data_final = max(data$Data)
  data_max_prev = max(prev$Data)
  
  # Risco e Retorno Atual
  risco = as.vector(sd.annualized(as.matrix(prev_r) %*% mshp$par, scale = 252))
  retorno = as.vector(Return.annualized(as.matrix(prev_r) %*% mshp$par, scale = 252))
  risco_al = as.vector(sd.annualized(as.matrix(prev_r) %*% ale, scale = 252))
  retorno_al = as.vector(Return.annualized(as.matrix(prev_r) %*% ale, scale = 252))
  risco_m = as.vector(sd.annualized(as.matrix(prev_r) %*% market, scale = 252))
  retorno_m = as.vector(Return.annualized(as.matrix(prev_r) %*% market, scale = 252))
  
  
  nlin = data.frame(data_inicial, data_final, data_max_prev,
                    retorno_a, retorno_m_a, retorno_a_a, retorno, retorno_m, retorno_al,
                    risco_a, risco_m_a, risco_a_a, risco, risco_m, risco_al)
  results = rbind(results, nlin)
    
  # Verificação de Resultados preliminares
  print(median(results$retorno))
  print(median(results$retorno_m))
  print(median(results$retorno) - median(results$retorno_m))
}
rownames(results) = 1:nrow(results)

# Calculando dias de previsão
results$dias_prev = as.numeric(results$data_final-results$data_inicial)

# Referência com o mercado
results$difret = results$retorno - results$retorno_m
results$difret_a = results$retorno_a - results$retorno_m_a
results$difris = results$risco - results$risco_m
results$difris_a = results$risco_a - results$risco_m_a

# Removendo outliers
for (j in c(4,7)) {
  outlier = boxplot(results[,j], plot=FALSE)$out
  results = results[-which(results[,j] %in% outlier),]
}

# Criando repartição
results$r1 = as.factor(ifelse(results$difret_a<=0, 1, 2))
results$r2 = as.factor(ifelse(results$difret<=0, 1, 2))

# Resultados
# Gráficos e Descrição
ggplot(results) +
  aes(x = retorno_a) +
  geom_histogram(bins = 50L, fill = "#00B050") +
  labs(
    x = "Retorno",
    y = "Frequência",
    title = "Distribuição dos Retornos de Markowitz",
    subtitle = "Período de análise",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14L, face = "bold"))

summary(results$retorno_a)

ggplot(results) +
  aes(x = retorno_m_a) +
  geom_histogram(bins = 50L, fill = "#00B050") +
  labs(
    x = "Retorno",
    y = "Frequência",
    title = "Distribuição dos Retornos de Mercado",
    subtitle = "Período de análise",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14L, face = "bold"))

summary(results$retorno_m_a)

ggplot(results) +
  aes(x = difret_a, fill = r1) +
  geom_histogram(bins = 80L) +
  scale_fill_manual(
    values = c(`1` = "#C00000",
               `2` = "#00B050")
  ) +
  labs(
    x = "Diferença do Retorno",
    y = "Frequência",
    title = "Diferença entre Retorno de Markowitz e Mercado",
    subtitle = "Período de análise",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14L,
                              face = "bold")
  )

summary(results$difret_a)
quantile(results$difret_a, 0.05)
quantile(results$difret_a, 0.95)
sum(results$difret_a>0)/nrow(results)
wilcox.test(results$difret_a, alternative = 'less')


ggplot(results) +
  aes(x = retorno) +
  geom_histogram(bins = 50L, fill = "#00B050") +
  labs(
    x = "Retorno",
    y = "Frequência",
    title = "Distribuição dos Retornos de Markowitz",
    subtitle = "Período de previsão",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14L, face = "bold"))

summary(results$retorno)

ggplot(results) +
  aes(x = retorno_m) +
  geom_histogram(bins = 50L, fill = "#00B050") +
  labs(
    x = "Retorno",
    y = "Frequência",
    title = "Distribuição dos Retornos de Mercado",
    subtitle = "Período de previsão",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14L, face = "bold"))

summary(results$retorno_m)

ggplot(results) +
  aes(x = difret, fill = r2) +
  geom_histogram(bins = 80L) +
  scale_fill_manual(
    values = c(`1` = "#C00000",
               `2` = "#00B050")
  ) +
  labs(
    x = "Diferença do Retorno",
    y = "Frequência",
    title = "Diferença entre Retorno de Markowitz e Mercado",
    subtitle = "Período de previsão",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14L,
                              face = "bold")
  )

summary(results$difret)
quantile(results$difret, 0.05)
quantile(results$difret, 0.95)
sum(results$difret>0)/nrow(results)

results$difal = results$retorno_al - results$retorno_m
results$r3 = as.factor(ifelse(results$difal<=0, 1, 2))
ggplot(results) +
  aes(x = difal, fill = r3) +
  geom_histogram(bins = 80L) +
  scale_fill_manual(
    values = c(`1` = "#C00000",
               `2` = "#00B050")
  ) +
  labs(
    x = "Diferença do Retorno",
    y = "Frequência",
    title = "Diferença entre Retorno Aleatório e Mercado",
    subtitle = "Período de previsão",
    caption = "Fonte: Yahoo Finance. Elaboração: Evânio Marques"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14L,
                              face = "bold")
  )

summary(results$difal)
sum(results$difal<=0)/nrow(results)