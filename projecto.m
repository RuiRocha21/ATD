clear all;
close all;
clc;

nomeFicheiro = 'dataset_ATD_PL8.csv';

%abertura do ficheiro
idFicheiro = fopen(nomeFicheiro);
%ler ficheiro
dados = textscan(idFicheiro,'%s %f','HeaderLines',1,'Delimiter',';');
%extracao de dados
xDados = dados{1,1}(:,1);
yDados = dados{1,2}(:,1);

%plot dos valores lidos
dias = 1:365;
N = length(yDados);
t=0 : N-1;
figure(1);

plot(t,yDados);
xlabel('numero de dias');
ylabel('consumo de energia diario (KW/h)');
title('serie temporal original');
%------------------------------------------------------

%tratar de valores NaN fazendo extrapolaçao de dados
if any(isnan(yDados))
    indice = find(isnan(yDados));
    valoresSerieTemporal = yDados;
    for k = 1 :length(indice)
        ttempos = t(indice(k)-4:indice(k)-1);
        valores = valoresSerieTemporal(indice(k)-4:indice(k)-1);
        valoresSerieTemporal(indice(k))=interp1(ttempos, valores, t(indice(k)), 'pchip', 'extrap');
    end
end
figure(2);
plot(t,valoresSerieTemporal);
xlabel('numero de dias');
ylabel('consumo de energia diario (KW/h)');
title('serie temporal sem valores NaN');
%----------------------------------------------------

%tratar outliers
media = mean(valoresSerieTemporal);
desvioPadrao = std(valoresSerieTemporal);
media_Atipica = repmat(media,N,1);
desvioPadrao_Atipico = repmat(desvioPadrao,N,1);
indice = find(abs(valoresSerieTemporal-media_Atipica) > 3*desvioPadrao);
valoresSerieTemporalCopia = valoresSerieTemporal;

if length(indice) > 1
    for k = 1 : length(indice)
        if valoresSerieTemporalCopia(indice(k))> media
            valoresSerieTemporalCopia(indice(k)) = media + 2.5*desvioPadrao;
        else
            valoresSerieTemporalCopia(indice(k)) = media - 2.5*desvioPadrao;
        end
    end
end

figure(3);
plot(t, valoresSerieTemporalCopia);
legend('Serie temporal sem outliers');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('serie temporal sem sem outliers');
%-------------------------------------------------------------
%Identificação das componentes da série temporal
%serie temporal = tendencia + sazonal + irregular
%separar a série em 3 partes e encontrar qual parte é cada modelo

%detrend ->retira a tendência
tamanho = length(valoresSerieTemporalCopia);
novo_t = 0:tamanho-1;
tendenciaConstante = detrend(valoresSerieTemporalCopia, 'constant');
tendenciaLinear = detrend(valoresSerieTemporalCopia, 'linear'); %aproximação por uma reta

figure(4);
subplot(211);
plot(novo_t, valoresSerieTemporalCopia, '-+', novo_t, valoresSerieTemporalCopia-tendenciaConstante, '-*') %valoresSerieTemporalCopia-tendenciaConstante -> tendência
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal com tendência grau zero');
subplot(212)
plot(novo_t, tendenciaConstante, '-o');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal sem tendência grau zero');
%---------------------------------------------------------------------
%obter o coeficientes do polinómio da serie

polinomio= polyfit(novo_t', valoresSerieTemporalCopia, 2); %ultimo parametro -> grau
%tendencia
tendencia= polyval(polinomio, novo_t');

figure(5);
subplot(211);
plot(novo_t, valoresSerieTemporalCopia, '-+', novo_t, tendencia, '-*');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal com tendência grau dois');
subplot(212);
plot(novo_t, tendenciaConstante - tendencia, '-o'); %tendenciaConstante - tendencia ->série sem tendência
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal sem tendência grau dois');

%---------------------------------------------------------------------
%sazonalidade -> ao longo do gráfica existir "ondas" que se adaptam a uma função sinusoidal
%Estimar a componente da sazonalidade da série temporal

diferenca = valoresSerieTemporalCopia - tendencia;
h0= (1:365)';
dv= repmat(h0, 1,1);
sazonalidade= dummyvar(dv); %matriz c 365 linhas e 365 colunas, diagonais a zeros

xTemp = sazonalidade\diferenca; %divisao a esquerda
sazonalidadeSerieTemp = sazonalidade*xTemp; %componente sazonal

figure(6);
subplot(211);
plot(novo_t, valoresSerieTemporalCopia, '-+', novo_t, valoresSerieTemporalCopia - sazonalidadeSerieTemp,'-*', novo_t, tendencia, '-o');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal sem sazonalidade');
subplot(212);
plot(novo_t, sazonalidadeSerieTemp, '-o');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Sazonalidade da série temporal');

%------------------------------------------------------------------
%componentes irregular da serie
%irregular -> se for significativa tenta-se arranjar um modelo para ela, se for menosprezável pode ser ignorada

componenteIrregular = valoresSerieTemporalCopia - sazonalidadeSerieTemp;

figure(7);
subplot(211);
plot(novo_t, valoresSerieTemporalCopia, '-+', novo_t, valoresSerieTemporalCopia-componenteIrregular,'-*');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal sem componente irregular');
subplot(212);
plot(novo_t,componenteIrregular, '-o');
xlabel('t[dias]');
ylabel('consumo de energia diario (KW/h)');
title('Componente irregualar da série temporal');

%-------------------------------------------------------------------
%modelo AR da serie

estacionaridadeDaserie = adftest(sazonalidadeSerieTemp);
FAC = autocorr(sazonalidadeSerieTemp);
FACP = parcorr(sazonalidadeSerieTemp);

Ts = 1;
objectoData = iddata(sazonalidadeSerieTemp,[],Ts,'TimeUnit','days');
figure(8); 
subplot(2,1,1);% daqui se tira valor 15
autocorr(sazonalidadeSerieTemp); %autocorrelaçao da serie
subplot(2,1,2);
parcorr(sazonalidadeSerieTemp); %autocorrelaçao parcial da serie
na = 15;
opt = arOptions('Approach','ls','Window','ppw');
modelo = ar(objectoData,na,'ls');
[A] = polydata(modelo);
optForecast = forecastOptions('InitialCondition','z');
figure(9);
ARforecast = forecast(modelo,sazonalidadeSerieTemp(1:na),30-na,optForecast);
ARForecast_Final=repmat([sazonalidadeSerieTemp(1:na); ARforecast],12,1);
plot(1:360,ARForecast_Final);
xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Previsão com o modelo AR (-o)')

previsaoAR = repmat(ARForecast_Final,2,1);
dobroDaSerie = repmat(valoresSerieTemporalCopia(1:360),2,1);
figure(16);
plot(1:720,dobroDaSerie,'-+',1:720,previsaoAR,'-o');

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal (-+) e Previsão com o modelo AR (-o)')

%------------------modelo ARMA
na = 15;
nc = 1;
opcoes  =armaxOptions('SearchMethod', 'auto');
modeloArma = armax(objectoData,[na nc],opcoes)
[pa_ARMA,pb_ARMA,pc_ARMA] = polydata(modeloArma);

ruidoBranco = randn(30,1); %ruido branco

ARMA=valoresSerieTemporalCopia(1:na);  % Simulação do Modelo ARMA

for k=na+1:30,ARMA(k)=sum(-pa_ARMA(2:end)'.*flip(ARMA(k-na:k-1)))+sum(pc_ARMA'.*flip(ruidoBranco(k-nc:k))); 

end

ARMA_fim=repmat(ARMA,12,1);

%simulacao do modelo ARMA com forecast

armaForecast=forecast(modeloArma,valoresSerieTemporalCopia(1:na),30-na);
armaForecast_Final=repmat([valoresSerieTemporalCopia(1:na); armaForecast],12,1);

figure(10);

plot(1:360,ARMA_fim,'-*',1:360,armaForecast_Final);

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Componente sazonal 1 (-+) e Estimação com o modelo ARMA (-o)')


figure(11)  % Compara a série com o modelo ARMA + tendência

plot(1:360,valoresSerieTemporalCopia(1:360),'-+',1:360,ARMA_fim+tendencia(1:360),'-o');

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal (-+) e Estimação com o modelo ARMA (-o)')


% dados de analise do modelo

estimativaModeloARMA = sum((valoresSerieTemporalCopia(1:360)-ARMA_fim(1:360)).^2);

figure(12)  % Faz a previsão para o dobro da serie

previsaoARMA = repmat(ARMA_fim,2,1)+repmat(tendencia(1:360),2,1);
dobroDaSerie = repmat(valoresSerieTemporalCopia(1:360),2,1);
plot(1:720,dobroDaSerie,'-+',1:720,previsaoARMA,'-o');

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal (-+) e Previsão com o modelo ARMA (-o)')

%-----------------MODELO ARIMA-------------------------------------
paramentro_P=30;  % Parâmetros do modelo ARIMA

parametro_D=1;

parametro_Q=1;

% estrutura do modelo ARIMA

modeloARIMA=arima(paramentro_P,parametro_D,parametro_Q);

% Estimação do modelo ARIMA

estimacaoARIMA = estimate(modeloARIMA,valoresSerieTemporalCopia(1:360),'Y0',valoresSerieTemporalCopia(1:paramentro_P+1));

% Simulação do modelo ARIMA

simulacaoARIMA = simulate(estimacaoARIMA,360);

figure(14) %serie temporal com a sua estimação

plot(1:360,valoresSerieTemporalCopia(1:360),'-+',1:360,simulacaoARIMA,'-o');

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal (-+) e Estimação do o modelo ARIMA (-o)')

% análise do modelo

analiseModeloArima=sum((valoresSerieTemporalCopia(1:360)-simulacaoARIMA(1:360)).^2);

% Simulação do modelo ARIMA para 2N

simulacaoDrobroSerie = simulate(estimacaoARIMA,720);

figure(15)  % Faz a previsão para 2N

plot(1:720,dobroDaSerie,'-+',1:720,simulacaoDrobroSerie,'-o');

xlabel('t [dias]');
ylabel('consumo de energia diario (KW/h)');
title('Série temporal (-+) e Previsão com o modelo ARIMA (-o)')