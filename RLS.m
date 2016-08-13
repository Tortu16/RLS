%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         Proyecto 2
%
%                        Luis Jimenez
%
%   Cancelacion activa de ruido auditivo para mejora en la 
%   captura de voz utilizando un algoritmo adaptativo RLS.   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
file = csvread('C:\Users\luisd\Dropbox\Maestria\Procesamiento Adaptativo\Proyecto\Adaptive Noise Canceling\FanNoise2.txt',1);
%file = csvread('C:\Users\ljimenez\Documents\MATLAB\FanNoise2.txt',1);

% Longitud de la senal en muestras
N = size(file);
N = N(1);

% Varianzas de dos procesos estocasticos artificiales
var1 = 0.27;
var2 = 0.1;

% Generacion de los procesos estocasticos artificiales
%va1 = 2*file(1:N,1);
va1 = sqrt(var1)*randn(N,1);
va2 = sqrt(var2)*randn(N,1);

% Calculo de la varianza del ruido correlacionado d(n)
% Esta senal es el resultado del efecto de un filtro IIR sobre ruido blanco
% va1

vard =  var1/(1-0.8458*0.8458);

% Construye la senal de ruido correlacionado al de referencia

d(N) = 0;
d(1) = va1(1);
for k = 2:N
    d(k) = va1(k) - 0.8458*d(k-1);
end

% Construye la senal de entrada del filtro u = senal + v1

dn(N)=0;
for n=1:N;
dn(n)=sin((0.03*pi).*n+0.2*pi);
end;

u = dn' + va1;

Rsignal = 2*file(:,1);
Rnoise = 2*file(:,2);

Asignal = u;
Anoise = d;

varR = var(Rsignal);

%%%%%%%%%%%%%%      Filtro      %%%%%%%%%%%%%%
% d => signal
% x => noise

filterOrder = 16;

delta = 0.001*varR;

P0 = (1/delta)*eye(filterOrder);
Pn = P0;

factorOlvido = 0.3;
lambdaI = 1/factorOlvido;

phi = zeros(filterOrder,filterOrder);
cita = zeros(filterOrder,1);
Kn1 = zeros(filterOrder,1);
an1(N) = 0; % alfa(n+1)

AwRLS = zeros(filterOrder,1);
RwRLS = zeros(filterOrder,1);

Ae(N) = 0;
Re(N) = 0;

Awtrls(N,:) = AwRLS;
Rwtrls(N,:) = RwRLS;

AtempRLS = zeros(filterOrder,1);
RtempRLS = zeros(filterOrder,1);

AJrls (N,:) = AtempRLS;
RJrls (N,:) = RtempRLS;

for k = filterOrder:N-1
    % Xn+1 = temp = [noise(k+1); noise(k);noise(k-1); ... ; noise(k-m+2)]
    for m = 1:filterOrder
        AtempRLS(m) = Anoise(k-m+2);
        RtempRLS(m) = Rnoise(k-m+2);
    end
    
    %%%%%%%%%%%%%%      RLS      %%%%%%%%%%%%%%
    
    % Opera en la senal Artificial
    Kn1 = (lambdaI*Pn*AtempRLS)/(1+lambdaI*AtempRLS'*Pn*AtempRLS);
    an1(k+1) = Asignal(k+1) - AwRLS'*AtempRLS;
    AwRLS = AwRLS + Kn1*an1(k+1);
    Pn = lambdaI*(Pn-Kn1*AtempRLS'*Pn);
    
    
    %%%%%%%%%%%%%%      Evolucion de coeficientes      %%%%%%%%%%%%%%
    
    Awtrls(k,:) = AwRLS;
    
    %%%%%%%%%%%%%%      Funcion de Coste: J      %%%%%%%%%%%%%%
    %AJrls(k,:) = -Ae(k)*AtempRLS;

end

%%%%%%%%%%%%%%      Presentacion de resultados      %%%%%%%%%%%%%%

AwRLS;
RwRLS;

Asnr = [snr(Asignal, 51100, 6) snr(Ae, 51100, 6)]
Athd = [thd(Asignal) thd(Ae)]
Rsnr = [snr(Rsignal, 51100, 6) snr(Re, 51100, 6)]
Rthd = [thd(Rsignal) thd(Re)]

%%%%%%%%%%%%%%      Graficacion de la senal      %%%%%%%%%%%%%%

Afigure = figure('Name','Senales Simuladas');
subplot(2,2,1)
plot(Anoise(N/2:floor(2*N/3)),'b')
title('Senal de ruido artificial');
subplot(2,2,2)
plot(Asignal(N/2:floor(2*N/3)),'r')
title('Senal sinusoidal ruidosa artificial');
subplot(2,2,3)
plot(dn(N/2:floor(2*N/3)),'g')
title('Senal de sinusoidal artificial');
subplot(2,2,4)
plot(Ae(N/2:floor(2*N/3)),'k')
title('Senal de error artificial');


Rfigure = figure('Name','Senales Reales');
subplot(2,2,1)
plot(Rnoise(N/2:floor(2*N/3)),'b')
title('Senal de ruido adquirida');
subplot(2,2,2)
plot(Rsignal(N/2:floor(2*N/3)),'r')
title('Senal de voz ruidosa adquirida');
subplot(2,2,3)
plot(Re(N/2:floor(2*N/3)),'k')
title('Senal de error adquirida');

wfigure = figure('Name', 'Evolucion de coeficientes');
subplot(2,2,1);
plot(Awtrls(:,1));
hold
plot(Awtrls(:,2));
plot(Awtrls(:,3));
plot(Awtrls(:,4));
hold
title('Evolucion de w simulado');
subplot(2,2,2);
plot(Rwtrls(:,1));
hold
plot(Rwtrls(:,2));
plot(Rwtrls(:,3));
plot(Rwtrls(:,4));
hold
title('Evolucion de w adquirido');
subplot(2,2,3);
plot(Awtrls(1:50000,1));
title('w[1] primer segundo');
subplot(2,2,4);
plot(Rwtrls(1:50000,1));
title('w[1] primer segundo');

Jfigure = figure('Name', 'Funcion de coste');
subplot(2,1,1);
plot(AJrls(1:50000,1));
title('Funcion de coste para w[1] simulado');
subplot(2,1,2);
plot(RJrls(1:50000,1));
title('Funcion de coste para w[1] adquirido');
%sound(Rsignal(1:248574),51100)
%sound(Re(1:248574),51100)