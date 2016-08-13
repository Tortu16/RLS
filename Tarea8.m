clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         Tarea 8
%
%                       Luis Jimenez
%
% Basado en ejercico 2.11 del Haykin
%
%   - Filtro de Wiener
%
%   - Algoritmo de Descenso Maximo
%
%   - Algoritmo de Gradiente Estocastico LMS
%
%   - Variante de algoritmo de Gradiente Estocastico LMS
%
%   - RLS
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Definicion de variables base

N = 1000; % Tamano de los vectores de senales

% Varianzas de los procesos estocasticos involucrados
var1 = 0.27;
var2 = 0.1;

% Generacion de los procesos estocasticos involucrados
v1 = sqrt(var1)*randn(N,1);
v2 = sqrt(var2)*randn(N,1);

% Calculo de la varianza de la senal deseada d(n)
% Esta senal es el resultado del efecto de un filtro IIR sobre ruido blanco
% v1

vard =  var1/(1-0.8458*0.8458);

% Construye la senal deseada

d(1) = v1(1);
for k = 2:N
    d(k) = v1(k) - 0.8458*d(k-1);
end

% Construye la senal de entrada del filtro u = x + v2

x(1) = d(1);
for k = 2:N
    x(k) = d(k) + 0.9458*x(k-1);
end

x=x';
u = x + v2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    Filtro Wiener
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% E{u(n)u'(n)} = Ru = E{[x(n)+v2(n)][x'(n)+v2'(n)]} = Rx + Rv2
% Se asume que x(n) y v2(n) son procesos aleatorios independientes de media cero.
% rx(m) + 0.1rx(m-1) + 0.8rx(m-2) = rv1x = E{v1(n)x(n-m)} = 0

% Solucion de Ecuaciones Yule-Walker para m = 1 y m = 2 nos permite
% calcular r0 y r1:
r0 = 0.27/(1+(-0.01/1.8)+(0.008/1.8)+(-0.0016/1.8))
r1 = 0.27/(1+(-0.01/1.8)+(0.008/1.8)+(-0.0016/1.8))*(-0.1/1.8)

% Con estos datos calculamos matriz R y vector p

R = [r0 r1;r1 r0]+[0.1 0;0 0.1]

p = [r0-0.9458*r1; r1-0.9458*r0]

% Con estos datos calculamos los coeficientes del filtro de Wiener:

wopt = inv(R)*p

% Calculamos el valor de la funcion de coste minima Jmin:

Jmin = vard - p'*inv(R)*p

% Graficamos la superficie de la funcion de coste

[w0, w1] = meshgrid(-2:0.1:2);
J = vard + (w0.^2 + w1.^2)*r0 + 2*w0.*w1.*r1-2*w0*r0-2*w1*r1;

figure('Name','Superficie MSE','NumberTitle','off')
mesh (w0,w1,J)
figure('Name','Grafica de contornos MSE','NumberTitle','off')
contour(w0,w1,J);
hold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    Algoritmo Descenso Maximo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inicializamos un primer vector con nuestra suposicion inicial.
wDM = [-2; -2];

% Se elige un factor de aprendizaje arbitrario y pequeno.
mu = 0.05;

% Variable auxiliar para guardar la evolucion de nuestra estimacion.
wt = wDM;
% contour(w0,w1,J);
% hold
for k = 1:N
    wDM = wDM + mu*(p-R*wDM);
    wt = [wt wDM];
    scatter(wDM(1),wDM(2),1,'g');
end

% El vector con los valores estimados finales del filtro 
wDM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  Algoritmo Gradiente Estocastico LMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wLMS = [-2;-2];

mu = 0.01;

e(N) = 0;

e(1) = d(1) - wLMS'*[u(1); 0];
wLMS = wLMS + mu*e(1)*[x(1);0];
wt2 = wLMS;

for k = 2:N
    e(k) = d(k) - wLMS'*[u(k); u(k-1)];
    wLMS = wLMS + mu*e(k)*[x(k);x(k-1)];
    wt2 = [wt2 wLMS];
    scatter(wLMS(1),wLMS(2),1,'r');
end

wLMS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  Algoritmo Gradiente Estocastico LMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wLMA = [-2;-2];

mu = 0.01;

e(N) = 0;

e(1) = d(1) - wLMA'*[u(1); 0];
wLMA = wLMA + mu*e(1)*[x(1);0];
wt3 = wLMA;

for k = 2:N
    e(k) = d(k) - wLMA'*[u(k); u(k-1)];
    wLMA = wLMA + mu*sign(e(k))*[x(k);x(k-1)];
    wt3 = [wt3 wLMA];
    %scatter(wLMA(1),wLMA(2),1,'b');
end

wLMA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          Algoritmo RLS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wRLS = [-2;-2];

factorOlvido = 0.03;
lambdaI = 1/factorOlvido;

delta = 0.001*var(x);

P0 = (1/delta)*eye(2);
Pn = P0;

phi = zeros(2,2);
cita = zeros(2,1);
Kn1 = zeros(2,1);
an1(N) = 0; % alfa(n+1) error a priori
e(N) = 0;   % error a posteriori
tempRLS = zeros(2,1);
wt4 = wRLS;

for k = 2:N-1
    % Xn+1 = temp = [noise(k+1); noise(k);noise(k-1); ... ; noise(k-m+2)]
    for m = 1:2
        tempRLS(m) = x(k-m+2);
    end
    
    %%%%%%%%%%%%%%      RLS      %%%%%%%%%%%%%%
    
    % Opera en la senal Artificial
    Kn1 = (lambdaI*Pn*tempRLS)/(1+lambdaI*tempRLS'*Pn*tempRLS);
    an1(k+1) = d(k+1) - wRLS'*tempRLS;
    wRLS = wRLS + Kn1*an1(k+1);
    Pn = lambdaI*(Pn-Kn1*tempRLS'*Pn);
    wt4 = [wt4 wRLS];
    scatter(wRLS(1),wRLS(2),1,'b');
end

wRLS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              Graficas con los resultados de los algoritmos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wWiener = [wopt(1)*ones(N); wopt(2)*ones(N)];

figure('Name','Evolucion de la estimacion del filtro w','NumberTitle','off')
subplot(1,2,1);
plot(wWiener(1,:),'color','r');
hold
plot(wt(1,:),'color','g');
plot(wt2(1,:),'color','b');
plot(wt4(1,:),'color','k');
hold
wWiener = wopt(2)*ones(N);
subplot(1,2,2);
plot(wWiener(2,:),'color','r');
hold
plot(wt(2,:),'color','g');
plot(wt2(2,:),'color','b');
plot(wt4(2,:),'color','k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Gracias!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
