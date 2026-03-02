%Metodo de diferencas finitas em uma chapa 2D em regime permanente
clc; clear; close all;

%Dimensoes
Lx = input('Largura (m) da chapa:');
Ly = input('Comprimento(m):');

%Numero de nos, difere nos eixos para maior precisao a depender de Lx e Ly
%Lembrete para ajuste de malha (razao entre elas)
%O ideal seria uma malha em que dy seja proximo de dx (quadrada)
Nx = input('Numero de nós eixo X:');
Ny = input('Numero de nós eixo Y:');

%Arestas isoladas
I = input('Isolar alguma aresta?(n/inf/sup/esq/dir):','s');

%Temperatura de contorno
Ti = input('Temperatura (K) inferior:');
Ts = input('Temperatura (K) superior:');
Te = input('Temperatura (K) esquerda:');
Td = input('Temperatura (K) direita:');

%Espaco entre nos (m)
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
%Malha fina em uma direcao e grossa em outra pode fazer o calculo nao
%convergir

%Pre-processamento da matriz das temperaturas
T = zeros(Ny, Nx);

%Temperaturas de contorno fixadas
T(1,:) = Ti; %Inferior
T(Ny,:) = Ts; %Superior
T(:,1) = Te; %Esquerda
T(:,Nx) = Td; %Direita

%As pontas terão os valores médios entre elas
T(1,1) = (Ti + Te)/2;
T(1,Nx) = (Ti + Td)/2;
T(Ny,1) = (Ts + Te)/2;
T(Ny,Nx) = (Ts + Td)/2;

%Simplificando as variveis dx e dy na formula discretiada
A = (2/dx^2)+(2/dy^2);
B = (1/dx^2);
C = (1/dy^2);

%Iteracoes ate que se chegue no numero quase exato
erro = 100; 
while erro > 0.0001
t=T; %Para comparar o t antigo com T novo

%Duas formas, em que dx=dy (malha quadrada) e dx~=dy
if dx == dy
    for i = 2:Ny-1
        for j = 2:Nx-1
    T(i,j) = (T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1))/4;
    %Metodo de Liebmann (29.11 da apostila Chapra & Canale)
        end
    end
else 
    for i = 2:Ny-1
        for j = 2:Nx-1
    T(i,j) = B/A * (T(i+1,j) + T(i-1,j)) + C/A * (T(i,j+1) + T(i,j-1));
    %Equacao de Laplace discretizada
        end
    end
end

%Caso da borda isolada (fluxo de calor e gradiente termico = 0)
switch I
    case 'n' %sem aresta isolada

    case 'inf'
        for j = 2 : Nx-1 %Os valores do meio de j2 até o penultimo j
        T(1,j) = (T(1,j+1) + T(1,j-1) + 2*T(2,j))/4;
        %Formula 29.2 da apostila de Chapra & Canale
        %Quando ha aresta isolada e meio que o valor do no inverso *2
        end

    case 'sup'
        for j = 2 : Nx-1
        T(Ny,j) = (T(Ny,j+1) + T(Ny,j-1) + 2*T(Ny-1,j))/4;
        %O superior normalmente usa-se (1,j) pois a matriz comeca de cima
        %para baixo, porem ao plotar sera invertida para que Ny esteja no
        %topo, logo a superior tera o Ny fixada
        end

    case 'esq'
        for i = 2: Ny-1
        T(i,1) = (T(i+1,1) + T(i-1,1) + 2*T(i,2))/4;
        end

    case 'dir'
        for i = 2: Ny-1
        T(i,Nx) = (T(i+1,Nx) + T(i-1,Nx) + 2*T(i,Nx-1))/4;
        end
end


%Maior valor da matriz com os maiores valores de erro
erro = max(  max(   abs((T-t)./T)* (100 + eps)   )  );
%A formula acima sera executada e a iteracao so ira acabar se for <0.0001
end

figure;
contourf(T,20,'LineColor','none');
title('Mapa de Calor da Chapa')
yticks([]);xticks([]);
ylabel(sprintf('%g m', Ly));xlabel(sprintf('%g m', Lx));
colorbar;colormap hot;
axis equal;
exportgraphics(gcf, 'Mapa_Calor.png', 'BackgroundColor', 'none', 'Resolution', 300);


resp = input('Mostrar fluxo de calor? (s/n):','s');
if resp == 's'
    k = input('Coeficiente de condutividade térmica do material [W/m.K]:');
%Prata -> 420 W/m.K
%Aluminio -> 205 W/m.K
%Niquel -> 91 W/m.K
%Aço carbono -> 54 W/m.K
%Vidro -> 0.8 W/m.K

figure;
imagesc([0 Lx], [0 Ly], T);
title('Fluxo de Calor')
set(gca, 'YDir', 'normal');
colorbar; colormap hot; axis fill;
yticks([]);xticks([]);
ylabel(sprintf('%g m', Ly));xlabel(sprintf('%g m', Lx));
hold on;

exportgraphics(gcf, 'Mapa_Calor.png', 'BackgroundColor', 'none', 'Resolution', 300);

    %Lei de Fourier conducao termica (29.14 da apostila Chapra & Canale)
    [qx, qy] = gradient(T, dx, dy);
    qx = -k .* qx;
    qy = -k .* qy;
    
    %Sera criado um vetor apenas para identificar cada no da malha
    x = linspace(0, Lx, Nx); 
    y = linspace(0, Ly, Ny);

    %Usa-se o meshgrid para 'mapear' as coordenadas em cada ponto
    %De acordo exclusivamente de Lx Ly
    [X, Y] = meshgrid(x, y);
    %X e Y são matrizes ligadas da metragem da chapa 2D
    
    %Quiver X e Y sao as coordenadas e qx e qy os vetores para cada direcao
    %No livro de Chapra & Canale foi utilizado a formula em que 
    %Theta == arctan(qy/qx)
    %O comando quiver ja realiza o calculo da direcao
    
    quiver(X(Ny/10:Ny/5:end, Nx/10:Nx/5:end), ...
        Y(Ny/10:Ny/5:end, Nx/10:Nx/5:end), ...
        qx(Ny/10:Ny/5:end, Nx/10:Nx/5:end)*0.5, ...
        qy(Ny/10:Ny/5:end, Nx/10:Nx/5:end)*0.5, ...
        'AutoScale', 'off', 'Color', 'cyan', 'LineWidth', 1.2);

    exportgraphics(gcf, 'Fluxo_Calor_Niquel.png', 'BackgroundColor', 'none', 'Resolution', 300);
    %Os valores para se aplicar a seta foram definidos a partir da largura
    %da chapa, as primeiras em Ly/10 e seguindo Ly/5 sucessivamente
    %feito da mesma forma no eixo x
    %Escala das setas refletem os valores de qx e qy (fluxo de calor)
end