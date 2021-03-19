%% TRABALHO PR�TICO - PROCESSAMENTO DE SINAIS
% ELE 042
% V�tor Gabriel Reis Caitit� - 2016111849 \\
% Willian Braga da Silva - 201602762


%% QUEST�O 1:
% Abrir o arquivo Blind_intro.wav e gerar um gr�fico da 
% forma de onda em fun��o do tempo, similar ao mostrado na figura
% do tp.
%_____________________________________________________________________

clear all;
clc;
audioArquive = 'Blind_intro.wav';
%    [Y, FS]=audioread(FILENAME) reads an audio file specified by the 
%    string FILE, returning the sampled data in Y and the sample rate 
%    FS, in Hertz. 
[Y, FS]=audioread(audioArquive);

%Caso deseje ouvir o audio basta descomentar a linha abaixo
%sound(Y, FS)

%Gera��o da escala de tempo a partir de FS:
period = 1/FS;                              % Per�odo
L = length(Y);                              % Largura de Faixa
timeScale = linspace(0, L-1, L)*period;     % Escala de Frequ�ncia
%Plot:
figure(1)
plot(timeScale, Y);                         % Gr�fico do sinal no tempo
xlabel('Tempo (segundos)');                 % Eixo x
ylabel('Amplitude');                        % Eixo y   
title('Sinal no Dom�nio do Tempo');         % T�tulo

% Descri��o da Atividade: Nesse exerc�cio o primeiro passo foi ler o 
% arquivo de �udio com a fun��o audioread, a qual retorna os dados  
% amostrados e tamb�m a frequ�ncia de amostragem.

%_____________________________________________________________________
%% QUEST�O 2:
% Calcular e mostrar a composi��o espectral do sinal utilizando a 
% transformada r�pida de Fourier (FFT). Gere os gr�ficos de amplitude e
% fase do espectro de frequ�ncias.

% Transformada R�pida de Fourier Y:
fftY = fft(Y); 
% M�dulo da fft:
moduloY = abs(fftY);
% Fase da fft:
faseY = unwrap(angle(fftY)); 
%Gera��o de escala de frequ�ncia positiva:
frequencyScale = (linspace(0, L/2, L/2+1))*FS/L;
%Gera��o de escala de frequ�ncia completa:
aux = length(frequencyScale)+1;
for i=length(frequencyScale):-1:1
aux = aux-1;
completeScale(i) = -frequencyScale(aux);
end
aux=0;
for i=(2*length(frequencyScale)-1):-1:length(frequencyScale)
aux=aux+1;
completeScale(i) = frequencyScale(aux);
end
f = linspace(-FS/2,FS/2, L);
%Plots:
completeScale=completeScale';
figure(2)
plot(completeScale, moduloY);                           % Plot m�dulo
ylabel('|X(e^{jw})|');                                  % Eixo y
xlabel('Freq (kHz)');                                   % Eixo x
title('M�dulo de Y(e^{jw})');                          % T�tulo
figure(3)
plot(f, faseY);                                         % Plot fase
ylabel('$\theta(\omega)$', 'interpreter', 'latex');     % Eixo y 
xlabel('Frequ�ncia (kHz)');                             % Eixo x
title('Fase de $Y(e^{jw})$');    % T�tulo

% Descri��o da Atividade: Primeiramente gerou-se a transformada r�pida de
% Fourier usando a fun��o fft(), dpois disso p�de-se obter o modulo e fase
% com as fun��es abs() e angle() respectivamente. Ent�o bastou gerar as
% escalas de frequeu�ncia e plotar.

%_____________________________________________________________________
%% QUEST�O 3: 
% Fa�a o projeto do banco de filtros.
%   . Subgraves: entre 16 Hz to 60 Hz;
%   . Graves: entre 60 Hz to 250 Hz;
%   . M�dio-graves: entre 250 Hz to 2 kHz;
%   . M�dio-agudos: entre 2 kHz to 4 kHz;
%   . Agudos: entre 4 kHz to 6 kHz;
%   . Brilho: entre 6 kHz to 16 kHz;

% Vetor de Frequ�ncias Inferiores:
freqInferior = [16 60 250 2000 4000 6000];  
% Vetor de Frequ�ncias Superiores:
freqSuperior = [60 250 2000 4000 6000 16000];       
% Vetor de t�tulos dos filtros:
filterType = {'Filtro Subgrave'; 'Filtro Grave'; 'Filtro Medio Grave'
              'Filtro Medio Agudo'; 'Filtro Agudo'; 'Filtro Brilho'}; 

for count = 1:6
    % Definindo par�metros para a fun��o cheb1ord:
    % [N, Wp] = cheb1ord(Wp, Ws, Rp, Rs) returns the order N of the lowest
    % order digital Chebyshev Type I filter which has a passband ripple of 
    % no more than Rp dB and a stopband attenuation of at least Rs dB. Wp 
    % and Ws are the passband and stopband edge frequencies, normalized 
    % from 0 to 1 (where 1 corresponds to pi radians/sample). For example,
    %  .  Lowpass:    Wp = .1,      Ws = .2
    %  .  Highpass:   Wp = .2,      Ws = .1
    %  .  Bandpass:   Wp = [.2 .7], Ws = [.1 .8]
    %  .  Bandstop:   Wp = [.1 .8], Ws = [.2 .7]
    
    fpI = freqInferior(count)*2/FS;
    fpS = freqSuperior(count)*2/FS;
    RsInferior = fpI - 0.08;
    RsSuperior = fpS + 0.08;
    if RsInferior < 0
        RsInferior = 0.0000001;
    end
    if RsSuperior > 1
        RsSuperior = 0.999999;
    end
    
    % Gerando Ordem do Filtro Chebyshev tipo I:
    [order ~] = cheb1ord([fpI fpS], [RsInferior RsSuperior],  0.5, 40); 
    
    % Designs a bandpass filter: cheby1(N,R,Wp,'bandpass')
    % Cheby1 returns the filter coefficients in length N+1 vectors B 
    %(numerator) and A (denominator).
    [B A] = cheby1(order, 0.5, [fpI fpS], 'bandpass');
    
    %[NUMd,DENd] = bilinear(B,A,FS)
    % Gerandoo sinal passado pelo filtro: 
    y{count} = filter(B,A, Y);          
    %h = fvtool(B,A);      
      
    figure(count+3)
    % Gr�fico de M�dulo e Fase do Filtro:
    freqz(B, A); 
    title(filterType{count});
    
    figure(10) 
    subplot(2, 3, count)
    plot(timeScale, y{count}) % Sinal ap�s o filtro
    xlabel('Tempo(s)')
    ylabel('Amplitude')
    title(filterType{count})
    
    figure(11)
    % Diagrama de polos e zeros do filtro:
    subplot(3, 2, count)
    %Fun��o de Transfer�ncia: 
    H(count) = tf(B, A, period)
    pzmap(H(count))
    title(filterType{count});
    
end

% Descri��o da Atividade: Primeiramente gerou-se os filtros de Chebyshev 
% tipo I, utilizando as fun��es cheb1ord(), cheby1(), cada uma ja
% explicada. Logo bastou passar o sinal pelos filtros especificados,
% utilizando a fun��o filter(). Assim pod�-se gerar todos os gr�ficos, a
% fun��o de tranfer�ncia (usando a fun��o tf()), e o diagrama de polos e
% zeros (usando a fun��o pzmap()).

%_____________________________________________________________________
%% QUEST�O 4:
% Ap�s o ajuste dos ganhos, o sistema dever� realizar a filtragem e 
% reconstru��o do sinal por meio da soma das sa�das dos filtros 
% individuais.

finalY = y{1}+y{2}+y{3}+y{4}+y{5}+y{6};
figure(12)
plot(timeScale, finalY);                    % Gr�fico do sinal no tempo
xlabel('Tempo (segundos)');                 % Eixo x
ylabel('Amplitude');                        % Eixo y   
title('Sinal no Dom�nio do Tempo Filtrado');% T�tulo

% Descri��o da Atividade: Bastou somar todos os sinais, cada um passado por
% um tipo de filtro e plotar.
%% QUEST�O 5:
%Reproduza o sinal filtrado no sistema de �udio do computador.
sound(finalY, FS);

% Descri��o da Atividade: Simplemente utilizar a fun��o sound para
% reproduzir o audio.
%%

