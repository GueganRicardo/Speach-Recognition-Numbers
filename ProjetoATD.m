eixo_x = 0:9;
n_figura = 1;
%Ler os ficheiros de audio
wave = cell(50, 10);%contem todas as ondas sonoras
fs = cell(50, 10);
possibilidades = zeros(10,1);
for j = 0:9
    for i = 0:49
        local_fich = ("34\");
        local_fich = local_fich.append(int2str(j));
        local_fich = local_fich.append("_34_");
        local_fich = local_fich.append(int2str(i));
        local_fich = local_fich.append(".wav");
        [wave{i+1,j+1} ,fs{i+1,j+1}] =audioread(local_fich);
    end
end
%% filtrar frequencias de voz do ruido (1000 - 4000 hz)
filtered_waves = cell(50, 10);

tamJanela = 48000;%testar para várias janelas
freq = cell(10,1);%ver os picos em diferentes itervalos de freq
Q1 = cell(10,1);
Q3 = cell(10,1);
figure(n_figura);
n_figura = n_figura + 1;
for j = 1:10
    sinais = zeros(tamJanela/2+1, 50);
    for i = 1:50
        X = fft(wave{i,j}, tamJanela);
        X_positivo = X(1:tamJanela/2+1);
        X_normalizado = abs(X_positivo) / tamJanela;
        sinais(:,i) = X_normalizado;
        sinais (1:1000, i) = 0; % filtro
        filtered_waves{i, j} = ifft(X);
    end
    freq{j} = median(sinais, 2);
    Q1{j} = quantile(sinais, 0.25, 2);
    Q3{j} = quantile(sinais, 0.75, 2);
end
for j = 1:10
    subplot(2,5,j);
    plot(wave{i, j});
    hold on;
    title(num2str(j-1));
    xlabel('Frequência (Hz)');
    ylabel('Amplitude (normalizada)');
end
%% fazer a analise no tempo
figure(n_figura);
n_figura = n_figura + 1;
for j = 0:9
    subplot(2,5,j+1);
    for i =1:50
        t = 1:size(wave{i,j+1});
        plot(t,wave{i,j+1})
        hold on
    end
    title(int2str(j))
    xlabel("Tempo")
    ylabel("Amplitude")
end


%% analise de tempo / amplitude mas norm em x e y
figure(n_figura);
n_figura = n_figura + 1;
[k, l] = cellfun(@size, wave);
target_sz_matrix = max(k);
for j = 0:9
    subplot(2,5,j+1);
    for i =1:50
        %disp("dentro do ciclo")
        %this_wave = interp1(1:size(wave{i,j+1}), wave{i,j+1}, 1:target_sz_matrix(j + 1));
        this_wave = wave{i, j + 1};
        max_amp = max(abs(this_wave));
        norm_wave = this_wave / max_amp;
        t = linspace(0, 1, size(norm_wave, 1));
        plot(t,norm_wave);
        hold on
    end
    title(int2str(j))
    xlabel("Tempo")
    ylabel("Amplitude")
end


%% dar trim para alinhar, com base na energia de janelas constantes
trim_waves = cell(50, 10);
janela_sz = 100;
limiar = 0.025;

caracteristica = zeros(50,3);
figure(n_figura);
n_figura = n_figura + 1;
for j = 0:9
    subplot(2,5,j+1);
    for i =1:50
        %disp("dentro do ciclo")
        %this_wave = interp1(1:size(wave{i,j+1}), wave{i,j+1}, 1:target_sz_matrix(j + 1));
        this_wave = wave{i, j + 1};
        max_amp = max(abs(this_wave));
        norm_wave = this_wave / max_amp;

        %lower lim
        janela = -1;
        tolerancia = 0;
        while tolerancia < 10
            janela = janela + 1;
            aux = power(norm_wave((janela_sz * janela) + 1:(janela + 1) * janela_sz, 1), 2);
            energia_janela = sum(aux);
            if energia_janela > limiar
                tolerancia = tolerancia + 1;
                if tolerancia == 1
                    lower_lim = (janela * janela_sz) + 1;
                end
            else
                tolerancia = 0;
            end
        end
        
        %upper_lim
        energia_janela = 0;
        janela = -1;
        tolerancia = 0;
        while tolerancia < 5
            janela = janela + 1; %as janelas agora avançam para trás no tempo
            aux = power(norm_wave(size(this_wave) - ((janela + 1) * janela_sz) + 1 : size(this_wave) - (janela_sz * janela) , 1), 2);
            energia_janela = sum(aux);
            if energia_janela > limiar
                tolerancia = tolerancia + 1;
                if tolerancia == 1
                    upper_lim = size(this_wave) - (janela_sz * janela);
                end
            else
                tolerancia = 0;
            end
        end



        trim_wave = norm_wave(lower_lim:upper_lim);
        trim_waves{i, j + 1} = trim_wave;
        t = linspace(0, 1, size(trim_wave, 1));
        plot(t,trim_wave);
        hold on
    end
    title(int2str(j))
    xlabel("Tempo")
    ylabel("Amplitude")
end

%% plotar a enenrgia em boxplots
energia = zeros(50, 10);
figure(n_figura);
n_figura = n_figura + 1;
% cada digito
for e = 0:9
    % cada sinal
    for s = 1:50
        aux = power(wave{s, e + 1}, 2);
        energia(s, e + 1) = sum(aux);
    end
end
boxplot(energia, eixo_x);

%% plotar a energia em linhas mas em janelas
n_janelas = 50; %tamanho do vetor / numero de janelas = numero de samples em cada janela. assim todas as janelas têm o mesmo tamanho relativo. possivelmente sera necessario dar trim do vetor para detetar quando começa e acaba o som
new_wave = trim_waves;
energias_janelas = zeros(50, 10, n_janelas);

for j = 1:n_janelas
    % cada digito
    for e = 0:9
        % cada sinal
        for s = 1:50
            this_wave = new_wave{s, e + 1};
            max_amp = max(this_wave);
            %disp(max_amp);
            this_wave = this_wave / max_amp;
            window_sz = fix(size(this_wave, 1) / n_janelas);
            lower_lim = window_sz * (j - 1) + 1;
            upper_lim = window_sz * j;
            aux = power(this_wave(lower_lim:upper_lim, 1), 2);
            energia(s, e + 1) = sum(aux);
            %disp(energia(s, e + 1))
            energias_janelas(s, e + 1, j) = energia(s, e + 1);
        end 
    end
    %plot(energia, eixo_x);
    %ylim([0 200]);
end

figure(n_figura);
n_figura = n_figura + 1;
for e = 0:9
    subplot(2,5,e+1);
    for s = 1:50
        plot(1:n_janelas, squeeze(energias_janelas(s, e + 1, :)))
        hold on
    end
    title(int2str(e))
    xlabel("Tempo")
    ylabel("Energia")
    ylim([0 200])
end

% plotar as energias normalizadas (e_parcial_1 + e_parcial_2 + ... +
% e_parcial_n = e_total
figure(n_figura);
n_figura = n_figura + 1;
for e = 0:9
    subplot(2,5,e+1);
    for s = 1:50
        plot(1:n_janelas, squeeze(energias_janelas(s, e + 1, :)))
        hold on
    end
    title(int2str(e))
    xlabel("Tempo")
    ylabel("Energia Normalizada")
    ylim([0 200])
end

%% boxplot do racio
% figure(n_figura);
% n_figura = n_figura + 1;
% racios = energias_metades(1:50, 1:10, 1) ./ energias_metades(1:50, 1:10, 2);
% boxplot(racios, eixo_x);
% ylim([0 30]);
%% plotar a energia em boxplots mas em janelas
n_janelas = 10; %tamanho do vetor / numero de janelas = numero de samples em cada janela. assim todas as janelas têm o mesmo tamanho relativo. possivelmente sera necessario dar trim do vetor para detetar quando começa e acaba o som
new_wave = trim_waves;
energias_metades = zeros(50, 10, n_janelas);
for j = 1:n_janelas
    hit = zeros(10,1);
    figure(n_figura);
    n_figura = n_figura + 1;
    % cada digito
    for e = 0:9
        % cada sinal
        for s = 1:50
            this_wave = new_wave{s, e + 1};
            max_amp = max(this_wave);
            %disp(max_amp);
            this_wave = this_wave / max_amp;
            window_sz = fix(size(this_wave, 1) / n_janelas);
            lower_lim = window_sz * (j - 1) + 1;
            upper_lim = window_sz * j;
            aux = power(this_wave(lower_lim:upper_lim, 1), 2);
            energia(s, e + 1) = sum(aux);
            %disp(energia(s, e + 1))
            energias_metades(s, e + 1, j) = energia(s, e + 1);
            %fazer os if's das janelas
            % 90% rate em 0 3 100% 2 7 98% 6
            if j == 2 && energia(s,e+1)>100
                hit(e+1)= hit(e+1)+1;
                possibilidades(0 + 1) =-1;
                possibilidades(3 + 1) =-1;
                possibilidades(2 + 1) =-1;
                possibilidades(7 + 1) =-1;
                possibilidades(6 + 1) =-1;
            end
%             94/98% rate em 6 8
%             if j == 6 && energia(s,e+1)>110
%                 hit(e+1)= hit(e+1)+1;
%                 possibilidades(6 + 1) =-1;
%                 possibilidades(8 + 1) =-1;
%             end
            % 100% rate em 6 8
            if j == 7 && energia(s,e+1)>85
                hit(e+1)= hit(e+1)+1;
                possibilidades(6 + 1) =-1;
                possibilidades(8 + 1) =-1;
            end
        end
    end
    boxplot(energia, eixo_x);
    ylim([0 1000]);
    disp(hit);
end

%boxplot do racio
figure(n_figura);
n_figura = n_figura + 1;
racios = energias_metades(1:50, 1:10, 1) ./ energias_metades(1:50, 1:10, 2);
boxplot(racios, eixo_x);
ylim([0 30]);

%% decisão com base nos piccos negativos 
%permite concluir que não pode ser um 2 ou 6 ou 7
%com base nas amostras tem 2% de falha
treshHold = -0.80;
hit = zeros(10,1);
for sinal = 1:50
    for digito = 0:9
        count = sum(trim_waves{sinal,digito+1}<treshHold);
        if count<10
            hit(digito+1,1) = hit(digito+1,1) + 1;
            possibilidades(2+1) = -1;
            possibilidades(6+1)=-1;
            possibilidades(7+1)=-1;
        end
    end
end
disp(hit)
disp(possibilidades)
%% fazer a analise com a tranformada de fourier
tamJanela = 48000%testar para várias janelas
freq = cell(10,1);%ver os picos em diferentes itervalos de freq
Q1 = cell(10,1);
Q3 = cell(10,1);
figure(n_figura);
n_figura = n_figura + 1;
sinais_todos=zeros(10,tamJanela/2+1, 50);
for j = 1:10
    sinais = zeros(tamJanela/2+1, 50);
    for i = 1:50
        % Apply a highpass filter to the signal
        trimmed_signal = trim_waves{i, j};
        filtered_signal = highpass(trimmed_signal, 1000, 48000);
        
        % Apply the Hamming window to the filtered signal
        windowed_signal = filtered_signal;% .* blackman(length(filtered_signal));
        
        X = fft(windowed_signal, tamJanela);
        X_positivo = X(1:tamJanela/2+1);
        X_normalizado = abs(X_positivo) / tamJanela;
        sinais(:,i) = X_normalizado;
    end
    freq{j} = median(sinais, 2);
    Q1{j} = quantile(sinais, 0.25, 2);
    Q3{j} = quantile(sinais, 0.75, 2);
    sinais_todos(j,:,:)=sinais(:,:);
end

for d = 1:10
    subplot(2,5,d);
    plot(freq{d});
    hold on;
    plot(Q1{d},'r--');
    plot(Q3{d},'k--');
    ylim([0,0.005])
    hold off;
    legend('Median', 'Q1', 'Q3');
    title(num2str(d-1));
    xlabel('Frequência (Hz)');
    ylabel('Amplitude (normalizada)');
end
%% boxplot do spectrall roll off
figure(n_figura);
n_figura = n_figura + 1;
rolloff = zeros(50,10);
for d = 0:9
    
    for sinal = 1:50
        this_wave = sinais_todos(d+1,:,sinal);
        this_wave = this_wave.^2;
        this_wave = this_wave/(sum(this_wave));
        somatorio = cumsum(this_wave);
        rolloff(sinal,d+1) = find(somatorio >= 0.85,1);
    end
    
end
boxplot(rolloff,0:9)
%% encontrar os picos do sinal
figure(n_figura);
n_figura = n_figura + 1;
picos = zeros(50,10,3);
for d = 0:9
    for sinal = 1:50
        this_wave = sinais_todos(d+1,:,sinal);
        findpeaks(this_wave,'MinPeakHeight',0.001)
    end
end
%% fazer a analise com a STFT
%efetito das diferentes janelas sobreposição numerod e pontos
figure(n_figura);
n_figura = n_figura + 1;
% Define the signal and its parameters
fs = 48000;  % Sample rate (Hz)
t = 0:1/fs:1;  % Time vector (1 second)
f1 = 1000;  % Frequency of the signal (Hz)
% Define the STFT parameters
frameSize = 1024;
overlap = 0.75;  % 75% overlap
nfft = 1024;  % FFT size

for e = 0:9
    subplot(4, 3, e + 1)
    sinais = zeros(tamJanela/2+1, 50);
    for i = 1:50
        % Apply a highpass filter to the signal
        trimmed_signal = trim_waves{i, e + 1};
        filtered_signal = highpass(trimmed_signal, 1000, 48000);
        sinais(1:size(filtered_signal),i) = filtered_signal;
    end

    
    
    % Calculate the spectrogram
    spectrogram(median(sinais, 2), frameSize, round(overlap*frameSize), nfft, fs, 'yaxis');
    
    % Set the color map for the spectrogram
    colormap jet;
    
    % Add labels and title to the plot
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(num2str(e));
    end



