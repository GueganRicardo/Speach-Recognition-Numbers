eixo_x = 0:9;
n_figura = 1;
%Ler os ficheiros de audio
wave = cell(50, 10);%contem todas as ondas sonoras
fs = cell(50, 10);
%caracteristicas = cell(50,10);
for j = 0:9
    for i = 0:49
        local_fich = ["34\"];
        local_fich = local_fich.append(int2str(j));
        local_fich = local_fich.append("_34_");
        local_fich = local_fich.append(int2str(i));
        local_fich = local_fich.append(".wav");
        [wave{i+1,j+1} ,fs{i+1,j+1}] =audioread(local_fich);
    end
end
%% fazer a analise no tempo
caracteristica = zeros(50,3);
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

caracteristica = zeros(50,3);
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
%trim_wave = cell(50, 10);
janela_sz = 100;
limiar = 0.05;

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


        energia_janela = 0;
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
        trim_wave = norm_wave(lower_lim:size(norm_wave));
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

%% plotar a energia em boxplots mas em janelas

n_janelas = 50; %tamanho do vetor / numero de janelas = numero de samples em cada janela. assim todas as janelas têm o mesmo tamanho relativo. possivelmente sera necessario dar trim do vetor para detetar quando começa e acaba o som
new_wave = wave;
energias_metades = zeros(50, 10, n_janelas);
for j = 1:n_janelas
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
            disp(energia(s, e + 1))
            energias_metades(s, e + 1, j) = energia(s, e + 1);
        end 
    end
    boxplot(energia, eixo_x);
    ylim([0 200]);
end

%boxplot do racio
figure(n_figura);
n_figura = n_figura + 1;
racios = energias_metades(1:50, 1:10, 1) ./ energias_metades(1:50, 1:10, 2);
boxplot(racios, eixo_x);
ylim([0 30]);


%% fazer a analise com a tranformada de fourier
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
    end
    freq{j} = median(sinais, 2);
    Q1{j} = quantile(sinais, 0.25, 2);
    Q3{j} = quantile(sinais, 0.75, 2);
end

for d = 1:10
    subplot(2,5,d);
    plot(freq{d});
    hold on;
    plot(Q1{d},'r--');
    plot(Q3{d},'k--');
    hold off;
    legend('Median', 'Q1', 'Q3');
    title(num2str(d-1));
    xlabel('Frequência (Hz)');
    ylabel('Amplitude (normalizada)');
end
%% fazer a analise com a STFT
%efetito das diferentes janelas sobreposição numerod e pontos


