%Ler os ficheiros de audio
wave = cell(50);%contem todas as ondas sonoras
fs = cell(50);
%caracteristicas = cell(50,10);
figure(1)
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
for j = 0:9
    for i =1:50
        t = 1:size(wave{i,j+1});
        subplot(2,5,j+1);
        plot(t,wave{i,j+1})
    end
    title(int2str(j))
    xlabel("Tempo")
    ylabel("Amplitude")
end
%% fazer a analise com a tranformada de fourier
tamJanela = 48000;%testar para várias janelas
freq = cell(10,1);%ver os picos em diferentes itervalos de freq
Q1 = cell(10,1);
Q3 = cell(10,1);
figure(2)
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


