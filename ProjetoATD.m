%Ler os ficheiros de audio
wave = cell(50);
fs = cell(50);
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
    for i =1:10
        t = 1:size(wave{i,j+1});
        subplot(4,3,j+1);
        plot(t,wave{i,j+1})
    end
    title(int2str(j))
    xlabel("Tempo")
    ylabel("Amplitude")
end
%%