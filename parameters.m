%% 2) Separar em intervalos de 5 segundos %%

% Inicializar a célula que armazenará os intervalos
intervals = {};

% Definir a duração do bloco de 5 segundos em termos de amostras
fs = EEG.srate; % Taxa de amostragem
block_duration = 5; % Duração de cada bloco em segundos
samples_per_block = block_duration * fs;

% Loop pelos eventos, exceto o último
for i = 1:length(EEG.event)-1
    % Pegar o tipo de evento atual e próximo
    current_event = EEG.event(i).type;
    next_event = EEG.event(i+1).type;
    
    % Pegar a latência dos eventos
    start_latency = EEG.event(i).latency;
    end_latency = EEG.event(i+1).latency;
    
    % Verificação para garantir que as latências são diferentes
    if start_latency == end_latency
        disp(['Eventos com latência igual encontrados: ' current_event ' e ' next_event]);
        continue; % Ignorar este intervalo
    end
    
    % Processar o intervalo se as latências forem diferentes
    interval_data = EEG.data(:, round(start_latency):round(end_latency));
    interval_duration = (end_latency - start_latency) / fs; % Duração em segundos
    
    % Verificar se o intervalo tem pelo menos 5 segundos
    if interval_duration < block_duration
        disp(['Intervalo entre ' current_event ' e ' next_event ' é menor que 5 segundos.']);
        continue; % Ignorar este intervalo
    end
    
    % Dividir o intervalo em blocos de 5 segundos
    num_samples = size(interval_data, 2); % Número total de amostras no intervalo
    num_blocks = floor(num_samples / samples_per_block); % Número de blocos completos de 5s
    
    blocks = {}; % Inicializa célula para armazenar os blocos de 5 segundos
    for b = 1:num_blocks
        start_idx = (b-1)*samples_per_block + 1;
        end_idx = b*samples_per_block;
        blocks{b} = interval_data(:, start_idx:end_idx);
        
        % Plotar o bloco
        %figure;
        %plot((1:samples_per_block)/fs, blocks{b});
        %title(['Bloco ' num2str(b) ' - Intervalo: ' current_event ' até ' next_event]);
        %xlabel('Tempo (s)');
        %ylabel('Amplitude');
    end
    
    % Caso tenha amostras restantes após o último bloco completo de 5 segundos
    if num_samples > num_blocks * samples_per_block
        remaining_samples = num_samples - num_blocks * samples_per_block;
        blocks{end+1} = interval_data(:, end-remaining_samples+1:end);
        
        % Plotar o bloco restante
        %figure;
        %plot((1:remaining_samples)/fs, blocks{end});
        %title(['Bloco restante - Intervalo: ' current_event ' até ' next_event]);
        %xlabel('Tempo (s)');
        %ylabel('Amplitude');
    end
    
    % Armazenar o intervalo e seus blocos
    intervals{end+1} = struct(...
        'start_event', current_event, ...
        'end_event', next_event, ...
        'start_latency', start_latency, ...
        'end_latency', end_latency, ...
        'duration', interval_duration, ...
        'EEG_data', interval_data, ...
        'blocks', {blocks} ...
    );
end
disp('Processamento da parte 2 completo!');

%% 3) PSD - Band Power Analysis %%

% Número de amostras e frequência de amostragem (definir conforme o seu sinal EEG)
fs = 500;  % Exemplo de taxa de amostragem (250 Hz), ajustar conforme necessário

% Definir as bandas de frequência
freq_bands = [1 4; 4 8; 8 13; 13 30; 30 50];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};

% Inicializando a variável para armazenar a potência por banda
band_power = cell(1, length(intervals));

for i = 1:length(intervals)
    % Extrair os dados do EEG para o intervalo atual
    data_segment = intervals{i}.EEG_data;

    % Aplicar a FFT para calcular a densidade espectral de potência (PSD)
    [psd, freq] = pwelch(data_segment', [], [], [], EEG.srate);

    % Inicializar a célula para armazenar a potência de cada banda
    band_power{i} = struct('delta', [], 'theta', [], 'alpha', [], 'beta', [], 'gamma', []);

    % Definir os intervalos de frequência para cada banda
    freq_bands = [1 4; 4 8; 8 13; 13 30; 30 50];
    band_names = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

    % Para cada banda de frequência
    for k = 1:size(freq_bands, 1)
        % Encontrar os índices das frequências que estão dentro do intervalo da banda
        freq_idx = find(freq >= freq_bands(k, 1) & freq <= freq_bands(k, 2));

        % Calcular a potência da banda integrando o PSD
        band_power{i}.(band_names{k}) = trapz(freq(freq_idx), psd(freq_idx, :));
    end
end

% Para exibir os valores da potência das bandas de cada intervalo
% for i = 1:length(band_power)
%     disp(['Intervalo ' num2str(i)]);
%     disp(['Delta power: ' num2str(band_power{i}.delta)]);
%     disp(['Theta power: ' num2str(band_power{i}.theta)]);
%     disp(['Alpha power: ' num2str(band_power{i}.alpha)]);
%     disp(['Beta power: ' num2str(band_power{i}.beta)]);
%     disp(['Gamma power: ' num2str(band_power{i}.gamma)]);
%     disp('------------------------------');
% end

% Para cada intervalo principal
for i = 1:length(intervals)
    % Recuperar a duração do intervalo
    duration = intervals{i}.duration;
    
    % Calcular o número de blocos de 5 segundos no intervalo
    num_blocks = floor(duration / block_duration); 
    
    % Inicializar uma célula para armazenar as potências de banda de cada bloco de 5 segundos
    intervals{i}.band_power = cell(1, num_blocks);
    
    % Para cada bloco de 5 segundos dentro do intervalo
    for k = 1:num_blocks
        % Obter os dados EEG do bloco atual
        start_idx = (k - 1) * samples_per_block + 1;
        end_idx = k * samples_per_block;
        EEG_block = intervals{i}.EEG_data(:, start_idx:end_idx);
        
        % Calcular a PSD usando pwelch
        [psd, freq] = pwelch(EEG_block', [], [], [], fs);
        
        % Inicializar uma estrutura para armazenar as potências de banda
        band_power_block = struct();
        
        % Para cada banda de frequência, calcular a potência
        delta_idx = (freq >= 1 & freq <= 4);
        theta_idx = (freq >= 4 & freq <= 8);
        alpha_idx = (freq >= 8 & freq <= 13);
        beta_idx = (freq >= 13 & freq <= 30);
        gamma_idx = (freq >= 30 & freq <= 50);
        
        % Calcular a potência em cada banda integrando a PSD
        band_power_block.delta = trapz(freq(delta_idx), psd(delta_idx, :));
        band_power_block.theta = trapz(freq(theta_idx), psd(theta_idx, :));
        band_power_block.alpha = trapz(freq(alpha_idx), psd(alpha_idx, :));
        band_power_block.beta = trapz(freq(beta_idx), psd(beta_idx, :));
        band_power_block.gamma = trapz(freq(gamma_idx), psd(gamma_idx, :));
        band_power_block.total = band_power_block.delta + band_power_block.theta + band_power_block.alpha + band_power_block.beta + band_power_block.gamma;

        % Armazenar a potência de banda no campo 'band_power' do intervalo atual
        intervals{i}.band_power{k} = band_power_block;
    end
end

% Para plotar o PSD dos intervalos completos
for i = 1:length(intervals)
    %figure;  % Criar uma nova figura para cada intervalo
    [psd, freq] = pwelch(intervals{i}.EEG_data', [], [], [], fs);  % Calcular a PSD novamente, se necessário

    %Plotar o gráfico da PSD para o intervalo atual
    % plot(freq, 10*log10(psd)); 
    % title(['Intervalo ' num2str(i) ' - Power Spectral Density']);
    % xlabel('Frequência (Hz)');
    % ylabel('Densidade Espectral de Potência (dB/Hz)');
    % xlim([0 60]);  % Limitar o eixo x até 60 Hz
    % grid on;

    % Adicionar uma legenda indicando o intervalo
    % legend(band_names);
end

% Para plotar a PSD dos blocos de 5 segundos
for i = 1:length(intervals)
    num_blocks = length(intervals{i}.EEG_data) / (block_duration * fs);
    
    for k = 1:num_blocks
        % Obter os dados do bloco atual
        start_idx = (k - 1) * samples_per_block + 1;
        end_idx = k * samples_per_block;
        EEG_block = intervals{i}.EEG_data(:, start_idx:end_idx);
        
        % Calcular a PSD do bloco de 5 segundos
        [psd, freq] = pwelch(EEG_block', [], [], [], fs);

        % Criar um novo gráfico para cada bloco
        % figure;  
        % plot(freq, psd);
        % title(['Intervalo ' num2str(i) ' - Bloco ' num2str(k) ' de 5 segundos']);
        % xlabel('Frequência (Hz)');
        % ylabel('Densidade Espectral de Potência (dB/Hz)');
        % xlim([0 60]); %  Isso limita o eixo X para mostrar apenas até 60 Hz. Você pode ajustar para mostrar um intervalo mais apropriado.
        % grid on;
        % 
        % % Adicionar uma legenda com o nome das bandas
        % legend(band_names);
    end
end
disp('Processamento da parte 3 completo!');

%% 4) PSD - Relative Power Band  %%

fs = 500;  % Frequência de amostragem do EEG

% Definir as bandas de frequência
freq_bands = [1 4; 4 8; 8 13; 13 30; 30 50];
band_names = {'delta', 'theta', 'alpha', 'beta', 'gamma'};

for i = 1:length(intervals) %indo de 1 a 20, percorrendo todos os intervalos maiores.
    num_blocks = length(intervals{i}.band_power); 
    % Número de blocos de 5 segundos já existentes, no caso seu valor sera a quantidade de blocos de 5 segundos existentes em cada bloco maior.

    % Inicializando a célula para armazenar a RBP de cada bloco
    intervals{i}.rbp = cell(1, num_blocks);

    for k = 1:num_blocks % Valor que vai de 1 a quantidade de intervalos de 5 segundos que existe em cada bloco
        % Extrair o bloco de 5 segundos atual
        data_block = intervals{i}.band_power{k}; % percorre cada band power de todos os intervalos de 5 segundos.

        % Inicializar a estrutura para armazenar a RBP para o bloco atual
        rbp_block = struct();

        % Cálculo da Relative Band Power (RBP) para cada banda
        for b = 1:size(freq_bands, 1) % Percorre cada uma das frequency bands (alpha,beta,gamma,delta,theta)
            for a = 1:size(EEG_block, 1) % Era 63, porém pode alternar se removermos certos canais, por isso tem q ter um tamanho dinâmico
                rbp_block.(band_names{b})(a) = data_block.(band_names{b})(1,a) / data_block.total(a);
            end
        end

        % Armazenar o RBP calculado no bloco atual
        intervals{i}.rbp{k} = rbp_block;
    end
end

disp('Processamento da parte 4 completo!');

%% 5) Peak frequency and Power  %%

fs = 500;  % Frequência de amostragem, ajuste conforme necessário
alpha_band = [8, 13];  % Definindo a faixa de frequência Alpha

for i = 1:length(intervals)
    % Inicializar células para armazenar frequência e potência de pico
    intervals{i}.alpha_peak_frequency = cell(1, length(intervals{i}.band_power));
    intervals{i}.alpha_peak_power = cell(1, length(intervals{i}.band_power));

    % Loop pelos blocos de 5 segundos
    for k = 1:length(intervals{i}.band_power)
        EEG_block = intervals{i}.blocks{k};

        % Inicializar variáveis de pico para cada canal
        alpha_peak_frequency = zeros(1, size(EEG_block, 1));
        alpha_peak_power = zeros(1, size(EEG_block, 1));

        % Loop pelos canais
        for ch = 1:size(EEG_block, 1)
            % Calcular a PSD para o canal e bloco atuais
            [psd, freq] = pwelch(EEG_block(ch, :), [], [], [], fs);

            % Encontrar índices dentro da banda Alpha
            alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));

            % Encontrar a frequência e a potência máxima dentro da banda Alpha
            [alpha_peak_power(ch), max_idx] = max(psd(alpha_idx));
            alpha_peak_frequency(ch) = freq(alpha_idx(max_idx));
        end

        % Armazenar nos campos de intervals
        intervals{i}.alpha_peak_frequency{k} = alpha_peak_frequency;
        intervals{i}.alpha_peak_power{k} = alpha_peak_power;
    end
end

disp('Processamento da parte 5 completo!');

%% 6) Power Ratio %%

fs = 500;  % Frequência de amostragem do EEG

% Definir as bandas de frequência
freq_bands = [4 8; 8 13];
band_names = {'theta', 'alpha'};

for i = 1:length(intervals) %indo de 1 a 20, percorrendo todos os intervalos maiores.
    num_blocks = length(intervals{i}.band_power); 
    % Número de blocos de 5 segundos já existentes, no caso seu valor sera a quantidade de blocos de 5 segundos existentes em cada bloco maior.

    % Inicializando a célula para armazenar a RBP de cada bloco
    intervals{i}.power_ratio = cell(1, num_blocks);

    for k = 1:num_blocks % Valor que vai de 1 a quantidade de intervalos de 5 segundos que existe em cada bloco
        % Extrair o bloco de 5 segundos atual
        data_block = intervals{i}.band_power{k}; %percorre cada band power de todos os intervalos de 5 segundos.

        % Inicializar a estrutura para armazenar o Power Ratio para o bloco atual
        power_ratio_block = struct();

        % Cálculo do Power Ratio para cada banda
        for b = 1:size(freq_bands, 1) % Percorre cada uma das frequency bands
            for a = 1:size(EEG_block, 1) % Era 63, porém pode alternar se removermos certos canais, por isso tem q ter um tamanho dinâmico
                power_ratio_block.ratio(a) = data_block.theta(a) / data_block.alpha(a);
            end
        end

        % Armazenar o Power Ratio calculado no bloco atual
        intervals{i}.power_ratio{k} = power_ratio_block;
    end
end

disp('Processamento da parte 6 completo!');

%% 7) Spectral Entropy of PSD %%

fs = 500;  % Frequência de amostragem
freq_bands = [1 4; 4 8; 8 13; 13 30; 30 50];
band_names = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};

for i = 1:length(intervals)
    intervals{i}.psd_entropy = struct();

    % Inicializar células para armazenar entropia espectral para cada banda
    for b = 1:length(band_names)
        intervals{i}.psd_entropy.(['entropy_' band_names{b}]) = cell(1, length(intervals{i}.band_power)); 
    end

    % Loop pelos blocos de 5 segundos
    for k = 1:length(intervals{i}.band_power)
        EEG_block = intervals{i}.blocks{k};

        % Loop pelos canais
        for ch = 1:size(EEG_block, 1)
            % Calcular a PSD para o canal e bloco atuais
            [psd, freq] = pwelch(EEG_block(ch, :), [], [], [], fs);

            % Loop por cada banda de frequência para calcular a entropia espectral
            for b = 1:length(freq_bands)
                % Encontrar os índices das frequências dentro da banda atual
                band_idx = find(freq >= freq_bands(b, 1) & freq <= freq_bands(b, 2));
                
                % Obter o PSD na banda específica
                band_psd = psd(band_idx);

                % Normalizar o PSD da banda
                norm_psd = band_psd / sum(band_psd);

                % Calcular a entropia espectral (evitando log(0) com a função eps)
                spectral_entropy = -sum(norm_psd .* log2(norm_psd + eps));

                % Armazenar a entropia espectral na estrutura intervals
                intervals{i}.psd_entropy.(['entropy_' band_names{b}]){k}(ch) = spectral_entropy;
            end
        end
    end
end

% Exibir os resultados de entropia espectral para verificação
% for i = 1:length(intervals)
%     disp(['Intervalo ' num2str(i)]);
%     for k = 1:length(intervals{i}.blocks)
%         disp(['  Bloco ' num2str(k)]);
%         for b = 1:length(band_names)
%             entropy_values = intervals{i}.(['entropy_' band_names{b}]){k};
%             disp(['    ' band_names{b} ' Entropy: ' num2str(entropy_values)]);
%         end
%     end
% end

disp('Processamento da parte 7 completo!');

%% 8) Spectral Edge Frequency (SEF) %%

% Configurações de parâmetros
fs = 500;  % Frequência de amostragem, ajustar conforme necessário
sef_percent = 0.5;  % Porcentagem da SEF, por exemplo, 50%

% Loop pelos intervalos
for i = 1:length(intervals)
    % Inicializar a estrutura SEF dentro de intervals
    intervals{i}.SEF = struct();

    % Loop pelos blocos de 5 segundos
    for k = 1:length(intervals{i}.band_power)
        EEG_block = intervals{i}.blocks{k};

        % Inicializar uma estrutura para armazenar a SEF por canal
        block_sef = struct();

        % Loop pelos canais
        for ch = 1:size(EEG_block, 1)
            % Calcular a PSD para o canal e bloco atuais
            [psd, freq] = pwelch(EEG_block(ch, :), [], [], [], fs);

            % Calcular a potência acumulada
            cum_power = cumsum(psd) / sum(psd);

            % Encontrar o índice da frequência onde a potência acumulada atinge 50%
            sef_idx = find(cum_power >= sef_percent, 1);

            % Obter a SEF correspondente (frequência)
            sef_value = freq(sef_idx);

            % Armazenar a SEF na estrutura block_sef
            block_sef(ch).SEF_50 = sef_value;
        end

        % Armazenar a SEF do bloco atual na estrutura SEF do intervalo
        intervals{i}.SEF.(['block_' num2str(k)]) = block_sef;
    end
end

% Exibir os resultados do SEF para verificação
% for i = 1:length(intervals)
%     disp(['Intervalo ' num2str(i)]);
%     for k = 1:length(intervals{i}.blocks)
%         disp(['  Bloco ' num2str(k)]);
%         for ch = 1:size(intervals{i}.blocks{k}, 1)
%             sef_value = intervals{i}.SEF.(['block_' num2str(k)])(ch).SEF_50;
%             disp(['    Canal ' num2str(ch) ': SEF 50% = ' num2str(sef_value) ' Hz']);
%         end
%     end
% end

disp('Processamento da parte 8 completo!');
