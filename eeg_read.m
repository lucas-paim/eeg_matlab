% Original script made by Nahian July 2021 based on Andrew Furman scripts
% Sammy adapted for Raf's project (extension of MODULATE) Sept 2021
% max clean, PAF 9-11 Hz
% eyes closed, resting states 1, 2, 3, and 4 of MODULATE

% Sammy added a more visually pleasing trial rejection - Oct 2021
% Also changed the final PAF calc section to not need hard coding of the
% frequency ranges - an index calculation based on foi is used instead

% Sammy adapted to work on desktop computer for MODULATE - April/May 2022
% See EC_modulate branch of PAF_Script.m on gitlab for commit messages
%%
% Have a look at the read.ME in gitlab to see how to run this script

% 1) Set-up workspace
% 2) EEGLAB section
% 3) FIELDTRIP section
% 4) Epoch check (checks there are >24 epochs (i.e. >2 mins) clean data)
% 5) PAF and power calculation (plus outputs averages in longformat .csv)
% 6) Functions

clear;
close all;
%% 1) Set-up Workspace %%
restoredefaultpath

cd('C:\Users\lucas\source\matlab_eeg\dados_lucasp')
cwd = [cd filesep]; 

wpms = []; % pre-allocate a workspace variables in struct

wpms.DATAIN     = fullfile(cwd, 'Datain', filesep); 
wpms.DATAOUT    = [cwd 'Dataout' filesep];
wpms.FUNCTIONS  = fullfile(cwd, 'Functions', filesep);
wpms.software   = 'C:\Users\lucas\source\matlab_eeg\dados_lucasp\Functions\';

load([wpms.FUNCTIONS filesep 'chanlocs.mat'])
load([wpms.FUNCTIONS 'neighbour_template.mat'])
neighbours = neighboursCopy;

addpath([wpms.software 'eeglab_current\eeglab2024.2'])
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\functions']));
addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\xdfimport1.19']);
addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\bva-io']);
%addpath([wpms.software 'eeglab_current\eeglab2024.2\plugins']); % add all plugins 
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\PrepPipeline0.56.0']));
addpath(genpath([wpms.software 'eeglab_current\eeglab2024.2\plugins\firfilt']));

addpath([wpms.software 'fieldtrip-20241025']);
addpath([wpms.software 'fieldtrip-20241025\external\eeglab']);

% Structure and store participant (px) codes
% pxlist = dir([wpms.DATAIN 'px-*']); % might also need a separation to reach EEG
% test only on px-02
pxlist = dir([wpms.DATAIN 'px-01*']); % might also need a separation to reach EEG

mkdir(wpms.DATAOUT); % Create Output directory
for subout = 1:length(pxlist)
    mkdir([wpms.DATAOUT pxlist(subout).name]) % Create folder for each px
end

%% 2) EEGLAB Section %%
eeglab
close all

skipped = 1; % contador para participantes cujos dados não serão processados corretamente.

for px =1:length(pxlist) % Inicia um loop para processar os dados de cada participante listado em pxlist.
    clearvars all_eeg_files EEG n_triggers triggers rejected % Limpa as variáveis das iterações anteriores para evitar conflitos de dados.
    fprintf(['\n Analysing participant: ' pxlist(px).name '\n\n']); % Exibe o nome do participante atual no console.
    all_eeg_files  = dir([wpms.DATAIN pxlist(px).name filesep 'EEG' filesep 'Raw' filesep '*eeg_php*.vhdr']);  % changed to specifically be only resting data
    % Busca todos os arquivos de EEG no diretório do participante que correspondem ao padrão de arquivos com prefixo 'eeg_php'.
    
    % another loop within initial px loop, loops through each resting state file for the current px
    for file = 1:length(all_eeg_files)
    
        clearvars EEG end_of_window end_time n_triggers start_of_window triggers rejected filename
        filename = get_filename(all_eeg_files(file).name);
        
        % code below makes sure the correct recording name is printed in the command window
        %fprintf(['\n Analysing session: ' all_eeg_files(file).name '\n\n']);
        fprintf(['\n Analysing session: ' filename '\n\n']);
       
        % Loads the raw EEG file
        %  Carrega o arquivo de EEG no formato BrainVision usando o EEGLAB.
        [EEG, ~] = pop_loadbv(all_eeg_files(file).folder, [all_eeg_files(file).name]); % Do you need to have opened and cleared eeglab?
        
        % store all the triggers/events in a variable
        for i = 1:length(EEG.event)
            triggers{i} = EEG.event(1, i).type; % Coleta o tipo de cada trigger (evento). 
        end
        % Loop que percorre os eventos de EEG e armazena os gatilhos (triggers).
        
         triggers = string(triggers); % Converte os triggers para string.

         %LUCAS: triggers are markers in the EEG file that mark when something happened. In this case each marker means:
         % R5 - heat beginning to ramp up
         % R4 - heat reached 46 degrees (pain)
         % R6 - heat begining to ramp down
         % R15 - heat reached 32 degrees (neutral temp - no pain)
         % So at this point in the script you would need to check whether it have 5*each of the different triggers.
         % Below I have a manual check to see if the correct triggers were present for the resting state. 
         % Also double check the durations between R4 and R6 are 40 seconds long, and the durations between R15 and R5 are 20 seconds long - See experimental protocol

        % pause;
        % ask for input depending
        % if 1 entered, mark file as being processed and continue
        % if 0 entered, mark file as not being processed 
        prompt = 'Please indicate whether file has 4 codes indicating start and stop of resting states (1 = yes, 0 = no):'; 
        % note: clicking on the lines in the spectra plot will tell you in
        % the command window the channel number of that line.
        % O usuário decide manualmente se os triggers corretos estão presentes. O valor 1 indica que o arquivo será processado; 0 indica que o arquivo será ignorado.

        process = input(prompt); % prompt saved in this variable
        processed{px, file} = process;
        
        % n_triggers = length(triggers(triggers == 'Start/End'));
        n_triggers = length(triggers); % count how many triggers
        % number of triggers could be used for automatic check

        % Dealing with edge cases
        % Commented out code to manually remove buffer overflow triggers
        % test = EEG.event(1,1:4);
        % test(1,5) = EEG.event(1,13); % change the 13 to whatever other
        % triggers you want to keep. You'll need to open EEG.event to look
        % EEG.event = test;
        % look at EEG.event again to check it worked. Then manually run the rest of the code
        % if number of triggers is 6, that means it only contains artifacts so that could have been used for automatic resting state check

        if process == 1

            % downsample to 500
            EEG = pop_resample(EEG, 500); % this is the slowest part. Reduz a taxa de amostragem para 500 Hz

            EEG = pop_select(EEG, 'nochannel', {'GSR','HR','RESP'}); % Remove canais relacionados à resposta galvânica da pele (GSR), frequência cardíaca (HR) e respiração (RESP).
            %EEG = pop_reref(EEG, []); % common average % skip till after the removing channels
            EEG = pop_eegfiltnew(EEG, 2, 100, 826, 0, [], 1); % Aplica um filtro passa-banda para manter frequências entre 2 e 100 Hz.
            close all;

            % LUCAS: at this point I would recommend opening up the user interface and having a visual look at the data :)

            % You can uncomment this piece of code for example
            pop_eegplot( EEG, 1, 1, 1); % Abre uma interface gráfica para visualizar os dados e verificar se os triggers estão corretos.

            % you can write 'eeglab redraw' into the command window and
            % have a click around trying to plot things (data scroll or spectra and maps) - just to orient yourself to the data.
            % if you click around in the eeglab user interface, you can then type 'EEG.history' into the command window and get it to print the code for what you just did manually.

            % below is Nahian's code for dealing with different/changing triggers, use if you need.
            % you will need to change this part with the triggers and get it to extract 'heat_on' as between R4 and R6 and 'heat_off' as between R15 and R5
            % only for EC for now, will need to change/add more another
            if sum(count(triggers, "EC")) == 1
                % all good
            %elseif sum(count(triggers,"EC")) == 2
                % If there are accidentally 2 EC triggers, we should use the second one. But this would need to be manually checked first and then manually implemented
                %c = find( triggers == "EC", 2);
                %c = c(1);
                %EEG.event(1, c).type = 'EO';
                %c = find( triggers == "Start/End", 2);
                %c = c(1);
                %EEG.event(1, c).type = 'EC';
            
            else
                eo_count = sum(count(triggers, "EO Start/End"));
                se_count = sum(count(triggers, "Start/End"));
                if eo_count == 2 % So this works for if we had "EO Start/End" for eyes open, and "Start/End" for eyes closed
                    EEG.event(1, find( triggers == "Start/End", 1)).type = 'EC'; 
                elseif eo_count == 0 && se_count == 4 % So this works for if we had "Start/End" 4 times
                    c = find( triggers == "Start/End", 3);
                    c = c(3);
                    EEG.event(1, c).type = 'EC';
                end
            end
            
            %if n_triggers == 2;
            %    EEG = pop_rmdat(EEG, {'Start'},[3 303], 0);
            %else
            %    EEG.event(1, (length(EEG.event)+1)).type = 'End';
            %    EEG.event(1, end).latency = length(EEG.times);
            %    pop_eegplot(EEG, 1, 1, 1);
            %    prompt = 'Please identify the end time';
            %    end_time = input(prompt);
            %    end_of_window = end_time - length(EEG.times)/500;
            %    start_of_window = end_of_window - 300;
            %    EEG = pop_rmdat(EEG, {'End'},[start_of_window end_of_window], 0);
            %end
            
            % This cuts the file based around the trigger we want from 3 seconds after to 183 seconds after the eyes closed trigger
            % You'll need to edit this part
            EEG = pop_rmdat(EEG, {'EC'},[3 183], 0); % Corta a janela temporal relevante para os dados de olhos fechados (EC), removendo os dados fora deste intervalo.
            %EEG = pop_rmdat(EEG, {'EC'},[0 210], 0);
            %EEG = pop_rmdat(EEG, {'R  4'},[0 40], 0);
            %pop_eegplot( EEG, 1, 1, 1);

            close all 

            % rejecting bad channels/electrodes/sensors
            % opens up two figures to view the data
            figure; pop_spectopo(EEG, 1, [0 15000], 'EEG', 'freqrange', [2 100], 'electrodes','off'); % Mostra um gráfico de espectro de potência para ajudar a identificar canais ruins.
            pop_eegplot(EEG, 1, 1, 1);
            % prompt will appear in command window we type a number to indicate which channel to remove or just hit 'Enter' on the keyboard if there are none to remove.
            prompt = 'Please identify sensor # to remove in matrix format (e.g., [12,14]):'; 
            % note: clicking on the lines in the spectra plot will tell you in the command window the channel number of that line

            rejected = input(prompt); % prompt saved in this variable
            rejected_channels{px, file} = rejected; % saves all rejected channels 

            % if statement 
            % if the number in rejected is not zero, then remove the bad
            % channel storred in the variable 'rejected'
            if size(rejected,1)~=0
                EEG = pop_select(EEG, 'nochannel',rejected);  
            end

            % re-reference to the common average of all electrodes
            EEG = pop_reref(EEG, []);

            % pops up another figure to double check what the spectra look like
            figure; pop_spectopo(EEG, 1, [0 15000], 'EEG', 'freqrange', [2 100], 'electrodes','off');

            % Saves the data 
            save([wpms.DATAOUT pxlist(px).name filesep pxlist(px).name '_' filename '_preICA_maxclean.mat'], 'EEG');% Change to PCA?
            fprintf(['\n Finished analysing session: ' filename '\n\n']);
                
            
            %close all
        else
            processed_no{skipped, 1} = all_eeg_files(file).name;
            processed_no{skipped, 2} = process; % save whether 1 or 0 for sanity check
            % saves a list of files that were not processed and would need
            % to be checked later
            skipped = skipped + 1;
            fprintf(['\n File: ' filename ' did not contain correct resting state recording triggers, so was skipped \n\n']); % where does it save to? Output folder? or just struct for now, is struct saved later?
        end
    end
    fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);
 end
% saves the variable that had been storing a list of all the rejected
% channels for all participants - useful for write up

save([wpms.DATAOUT 'rejected_channels.mat'], 'rejected_channels');
save([wpms.DATAOUT 'not_processed_list.mat'], 'processed_no');
save([wpms.DATAOUT 'processed_codes.mat'], 'processed');

% %% 3) FIELDTRIP Section %%
% % next round of processing conducted in the fieldtrip toolbox
% 
% % LUCAS: you may want to adjust parts of this pipeline if you need different epoch lengths or anything for particular analyses we want. 
% % For now I think we can leave as 5 seconds - but we might need to change to 2 seconds, because we only have 20 or 40 second sections, and so if
% % there is noise within those seconds, we might loose too many epochs (i.e. sections of data)
% % epochs (sections of data - sometimes called 'trials' as well)
% % reject epochs containing artifacts
% % run ICA
% % reject ICA components (i.e. eye blinks and saccades)
% % run fourier transform
% 
% clearvars all_eeg_files EEG n_triggers prompt rejected subjlist triggers px file
% pxlist = dir([wpms.DATAOUT 'px-*']); % Lista todos os participantes na pasta de saída
% not_2min_count = 1; % Contador para arquivos com menos de 2 minutos de dados limpos
% % Limpa variáveis e gera uma lista dos participantes cujos dados estão armazenados na pasta wpms.DATAOUT.
% 
% for px = 8% length(pxlist) % Itera sobre os participantes
% 
%     clearvars all_eeg_files cfg data data_rejected data_freq data_pruned EEG freq paf power promt n_rejected rejected temp temp2 X y z
%     fprintf(['\n Analysing participant: ' pxlist(px).name '\n\n']);
%     % Para cada participante, o código carrega seus arquivos e faz um processamento em várias etapas, limpando variáveis desnecessárias ao início de cada iteração.
% 
%     % Lista todos os arquivos EEG pré-processados para o participante atual.
%     all_eeg_files  = dir([wpms.DATAOUT pxlist(px).name filesep '*_preICA_maxclean.mat']);
% 
%     % loop over each file for current px
%     for file = 1:length(all_eeg_files)
%         clearvars cfg data data_freq data_pruned EEG freq paf power promt rejected temp temp2 X y z n_rejected filename
%         filename = get_filename_ICA(all_eeg_files(file).name);
% 
%         % Carrega o arquivo EEG específico do participante.
%         load([all_eeg_files(file).folder filesep all_eeg_files(file).name], 'EEG');
%         fprintf(['\n Analysing participant and session: ' pxlist(px).name ' ' filename '\n\n']);
% 
%         % the two toolboxes (i.e. eeglab and fieldtrip) have functions that allow you to convert between the different structures used to store the data for each
%         % EEG is the data structure for eeglab, and here we convert EEG using the function 'eeglab2fieldtrip' to a different structre that we name data
%         % Converte os dados do formato EEGLAB para o formato FieldTrip para usar as funções do FieldTrip.
%         data = eeglab2fieldtrip(EEG, 'preprocessing');
%         data.label = {EEG.chanlocs.labels};
% 
%         % notice all fieldtrip functions start with 'ft_' and require you to set up a configuration variable, which we name 'cfg'
% 
%         % cut the signal we have into epochs of 5 seconds
%         cfg = []; % set up empty configuration 
%         cfg.length = 5; % Define épocas de 5 segundos
%         cfg.overlap = 0; % Sem sobreposição
%         data=ft_redefinetrial(cfg,data); % redefine the trials, and give the function the config and data (Redefine os dados em épocas de 5 segundos) 
% 
%         % Nahian's code to reject epochs (i.e. not pretty)
%         %cfg          = [];
%         %cfg.method   = 'trial';
%         %cfg.alim = 100;
%         %data_rejected = ft_rejectvisual(cfg,data);
%         %n_rejected = length(data.trial) - length(data_rejected.trial);
% 
%         % Sammy's code to reject epochs (i.e. more complex, but pretty)
%         % Trial rejection taken from chronic itch script - S.K.Millard
%         chanSel = 1:63; % number of channels we have (may not need)
%         cfg = [];
%         cfg.channel = 'all';
%         cfg.ylim = [0 20]; % auto set scale, but can be changed manually in window once it pops up
%         cfg.blocksize = 6;  % Tamanho do bloco de visualização
%         cfg.viewmode = 'vertical'; % Visualiza os dados verticalmente
% 
%         % opens the figure to browse data
%         % saves the artifacts you highlight
%         visualArtifacts = ft_databrowser(cfg,data); % Navega pelos dados para marcar artefatos
%         % artifact sample points are stored in visualArtifacts.artfctdef.visual.artifact
% 
%         % Read out the trials that contain artifacts
%         artifacts = visualArtifacts.artfctdef.visual.artifact; % Salva os artefatos marcados
% 
%         a_trials = []; % trials with artifacts will be numbered here
%         for n=1:size(artifacts,1)
%             a_trials = [a_trials;find(data.sampleinfo(:,1) <= artifacts(n,1) ...
%                 & data.sampleinfo(:,2) >= artifacts(n,1))]; % concatenate the trials that contain artifacts
%             % Para cada artefato identificado, o código encontra as épocas (trials) que contêm esses artefatos.
%         end
% 
%         % find out which trials are artefact free
%         trials = 1:numel(data.trial); % list all the trials
%         trials_visual = find(~ismember(trials,a_trials))';  % trials without visual artifacts
% 
%         % Add information about artifact free trials to data_reref.trialinfo
%         % add a new section to the data struct with trial numbers in 
%         data.trialinfo = trials'; 
%         % the ' makes it long format (i.e. a column rather than a row)
% 
%         % add column with artifact detection index, 1 is clean, 0 is artifact
%         noVisual = ismember(1:length(trials),trials_visual)';
%         data.trialinfo = [data.trialinfo noVisual];
% 
%         % update labels - I don't think this is needed
%         %data.trialinfo_labels = {'1. TrialNo' '2. Triggers' '3. No Visual Artifact'};
% 
%         % reject artifacts that were manually highlighted
%         cfg = [];
%         cfg.artfctdef.reject = 'complete'; % default to delete complete trial
%         cfg.artfctdef.visual.artifact = visualArtifacts.artfctdef.visual.artifact;
%         data_rejected = ft_rejectartifact(cfg, data);
%         % O código rejeita (remove) as épocas que contêm artefatos.
% 
%         % store the number of epoched rejected - useful for writeup
%         n_rejected = length(data.trial) - length(data_rejected.trial);
% 
%         rest_num = 1;
% 
%         %n_rejected = length(data.trial) - length(data_rejected.trial);
%         n_clean_epochs{px, rest_num} = length(data_rejected.trialinfo);
% 
% 
% 
%         if n_rejected > 12
%             not_2mins{not_2min_count, 1} = all_eeg_files(file).name;
%             not_2mins{not_2min_count, 2} = n_rejected; % save whether 1 or 0 for sanity check
%             % saves a list of files that were not processed and would need
%             % to be checked later
%             not_2min_count = not_2min_count + 1;
%         end
% 
% 
%         %Run ICA with 15 components to observe
%         cfg = [];
%         cfg.method = 'runica'; % the standard method
%         cfg.runica.pca = 15; % pca - which is slightly different to ica
%         % data struct is temporarily saved as a variable called 'X' here
%         X = ft_componentanalysis(cfg,data_rejected);
%         % Realiza a ICA com 15 componentes principais para identificar e separar fontes independentes de sinais (como piscadas de olho).
% 
%         %Create plots to look at the components
%         %Look for topology that look like eye artifacts first then see how the signal for those components change over time to double check they are eye artifacts
%         cfg.component = [1:15];
%         cfg.layout = 'EEG1010.lay';
%         ft_topoplotIC(cfg,X); % topoplots (i.e. heads) (Mostra o topo das componentes ICA)
%         ft_databrowser(cfg,X); % browse the signal across time (Navega pelos sinais temporais das componentes)
% 
% 
%         %Reject Visual Components
%         %prompt will appear in the command window
%         prompt = 'Please Identify Component ro reject in matrix form, (i.e. [1,4]) :';
%         rejected = input(prompt); % stores components you want to reject
%         rejected_components{px, file} = rejected;
% 
%         cfg.component = rejected; % add those components to the config
%         data_postICA = ft_rejectcomponent(cfg,X); % reject components 
% 
%         %interpolate missing channels
%         % We aren't going to do this for now
%         %cfg= [];
%         %cfg.neighbours = neighbours;
%         %cfg.method = 'nearest';
%         %cfg.layout = 'EEG1010.lay';
%         %cfg.missingchannel = {neighbours((~ismember({neighbours(:).label}, data_pruned.label))).label};
%         %data_repaired = ft_channelrepair(cfg, data_pruned);
% 
%         %Fourier Transform with Hanning taper
%         cfg = [];
%         cfg.method = 'mtmfft'; % Transformada de Fourier
%         cfg.taper = 'hanning';
%         cfg.foi = 2:.20:50; %freq of interest 2-50Hz with 0.2Hz bins
%         cfg.keeptrials = 'yes'; %do for each trial/epoch 
%         data_freq = ft_freqanalysis(cfg, data_postICA); % Análise de frequências
%         close all
% 
%         % Save the variables data, data_freq etc., in a new file  called 
%         % _postICA in the output folder
%         save([wpms.DATAOUT pxlist(px).name filesep pxlist(px).name '_' filename '_postICA_maxclean.mat'],'data','data_freq', 'data_postICA','data_rejected', 'n_rejected');% Nahian - Adjust based on what you th
%         fprintf(['\n Finishined analysing participant and timepoint: ' pxlist(px).name ' ' filename '\n\n']);
% 
% EEG = fieldtrip2eeglab(data_postICA);
%     end
%     fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);
% 
%     save([wpms.DATAOUT 'n_clean_epochs.mat'], 'n_clean_epochs');
%     save([wpms.DATAOUT 'rejected_components.mat'], 'rejected_components');
% %     
% end
% save([wpms.DATAOUT 'n_clean_epochs.mat'], 'n_clean_epochs');
% save([wpms.DATAOUT 'rejected_components.mat'], 'rejected_components');
% fprintf('\n Finished ICA part of the script. \n\n');
% 
% %% 4) Count the number of trials
% % to check if any need re-doing with a longer length of resting state rather than +3 - +183 seconds after EC
% 
% clearvars all_eeg_files EEG n_triggers prompt rejected subjlist triggers px file data_postICA file filename px n_clean_epochs n_clean_epochs_check not_2mins data_freq artifacts
% clearvars not_2min_count not_2mins_check
% 
% pxlist = dir([wpms.DATAOUT 'px-*']); % list of all px in output folder
% not_2min_count = 1;
% % gera uma lista de diretórios (ou arquivos) para cada participante (px-*), que já tiveram seus dados pré-processados e salvos no diretório wpms.DATAOUT
% % A variável not_2min_count serve para contar quantos participantes tiveram menos do que 24 épocas limpas.
% 
% for px = 1:length(pxlist) % loop through participant list
%     all_eeg_files  = dir([wpms.DATAOUT pxlist(px).name filesep '*postICA_maxclean.mat']);
%     % percorre todos os participantes listados no diretório de saída e, para cada participante, carrega os arquivos que passaram pela etapa de ICA (análise de componentes independentes) e limpeza de artefatos.
% 
%     % loop through each file for the current participant
%     for file = 1:length(all_eeg_files)
%         % load file saved after fourier transform in fieldtrip
%         load([all_eeg_files(file).folder filesep all_eeg_files(file).name], 'data_freq');
% 
% 
%         rest_num = get_rest_num(all_eeg_files(file).name);
%         if rest_num == 11
%             rest_num = 1;
%         end
%         filename = get_filename_ICA(all_eeg_files(file).name);
%         fprintf(['\n Analysing participant and session: ' pxlist(px).name ' ' filename '\n\n']);
% 
%         %n_rejected = length(data.trial) - length(data_rejected.trial);
%         n_clean_epochs_check{px, rest_num} = length(data_freq.trialinfo);
% 
% 
%         if n_clean_epochs_check{px, rest_num} < 24
%             not_2mins_check{not_2min_count, 1} = all_eeg_files(rest_num).name;
%             not_2mins_check{not_2min_count, 2} = n_clean_epochs_check{px, rest_num}; % save whether 1 or 0 for sanity check
%             % saves a list of files that were not processed and would need
%             % to be checked later
%             not_2min_count = not_2min_count + 1;
%         end
%         fprintf(['\n Finishined analysing participant and timepoint: ' pxlist(px).name ' ' filename '\n\n']);
%     end
%     fprintf(['\n Finished analysing participant: ' pxlist(px).name '\n\n']);
%     save([wpms.DATAOUT 'n_clean_epochs_check.mat'], 'n_clean_epochs_check');
% 
% 
% end
% save([wpms.DATAOUT 'not_2mins_check_02052022.mat'], 'not_2mins_check');
% 
% not_2mins_check;
% %pause;
% 
% fprintf('\n Finished check of participants requiring manual restart to get more epochs. \n\n');
% 
% % Acessar a matriz 'data' dentro da estrutura
% data = EEG.data;
% 
% % Plotar o canal 1
% plot(data(1, :));
% xlabel('Amostras');
% ylabel('Amplitude');
% title('Sinal do Canal 1');
    
%% 5) Functions %%

function filename = get_filename(full_filename)
    filename = split(full_filename, '_');
    filename = filename(length(filename)); 
    filename = split(filename, '.'); 
    filename = filename{1};
end

function filename = get_filename_ICA(full_filename)
    filename = split(full_filename, '_');
    filename = filename{2}; 
    %filename = split(filename, '.'); 
    %filename = filename{1};
end

function rest_num = get_rest_num(full_filename)
    rest_num = split(full_filename, '_');
    rest_num = split(rest_num(2), 'rest'); 
    rest_num = rest_num{2};
    rest_num = str2double(rest_num);
end
