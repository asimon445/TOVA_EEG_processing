% This script will epoch data, reject epochs containing artifacts (above
% threshold of 7.5 µV), and baseline correct. 

clear;clc;

eeglab;

%% Parameters
origfilepath='/Volumes/Seagate Backup Plus Drive/Cortica_EEG/AJ_pipeline/PediBBT/TOVA/preprocessed/Post/impulsive/'; %Folder where all the data is stored
replacefilepath='/Volumes/Seagate Backup Plus Drive/Cortica_EEG/AJ_pipeline/PediBBT/TOVA/epcohed_DataRecovery/Post/impulsive/';
filetype='_preprocessed.set'; %This is based on whatever is added to the end of the PreProcessEEG script default: _cleaned.set
ix=0;

fileList = getAllFiles(origfilepath);
Fileidx=strfind(fileList,filetype); Fileidx=find(~cellfun(@isempty,Fileidx)); %Find the cleaned set Files
fileList=fileList(Fileidx); %Extract just these files


allevents={'32','128'};  %stimulus onset markers for TOVA
lim=[-1 1]; %Set the epoch size, about stimulus, pre and post stimulus
timerange = [-200 0];  %timerange to baseline correct to

IndRej=zeros(length(fileList),250);

for fn=1:length(fileList)
    
    pathsep=strfind(fileList{fn},filesep);pathsep=pathsep(end); %find where the directory ends and filename starts
    filename=fileList{fn}(pathsep+1:end-length(filetype));
    
    filepath=fileList{fn}(1:pathsep);
    
    newfilepath=strrep(filepath, origfilepath, replacefilepath);
    
    if ~exist([newfilepath filename '_epoched.set'],'file')
        if ~exist(newfilepath,'dir');
            mkdir(newfilepath);
        end
        
        EEG=pop_loadset([filename filetype],filepath);
        
        % This section will check if any photodiodes have not been
        % converted. They should have been by this point, but if not,
        % this will make sure they are
        
%         nonconv=zeros(1,length(EEG.event)); conv=[];
%         num=zeros(1,length(EEG.event));
%         
%         for e=1:length(EEG.event)
%             concheck=strfind(EEG.event(1,e).type,'8192');
%             if isempty(concheck)
%                 nonconv(1,e)=0;
%             else
%                 nonconv(1,e)=concheck;
%             end
%         end
%         
%         conv=find(nonconv==1);
%         if length(conv)>1
%             EEG=PD_to_numeric(EEG);
%         end
%         
%         z_ix=0;
%         for e=1:length(EEG.event)
%             num=isnumeric(EEG.event(1,e).type);
%             if num==1
%                 EEG.event(1,e).type=num2str(EEG.event(1,e).type);
%             end
%             zro=strmatch(EEG.event(1,e).type,'0');
%             if ~isempty(zro)
%                 z_ix=z_ix+1;
%             end
%             
%             %recode mismatches
%             one60=strmatch(EEG.event(1,e).type,'160');
%             if ~isempty(one60)
%                 EEG.event(1,e).type='32';
%             end
%             one44=strmatch(EEG.event(1,e).type,'144');
%             if ~isempty(one44)
%                 EEG.event(1,e).type='128';
%             end
%             
%         end
%         
%         if z_ix<5
            if length(EEG.event)>20
                ix=ix+1;
                
                EEG=pop_epoch(EEG,allevents,lim);   %epoch
                
                OriginalEpochs=length(EEG.epoch);  %index number of epochs generated before artifact rejection
                
                [EEG, Irej] = pop_eegthresh(EEG,1,[3 23:27 30 33 36 37 38 47 60 61:64] ,-150,150,-1,0.99902,0,1);  %reject artifactual epochs (fluctuations over +/- 7.5 µV) from selected channels
                
                %remove first and last trial
                en=find(Irej==OriginalEpochs);
                pr=find(Irej==1);
                
                if isempty(en) && isempty(pr)   %if neither the first nor the last events were rejected
                    EEG.data=EEG.data(:,:,2:end-1);
                    EEG.epoch=EEG.epoch(:,2:end-1);
                    EEG=eeg_checkset(EEG);
                    Irej(1,end+1)=OriginalEpochs;  %add the last event to the list of rejected events
                    Irej(2:length(Irej)+1)=Irej; Irej(1,1)=1;  %add the first event to the list of rejected events
                elseif ~isempty(pr)  %if only the first event was rejected
                    EEG.data=EEG.data(:,:,1:end-1);
                    EEG.epoch=EEG.epoch(:,1:end-1);
                    EEG=eeg_checkset(EEG);
                    Irej(1,end+1)=OriginalEpochs;  %add the last event to the list of rejected events
                elseif ~isempty(en)  %if only the last event was rejected
                    EEG.data=EEG.data(:,:,2:end);
                    EEG.epoch=EEG.epoch(:,2:end);
                    EEG=eeg_checkset(EEG);
                    Irej(2:length(Irej)+1)=Irej; Irej(1,1)=1;  %add the first event to the list of rejected events
                end
                
                IndRej(fn,1:length(Irej))=Irej(1,:);
                
                NewEpochs=length(EEG.epoch);
                RemovedEpochs=OriginalEpochs-NewEpochs;
                
                EEG=pop_rmbase(EEG,timerange);   %remove baseline
                
                ArtifactRejInfo(1,1)={'Filename'};
                ArtifactRejInfo(1,2)={'Original Epoch Count'};
                ArtifactRejInfo(1,3)={'New Epoch Count'};
                ArtifactRejInfo(1,4)={'Removed Epochs'};
                
                ArtifactRejInfo(ix+1,1)={filename};
                ArtifactRejInfo(ix+1,2)={OriginalEpochs};
                ArtifactRejInfo(ix+1,3)={NewEpochs};
                ArtifactRejInfo(ix+1,4)={RemovedEpochs};
                
                pop_saveset(EEG,'filename',[filename '_epoched'],'filepath',newfilepath);
            end
        %end
    end
end

str=date; 
outfile=sprintf('Artifact_rejection_info_TOVA_%s_01.mat',str);
save([replacefilepath outfile],'ArtifactRejInfo');  %save information about how many epochs were rejected per dataset

outfile1=sprintf('Indexed_epochs_rejected_%s_01.mat',str);
save([replacefilepath outfile1],'IndRej'); 