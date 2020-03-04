%%% This will get ERPs and/or ERSPs

clc;clear;

%% Parameters
origfilepath={'/Users/ajsimon/Documents/Data/Cortica_EEG/Evo/Merged/Post/NR/'}; %Folder where all the data is stored
replacefilepath={'/Users/ajsimon/Documents/Data/Cortica_EEG/Evo/ERP_ERSP/Post/NR/'};
filetype='.set'; %Ending- all raw data should be labelled ending from 1 up to 9. The _processed is output from FilterEpochEEG

Run_ERP=true;
Run_ERSP=true;

%allevents={'128' , '32'};
allevents = {'240'};
exptypes={'TOVA'}; %This should be the same length as allevents

%eventnames={'Targets' , 'Non-Targets'};
eventnames = {'all'};

eeglab; %Start EEGLAB
pop_editoptions( 'option_single', false);
%% Now let's Process!

for fp = 1:length(origfilepath)
fileList = getAllFiles(origfilepath{fp});
Fileidx=strfind(fileList,filetype); Fileidx=find(~cellfun(@isempty,Fileidx)); %Find the cleaned set Files
fileList=fileList(Fileidx); %Extract just these files

for fn=1:length(fileList)
    
    pathsep=strfind(fileList{fn},filesep);pathsep=pathsep(end); %find where the directory ends and filename starts
    filename=fileList{fn}(pathsep+1:end-length(filetype));
    
    filepath=fileList{fn}(1:pathsep);
    newfilepath=strrep(filepath, origfilepath{fp}, replacefilepath{fp});
    
    if ~exist(newfilepath,'dir');
        mkdir(newfilepath);
    end
    
    if ~exist([newfilepath filename '_Targets_ERPs.mat'],'file') && ~exist([newfilepath filename '_Non-Targets_ERPs.mat'],'file')
        EEG=pop_loadset([filename filetype],filepath);
        %% ERP section
        for eco=1:length(allevents)
            ERP=struct;
            indx=0;
            for e=1:length(EEG.epoch)
                if length(EEG.epoch(1,e).event)>1
                    match=strmatch(allevents{eco},EEG.epoch(1,e).eventtype{1,1});
                    if match==1
                        indx=indx+1;
                        tempdata(:,:,indx)=EEG.data(:,:,e);
                    end
                elseif length(EEG.epoch(1,e).event)==1
                    match=strmatch(allevents{eco},EEG.epoch(1,e).eventtype);
                    if match==1
                        indx=indx+1;
                        tempdata(:,:,indx)=EEG.data(:,:,e);
                    end
                end
            end
        
        %% compute ERP
        if Run_ERP==true
            
            if indx>0
                
                ERP.data=mean(tempdata(:,:,:),3);
                ERP.numtrials=indx;
                
                if ~exist([newfilepath filename '_' eventnames{eco} '_ERPs.mat'],'file');
                    
                    [ERP.P1amp, ERP.P1lat] = fPickPeakWin(ERP.data,EEG.times,[65 140],'max',5);   %65 to 140 msec
                    [ERP.N1amp, ERP.N1lat] = fPickPeakWin(ERP.data,EEG.times,[130 230],'min',5);   %130 to 230 msec
                    [ERP.P2amp, ERP.P2lat] = fPickPeakWin(ERP.data,EEG.times,[150 270],'max',5);   %150 to 270 msec
                    [ERP.P300amp, ERP.P300lat] = fPickPeakWin(ERP.data,EEG.times,[250 500],'max',5);  %250 to 500 msec
                    %
                    ERP.chanlocs=EEG.chanlocs;
                    ERP.channels={ERP.chanlocs.labels};
                    ERP.times=EEG.times;
                    
                    %get P3b stats from Pz
                    tindx = find(ERP.data(31,1281:1589) > 0 );   %200-500 msec
                    col_ix=1; row_ix=1;
                    if ~isempty(tindx)
                        new_tindx(row_ix,col_ix)=tindx(1,1);
                    end
                    for ti=2:length(tindx)
                        if tindx(ti)==tindx(ti-1)+1
                            col_ix=col_ix+1;
                            new_tindx(row_ix,col_ix)=tindx(ti);
                        else
                            row_ix=row_ix+1; col_ix=1;
                            new_tindx(row_ix,col_ix)=tindx(ti);
                        end
                    end
                    if ~isempty(tindx)
                        for p=1:size(new_tindx,1)
                            tbeg=new_tindx(p,1)+1281;
                            tend=max(new_tindx(p,:))+1281;
                            tempAUC(1,p) = trapz(ERP.data(31,tbeg:tend));
                        end
                        ERP.P3bAUC(1,1)=sum(tempAUC);
                    else
                        ERP.P3bAUC(1,1) = 0;
                    end
                    
                    save([newfilepath filename '_' eventnames{eco} '_ERPs.mat'],'ERP');
                    clear new_tindx tindx tempAUC tbeg tend
                    numevents{fn,1}=filename(1:3);
                    numevents{fn,2}=size(tempdata,3);
                end
            else
                fprintf(['No ERPs generated for: '  filename '_' eventnames{eco} '_ERPs.mat, \n'  eventnames{eco} ' stimulus onset event markers do not exist.\n']);
            end
        end
        
        %% Compute ERSP
        if Run_ERSP==true
            
            if ~exist([newfilepath filename '_' eventnames{eco} '_ERSP.mat'],'file') && indx>0
                frames=size(EEG.data,2);
                tlimits=[EEG.times(1), EEG.times(end)];
                ERSPs=struct('data',[],'itc',[],'powbase',[],'times',[],'freqs',[],'itcphase',[]);
                ERSPs(size(EEG.data,1)).data=[];
                
                % This commented portion is the normal chunk of
                % code to use in this script that will compute the
                % ERSP on all elecs. I commented it out on 2/13/18
                % so I could just get the ERSP in two electrodes
                for elec=1:size(EEG.data,1)
                    disp(['Now Processing Electrode ' num2str(elec)]);
                    data=squeeze(tempdata(elec,:,:)); data=data(:)';
                    [ERSPs(elec).data,ERSPs(elec).itc,ERSPs(elec).powbase,ERSPs(elec).times,ERSPs(elec).freqs,~,~,ERSPs(elec).itcphase]= ...
                        timef(data,frames,tlimits,EEG.srate,0,'timesout',400,'maxfreq',40,'baseline',-700,'plotersp','off','itctype','coher','plotitc','off','freqscale','linear');
                    ERSPs(elec).data=ERSPs(elec).data(1:20,:);
                    ERSPs(elec).freqs=ERSPs(elec).freqs(1,1:20);
                    ERSPs(elec).itc=ERSPs(elec).itc(1:20,:);
                    ERSPs(elec).itcphase=ERSPs(elec).itcphase(1:20,:);
                end
                
                %                         data=squeeze(tempdata(38,:,:)); data=data(:)';
                %                         [ERSPs(1).data,ERSPs(1).itc,ERSPs(1).powbase,ERSPs(1).times,ERSPs(1).freqs,~,~,ERSPs(1).itcphase]= ...
                %                             timef(data,frames,tlimits,EEG.srate,0,'timesout',400,'maxfreq',40,'baseline',-700,'plotersp','off','itctype','coher','plotitc','off','freqscale','linear');
                %                         ERSPs(1).data=ERSPs(1).data(1:20,:);
                %                         ERSPs(1).freqs=ERSPs(1).freqs(1,1:20);
                %                         ERSPs(1).itc=ERSPs(1).itc(1:20,:);
                %                         ERSPs(1).itcphase=ERSPs(1).itcphase(1:20,:);
                %
                %                         ERSP=struct;
                %                         ERSP.chanlocs=EEG.chanlocs;
                %                         ERSP.channels={ERSP.chanlocs.labels};
                %                         ERSP.times=ERSPs(1).times;
                %                         ERSP.freqs=ERSPs(1).freqs;
                %
                %                         ERSP.data(1,:,:)=ERSPs(1).data;
                %                         ERSP.itc(1,:,:)=ERSPs(1).itc;
                
                for el=1:length(ERSPs)
                    ERSP.data(el,:,:)=ERSPs(el).data;
                    ERSP.itc(el,:,:)=ERSPs(el).itc;
                end
                
                ERSP.chanlocs=EEG.chanlocs;
                ERSP.channels={ERSP.chanlocs.labels};
                ERSP.times=ERSPs(1).times;
                ERSP.freqs=ERSPs(1).freqs;
                
                save([newfilepath filename '_' eventnames{eco} '_ERSP.mat'],'ERSP');
                clear ERSP;
            end
        end
        clear indx tempdata
        end
    end
end
end