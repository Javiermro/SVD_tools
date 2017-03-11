function FuncCreateStrainSnapShotsHardeningFromMatFile(snpFolder, nE, nPos)

cd(snpFolder)

% if exist([snpFolder,'Strains.mat']) && exist([snpFolder,'Inelas.mat'])
%     load('Strains.mat') % File that stores all Strains at every GP's
%     strainSnp=STRAINS; clear STRAINS
%     load('Inelas.mat') % File that stores the flag for inelastic steps
%     macroEnergySnp=INELAS; clear INELAS
if exist([snpFolder,'SnapStrain.mat']) && exist([snpFolder,'SnapEnergy.mat']) && exist([snpFolder,'Snapflag.mat'])
    load('SnapStrain.mat') % File that stores all Strains at every GP's
    strainSnp=SnapStrain; clear SnapStrain
%     load('SnapEnergy.mat') % File that stores the flag for inelastic steps
%     macroEnergySnp=SnapEnergy; clear SnapEnergy
    load('Snapflag.mat') % File that stores the flag for inelastic steps
    flagSnp = Snapflag ; clear Snapflag % uso la variable flag que viene de Inelastic.mat
else
    error('binary files not detected, please check if the mat files are already created!')
end

display(' ')
display('Size of Training Matrix')
% display(['-> Energy:       ' num2str(size(macroEnergySnp,1)) 'x' num2str(size(macroEnergySnp,2))])
display(['-> flag   :     ' num2str(size(flagSnp,1)) 'x' num2str(size(flagSnp,2))])
display(['-> strain :     ' num2str(size(strainSnp,1)) 'x' num2str(size(strainSnp,2))])

%% Sizes
% make test on sizes
nTraj=1;
nSnp =size(strainSnp,2);
nTS  =nSnp/nTraj;
nStrainS=size(strainSnp,1);

%TestGP=1;
%macroEnergySnp = macroEnergySnp(TestGP,:);

%% Create Elastic Snapshot matrix
ElastTraj=1:nTraj;
% number of mode taken for each trajectory (starting at the first time step)
nElastModesPerTraj=nE;
%if(ElastTraj*nElastModesPerTraj<3)
%  warning('Not enough elastic snapshots. Changed to 3.')
%  nElastModesPerTraj=3;
%end

elasticStrainSnp=zeros(nStrainS, size(ElastTraj,2)*nElastModesPerTraj);

for i=1:size(ElastTraj,2)
    it=ElastTraj(i);
    a=(it-1)*nTS+1;
    b=(i-1)*nElastModesPerTraj+1;
    
    elasticStrainSnp(:,b:(b+nElastModesPerTraj-1))=strainSnp(:,a:(a+nElastModesPerTraj-1));
    
    aCheckDissip=a+nElastModesPerTraj-1;
%     if(macroEnergySnp(aCheckDissip)>0.0)
%     if(flagSnp(aCheckDissip)>0.0)
%         warning(['You dissipate energy in your elastic modes you twat !!!!! in traj: '  num2str(it)  '  value: ' num2str(dissipSnp(aCheckDissip))])        
%     end
    
end

display('Elastic Snapshot matrix')
display(['-> Number of trajectory                : ' num2str(size(ElastTraj,2))])
display(['-> Number of time steps per trajectory : ' num2str(nElastModesPerTraj)])
display(['-> Number of snapshots                 : ' num2str(size(ElastTraj,2)*nElastModesPerTraj)])
display(['-> Size of the strain matrix           : ' num2str(size(elasticStrainSnp,1)) 'x' num2str(size(elasticStrainSnp,2))])

%% INELASTIC SNAPSHOTS   
if nPos~=0 % only elastic modes are considered 
    display('Inelastic Snapshot matrix')
    display(['-> Select ' num2str(nPos) ' Inelastic snapshots (POS)'])
    inelasPosStrainSnp=[];
    
    for i=1:nTraj

    %     iBif=min(find(macroEnergySnp((i-1)*nTS+1:i*nTS)==1))+(i-1)*nTS;
        iBif=min(find(flagSnp((i-1)*nTS+1:i*nTS)==1))+(i-1)*nTS;
        %iBifs(1,iBif)=dissipSnp(1,iBif);
        %iDis=min(find(dissipSnp((i-1)*nTS+1:i*nTS)>1e-18))+(i-1)*nTS;
        lastTS=i*nTS;
        if((lastTS-iBif)<nPos)
          warning(['Traj ' num2str(i) ': ' num2str(lastTS-iBif+1) ' TS for ' num2str(nPos) ' SNP of POST'])
        else 
          display(['---> Traj ' num2str(i) ': ' num2str(nPos/(lastTS-iBif)*100.0) '% of POST'])
        end

        % In the case of have less TS than the imposed in the previous routine
        if (lastTS-iBif)<nPos
            %tabPos=(iBif+1:lastTS) ;
            tabPos=(iBif:lastTS) ;
            warning(['TAKING ' num2str(length(tabPos)) ' POST Snapshots']);
        else
            tabPos=ceil([1:nPos]*(lastTS-iBif)/nPos+iBif);
        end

        for iPos=1:length(tabPos)
          %a=nPos-iPos;
          b=tabPos(iPos);
          %inelasPosEnergySnp(:,nPos*i-a)=EnergySnp(:,b);
          inelasPosStrainSnp = [inelasPosStrainSnp strainSnp(:,b)];
        end     

    end   

    inelasStrainSnp=[inelasPosStrainSnp];
    save(['inelasPosStrainSnp_' num2str(nPos)], 'inelasPosStrainSnp')
    save('inelasStrainSnp', 'inelasStrainSnp')
    
end

save('elasticStrainSnp', 'elasticStrainSnp')


