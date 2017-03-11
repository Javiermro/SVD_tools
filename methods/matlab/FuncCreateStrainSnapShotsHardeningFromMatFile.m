function FuncCreateStrainSnapShotsHardeningFromMatFile(snpFolder, nE, nPos)

cd(snpFolder)

% if exist([snpFolder,'Strains.mat']) && exist([snpFolder,'Inelas.mat'])
%     load('Strains.mat') % File that stores all Strains at every GP's
%     strainSnp=STRAINS; clear STRAINS
%     load('Inelas.mat') % File that stores the flag for inelastic steps
%     macroEnergySnp=INELAS; clear INELAS
% else
%     error('binary files not detected, please check if the mat files are already created!')
% end
if exist([snpFolder,'SnapStrain.mat']) && exist([snpFolder,'SnapEnergy.mat']) && exist([snpFolder,'Snapflag.mat'])
    load('SnapStrain.mat') % File that stores all Strains at every GP's
    strainSnp=SnapStrain; clear SnapStrain
    load('SnapEnergy.mat') % File that stores the flag for inelastic steps
    energySnp_e = SnapEnergy_e     ; clear SnapEnergy_e
    energySnp_v = SnapEnergy_e_vol ; clear SnapEnergy_e_vol
    energySnp_d = SnapEnergy_e_dev ; clear SnapEnergy_e_dev
    energySnp_p = SnapEnergy_p     ; clear SnapEnergy_p
    energySnp_t = energySnp_e + energySnp_p ;
    load('Snapflag.mat') % File that stores the flag for inelastic steps
    flagSnp = Snapflag ; clear Snapflag % uso la variable flag que viene de Inelastic.mat
else
    error('binary files not detected, please check if the mat files are already created!')
end

display('Size of Training Matrix')
display(['-> Strain:       ' num2str(size(strainSnp,1)) 'x' num2str(size(strainSnp,2))])
display(['-> Energy: Psi_e ' num2str(size(energySnp_e,1)) 'x' num2str(size(energySnp_e,2))])
display(['-> Energy: Psi_v ' num2str(size(energySnp_v,1)) 'x' num2str(size(energySnp_v,2))])
display(['-> Energy: Psi_d ' num2str(size(energySnp_d,1)) 'x' num2str(size(energySnp_d,2))])
display(['-> Energy: Psi_p ' num2str(size(energySnp_p,1)) 'x' num2str(size(energySnp_p,2))])
display(['-> Energy: Psi_t ' num2str(size(energySnp_t,1)) 'x' num2str(size(energySnp_t,2))])

%% Sizes
% make test on sizes
nTraj=1;
nSnp =size(strainSnp,2);
nTS  =nSnp/nTraj;
nStrainS=size(strainSnp,1);
nEnegyS_e =size(energySnp_e,1);
nEnegyS_v =size(energySnp_v,1);
nEnegyS_d =size(energySnp_d,1);
nEnegyS_p =size(energySnp_p,1);
nEnegyS_t =size(energySnp_t,1);

%TestGP=1;
%macroEnergySnp = macroEnergySnp(TestGP,:);

%% Create Elastic Snapshot matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ElastTraj=1:nTraj;
% number of mode taken for each trajectory (starting at the first time step)
% nElastModesPerTraj=nE;
%if(ElastTraj*nElastModesPerTraj<3)
%  warning('Not enough elastic snapshots. Changed to 3.')
%  nElastModesPerTraj=3;
%end

elasticStrainSnp   = zeros(nStrainS , size(ElastTraj,2)*nE);
elasticEnergySnp_e = zeros(nEnegyS_e, size(ElastTraj,2)*nE); %JLM
elasticEnergySnp_v = zeros(nEnegyS_v, size(ElastTraj,2)*nE); %JLM
elasticEnergySnp_d = zeros(nEnegyS_d, size(ElastTraj,2)*nE); %JLM
elasticEnergySnp_p = zeros(nEnegyS_p, size(ElastTraj,2)*nE); %JLM
elasticEnergySnp_t = zeros(nEnegyS_t, size(ElastTraj,2)*nE); %JLM

% for i=1:size(ElastTraj,2)
%     it=ElastTraj(i);
%     a=(it-1)*nTS+1;
%     b=(i-1)*nElastModesPerTraj+1;     
%     elasticStrainSnp(:,b:(b+nElastModesPerTraj-1))   = strainSnp(:,a:(a+nElastModesPerTraj-1));
%     elasticEnergySnp_e(:,b:(b+nElastModesPerTraj-1)) = energySnp_e(:,a:(a+nElastModesPerTraj-1));  %JLM
%     elasticEnergySnp_v(:,b:(b+nElastModesPerTraj-1)) = energySnp_v(:,a:(a+nElastModesPerTraj-1));  %JLM
%     elasticEnergySnp_d(:,b:(b+nElastModesPerTraj-1)) = energySnp_d(:,a:(a+nElastModesPerTraj-1));  %JLM
%     elasticEnergySnp_p(:,b:(b+nElastModesPerTraj-1)) = energySnp_p(:,a:(a+nElastModesPerTraj-1));  %JLM     
%     aCheckDissip=a+nElastModesPerTraj-1;
%     if(nElastModesPerTraj > size(find(flagSnp==0),2))
%         warning(['You dissipate energy in your elsatic modes you twat !!!!! in traj: '  num2str(it)  '  value: ' num2str(dissipSnp(aCheckDissip))])        
%     end    
% end

if(nE > size(find(flagSnp==0),2))
    warning(['You dissipate energy in your elsatic modes you twat !!!!! in traj: '  num2str(it)  '  value: ' num2str(dissipSnp(aCheckDissip))])        
end   
nElasSnp = size(find(flagSnp==0),2);
SelSnp   = round([1:nE]*nElasSnp/nE);
elasticStrainSnp(:,[1:nE]) = strainSnp(:,SelSnp);
elasticEnergySnp_e(:,[1:nE]) = energySnp_e(:,SelSnp);
elasticEnergySnp_v(:,[1:nE]) = energySnp_v(:,SelSnp);
elasticEnergySnp_d(:,[1:nE]) = energySnp_d(:,SelSnp);
elasticEnergySnp_p(:,[1:nE]) = energySnp_p(:,SelSnp);
elasticEnergySnp_t(:,[1:nE]) = energySnp_t(:,SelSnp);

display('Elastic Snapshot matrix')
display(['-> Number of trajectory                : ' num2str(size(ElastTraj,2))])
display(['-> Number of time steps per trajectory : ' num2str(nE)])
display(['-> Number of snapshots                 : ' num2str(size(ElastTraj,2)*nE)])
display(['-> Size of the strain matrix           : ' num2str(size(elasticStrainSnp,1)) 'x' num2str(size(elasticStrainSnp,  2))])
display(['-> Size of the energy matrix, Psi_e    : ' num2str(size(elasticEnergySnp_e,1)) 'x' num2str(size(elasticEnergySnp_e,2))])
display(['-> Size of the energy matrix, Psi_v    : ' num2str(size(elasticEnergySnp_v,1)) 'x' num2str(size(elasticEnergySnp_v,2))])
display(['-> Size of the energy matrix, Psi_d    : ' num2str(size(elasticEnergySnp_d,1)) 'x' num2str(size(elasticEnergySnp_d,2))])
display(['-> Size of the energy matrix, Psi_p    : ' num2str(size(elasticEnergySnp_p,1)) 'x' num2str(size(elasticEnergySnp_p,2))])
display(['-> Size of the energy matrix, Psi_t    : ' num2str(size(elasticEnergySnp_t,1)) 'x' num2str(size(elasticEnergySnp_t,2))])


%% Create INELASTIC snapshots matrix
display('Inelastic Snapshot matrix')
display(['-> Select ' num2str(nPos) ' Inelastic snapshots (POS)'])

% inelasStrainSnp   = [];
% inelasEnergySnp_e = [];
% inelasEnergySnp_v = [];
% inelasEnergySnp_d = [];
% inelasEnergySnp_p = [];
% for i=1:nTraj    
%     iPlas=min(find(flagSnp((i-1)*nTS+1:i*nTS)==1))+(i-1)*nTS;  % identify the begining of the plastification 
%     if nPos==0 ; iPlas=0 ; end %JLM
%     %iBifs(1,iBif)=dissipSnp(1,iBif);
%     %iDis=min(find(dissipSnp((i-1)*nTS+1:i*nTS)>1e-18))+(i-1)*nTS;
%     lastTS=i*nTS; %JLM -1
%     if((lastTS-iPlas)<nPos)
%       warning(['Traj ' num2str(i) ': ' num2str(lastTS-iPlas+1) ' TS for ' num2str(nPos) ' SNP of POST'])
%     else 
%       display(['---> Traj ' num2str(i) ': ' num2str(nPos/(lastTS-iPlas)*100.0) '% of POST'])
%     end
%     % In the case of have less TS than the imposed in the previous routine
%     if (lastTS-iPlas)<nPos
%         %tabPos=(iBif+1:lastTS) ;
%         tabPos=(iPlas:lastTS) ;
%         warning(['TAKING ' num2str(length(tabPos)) ' POST Snapshots']);
%     else
% %         tabPos=ceil([1:nPos]*(lastTS-iBif)/nPos+iBif); %JLM lo cambie 
%         tabPos=(iPlas:nPos+iPlas-1) ;
%     end
%     
%     
%     for iPos=1:length(tabPos)
%       %a=nPos-iPos;
%       b=tabPos(iPos);
%       %inelasPosEnergySnp(:,nPos*i-a)=EnergySnp(:,b);      
%       inelasStrainSnp   = [inelasStrainSnp   strainSnp(:,b)];
%       inelasEnergySnp_e = [inelasEnergySnp_e energySnp_e(:,b)]; %JLM
%       inelasEnergySnp_v = [inelasEnergySnp_v energySnp_v(:,b)]; %JLM
%       inelasEnergySnp_d = [inelasEnergySnp_d energySnp_d(:,b)]; %JLM
%       inelasEnergySnp_p = [inelasEnergySnp_p energySnp_p(:,b)]; %JLM
%     end     
%     
% end   


InelasPos  = find(flagSnp==1) ;
nInelasSnp = size(InelasPos,2);
if nInelasSnp<nPos
    warning(['---> there are less inelastic snapshots than requiered'])
else
    InelasPos(1:nPos) ;
    SelSnp   = InelasPos*nInelasSnp/nPos;
    inelasStrainSnp   = strainSnp(:,InelasPos(1:nPos)) ;
    inelasEnergySnp_e = energySnp_e(:,InelasPos(1:nPos));
    inelasEnergySnp_v = energySnp_v(:,InelasPos(1:nPos)); 
    inelasEnergySnp_d = energySnp_d(:,InelasPos(1:nPos)); 
    inelasEnergySnp_p = energySnp_p(:,InelasPos(1:nPos)); 
    inelasEnergySnp_t = energySnp_t(:,InelasPos(1:nPos));     
end



%% STRAIN SNAPSHOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inelasStrainSnp=[inelasPosStrainSnp];
% save(['inelasPosStrainSnp_' num2str(nPos)], 'inelasPosStrainSnp')
save('elasticStrainSnp', 'elasticStrainSnp')
save('inelasStrainSnp', 'inelasStrainSnp')

%% ENERGY SNAPSHOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inelasEnergySnp=[inelasPosEnergySnp];
% save(['inelasPosEnergySnp_' num2str(nPos)], 'inelasPosEnergySnp')
% save('inelasEnergySnp', 'inelasEnergySnp')
% save('elasticEnergySnp', 'elasticEnergySnp')
elasticEnergySnp = elasticEnergySnp_e; save('Phi_ela_elastic', 'elasticEnergySnp'); clear elasticEnergySnp;
elasticEnergySnp = elasticEnergySnp_v; save('Phi_vol_elastic', 'elasticEnergySnp'); clear elasticEnergySnp;
elasticEnergySnp = elasticEnergySnp_d; save('Phi_dev_elastic', 'elasticEnergySnp'); clear elasticEnergySnp;
elasticEnergySnp = elasticEnergySnp_p; save('Phi_pla_elastic', 'elasticEnergySnp'); clear elasticEnergySnp;
elasticEnergySnp = elasticEnergySnp_t; save('Phi_tot_elastic', 'elasticEnergySnp'); clear elasticEnergySnp;

inelasEnergySnp = inelasEnergySnp_e; save('Phi_ela_inelastic', 'inelasEnergySnp'); clear inelasEnergySnp ;
inelasEnergySnp = inelasEnergySnp_v; save('Phi_vol_inelastic', 'inelasEnergySnp'); clear inelasEnergySnp ;
inelasEnergySnp = inelasEnergySnp_d; save('Phi_dev_inelastic', 'inelasEnergySnp'); clear inelasEnergySnp ;
inelasEnergySnp = inelasEnergySnp_p; save('Phi_pla_inelastic', 'inelasEnergySnp'); clear inelasEnergySnp ;
inelasEnergySnp = inelasEnergySnp_t; save('Phi_tot_inelastic', 'inelasEnergySnp'); clear inelasEnergySnp ;

 display('-')
