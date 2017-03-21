clc; clear all; tic;
display('------------------------------------------------------------')
display('-                       START (Step 03)                    -')
display('------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%% USERS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/javiermro/Projects/SVD_tools/methods/matlab')
offlineFolder='/home/javiermro/Projects/SVD_tools/Offline/'; cd(offlineFolder); 
ortSplit=1; % orth Singular and Regular snapshots
normaliz=0; % normalize snapshots before SVD
tol_nE = 1e-5 ; % tolerancia para determinar el nro de modos elasticos
% tol_nI = 1e-5 ; % tolerancia para determinar el nro de modos inelasticos
nE_max = 6; % number of elastic snapshots (Maximo) JLM
nPos = 83; % number of Inelastic snapshots (DEBE COINCIDIR CON STEP 02)
Dom_Descomp = 1 ; % 1= yes; 0 = no
nameExtensionW=['RVE_StructHole10_nPos' num2str(nPos)]; % Directory of Modes 
load('/home/javiermro/Projects/Examples/StructHole10/DomainPointers.mat'); % Pointer for Domain Decomposition 

nModes = 26; %number of modes or trajectories
snpFolder0 = '/home/javiermro/Projects/Examples/StructHole10/Modo';
snpFile='SNAPSHOTS_RVE_StructHole10.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

snpFolder = [];
for i=1:nModes
    if i<10
        snpFolder =[ snpFolder ; [snpFolder0 '0' num2str( i )]] ;
    else
        snpFolder =[ snpFolder ; [snpFolder0 num2str( i )]] ;
    end
end

runFolder=[offlineFolder 'Modes/']; if(~isdir(runFolder)); mkdir(runFolder); end; cd(runFolder);W=1;

if(~isdir(nameExtensionW)); mkdir(nameExtensionW); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELASTIC FULCTUATION MATRIX                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step=1;
display(['STEP ' num2str(step) ': Creation of the elastic fluctuation snapshots matrix']); step=step+1;
if(size(snpFolder,1)==1)
    cd(snpFolder);
    load('elasticStrainSnp.mat');
%     load('elasticEnergySnp.mat');
else
    nSnpFolder=size(snpFolder,1);
    cd(snpFolder(1,:));
    load('elasticStrainSnp.mat');  tmpS=elasticStrainSnp;
%     load('elasticEnergySnp.mat');  tmpE=elasticEnergySnp;
    for ifold=2:nSnpFolder
        cd(snpFolder(ifold,:));
        load('elasticStrainSnp.mat');  tmpS=[tmpS elasticStrainSnp];
%         load('elasticEnergySnp.mat');  tmpE=[tmpE elasticEnergySnp];
    end
    elasticStrainSnp = tmpS; clear tmpS;
%     elasticEnergySnp = tmpE; clear tmpE ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC SNP                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Creation of the inelastic fluctuation snapshots matrix']); step=step+1;

if(size(snpFolder,1)==1)
    cd(snpFolder);
    load(['inelasStrainSnp.mat']);
%     inelasStrainSnp=[inelasPosStrainSnp];
%     load(['inelasPosEnergySnp_' num2str(nPos)]);
%     inelasEnergySnp=[inelasPosEnergySnp];
else
    nSnpFolder=size(snpFolder,1);
    cd(snpFolder(1,:));
    load(['inelasStrainSnp.mat']);
%     inelasStrainSnp=[inelasStrainSnp];
    tmpS=inelasStrainSnp;
%     load(['inelasPosEnergySnp_' num2str(nPos)]);
%     inelasEnergySnp=[inelasPosEnergySnp];
%     tmpE=inelasEnergySnp;
    for ifold=2:nSnpFolder
        cd(snpFolder(ifold,:));
        load(['inelasStrainSnp.mat']);
%         inelasStrainSnp=[inelasPosStrainSnp];
        tmpS=[tmpS inelasStrainSnp];
%         load(['inelasPosEnergySnp_' num2str(nPos)]);
%         inelasEnergySnp=[inelasPosEnergySnp];
%         tmpE=[tmpE inelasEnergySnp];
    end
    inelasStrainSnp = tmpS; clear tmpS;
%     inelasEnergySnp = tmpE; clear tmpE;
end
%inelasStrainSnp=STRAIN_MODES_UPDATING(inelasStrainSnp, 'GidToFE');
%inelasStrainSnp=STRAIN_MODES_UPDATING(strainSnp(:,1:50), 'GidToFE');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fusion inelastic/elastic SNP                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Fusion of elastic and inelastic snapshots']); step=step+1;
nElast=size(elasticStrainSnp,2);
nInelast=size(inelasStrainSnp,2);
allStrainSnp=[elasticStrainSnp inelasStrainSnp];
% allEnergySnp=[elasticEnergySnp inelasEnergySnp];

cd(runFolder)
if(~isdir([nameExtensionW '/all'])); mkdir([nameExtensionW '/all']); end;
plot_cov_matrix(allStrainSnp, W, 1)
% plot_cov_matrix(allEnergySnp, W, 2)
% saveas(gca, [nameExtensionW '/all/elasticStrainCOV_' nameExtensionW '.eps'],'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Splitting                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Splitting of all snapshots in Regular/Singular']); step=step+1;

display('        ---> splitting')
allStrainSnp_REG=zeros(size(allStrainSnp)); allStrainSnp_SIN=zeros(size(allStrainSnp));
% allStrainSnp_REG(pointersToMatrixElem,:)=allStrainSnp(pointersToMatrixElem,:);
% allStrainSnp_SIN(pointersToCBElem,:)    =allStrainSnp(pointersToCBElem,:);
allStrainSnp_REG(PointersToSet1,:) =allStrainSnp(PointersToSet1,:);
% allStrainSnp_REG = allStrainSnp(PointersToSet1,:);% JLM
allStrainSnp_SIN(PointersToSet2,:) =allStrainSnp(PointersToSet2,:);
% allStrainSnp_SIN = allStrainSnp(PointersToSet2,:);% JLM 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELASTIC FULCTUATION SVD                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': SVD of the Regular elastic snp']); step=step+1;
Y=allStrainSnp_REG(:,1:nElast);% E=allEnergySnp_REG(:,1:nElast);
plot_cov_matrix(allStrainSnp_REG(:,1:nElast), W, 2)
% plot_cov_matrix(allEnergySnp_REG(:,1:nElast), W, 2)
display(['        ---> rank of elastic strain snapshots: ' num2str(rank(Y'*W*Y))])
% display(['        ---> rank of elastic energy snapshots: ' num2str(rank(E'*W*E))])

[V,eigElast_REG,~] = svd(Y,0);

tmp = diag(eigElast_REG); eigElast_REG=tmp(find(tmp>1e-14)); %JLM
phiElast_REG=zeros(size(Y,1),size(eigElast_REG,1));

if size(eigElast_REG,1)==0 ; display(' THERE ARE NO ELASTIC MODES ON REGULAR DOMAIN!!') ; return ;end %JLM
for i=1:size(eigElast_REG,1) ;
    phiElast_REG(:,i)=V(:,i); 
end %Y*V(:,i)/sqrt(eigElast_REG(i)); end

plot_cov_matrix(phiElast_REG, W, 3);
figure(4); semilogy(eigElast_REG, '-*'); title('Singular values - Regular domain (Elastic) (E_{SNP})');

nE = min(nE_max,size(find(eigElast_REG>tol_nE),1)) ; %JLM
display(['        ---> keep only the first ' num2str(nE) ' elastic modes']); %JLM
eigElast_REG = eigElast_REG(1:min(nE,size(eigElast_REG,1))); %JLM
phiElast_REG = phiElast_REG(:,1:min(nE,size(eigElast_REG,1))); %JLM
if(~isdir([nameExtensionW '/Regular'])); mkdir([nameExtensionW '/Regular']); end;
fid=fopen([runFolder nameExtensionW '/Regular/elasticStrainEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigElast_REG, 'descend')'./max(eigElast_REG));   fclose(fid);

if Dom_Descomp % Domain decomposition
    display(['STEP ' num2str(step) ': SVD of the Singular elastic snp']); step=step+1;
    Y=allStrainSnp_SIN(:,1:nElast);
    plot_cov_matrix(allStrainSnp_SIN(:,1:nElast), W, 5)
    display(['        ---> rank of elastic snapshots: ' num2str(rank(Y'*W*Y))])

    [V,eigElast_SIN,~] = svd(Y,0);

    tmp = diag(eigElast_SIN); eigElast_SIN=tmp(find(tmp>1e-14)); %JLM
    phiElast_SIN=zeros(size(Y,1),size(eigElast_SIN,1));

    if size(eigElast_SIN,1)==0 ; display(' THERE ARE NO ELASTIC MODES ON SINGULAR DOMAIN !!') ; return ;end %JLM
    for i=1:size(eigElast_SIN,1);
        phiElast_SIN(:,i)=V(:,i); 
    end %Y*V(:,i)/sqrt(eigElast_SIN(i)); end

    %check_compatibility(phiElast_SIN, 2e-3, 'elastic mode');
    plot_cov_matrix(phiElast_SIN, W, 6);
    figure(7); semilogy(eigElast_SIN, '-*'); title('Singular values - Discontinuous domain (Elastic) (E_{SNP})');

    nE = min(nE_max,size(find(eigElast_SIN>tol_nE),1)) ; %JLM
    display(['        ---> keep only the first ' num2str(nE) ' elastic modes']); %JLM

    eigElast_SIN = eigElast_SIN(1:min(nE,size(eigElast_SIN,1))); %JLM
    phiElast_SIN = phiElast_SIN(:,1:min(nE,size(eigElast_SIN,1))); %JLM
    % eigElast_SIN = eigElast_SIN(1:min(3 ,size(eigElast_SIN,1))); phiElast_SIN = phiElast_SIN( :,1:min(3,size(eigElast_SIN,1)));
    if(~isdir([nameExtensionW '/Singular'])); mkdir([nameExtensionW '/Singular']); end;
    fid=fopen([runFolder nameExtensionW '/Singular/elasticStrainEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigElast_SIN, 'descend')'./max(eigElast_SIN));   fclose(fid);

else % No Domain decomposition
    phiElast_SIN = [] ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orthogonalisation of inelastic snp                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPos>0
display(['STEP ' num2str(step) ': ortho of the Regular and Singular snp']); step=step+1;
inelasStrainSnp_REG=zeros(size(inelasStrainSnp));
inelasStrainSnp_SIN=zeros(size(inelasStrainSnp));
% inelasStrainSnp_REG = allStrainSnp(PointersToSet1,nElast+1:end); % JLM 
inelasStrainSnp_REG(PointersToSet1,:) = allStrainSnp(PointersToSet1,nElast+1:end);
% inelasStrainSnp_SIN = allStrainSnp(PointersToSet2,nElast+1:end); % JLM 
inelasStrainSnp_SIN(PointersToSet2,:) = allStrainSnp(PointersToSet2,nElast+1:end);

if ortSplit
    display('        ---> orthogonalization of REGULAR snapshots')
    ortInelasticStrainSnp_REG = inelasStrainSnp_REG;
    for k=1:size(inelasStrainSnp_REG,2)
        for i=1:size(phiElast_REG,2)
            ortInelasticStrainSnp_REG(:,k) = ortInelasticStrainSnp_REG(:,k) - ...
                phiElast_REG(:,i)'*inelasStrainSnp_REG(:,k)*phiElast_REG(:,i);
        end
    end
    inelasStrainSnp_REG=ortInelasticStrainSnp_REG; clear ortInelasticStrainSnp_REG;
    
    if Dom_Descomp % Domain decomposition
        display('        ---> orthogonalization of SINGULAR bands snapshots')
        ortInelasticStrainSnp_SIN = inelasStrainSnp_SIN;
        for k=1:size(inelasStrainSnp_SIN,2)
            for i=1:size(phiElast_SIN,2)
                ortInelasticStrainSnp_SIN(:,k) = ortInelasticStrainSnp_SIN(:,k) - ...
                    phiElast_SIN(:,i)'*inelasStrainSnp_SIN(:,k)*phiElast_SIN(:,i);
            end
        end
        inelasStrainSnp_SIN=ortInelasticStrainSnp_SIN;     clear ortInelasticStrainSnp_SIN;

        %check_compatibility(inelasStrainSnp_SIN, 2e-3, 'orthogonal inelastic snapshot in Singular');
        %check_compatibility(inelasStrainSnp_REG, 2e-3, 'orthgonal inelastic snapshot in Regular');
    else % No Domain decomposition
        inelasStrainSnp_SIN = [] ;
    end
end

if normaliz
    display('        ---> normalization of each')    
    if Dom_Descomp % Domain decomposition
        for i=1:size(inelasStrainSnp_REG,2)
            display(['Percents: ' num2str(i/size(inelasStrainSnp_REG,2)*100.0)])
            inelasStrainSnp_REG(:,i)=inelasStrainSnp_REG(:,i)./sqrt(inelasStrainSnp_REG(:,i)'*W*inelasStrainSnp_REG(:,i));
            inelasStrainSnp_SIN(:,i)=inelasStrainSnp_SIN(:,i)./sqrt(inelasStrainSnp_SIN(:,i)'*W*inelasStrainSnp_SIN    (:,i));
        end
        if(~isdir([nameExtensionW '/Regular'])); mkdir([nameExtensionW '/Regular']); end;
    else % No Domain decomposition
        for i=1:size(inelasStrainSnp_REG,2)
            display(['Percents: ' num2str(i/size(inelasStrainSnp_REG,2)*100.0)])
            inelasStrainSnp_REG(:,i)=inelasStrainSnp_REG(:,i)./sqrt(inelasStrainSnp_REG(:,i)'*W*inelasStrainSnp_REG(:,i));        
        end
        if(~isdir([nameExtensionW '/Singular'])); mkdir([nameExtensionW '/Singular']); end;
    end
end

display('        ---> plotting')
clear inelasStrainSnp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC FULCTUATION SVD - Regular                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': SVD of the regular inelastic snp']); step=step+1;
plot_cov_matrix(inelasStrainSnp_REG, W, 6);
Y=inelasStrainSnp_REG;
display(['        ---> rank of inelastic Regular snapshots: ' num2str(rank(Y'*W*Y,1e-12))])

[V,eigInelast_REG,~] = svd(Y,0);

tmp = real(diag(eigInelast_REG)); eigInelast_REG=tmp(find(tmp>1e-12*tmp(1))); % JLM para nPos=0 hay que comentar esta linea
phiInelast_REG=zeros(size(Y,1),size(eigInelast_REG,1)); %JLM size(eigInelast_REG,2) por size(eigInelast_REG,1)

for i=1:size(eigInelast_REG,1) ;
    phiInelast_REG(:,i)=V(:,i); 
end %Y*V(:,i)/sqrt(eigInelast_REG(i)); end

clear Y; clear V;
plot_cov_matrix(phiInelast_REG, W, 7);
figure(8); semilogy(sort(eigInelast_REG, 'descend')'./max(eigInelast_REG), '-*'); title('Singular values (E_{SNP})');
fid=fopen([runFolder nameExtensionW '/Regular/inelasStrainEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigInelast_REG, 'descend')'./max(eigInelast_REG));   fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC FULCTUATION SVD - Singular                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Dom_Descomp % Domain decomposition
    display(['STEP ' num2str(step) ': SVD of the singular inelastic fluctuation snp matrix']); step=step+1;
    plot_cov_matrix(inelasStrainSnp_SIN, W, 9);
    Y=inelasStrainSnp_SIN;
    display(['        ---> rank of inelastic Singular snapshots: ' num2str(rank(Y'*W*Y))])

    [V,eigInelast_SIN,~] = svd(Y,0);

    tmp = real(diag(eigInelast_SIN)); eigInelast_SIN=tmp(find(tmp>1e-12*tmp(1))); % JLM para nPos=0 hay que comentar esta linea
    phiInelast_SIN=zeros(size(Y,1),size(eigInelast_SIN,1));%JLM size(eigInelast_SIN,2) por size(eigInelast_SIN,1)


    for i=1:size(eigInelast_SIN,1) ;
        phiInelast_SIN(:,i)=V(:,i); 
    end %/sqrt(eigInelast_SIN(i)); end

    clear Y; clear V;
    plot_cov_matrix(phiInelast_SIN, W, 10);
    figure(11); semilogy(sort(eigInelast_SIN, 'descend')'./max(eigInelast_SIN), '-*'); title('Singular values (E_{SNP})');
    fid=fopen([runFolder nameExtensionW '/Singular/inelasStrainEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigInelast_SIN, 'descend')'./max(eigInelast_SIN));   fclose(fid);
else % No Domain decomposition
    phiInelast_SIN = []
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUSION OF ALL MODES                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPos==0  %JLM
    PHI_EPS_REG=[phiElast_REG ]; % phi_bar
    PHI_EPS_DIS=[phiElast_SIN ]; % phi_barbar
else
    PHI_EPS_REG=[phiElast_REG phiInelast_REG]; % phi_bar
    PHI_EPS_DIS=[phiElast_SIN phiInelast_SIN]; % phi_barbar
end

save([runFolder nameExtensionW '/allStrainModes_' nameExtensionW '.mat'], 'PHI_EPS_REG', 'PHI_EPS_DIS')

display(['---> Number of Regular modes : ' num2str(size(PHI_EPS_REG,2))])
display(['---> Number of Singular modes: ' num2str(size(PHI_EPS_DIS,2))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END                                                 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('/CIMNE/Codes/toolsSVD/Step03_StrainSVD/')

display('------------------------------------------------------------')
display('-                         THE END                          -')
display('------------------------------------------------------------')
%toc; clear all
%exit
