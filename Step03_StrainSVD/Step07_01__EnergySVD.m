clc; clear all; close all; tic;
display('------------------------------------------------------------')
display('-                      START (Step07)                      -')
display('------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%% USERS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/javiermro/Projects/SVD_tools/methods/matlab')
StrainModesFolder='/home/javiermro/Projects/SVD_tools/Offline/Modes/'; 
offlineFolder='/home/javiermro/Projects/SVD_tools/Offline/'; cd(offlineFolder);
ortSplit=1; % orth Singular and Regular snapshots
normaliz=0; % normalize snapshots before SVD
conStres=0;
ntens = 5 ; % stress/strain component for gauss point, 4 = PD; 5 = GD
tol_nE = 1e-5 ; % tolerancia para determinar el nro de modos elasticos
% tol_nI = 1e-5 ; % tolerancia para determinar el nro de modos inelasticos
nE_max = 6; % number of elastic snapshots (Maximo) JLM
nPos = 0; % number of Inelastic snapshots (DEBE COINCIDIR CON STEP 02)
nameEnergySnap = ['Phi_ela'] ; %Phi_ela  Phi_vol  Phi_dev  Phi_pla  Phi_tot % name of Snapshot energy in the folder
nameExtensionW = ['RVE_Hole01_nPos' num2str(nPos)]; % Directory of Modes (debe coincidir con Step03_01
load('/home/javiermro/Projects/Examples/RVE_Hole01/DomainPointers.mat'); % Pointer for Domain Decomposition 

nModes = 6; %number of modes or trajectories
snpFolder0 = '/home/javiermro/Projects/Examples/RVE_Hole01/Modo';
snpFile='SNAPSHOTS_RVE_Hole01.mat';
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
% Pointers corrections for energy domain
PointersToSet1 = PointersToSet1(1:ntens:size(PointersToSet1,1)) ;
PointersToSet2 = PointersToSet2(1:ntens:size(PointersToSet2,1)) ;

if(~isdir(nameExtensionW)); mkdir(nameExtensionW); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELASTIC FULCTUATION MATRIX                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step=1;
display(['STEP ' num2str(step) ': Creation of the elastic fluctuation snapshots matrix']); step=step+1;
if(size(snpFolder,1)==1)
    cd(snpFolder);
    load([nameEnergySnap '_elastic.mat']);
else
    nSnpFolder=size(snpFolder,1);
    cd(snpFolder(1,:));
    load([nameEnergySnap '_elastic.mat']);
    tmp=elasticEnergySnp;
    for ifold=2:nSnpFolder
        cd(snpFolder(ifold,:));
        load([nameEnergySnap '_elastic.mat']);
        tmp=[tmp elasticEnergySnp];
    end
    elasticEnergySnp = tmp; clear tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC SNP                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Creation of the inelastic fluctuation snapshots matrix']); step=step+1;
if(size(snpFolder,1)==1)
    cd(snpFolder);
    load([nameEnergySnap '_inelastic.mat']);
else
    nSnpFolder=size(snpFolder,1);
    cd(snpFolder(1,:));
    load([nameEnergySnap '_inelastic.mat']);
    tmp=inelasEnergySnp;
    for ifold=2:nSnpFolder
        cd(snpFolder(ifold,:));
        load([nameEnergySnap '_inelastic.mat']);
        tmp=[tmp inelasEnergySnp];
    end
    inelasEnergySnp = tmp; clear tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fusion inelastic/elastic SNP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Fusion of elastic and inelastic snapshots']); step=step+1;
nElast=size(elasticEnergySnp,2);
nInelast=size(inelasEnergySnp,2);
allEnergySnp=[elasticEnergySnp inelasEnergySnp];

cd(runFolder)
if(~isdir([nameExtensionW '/all'])); mkdir([nameExtensionW '/all']); end;
%plot_modes(get_disp_from_vec(allStress_sSnp), runFolder, [nameExtensionW '/all/elasticDispFromStress_sSNP_' nameExtensionW])
%plot_modes_mat(STRAIN_MODES_UPDATING(allStress_sSnp, 'FEToGid')', runFolder, [nameExtensionW '/all/elasticStress_sSNP_' nameExtensionW])
plot_cov_matrix(allEnergySnp, W, 1)
saveas(gca, [nameExtensionW '/all/elasticEnergyCOV_' nameExtensionW '.eps'],'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Splitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': Splitting of all snapshots in Regular/singular']); step=step+1;

% For hardening
% *************
% load('/home/mcaicedo/toolsemmanuel/Step07_StressSVD/src/pointersForInelasticDomain&AggeregatesHardening_ForEnergy.mat');

display('        ---> splitting')
allEnergySnp_REG=zeros(size(allEnergySnp)); allEnergySnp_SIN=zeros(size(allEnergySnp));
% allEnergySnp_REG=allEnergySnp(pointersToMatrixElem_ForEnergy,:);
% allEnergySnp_SIN=allEnergySnp(pointersToCBElem_ForEnergy,:);
allEnergySnp_REG = allEnergySnp(PointersToSet1,:);% JLM 
% allEnergySnp_REG(PointersToSet1,:) =allEnergySnp(PointersToSet1,:);
allEnergySnp_SIN = allEnergySnp(PointersToSet2,:);% JLM 
% allEnergySnp_SIN(PointersToSet2,:) =allEnergySnp(PointersToSet2,:);


% allEnergySnp_REG(PointersToSet1(1:ntens:size(PointersToSet1,1)),:) = allEnergySnp(PointersToSet1(1:ntens:size(PointersToSet1,1)),:);
% allEnergySnp_SIN(PointersToSet2(1:ntens:size(PointersToSet2,1)),:) = allEnergySnp(PointersToSet2(1:ntens:size(PointersToSet2,1)),:);

if conStres
ne=(size(allEnergySnp_SIN,1)/16);
load('/home/mcaicedo/toolsemmanuel/Step07_Stress_sSVD/src/COHESIVE_BANDS_NORMALS.mat')

for i=1:ne
    n=COHESIVE_BANDS_NORMALS(i,2:3);
    t=[-n(2) n(1)];
    theta=atan2(n(2),n(1));
    PASS=[cos(theta) -sin(theta); sin(theta) cos(theta)];
    for j=1:size(allEnergySnp_SIN,2)
        display('        ---> constant stress_s')
        pt=16*(i-1);
        sigXX=(allEnergySnp_SIN(pt+1,j)+allEnergySnp_SIN(pt+5,j)+allEnergySnp_SIN(pt+9,j) +allEnergySnp_SIN(pt+13,j))./4.0;
        sigYY=(allEnergySnp_SIN(pt+2,j)+allEnergySnp_SIN(pt+6,j)+allEnergySnp_SIN(pt+10,j)+allEnergySnp_SIN(pt+14,j))./4.0;
        sigZZ=(allEnergySnp_SIN(pt+3,j)+allEnergySnp_SIN(pt+7,j)+allEnergySnp_SIN(pt+11,j)+allEnergySnp_SIN(pt+15,j))./4.0;
        sigXY=(allEnergySnp_SIN(pt+4,j)+allEnergySnp_SIN(pt+8,j)+allEnergySnp_SIN(pt+12,j)+allEnergySnp_SIN(pt+16,j))./4.0;
        
        SIG_XY=[sigXX sigXY; sigXY sigYY];
        SIG_NT=transpose(PASS)*SIG_XY*PASS;
        SIG_NT(2,2)=0;
        SIG_XY=PASS*SIG_NT*transpose(PASS);

        allEnergySnp_SIN(pt+1,j)=SIG_XY(1,1);
        allEnergySnp_SIN(pt+2,j)=SIG_XY(2,2);
        allEnergySnp_SIN(pt+3,j)=sigZZ;
        allEnergySnp_SIN(pt+4,j)=SIG_XY(1,2);
        allEnergySnp_SIN(pt+5,j)=SIG_XY(1,1);
        allEnergySnp_SIN(pt+6,j)=SIG_XY(2,2);
        allEnergySnp_SIN(pt+7,j)=sigZZ;
        allEnergySnp_SIN(pt+8,j)=SIG_XY(1,2);
        allEnergySnp_SIN(pt+9,j)=SIG_XY(1,1);
        allEnergySnp_SIN(pt+10,j)=SIG_XY(2,2);
        allEnergySnp_SIN(pt+11,j)=sigZZ;
        allEnergySnp_SIN(pt+12,j)=SIG_XY(1,2);
        allEnergySnp_SIN(pt+13,j)=SIG_XY(1,1);
        allEnergySnp_SIN(pt+14,j)=SIG_XY(2,2);
        allEnergySnp_SIN(pt+15,j)=sigZZ;
        allEnergySnp_SIN(pt+16,j)=SIG_XY(1,2);

    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ELASTIC FULCTUATION SVD                             %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': SVD of the Regular elastic snp']); step=step+1;
Y=allEnergySnp_REG(:,1:nElast);
plot_cov_matrix(allEnergySnp_REG(:,1:nElast), W, 2)
display(['        ---> rank of elastic snapshots: ' num2str(rank(Y'*W*Y))])

%[V, eigElast_REG] = eigs(Y'*W*Y, size(Y,2)-2);
[V,eigElast_REG,~]=svd(Y,0);


tmp = diag(eigElast_REG); eigElast_REG=tmp(find(tmp>1e-14));
phiElast_REG=zeros(size(Y,1),size(eigElast_REG,1)); % JLM cambie size(eigElast_REG,2)) por size(eigElast_REG,1))

for i=1:size(eigElast_REG,1) ;
    phiElast_REG(:,i)=V(:,i); 
end % Y*V(:,i)/sqrt(eigElast_REG(i)); end

%check_compatibility(phiElast_REG, 2e-3, 'elastic mode');
plot_cov_matrix(phiElast_REG, W, 3);
figure(4); semilogy(eigElast_REG, '-*'); title('Singular values (E_{SNP})');
% display('        ---> keep only the first six elastic modes')
% eigElast_REG = eigElast_REG(1:min(6,size(eigElast_REG,1))); 
% phiElast_REG = phiElast_REG(:,1:min(6,size(eigElast_REG,1)));

nE = min(nE_max,size(find(eigElast_REG>tol_nE),1)) ; %JLM
display(['        ---> keep only the first ' num2str(nE) ' elastic modes']); %JLM
eigElast_REG = eigElast_REG(1:min(nE,size(eigElast_REG,1))); %JLM
phiElast_REG = phiElast_REG(:,1:min(nE,size(eigElast_REG,1))); %JLM
if(~isdir([nameExtensionW '/Regular'])); mkdir([nameExtensionW '/Regular']); end;

fid=fopen([runFolder nameExtensionW '/Regular/elasticEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', eigElast_REG');   fclose(fid);

display(['STEP ' num2str(step) ': SVD of the Singular elastic snp']); step=step+1;
Y=allEnergySnp_SIN(:,1:nElast);
plot_cov_matrix(allEnergySnp_SIN(:,1:nElast), W, 5)
display(['        ---> rank of elastic snapshots: ' num2str(rank(Y'*W*Y))])

%[V, eigElast_SIN] = eigs(Y'*W*Y, size(Y,2)-2);
[V,eigElast_SIN,~]=svd(Y,0);

tmp = diag(eigElast_SIN); eigElast_SIN=tmp(find(tmp>1e-14));
phiElast_SIN=zeros(size(Y,1),size(eigElast_SIN,1)); % JLM cambie size(eigElast_SIN,2)) por size(eigElast_SIN,1))

for i=1:size(eigElast_SIN,1) ;
    phiElast_SIN(:,i)=V(:,i);
end

%check_compatibility(phiElast_SIN, 2e-3, 'elastic mode');
plot_cov_matrix(phiElast_SIN, W, 6);
figure(7); semilogy(eigElast_SIN, '-*'); title('Singular values (E_{SNP})');
% display('        ---> keep only the first six elastic modes')
% eigElast_SIN = eigElast_SIN(1:min(6,size(eigElast_SIN,1))); 
% phiElast_SIN = phiElast_SIN(:,1:min(6,size(eigElast_SIN,1)));

nE = min(nE_max,size(find(eigElast_SIN>tol_nE),1)) ; %JLM
display(['        ---> keep only the first ' num2str(nE) ' elastic modes']); %JLM
eigElast_SIN = eigElast_SIN(1:min(nE,size(eigElast_SIN,1))); %JLM
phiElast_SIN = phiElast_SIN(:,1:min(nE,size(eigElast_SIN,1))); %JLM

if(~isdir([nameExtensionW '/Singular'])); mkdir([nameExtensionW '/Singular']); end;
%%%plot_modes(get_disp_from_vec(phiElast_SIN), runFolder, [nameExtensionW '/Singular/elasticDispFromStress_sModes_' nameExtensionW])
%%%plot_modes_mat(STRAIN_MODES_UPDATING(phiElast_SIN, 'FEToGid'), runFolder, [nameExtensionW '/Singular/elasticStress_sModes_' nameExtensionW])
%%fid=fopen([runFolder nameExtensionW '/Singular/elasticStress_sEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigElast_SIN, 'descend')'./max(eigElast_SIN));   fclose(fid);

fid=fopen([runFolder nameExtensionW '/Singular/elasticEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', eigElast_SIN');   fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orthogonalisation of inelastic snp                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': ortho of the Regular and Singular snp']); step=step+1;

% inelasEnergySnp_REG=allEnergySnp(PointersToSet1(1:ntens:size(PointersToSet1,1)),nElast+1:end);
% inelasEnergySnp_SIN=allEnergySnp(PointersToSet1(1:ntens:size(PointersToSet2,1)),nElast+1:end);

inelasEnergySnp_REG=zeros(size(inelasEnergySnp));
inelasEnergySnp_SIN=zeros(size(inelasEnergySnp));
inelasEnergySnp_REG = allEnergySnp(PointersToSet1,nElast+1:end); % JLM
% inelasEnergySnp_REG(PointersToSet1,:) = allEnergySnp(PointersToSet1,nElast+1:end); %
inelasEnergySnp_SIN = allEnergySnp(PointersToSet2,nElast+1:end); % JLM
% inelasEnergySnp_SIN(PointersToSet2,:) = allEnergySnp(PointersToSet2,nElast+1:end); % 

if ortSplit
    display('        ---> orthogonalization of REGULAR snapshots')
    ortInelasticEnergySnp_REG = inelasEnergySnp_REG;
    for k=1:size(inelasEnergySnp_REG,2)
        for i=1:size(phiElast_REG,2)
            ortInelasticEnergySnp_REG(:,k) = ortInelasticEnergySnp_REG(:,k) - ...
                phiElast_REG(:,i)'*inelasEnergySnp_REG(:,k)*phiElast_REG(:,i);
        end
    end
    inelasEnergySnp_REG=ortInelasticEnergySnp_REG; clear ortInelasticEnergySnp_REG;

    display('        ---> orthogonalization of SINGULAR snapshots')
    ortInelasticEnergySnp_SIN = inelasEnergySnp_SIN;
    for k=1:size(inelasEnergySnp_SIN,2)
        for i=1:size(phiElast_SIN,2)
            ortInelasticEnergySnp_SIN(:,k) = ortInelasticEnergySnp_SIN(:,k) - phiElast_SIN(:,i)'*inelasEnergySnp_SIN(:,k)*phiElast_SIN(:,i);
        end
    end
    inelasEnergySnp_SIN=ortInelasticEnergySnp_SIN;     clear ortInelasticEnergySnp_SIN;
    
    %check_compatibility(inelasStress_sSnp_SIN, 2e-3, 'orthogonal inelastic snapshot in Singular');
    %check_compatibility(inelasStress_sSnp_REG, 2e-3, 'orthgonal inelastic snapshot in Regular');
end

if normaliz
    display('        ---> normalization of each')
    for i=1:size(inelasEnergySnp_REG,2)
        display(['Percents: ' num2str(i/size(inelasEnergySnp_REG,2)*100.0)])
        inelasEnergySnp_REG(:,i)=inelasEnergySnp_REG(:,i)./sqrt(inelasEnergySnp_REG(:,i)'*W*inelasEnergySnp_REG(:,i));
        inelasEnergySnp_SIN(:,i)=inelasEnergySnp_SIN(:,i)./sqrt(inelasEnergySnp_SIN(:,i)'*W*inelasEnergySnp_SIN    (:,i));
    end
end

if(~isdir([nameExtensionW '/Regular'])); mkdir([nameExtensionW '/Regular']); end; if(~isdir([nameExtensionW '/Singular'])); mkdir([nameExtensionW '/Singular']); end;
display('        ---> plotting')
%%plot_modes(get_disp_from_vec(inelasStress_sSnp_REG), runFolder, [nameExtensionW '/Regular/inelastDispFromStress_sSNP_REG_' nameExtensionW])
%%plot_modes(get_disp_from_vec(inelasStress_sSnp_SIN), runFolder, [nameExtensionW '/Singular/inelastDispFromStress_sSNP_SIN_' nameExtensionW])
%%plot_modes_mat(STRAIN_MODES_UPDATING(inelasStress_sSnp_REG, 'FEToGid'), runFolder, [nameExtensionW '/Regular/inelastStress_sSNP_REG_' nameExtensionW])
%%plot_modes_mat(STRAIN_MODES_UPDATING(inelasStress_sSnp_SIN,     'FEToGid'), runFolder, [nameExtensionW '/Singular/inelastStress_sSNP_SIN_' nameExtensionW])
clear inelasEnergySnp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC FULCTUATION SVD - Regular                  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPos>0
display(['STEP ' num2str(step) ': SVD of the regular inelastic snp']); step=step+1;
plot_cov_matrix(inelasEnergySnp_REG, W, 6);
Y=inelasEnergySnp_REG;
display(['        ---> rank of inelastic Regular snapshots: ' num2str(rank(Y'*W*Y,1e-12))])

[V,eigInelast_REG,~]=svd(Y,0);
%[V, eigInelast_REG] = eigs(Y'*W*Y, size(Y,2)-2);

tmp = diag(eigInelast_REG); eigInelast_REG=tmp(find(tmp>1e-12*tmp(1)));
phiInelast_REG=zeros(size(Y,1),size(eigInelast_REG,1)); %JLM size(eigInelast_REG,2) por size(eigInelast_REG,1)

for i=1:size(eigInelast_REG,1) ;
    phiInelast_REG(:,i)=V(:,i); 
end

%tmp = real(diag(eigInelast_REG)); eigInelast_REG=tmp(find(tmp>1e-12*tmp(1)));
%phiInelast_REG=zeros(size(Y,1),size(eigInelast_REG,2));
%for i=1:size(eigInelast_REG,1) phiInelast_REG(:,i)=Y*V(:,i)/sqrt(eigInelast_REG(i)); end
clear Y; clear V;
plot_cov_matrix(phiInelast_REG, W, 7);
figure(8); semilogy(sort(eigInelast_REG, 'descend')'./max(eigInelast_REG), '-*'); title('Singular values (E_{SNP})');
%plot_modes(get_disp_from_vec(phiInelast_REG), runFolder, [nameExtensionW '/Regular/inelastDispFromStress_sModes_' nameExtensionW])
%%plot_modes_mat(STRAIN_MODES_UPDATING(phiInelast_REG, 'FEToGid'), runFolder, [nameExtensionW '/Regular/inelastStress_sModes_' nameExtensionW])
%%fid=fopen([runFolder nameExtensionW '/Regular/inelasEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigInelast_REG, 'descend')'./max(eigInelast_REG));   fclose(fid);

fid=fopen([runFolder nameExtensionW '/Regular/inelasEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', eigInelast_REG');   fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INELASTIC FULCTUATION SVD - Singular                      %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['STEP ' num2str(step) ': SVD of the singular inelastic fluctuation snp matrix']); step=step+1;
plot_cov_matrix(inelasEnergySnp_SIN, W, 9);
Y=inelasEnergySnp_SIN;
display(['        ---> rank of inelastic Singular snapshots: ' num2str(rank(Y'*W*Y))])

[V,eigInelast_SIN,~]=svd(Y,0);
%[V, eigInelast_SIN] = eigs(Y'*W*Y, size(Y,2)-2);

tmp = diag(eigInelast_SIN); eigInelast_SIN=tmp(find(tmp>1e-12*tmp(1)));
phiInelast_SIN=zeros(size(Y,1),size(eigInelast_SIN,1)); %JLM size(eigInelast_SIN,2) por size(eigInelast_SIN,1)

for i=1:size(eigInelast_SIN,1) ;
    phiInelast_SIN(:,i)=V(:,i); 
end

%tmp = real(diag(eigInelast_SIN)); eigInelast_SIN=tmp(find(tmp>1e-12*tmp(1)));
%phiInelast_SIN=zeros(size(Y,1),size(eigInelast_SIN,2));
%for i=1:size(eigInelast_SIN,1) phiInelast_SIN(:,i)=Y*V(:,i)/sqrt(eigInelast_SIN(i)); end
clear Y; clear V;
plot_cov_matrix(phiInelast_SIN, W, 10);
figure(11); semilogy(sort(eigInelast_SIN, 'descend')'./max(eigInelast_SIN), '-*'); title('Singular values (E_{SNP})');
%plot_modes(get_disp_from_vec(phiInelast_SIN), runFolder, [nameExtensionW '/Singular/inelastDispFromStress_sModes_' nameExtensionW])
%plot_modes_mat(STRAIN_MODES_UPDATING(phiInelast_SIN, 'FEToGid'), runFolder, [nameExtensionW '/Singular/inelastStress_sModes_' nameExtensionW])
%%fid=fopen([runFolder nameExtensionW '/Singular/inelasEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', sort(eigInelast_SIN, 'descend')'./max(eigInelast_SIN));   fclose(fid);

fid=fopen([runFolder nameExtensionW '/Singular/inelasEnergyEIG_' nameExtensionW '.dat'],'w');   fprintf(fid, '%.25f\n', eigInelast_SIN');   fclose(fid);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUSION OF ALL STRAIN AND ENERGY MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nPos==0  %JLM
    PHI_ENER_REG=[phiElast_REG ]; % phi_bar
    PHI_ENER_DIS=[phiElast_SIN ]; % phi_barbar
    
    Sing_Val_ELAS_REG=eigElast_REG;
    Sing_Val_INELAS_REG=[];
    Sing_Val_ELAS_DIS=eigElast_SIN;
    Sing_Val_INELAS_DIS=[];
else
    PHI_ENER_REG=[phiElast_REG phiInelast_REG]; % phi_bar
    PHI_ENER_DIS=[phiElast_SIN phiInelast_SIN]; % phi_barbar
    
    Sing_Val_ELAS_REG=eigElast_REG;
    Sing_Val_INELAS_REG=eigInelast_REG;
    Sing_Val_ELAS_DIS=eigElast_SIN;
    Sing_Val_INELAS_DIS=eigInelast_SIN;
end

load([StrainModesFolder nameExtensionW '/allStrainModes_' nameExtensionW '.mat']);

% save([runFolder nameExtensionW '/allEnergyModes4ROMI_' nameExtensionW '.mat'], 'PHI_ENER_REG', 'Sing_Val_ELAS_REG', 'Sing_Val_INELAS_REG', 'PHI_ENER_DIS', 'Sing_Val_ELAS_DIS', 'Sing_Val_INELAS_DIS')
save([runFolder nameExtensionW '/allStrainEnergyModes_' nameExtensionW '_' nameEnergySnap '.mat'], ...
    'PHI_EPS_REG', 'PHI_EPS_DIS', 'PHI_ENER_REG', 'Sing_Val_ELAS_REG', 'Sing_Val_INELAS_REG', 'PHI_ENER_DIS', ...
    'Sing_Val_ELAS_DIS', 'Sing_Val_INELAS_DIS')

display(['---> Number of Regular modes : ' num2str(size(PHI_ENER_REG,2))])
display(['---> Number of Singular modes: ' num2str(size(PHI_ENER_DIS,2))])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cd('/CIMNE/Codes/toolsSVD/Step03_StrainSVD/')

display('------------------------------------------------------------')
display('-                         THE END                          -')
display('------------------------------------------------------------')
toc; clear all
%exit
