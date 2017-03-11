clc; clear all
display('------------------------------------------------------------')
display('-                          START                           -')
display('------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%% USERS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('C:/CIMNE/Codes/toolsSVD/methods/matlab')
% load('C:/CIMNE/Codes/MultiScale_SantaFe/Examples/REV_Laminate_01/DomainPointers.mat'); % Pointer for Domain Decomposition 
load('C:/CIMNE/Codes/MultiScale_SantaFe/Examples/RVE_Hole_05/DomainPointers_Hole05.mat'); % Pointer for Domain Decomposition 
% load('C:/CIMNE/Codes/MultiScale_SantaFe/Examples/TEST13/DomainPointers_TEST13.mat'); % Pointer for Domain Decomposition 
femSolution='HF_Snapshot';% Mode_01';  % LOAD FEM SOLUTION (en esa carpeta se copia los Snapshot de deformacion del HF correspondientes al modo que se esta comparando)
femRoot='C:/CIMNE/Codes/toolsSVD/Offline/Modes/'; 
phiRoot='C:/CIMNE/Codes/toolsSVD/Offline/Modes/'; % LOAD BASIS
phiFolders=dir([phiRoot 'HOLE05_nPos74_GDJ2_PSI_EP*']); 
ntens = 5 ; % stress/strain component for gauss point, 4 = PD; 5 = GD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
femFolder=[femRoot femSolution];

% Pointers corrections for energy domain
PointersToSet1 = PointersToSet1(1:ntens:size(PointersToSet1,1)) ;
PointersToSet2 = PointersToSet2(1:ntens:size(PointersToSet2,1)) ;

cd(femFolder);
load('SnapEnergy.mat'); 
fem_reg=SnapEnergy(PointersToSet1,:);
fem_sin=SnapEnergy(PointersToSet2,:); 
clear SnapEnergy;

%     cd(phiFolders);
%     load(['inelasPosEnergySnp_' num2str(nPos)])
%     inelasEnergySnp=[inelasPosEnergySnp];
% 
% 
% EnergySnp=[elasticEnergySnp inelasPosEnergySnp];
% 
% fem_dis=EnergySnp(pointersToCBElem,:); fem_reg=EnergySnp(pointersToMatrixElem,:);
% clear EnergySnp;
% 
% %% LOAD BASIS
% phiRoot='/home/mcaicedo/toolsemmanuel/Offline/Modes/';


for n=1:size(phiFolders,1)    
    modeName=phiFolders(n).name;
    phiFolder=[phiRoot modeName];
    %cd(phiFolder); load(['allEnergyModes4ROMI_' modeName '.mat']); phi_dis=PHI_ENER_DIS(pointersToCBElem,:); phi_reg=PHI_ENER_REG(pointersToMatrixElem,:);
    cd(phiFolder); load(['allStrainEnergyModes4ROMI_' modeName '.mat']); phi_dis=PHI_ENER_DIS; phi_reg=PHI_ENER_REG;
    clear PHI_ENER_DIS PHI_ENER_REG; 
    
    %% ERROR REG
    display(['*** ' modeName ' regular ***'])
    normType='fro'; phiAll=phi_reg; fem=fem_reg; FileName='reg';
    nModes=size(phiAll,2); err=[];
    for m = 1:nModes
        phi = phiAll(:,1:m);
        TrialPhi = phi*(phi'*fem);
        TrialOrt = fem - TrialPhi ;
        nTrialOrt = norm(TrialOrt, normType);
        nfem      = norm(fem,      normType);
        err = [err nTrialOrt/nfem*100];
        disp([' --> err(' num2str(m) '/' num2str(nModes) ') = ', num2str(err(m)) '%']);
        clear phi TrialOrt TrialPhi nTrialOrt nfem orthogonalError;
        if(m>1)
            if(err(m)<0.1 || err(m)>err(m-1))
                break
            end
        end
    end
    
    fid=fopen([phiFolder '/' femSolution '_' modeName '_AprioriError_' FileName '.dat'], 'w');
    index=1:size(err,2);
    fprintf(fid, '%1.0f %1.5f\n', [index; err]);
    fclose(fid);
    
    %% ERROR DIS
    display(['*** ' modeName ' singular ***'])
    normType='fro'; phiAll=phi_dis; fem=fem_sin; FileName='dis';
    nModes=size(phiAll,2); err=[];
    for m = 1:nModes
        phi = phiAll(:,1:m);
        TrialPhi = phi*(phi'*fem);
        TrialOrt = fem - TrialPhi ;
        nTrialOrt = norm(TrialOrt, normType);
        nfem      = norm(fem,      normType);
        err = [err nTrialOrt/nfem*100];
        disp([' --> err(' num2str(m) '/' num2str(nModes) ') = ', num2str(err(m)) '%']);
        clear phi TrialOrt TrialPhi nTrialOrt nfem orthogonalError;
        if(m>1)
        if(err(m)<0.1 || err(m)>err(m-1))
            break
        end
        end
    end
    
    fid=fopen([phiFolder '/' femSolution '_' modeName '_AprioriError_' FileName '.dat'], 'w');
    index=1:size(err,2);
    fprintf(fid, '%1.0f %1.5f\n', [index; err]);
    fclose(fid);
end


%exit
