clc; clear all
display('------------------------------------------------------------')
display('-                          START                           -')
display('------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%% USERS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/javiermro/Projects/SVD_tools/methods/matlab')
load('/home/javiermro/Projects/Examples/StructHole10/DomainPointers.mat'); % Pointer for Domain Decomposition 
femSolFolder='HF_Snapshot';% Mode_01';  % LOAD FEM SOLUTION (en esa carpeta se copia los Snapshot de deformacion del HF correspondientes al modo que se esta comparando)
femSolFile='SNAPSHOTS_RVE_StructHole10.mat';
phiRoot='/home/javiermro/Projects/SVD_tools/Offline/Modes/'; 
phiFolders=dir([phiRoot 'RVE_StructHole10_nPos83*']);                          
ntens = 5 ; % stress/strain component for gauss point, 4 = PD; 5 = GD
nameEnergySnap = ['Phi_ela'] ; %Phi_ela  Phi_vol  Phi_dev  Phi_pla  Phi_tot % name of Snapshot energy in the folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phiFolders=dir([phiRoot 'RVE_StructHole10_nPos83*']);  
femFolder=[phiRoot femSolFolder];
cd(femFolder);

% Pointers corrections for energy domain
PointersToSet1 = PointersToSet1(1:ntens:size(PointersToSet1,1)) ;
PointersToSet2 = PointersToSet2(1:ntens:size(PointersToSet2,1)) ;

if exist( [femFolder '/' femSolFile])
    switch nameEnergySnap
        case 'Phi_ela'
            load(femSolFile,'SnapEnergy_e') ; SnapEnergy = SnapEnergy_e ; clear SnapEnergy_e
        case 'Phi_vol' 
            load(femSolFile,'SnapEnergy_e_vol') ; SnapEnergy = SnapEnergy_e_vol ; clear SnapEnergy_e_vol
        case 'Phi_dev' 
            load(femSolFile,'SnapEnergy_e_dev') ; SnapEnergy = SnapEnergy_e_dev ; clear SnapEnergy_e_dev
        case 'Phi_pla' 
            load(femSolFile,'SnapEnergy_p') ; SnapEnergy = SnapEnergy_p ; clear SnapEnergy_p
        case 'Phi_tot' 
            load(femSolFile,'SnapEnergy_t') ; SnapEnergy = SnapEnergy_t ; clear SnapEnergy_t
        otherwise
            warning('Unexpected energy snapshot definition')
    end
    fem_reg = SnapEnergy(PointersToSet1,:);
    fem_sin = SnapEnergy(PointersToSet2,:); 
    load(femSolFile,'Snapflag') ; flagSnp = Snapflag ; clear Snapflag 
else
    error('binary files not detected, please check if the mat files are already created!')
end

for n=1:size(phiFolders,1)    
    modeName=phiFolders(n).name;
    phiFolder=[phiRoot modeName]; 
    cd(phiFolder); load(['allStrainEnergyModes_' modeName '_' nameEnergySnap '.mat']); 
    phi_dis = PHI_ENER_DIS; phi_reg = PHI_ENER_REG;
    clear PHI_ENER_DIS  PHI_ENER_REG  PHI_EPS_DIS  PHI_EPS_REG
    clear Sing_Val_ELAS_REG  Sing_Val_INELAS_REG  Sing_Val_ELAS_DIS  Sing_Val_INELAS_DIS
    
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
    
    fid=fopen([phiFolder '/' femSolFolder '_' modeName '_AprioriError_' FileName '.dat'], 'w');
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
    
    fid=fopen([phiFolder '/' femSolFolder '_' modeName '_AprioriError_' FileName '.dat'], 'w');
    index=1:size(err,2);
    fprintf(fid, '%1.0f %1.5f\n', [index; err]);
    fclose(fid);
end


%exit
