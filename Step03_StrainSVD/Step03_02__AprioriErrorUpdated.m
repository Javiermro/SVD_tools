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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phiFolders=dir([phiRoot 'RVE_StructHole10_nPos83*']);  

femFolder=[phiRoot femSolFolder];
cd(femFolder);

if exist( [femFolder '/' femSolFile])
    load(femSolFile,'SnapStrain', 'Snapflag') ;
    strainSnp = SnapStrain ; 
    fem_reg   = SnapStrain(PointersToSet1,:);
    fem_sin   = SnapStrain(PointersToSet2,:); 
    flagSnp   = Snapflag   ;
    clear SnapStrain Snapflag     
else
    error('binary files not detected, please check if the mat files are already created!')
end



for n=1:size(phiFolders,1)    
    modeName=phiFolders(n).name;
    phiFolder=[phiRoot modeName];
    cd(phiFolder); load(['allStrainModes_' modeName '.mat']); 
    phi_reg=PHI_EPS_REG(PointersToSet1,:);
    phi_sin=PHI_EPS_DIS(PointersToSet2,:); 
    clear PHI_EPS_DIS PHI_EPS_REG;
    
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
    
    %% ERROR SIN
    display(['*** ' modeName ' singular ***'])
    normType='fro'; phiAll=phi_sin; fem=fem_sin; FileName='dis';
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
display('------------------------------------------------------------')
display('-                         THE END                          -')
display('------------------------------------------------------------')

%exit
