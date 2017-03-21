clc; clear all
display('------------------------------------------------------------')
display('-                       START (Step 02)                    -')
display('------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%% USERS DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/javiermro/Projects/SVD_tools/methods/matlab')
nE     = 6; % number of elastic snapshots (el menor de todos los modos)
nPos   = 83; % number of Inelastic snapshots (diferencia entre la cant. total y en el nE mayor)
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

%try
%    matlabpool close
%catch
%    n_workers = size(snpFolder,1);
%    matlabpool('open', n_workers)
%end  

for i=1:size(snpFolder,1)
   display(['** Trajectory : ' num2str(i) ' *******'])
   FuncCreateStrainSnapShotsHardeningFromMatFile(snpFolder(i,:), snpFile, nE, nPos);
end

