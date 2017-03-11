function [snp, nFile] = makesnapshots( files )

display('------------')
nFile=size(files,1);
display(['  -> There is ', num2str(nFile) ' trajectories'])

%% First file in order to set the variables
f=files(1).name;
fid=fopen(f, 'r');
display([' f = ' f])

% size of the file
str = fileread( files(1).name );
%ix  = strfind( str, char(10) );
%mx  = max( ix );
tbsz = size( str );
if(tbsz(2)<100)
    sz=100;
else
    sz=tbsz(2);
end
% brut reading
txt=textscan(fid, '%s', 'delimiter', '\n', 'bufsize', sz); fclose(fid);
% Determine number of time step
nTS=size(txt{1},1);
display(['  -> nTS  = ', num2str(nTS)])

firstLine=textscan(txt{1}{1},'%f','delimiter',' ');
sDof=size(firstLine{1},1);
display(['  -> sDof = ', num2str(sDof)])

snp=zeros(sDof,nTS*nFile);
parfor k = 1:nTS
    firstLine=textscan(txt{1}{k},'%f','delimiter',' ');
    snp(:,k)=firstLine{1};
end




%% Other Files
for i=2:nFile
    f=files(i).name;
    fid=fopen(f, 'r');
    display([' f = ' f])
    % brut reading
    txt=textscan(fid, '%s', 'delimiter', '\n', 'bufsize', sz);  fclose(fid); 
    for k = 1:nTS
        firstLine=textscan(txt{1}{k},'%f','delimiter',' ', 'bufsize', sz);
        snp(:,(i-1)*nTS+k)=firstLine{1};
    end
end

end

