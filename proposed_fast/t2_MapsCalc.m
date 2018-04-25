
function [t2Maps,S0Maps] = t2_MapsCalc(vol,TEs)

[Nx, Ny, Ns, Nz] = size(vol);
numTEs=size(TEs,1);
sets = floor(Ns/(numTEs));
disp('T2 map generation...');
disp(['--> imaging matrix = ' num2str(Nx) ' x ' num2str(Ny) ' x ' num2str(Nz)]);
disp(['--> TE values = ' num2str(TEs') ' ms']);

tic;
toc0 = toc;

% calculate T1rho maps
t2Maps = zeros(Nx,Ny,sets,Nz);  % one map for each set of TE values
S0Maps = zeros(Nx,Ny,sets,Nz); 
NxNy = Nx*Ny;
NxNyTEs = NxNy*numTEs;   %actually not needed... multiplied by 0 ahead..

stepInc=(0:99)'*NxNy;
te=[TEs];te=te(:);
A=[ones(length(te),1),te];

for nz = 1:Nz
    disp(['Analyzing slice ' num2str(nz) ' / ' num2str(Nz) '... (t = ' num2str(toc-toc0) ' sec)']);
    % work with smaller volumes for increased performance
    volBuffer = vol(:,:,:,nz);
    volMaskBuffer = squeeze(vol(:,:,1,nz));
    t2MapsBuffer = zeros(Nx,Ny,sets);
    S0Buffer=zeros(Nx,Ny,sets);
     nonZeroIndices = find(abs(volMaskBuffer(:))>0);
%     nonZeroIndices = find(volMaskBuffer(:)>0);
    for ns = 0:sets-1  % index shifted by 1 for faster processing
        steps = stepInc + NxNyTEs * ns;
        NxNyns = NxNy * ns;
        for ni = 1:size(nonZeroIndices,1)
            Svalues = volBuffer(nonZeroIndices(ni) + steps(1:numTEs));
            P=A\log(Svalues);
            %P = coeffvalues(fit(TEs,Svalues,'exp1'));
            t2MapsBuffer(nonZeroIndices(ni) + NxNyns)=-1/P(2);
            S0Buffer(nonZeroIndices(ni) + NxNyns)=exp(P(1));
        end
    end
S0Maps(:,:,:,nz)=S0Buffer;
t2Maps(:,:,:,nz) = t2MapsBuffer;
end

% clean up T2 maps
t2Maps(~isfinite(t2Maps)) = 0;  % set any NANs or INFs to zero
t2Maps(t2Maps<0) = 0;  % set any negative T1rho values to zero
t2Maps(t2Maps>(median(t2Maps(t2Maps(:)>0))*10)) = 0;  % set high noise values to zero
% 
% clean up S0 maps
S0Maps(~isfinite(S0Maps)) = 0;  % set any NANs or INFs to zero
S0Maps(S0Maps<0) = 0;  % set any negative T1rho values to zero
S0Maps(S0Maps>(median(S0Maps(S0Maps(:)>0))*10)) = 0;  % set high noise values to zero
clc;
