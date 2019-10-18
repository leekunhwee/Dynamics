function FRF = FRFGeneration(Dynamic,Simulation,flag)

if flag == 1 %% Residues type
    %% Define modal parameters
    wnX = Dynamic.f.x*2*pi;
    wnY = Dynamic.f.y*2*pi;
    zetaX = Dynamic.zeta.x;
    zetaY = Dynamic.zeta.y;    
    resiX = Dynamic.residues.x;
    resiY = Dynamic.residues.y;
    %% Simulation parameter
    [f,w] = FrequencySequence(Simulation);

    %% Define cells
    FRFiX={1,length(resiX)};
    FRFiY={1,length(resiY)};
    %% Initialization
    FRFX   = zeros(length(w),1);
    FRFY   = zeros(length(w),1);
    alphaX = zeros(1,length(resiX));
    alphaY = zeros(1,length(resiY));
    betaX  = zeros(1,length(resiX));
    betaY  = zeros(1,length(resiY));
    wdX    = zeros(1,length(resiX));
    wdY    = zeros(1,length(resiY));
    sigmaX = real(resiX);
    sigmaY = real(resiY);
    nuX    = imag(resiX);
    nuY    = imag(resiY);
    %% Modes superposition
    for cnt = 1:length(resiX)
        wdX(cnt)   = wnX(cnt)*sqrt(1-zetaX(cnt)^2);
        alphaX(cnt)= 2*(zetaX(cnt)*wnX(cnt)*sigmaX(cnt)-wdX(cnt)*nuX(cnt));
        betaX(cnt) = 2*sigmaX(cnt);
        FRFiX{cnt} = (alphaX(cnt).*ones(length(w),1)+1j*betaX(cnt).*w)./(wnX(cnt)^2-w.^2 +1j* (2*zetaX(cnt)*wnX(cnt).*w));
        FRFX= FRFX + FRFiX{cnt};
    end
    for cnt = 1:length(resiY)
        wdY(cnt)   = wnY(cnt)*sqrt(1-zetaY(cnt)^2);
        alphaY(cnt)= 2*(zetaY(cnt)*wnY(cnt)*sigmaY(cnt)-wdY(cnt)*nuY(cnt));
        betaY(cnt) = 2*sigmaY(cnt);
        FRFiY{cnt} = (alphaY(cnt).*ones(length(w),1)+1j*betaY(cnt).*w)./(wnY(cnt)^2-w.^2 +1j* (2*zetaY(cnt)*wnY(cnt).*w));
        FRFY= FRFY + FRFiY{cnt};
    end
    % Dividing to real and imag
    HX(:,1) = real(FRFX);
    HX(:,2) = imag(FRFX);
    HY(:,1) = real(FRFY);
    HY(:,2) = imag(FRFY);
    % Save
    FRF = struct;
    FRF.X = FRFX;% Feeding direction
    FRF.Y = FRFY;% Normal direction
    FRF.w = w;
    %% Print to .frf for CutPro
    Matlab2CutPro(f,HX,HY)
    
elseif flag == 2 %% Stiffness type
    %% Define modal parameters
    wnX = Dynamic.f.x*2*pi;
    wnY = Dynamic.f.y*2*pi;
    zetaX = Dynamic.zeta.x;
    zetaY = Dynamic.zeta.y;    
    kX = Dynamic.k.x;
    kY = Dynamic.k.y;
    %% Simulation parameter
    f_start = Simulation.f_start;
    f_end   = Simulation.f_end;
    df      = Simulation.df;
    f = (f_start:df:f_end)';   % frequency, Hz
    w = f*2*pi;                % frequency, rad
    %% Define cells
    FRFiX={1,length(wnX)};
    FRFiY={1,length(wnY)};
    %% Initialization
    FRFX = zeros(length(w),1);
    FRFY = zeros(length(w),1);
    %% Modes superposition
    for cnt = 1:length(wnX)
        FRFiX{cnt} = (wnX(cnt)^2/kX(cnt))./(wnX(cnt)^2 - w.^2 + 1i*2*zetaX(cnt)*wnX(cnt).*w);
        FRFX = FRFX + FRFiX{cnt};
    end
    for cnt = 1:length(wnY)
        FRFiY{cnt} = (wnY(cnt)^2/kY(cnt))./(wnY(cnt)^2 - w.^2 + 1i*2*zetaY(cnt)*wnY(cnt).*w);
        FRFY = FRFY + FRFiY{cnt};
    end
    % Dividing to real and imag
    HX(:,1) = real(FRFX);
    HX(:,2) = imag(FRFX);
    HY(:,1) = real(FRFY);
    HY(:,2) = imag(FRFY);
    % Save
    FRF = struct;
    FRF.X = FRFX;% Feeding direction
    FRF.Y = FRFY;% Normal direction
    FRF.w = w;
    %% Print to .frf for CutPro
    Matlab2CutPro(f,HX,HY)
end