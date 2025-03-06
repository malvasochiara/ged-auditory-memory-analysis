function [O] = GED_single_subject_singletrial_outlier(S)
O = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads a single subject, computed GED on it and store the data in a .mat file
% the covariance matrix is computed on each trial independently. Outliers
% are removed both from covS and covR. The covraince matrices for GED are
% the one averaged over trials.

% INPUTS :  -S.filenameNarrow  = name of the narrowband file
%           -S.pathNarrow      = path to the narrowband file
%           -S.filenameBroad   = name of the broadband file
%           -S.pathBroad       = path to the broadband file
%           -S.destpath        = path to the folder where you store the results      
%           -S.ncomps          = numper of components you want to save
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filenameNarrow = S.filenameNarrow;
    pathNarrow = S.pathNarrow;
    filenameBroad = S.filenameBroad;
    pathBroad = S.pathBroad;
    destinationpath = S.destpath;
    ncomps = S.ncomps;
%    comps2see = S.comp2see;
    shr   = 0.01; %fixed
%%%%%  DATA LOADING  %%%%%
    dataN = load([pathNarrow '/' filenameNarrow]);
    dataN = dataN.OUT.sources_ERFs{1,1};%selecting the desired condition
    %broadband data
    dataB = load([pathBroad '/' filenameBroad]); 
    dataB = dataB.OUT.sources_ERFs{1,1};%selecting the desired condition

    sizebroad = size(dataB);
    sizenarrow = size(dataN);
%%%%% GED COMPUTATION %%%%%
% Initialize the covariance matrix with a third dimension for trials
    covS_trials = zeros(sizenarrow(1),sizenarrow(1), sizenarrow(3));
    covR_trials = zeros(sizebroad(1),sizebroad(1), sizebroad(3));
% Compute covariance matrices
    for tt = 1:sizebroad(3) %over trials
        covS_trials(:,:,tt) = cov(dataN(:,:,tt)');
        covR_trials(:,:,tt) = cov(dataB(:,:,tt)');

        % Regularisation (to bring covariance matrices to full rank)
        % narrow data
        covS_trials(:,:,tt) = covS_trials(:,:,tt)  + 1e-6*eye(size(covS_trials(:,:,tt) )); %making matrix full rank again by adding a small perturbation/noise (matrix regularisation)
        evalsR = eig(covR_trials(:,:,tt) );  % for regularization of R covariance matrix
        covR_trials(:,:,tt)  = (1-shr)*covR_trials(:,:,tt)  + shr * mean(evalsR) * eye(size(covR_trials(:,:,tt) )); % regularization
        clear evalsR
    end
    %averaging over trials
    covS_avg = mean(covS_trials(:,:,:),3);
    covR_avg = mean(covR_trials(:,:,:),3);
    
    % Removing bad trialsbased on distance
    for tt = 1:sizebroad(3)
        % computing the distance
        distS(tt) = sqrt(trace(covS_trials(:,:,tt)*covS_avg));
        distR(tt) = sqrt(trace(covR_trials(:,:,tt)*covR_avg));
    end
    
    %Normalize distance (temporary variable)
    distS_z = (distS - mean(distS)) / std(distS);
    distR_z = (distR - mean(distR)) / std(distR);
    z_thresh = 2.3;  %~.01
    
    %find indexes outliers
    toofarS = abs(distS_z) > z_thresh;
    toofarR = abs(distR_z) > z_thresh;
 
    %remove distances just for visualization, using logical indexing
    distS_z(toofarS) = NaN;
    distR_z(toofarR) = NaN;
    
    %Remove actual outlier matrices, using logical indexing
    covS_trials(:,:,toofarS) = NaN;
    covR_trials(:,:,toofarR) = NaN;
    
    badtrialS = [];
    badtrialR = [];
    countR = 0;
    countS =0;
    % finding the index of "bad" trials to save them for later
    for tt = 1:sizebroad(3) % over trials
        if toofarR(tt)==1
            countR = countR+1;
            badtrialR(countR) =tt;
        end
         
        if toofarS(tt)==1
            countS = countS+1;
            badtrialS(countS) =tt;
         end
    end
    
    
    
    %Re-compute grand-averages (whatch out the NaNs, use 'omitnan')
    ntrials = size(dataB,3) - length(badtrialR); % getting the effective number of trials, to beused later on in the time series
    covS = squeeze(mean(covS_trials , 3 , 'omitnan')); %along 3D dimension
    covR= squeeze(mean(covR_trials, 3 , 'omitnan'));
    % Eigendecomposition 
    [evecs,evals] = eig(covS ,covR);
    [evals,sidx]  = sort( diag(evals),'descend' ); % first output returns sorted evals extracted from diagonal
    evecs = evecs(:,sidx);          % sort eigenvectors
    evals = evals.*100./sum(evals); % normalize eigenvalues to percent variance explained   

% compute filter forward model and flip sign
     GEDmap_abs = zeros(ncomps,sizebroad(1));
%     GEDts_abs  = zeros(ncomps,size(broadData,2), size(dataB,3));
    GEDmap = zeros(ncomps,sizebroad(1));
    GEDts  = zeros(ncomps,sizebroad(2), ntrials);
    
    dataB(:,:,badtrialR) = []; %removing bad trials
    for compi = 1:ncomps
        
         % Absolute value
%         GEDmap_abs(compi,:) = evecs(:,compi)'*covS; % get component
%         [~,idxmax] = max(abs(GEDmap_abs(:,1)));     % find max magnitude
%         GEDmap_abs(compi,:)  = abs(GEDmap_abs(compi,:)*sign(GEDmap_abs(compi,idxmax))); % possible sign flip
%             % compute time series (projections) USING BROADBAND DATA 
%             GEDts_abs(compi,:) = abs(evecs(:,compi)'*dataB(:,:,tt));


        % NO absolute value used: abs value before plotting the brain images
        GEDmap(compi,:) = evecs(:,compi)'*covS ; % get component
        [~,idxmax] = max(abs(GEDmap(:,1)));     % find max magnitude
        GEDmap(compi,:)  = GEDmap(compi,:)*sign(GEDmap(compi,idxmax)); % possible sign flip

        for tt = 1:ntrials % over good trials, reconstructing a timeserie for each trial
            % compute time series (projections) USING BROADBAND DATA 
            GEDts(compi,:,tt) = evecs(:,compi)'*dataB(:,:,tt);
        end
    end
    
%%%%%  SAVING  %%%%%
        filename = [filenameNarrow(1:10) 'GEDresults_singlesub'];
        %savingpath = [destinationpath '/' pathNarrow(90:112) '/'];
        save([destinationpath '/' filename '.mat'], 'covS','covR','evecs','evals','GEDmap','GEDts','sizebroad','sizenarrow','badtrialS','badtrialR','-v7.3')
