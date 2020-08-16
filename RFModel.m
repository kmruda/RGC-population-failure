
% Code for model in Figure 7
% Simpel RF model to see how correlation structure impacts the consequences of ignoring noise correlations
% For both white noise stimuli and natural images

% Uses functions:
% - NaturalImageCov: script in same repository
% - L2_distance: http://web.mit.edu/cocosci/isomap/code/L2_distance.m
% - hline: https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline

% Written by JZ, organized by KR
% August 2020

%% Define parameters and define the pixel grid

% Add the path for the code & data for natural scenes stimuli
addpath(genpath('/Downloads/RGC-pop-failure/'))

% Define parameters
rmaxlist = [0 0.25 0.5 0.75];   % Max noise correlation, for nearby neurons
lambdalist = [100:100:500]; % Scale over which correlations decay (exponentially, um)
pixSize = 250; % Size of pixel (um)
RFSize = 250;  % Receptive field diameter (um). We'll model RFs as circles.
overlap = 0.1; % Amount of linear overlap between any two neighboring cells
maxFR = 30; % Max firing rate (Hz). This FR is achieved for a cell with RF fully covered by a bright pixel.
nPix = 25; % Number of pixels. Should be the square of an odd number so there is a central pizel to be decoded.
sz = sqrt(nPix); % Edgelength of the pixel grid.

uncorr = 0; % Boolean flag for using a binary white noise stimulus
natscenes = 1; % Boolean flag for using a natural scene stimulus

% Define the pixel grid
edgeLength = pixSize*sqrt(nPix);
centsX = linspace(0,edgeLength,sqrt(nPix)+1);
centsX = centsX(1:end-1) + pixSize/2;
centsY = centsX;
pixCX = centsX; % Save the centers for later
pixCY = centsY;
[pixX pixY] = meshgrid(centsX,centsY);

%% Define the neuron grid and compute overlap with pixels

% Define the neuron grid
central_center = edgeLength/2; % Center location in x,y of center neuron. Assume it is centered under the central pixel
n_extra = ceil(central_center/(RFSize*(1-overlap)));% Number of extra neurons on each side to get past edge of pixel grid. Neuron grid is 2*n_extra + 1 x 2*n_extra + 1
nNeur = (2*n_extra + 1)^2; % Total number of neurons. Note, some are extending past pixel grid, but they will be removed later
orig_nNeur = nNeur;
% Plop neurons down that tile the space: define their center positions here. Use triangles for a hexagonal grid.
xDist = (RFSize*(1-overlap)*sqrt(3)/2); % Distance between columns of cells
c1x = fliplr(central_center - xDist*(0:n_extra));
c2x = central_center + xDist*(0:n_extra);
centsX = [c1x c2x(2:end)];
c1y = fliplr(central_center - RFSize*(1-overlap)*(0:n_extra));
c2y = central_center + RFSize*(1-overlap)*(0:n_extra);
centsY = [c1y c2y(2:end)];
[neurX neurY] = meshgrid(centsX,centsY);
cntr_col = ceil((2*n_extra+1)/2);
if mod(cntr_col,2) == 1 
    cls = 2:2:(2*n_extra + 1); % Offset even columns by (RFSize*(1-overlap))/2
elseif mod(cntr_col,2) == 0
    cls = 1:2:(2*n_extra + 1); % Offset odd columns by (RFSize*(1-overlap))/2
end

% Compute distances between RF centers, which partially defines the strength of noise correlations
neurY(:,cls) = neurY(:,cls) + (RFSize*(1-overlap))/2;
neurXvec = reshape(neurX,1,nNeur);
neurYvec = reshape(neurY,1,nNeur);
distanceMat = sqrt((L2_distance(neurXvec,neurXvec,1)).^2 + (L2_distance(neurYvec,neurYvec,1)).^2);

% Compute overlaps between RFs and pixels
[posX posY] = meshgrid(1:edgeLength,1:edgeLength); % Use this to write out RFs and pixels on a fine grid. Compute overlaps numerically
% Get RF size of center neuron (one that is fully visible on the grid); used for overlapfrac
nn = (nNeur+1)/2;
RF = zeros(edgeLength,edgeLength); 
RF((posX-neurXvec(nn)).^2 + (posY-neurYvec(nn)).^2  <= (RFSize/2)^2) = 1; % Fill in with ones for the part in the RF in the circle
RF_area = sum(sum(RF));

% Create  a nNeur x nPix matrix of RF-pixel overlaps.
neurPixOverlap = zeros(nNeur,nPix);
all_rfs = zeros(nNeur,edgeLength,edgeLength);
for nn = 1:nNeur
    for pp = 1:nPix
        % Make the neuron's RF, on a edgeLength x edgeLength grid. 
        RF = zeros(edgeLength,edgeLength);
        RF((posX-neurXvec(nn)).^2 + (posY-neurYvec(nn)).^2  <= (RFSize/2)^2) = 1; % Fill in with ones for the part in the RF in the circle

        % Now put the pixel on the grid
        pix = zeros(edgeLength,edgeLength);
        pix( (abs(posX - pixX(pp)) < pixSize/2) & (abs(posY - pixY(pp)) < pixSize/2) ) = 1;

        % Multiply elementwise to get overlap
        prod = pix.*RF;

        % Compute fraction of RF covered by pixel
        overlapfrac = sum(sum(prod))/RF_area;
        neurPixOverlap(nn,pp) = overlapfrac;
    end
    all_rfs(nn,:,:) = RF;
end

% Get rid of neurons that don't overlap with the pixel grid
sumz = sum(neurPixOverlap,2);
sumz(sumz==0)=nan;
nan_inds = isnan(sumz);
for nn = 1:nNeur
    if nan_inds(nn)==1
        neurPixOverlap(nn,:) = nan(1,nPix);
    end
end

% Update distanceMat variable: omit the nan rows and columns 
i = 1;
for rr = 1:nNeur
    if nan_inds(rr) == 1
        % Do nothing
    else
        temp0_distanceMat(i,:) = distanceMat(rr,:);
        i = i+1;
    end
end
i = 1;
for cc = 1:nNeur
    if nan_inds(cc) == 1
        % Do nothing
    else
        temp1_distanceMat(:,i) = temp0_distanceMat(:,cc);
        i = i+1;
    end
end
distanceMat = temp1_distanceMat;
figure;imagesc(distanceMat)
    
% Update neurPixOverlap
k = 1;
for nn = 1:nNeur
    if nan_inds(nn)==1
        % Do nothing
    else
        temp0_neurPixOVerlap(k,:) = neurPixOverlap(nn,:);
        k = k+1;
    end
end
figure;imagesc(temp0_neurPixOVerlap);
neurPixOverlap = temp0_neurPixOVerlap;

% Update nNeur
nNeur = k-1;

% Plot RFs and pixels
running_sum = zeros(edgeLength,edgeLength);
for nn = 1:orig_nNeur
    if ~isnan(nan_inds(nn))
        running_sum = running_sum + squeeze(all_rfs(nn,:,:));
    end
end
figure();imagesc(running_sum);hold on;colorbar;
vline(pixCX + pixSize/2); vline(pixCX - pixSize/2);
hline(pixCY + pixSize/2); hline(pixCY-pixSize/2);
hold off

%% Compute information calculations for different rmax and lambda values

for ll = 1:length(lambdalist)
    for rr = 1:length(rmaxlist)
        lambda = lambdalist(ll);
        rmax = rmaxlist(rr);
        neurNoiseCorrs = rmax*exp(-distanceMat/lambda).*[ones(nNeur) - eye(nNeur)] + eye(nNeur); % Pairwise noise correlations

        % Compute pixel-pixel covariances
        if(uncorr ==1) % Binary white noise stimulus
            % Always 0 (except for diagonal entries) for white noise
            pixpixcov_on = 0.25*eye(nPix); %when center pixel is on or off
            pixpixcov_off = 0.25*eye(nPix); 

            % Mean pixel pattern when the central pixel is on or off
            mean_on = 0.5*ones(sz,sz);
            mean_on((sz-1)/2 + 1,(sz-1)/2 + 1) = 1;
            mean_on = reshape(mean_on, nPix,1);

            mean_off = 0.5*ones(sz,sz);
            mean_off((sz-1)/2 + 1,(sz-1)/2 + 1) = 0;
            mean_off = reshape(mean_off, nPix,1);
        end
        if(natscenes == 1) % Natural scene stimulus
            NaturalImageCov % Get the natural scenes covariance matrix and mean image patches from binarized image patches
        end

        % Compute neuron-pixel covariances, conditioned on when the central pixel is on or off
        neurPixCov_on = neurPixOverlap*maxFR*pixpixcov_on; % All neurons to all pixels when central pixel is on
        neurPixCov_off = neurPixOverlap*maxFR*pixpixcov_off; % All neurons to all pixels when central pixel is off

        % Chop out entries for the central pixel
        neurPixCov_on_noCenter = neurPixCov_on;
        neurPixCov_on_noCenter(:,(nPix-1)/2 + 1) = [];
        neurPixCov_off_noCenter = neurPixCov_off;
        neurPixCov_off_noCenter(:,(nPix-1)/2 + 1) = [];
        % Same for the neuron-pixel overlap-- make one excluding center pixel
        neurPixOverlap_noCenter = neurPixOverlap;
        neurPixOverlap_noCenter(:,(nPix-1)/2 + 1) = [];
        % and the pixel-pixel covariance-- make one excluding center pixel
        pixpixcov_on_noCenter = pixpixcov_on;
        pixpixcov_on_noCenter(:,(nPix-1)/2 + 1) = [];
        pixpixcov_on_noCenter((nPix-1)/2 + 1,:) = [];
        pixpixcov_off_noCenter = pixpixcov_off;
        pixpixcov_off_noCenter(:,(nPix-1)/2 + 1) = [];
        pixpixcov_off_noCenter((nPix-1)/2 + 1,:) = [];

        % Compute information values from Averbeck and Lee J Neurophys 2006
        % The d2 values are similar to Fisher information, but for discrete stimuli

        % du is the vector of differences in mean neural firing rates (or pixel values 
        % for the outer, non-decoded pixels used in estimation) between central pixel on vs off. 
        du_pix = mean_on - mean_off;
        du_pix((nPix-1)/2+1) = [];
        % Assume linear neuron, with firing rate r = RF*image + noise and define pixels in binary image to be +1 (on) or 0 (off).
        du = [maxFR*neurPixOverlap*(mean_on - mean_off) ; du_pix];
        % The first nNeur values are the neurons' differences in mean response for
        % central pixel being on vs off. The next nPix - 1 values are differences in mean pixel values
        % of outer, non-decoded pixels for central pixel being on vs off. 

        % Q is the mean covariance matrix (averaged over having central pixel on or off)
        % First calculate it for central pixel on
        Qon = zeros(nNeur + nPix-1, nNeur + nPix-1); 
        % Fill in Poisson noise for the neuron-neuron covariance, with the given
        % noise correlation structure
        % Use law of Total covariance: 
        % cov(n1,n2) = E[cov(n1,n2|outer_pix)] + cov(E[n1|outer_pix],E[n2|outer_pix])
        % First term is the usual noise covariance, second term is the
        % signal covariance driven by variations in outer pixels. With
        % regards to central pixel, that's effective noise.
        Qon(1:nNeur,1:nNeur) = maxFR*sqrt(diag(neurPixOverlap*mean_on))*neurNoiseCorrs*sqrt(diag(neurPixOverlap*mean_on)) ...
                               +maxFR^2*neurPixOverlap_noCenter*pixpixcov_on_noCenter*neurPixOverlap_noCenter'; 
        % Add pixel-pixel covariance
        Qon(nNeur+1:end,nNeur+1:end) = pixpixcov_on_noCenter;
        % neuron-pixel covariance as computed above, for non-central pixels
        Qon(1:nNeur,nNeur+1:end) = neurPixCov_on_noCenter; 
        Qon(nNeur+1:end,1:nNeur) = neurPixCov_on_noCenter';

        % Repeat for when the central pixel is off
        Qoff = zeros(nNeur + nPix-1, nNeur + nPix-1); 
        Qoff(1:nNeur,1:nNeur) = maxFR*sqrt(diag(neurPixOverlap*mean_off))*neurNoiseCorrs*sqrt(diag(neurPixOverlap*mean_off)) ...
                                +maxFR^2*neurPixOverlap_noCenter*pixpixcov_off_noCenter*neurPixOverlap_noCenter';
        Qoff(nNeur+1:end,nNeur+1:end) = pixpixcov_off_noCenter;
        Qoff(1:nNeur,nNeur+1:end) = neurPixCov_off_noCenter; 
        Qoff(nNeur+1:end,1:nNeur) = neurPixCov_off_noCenter'; 

        % Average over ON and OFF central pixel to get Q
        Q = 0.5*Qon + 0.5*Qoff;

        % Check that covariance is sensible
        if min(eig(Q))<0
            disp('ERROR: non-PSD covariance')
            rmax
            lambda
            break
        end

        % Make the diagonal version of Q to ignore correlations
        % Remove only neuron-neuron noise correlations, keeping the neuron-neuron covariance due to changing outer pixels.
        Qd = Q;
        Qd(1:nNeur,1:nNeur) = diag(diag(Q(1:nNeur,1:nNeur))) ...
            +(0.5*+maxFR^2*neurPixOverlap_noCenter*pixpixcov_on_noCenter*neurPixOverlap_noCenter' ...
            +0.5*maxFR^2*neurPixOverlap_noCenter*pixpixcov_off_noCenter*neurPixOverlap_noCenter').*(1-eye(nNeur,nNeur));

        % d2_diag is the info that would be obtained on the real data, using a
        % decoder that ignores correlations
        d2_diag(rr,ll) = (du'*inv(Qd)*du)^2/(du'*inv(Qd)*Q*inv(Qd)*du);

        % d2 is the info that is obtained on the real data using a decoder that 
        % accounts for correlations
        d2(rr,ll) = du'*inv(Q)*du;

    end
end

% Percent difference between ignoring correlations and accounting for them
pd = 100*(d2-d2_diag)./d2_diag;

% For comparison, compute single-neuron info, including the outer, non-decoded pixels
central_neuron = ceil(nNeur/2);
Q_single = zeros(nPix,nPix);
Q_single(1,1) = Q(central_neuron,central_neuron);
Q_single(1,2:end) = Q(central_neuron,nNeur+1:end);
Q_single(2:end,1) = Q(nNeur+1:end,central_neuron);
Q_single(2:end,2:end) = Q(nNeur+1:end,nNeur+1:end);
dprime_single = [du(central_neuron) ; du_pix]'*inv(Q_single)*[du(central_neuron) ; du_pix];

%% Plot results

% Discriminability
colors = [0.9 0.9 0.9;0.6 0.6 0.6;0.3 0.3 0.3;0 0 0];
figure;
for rr = 1:length(rmaxlist)
    plot(lambdalist,d2(rr,:),'o-','Color',colors(rr,:));hold on
    p = plot(lambdalist,d2_diag(rr,:),'o--','Color',colors(rr,:));
    p.HandleVisibility = 'off';
end
hline(dprime_single,'r-')
xlabel('Correlation spatial scale (um)');ylabel('Discriminability')
legend(num2str(rmaxlist'))

% Percent difference between the independent and coupled decoders
figure;
for rr = 1:length(rmaxlist)
    plot(lambdalist,pd(rr,:),'o-','Color',colors(rr,:));hold on
end
xlabel('Correlation spatial scale (um)');ylabel('Percent difference')
legend(num2str(rmaxlist'))

%%



