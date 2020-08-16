
% Code for model in Figure 6 
% 2-cell model for population failure

% Uses function error_ellipse: https://www.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse

% Written by JZ, organized by KR
% August 2020

%%

% Set the firing rate response of each cell to the 2 stimuli. Each row is a neuron and each column in a stimulus response.
FRs = [15 20; 20 23]; % Positive signal correlations. % For negative signal correlations, try [15 20; 23 20] 

% Set the noise correlation value between the 2 cells
r = 0.8; % Any value between -0.9 and 0.9
corrmat = [1 r; r 1];

% Make a vector of possible firing rates for plotting distributions
FRlist = linspace(0,30,31);

% For each stimulus, compute covariance, plot ellipses, and get distribution of single neuron firing rates
figure();
for ii = 1:length(FRs(1,:)) 
    mean = FRs(:,ii);      
    covmat = diag(sqrt(mean))*corrmat*diag(sqrt(mean)); % Covariance matrix with Poisson noise and given noise correlation
    subplot(221);error_ellipse(covmat,mean); hold on; % Plot the full error ellipse
    e = error_ellipse(diag(diag(covmat)),mean); % Plot the uncorrelated ellipse
    e.LineStyle = '--';

    prob1(ii,:) = (1/sqrt(2*pi*mean(1,:)))*exp(-((FRlist-mean(1,:)).^2)/(2*mean(1,:))); % Distribution of firing for cell 1 (Poisson with mean = var)
    prob2(ii,:) = (1/sqrt(2*pi*mean(2,:)))*exp(-((FRlist-mean(2,:)).^2)/(2*mean(2,:))); % Distribution of firing for cell 2 (Poisson with mean = var)

end
xlabel('Firing Rate Cell 1 (Hz)')
ylabel('Firing Rate Cell 2 (Hz)')
xlim([5 30])
hold off

% Compute d prime squared (following Averbeck and Lee J Neurophys 2006)
du = FRs(:,1) - FRs(:,2); % Difference in firing rates for the two stimuli
covmat1 = diag(sqrt(FRs(:,1)))*corrmat*diag(sqrt(FRs(:,1))); % Covariance for stimulus 1
covmat2 = diag(sqrt(FRs(:,2)))*corrmat*diag(sqrt(FRs(:,2))); % Covariance for stimulus 2
Q = 0.5*(covmat1 + covmat2); % Average covariance
Qd = diag(diag(Q)); % Covariance for independent model (no noise correlations)

dprime_single1 = du(1)'*inv(Q(1,1))*du(1); % Single neuron discriminability for cell 1
dprime_single2 = du(2)'*inv(Q(2,2))*du(2); % Single neuron discriminability for cell 2

d2_diag = (du'*inv(Qd)*du)^2/(du'*inv(Qd)*Q*inv(Qd)*du); % Discriminability of the real data using a decoder that ignores noise correlations
d2 = du'*inv(Q)*du; % Discriminability of the real data using a decoder that takes the whole noise covariance matrix Q into account


% Plot decoding performance
subplot(222);
hold on
x = [1, 11, 21];
bar(x,[max(dprime_single1,dprime_single2),d2_diag,d2])
hline(max(dprime_single1,dprime_single2))
set(gca,'XTick',x)
set(gca,'XTickLabel',{'single cell','ind.','cpl.'});
ylabel('discriminability')


% Plot response distributions for the best single cell for the single best cell
cell1best = dprime_single1 > dprime_single2; % boolean: 1 if cell 1 is best, 0 if cell 2 is best
subplot(223);
hold on
if(cell1best)
    plot(FRlist,prob1);
else
    plot(FRlist,prob2);
end
xlabel('Firing Rate (Hz)')
ylabel('Probability density')
hold off





%%






