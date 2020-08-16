
% Estimate pixel covariance and mean values in natural images by sampling random patches for use with RFModel.m
% Conditioned on when a center pixel is on or off


% Load the natural images from Olshausen and Field 1996 
load('IMAGES_RAW.mat')

sz = sqrt(nPix); % Size of image patches
BUFF = 10; % Buffer zone to not get too close to the edge
batch_size = 10000; % Number patches per image 
imsize = 512;

% Go through the images one by one
for i = 1:10
    img = squeeze(IMAGESr(:,:,i));
    % Make a binarized image, with roughly half ones and half zeros
    med = median(median(img));
    img = (img>med);
    
    ctr_on = 1;
    ctr_off = 1;
    % Draw a bunch of random patches to estimate the covariance 
    for j=1:batch_size
		r=BUFF+ceil((imsize-sz-2*BUFF)*rand);
		c=BUFF+ceil((imsize-sz-2*BUFF)*rand);
        patch = reshape(img(r:r+sz-1,c:c+sz-1),sz^2,1);
        
        if(patch((nPix-1)/2 + 1) == 1)
            patches_on(ctr_on,:) = patch;
            ctr_on = ctr_on + 1;
        end
        
        if(patch((nPix-1)/2 + 1) == 0)
            patches_off(ctr_off,:) = patch;
            ctr_off = ctr_off + 1;
        end
    end
    
    covs_on(i,:,:) = cov(patches_on);
    covs_off(i,:,:) = cov(patches_off);
    
    meanpatch_on(i,:) = mean(patches_on)';
    meanpatch_off(i,:) = mean(patches_off)';

end

pixpixcov_on  = squeeze(mean(covs_on));
pixpixcov_off  = squeeze(mean(covs_off));

mean_on = squeeze(mean(meanpatch_on))';
mean_off = squeeze(mean(meanpatch_off))';


        
    