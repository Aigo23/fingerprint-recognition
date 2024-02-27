
function [newim, binim, mask, orientim, reliability] =  testfin(im)
    
    %if nargin == 0
	%im = imread('finger.png');
   % end
    
    % Identify ridge-like regions and normalise image
    blksze = 16; thresh = 0.2;
    [normim, mask] = ridgesegment(im, blksze, thresh);
    %figure(11),subplot(2,3,1),imshow(mask),title("mask");%add
    %show(normim,1);
    
    % Determine ridge orientations
    [orientim, reliability] = ridgeorient(normim, 1, 5, 5);
    %plotridgeorient(orientim, 20, im, 2)
    %show(reliability,6)
    
    % Determine ridge frequency values across the image
    blksze = 32; 
    [freq, medfreq] = ridgefreq(normim, mask, orientim, blksze, 5, 5, 15);
    %show(freq,3) 
    
    % Actually I find the median frequency value used across the whole
    % fingerprint gives a more satisfactory result...
    

    mask = freq > 0;
    %figure(11),subplot(2,3,2),imshow(mask),title("mask_freqed");
    mask = bwareaopen(mask,30*16*16,4);
    mask = ~bwareaopen(~mask,9*16*16,4);
    %figure(11),subplot(2,3,3),imshow(mask),title("mask_areaed");

    disk = strel('disk', blksze/2);
    mask = imopen(mask,disk);
    mask = imclose(mask,disk);
    mask = bwareaopen(mask,30*16*16,4);
    mask = ~bwareaopen(~mask,9*16*16,4);
    
    %figure(11),subplot(2,3,4),imshow(mask),title("mask_closed");

    freq = medfreq.*mask;
    
    % Now apply filters to enhance the ridge pattern
    newim = ridgefilter(normim, orientim, freq, 0.5, 0.5, 0);
    %show(newim,4);
    
    % Binarise, ridge/valley threshold is 0
    binim = newim > 0;
    show(binim,5);
    % Display binary image for where the mask values are one and where
    % the orientation reliability is greater than 0.5
   %show(binim.*mask.*(reliability>0.5), 7)