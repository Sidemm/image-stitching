   
    
   
    N = size(keyPt, 1);
    patch_size = 16;
    grid_size = 4; 
    pixelsPerCell= patch_size / grid_size;  
    num_bins = 8;   
    descriptors = zeros(N, grid_size*grid_size*num_bins);
    
 
   
    
    %   ComputeGradient() will take the pyramid images as input
    %   and output the gradient magnitudes and angles in two cell arrays
    %   of size equal to number of pyramid images. 
    
 
        grad_mag = cell(length(pyramid),1);
        grad_theta = cell(length(pyramid),1);
        
        [grad_mag,grad_theta] = ComputeGradient(pyramid);
    
	
	
  
    % Iterate over all keypoints
    for i = 1 : N



        % Use the function Normalize_Orientation() to normalize the  
        % gradient directions relative to the dominant gradient direction. 
  
        scale = round(keyPtScale(i)); 
        magnitudes = grad_mag{scale};
        thetas = grad_theta{scale};
        
        [patch_mag,patch_theta] = Extract_Patch(keyPt(i,:),patch_size,magnitudes,thetas);
        
        if( isempty(patch_mag))
            continue;
        end
        
        patch_theta = Normalize_Orientation(patch_mag, patch_theta);
        

        patch_mag = patch_mag .* fspecial('gaussian', patch_size, patch_size / 2);
        
        feature = ComputeSIFTDescriptor...
            (patch_mag,patch_theta,grid_size,pixelsPerCell,num_bins);
  
       % Add the feature vector we just computed to our matrix of SIFT descriptors.
        descriptors(i, :) = feature;
    end
    
    % Normalize the descriptors.
    descriptors = NormalizeDescriptors(descriptors);
end


function [grad_mag,grad_theta] = ComputeGradient(pyramid)

    
    grad_theta = cell(length(pyramid),1);
    grad_mag = cell(length(pyramid),1);
  
    for scale = 1:length(pyramid)
       
        currentImage = pyramid{scale}; 
        grad_mag{scale} = zeros(size(currentImage));
        grad_theta{scale} = zeros(size(currentImage));

        

    img_dx = filter2([-1 0 1], currentImage); %gradient in x direction(looks at pixels x+1 and x-1)
    img_dy = filter2([-1;0;1], currentImage); %gradient in y direction(looks at pixels y+1 and y-1
    
    grad_theta{scale} = atan2(img_dy,img_dx) - 2*pi/4; %atan2 returns radian and 90° decreased in radians 
    grad_mag{scale} = sqrt(img_dx.^2 + img_dy.^2);
    

        % atan2 gives angles from -pi to pi. To make the histogram code
        % easier, we'll change that to 0 to 2*pi.
        grad_theta{scale} = mod(grad_theta{scale}, 2*pi);
        
    end
    
end
        
        
function norm_angles = Normalize_Orientation(gradient_magnitudes, gradient_angles)
% Computes the dominant gradient direction for the region around a keypoint
% given the gradient magnitudes and gradient angles of the pixels in the 
% region surrounding the keypoint.


 
    
     num_bins = 36;     
     

    [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles);
    max_value = max(histogram);
    max_index = find (histogram == max_value); %finds index according to calculated max histogram value
    direction = angles(max_index);
    norm_angles = gradient_angles - direction; %normalizes the gradient orientation 

 

   % This line will re-map norm_theta into the range 0 to 2*pi
   norm_angles = mod(norm_angles, 2*pi);
        
end


function [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles)
% Compute a gradient histogram using gradient magnitudes and directions.
% Each point is assigned to one of num_bins depending on its gradient
% direction; the gradient magnitude of that point is added to its bin.

    angle_step = 2 * pi / num_bins;
    angles = 0 : angle_step : (2*pi-angle_step);
    histogram = zeros(1, num_bins);
    
    
    for i=1:size(gradient_magnitudes,1) %rows
        for j=1:length(gradient_magnitudes) %columns
            top = ceil(gradient_angles(i,j) / angle_step);
            down = floor(gradient_angles(i,j) / angle_step +1);
            % finds the angle range and chooses closest orientation
            a = min((abs(top - gradient_angles(i,j) / angle_step)),(abs(down - gradient_angles(i,j) / angle_step)));
            if (a == (abs(top-gradient_angles(i,j) / angle_step)))
                bin = top; 
            else 
                bin = down;
            end
            %solves matlab index problem with 0 
            if bin == 0
                bin = 1;
            end
            %adds the magnitudes with the same orientation
            histogram(bin) = histogram(bin) + gradient_magnitudes(i,j);
            
        end
    end

end

function descriptor = ComputeSIFTDescriptor...
            (patch_mag,patch_theta,grid_size,pixelsPerCell,num_bins)

                                                                     %
%         Compute the gradient histograms and concatenate them in the          %
%  feature variable to form a size 1x128 SIFT descriptor for this keypoint.    %

    feature = [];
    for i=1:grid_size:length(patch_mag) %proceeds 4 by 4 in rows on the 16 x 16 patch_mag 
        for j=1:grid_size:length(patch_mag) %proceeds 4 by 4 in columns on the 16 x 16 patch_mag 
            gradient_magnitudes = patch_mag(i:i+pixelsPerCell-1,j:j+pixelsPerCell-1);
            gradient_angles = patch_theta(i:i+pixelsPerCell-1,j:j+pixelsPerCell-1);
            %calculates histogram and angles for each 4 x 4 grid 
            [histogram, angles] = ComputeGradientHistogram(num_bins, gradient_magnitudes, gradient_angles); 
            feature = [feature,histogram];
        end
    end
    
    descriptor = feature;


end

    
function [patch_mag,patch_theta] = Extract_Patch(keyPt,patch_size,magnitudes,thetas)
% Finds the window of pixels that contributes to the descriptor for the
        % current keypoint.
        xAtScale = keyPt(1);%center of the DoG keypoint in the pyramid{2} image
        yAtScale = keyPt(2);
        x_lo = round(xAtScale - patch_size / 2);
        x_hi = x_lo+patch_size-1;
        y_lo = round(yAtScale - patch_size / 2);
        y_hi = y_lo+patch_size-1;

        % These are the gradient magnitude and angle images from the 
        % correct scale level. You computed these above.
       
        patch_mag = [];
        patch_theta = [];
        try    
            % Extract the patch from that window around the keypoint
            patch_mag = magnitudes(y_lo:y_hi,x_lo:x_hi);
            patch_theta = thetas(y_lo:y_hi,x_lo:x_hi);
        catch err
            % If any keypoint is too close to the boundary of the image
            % then we just skip it.
            
        end
end


function descriptors = NormalizeDescriptors(descriptors)
% Normalizes SIFT descriptors so they're unit vectors. 

    % normalize all descriptors so they become unit vectors
    lengths = sqrt(sum(descriptors.^2, 2));
    nonZeroIndices = find(lengths);
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

    % suppress large entries
    descriptors(descriptors > 0.2) = 0.2;

    % finally, renormalize to unit length
    lengths = sqrt(sum(descriptors.^2, 2));
    lengths(lengths == 0) = 1;
    descriptors = descriptors ./ repmat(lengths, [1 size(descriptors,2)]);

end