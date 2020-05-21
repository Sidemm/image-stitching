function H_best = RANSACFit(p1, p2, match, seedSampleSize, maxInlierError, goodFitThresh )



    N = size(match, 1);
    if N<3
        error('not enough matches to produce a transformation matrix')
    end
    
    if ~exist('seedSetSize', 'var'),
        seedSampleSize = ceil(0.2 * N);
    end
     
    if ~exist('maxInlierError', 'var'),
        maxInlierError = 30;
    end
    
    if ~exist('goodFitThresh', 'var'),
        goodFitThresh = floor(0.7 * N);
    end



    
    %probability that all samples fail
    p_fail = 0.001;
    
    %fraction of inliers
    w = 0.7;
    
    %seed sample size
    n = seedSampleSize;
    
    %number of iterations
    maxIter=[];
    
   
       maxIter = log(p_fail)/log(1-w^n);
   
    H_best = eye(3);
    min_error = inf;
        
    for i = 1 : maxIter

        % Randomly select a seed group 
        idx = randperm(size(match, 1));
        seed_group = match(idx(1:seedSampleSize), :);

        % Select remaining as the non-seed group
        non_seed_group = match(idx(seedSampleSize+1:end), :);


        H = ComputeAffineMatrix( p1(seed_group(:, 1), :), p2(seed_group(:, 2), :) );
        

        err = ComputeError(H, p1(non_seed_group(:, 1), :), p2(non_seed_group(:, 2),:));
                                                 
        
        inliers = [];
       
        inliers = (err < maxInlierError);                                                                 
      
   
        number_of_inliers = size(inliers,1) + size(seed_group,1); 
        if( number_of_inliers > goodFitThresh )

	
         
            inliers_seed = [seed_group ; non_seed_group(inliers,:)]; % accumulates all the inliners adding seed group and 
                                                                     % the members of non_seed_group with low error(inliers)                                                                
            H= ComputeAffineMatrix(p1(inliers_seed(:,1),:),p2(inliers_seed(:,2),:));
            err = ComputeError(H,p1(inliers_seed(:,1),:),p2(inliers_seed(:,2),:));
            if sum(err)<min_error 
                H_best = H;
                min_error = sum(err); %saves the new min_error in every iteration
            end
          

        end

    end
    
    if sum(sum((H_best - eye(3)).^2)) == 0,
        disp('No RANSAC fit was found.')
    end
end


function H = ComputeAffineMatrix( Pt1, Pt2 )


    N = size(Pt1,1);
    if size(Pt1, 1) ~= size(Pt2, 1),
        error('Dimensions unmatched.');
    elseif N<3
        error('At least 3 points are required.');
    end
    
    % Convert the input points to homogeneous coordintes.
    P1 = [Pt1';ones(1,N)];
    P2 = [Pt2';ones(1,N)];

    
    
    H = [];

    H = P1'\P2'; %solves for P1'*H'=P2'
    H = H';
    

    H(3,:) = [0 0 1];
end


function dists = ComputeError(H, pt1, pt2)


   dists = zeros(size(pt1,1),1);


    
     point1 = H*[pt1,ones(length(pt1),1)]'; %adds a third dimension filling it with ones and 
                                            %transforms it to pt2's coordinate system                                           %coordinate system
     point2 = [pt2,ones(length(pt2),1)];    
     point1 = point1';
     for i = 1 : length(dists)
     dists(i) = norm(point1(i)-point2(i));  %fills dists array with euclidean distances
     end


    if size(dists,1) ~= size(pt1,1) || size(dists,2) ~= 1
        error('wrong format');
    end
end





