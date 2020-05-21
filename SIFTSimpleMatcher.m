function match = SIFTSimpleMatcher(descriptor1, descriptor2, thresh)

    if ~exist('thresh', 'var'),
        thresh = 0.7;
    end

    match = [];
                             
    for i=1:length(descriptor1)
        for j=1:length(descriptor2)
            euclidean(j) = norm(descriptor1(i,:) - descriptor2(j,:));
        end
            distance = sort(euclidean);
            if distance(1)<distance(2)*thresh % compares the best match and the second best match
            match = [match; [i, find(euclidean == distance(1))]] %concatenates matched features
            end
    end
    

end
