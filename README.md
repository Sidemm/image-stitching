# image-stitching
Image stitching using SIFT descriptor and RANSAC.

To determine the gradient angles and magnitude, the function ComputeGradient() is used.
ComputeGradient() takes the pyramid images as inputand output the gradient magnitudes and angles 
in two cell arrays of size equal to number of pyramid images. 

Patch of magnitudes and directions (16x16) around the  keypoint are extracted using the function Extract_Patch().

Normalize_Orientation() is used to normalize the gradient directions relative to the dominant gradient direction. 

The patch extracted is subdivided intogrid_size x grid_size cells, each of which is sizepixelsPerCell x pixelsPerCell. 
        
A gradient histogram is computed for each cell, and the histograms are concataneted into a single feature vector of length 128.
        
The patch is traversed row by row, starting in the top left, in order to match the given solution. 

Lastly Ransac algorithm is used to find a transformation from points from one image to points from another image which 
I found as good matches earlier. 
      
