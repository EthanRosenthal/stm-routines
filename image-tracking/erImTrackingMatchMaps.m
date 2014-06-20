function [movingMapOut, movingTopoOut, fixedTopoOut] = erImTrackingMatchMaps(movingMap,movingTopo,referenceTopo,cropSize)
% Use this program to warp a dI/dV map and topography to match a reference
% topography.

%Inputs
%  movingMap = dI/dV map to warp (3-D matrix where last dimension corresponds to energy slice)
%  movingTopo = topography coincident to movingMap
%  referenceTopo = topography to which one wants to map movingTopo
%  cropSize = number of pixels to which one wants to crop everything.

%Outputs
%  movingMapOut = warped dI/dV map
%  movingTopoOut = warped topography coincident to movingMapOut
%  fixedTopoOut = cropped referenceTopo


% Get tforms using topographies.
[tformSimilarity,tformAffine,fixed_out] = erImTrackingTforms(movingTopo,referenceTopo,cropSize);
%Warp the dI/dV map using the tforms from above
[movingMapOut]=erImTrackingTransform(movingMap,fixed_out,tformSimilarity,tformAffine,1);

%Warp the moving topography using the tforms from above
movingTopoOut = erImTrackingTransform(movingTopo,fixed_out,tformSimilarity,tformAffine,0);
%Crop the reference topography to match moving outputs
fixedTopoOut=ercrop(referenceTopo,cropSize);

end

%% 
function [tformSimilarity,tformAffine,fixed] = erImTrackingTforms(moving_input,fixed_input,cropSize)
% This routine calculates the two transormation matrices (similarity and affine) that map the moving_input 
% image to the fixed_input image.
% cropSize corresponds to the desired output number of pixels. Useful if the edges are highly skewed.

% Remove line glitches. Smooth, crop.
moving=ercrop(ersmooth(erdescratch(moving_input,2),.75),cropSize);
fixed=ercrop(ersmooth(erdescratch(fixed_input,2),.75),cropSize);

% Convert to grayscale.
moving=mat2gray(moving); 
fixed=mat2gray(fixed);

% Define image registration objects that work for data acquired under
% slightly different conditions (multimodal vs monomodal). This takes care
% of slight tip changes and other changes between temperatures.

[optimizer, metric]=imregconfig('multimodal');

% Helps with convergence
optimizer.InitialRadius=optimizer.InitialRadius/5.5; 
% For better image registrations (sacrifices convergence time)
optimizer.MaximumIterations=1300;

% Get similarity transform corresponding to translation and scaling of the
% matrix.
tformSimilarity=imregtform(moving,fixed,'Similarity',optimizer,metric);

% Warp the moving image using tformSimilarity.
movingSimilarity=imwarp(moving,tformSimilarity,'OutputView',imref2d(size(fixed)));

% Get affine transform corresponding to similarity transform + shearing.
tformAffine=imregtform(movingSimilarity,fixed,'affine',optimizer,metric);

end

%%
function out = erImTrackingTransform(moving,fixed,tformSimilarity,tformAffine,f)
% This routine uses two tforms (tformSimilarity and tformAffine) from
% erImTrackingTforms to warp the matrix "moving."


cropSize=size(fixed,1);


%Flag for cleaning up the data. (Normalize lines by mean and deglitch).
if f==1 
    for i=1:size(moving,3)
        for j=1:size(moving,2)
            moving(:,j,i)=moving(:,j,i)./mean(squeeze(moving(:,j,i)));
        end
        moving(:,:,i)=erdeglitch(moving(:,:,i),2);        
    end
end


% Crop data, warp the moving image using tformSimilarity and tformAffine.
for i=1:size(moving,3)
    tmp=ercrop(moving(:,:,i),cropSize);
    tmp=imwarp(tmp,tformSimilarity,'OutputView',imref2d(size(fixed)));
    out(:,:,i)=imwarp(tmp,tformAffine,'OutputView',imref2d(size(fixed)));
end


end

%% 
function out=ercrop(data,L)
% Crops a matrix to size L while keeping the same center point.
% L can be one number for a square crop or an array [L1 L2].

% Program is clunky. Error if size(data) is even and L is odd and vice versa.
% Need to clean up.

if size(L,2)==2 %For non-square cropping
            L1=round(L(1)/2);
            L2=round(L(2)/2);
    sz1=round(size(data,1)/2)+1;
sz2=round(size(data,2)/2)+1;

for i=sz1-L1:sz1+L1-1
    for j=sz2-L2:sz2+L2-1
        out(i-sz1+L1+1,j-sz2+L2+1)=data(i,j);
    end
end

elseif size(data,1)==size(data,2) && size(L,2)==1

if mod(size(data,1),2)==0 && mod(size(data,2),2)==0 % data is even for both dims.
L=L/2;
sz1=round(size(data,1)/2)+1;
sz2=round(size(data,2)/2)+1;

for i=sz1-L:sz1+L-1
    for j=sz2-L:sz2+L-1
        out(i-sz1+L+1,j-sz2+L+1)=data(i,j);
    end
end
elseif mod(size(data,1),2)==1 && mod(size(data,2),2)==1 % data is odd for both dims.
    L=(L-1)/2;
    sz=(size(data,1)-1)/2;
    for i=sz+1-L:sz+L+1
        for j=sz+1-L:sz+L+1
            out(i-sz+L,j-sz+L)=data(i,j);
        end
    end
else
    L=round(L/2);
    sz1=round(size(data,1)/2)+1;
sz2=round(size(data,2)/2)+1;

for i=sz1-L:sz1+L-1
    for j=sz2-L:sz2+L-1
        out(i-sz1+L+1,j-sz2+L+1)=data(i,j);
    end
end
end
else
        L=round(L/2);
    sz1=round(size(data,1)/2)+1;
sz2=round(size(data,2)/2)+1;

for i=sz1-L:sz1+L-1
    for j=sz2-L:sz2+L-1
        out(i-sz1+L+1,j-sz2+L+1)=data(i,j);
    end
end
end



end

%%
function out = ersmooth(data,sigma)
%Applies low-pass filter gaussian of width sigma to input data.

sz=size(data,1);

[nx,ny]=ndgrid(1:sz,1:sz);

% Make gaussian centered at origin.
gaussmat=exp(-((nx-sz/2-1).^2+(ny-sz/2-1).^2)/2/(sigma)^2/2);

gaussmat=gaussmat./sum(sum(gaussmat)); %Normalize

out = conv2(data,gaussmat,'same'); % Convolve input data with gaussian.

end
%%
function out = erdeglitch(data,sigma)
% Removes glitches sigma standard deviations above or below the mean in an
% image. Replaces glitches with average of four nearest neighbors.
% Does not operate along edges.

mn = mean2(data);
stdev = std2(data);
limplus = mn+sigma*stdev;
limmin = mn-sigma*stdev;



for i = 2:size(data,1)-1
    for j = 2:size(data,2)-1
        if data(i,j) > limplus
            data(i,j) = (data(i+1,j)+data(i-1,j)+data(i,j+1)+data(i,j-1))/4;
        end
        if data(i,j) < limmin
            data(i,j) = (data(i+1,j)+data(i-1,j)+data(i,j+1)+data(i,j-1))/4;
        end
    end
end
out = data;
end

%%
function out = erdescratch(data,sigma)
% Removes line-glitches or "scratches" that are sigma standard deviations
% above or below the mean. The scratches are replaced by the average of the
% surrounding points.

mn = mean2(data);
stdev = std2(data);
limplus = mn+sigma*stdev;
limmin = mn-sigma*stdev;



for i = 2:size(data,1)-1
    for j = 2:size(data,2)-1
        if data(i,j)>limplus || data(i,j)<limmin % is point on a scratch?
            if data(i+1,j)>limplus || data(i+1,j)<limmin %is next point part of scratch?
                if data(i-1,j)>limplus || data(i-1,j)<limmin %was previous point part of scratch?
                    data(i,j)=(data(i,j-1)+data(i,j+1))/2; %replace by upper and lower lines
                else % currently at start of scratch - replace by previous + upper + lower
                    data(i,j)=(data(i-1,j)+data(i,j+1)+data(i,j-1))/3; 
                end
            elseif data(i-1,j)>limplus || data(i-1,j)<limmin % currently at end of scratch
                data(i,j)=(data(i+1,j)+data(i,j+1)+data(i,j-1))/3; 
            else % Just a point glitch - replace by 4 nearest neighbors.
                data(i,j)=(data(i-1,j)+data(i+1,j)+data(i,j-1)+...
                    data(i,j+1))/4;
            end
        end
    end
end

out=data;

end



