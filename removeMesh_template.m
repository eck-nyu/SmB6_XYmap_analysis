function [newImage,allHotInds,finalVectList] =removeMesh_template(theData,threshVal,finalVectList)
% 1st order correction sample call:
% theInd=8; newImage=removeMesh(Idata_3d(:,:,theInd),[1e-3, 0*1e-3]);
% 2nd order correction:
% theInd=8; newImage=removeMesh(Idata_3d(:,:,theInd),[1e-3, 0*1e-3; 5e-4, 0*1e-3]);

% threshVal=[intensity threshold for identifying primary mesh peak, new
% intensity for the mesh feature; ...same parameters for secondary features...]
ptRadius=15;
plot_fig=0;
%Enter the 2 inequivalent Fourier coordinates here:
% theVectList=[85,50; 67,-65]; %inverse frame units
% theVectList = [429,31; 44,42];

theVectList = [300 - round(485/2), 165 - round(389/2);...
                286- round(485/2),236- round(389/2)]; % for 2020 11 ESM data

theVectList = [55,43;  72, 360-size(theData,2)]; % for 2021 Aug ESM data 
theVectList=[theVectList;-theVectList];
FFTimg=fft2(theData);

%set the threshold itself as the maximum value, if no max val was asigned
threshVal=threshVal*FFTimg(1,1); % set threshold to a fraction of average intensity
if size(threshVal,2)==1
    threshVal=[threshVal,threshVal];
end

if isempty(finalVectList)
    %This identifies the Fourier coordinates in the full image
    latticeVects=[0,0; size(FFTimg,1)+1, 0; 0, size(FFTimg,2)+1; size(FFTimg,1)+1, size(FFTimg,2)+1];
    finalVectList=[];
    for latticePl=1:size(latticeVects,1)
        for basisPl=1:size(theVectList,1)
            newPl=theVectList(basisPl,:)+latticeVects(latticePl,:);

            %if the place is in range:
            if newPl(1)>1 & newPl(1)<size(theData,1)
                if newPl(2)>1 & newPl(2)<size(theData,2)
                    finalVectList=[finalVectList; [newPl,threshVal(1,:)]];
                end
            end
        end
    end

    %now add in the 2nd order points
    % if size(threshVal,1)>1 %if there are parameters for the 2nd order
    if size(threshVal,1)>1
        if plot_fig==1, disp(['Doing second order']); end
        for primaryPl=1:size(finalVectList,1)
            for basisPl=1:size(theVectList,1)
                newPl=finalVectList(primaryPl,1:2)+theVectList(basisPl,:);
                if newPl(1)>1 & newPl(1)<size(theData,1)
                    if newPl(2)>1 & newPl(2)<size(theData,2)
                        distMat = sqrt( (newPl(1)-finalVectList(:,1)).^2 + (newPl(2)-finalVectList(:,2)).^2);
                        if min(distMat) ~= 0
                            finalVectList=[finalVectList; [newPl,threshVal(2,:)]];
                        end
                    end
                end
            end
        end
    end
    % finalVectList

    if size(threshVal,1)>2
    %     disp(['Doing third order'])
        for primaryPl=1:size(finalVectList,1)
            for basisPl=1:size(theVectList,1)
                newPl=finalVectList(primaryPl,1:2)+theVectList(basisPl,:);
                if newPl(1)>1 & newPl(1)<size(theData,1)
                    if newPl(2)>1 & newPl(2)<size(theData,2)
                        distMat = sqrt( (newPl(1)-finalVectList(:,1)).^2 + (newPl(2)-finalVectList(:,2)).^2);
                        if min(distMat) ~= 0
                            finalVectList=[finalVectList; [newPl,threshVal(3,:)]];
                        end
                    end
                end
            end
        end
    end
    if plot_fig==1, finalVectList, end
end

%Now apply the filter on the FFTimage
newImage=FFTimg;
allHotInds = [];
% figure, imagesc(abs(newImage)), caxis([0,1e6]) %***********
for theVectPl=1:size(finalVectList,1)
    theVect=finalVectList(theVectPl,:);
    
    %*** slow implementation for distance, and fails across boundary; should use template:
    Xcoords=[1:size(theData,1)]' * ones(1,size(theData,2));
    Ycoords=ones(1,size(theData,1))' * [1:size(theData,2)];
    distMat=sqrt((Xcoords-theVect(1)).^2+(Ycoords-theVect(2)).^2);

    actInds=find(distMat<ptRadius);
    hotInds=find(abs(newImage(actInds))>theVect(3)); hotInds=actInds(hotInds);
    newImage(hotInds)=theVect(4)*(newImage(hotInds)./abs(newImage(hotInds)));
    allHotInds = [allHotInds;hotInds];
%     
end
allHotInds = sort(unique(allHotInds));

if plot_fig==1
figure, 
cax = [3,15];
plot_colors = {'r+','g+','b+'};
subplot(221); imagesc(log(log(1+abs(FFTimg)))), axis xy; 
for tv_i = 1:size(threshVal,1)
    tv = threshVal(tv_i,1);
    tv_ii = find(finalVectList(:,3)==tv);
    for tv_iii = 1:numel(tv_ii)
        fv_i = tv_ii(tv_iii);
        hold on, plot(finalVectList(fv_i,2),finalVectList(fv_i,1),plot_colors{tv_i})
    end
end
subplot(222), imagesc(log(1+abs(newImage))), axis xy; 
%caxis([0,mean(abs(newImage(:)))])
end

newImage=abs(fft2(newImage));
newImage=newImage(end:-1:1, end:-1:1);

if plot_fig==1
subplot(223), imagesc(theData), axis xy; 
subplot(224), imagesc(newImage), axis xy;
sgtitle(mat2str(threshVal/FFTimg(1,1)))
end
end

