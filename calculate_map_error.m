function mapVars = calculate_map_error( maps, mapStrs )

maps = rot90(maps); % Rotate so that strip-like features from hysteresis run along vertical axis (seen with figure, imagesc(map))
NaNoutStd = 2;

Eu_or_Ce = 1; % 1=Eu, 0=Ce 
Xsz = size(maps,2); Ysz = size(maps,1);

mapVars = zeros(1,numel(mapStrs)); 
    
figure; 

if Eu_or_Ce == 1, mapPts = 521; % Eu05 has uncomplete xy-scan. Otherwise, data pnts = Xsz*Ysz
else, mapPts = size(maps,1)*size(maps,2);
end
botEdgeIdx = 1 : Ysz : mapPts; 
topEdgeIdx = Ysz : Ysz : mapPts; % Bottom edge idx, for throwing away pixel groups that are split into 2 columns

for map_i = 1:numel(mapStrs)
    mapStr = mapStrs{map_i};

    mapVals = maps(:,:,map_i);
    mapVals(find(abs(mapVals-nanmean(mapVals))>NaNoutStd*nanstd(mapVals))) = NaN;

    subplot(4,size(maps,3),map_i);
    imagesc(mapVals), axis xy, colormap hot, colorbar();
    title([mapStr]); 

    subplot(4,size(maps,3),map_i+3*size(maps,3));

    meanUltVariance = 0; % Converged value of variance
    meanVarianceAt50 = 0; % Check convergence of variance at 50 smaples
    meanVarianceAt20 = 0; % Check convergence of variance at 20 samples
    meanRandUltVariance = 0;  % Compare to variance of map of random values      
    iterCount = 0;

    for iter = 1:10   

    randMapVals = nanmin(mapVals(:))+(nanmax(mapVals(:))-nanmin(mapVals(:)))*rand(1,mapPts);%NaN*ones(size(mapVals));

    Ns = [10,20,50,100,200,500,1000,2000,5000];
    Variances = NaN*ones(size(Ns));
    randVariances = NaN*ones(size(Ns));

    % Create the array of data where each entry is the mean value of a strip of
    % map values (1x3 vertical column of map coors)
    for Ni = 1:numel(Ns)
        N = Ns(Ni);

        snipVariances = NaN*ones(1,N);
        randSnipVariances = NaN*ones(1,N);
        for i = 1:N     
            randIdx = round((2 + (mapPts-1-2)*rand()));

            randIdx = round((2 + (mapPts -2)*rand()));
            while ismember(randIdx-1,topEdgeIdx) || ismember(randIdx,botEdgeIdx)
                randIdx = round((2+(mapPts-2)*rand()));
            end
            snipVals = mapVals(randIdx-1 : randIdx);
            snipVals = snipVals(~isnan(snipVals));

            if numel(snipVals)==2
                randSnipVals = randMapVals(randIdx-1 : randIdx);

                snipVariance = sum((snipVals - mean(snipVals)).^2);
                randSnipVariance = sum((randSnipVals - mean(randSnipVals)).^2);

                snipVariances(i) = snipVariance;
                randSnipVariances(i) = randSnipVariance;
            end
        end

        Variances(Ni) = sqrt( nanmean( snipVariances ) );
        randVariances(Ni) = sqrt( nanmean( randSnipVariances ) );

    end
    log10Ns = log10(Ns);

    meanVarianceAt20 = meanVarianceAt20 + Variances(2);
    meanVarianceAt50 = meanVarianceAt50 + Variances(3);
    meanUltVariance = meanUltVariance + Variances(end);

    meanRandUltVariance = meanRandUltVariance + randVariances(end);

    iterCount = iterCount + 1;


    hold on, 
    plot(log10Ns, Variances, 'rx','HandleVisibility','off')
    hold on, plot(log10Ns, Variances, 'r'); 

    plot(log10Ns, randVariances, 'bo','HandleVisibility','off'); 
    hold on, plot(log10Ns, randVariances, 'b--');

    xlabel('log10 N'); ylabel('sqrt ( mean ( stripVariances ) )');

    end        
    meanUltVariance = meanUltVariance / iterCount;
    meanVarianceAt50 = meanVarianceAt50 / iterCount;
    meanVarianceAt20 = meanVarianceAt20 / iterCount;

    meanRandUltVariance = meanRandUltVariance / iterCount;

    if map_i==1, legend('Data','Random'); end
    title({['At N=20 ',num2str(meanVarianceAt20,'%.1e')];['At N=50 ',num2str(meanVarianceAt50,'%.1e')];['Ninf to ',num2str(meanUltVariance,'%.1e')]},'FontSize',8);

    % Plot histogram setting binWidth to variance
    subplot(4,size(maps,3),[map_i+1*size(maps,3),map_i+2*size(maps,3)]);
    histogram(mapVals,'FaceColor',[1,0,0]);
    
    % [N,x]=histcounts(mapVals,'BinWidth',meanUltVariance);
    hold on, histogram(mapVals, 'BinWidth',meanUltVariance,'FaceColor',[0,0,1]);

    if map_i == 1, legend('w/o OL','width=var'), end

    mapVars(map_i) = meanUltVariance; 
end

pause(.05);
% end

end