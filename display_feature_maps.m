
function [mapGrids,Xq,Yq] = display_feature_maps( filename, FL_E, X_list, Y_list, X_step, Y_step, feature_results, feature_names )


numMaps = size(feature_results,2);

[XQ,YQ,Xq,Yq] = get_XY_grid( X_list, Y_list, X_step, Y_step );
mapGrids = NaN*ones(size(XQ,1),size(XQ,2),size(feature_results,2));
figure,

for map_i = 1:numMaps
    
    map_list = feature_results(:,map_i);        
    mapName = feature_names{map_i};
    if contains(mapName, 'E'), map_list = map_list - FL_E; end
    
    notOutIdx = find(abs(nanmean(map_list)-map_list) < 2*nanstd(map_list));
    
    mapGrid = griddata(X_list, Y_list, map_list, XQ, YQ);
    mapGrids(:,:,map_i) = mapGrid;
    
    subplot(2,numMaps,map_i)
    imagesc(Xq,Yq,mapGrid);  daspect([1,1,1]); axis xy;
    caxis([ max([mean(map_list(notOutIdx))-2*std(map_list(notOutIdx)), min(map_list(notOutIdx))]),...
            min([nanmean(map_list(notOutIdx))+2*nanstd(map_list(notOutIdx)), max(map_list(notOutIdx))]) ]);
    
        title([num2str(map_i),'. ',feature_names{map_i}], 'Interpreter','None','FontSize',8 )
    
    subplot(2,numMaps, numMaps+map_i)
    histogram(map_list(notOutIdx));
    
end
colormap hot
sgtitle(filename, 'Interpreter','None', 'FontSize',8)

if contains(filename, 'Eu05')
    assignin('base','mapGrids_Eu',mapGrids); assignin('base','X_Eu',Xq); assignin('base','Y_Eu',Yq);
elseif contains(filename, 'Ce30')
    assignin('base','mapGrids_Ce',mapGrids); assignin('base','X_Ce',Xq); assignin('base','Y_Ce',Yq);
else 
    assignin('base','mapGrids',mapGrids); assignin('base','X',Xq); assignin('base','Y',Yq);
end


end







