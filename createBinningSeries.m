function createBinningSeries( data )

maps = data.feature_results; 
mapNames = data.feature_names; 

% First, get outlier idx for all maps
NaNoutStd = 2;
badIdxs = [];
for map_i = 1:size(maps,2)
    badIdx = find( abs(maps(:,map_i)-nanmean(maps(:,map_i)))>NaNoutStd*nanstd(maps(:,map_i)) );
    badIdxs = [badIdxs; badIdx];
end
badIdxs = sort(unique(badIdxs));
data.badIdxs = badIdxs;

data.Eax = data.E_list(:,1) - data.FL_E; 
data.thax = data.th_list(:,1);
data.kax = nanmean(data.k_list,2);

data.hName = mapNames{data.hMapIdx};  data.vName = mapNames{data.vMapIdx};
hMap = maps(:,data.hMapIdx);          vMap = maps(:,data.vMapIdx);
% Put energies relative to FL 
if contains(data.hName, 'E'), hMap = hMap - data.FL_E; end 
if contains(data.vName, 'E'), vMap = vMap - data.FL_E; end 

% NaN out map outliers 
% vMap(find(abs(vMap-nanmean(vMap))>NaNoutStd*nanstd(vMap))) = NaN;
% hMap(find(abs(hMap-nanmean(hMap))>NaNoutStd*nanstd(hMap))) = NaN;
hMap(badIdxs) = NaN;
vMap(badIdxs) = NaN;
% Save processed maps to data
data.hMap = hMap; 
data.vMap = vMap;

if contains(data.mat_filename, 'Eu05')
    assignin('base','hMap_Eu',data.hMap); assignin('base','vMap_Eu',data.vMap); assignin('base','badIdxs_Eu',data.badIdxs);
elseif contains(data.mat_filename, 'Ce30')
    assignin('base','hMap_Ce',data.hMap); assignin('base','vMap_Ce',data.vMap); assignin('base','badIdxs_Ce',data.badIdxs);
end

% Assign rest of new gui params in data structure
data.binSpecs = [];         
data.cax1 = 0;                   data.cax2 = 1;
data.xlim1 = min(data.kax);      data.xlim2 = max(data.kax);
data.ylim1 = min(data.Eax);      data.ylim2 = max(data.Eax);
data.mdcE = -0.5;                data.mdcEw = 0.01; 
data.plotMdc = 0;                data.calcV5d = 0; 
data.v5d_E1 = -0.6;              data.v5d_E2 = -0.8; 

data.plotFig3 = 0;               data.plotFig4 = 0;
data.kF_E = 0.000; % Energy to take MDC for determining kF_SS
data.kF_Ew = 0.005; % Energy half-width for summing kF_SS MDC
data.kF_kw = 0.02; % k-width of gaussian blurring kF_SS MDC
data.E5d_E1 = -.50; % Range of energies to take MDCs for determining E5d
data.E5d_E2 = -.35; 
data.E5d_EN = 15; % Number of MDCs to take within range 
data.symmE5d = 0; % Symmetrize specs before calculating E5d

% For initial placeholder, make bin edges such that there's ~equal specs per bin
hSorted = sort(hMap); 
vSorted = sort(vMap);
hPerBin = ceil(numel(hMap(~isnan(hMap)))/data.hNumBins);
vPerBin = ceil(numel(vMap(~isnan(vMap)))/data.vNumBins);
hBinEdges = hSorted(1); 
vBinEdges = vSorted(1); 

for hBin_i = 1:data.hNumBins
    hBinEdges(end+1) = hSorted( min([find(isnan(hSorted),1,'first')-1, hBin_i*hPerBin]) );
end
for vBin_i = 1:data.vNumBins
    vBinEdges(end+1) =  vSorted( min([find(isnan(vSorted),1,'first')-1, vBin_i*vPerBin]) );
end
data.hBinWidth = data.mapVariances(data.hMapIdx);
data.vBinWidth = data.mapVariances(data.vMapIdx);

data.hBinEdges = hBinEdges; data.vBinEdges = vBinEdges;
for hBin_i = 1:numel(hBinEdges)
    eval(['data.hBinEdge',num2str(hBin_i),' = ',num2str(hBinEdges(hBin_i)),';']);
end
for vBin_i = 1:numel(vBinEdges)
    eval(['data.vBinEdge',num2str(vBin_i),' = ',num2str(vBinEdges(vBin_i)),';']);
end


ff = figure('Units','normalized','position',[.1,.5,.3,.4]);
guidata(ff,data);
pp1 = uipanel(ff, 'position',[.01,.01,.8,.8]); % Plot specs
pp2 = uipanel(ff, 'position',[.01,.85,.8,.15]); % Plot hMap distribution
pp3 = uipanel(ff, 'position',[.81,.01,.19,.8]); % Plot vMap distribution
pp4 = uipanel(ff, 'position',[.81,.81,.19,.19]); % More panel controls and plotting 
pp5 = uipanel(ff, 'position',[.01,.81,.8,.04]); % Panel for hBinEdges

ax3 = axes(pp3);
ax2 = axes(pp2); 


% Bin the specs with current params
uicontrol('parent',pp1,'style','pushbutton','string','Plot it','units','normalized','position',[.9,.95,.1,.05],'callback',{@binTheSpecs});

% Set color axis limits 
uicontrol('parent',pp1,'style','text','string','caxMin', 'units','normalized','position',[0.65,0.03, 0.1, 0.02]);
cax1 = uicontrol('parent',pp1,'style','slider', 'units','normalized','position',[0.6,0,0.15,0.03]);
cax1.Value = 0; cax1.Max = 1; cax1.Min = 0; cax1.Callback = {@editPlotVal,'cax1'};

uicontrol('parent',pp1,'style','text','string','caxMax', 'units','normalized','position',[0.85,0.03, 0.1, 0.02]);
cax2 = uicontrol('parent',pp1,'style','slider','units','normalized','position',[0.8,0,0.15,0.03]);
cax2.Value = 1; cax2.Max = 1; cax2.Min = 0; cax2.Callback = {@editPlotVal,'cax2'};

% Set x-axis limits 
uicontrol('parent',pp1,'style','text','string','xLim', 'units','normalized','position',[0.05,0.03, 0.1, 0.02]);
xlim1 = uicontrol('parent',pp1,'style','edit', 'units','normalized','position',[0,0,0.1,0.03]);
xlim1.String = min(data.kax); xlim1.Callback = {@editPlotStr, 'xlim1'};

xlim2 = uicontrol('parent',pp1,'style','edit', 'units','normalized','position',[0.11,0,0.1,0.03]);
xlim2.String = max(data.kax); xlim2.Callback = {@editPlotStr, 'xlim2'};

% Set y-axis limits 
uicontrol('parent',pp1,'style','text','string','yLim', 'units','normalized','position',[0.3,0.03, 0.1, 0.02]);
ylim1 = uicontrol('parent',pp1,'style','edit', 'units','normalized','position',[0.25,0,0.1,0.03]);
ylim1.String = min(data.Eax); ylim1.Callback = {@editPlotStr, 'ylim1'};

ylim2 = uicontrol('parent',pp1,'style','edit', 'units','normalized','position',[0.36,0,0.1,0.03]);
ylim2.String = max(data.Eax); ylim2.Callback = {@editPlotStr, 'ylim2'};

% Plot MDCs
mdcE_ui = uicontrol('parent',pp4,'style','edit','units','normalized','position',[.5, 0,.25,.2],'string',data.mdcE,'callback',{@editPlotStr,'mdcE'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.75, 0, 0.25, .2], 'string',data.mdcEw,'callback',{@editPlotStr,'mdcEw'});
uicontrol('parent',pp4,'style','pushbutton','units','normalized','position',[0,0,.5,.2],'string','Plot MDCs','callback',{@makeMDCs,'plotMdc'});

% Calculate v5d
uicontrol('parent',pp4,'style','pushbutton','units','normalized','position',[0,0.2,0.5,.2],'string','Calc v5d','callback',{@calc_v5d,'calcV5d'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.5, 0.2, 0.25, .2], 'string',data.v5d_E1,'callback',{@editPlotStr,'v5d_E1'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.75,0.2,0.25,0.2], 'string',data.v5d_E2,'callback',{@editPlotStr,'v5d_E2'});

% Plot E5ds
uicontrol('parent',pp4,'style','pushbutton','units','normalized','position',[0,0.4,0.15,0.2],'string','E5d','callback',{@plot_E5d});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.15,0.4,0.25,0.2],'string',data.E5d_E1,'callback',{@editPlotStr,'E5d_E1'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.40,0.4,0.25,0.2],'string',data.E5d_E2,'callback',{@editPlotStr,'E5d_E2'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.65,0.4,0.25,0.2],'string',data.E5d_EN,'callback',{@editPlotStr,'E5d_EN'});
uicontrol('parent',pp4,'style','checkbox','units','normalized','position',[0.90,0.4,0.1,0.2],'string','Symm','callback',{@editPlotVal,'symmE5d'});

% Plot SS kF
uicontrol('parent',pp4,'style','pushbutton','units','normalized','position',[0,0.6,0.15,0.2],'string','kSS','callback',{@plot_SSkF});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.15,0.6,0.25,0.2],'string',data.kF_E,'callback',{@editPlotStr,'kF_E'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.40,0.6,0.25,0.2],'string',data.kF_Ew,'callback',{@editPlotStr,'kF_Ew'});
uicontrol('parent',pp4,'style','edit','units','normalized','position',[0.65,0.6,0.25,0.2],'string',data.kF_kw,'callback',{@editPlotStr,'kF_kw'});


% hBinWidth 
uicontrol('parent', pp5, 'style','text','units','normalized','position',[0, 0.6, 0.1, 0.4],'string','bin width');
uicontrol('parent', pp5, 'style','edit','units','normalized','position',[0, 0.0, 0.1, 0.6],'string',data.hBinWidth, 'callback',{@editPlotStr,'hBinWidth'});

% hBinEdges
for hBin_i = 1:data.hNumBins+1
    xPos = (hBin_i)*(1/(data.hNumBins+2));
    uicontrol('parent',pp5,'style','edit','units','normalized','position',[xPos,0,0.1,1.0],'string',data.hBinEdges(hBin_i), 'callback',{@editPlotStr,['hBinEdge',num2str(hBin_i)]});
end

% vBinWidth 
uicontrol('parent', pp3, 'style','text','units','normalized','position',[0, 0.95, 0.3, 0.02],'string','bin width');
uicontrol('parent', pp3, 'style','edit','units','normalized','position',[0, 0.92, 0.3, 0.02],'string',data.vBinWidth, 'callback',{@editPlotStr,'vBinWidth'});

%vBinEdges 
for vBin_i = 1:data.vNumBins+1
    yPos = 1-(vBin_i)*(1/(data.vNumBins+2));
    uicontrol('parent',pp3, 'style','edit','units','normalized','position',[0,yPos,0.3,0.02],'string',data.vBinEdges(vBin_i), 'callback',{@editPlotStr,['vBinEdge',num2str(vBin_i)]});
end
  

%% Interactive functions
function editPlotStr(src,event,param)
    data = guidata(src); 
    data.(param) = str2num(get(src,'String'));
        
    guidata(src,data);
        
    if contains(param, 'Bin')
        % If changed some binning params, re-bin the spectra 
        plotHbins(data);
        plotVbins(data);
        
        guidata(src,data);
        binTheSpecs(src,event,param); 
    else % If not, just crop/edit the already-binned spectra  
        plotBinnedSpecs(data);
    end
end

function editPlotVal(src,event,param)
    data = guidata(src); 
    data.(param) = get(src,'Value');
    guidata(src,data);
    
    plotBinnedSpecs(data);
end
       
function calc_v5d(src, event, param)
    data = guidata(src);
    data.(param) = get(src,'Value');
        
    v5ds = calculate_v5d(data);
    data.v5ds = v5ds;
    guidata(src,data);
end

function makeMDCs(src, event, param)
    data = guidata(src);
    data.(param) = get(src,'Value'); 
    
    plotMDCs(data);
    guidata(src,data);
end    


function binTheSpecs( src, event, param )
    disp('Plotting the specs with params...');
    data = guidata(src);
    
    spec_list = data.spec_list; kShift = data.kShift; 
    kax = data.kax;             thax = data.thax;        
    hNumBins = data.hNumBins;   vNumBins = data.vNumBins;
    
    
    hBinEdges = zeros(1,hNumBins+1);
    for hi = 1:hNumBins+1
        hBinEdges(hi) = eval(['data.hBinEdge',num2str(hi)]);
    end
    vBinEdges = zeros(1,vNumBins+1);
    for vi = 1:vNumBins+1
        vBinEdges(vi) = eval(['data.vBinEdge',num2str(vi)]);
    end
    
    % Get spec idxs for each horizontal superbin 
    for hBin_i = 1:hNumBins
        hBinIdx_list = find(hMap >= hBinEdges(hBin_i) & hMap < hBinEdges(hBin_i+1));
        if hBin_i == hNumBins % For last superbin, include idxs == max bin edge 
            hBinIdx_list = find(hMap >= hBinEdges(hBin_i) & hMap <= hBinEdges(hBin_i+1));
        end
        hBinIdx{hBin_i} = hBinIdx_list;
    end
    % Get spec idxs for each vertical superbin 
    for vBin_i = 1:vNumBins
        vBinIdx_list = find(vMap >= vBinEdges(vBin_i) & vMap < vBinEdges(vBin_i+1));
        if vBin_i == vNumBins % For last superbin, include idxs == max bin edge 
            vBinIdx_list = find(vMap >= vBinEdges(vBin_i) & vMap <= vBinEdges(vBin_i+1));
        end
        vBinIdx{vBin_i} = vBinIdx_list;
    end
    
    
    if kShift == 1, [Kq,Eq] = meshgrid(kax,[1:size(spec_list,1)]); end
    
    binIdxs = cell(1,hNumBins*vNumBins);
    binSpecs = zeros(size(spec_list,1),size(spec_list,2), hNumBins*vNumBins);
    for vBin_i = 1:numel(vBinIdx)
        for hBin_i = 1:numel(hBinIdx)
            bin_i = (vBin_i-1)*numel(hBinIdx) + hBin_i;
            binIdx = intersect( vBinIdx{vBin_i}, hBinIdx{hBin_i} );
            binIdx = binIdx(~ismember(binIdx, data.badIdxs));
            binIdxs{bin_i} = binIdx;
            
            binSpec = zeros(size(spec_list,1),size(spec_list,2));
            for bin_ii = 1:numel(binIdx)
                spec_i = binIdx(bin_ii);
                spec = mat2gray( spec_list(:,:,spec_i) );
                spec(:,end-2:end) = 0;
%                 spec(spec>10*nanmean(spec(:))) = 0;
                if kShift==1
                    
                    theKax = data.k_list(:,spec_i);                
                    [K,E] = meshgrid(theKax, [1:size(spec_list,1)]);
                    theBinSpec = interp2(K, E, spec, Kq, Eq);
                    
                else, theBinSpec = spec;
                end
                binSpec = binSpec + mat2gray(theBinSpec);
            end
            binSpecs(:,:,bin_i) = (binSpec) / numel(binIdx);
        end
    end  
    
    data.hBinEdges = hBinEdges; 
    data.vBinEdges = vBinEdges;
    data.hBinIdx = hBinIdx;
    data.vBinIdx = vBinIdx;
    
    data.binSpecs = binSpecs; 
    data.binIdxs = binIdxs;
    
    guidata(src,data);
    
    plotHbins(data);
    plotVbins(data);
    plotBinnedSpecs(data);
    disp('Done binning.')
                
end

function plot_SSkF( src, event, param )
    % Plot SS kF versus E4f
    data = guidata(src);
    
    % First create the interpolated binned spectra
    if ~isfield(data, 'ex_binSpecs')
        binSpecs = data.binSpecs;
        kax = data.kax; Eax = data.Eax;

        kq = linspace(-0.65, 0.65, 500); 
        eq = linspace(-0.75, 0.10, 1000);
        [Kq, Eq] = meshgrid( kq, eq );
    
        ex_binSpecs = zeros( numel(eq), numel(kq), size(binSpecs,3) );
        for bi = 1:size(binSpecs,3)
            [K,E] = meshgrid( kax, Eax );
            ex_binSpecs(:,:,bi) = interp2( K, E, binSpecs(:,:,bi), Kq, Eq);
        end

        assignin('base','binSpecs',ex_binSpecs); assignin('base','binIdxs',data.binIdxs);
        assignin('base','binEax', eq); assignin('base','binKax',kq); assignin('base','f1pks',data.hMap);
        data.ex_binSpecs = ex_binSpecs; data.kq = kq; data.eq = eq;
    else
        ex_binSpecs = data.ex_binSpecs; kq = data.kq; eq = data.eq; 
        assignin('base','binSpecs',ex_binSpecs); assignin('base','binIdxs',data.binIdxs);
        assignin('base','binEax', eq); assignin('base','binKax',kq); assignin('base','f1pks',data.hMap);
    end
    
    if contains(data.mat_filename, 'Eu05'), Eu_or_Ce = 1; 
    elseif contains(data.mat_filename, 'Ce30'), Eu_or_Ce = 0;
    end
    
    kF_E = data.kF_E; symmetrize = 1; kF_conf_level = 0.68; 
    f_E1 = data.f_E1 - data.FL_E; f_E2 = data.f_E2-data.FL_E; f2_E1 = data.f2_E1-data.FL_E; f2_E2 = data.f2_E2-data.FL_E; 
    topX = 0.1; interpDiv = 100; 

    [~,~,~,~,~] = fit_SS_kF(Eu_or_Ce, data.hNumBins, data.vNumBins, ex_binSpecs, data.binIdxs, data.hMap, kq, eq, kF_E, symmetrize, kF_conf_level,...
        data.kF_Ew, data.kF_kw, f_E1, f_E2, f2_E1, f2_E2, data.f_k1, data.f_k2, data.fS_k1, data.fS_k2, topX, data.fS_X, interpDiv );
    
end


function plot_E5d( src, event, param )
    % Plot E5d versus E4f
    data = guidata(src); 
    
    % First create the interpolated binned spectra
    if ~isfield(data, 'ex_binSpecs')
        binSpecs = data.binSpecs;
        kax = data.kax; Eax = data.Eax;

        kq = linspace(-0.65, 0.65, 500); 
        eq = linspace(-0.75, 0.10, 1000);
        [Kq, Eq] = meshgrid( kq, eq );
    
        ex_binSpecs = zeros( numel(eq), numel(kq), size(binSpecs,3) );
        for bi = 1:size(binSpecs,3)
            [K,E] = meshgrid( kax, Eax );
            ex_binSpecs(:,:,bi) = interp2( K, E, binSpecs(:,:,bi), Kq, Eq);
        end

        assignin('base','binSpecs',ex_binSpecs); assignin('base','binIdxs',data.binIdxs);
        assignin('base','binEax', eq); assignin('base','binKax',kq); assignin('base','f1pks',data.hMap);
        data.ex_binSpecs = ex_binSpecs; data.kq = kq; data.eq = eq;
    else
        ex_binSpecs = data.ex_binSpecs; kq = data.kq; eq = data.eq; 
    end
    
    if contains(data.mat_filename, 'Eu05'), Eu_or_Ce = 1; 
    elseif contains(data.mat_filename, 'Ce30'), Eu_or_Ce = 0;
    end

    dBand_mdc_Es =  linspace( data.E5d_E1, data.E5d_E2, data.E5d_EN ) ; 

    [E4f, E4f_err, E5d, E5d_err, k5d, k5d_err, v5dL, v5dR] = track_dBand_Energies( Eu_or_Ce, dBand_mdc_Es, eq, kq, ex_binSpecs, ...
                    data.binIdxs, data.vNumBins, data.hNumBins, data.hMap, data.symmE5d);          
    % Save variables to base workspace to make combined Eu/Ce Fig.4a 
    if Eu_or_Ce == 1
        assignin('base','E4f_Eu',E4f);    assignin('base','E4f_err_Eu',E4f_err);
        assignin('base','E5d_Eu',E5d);    assignin('base','E5d_err_Eu',E5d_err);
        assignin('base','k5d_Eu',k5d);    assignin('base','k5d_err_Eu',k5d_err);
        assignin('base','v5dL_Eu',v5dL);  assignin('base','v5dR_Eu',v5dR);
    elseif Eu_or_Ce == 0
        assignin('base','E4f_Ce',E4f);    assignin('base','E4f_err_Ce',E4f_err);
        assignin('base','E5d_Ce',E5d);    assignin('base','E5d_err_Ce',E5d_err);
        assignin('base','k5d_Ce',k5d);    assignin('base','k5d_err_Ce',k5d_err);
        assignin('base','v5dL_Ce',v5dL);  assignin('base','v5dR_Ce',v5dR);
    end
    
    guidata(src,data);    
end



%% Functions to display what you have
% Plot the binned spectra
function plotBinnedSpecs(data)
    
    binSpecs = data.binSpecs; 
    binIdxs = data.binIdxs; 
    kax = data.kax; Eax = data.Eax;
    hNumBins = data.hNumBins; vNumBins = data.vNumBins; 
    
    ax1 = tiledlayout(pp1, vNumBins,hNumBins,'TileSpacing','compact');
    for bin_i = 1:vNumBins*hNumBins
        theAx = nexttile(ax1);

        imagesc(theAx, kax, Eax, binSpecs(:,:,bin_i)), set(theAx, 'YDir','normal');
        hold on, plot(theAx, [min(kax), max(kax)], data.mdcE-0.5*data.mdcEw*[1,1], 'w:');
        hold on, plot(theAx, [min(kax), max(kax)], data.mdcE+0.5*data.mdcEw*[1,1], 'w:');
        hold on, plot(theAx, [min(kax), max(kax)], data.v5d_E1*[1,1],'y:');
        hold on, plot(theAx, [min(kax), max(kax)], data.v5d_E2*[1,1],'y:');
        title(theAx, num2str(numel(binIdxs{bin_i})))
        
        if bin_i < (vNumBins-1)*hNumBins+1, set(theAx,'xTickLabel',[]); end
        if rem(bin_i,hNumBins) ~= 1, set(theAx,'yTickLabel',[]); end
        set(theAx,'TickDir','out');

        caxis(theAx, [data.cax1,data.cax2]);
        xlim(theAx, [data.xlim1,data.xlim2]); ylim(theAx, [data.ylim1,data.ylim2]);
    end
    colormap turbo 
    sgtitle(ax1, data.mat_filename, 'Interpreter','None');
end

function plotHbins( data )
    hMap = data.hMap; 
    hBinEdges = data.hBinEdges;
    hBinWidth = data.hBinWidth;
    hName = data.hName;
    
    if ~isempty(hBinWidth) % If specific binwidth provided
        [hN,hx] = histcounts(hMap,'BinWidth',hBinWidth);
    else % Else, just use a small binwidth, but large enough that there aren't any bins with counts less than 2
        hDispN = 11; [hN,~] = histcounts(hMap,'BinEdges',linspace(hBinEdges(1),hBinEdges(end),hDispN)); 
        while min(hN) > 1
            hDispN = hDispN+10; 
            [hN,~] = histcounts(hMap,'BinEdges',linspace(hBinEdges(1),hBinEdges(end),hDispN));
        end
        if hDispN == 11, hDispN = 21; end
        [hN,hx] = histcounts(hMap,'BinEdges',linspace(hBinEdges(1),hBinEdges(end),hDispN-10));
    end
    
    bar(ax2, 0.5*(hx(2:end)+hx(1:end-1)), hN, 1, 'FaceColor',[1,0,1]); 
    set(ax2,'TickDir','out');
    for hBin_i = 1:data.hNumBins+1
        hold(ax2,'on'); 
        plot(ax2, hBinEdges(hBin_i)*[1,1],[0,max(hN)*1.15],'k-','LineWidth',1); 
        ylim(ax2, [0,1.15*max(hN)]);
    end
    title(ax2, hName);
    xlim(ax2, [nanmin(hMap),nanmax(hMap)]);
    hold(ax2,'off');
    
    if contains(data.mat_filename, 'Eu05')
        assignin('base','hBinEdges_Eu',data.hBinEdges); 
    elseif contains(data.mat_filename, 'Ce30')
        assignin('base','hBinEdges_Ce',data.hBinEdges); 
    end
    
end

function plotVbins( data )
    % Plot the vertical bin distribution, superbin edges 
    vMap = data.vMap; 
    vBinEdges = data.vBinEdges;
    vBinWidth = data.vBinWidth;
    vName = data.vName;
    if ~isempty(vBinWidth) 
        [vN,vx] = histcounts(vMap,'BinWidth',vBinWidth);
    else
        % Start with 11 bins, increase by 10 until bin counts start getting too small
        vDispN = 11; [vN,~] = histcounts(vMap,'BinEdges',linspace(vBinEdges(1),vBinEdges(end),vDispN)); 
        while min(vN) > 1
            vDispN = vDispN+10; 
            [vN,~] = histcounts(vMap,'BinEdges',linspace(vBinEdges(1),vBinEdges(end),vDispN));
        end
        if vDispN== 11, vDispN = 21; end

        [vN,vx] = histcounts(vMap,'BinEdges',linspace(vBinEdges(1),vBinEdges(end),vDispN-10));
    end
    
    barh(ax3, 0.5*(vx(2:end)+vx(1:end-1)), vN, 1, 'FaceColor',[1,0,0]); 
    set(ax3,'TickDir','out');
    for vBin_i = 1:data.vNumBins+1
        hold(ax3,'on'); plot(ax3, [0,max(vN)*1.15],vBinEdges(vBin_i)*[1,1],'k-','LineWidth',1); 
        xlim(ax3, [0,1.15*max(vN)]);
    end
    set(ax3,'Ydir','reverse');
    title(ax3, vName);
    hold(ax3,'off');
    
    if contains(data.mat_filename, 'Eu05')
        assignin('base','vBinEdges_Eu',data.vBinEdges); 
    elseif contains(data.mat_filename, 'Ce30')
        assignin('base','vBinEdges_Ce',data.vBinEdges); 
    end
end
end