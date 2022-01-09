function XYmap_gui()

%% NOTE: For some reason, callback functions work better with input arguments src, event, param (even if unused)

% Get name of available .mat files in the directory for gui drop-down list 
matdata_list = dir('*_data*.mat');
matdata_names = cell(1,numel(matdata_list));
for i = 1:numel(matdata_names)
    matdata_names{i} = matdata_list(i).name; 
end
data.matdata_names = matdata_names;

%% Itx folders have set of "raw" itx files from XY scan at ESM
% If converting from itx folder, push "Convert itx" then re-run gui to access .mat file
data.itx_foldername = '/home/data/eck325/SmB6/BNL_SmB6_20210808_all_the_data/XY120KSmB6_itx/';
% Each XY scan should have corresponding csv file with manipulator XY position data 
data.csv_filename = '/home/data/eck325/SmB6/BNL_SmB6_20210808_all_the_data/meta_XY120KSmB6.csv';
        % name_col_idx = 2;  X_col_idx = 6;  Y_col_idx = 7;
      
data.mat_filename = data.matdata_names{1}; % Placeholder for current mat file 
data.convert_itx = 0;  % Convert folder of itx to a mat file
data.import_mat = 0; % Import a new selected mat file
data.run_params = 0; % Calculate maps with the current set of params 
data.save_params = 0; % Option to save set of params to mat file
data.import_params = 1; % Import the params saved to mat file with mat file

% Placeholder params: 
data.hMapIdx = {'(No maps yet)'};
data.vMapIdx = {'(No maps yet)'};
data.hNumBins = 1; % Number of horizontal bins for binning series 
data.vNumBins = 1; % Number of vertical bins for binning series 
data.kShift = 0; % 1=Add k-correction using tilt map
data.convert2k = 0; % Convert theta to k

data.feature_results = []; % 2d-array of feature map results 
data.feature_names = {}; % string of feature names

% All energies in this demo are kinetic energy 
data.FL_E = 65.71; % Energy of Fermi level
data.bg_E1 = 65.1;   data.bg_E2 = 65.5; % Energy window to find background intensity
data.bg_th1 = -5;    data.bg_th2 = 8; % Theta window to find background intensity
data.bg_k1 = 0.5;    data.bg_k2 = 0.65; % If using converted k-axis, k window to find background intensity
data.f_E1 = 65.64;   data.f_E2 = 65.74; % Energy, theta, k windows to search for J=5/2 f-band peak in EDC 
data.f_th1 = -100;   data.f_th2 = 100; 
data.f_k1 = -0.6;    data.f_k2 = 0.6; 
data.fS_th1 = -100;  data.fS_th2 = 100; % theta, k windows to sum EDC for f-band peak width 
data.fS_k1 = -0.1;   data.fS_k2 = 0.1; 
data.f2_E1 = 65.45;  data.f2_E2 = 65.63; % Energy window to search for J=7/2 f-band peak
data.Ism_E1 = 64.9;  data.Ism_E2 = 65.23; % Energy, theta, k windows to find Sm-multiplet feature intensity
data.Ism_th1 = -13;  data.Ism_th2 = -6; 
data.Ism_k1 = 0.5;   data.Ism_k2 = 0.7; 
data.d_E1 = 65.2;    data.d_E2 = 65.45; % Energy window to track d-bands
data.dL_th1 = -5;    data.dL_th2 = 0; % theta window to search for left d-band crossing
data.dR_th1 = 3;     data.dR_th2 = 10; % theta window to search for right d-band crossing
data.fS_X = 0.5; % X in FWXM to calculate f-band peak width (self-energy)

% ARPES geometry params (assume constant throughout dataset)
data.azi0 = 0; % sample azimuth 
data.pol0 = 0; % sample polar 
data.hv = 70; % photon energy 
data.WF = 4.4; % work function 
data.EB = 0; % binding energy
data.V0 = 10; % inner potential

data.spec_list = []; % 3d-array of ARPES images (1/row=energy, 2/column=theta, 3=image)
data.E_list = []; % 2d-array of ARPES image energy axis (1/row=energy, 2/column=image)
data.th_list = []; % 2d-array of ARPES image theta axis (1/row=theta, 2/column=image)
data.X_list = []; % 1d-array of manipulator X positions
data.Y_list = []; % 1d-array of manipulator Y positions
data.X_step = 0.01; % Manipulator X step
data.Y_step = 0.01; % Manipulator Y step
data.show_samples = 10; % Show random sampling of spectra

f = figure('Units','normalized','position',[.1,.1,.3,.5]);
guidata(f, data)

p1 = uipanel(f, 'position',[0.0, 0.9, 1.0, 0.1]); % Panel for folder/data info 
p2 = uipanel(f, 'position',[0.0, 0.0, 0.2, 0.9]); % Panel for all the calculation params 
p3 = uipanel(f, 'position',[0.2, 0.85, 0.80, 0.05]); % Panel for experimental setup/geometry params 
p4 = uipanel(f, 'position',[0.2, 0.1, 0.80, 0.75]); % Panel for displaying spec and current params
p5 = uipanel(f, 'position',[0.2, 0.0, 0.80, p4.Position(2)]); % Panel for binning series 

p4_t = tiledlayout(p4,1,1);

ax4_2 = axes(p4_t);
ax4 = axes(p4_t);

%%%%%% FOLDER DATA STUFF %%%%%%
uicontrol('parent',p1,'style','edit','units','normalized',      'position',[.02,.5,.8,.4],...
    'string',data.itx_foldername, 'callback', {@gotFoldername, 'itx_foldername'} );

uicontrol('parent',p1,'style','pushbutton','units','normalized','position',[.82,.5,.1,.4],...
    'string','Convert itx', 'callback',{@convertItx2Mat,'convert_itx'});

% Import spectra from specified saved .mat file 
uicontrol('parent',p1,'style','popupmenu', 'units','normalized','position',[.10,.02,.50,.40],...
    'string',data.matdata_names,'callback',{@gotMenuSelection, 'mat_filename'} );
uicontrol('parent',p1,'style','pushbutton','units','normalized','position',[.60,.02,.20,.40],...
    'string','Import mat',  'callback',{@importMatData, 'import_mat'});
% Import pre-saved parameters with mat file
uicontrol('parent',p1,'style','pushbutton','units','normalized', 'position',[.80,.02,.20,.40],...
    'String','Import params','callback',{@importParams,'import_params'});
 
%%%%%%% EXPERIMENTAL SETUP GEOMETRY %%%%%%%%%
% azi0 - sample azimuthal rotation (degrees)
uicontrol('parent',p3, 'style','text','string','azi0','units','normalized','position',[0.0, 0.5, 0.1, 0.5]);
uicontrol('parent',p3, 'style','edit','string',data.azi0, 'units','normalized','position',[0.0, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'azi0'});
% pol0 - DA30 theta-Y offset (degrees)
uicontrol('parent',p3, 'style','text','string','pol0','units','normalized','position',[0.1, 0.5, 0.1, 0.5]);
uicontrol('parent',p3,'style','edit','string',data.pol0, 'units','normalized','position',[0.1, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'pol0'});
% hv - photon energy (eV)
uicontrol('parent',p3, 'style','text','string','hv','units','normalized','position',[0.2, 0.5, 0.1, 0.5]);
uicontrol('parent',p3,'style','edit','string',data.hv, 'units','normalized','position',[0.2, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'hv'});
% WF - work function (eV)
uicontrol('parent',p3, 'style','text','string','WF','units','normalized','position',[0.3, 0.5, 0.1, 0.5]);
uicontrol('parent',p3,'style','edit','string',data.WF, 'units','normalized','position',[0.3, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'WF'});
% EB - binding energy (approximating a single value for spectra) (eV)
uicontrol('parent',p3, 'style','text','string','EB','units','normalized','position',[0.4, 0.5, 0.1, 0.5]);
uicontrol('parent',p3,'style','edit','string',data.EB, 'units','normalized','position',[0.4, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'EB'});
% V0 - inner potential (eV)
uicontrol('parent',p3, 'style','text','string','V0','units','normalized','position',[0.5, 0.5, 0.1, 0.5]);
uicontrol('parent',p3,'style','edit','string',data.V0, 'units','normalized','position',[0.5, 0, 0.1, 0.5], 'callback',{@updateParameterValue,'V0'});
% Convert axes with current geometry / calc parameters
uicontrol('parent',p3, 'style','pushbutton','string','Convert axes','units','normalized','position',[.7,0,.3,1],'callback',{@convertAxes, 'convert2k'});


%%%% PARAMETERS %%%%

p2N = 20;  p2h = 1/p2N; % Set height of param boxes 

% FL energy - set manually, since setting max slope of edc might be too low? 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-1)*p2h, .3, p2h], 'string','FL E');
FL_E_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-1)*p2h, .59, p2h], ...
                    'string',data.FL_E, 'callback',{@updateParameterValue,'FL_E'});

% Energy, theta range for calculating background intensity (used in Ism) 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-2)*p2h, .3, p2h],  'string','BG E');
bg_E1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-2)*p2h, .3, p2h], ...
                    'string',data.bg_E1,'callback',{@updateParameterValue,'bg_E1'});
bg_E2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-2)*p2h, .3, p2h],  ...
                    'string',data.bg_E2,'callback',{@updateParameterValue,'bg_E2'});

uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-3)*p2h, .3, p2h],  'string','BG th');
bg_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-3)*p2h, .3, p2h], ...
                    'string',data.bg_th1,'callback',{@updateParameterValue,'bg_th1'});
bg_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-3)*p2h, .3, p2h], ...
                    'string',data.bg_th2,'callback',{@updateParameterValue,'bg_th2'});

% option to integrate k-range for bacground
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-4)*p2h, .3, p2h],  'string','BG k');
bg_k1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-4)*p2h, .3, p2h], ...
                    'string',data.bg_k1,'callback',{@updateParameterValue,'bg_k1'});
bg_k2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-4)*p2h, .3, p2h], ...
                    'string',data.bg_k2,'callback',{@updateParameterValue,'bg_k2'});

                
                
% Energy range to search for f (J=5/2) peak
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-5)*p2h, .3, p2h],  'string','f1 E');
f_E1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-5)*p2h, .3, p2h], ...
                    'string',data.f_E1,'callback',{@updateParameterValue,'f_E1'});
f_E2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-5)*p2h, .3, p2h], ...
                    'string',data.f_E2,'callback',{@updateParameterValue,'f_E2'});

% theta range for integrating EDC to calculate f-band energies
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-6)*p2h, .3, p2h],  'string','f th');
f_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-6)*p2h, .3, p2h],...
                    'string',data.f_th1,'callback',{@updateParameterValue,'f_th1'});
f_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-6)*p2h, .3, p2h],...
                    'string',data.f_th2,'callback',{@updateParameterValue,'f_th2'});

% Option: k range for integrating EDC to calculate f-band energies
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-7)*p2h, .3, p2h],  'string','f k');
f_k1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-7)*p2h, .3, p2h],...
                    'string',data.f_k1,'callback',{@updateParameterValue,'f_k1'});
f_k2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-7)*p2h, .3, p2h],...
                    'string',data.f_k2,'callback',{@updateParameterValue,'f_k2'});

% Energy range to search for f2 (J=7/2) peak
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-8)*p2h, .3, p2h], 'string','f2 E');
f2_E1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-8)*p2h, .3, p2h], ...
                    'string',data.f2_E1,'callback',{@updateParameterValue,'f2_E1'});
f2_E2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-8)*p2h, .3, p2h], ...
                    'string',data.f2_E2,'callback',{@updateParameterValue,'f2_E2'});

% theta range for integrating EDC to calculate f-peak widths 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-9)*p2h, .3, p2h],  'string','fS th');
fS_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-9)*p2h, .3, p2h], ...
                    'string',data.fS_th1,'callback',{@updateParameterValue,'fS_th1'});
fS_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-9)*p2h, .3, p2h], ...
                    'string',data.fS_th2,'callback',{@updateParameterValue,'fS_th2'});
                
% option: k range for integrating EDC to calculate f-peak widths 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-10)*p2h, .3, p2h],  'string','fS k');
fS_k1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-10)*p2h, .3, p2h], ...
                    'string',data.fS_k1,'callback',{@updateParameterValue,'fS_k1'});
fS_k2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-10)*p2h, .3, p2h], ...
                    'string',data.fS_k2,'callback',{@updateParameterValue,'fS_k2'});

% E, th windows for calculating Ism intensity 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-11)*p2h, .3, p2h],  'string','Ism E');
Ism_E1_ui = uicontrol('parent',p2,'style','edit','units','normalized', 'position',[.31, (p2N-11)*p2h, .3, p2h],...
                    'string',data.Ism_E1,'callback',{@updateParameterValue,'Ism_E1'});
Ism_E2_ui = uicontrol('parent',p2,'style','edit','units','normalized', 'position',[.61, (p2N-11)*p2h, .3, p2h],...
                    'string',data.Ism_E2,'callback',{@updateParameterValue,'Ism_E2'});

uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-12)*p2h, .3, p2h],  'string','Ism th');
Ism_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-12)*p2h, .3, p2h],...
                    'string',data.Ism_th1,'callback',{@updateParameterValue,'Ism_th1'});
Ism_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-12)*p2h, .3, p2h],...
                    'string',data.Ism_th2,'callback',{@updateParameterValue,'Ism_th2'});

% option: k range for integrating Ism window             
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-13)*p2h, .3, p2h],  'string','Ism k');
Ism_k1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-13)*p2h, .3, p2h],...
                    'string',data.Ism_k1,'callback',{@updateParameterValue,'Ism_k1'});
Ism_k2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-13)*p2h, .3, p2h],...
                    'string',data.Ism_k2,'callback',{@updateParameterValue,'Ism_k2'});

                
                
% E range for finding d-band peaks                 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-14)*p2h, .3, p2h],  'string','d E');
d_E1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-14)*p2h, .3, p2h],...
                    'string',data.d_E1,'callback',{@updateParameterValue,'d_E1'});
d_E2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-14)*p2h, .3, p2h],...
                    'string',data.d_E2,'callback',{@updateParameterValue,'d_E2'});

% th range for finding LHS d-band peaks                 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-15)*p2h, .3, p2h],  'string','dL th');
dL_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-15)*p2h, .3, p2h],...
                    'string',data.dL_th1,'callback',{@updateParameterValue,'dL_th1'});
dL_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-15)*p2h, .3, p2h],...
                    'string',data.dL_th2,'callback',{@updateParameterValue,'dL_th2'});

% th range for finding RHS d-band peaks                 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-16)*p2h, .3, p2h],  'string','dR th');
dR_th1_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-16)*p2h, .3, p2h],...
                    'string',data.dR_th1,'callback',{@updateParameterValue,'dR_th1'});
dR_th2_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-16)*p2h, .3, p2h],...
                    'string',data.dR_th2,'callback',{@updateParameterValue,'dR_th2'});

% X,Y step used in XY manipulator scan
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-17)*p2h, .3, p2h], 'string','XYdel');
X_step_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.31, (p2N-17)*p2h, .3, p2h],...
                    'string',data.X_step,'callback',{@updateParameterValue,'X_step'});
Y_step_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-17)*p2h, .3, p2h],...
                    'string',data.Y_step,'callback',{@updateParameterValue,'Y_step'});

% Fraction X in FWXM used to find f-peak widths 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-18)*p2h, .3, p2h], 'string','fS X');
fS_X_ui = uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-18)*p2h, .3, p2h], ...  
                                'string',data.fS_X,'callback',{@updateParameterValue,'fS_X'});

% Number of sample spectra/calcs to display 
uicontrol('parent',p2,'style','text','units','normalized','position',[.01, (p2N-19)*p2h, .6, p2h], 'string','Show samp.');
uicontrol('parent',p2,'style','edit','units','normalized','position',[.61, (p2N-19)*p2h, .3, p2h], 'string',data.show_samples,'callback',{@updateParameterValue,'show_samples'});

% Run calculations with current parameter set
uicontrol('parent',p2,'style','pushbutton','units','normalized', 'position',[.6,.0,.4,p2h], 'string','Run', 'callback',{@runParams,'run_params'});
% Save parameter set 
uicontrol('parent',p2,'style','pushbutton','units','normalized', 'position',[.1,.0,.4,p2h], 'string','Save', 'callback',{@saveParams,'save_params'});


%%%%%% Panel for plotting spec, params %%%%%%%%


% Remove mesh (edit threshold parameters in function code) 
uicontrol('parent',p4,'style','pushbutton','units','normalized','position',[0.8,0.95,.2,.05],...
    'String','Remove mesh', 'callback',{@removeTheMesh});

% Replot the feature maps 
uicontrol('parent',p4,'style','pushbutton','units','normalized','position',[0.8, 0.0, .2, .05],...
    'String','Replot maps', 'callback',{@replotMaps});

% Calculate feature map error 
uicontrol('parent',p4,'style','pushbutton','units','normalized','position',[0.6, 0.0, 0.2, 0.05], ...
    'string','Map error','callback',{@calculateMapError});

%%%%%%%%%% BINNING SERIES %%%%%%%%%%

% Choose horizontal bins map
hMap_ui = uicontrol('parent',p5,'style','popupmenu', 'units','normalized','position',[0,0,.3,.6], 'String',data.hMapIdx,'callback',{@gotValueSelection,'hMapIdx'});
             uicontrol('parent',p5,'style','text','units','normalized','position',[0,0.6,.3,.4],'String','H axis');
uicontrol('parent',p5,'style','edit', 'units','normalized','position',[.3,0,.12,.6], 'String',data.hNumBins,'callback',{@gotString2ValueSelection,'hNumBins'});
           uicontrol('parent',p5,'style','text','units','normalized','position',[.3,.6,.12,.4],'String','Num bins');

% Choose vertical bins map
vMap_ui = uicontrol('parent',p5,'style','popupmenu', 'units','normalized','position',[.425,0,.3,.6], 'String',data.vMapIdx,'callback',{@gotValueSelection,'vMapIdx'});
             uicontrol('parent',p5,'style','text','units','normalized','position',[.425,0.6,.3,.4],'String','V axis');
uicontrol('parent',p5,'style','edit', 'units','normalized','position',[.725,0,.12,.6], 'String',data.vNumBins, 'callback',{@gotString2ValueSelection,'vNumBins'});
            uicontrol('parent',p5,'style','text','units','normalized','position',[.725,.6,.12,.4],'String','Num bins');
            
% Make the binning series
uicontrol('parent',p5,'style','pushbutton','units','normalized', 'position',[.85,0,.15,.6], 'String','Make BS',...
            'callback',{@makeBinningSeries_call});
        
% k-correct individual spectra while making binning series
uicontrol('parent',p5,'style','checkbox','units','normalized','position',[.85,.7, .15, .3], 'String','k-shift',...
            'callback',{@gotValueSelection,'kShift'});


%%%%%%%%%% FUNCTIONS %%%%%%%%%%

    function [feature_results, feature_names, k_list ] = runParams(src, event, param)
    % With updated params, produce the new feature maps    
        data = guidata(src);
        data.(param) = get(src,'Value');
        
        disp('Calculating the feature maps...');
        [ feature_results, feature_names, k_list ] = run_the_params( data );

        pause(.1);
                                       
        data.feature_results = feature_results; 
        data.feature_names = feature_names;
        data.k_list = k_list;
        
        % Update the popup menu for making binning series, assign results in base workspace 
        set(hMap_ui,'String',feature_names); set(vMap_ui,'String',feature_names);     
        assignin('base','maps',feature_results); assignin('base','mapNames',feature_names); assignin('base','k_list',k_list);
            disp('Done.');
         
        [mapGrids,X,Y] = display_feature_maps( data.mat_filename, data.FL_E, data.X_list, data.Y_list, data.X_step, data.Y_step, data.feature_results, data.feature_names );
        
        data.mapGrids = mapGrids; data.X = X; data.Y = Y;
        
        
        
        guidata(src,data);      
    end
    
    function convertAxes(src, event, param)
        data = guidata(src);
        
        % Make spec_list contain only one spatially-averaged spectra, feed
        % into run_the_params to get the output kax, then plot on top of 
        spec_list = data.spec_list; 
        
        mn_spec = nanmean(data.spec_list,3);
        
        data.spec_list = mn_spec;
        data.show_samples = 0; 
        
        [~,~,mn_kax] = run_the_params( data );
        
        mn_Eax = data.E_list(:,1) - data.FL_E;
        
        xlim(ax4_2, [(mn_kax(1)),(mn_kax(end))]); 
        ylim(ax4_2, [mn_Eax(1),mn_Eax(end)]);
        xlabel(ax4_2,'k (invA)'); ylabel(ax4_2,'E-EF (eV)');
        ax4_2.XAxisLocation = 'top';
        ax4_2.YAxisLocation = 'right';
        ax4_2.YDir = 'normal';
        
        % Without saving data, run updateDisplaySpec 
        data.spec_list = spec_list; 
        data.kax = mn_kax; 
        guidata(src, data);
        
        updateDisplaySpec(src, event, param); 
    end

    function updateDisplaySpec(src, event, param)      
        data = guidata(src);

        mn_spec = mat2gray(nanmean(data.spec_list,3)); thax = data.th_list(:,1); Eax = data.E_list(:,1);
        edc = sum(mn_spec,2); edc = (edc-min(edc))/(max(edc)-min(edc));
        edc = edc * (max(thax)-min(thax)); edc = edc + min(thax);
        eax = linspace(Eax(1),Eax(end),numel(Eax)*1); edc = interp1(Eax, edc, eax,'spline');
        edc_diff = diff(edc); edc_diff = (edc_diff-min(edc_diff))/(max(edc_diff)-min(edc_diff)) * range(edc) + min(edc); 
        
        imagesc(ax4, thax, Eax, mn_spec); ax4.YDir = 'normal'; hold(ax4, 'on');
        plot(ax4, edc, eax, 'w'); 
        plot(ax4, edc_diff, eax(1:end-1), 'r');
        plot(ax4, [min(thax),max(thax)], data.FL_E*[1,1],'w:');  text(ax4, min(thax),data.FL_E, 'FL','Color','w');
        
        if data.bg_k2-data.bg_k1 > 0 && isfield(data,'kax')
            kax = data.kax;
            bg_k1_idx = find(kax >= data.bg_k1, 1,'first'); 
            bg_k2_idx = find(kax <= data.bg_k2, 1,'last');
            rectangle(ax4,'Position',[thax(bg_k1_idx), data.bg_E1, thax(bg_k2_idx)-thax(bg_k1_idx), data.bg_E2-data.bg_E1],'EdgeColor',[.8,.8,.8]);
            text(ax4, thax(bg_k1_idx), data.bg_E1, 'bg','Color',[.8,.8,.8]);
        else       
            rectangle(ax4, 'Position',[data.bg_th1,data.bg_E1,data.bg_th2-data.bg_th1,data.bg_E2-data.bg_E1],'EdgeColor',[.8,.8,.8]);
            text(ax4, data.bg_th1, data.bg_E1, 'bg','Color',[.8,.8,.8]);
        end
        if data.Ism_k2-data.Ism_k1 > 0 && isfield(data,'kax')
            Ism_k1_idx = find(kax >= data.Ism_k1,1,'first');
            Ism_k2_idx = find(kax <= data.Ism_k2,1,'last');
            rectangle(ax4,'Position',[thax(Ism_k1_idx),data.Ism_E1,thax(Ism_k2_idx)-thax(Ism_k1_idx),data.Ism_E2-data.Ism_E1],'EdgeColor',[.1,.1,.1]);
            text(ax4, thax(Ism_k1_idx), data.Ism_E1, 'Ism','Color',[.1,.1,.1]);
        else 
            rectangle(ax4, 'Position', [data.Ism_th1,data.Ism_E1,data.Ism_th2-data.Ism_th1,data.Ism_E2-data.Ism_E1],'EdgeColor',[.1,.1,.1]);
            text(ax4, data.Ism_th1, data.Ism_E1, 'Ism','Color',[.1,.1,.1]);
        end
        if data.f_k2-data.f_k1 > 0 && isfield(data,'kax')
            f_k1_idx = find(kax >= data.f_k1, 1,'first'); 
            f_k2_idx = find(kax <= data.f_k2, 1,'last');
            rectangle(ax4,'Position',[thax(f_k1_idx), data.f_E1, thax(f_k2_idx)-thax(f_k1_idx), data.f_E2-data.f_E1],'EdgeColor',[1,1,0]);
            text(ax4, thax(f_k1_idx), data.f_E1, 'f1','Color',[1,1,0]);
            rectangle(ax4,'Position',[thax(f_k1_idx), data.f2_E1, thax(f_k2_idx)-thax(f_k1_idx), data.f2_E2-data.f2_E1],'EdgeColor',[1,.5,0]);
            text(ax4, thax(f_k1_idx), data.f2_E1, 'f2','Color',[.8,.8,.8]);
        else       
            rectangle(ax4, 'Position',[data.f_th1,data.f_E1,data.f_th2-data.f_th1,data.f_E2-data.f_E1],'EdgeColor',[1,1,0]);
            text(ax4, data.f_th1, data.f_E1, 'f1','Color',[1,1,0]);
            rectangle(ax4, 'Position',[data.f_th1,data.f_E1,data.f_th2-data.f_th1,data.f_E2-data.f_E1],'EdgeColor',[1,.5,0]);
            text(ax4, max([min(thax),data.f_th1]), data.f_E1, 'f2','Color',[1,.5,0]);
        end 
         
        if data.fS_k2-data.fS_k1 > 0 && isfield(data,'kax')
            fS_k1_idx = find(kax >= data.fS_k1, 1,'first'); 
            fS_k2_idx = find(kax <= data.fS_k2, 1,'last');
            rectangle(ax4,'Position',[thax(fS_k1_idx), data.f_E1, thax(fS_k2_idx)-thax(fS_k1_idx), data.f_E2-data.f_E1],'EdgeColor',[1,1,0]);
            text(ax4, thax(fS_k1_idx), data.f_E1, 'f1,2S','Color',[1,1,0]);
        else       
            rectangle(ax4, 'Position',[data.fS_th1,data.f_E1,data.fS_th2-data.fS_th1,data.f_E2-data.f_E1],'EdgeColor',[1,1,0]);
            text(ax4, max([min(thax),data.fS_th1]), data.f_E1, 'f1,2S','Color',[1,1,0]);
        end   
        rectangle(ax4, 'Position', [data.dL_th1,data.d_E1,data.dL_th2-data.dL_th1,data.d_E2-data.d_E1],'EdgeColor',[0,0.9,0]);
        rectangle(ax4, 'Position', [data.dR_th1,data.d_E1,data.dR_th2-data.dR_th1,data.d_E2-data.d_E1],'EdgeColor',[0,0.9,0]);
        
        text(ax4, data.dL_th1, data.d_E1, 'dL', 'Color',[0,0.9,0]);
        text(ax4, data.dR_th1, data.d_E1, 'dR', 'Color',[0,0.9,0]);
        
        colormap(ax4, turbo)
        
        caxis([0.1,0.65]);
        xlabel('theta (deg)'); ylabel('KE (eV)');
    end



    function makeBinningSeries_call(src, event, param)
        data = guidata(src);
        data.sameBinWeightOrWidth = 1; 
        
        guidata(src,data);

        createBinningSeries( data );
    end

    function saveParams(src, event, param)
        % Store current set of parameters with the current mat file
        data = guidata(src);
        data.(param) = get(src,'Value');
        if data.save_params == 1
            mat_filename = data.mat_filename;
            param_filename = [mat_filename(1:end-8),'params.mat'];
            
            FL_E = data.FL_E; 
            bg_E1 = data.bg_E1;     bg_E2 = data.bg_E2;
            bg_th1 = data.bg_th1;   bg_th2 = data.bg_th2;
            bg_k1 = data.bg_k1;     bg_k2 = data.bg_k2; 
            f_E1 = data.f_E1;       f_E2 = data.f_E2;
            f_th1 = data.f_th1;     f_th2 = data.f_th2;
            f_k1 = data.f_k1;       f_k2 = data.f_k2; 
            fS_th1 = data.fS_th1;   fS_th2 = data.fS_th2;
            fS_k1 = data.fS_k1;     fS_k2 = data.fS_k2; 
            f2_E1 = data.f2_E1;     f2_E2 = data.f2_E2; 
            Ism_E1 = data.Ism_E1;   Ism_E2 = data.Ism_E2;
            Ism_th1 = data.Ism_th1; Ism_th2 = data.Ism_th2; 
            Ism_k1 = data.Ism_k1;   Ism_k2 = data.Ism_k2; 
            d_E1 = data.d_E1;       d_E2 = data.d_E2; 
            dL_th1 = data.dL_th1;   dL_th2 = data.dL_th2; 
            dR_th1 = data.dR_th1;   dR_th2 = data.dR_th2;
            fS_X = data.fS_X;
            
            disp(['Saving new params to ',param_filename,'...']);
            save(param_filename, 'FL_E','bg_E1','bg_E2','bg_th1','bg_th2','bg_k1','bg_k2','f_E1','f_E2','f_th1','f_th2','f_k1','f_k2','fS_th1','fS_th2','fS_k1','fS_k2','f2_E1','f2_E2',...
                    'Ism_E1','Ism_E2','Ism_th1','Ism_th2','Ism_k1','Ism_k2','d_E1','d_E2','dL_th1','dL_th2','dR_th1','dR_th2','fS_X');
                           disp('Done.');
        end
        guidata(src,data)
    end 
    
    function importMatData(src, event, param)
        % Import new XY set .mat file
        data = guidata(src);
        data.(param) = get(src,'Value');
        disp(['Loading data from ',data.mat_filename,'...'])
        
        newdata = load( data.mat_filename );
        data.spec_list = newdata.spec_list; 
        data.E_list = newdata.E_list; 
        data.th_list = newdata.th_list;
        data.X_list = newdata.X_list;
        data.Y_list = newdata.Y_list;
        
        guidata(src,data);
        updateDisplaySpec(src, event, param);
        disp('Done.');
    end
    
    function importParams(src, event, param)
        % Import the parameters from associated param mat file
        data = guidata(src); 
        data.(param) = get(src,'Value');
        mat_filename = data.mat_filename; 
        param_filename = [mat_filename(1:end-8),'params.mat'];
        disp(['Loading params from ',param_filename,'...']);
        
        newdata = load( param_filename );

        varList = who('-file',param_filename);
        data.FL_E = newdata.FL_E;
        data.bg_E1 = newdata.bg_E1;     data.bg_E2 = newdata.bg_E2; 
        data.bg_th1 = newdata.bg_th1;   data.bg_th2 = newdata.bg_th2;
        data.f_E1 = newdata.f_E1;       data.f_E2 = newdata.f_E2;
        data.f_th1 = newdata.f_th1;     data.f_th2 = newdata.f_th2;        set(Ism_th1_ui,'String',num2str(data.Ism_th1)); set(Ism_th2_ui,'String',num2str(data.Ism_th2));

        data.fS_th1 = newdata.fS_th1;   data.fS_th2 = newdata.fS_th2;
        data.f2_E1 = newdata.f2_E1;     data.f2_E2 = newdata.f2_E2; 
        data.Ism_E1 = newdata.Ism_E1;   data.Ism_E2 = newdata.Ism_E2;
        data.Ism_th1 = newdata.Ism_th1; data.Ism_th2 = newdata.Ism_th2;
        data.fS_X = newdata.fS_X; 
        if ismember('d_E1',varList)
            data.d_E1 = newdata.d_E1;     data.d_E2 = newdata.d_E2; 
            data.dL_th1 = newdata.dL_th1; data.dL_th2 = newdata.dL_th2; 
            data.dR_th1 = newdata.dR_th1; data.dR_th2 = newdata.dR_th2;
        end
        if ismember('f_k1',varList)
            data.f_k1 = newdata.f_k1;   data.f_k2 = newdata.f_k2; 
            data.fS_k1 = newdata.fS_k1; data.fS_k2 = newdata.fS_k2; 
        end
        if ismember('bg_k1',varList)
            data.bg_k1 = newdata.bg_k1;   data.bg_k2 = newdata.bg_k2; 
            data.Ism_k1 = newdata.Ism_k1; data.Ism_k2 = newdata.Ism_k2; 
        end
        if ismember('maps',varList)
            data.feature_results = newdata.maps; data.feature_names = newdata.mapNames; 
        end
        
        % Then update the display panel
        set(FL_E_ui,'String',num2str(data.FL_E));
        set(bg_E1_ui,'String',num2str(data.bg_E1));     set(bg_E2_ui,'String',num2str(data.bg_E2));
        set(bg_th1_ui,'String',num2str(data.bg_th1));   set(bg_th2_ui,'String',num2str(data.bg_th2));
        set(bg_k1_ui,'String',num2str(data.bg_k1));     set(bg_k2_ui,'String',num2str(data.bg_k2));
        set(f_E1_ui,'String',num2str(data.f_E1));       set(f_E2_ui,'String',num2str(data.f_E2));
        set(f_th1_ui,'String',num2str(data.f_th1));     set(f_th2_ui,'String',num2str(data.f_th2));
        set(f_k1_ui,'String',num2str(data.f_k1));       set(f_k2_ui,'String',num2str(data.f_k2));
        set(fS_th1_ui,'String',num2str(data.fS_th1));   set(fS_th2_ui,'String',num2str(data.fS_th2));
        set(fS_k1_ui,'String',num2str(data.fS_k1));     set(fS_k2_ui,'String',num2str(data.fS_k2));
        set(f2_E1_ui,'String',num2str(data.f2_E1));     set(f2_E2_ui,'String',num2str(data.f2_E2));
        set(Ism_E1_ui,'String',num2str(data.Ism_E1));   set(Ism_E2_ui,'String',num2str(data.Ism_E2));
        set(Ism_th1_ui,'String',num2str(data.Ism_th1)); set(Ism_th2_ui,'String',num2str(data.Ism_th2));
        set(Ism_k1_ui,'String',num2str(data.Ism_k1));   set(Ism_k2_ui,'String',num2str(data.Ism_k2));
        set(X_step_ui,'String',num2str(data.X_step));   set(Y_step_ui,'String',num2str(data.Y_step));
        set(d_E1_ui, 'String',num2str(data.d_E1));      set(d_E2_ui,'String',num2str(data.d_E2));
        set(dL_th1_ui,'String',num2str(data.dL_th1));   set(dL_th2_ui,'String',num2str(data.dL_th2));
        set(dR_th1_ui,'String',num2str(data.dR_th1));   set(dR_th2_ui,'String',num2str(data.dR_th2));        
        set(fS_X_ui,'String',num2str(data.fS_X)); 

        
        guidata(src,data)
        updateDisplaySpec(src, event, param)
        disp('Done.');
    end
    
    
    function replotMaps(src,event,param)
        data = guidata(src);

        [mapGrids,X,Y] = display_feature_maps( data.mat_filename, data.FL_E, data.X_list, data.Y_list, data.X_step, data.Y_step, data.feature_results, data.feature_names );
        
                
        data.mapGrids = mapGrids; data.X = X; data.Y = Y;
        guidata(src, data); 
    end
    
    function calculateMapError(src,event,param)
        data = guidata(src);
        
        disp('Calculating map error...');
        [mapVariances] = calculate_map_error( data.mapGrids, data.feature_names );
        
        data.mapVariances = mapVariances;
        if contains(data.mat_filename, 'Eu05')
            assignin('base','mapVars_Eu',mapVariances); 
        elseif contains(data.mat_filename, 'Ce30')
            assignin('base','mapVars_Ce',mapVariances); 
        else 
            assignin('base','mapVars',mapVariances);
        end
        guidata(src,data); disp('Done.');
    end
        

    function removeTheMesh(src,event,param)
        data = guidata(src);
        mat_filename = data.mat_filename;
        fS_X = data.fS_X;
        % Check if mesh has been removed before. Always use the original raw data for mesh removal. 
        variableList = who('-file',mat_filename);
        if ismember('spec_list_withMesh',variableList)
            theVar = load(mat_filename, 'spec_list_withMesh');
            spec_list_withMesh = theVar.spec_list_withMesh;
        else % If it hasn't been done yet, spec_list itself is the raw data.
            theVar = load(mat_filename, 'spec_list');
            spec_list_withMesh = theVar.spec_list;
        end
        % Remove the mesh 
        disp(['Removing mesh from ',num2str(size(spec_list_withMesh,3)),' spectra...']);
        spec_list_noMesh = zeros(size(spec_list_withMesh));
        threshVals = [2e-5; 1e-5; 1e-5];
        for i=1:size(spec_list_withMesh,3)
            [spec,~] = removeMesh( spec_list_withMesh(:,:,i), threshVals );
            spec_list_noMesh(:,:,i) = mat2gray(spec);
        end
        spec_list = spec_list_noMesh; % Update spec_list with the mesh-removed spectra.
        disp(['Saving mesh-removed spectra as new spec_list...']);
        save(mat_filename, 'spec_list', 'spec_list_withMesh','fS_X', '-append');   disp('Done.');
        data.spec_list = spec_list_noMesh;
        
        guidata(src,data);
        updateDisplaySpec(src)     
    end

    function gotFoldername(src, event, param)
        data = guidata(src);
        data.(param) = get(src, 'String');
        disp(['In itx folder: ', data.(param) ])
        guidata(src, data);
    end

    function convertItx2Mat(src,event, param)
        % Convert folder of itx files into 3d array stored as .mat file
        data = guidata(src);
        data.(param) = get(src, 'Value');
        if data.convert_itx == 1
            
            itx_foldername = data.itx_foldername;
            idx = find(itx_foldername=='/'); idx(idx==length(itx_foldername))=[]; idx=idx(end);
            csv_filename = itx_foldername; csv_filename(end+1:end+4) = '    ';
            csv_filename(idx+1:idx+5) = 'meta_'; csv_filename(idx+6:end-4) = itx_foldername(idx+1:end-5); csv_filename(end-3:end) = '.csv';
                name_col_idx = 2; X_col_idx = 6;  Y_col_idx = 7;
                
            disp(['Converting itx from folder: ', data.itx_foldername ]);
            disp([' getting metadata from ', csv_filename ]);
            
            [ spec_list, E_list, th_list, X_list, Y_list ] = itx_csv_to_spec( itx_foldername, data.convert_itx, csv_filename, name_col_idx, X_col_idx, Y_col_idx );
            data.spec_list = spec_list; 
            data.E_list = E_list; 
            data.th_list = th_list; 
            data.X_list = X_list;
            data.Y_list = Y_list;
            
            % assignin('base','spec_list',spec_list)
            % assignin('base','E_list',E_list)
            % assignin('base','th_list',th_list)
            % assignin('base','X_list',X_list)
            % assignin('base','Y_list',Y_list)
            % save('your_mat_filename.mat','spec_list','E_list','th_list','X_list','Y_list')
             
        end
    end

    function updateParameterValue(src, event, param)
    % Update the gui data with updated params. Need str2num since using "edit" uicontrol which defaults to string inputs
        data = guidata(src);
        data.(param) = str2num(get(src,'String'));
        
        guidata(src,data);
        updateDisplaySpec(src)
    end

    function gotString2ValueSelection(src, event, param)
        % Update value after converting from string
        data = guidata(src);
        data.(param) = str2double(get(src,'String'));
        guidata(src,data);
    end
    
    function gotValueSelection(src, event, param)
        % Update logic-value selections
        data = guidata(src);
        data.(param) = get(src,'Value');
        guidata(src,data);
    end

    function gotMenuSelection(src, event, param)
        % Take popupmenu selection, get the corresponding string and update gui data
        data = guidata(src);
        selectVal = get(src,'Value');
        selectStr = get(src,'String');
        data.(param) = selectStr{selectVal};
        guidata(src,data);
    end

       


end