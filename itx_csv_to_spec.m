% itx_foldername = '/home/data/eck325/SmB6/BNL_SmB6_20210808_all_the_data/XYSmB6_itx/';
% csv_filename = '/home/data/eck325/SmB6/BNL_SmB6_20210808_all_the_data/meta_XYSmB6.csv';
%     name_col_idx = 2;  X_col_idx = 6;  Y_col_idx = 7;
% 
% [spec_list, E_list, th_list, X_list, Y_list] = itx_to_spec( itx_foldername, csv_filename, 2, 6, 7);

function [ spec_list, E_list, th_list, X_list, Y_list ] = itx_csv_to_spec( itx_foldername, convert_itx, csv_filename, name_col_idx, X_col_idx, Y_col_idx )
%     function [X_list, Y_list] = itx_to_spec( itx_foldername, csv_filename, name_col_idx, X_col_idx, Y_col_idx)
% Convert folder of pxt files to a folder of itx files using Igor.
% In Igor, open > file > and open Desktop/igor_ptx_to_itx.ipf 
% 
% Inputs: 
% folder name containing list of itx files of spectra. 
% csv_filename is name of metadata file. 
% name_col_idx is column idx of filenames
% X(Y)_col_idx is col idx of X(Y) coordinates. 

% Outputs: 
% spec_list is 3d matrix Edim x thdim x numSpec 
% E(th)_list is 2d matrix E(th)dim x numSpec
% X(Y)_list is 1d column of sample-X(Y) coordinates


fileList = dir(fullfile(itx_foldername, '*.itx'));
fileNum = numel(fileList);

%% Import the XY coordinates from csv file

[X_list, Y_list] = csv_to_XY_list( csv_filename, fileList, fileNum, name_col_idx, X_col_idx, Y_col_idx );
% if isfile(csv_filename)
%     csv_data = readtable(csv_filename);
%     
%     X_list = zeros(fileNum,1);
%     Y_list = zeros(fileNum,1);
%     
%     varNames = csv_data.Properties.VariableNames;
%     csv_filename_list = eval(['csv_data.',varNames{name_col_idx},';']);
%     csv_X_list = eval(['csv_data.',varNames{X_col_idx},';']);
%     csv_Y_list = eval(['csv_data.',varNames{Y_col_idx},';']);
%     
%     for file_i = 1:fileNum
%         filename = fileList(file_i).name;
%         filename = erase(filename, '.itx');
%         % Sometimes csv doesn't list in order of filename number, so
%         % manually find the row with the filename for each one. 
%         csv_idx = find( contains( csv_filename_list, filename ) );
%         X_list(file_i) = csv_X_list( csv_idx );
%         Y_list(file_i) = csv_Y_list( csv_idx ); 
%     end
%     disp(['Loaded XY coordinates from ',csv_filename]);
% else
%     disp('csv file not found, assuming XY coordinates.');
%     sz = ceil(sqrt(fileNum));
%     X_list = repmat( (1:sz)', sz,1);
%     Y_list = repelem( (1:sz)',sz);
%     
%     X_list = X_list(1:fileNum);
%     Y_list = Y_list(1:fileNum);
% end

%% Now import the spectra from the itx folder. 

if convert_itx == 1
    disp(['Loading files from ',itx_foldername,'...']);
    for file_i = 1:fileNum
        filename = fileList(file_i).name; % name is generic field used in fileList function
        [spec, Eaxis, angleAxis] = load_itx([itx_foldername, filename]);
        if file_i ==1
            spec_list = spec;
            E_list = reshape( Eaxis,[],1);
            th_list = reshape( angleAxis, [],1);
        else
            if size(spec,1) ~= size(spec_list,1) || size(spec,2) ~= size(spec_list,2)
                spec = NaN*ones(size(spec_list,1),size(spec_list,2));
                Eaxis = NaN*ones(size(E_list(:,end)));
                angleAxis = NaN*ones(size(th_list(:,end)));
            end
            spec_list = cat(3, spec_list, spec);
            E_list = cat(2, E_list, reshape( Eaxis, [],1));
            th_list = cat(2, th_list, reshape( angleAxis, [],1));
        end
    end
    disp(['Finished loading ',num2str(fileNum),' spectra.'])
end
% fig = figure;
% imagesc(th_list(:,1),E_list(:,1),nanmean(spec_list,3)); axis xy;
% xlabel('theta (deg)'); ylabel('Energy (eV)');
% title(itx_foldername,'Interpreter','None');

end