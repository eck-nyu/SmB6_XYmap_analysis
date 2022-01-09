function [X_list, Y_list] = csv_to_XY_list( csv_filename, fileList, fileNum, name_col_idx, X_col_idx, Y_col_idx )
% Input the location of csv file, and the corresponding column index of the
% file/scan no., X-positions, and Y-positions.
% Historically from ESM, indices are 2,5,6 respectively

if isfile(csv_filename)
    csv_data = readtable(csv_filename);
    
    X_list = zeros(fileNum,1);
    Y_list = zeros(fileNum,1);
    
    varNames = csv_data.Properties.VariableNames;
    csv_filename_list = eval(['csv_data.',varNames{name_col_idx},';']);
    csv_X_list = eval(['csv_data.',varNames{X_col_idx},';']);
    csv_Y_list = eval(['csv_data.',varNames{Y_col_idx},';']);
    
    for file_i = 1:fileNum
        filename = fileList(file_i).name;
        filename = erase(filename, '.itx');
        % Sometimes csv doesn't list in order of filename number, so
        % manually find the row with the filename for each one. 
        csv_idx = find( contains( csv_filename_list, filename ) );
        X_list(file_i) = csv_X_list( csv_idx );
        Y_list(file_i) = csv_Y_list( csv_idx ); 
    end
    disp(['Loaded XY coordinates from ',csv_filename]);
else
    disp('csv file not found, assuming XY coordinates.');
    sz = ceil(sqrt(fileNum));
    X_list = repmat( (1:sz)', sz,1);
    Y_list = repelem( (1:sz)',sz);
    
    X_list = X_list(1:fileNum);
    Y_list = Y_list(1:fileNum);
end


end