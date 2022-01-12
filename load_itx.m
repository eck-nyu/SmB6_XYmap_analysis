function [ spec, Eaxis,angleAxis ] = load_itx( fileName )
%LOAD_ITX Summary of this function goes here
%   Detailed explanation goes here
%   laod itx format data 
    
    fdata=importdata(fileName);
    
    spec=fdata.data;
    
    % now load the energy and angle axis
    
    f1=fopen(fileName);
    
    for i=1:2
        fline=fgetl(f1);
        
        % now import data size
        t_size=cell2mat(textscan(fline,'WAVES/D/N=(%f, %f)'));
    end
    
    while ~strcmp(fline,'END')
        fline=fgetl(f1);
    end
    % now read in the next line contain the energy and angle axies
    fline=fgetl(f1);
    
    % split the line
    line_string=split(fline(2:end),';');
    xinfo=cell2mat(textscan(line_string{1},' SetScale/P x %f, %f'));
    yinfo=cell2mat(textscan(line_string{2},' SetScale/P y %f, %f'));
    
    
    
    Eaxis=xinfo(1)+[0:t_size(1)-1]*xinfo(2);
    angleAxis=yinfo(1)+[0:t_size(2)-1]*yinfo(2);

    spec=spec(:,~all(isnan(spec)));
    % close file
    fclose(f1);
    
end

