
function [tilt0_idx, dIdxs] = d_band_center( spec, Eax, thax, d_E1, d_E2, dL_th1, dL_th2, dR_th1, dR_th2 )

topX = 0.1; % Set d-band as com of topX percent surrounding peak mdc intensity 

mdc_Es = linspace(d_E1, d_E2, 25); 
mdc_rows = round( interp1( Eax, (1:numel(Eax)), mdc_Es)); 
mdc_width = round( .05 / abs(Eax(2)-Eax(1)));
mdc_conv_sigma = 0.5;

mdc_colRange = [round(interp1( thax, (1:numel(thax)), [dL_th1, dL_th2])); ...
                round(interp1( thax, (1:numel(thax)), [dR_th1, dR_th2]))];

mdc_spec = spec - 0.5*nanmean(spec(:));
mdc_spec(mdc_spec<0)=0;

dIdxs = NaN*ones(numel(mdc_rows),2);

for row_i =1:numel(mdc_rows)
    mdc_row = mdc_rows(row_i);
    mdc = sum(mdc_spec(mdc_row-mdc_width:mdc_row+mdc_width,:),1);
    mdc_conv = gaussmf( thax, [mdc_conv_sigma, 0] ); %Convolute w gauss width mdc_conv_sigma (centered @ 0deg)
    mdc = conv(mdc, mdc_conv, 'same');
    mdc = (mdc-min(mdc))/(max(mdc)-min(mdc));

    try
        [Lht,Lidx] = max(mdc(mdc_colRange(1,1):mdc_colRange(1,2)));
        [Rht,Ridx] = max(mdc(mdc_colRange(2,1):mdc_colRange(2,2)));

        Lidx = Lidx + mdc_colRange(1,1) - 1;
        Ridx = Ridx + mdc_colRange(2,1) - 1;

        Lcom_idx = mdc_com_finder( mdc, Lidx, Lht, topX );
        Rcom_idx = mdc_com_finder( mdc, Ridx, Rht, topX );

        dIdxs(row_i,:) = [Lcom_idx, Rcom_idx];
    end    

end

% NaN out dIdx where they fall out of the col windows specified 
out_idx = (dIdxs(:,1) < mdc_colRange(1,1) | dIdxs(:,1) > mdc_colRange(1,2) ...
                | dIdxs(:,2) < mdc_colRange(2,1) | dIdxs(:,2) > mdc_colRange(2,2) );

% Only keep pairs of dIdxs, throw away any with one leg missing 
nan_idx = any(isnan(dIdxs),2);
dIdxs(nan_idx,:) = NaN;
dIdxs(out_idx,:) = NaN; 

tilt0_idx = round(nanmean(dIdxs(:)));

% figure, imagesc(mdc_spec), axis xy; hold on, 
% plot( dIdxs, [mdc_rows', mdc_rows'],  'yx');
% plot(tilt0_idx*[1,1],[1,size(spec,1)],'w'); 
end
