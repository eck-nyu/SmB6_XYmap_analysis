function com_idx = mdc_com_finder( mdc, idx, ht, topX )

    topX_start = idx + 1 - find(fliplr(mdc(1:idx)) <= ht*(1-topX), 1,'first');
    topX_end = idx - 1 + find(mdc(idx:end) <= ht*(1-topX), 1,'first');
    topX_idx = topX_start : topX_end;

    com_idx = round( nansum(topX_idx .* mdc(topX_idx)) / nansum(mdc(topX_idx)) ); 

end