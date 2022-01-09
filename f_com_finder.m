function [edcq, eaxq, f_E, f_S, f2_E, f2_S, f21_ratio] = f_com_finder( spec, Eax, thax, ...
                                f_E1, f_E2, f2_E1, f2_E2, f_th1, f_th2, fS_th1, fS_th2, topX, sigX, interpDiv)
    % Calculates f-pk as the center of mass of highest topX % of edc peak
    % Sets width of f-pk as full-width-sigX-max 
    
    % Take two EDCs, one using th range for finding f Es, another for calculating f widths 
    edc = sum(spec(:,find(thax>f_th1,1,'first') : find(thax<=f_th2,1,'last')),2);
    edcS = sum(spec(:,find(thax>fS_th1,1,'first') : find(thax<=fS_th2,1,'last')),2);
    
    % Interp the energy axis, use this to interp EDCs
    eaxq = linspace(Eax(1),Eax(end), numel(Eax)*interpDiv);
    
    edcq = interp1(Eax, edc, eaxq, 'spline');
    edcqS = interp1(Eax, edcS, eaxq, 'spline');
    
    % Obtain E windows within to search for f,f2 
    f_rng_Eq = find(eaxq > f_E1, 1,'first') : find(eaxq <= f_E2, 1,'last');
    f2_rng_Eq = find(eaxq > f2_E1, 1,'first') : find(eaxq <= f2_E2, 1,'last');
    
    % Get cropped EDCs only within these windows for each f,f2 and their width EDCs
    edcq_f = edcq(f_rng_Eq); 
    edcq_fS = edcqS(f_rng_Eq); 
    
    edcq_f2 = edcq(f2_rng_Eq);
    edcq_f2S = edcqS(f2_rng_Eq);
    
    % Find the pk positions corresponding to f,fS, f2,f2S
    [fpk_ht, fpk_idx] = max(edcq_f);
    if numel(fpk_idx) > 1, fpk_idx = round(mean(fpk_idx)); end
    fpk_idx = f_rng_Eq(fpk_idx);
    
    [f2pk_ht, f2pk_idx] = max(edcq_f2); 
    if numel(f2pk_idx) > 1, f2pk_idx = round(mean(f2pk_idx)); end
    f2pk_idx = f2_rng_Eq(f2pk_idx);
    
    [fSpk_ht, fSpk_idx] = max(edcq_fS);
    if numel(fSpk_idx) > 1, fSpk_idx = round(mean(fSpk_idx)); end
    fSpk_idx = f_rng_Eq(fSpk_idx);
    
    [f2Spk_ht, f2Spk_idx] = max(edcq_f2S);
    if numel(f2Spk_idx) > 1, f2Spk_idx = round(mean(f2Spk_idx)); end
    f2Spk_idx = f2_rng_Eq(f2Spk_idx);
    
    % Get the top X-percent EDC peaks for f, f2 
    topXend = fpk_idx - 1 + find(edcq(fpk_idx:end) <= (1-topX) * fpk_ht, 1, 'first');
    topXstart = fpk_idx +1 - find(fliplr(edcq(1:fpk_idx)) <= (1-topX) * fpk_ht, 1, 'first');
    topXidx = topXstart : topXend;
    
    topXend2 = f2pk_idx - 1 + find(edcq(f2pk_idx:end) <= (1-topX) * f2pk_ht, 1, 'first');
    topXstart2 = f2pk_idx +1 - find(fliplr(edcq(1:f2pk_idx)) <= (1-topX) * f2pk_ht, 1, 'first');
    topXidx2 = topXstart2 : topXend2;
    
    % Calculate f,f2 E as center-of-mass of top X-percent peaks
    f_E = sum( eaxq(topXidx).*edcq(topXidx)) / sum(edcq(topXidx));   
    f_E_idx = round( sum(topXidx .* edcq(topXidx)) / sum(edcq(topXidx)) );
    
    f2_E = sum( eaxq(topXidx2).*edcq(topXidx2) ) / sum(edcq(topXidx2));
    f2_E_idx = round( sum(topXidx2 .* edcq(topXidx2)) / sum(edcq(topXidx2)) );
    
    % Get the surrounding X-fraction surrounding fS, f2S 
    sigXend = fSpk_idx - 1 + find(edcqS(fSpk_idx:end) <= (sigX * fSpk_ht), 1, 'first');
    sigXstart = fSpk_idx + 1 - find(fliplr(edcqS(1:fSpk_idx)) <= (sigX * fSpk_ht), 1, 'first');
    
    sigXend2 = f2Spk_idx - 1 + find(edcqS(f2Spk_idx:end) <= (sigX * f2Spk_ht), 1, 'first');
    sigXstart2 = f2Spk_idx + 1 - find(fliplr(edcqS(1:f2Spk_idx)) <= (sigX * f2Spk_ht), 1, 'first');

    
    % Calculate fS,f2S as full-width at X-maximum 
    f_S = eaxq(sigXend) - eaxq(sigXstart);   
    if isempty(f_S), f_S = NaN; end
    if sigXstart < f_rng_Eq(1) || sigXend > f_rng_Eq(end)
        sigXstart = NaN;
        sigXend = NaN; 
        f_S = NaN; 
    end
    
    f2_S = eaxq(sigXend2) - eaxq(sigXstart2);   
    if isempty(f2_S), f2_S = NaN; end
    if sigXstart2 < f2_rng_Eq(1) || sigXend2 > f2_rng_Eq(end), f2_S = NaN; end
    
    if sum(isempty([fpk_idx, f2pk_idx])) == 0
        f21_ratio = edcq(f2pk_idx) / edcq(fpk_idx); 
    else, f21_ratio = NaN;
    end
    
    % Plot everything
%     figure, plot(Eax, edc, 'k:'), hold on, 
%     plot(eaxq, edcq, 'k')
%     plot(eaxq(topXidx),edcq(topXidx),'r','LineWidth',2);
%     plot(eaxq(f_E_idx),edcq(f_E_idx),'r+');
%     plot(eaxq(fpk_idx),edcq(fpk_idx),'ro');
%     
%     plot(eaxq(topXidx2),edcq(topXidx2),'m','LineWidth',2);
%     plot(eaxq(f2_E_idx),edcq(f2_E_idx),'m+');
%     plot(eaxq(f2pk_idx),edcq(f2pk_idx),'mo');
%     
%     plot(Eax, edcS, 'k:'), hold on, 
%     plot(eaxq, edcqS, 'k')
%     plot(eaxq(sigXstart:sigXend),edcqS(sigXstart:sigXend),'r--');
%     plot(eaxq(fSpk_idx),edcqS(fSpk_idx),'rx');  
%     plot([eaxq(sigXstart),eaxq(sigXend)],[edcqS(sigXstart),edcqS(sigXend)],'b')
%     
%     plot(eaxq(sigXstart2:sigXend2),edcqS(sigXstart2:sigXend2),'m--');
%     plot(eaxq(f2Spk_idx),edcqS(f2Spk_idx),'mx');  
%     plot([eaxq(sigXstart2),eaxq(sigXend2)],[edcqS(sigXstart2),edcqS(sigXend2)],'c')
end