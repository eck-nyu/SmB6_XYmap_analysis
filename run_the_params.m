function [ feature_results, feature_names, k_list ] = run_the_params( data )
    tic;    
    %%%%%%%%%%%%%% Get variables from the data structure %%%%%%%%%%%%%%%%%%%%%%%%%
    mat_filename = data.mat_filename; 
    spec_list = data.spec_list;   E_list = data.E_list;  th_list = data.th_list; 
    X_list = data.X_list;         Y_list = data.Y_list;  
    FL_E = data.FL_E;             bg_E1 = data.bg_E1;     bg_E2 = data.bg_E2; 
    bg_th1 = data.bg_th1;         bg_th2 = data.bg_th2;  
    bg_k1 = data.bg_k1;           bg_k2 = data.bg_k2; 
    f_E1 = data.f_E1;             f_E2 = data.f_E2;   
    f_th1 = data.f_th1;           f_th2 = data.f_th2; 
    f_k1 = data.f_k1;             f_k2 = data.f_k2;
    fS_th1 = data.fS_th1;         fS_th2 = data.fS_th2;   
    fS_k1 = data.fS_k1;           fS_k2 = data.fS_k2;  
    fS_X = data.fS_X; 
    f2_E1 = data.f2_E1;           f2_E2 = data.f2_E2;   
    Ism_E1 = data.Ism_E1;         Ism_E2 = data.Ism_E2; 
    Ism_th1 = data.Ism_th1;       Ism_th2 = data.Ism_th2;  
    Ism_k1 = data.Ism_k1;         Ism_k2 = data.Ism_k2; 
    d_E1 = data.d_E1;             d_E2 = data.d_E2;  
    dL_th1 = data.dL_th1;         dL_th2 = data.dL_th2; 
    dR_th1 = data.dR_th1;         dR_th2 = data.dR_th2;  
    show_samples = data.show_samples;  
    azi0 = data.azi0;             pol0 = data.pol0;        hv = data.hv; 
    WF = data.WF;                 EB = data.EB;            innerPot = data.V0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    num_spec = size(spec_list, 3);
    Eax = E_list(:,1);    Estep = mean(abs(Eax(2:end)-Eax(1:end-1)));
    thax = th_list(:,1);  thstep = mean(abs(thax(2:end)-thax(1:end-1)));
    
    % E, th window size to search for min/max pixel intensities for bgnd/Ism
    Ewindow = round( 0.1 / Estep ); 
    thwindow = round( 3 / thstep ); 
    

    % Initialize arrays for feature maps
    tInts = NaN*ones(num_spec,1);
    tilts = NaN*ones(num_spec,1);
    f_Es = NaN*ones(num_spec,1);
    f_Ss = NaN*ones(num_spec,1);
    f2_Es = NaN*ones(num_spec,1);
    f2_Ss = NaN*ones(num_spec,1);
    f21_rats = NaN*ones(num_spec,1);
    Ibgs = NaN*ones(num_spec,1);
    
%     Isms = NaN*ones(num_spec,1);
    Ism_rat = NaN*ones(num_spec,1);
    Ism_abs = NaN*ones(num_spec,1);
    Ism_diff = NaN*ones(num_spec,1);
    Ism_ratInc = NaN*ones(num_spec,1);
    
    Ism_Es = NaN*ones(num_spec,1);
    
    k_list = NaN*ones(size(th_list,1),num_spec);
    
    for i = 1 : num_spec
        spec = mat2gray( spec_list(:,:,i) );
                
        if ~isnan(sum(sum(spec)))
            % Total intensity
            tInts(i) = sum(sum(spec));  
       
            %% Tilt from d-band peaks 
            [tilt0_idx, ~] = d_band_center( spec, Eax, thax, d_E1, d_E2, dL_th1, dL_th2, dR_th1, dR_th2);
            try
                tilts(i) = -thax(tilt0_idx); 
            end
            
            tParams = [azi0, pol0, hv, WF, EB, innerPot]; % azi0, pol0, hv, WF, EB, innerPot
            kax = theta2kx( thax, tParams(1), tilts(i), tParams(2), tParams(3), tParams(4), tParams(5), tParams(6) );
            k_list(:,i) = kax;
                       
            %% 4f energies, widths, amplitude ratio
            top_X = 0.1; interpDiv = 100; 
            
            % Use either raw-data theta ranges or converted-k axis ranges for taking EDCs for calculating f pk/widths 
            if f_k2 - f_k1 > 0
                [~, ~, f_E, f_S, f2_E, f2_S, f21_ratio] = f_com_finder( spec, Eax, kax, f_E1, f_E2, f2_E1, f2_E2, f_k1, f_k2, fS_k1, fS_k2, top_X, fS_X, interpDiv); 
            else           
                [~, ~, f_E, f_S, f2_E, f2_S, f21_ratio] = f_com_finder( spec, Eax, thax, f_E1, f_E2, f2_E1, f2_E2, f_th1, f_th2, fS_th1, fS_th2, top_X, fS_X, interpDiv); 
            end
            
            f_Es(i) = f_E;      
            f_Ss(i) = f_S;
            
            f2_Es(i) = f2_E;
            f2_Ss(i) = f2_S;
            
            f21_rats(i) = f21_ratio;
            

            %% Bgnd intensity
            % Convert to pixel idx 
            bg_E_srch = find(Eax > bg_E1, 1,'first') : find(Eax <= bg_E2, 1,'last');
            if bg_k2 - bg_k1 > 0
                bg_k_srch = find(kax > bg_k1,1,'first') : find(kax <= bg_k2, 1,'last'); 
                Ibg_window = spec(bg_E_srch, bg_k_srch); 
            else % If range not given in k, use raw-theta values 
                bg_th_srch = find(thax > bg_th1, 1,'first') : find(thax <= bg_th2, 1,'last');
                Ibg_window = spec( bg_E_srch, bg_th_srch );   
            end
            Ibg = mean(Ibg_window(:));
            Ibgs(i) = Ibg;
            
            %% Ism
            % Convert search-window ranges to pixel idx
            Ism_E_srch = find(Eax > Ism_E1, 1,'first') : find(Eax <= Ism_E2, 1,'last');           
            if Ism_k2 - Ism_k1 > 0
                Ism_k_srch = find(kax > Ism_k1,1,'first') : find(kax <= Ism_k2, 1,'last');
                Ism_window = spec(Ism_E_srch, Ism_k_srch);
            else
                Ism_th_srch = find(thax > Ism_th1, 1,'first') : find(thax <= Ism_th2, 1,'last');
                Ism_window = spec( Ism_E_srch, Ism_th_srch );  
            end
            % One method of evaluating Ism:
            Ism = mean(Ism_window(:)); 
                       
            % Another method of evaluating Ism: 
            Ism_edc = sum(spec,2);
            [~, Ism_Eidx] = max( Ism_edc(Ism_E_srch) );
%             [Ism, Ism_Eidx] = max( Ism_edc(Ism_E_srch) );
            Ism_Eidx = Ism_Eidx + Ism_E_srch(1) - 1;
            Ism_Es(i) = Eax(Ism_Eidx);
%             Ibg = min( Ism_edc(bg_E_srch) ); 

            
            Ism_abs(i) = Ism;
            Ism_diff(i) = Ism-Ibg;
            Ism_rat(i) = Ism / Ibg; 
            Ism_ratInc(i) = (Ism-Ibg) / Ism; 
            
        end

    end
    feature_names = {'tot int','tilt','f E','f S','f2_E','f2 S','f21_rat', 'Ism abs', 'Ism diff', 'Ism rat', 'Ism ratInc', 'Ism E'};
    feature_results = [ tInts, tilts, f_Es, f_Ss, f2_Es, f2_Ss, f21_rats,   Ism_abs,   Ism_diff,  Ism_rat, Ism_ratInc , Ism_Es];
    
    toc
    
    %% Plot figure showing input parameters 
    if show_samples > 0

        display_sample_spectra(mat_filename, show_samples, spec_list, thax, Eax, ...
                FL_E, f_Es, f2_Es, bg_E1,bg_E2, bg_th1,bg_th2, Ism_E1,Ism_E2, Ism_th1,Ism_th2, f_E1,f_E2, f_th1,f_th2, f2_E1,f2_E2 )
    pause(.1)
    end
    
    %% Plot feature maps
%     display_feature_maps( mat_filename, FL_E, X_list, Y_list, X_step, Y_step, feature_results, feature_names );
   
end



