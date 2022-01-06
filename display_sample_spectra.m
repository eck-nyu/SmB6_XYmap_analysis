function display_sample_spectra(mat_filename, show_num, spec_list, thax, Eax, ...
    FL_E, f_Es, f2_Es, bg_E1,bg_E2, bg_th1,bg_th2, Ism_E1,Ism_E2, Ism_th1,Ism_th2, f_E1,f_E2, f_th1,f_th2, f2_E1,f2_E2 )

fig = figure;
        
mn_spec = nanmean(spec_list,3);
[mn_spec,~,~] = removeMesh_template(mat2gray(mn_spec),[2e-5;1e-5],[]);
edc = sum(mn_spec,2); edc = (edc-min(edc))/(max(edc)-min(edc));
edc = edc * (max(thax)-min(thax)); edc = edc + min(thax);

pltSz = factor(show_num);
if numel(pltSz) == 1, pltSz = [1, pltSz]; 
elseif numel(pltSz) > 2
    idx = [1:numel(pltSz)];
    idx1 = round(numel(pltSz)/2) : round(numel(pltSz)/2)+1;
    idx2 = idx(idx~=idx1(1) & idx~=idx1(2));

    pltSz = [pltSz(idx1(1))*pltSz(idx1(2)), pltSz(idx2(1))*pltSz(idx2(2))];
end

for rand_smp = 1:show_num
    if rand_smp == 1

        ax = subplot(pltSz(1), pltSz(2), 1);
        imagesc(thax,Eax,mat2gray(mn_spec)); axis xy; hold on;
    else
        idx = round(1+(size(spec_list,3)-2)*rand()); 

        ax = subplot(pltSz(1), pltSz(2), rand_smp);
        [spec,~,~] = removeMesh_template(spec_list(:,:,idx),[2e-5;1e-5;1e-5],[]);
        spec = mat2gray(spec);

        imagesc(ax, thax,Eax,spec); set(ax,'YDir','normal'); hold(ax,'on'); 
    end

    plot(ax, [min(thax),max(thax)],FL_E*[1,1],'w:');
    plot(ax, edc, Eax, 'w');
    if rand_smp > 1
        f_Eidx = find(Eax >= f_Es(idx),1,'first');
        f2_Eidx = find(Eax >= f2_Es(idx),1,'first');

        plot(ax, edc(f_Eidx),Eax(f_Eidx), 'r+');
        plot(ax, edc(f2_Eidx), Eax(f2_Eidx), 'm+');
    end    

    rectangle(ax, 'Position', [bg_th1,bg_E1,bg_th2-bg_th1,bg_E2-bg_E1],'EdgeColor',[.5,.5,.5]);
    rectangle(ax, 'Position', [Ism_th1,Ism_E1,Ism_th2-Ism_th1,Ism_E2-Ism_E1],'EdgeColor',[.8,.8,.8]);
    rectangle(ax, 'Position', [f_th1,f_E1,f_th2-f_th1,f_E2-f_E1],'EdgeColor',[1,0,0]);
    rectangle(ax, 'Position', [f_th1,f2_E1,f_th2-f_th1,f2_E2-f2_E1],'EdgeColor',[1,.5,0]);
    
    caxis([0.05,0.7])
    xlim(ax, [min(thax),max(thax)]); ylim(ax,[min(Eax),max(Eax)]);
    if rand_smp == 1, title('mean','FontSize',8);
    else, title(['idx=',num2str(idx)],'FontSize',8,'Interpreter','None');
    end
    pause(.1); 
end
sgtitle(mat_filename,'Interpreter','None', 'FontSize',8);
end