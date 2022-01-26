function plot_qti_invariants(s)
% function to plot QTI invariants for different fits/protocols
%

n_fits = max(size(s));

% by default, we take the middle slice to be shown
slice = round(size(s{1}.s0,3)/2);

% open a figure
figure()

% for each fit (invariant structure), we plot the QTI main invariants

if license('test', 'image_toolbox')
    % Do some rescaling BASED ON DIB Data results!!
    s0_max = 1500;
    md_max = 3.2;
    fa_max = 1;
    cm_max = 1;
    cmu_max = 1;
    ufa_max = 1;
    cmd_max = 0.25;
    cc_max = 1;
    tmp = [];
    for i = 1:n_fits
        
        % for each fit (invariant structure), we plot the QTI main invariants
        s_tmp = s{i};
        s0 =  (fliplr(s_tmp.s0(:,:,slice)))' ./ s0_max;
        md =  (fliplr(s_tmp.MD(:,:,slice)))' ./ md_max;
        fa =  (fliplr(s_tmp.FA(:,:,slice)))' ./ fa_max;
        cm =  (fliplr(s_tmp.C_M(:,:,slice)))' ./ cm_max;
        cmu = (fliplr(s_tmp.C_mu(:,:,slice)))' ./ cmu_max;
        ufa = (fliplr(s_tmp.uFA(:,:,slice)))' ./ ufa_max;
        cmd = (fliplr(s_tmp.C_MD(:,:,slice)))'./ cmd_max;
        cc =  (fliplr(s_tmp.C_c(:,:,slice)))' ./ cc_max;
        tmp = cat(3,tmp,s0,md,fa,cm,cmu,ufa,cmd,cc);
        
    end
    montage(tmp, 'size', [n_fits, 8])
    
else % use subplot
    
    params = ["$S_0$", "$MD$", "$FA$", "$C_M$", "$C_{\mu}$", "$\mu FA$", "$C_{MD}$", "$C_c$"];
    
    for i = 1:n_fits
        
        % for each fit (invariant structure), we plot the QTI main invariants
        s_tmp = s{i};
        s0 =  (fliplr(s_tmp.s0(:,:,slice)))';
        md =  (fliplr(s_tmp.MD(:,:,slice)))';
        fa =  (fliplr(s_tmp.FA(:,:,slice)))';
        cm =  (fliplr(s_tmp.C_M(:,:,slice)))';
        cmu = (fliplr(s_tmp.C_mu(:,:,slice)))';
        ufa = (fliplr(s_tmp.uFA(:,:,slice)))';
        cmd = (fliplr(s_tmp.C_MD(:,:,slice)))';
        cc =  (fliplr(s_tmp.C_c(:,:,slice)))';
        tmp = cat(3,s0,md,fa,cm,cmu,ufa,cmd,cc);
        
        for j = 1:8
            
            subplot(n_fits,8, j+8*(i-1))
            imagesc(imresize(tmp(:,:,j),4)), colormap gray, axis square
            ax = gca;
            ax.YTick = [];
            ax.YTickLabel = [];
            ax.XTick = [];
            ax.XTickLabel = [];
            ax.FontSize = 20;
            ax.Title.Interpreter = 'latex';
            
            % adjust scaling
            switch j
                case 1
                    ax.CLim = [0 1200];
                case 2
                    ax.CLim = [0 3.2];
                case 7
                    ax.CLim = [0 0.3];
                otherwise
                    ax.CLim = [0 1];
            end
            
            % set title
            if i == 1
                
                ax.Title.String = params(j);
                
            end
            
        end
    end
end

end