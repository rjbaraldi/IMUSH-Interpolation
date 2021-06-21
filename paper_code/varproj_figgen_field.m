function varproj_figgen_field(X,tstring, form, params, savefile)
        if contains(savefile, 'difference')
            ca = [0, 1.2]; 
        elseif form==3
            ca = [-.2, 1.2];
        elseif form==2
            ca = [-.5, 1.2]; 
        end
        ax = params.Gx(1):params.resolution:params.Gx(2); 
        if form==2
            [n,m] = size(X); 
            if n==m
                X = imrotate(X,90);
            end
            figure; 
            imagesc(flipud(X));
            caxis(ca)
        elseif form==4
            [n,m] = size(X); 
            if n==m
                X = imrotate(X,90);
            end
            figure; %subsampled tt's
            X(X~=0) = 1; 
            imagesc(flipud(X));colormap gray; caxis([0 1]*1);
        elseif form==3
            X = imrotate(X(1:params.Rx, 1:params.Ry), 90); 
            figure; %subsampled tt's
            imagesc(ax, ax, flipud(X));
            caxis(ca)
        elseif form==1
            figure;
            imagesc(flipud(X));
            caxis(ca)
        end
         title(tstring, 'Fontsize', 22, 'Fontweight', 'Bold');colorbar;
        if form==1
           ylabel('(R_x,R_y)', 'Fontsize', 22, 'Fontweight', 'Bold')
           xlabel('S_xS_y', 'Fontsize', 22, 'Fontweight', 'Bold')
        elseif form==2 || (form==4 && n==m)
            xlabel('(R_x, S_x)', 'Fontsize', 22, 'Fontweight', 'Bold');
            ylabel('(R_y, S_y)', 'Fontsize', 22, 'Fontweight', 'Bold');
        elseif form==4 && n~=m
            ylabel('(R_x, R_y)', 'Fontsize', 22, 'Fontweight', 'Bold');
            xlabel('(S_{x}S_{y})', 'Fontsize', 22, 'Fontweight', 'Bold');
        elseif form==3
            xlabel('R_x', 'Fontsize', 22, 'Fontweight', 'Bold');
            ylabel('R_y', 'Fontsize', 22, 'Fontweight', 'Bold');
        end
%         set(gca, 'Fontsize', 22, 'Fontweight', 'Bold','xtick',[], 'xticklabel', [], 'ytick',[], 'yticklabel', [],'YDir','normal');
        set(gca, 'Fontsize', 22, 'Fontweight', 'Bold','YDir','normal');
        path = '../figs/';
        saveas(gcf,strcat(path, savefile, '.jpg'));
        close;
end