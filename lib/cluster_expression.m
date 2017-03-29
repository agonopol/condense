function cluster_expression(samples, clusterassigments)
    % Make expression matrices for EB data
    load nataliaBranchPoints

    % Plot the data

    % figure
    % scatter(Ynonmds(:,1),Ynonmds(:,2),5,'filled')
    % hold on
    % for m=1:length(branch)
    %     scatter(Ynonmds(branch{m},1),Ynonmds(branch{m},2),5,'filled')
    %     
    % end
    % 
    % hold off

    branch_labels=zeros(length(Ynonmds),1);
    cbranch=[];
    for m=1:length(branch)
        branch_labels(branch{m})=m;
        cbranch=[cbranch; m*ones(length(branch{m}),1)];
    end

    load nataliaMagic_fewer2_a10_k5_t25

    load nataliaMagic_fewer_a10_k5_t25

    data=[data_markers data_markers2];

    markers=[markers markers2];
    Y_max=max(data);
    Y_min=min(data);

    [N,d]=size(data);
    M=length(branch);

    mean_data=zeros(M,d);
    dremi=zeros(M,d);
    data_sorted=[];
    drevi_opt = struct('DrawPlot','false');

    z=zscore(data);
    bins=zeros(M,1);

    % cbranch=sort(branch_labels);
    % cbranch(cbranch==0)=[];
    for m=1:M
    %     Sorting along branches
        x1=Ynonmds(points(mainPoints(edges(m,1))),:);
        x2=Ynonmds(points(mainPoints(edges(m,2))),:);
        xnew=x2-x1;
        bpoints=branch{m};

        if m==M % make the last cluster bigger, just for this dataset
            Nrep=6;
            cbranch=[cbranch; M*ones((Nrep-1)*length(bpoints),1)];
            bpoints=repmat(bpoints,Nrep,1);

        end
        B=length(bpoints);
        Ytemp=Ynonmds(bpoints,:);
        Ytemp=bsxfun(@minus,Ytemp,x1);
        xdot=dot(Ytemp,repmat(xnew,B,1),2);
        [y,I]=sort(xdot);
        ztemp=z(bpoints(I),:);
        data_sorted=[data_sorted; ztemp];

    %     compute DREMI
        temp=data(bpoints(I),:);
        n=1:B;
        for l=1:d
            drevi_opt.MaxY=Y_max(l);
            drevi_opt.MinY=Y_min(l);
            dremi_obj=visualize_2D_DREVI(n', temp(:,l), drevi_opt);
            dremi(m,l)=dremi_obj.dremi_score;
        end
        mean_data(m,:)=mean(ztemp);
        bins(m)=B;



    end
    bins=cumsum(bins);

    D=pdist(data_sorted');  
    Z = linkage(D, 'average');
    ind_cols = optimalleaforder(Z, D);
    marker_names=[];
    for m=1:length(ind_cols)
       marker_names{m}=markers{ind_cols(m)}; 

    end

    % Make the figure
    figure
    imagesc(data_sorted(:,ind_cols)')
    zd=z(:);
    caxis([prctile(zd,5) prctile(zd,95)])

    set(gca,'ytick',1:38)
    set(gca,'ytickLabels',marker_names)
    set(gca,'xtick',[])
    set(gca,'Ticklength',[0 0])
    h=colorbar;
    set(h,'xtick',[])

    line([bins bins]', repmat(ylim, length(bins), 1)', 'color', 'k','Linewidth',2);


    % Cluster labels
    cmp=cbrewer('div','Spectral',M+1);
    cmp(1,:)=[.5 .5 .5];

    figure
    imagesc(cbranch')
    axis off
    % set(gcf,'paperposition',[0 0 8 0.5]);
    caxis([0 M])
    colormap(cmp)
    line([bins bins]', repmat(ylim, length(bins), 1)', 'color', 'k','Linewidth',2);

    % PHATE embedding
    figure
    scatter(Ynonmds(:,1),Ynonmds(:,2),5,branch_labels,'filled');
    colormap(cmp)

    % Dremi expression matrix
    figure
    D=pdist(dremi');  
    Z = linkage(D, 'average');
    ind_cols_dremi = optimalleaforder(Z, D);
    marker_names_dremi=[];
    for m=1:length(ind_cols_dremi)
       marker_names_dremi{m}=markers{ind_cols_dremi(m)}; 

    end

    imagesc(dremi(:,ind_cols_dremi))
    xticklabel_rotate(0:37,90,marker_names_dremi,'fontSize',12)

    % Mean expression matrix
    figure
    D=pdist(mean_data');  
    Z = linkage(D, 'average');
    ind_cols_mean = optimalleaforder(Z, D);
    marker_names_mean=[];
    for m=1:length(ind_cols_dremi)
       marker_names_mean{m}=markers{ind_cols_mean(m)}; 

    end

    imagesc(mean_data(:,ind_cols_mean))
    xticklabel_rotate(0:37,90,marker_names_mean,'fontSize',12)
end

