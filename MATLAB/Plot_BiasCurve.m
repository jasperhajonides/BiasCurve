function [y_all, theta] = Plot_BiasCurve(biases_all,angmus_all,config);

%%(c) Nick Myers 2018
addpath(genpath('/Users/jasperhvdm/Documents/Scripts/colourmap_BREWERMAP'));

figure; set(gcf,'color','white');
set(gcf,'Position',[10 10 300 300*max(config.same_diff_subplot)]);
hold on

areas = [];
psubs = config.psubs;
xlims = [-pi +pi]*0.5;
ylims = [-1 1]*.075;

% colours
[colors{1}] = colorcube(sum(config.same_diff_subplot==1));
[colors{2}] = flag(sum(config.same_diff_subplot==2)+1);
[colors{3}] = hsv(sum(config.same_diff_subplot==3)+1);
pcolors = [colors{1}(1:end,:); colors{2}(1:end,:) ; colors{3}];


%% format data

% if there are multiple conditions plotted we loop over these.
for icond = 1:length(config.same_diff_subplot)
    
    biases = biases_all(:,:,icond);
    angmus = angmus_all(:,:,icond);
    
    %% get area under the bias curve
    % average bias within the predefined window.
    y    = biases;
    y    = bsxfun(@minus,y,mean(y,2));
    area = [];
    nsub = size(y,1);
    
    for isub = 1:nsub
        theta         = angmus(isub,:);
        [theta,isort] = sort(theta); % sort all target-distractor distances.
        y(isub,:) = y(isub,isort);
        angmus(isub,:)   = theta;
        
        itheta = theta > 0; jtheta = theta < 0;
        area(isub,1) = squeeze(trapz(theta( itheta),y(isub, itheta),2)  ...
            - trapz(theta(~itheta),y(isub,~itheta),2));
    end;
    
    y_all(:,:,icond) = y;
    theta = mean(angmus,1);
    %% Plot Bias Curves
    % plot the biases using the 
    shading = 0.275;    
    
    %which plot
    subplot(max(config.same_diff_subplot),1,config.same_diff_subplot(icond))
    condnames    = {'Bias'};
    pcolor = pcolors(icond,:);
    
    hold on
    plotdat   = squeeze(y(psubs,:))/2;
    npts      = size(plotdat,2);
    
    theta     = circ_mean(angmus(psubs,:),[],1)/2;
    
    plot(xlims,[0 0],'k--')
    plot([0 0],ylims,'k--')
    
    [h parea c sarea] = ttest(area(psubs));
    
    [h p c s] = ttest(mean(plotdat,3));
    yheight   = min(ylims) + 0.50*range(ylims);
    pthresh   = 0.05;
    b         = bwconncomp(p < pthresh);
    pareathresh = 0.05;
    if parea < pareathresh
        for ib = 1:b.NumObjects
            icur = b.PixelIdxList{ib};
            %equal spacing significance markings
            max_val = sum(config.same_diff_subplot==config.same_diff_subplot(icond));
            sp = zeros(length(config.same_diff_subplot),1);
            if max_val == 1
                spacing = 0;
                sp(config.same_diff_subplot==config.same_diff_subplot(icond)) = spacing;
            elseif max_val == 2
                spacing = [ -0.003 0.003];
                sp(config.same_diff_subplot==config.same_diff_subplot(icond)) = spacing;
             elseif max_val == 3
                spacing = [ -0.003 0 0.003];
                sp(config.same_diff_subplot==config.same_diff_subplot(icond)) = spacing;
            elseif max_val == 4
                spacing = [ -0.003 -0.0015 0.0015 0.003];
                sp(config.same_diff_subplot==config.same_diff_subplot(icond)) = spacing;
            end
            y_val = ones(1,nnz(icur))*yheight+ sp(icond);
            plot(theta(icur),y_val,'color',pcolor,'LineStyle','-','linewidth',4)
        end
    end
    
    cfg = [];
    cfg.patchcolor = pcolor;
    cfg.alpha = 0.5;
    cfg.shading = true;
    [~,hpatch] = plotpatch(mean(plotdat,3),theta,pcolor, cfg);
    
    ylim(ylims) , xlim(xlims)
    set(gca,'ytick',unique([0 ylims]),'xtick',unique([0 xlims]),'xticklabel',{'-1/2 pi' '' '+1/2 pi'})
    ylabel(sprintf('Bias'),'fontsize',14,'fontname','helvetica'); 
    set(gca, 'Fontsize', 14);
    
end


if config.do_print
    set(gcf, 'PaperPositionMode', 'auto');
    figfolder = sprintf('%s/figure',cd); if ~exist(figfolder,'dir'), mkdir(figfolder); end
    figname   = sprintf('%s/BiasCurve_Avg',figfolder);
    print(gcf,figname,'-painters','-depsc2'); print(gcf,figname,'-dpng')
end

