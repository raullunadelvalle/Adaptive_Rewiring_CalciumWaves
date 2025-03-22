clear all
close all
clc


% Use this script to generate figures similar to those in Figure 7 and Figure 8 of:
%Luna, R., Li, J., Bauer, R., & van Leeuwen, C. (2023). Retinal waves in adaptive 
%rewiring networks orchestrate convergence and divergence in the visual system.
%Network Neuroscience

%Figures 1 to 9 correspond to Figure 7 in the article
%Figures 10 to 14 correspond to Figure 8 in the article



neighbours = 0; 
%0: neighbours are nodes targeted by initiator's outlinks 
%1: neighbours are nodes targeted by initiator's outlinks and nodes targeting initiators

%Here the neighbourhood of a ganglion cell are only nodes targeted by the
%out-links of the ganglion cell


cd '..\scripts\Output'

%% load evolved networks
fid = py.open('boost_nor_MatPlusCoords.pckl','rb'); 
A = py.pickle.load(fid);
coords = A(1);
adj = A(2);

num_nodes = size(double(coords{1}{1}),1);
node_list = [1:1:num_nodes];

fid=py.open('simulation_params.pckl','rb'); 
params = py.pickle.load(fid);

p_dist = params(3);
p_dist = cell2mat(cell(p_dist{1}));

p_boost = params(4);
p_boost = cell2mat(cell(p_boost{1}));

tau = params(1);
a = cell(tau{1});
tau = double(a{1});


K = params(6);
K = (double(K{1}));
R = params(7);
R = double(R{1});


fid = py.open('boosted_nodes.pckl','rb');
boosted_nodes = py.pickle.load(fid);




cd '..\..\Data_analysis\Metrics'

fid = py.open('all_hubs_15.pckl','rb');
Hubs = py.pickle.load(fid);

fid = py.open('all_hubs_nodeNo_15.pckl','rb');
Hubs_nodeNo = py.pickle.load(fid);

fid = py.open('all_cd_pairs_15.pckl','rb');  
cd_pairs = py.pickle.load(fid);





outdeg_neighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
outdeg_neighpacemaker_NoOutput_network = zeros(length(p_dist),length(p_boost));
outdeg_NoNneighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
outdeg_NoNneighpacemaker_NoOutput_network = zeros(length(p_dist),length(p_boost));
indeg_neighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
indeg_NoNneighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
inHubs_neighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
inHubs_NoNneighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
inHubsCD_neighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));
inHubsCD_NoNneighpacemaker_mean_network = zeros(length(p_dist),length(p_boost));

indeg_neighpacemaker_NoInput_mean_network = zeros(length(p_dist),length(p_boost));
indeg_NoNneighpacemaker_NoInput_mean_network = zeros(length(p_dist),length(p_boost));


ii = 0;
for i = p_dist
    ii = ii+1;
    jj = 0;
    for j = p_boost
        jj = jj+1;

        if (i+j) <= 1

            kk = 0;
            for k = 0:K-1
                kk = kk+1;

                rr = 0;
                for r = 0:R-1
                    rr = rr+1;

                    node = boosted_nodes{k};
                    node = double(node{1})+1; %+1 because in python vectors start from position 0, but in matlab from 1
                        
                    COORDS = double( coords{1} {k} );
                    xdata = COORDS(:,1)';
                    ydata = COORDS(:,2)';
                    ADJ_mat = double( adj{1} {tau}{0.5}{i}{j}{k}{r} ); %A_matrices[tau][p_in][pDist][pBoost][pBoostEach][random_every][k][r]


                    G = digraph(ADJ_mat');    
                    if neighbours == 0
                        neighPacemaker = successors(G,node);%neighbors are nodes targeted by initiator,
                    elseif neighbours == 1
                        neighPacemaker = sort([successors(G,node)' predecessors(G,node)'])'; %neighbors are nodes targeted by initiator, and nodes targeting initiator
                    end
                    NoNneighPacemaker = node_list( not(ismember(node_list,neighPacemaker')) );
                    
                
%% OUT-DEGREE                
                    %Neighbors
                    if isempty(neighPacemaker)==1
                        outdeg_neighpacemaker = NaN;
                        neighpacemaker_NoOutput = NaN;
                    else
                        for iii = 1:length(neighPacemaker)
                            outdeg_neighpacemaker(iii) = length(successors(G,neighPacemaker(iii)));
                        end
                        neighpacemaker_NoOutput = sum(outdeg_neighpacemaker==0) / length(outdeg_neighpacemaker);
                    end


                    outdeg_neighpacemaker_mean_r(rr) = mean(outdeg_neighpacemaker);
                    outdeg_neighpacemaker_NoOutput_r(rr) = neighpacemaker_NoOutput;


                    %non-Neighbors
                    if isempty(NoNneighPacemaker)==1
                        outdeg_NoNneighpacemaker = NaN;
                        NoNneighpacemaker_NoOutput = NaN;
                    else
                        for iii = 1:length(NoNneighPacemaker)
                            outdeg_NoNneighpacemaker(iii) = length(successors(G,NoNneighPacemaker(iii)));
                        end
                        NoNneighpacemaker_NoOutput = sum(outdeg_NoNneighpacemaker==0) / length(outdeg_NoNneighpacemaker);
                    end
                    outdeg_NoNneighpacemaker_mean_r(rr) = mean(outdeg_NoNneighpacemaker);
                    outdeg_NoNneighpacemaker_NoOutput_r(rr) = NoNneighpacemaker_NoOutput;
                  
%% IN-DEGREE
                    %Neighbors
                    if isempty(neighPacemaker)==1
                        indeg_neighpacemaker = NaN;
                        neighpacemaker_NoInput = NaN;
                    else
                        for iii = 1:length(neighPacemaker)
                            indeg_neighpacemaker(iii) = length(predecessors(G,neighPacemaker(iii)));
                        end
                        neighpacemaker_NoInput = sum(indeg_neighpacemaker==0) / length(indeg_neighpacemaker);
                    end
                    indeg_neighpacemaker_mean_r(rr) = mean(indeg_neighpacemaker);
                    indeg_neighpacemaker_NoInput_r(rr) = neighpacemaker_NoInput;


                    %non-Neighbors
                    if isempty(NoNneighPacemaker)==1
                        indeg_NoNneighpacemaker = NaN;
                        NoNneighpacemaker_NoInput = NaN;
                    else
                        for iii = 1:length(NoNneighPacemaker)
                            indeg_NoNneighpacemaker(iii) = length(predecessors(G,NoNneighPacemaker(iii)));
                        end
                        NoNneighpacemaker_NoInput = sum(indeg_NoNneighpacemaker==0) / length(indeg_NoNneighpacemaker);
                    end
                    indeg_NoNneighpacemaker_mean_r(rr) = mean(indeg_NoNneighpacemaker);  
                    indeg_NoNneighpacemaker_NoInput_r(rr) = NoNneighpacemaker_NoInput;


                    %% Neighbors and non-neighbors that form in-hubs
                    Hubs_nodeNo_in = int64(Hubs_nodeNo{'in'}{tau}{0.5}{i}{j}{k}{r}) + 1;
                    %Neighbors
                    if isempty(neighPacemaker)==1
                        inHubs_neighpacemaker = NaN;
                    else
                        inHubs_neighpacemaker = sum( ismember(neighPacemaker,Hubs_nodeNo_in) )  / length(neighPacemaker);                 
                    end
                    inHubs_neighpacemaker_mean_r(rr) = inHubs_neighpacemaker;
                    %Non-Neighbors
                    if isempty(NoNneighPacemaker)==1
                        inHubs_NoNneighpacemaker = NaN;
                    else
                        inHubs_NoNneighpacemaker = sum( ismember(NoNneighPacemaker,Hubs_nodeNo_in) )  / length(NoNneighPacemaker);                
                    end
                    inHubs_NoNneighpacemaker_mean_r(rr) = inHubs_NoNneighpacemaker;

                    %% Neighbors and non-neighbors that are convergent hubs in cd units
                    conv_div = int64(cd_pairs{tau}{0.5}{i}{j}{k}{r}) + 1;
                    convergent = conv_div(1,:);
                    convergent = unique(convergent(:).');
                    %Neighbors
                    if isempty(neighPacemaker)==1
                        inHubsCD_neighpacemaker = NaN;
                    else
                        inHubsCD_neighpacemaker = sum( ismember(neighPacemaker,convergent) )  / length(neighPacemaker);                 
                    end
                    inHubsCD_neighpacemaker_mean_r(rr) = inHubsCD_neighpacemaker;
                    %Non-Neighbors
                    if isempty(NoNneighPacemaker)==1
                        inHubsCD_NoNneighpacemaker = NaN;
                    else
                        inHubsCD_NoNneighpacemaker = sum( ismember(NoNneighPacemaker,convergent) )  / length(NoNneighPacemaker);                 
                    end
                    inHubsCD_NoNneighpacemaker_mean_r(rr) = inHubsCD_NoNneighpacemaker;

    
                end

                outdeg_neighpacemaker_mean_k(kk) = nanmean(outdeg_neighpacemaker_mean_r);
                outdeg_neighpacemaker_NoOutput_k(kk) = nanmean(outdeg_neighpacemaker_NoOutput_r);

                outdeg_NoNneighpacemaker_mean_k(kk) = nanmean(outdeg_NoNneighpacemaker_mean_r);
                outdeg_NoNneighpacemaker_NoOutput_k(kk) = nanmean(outdeg_NoNneighpacemaker_NoOutput_r);

                indeg_neighpacemaker_mean_k(kk) = nanmean(indeg_neighpacemaker_mean_r);

                indeg_NoNneighpacemaker_mean_k(kk) = nanmean(indeg_NoNneighpacemaker_mean_r);

                inHubs_neighpacemaker_mean_k(kk) = nanmean(inHubs_neighpacemaker_mean_r);
                inHubs_NoNneighpacemaker_mean_k(kk) = nanmean(inHubs_NoNneighpacemaker_mean_r);

                inHubsCD_neighpacemaker_mean_k(kk) = nanmean(inHubsCD_neighpacemaker_mean_r);
                inHubsCD_NoNneighpacemaker_mean_k(kk) = nanmean(inHubsCD_NoNneighpacemaker_mean_r);


                indeg_neighpacemaker_NoInput_k(kk) = nanmean(indeg_neighpacemaker_NoInput_r);
                indeg_NoNneighpacemaker_NoInput_k(kk) = nanmean(indeg_NoNneighpacemaker_NoInput_r);

            end

            outdeg_neighpacemaker_mean_network(ii,jj) = nanmean(outdeg_neighpacemaker_mean_k);
            outdeg_neighpacemaker_NoOutput_network(ii,jj) = nanmean(outdeg_neighpacemaker_NoOutput_k);

            outdeg_NoNneighpacemaker_mean_network(ii,jj) = nanmean(outdeg_NoNneighpacemaker_mean_k);
            outdeg_NoNneighpacemaker_NoOutput_network(ii,jj) = nanmean(outdeg_NoNneighpacemaker_NoOutput_k);

            indeg_neighpacemaker_mean_network(ii,jj) = nanmean(indeg_neighpacemaker_mean_k);

            indeg_NoNneighpacemaker_mean_network(ii,jj) = nanmean(indeg_NoNneighpacemaker_mean_k);

            inHubs_neighpacemaker_mean_network(ii,jj) = nanmean(inHubs_neighpacemaker_mean_k);
            inHubs_NoNneighpacemaker_mean_network(ii,jj) = nanmean(inHubs_NoNneighpacemaker_mean_k);

            inHubsCD_neighpacemaker_mean_network(ii,jj) = nanmean(inHubsCD_neighpacemaker_mean_k);
            inHubsCD_NoNneighpacemaker_mean_network(ii,jj) = nanmean(inHubsCD_NoNneighpacemaker_mean_k);

            indeg_neighpacemaker_NoInput_mean_network(ii,jj) = nanmean(indeg_neighpacemaker_NoInput_k);
            indeg_NoNneighpacemaker_NoInput_mean_network(ii,jj) = nanmean(indeg_NoNneighpacemaker_NoInput_k);

        end

    end
end



Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;


outdeg_neighpacemaker_mean_network = flipud(outdeg_neighpacemaker_mean_network).*Zero_diag;

ax1 = figure(1)
b = imagesc(outdeg_neighpacemaker_mean_network);
set(b,'AlphaData',~isnan(outdeg_neighpacemaker_mean_network)) %~isnan(idx)
if neighbours == 0
    caxis manual
    caxis([0 14]); 
elseif neighbours == 1
    caxis manual
    caxis([0 14]); %caxis([3 13.5])
end

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_neighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_neighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Out-degree neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




outdeg_NoNneighpacemaker_mean_network = flipud(outdeg_NoNneighpacemaker_mean_network).*Zero_diag;

ax1 = figure(2)
b = imagesc(outdeg_NoNneighpacemaker_mean_network);
set(b,'AlphaData',~isnan(outdeg_NoNneighpacemaker_mean_network)) %~isnan(idx)
if neighbours == 0
    caxis manual
    caxis([0 14]); 
elseif neighbours == 1
    caxis manual
    caxis([0 14]); %caxis([3 13.5])
end

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_NoNneighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_NoNneighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Out-degree non-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



ratio_outdeg = outdeg_neighpacemaker_mean_network./outdeg_NoNneighpacemaker_mean_network;

ax1 = figure(3)
b = imagesc(ratio_outdeg );
set(b,'AlphaData',~isnan(ratio_outdeg)) %~isnan(idx)
% caxis manual
% caxis([0 14]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(ratio_outdeg,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(ratio_outdeg,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Raio: out-deg neigh / out-deg non-neigh')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





indeg_neighpacemaker_mean_network = flipud(indeg_neighpacemaker_mean_network).*Zero_diag;

ax1 = figure(4)
b = imagesc(indeg_neighpacemaker_mean_network);
set(b,'AlphaData',~isnan(indeg_neighpacemaker_mean_network)) %~isnan(idx)
if neighbours == 0
    caxis manual
    caxis([0 30]); 
elseif neighbours == 1
    caxis manual
    caxis([0 16]); %caxis([3 13.5])
end

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(indeg_neighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(indeg_neighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('In-degree neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




indeg_NoNneighpacemaker_mean_network = flipud(indeg_NoNneighpacemaker_mean_network).*Zero_diag;

ax1 = figure(5)
b = imagesc(indeg_NoNneighpacemaker_mean_network);
set(b,'AlphaData',~isnan(indeg_NoNneighpacemaker_mean_network)) %~isnan(idx)
if neighbours == 0
    caxis manual
    caxis([0 30]); 
elseif neighbours == 1
    caxis manual
    caxis([0 16]); %caxis([3 13.5])
end

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(indeg_NoNneighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(indeg_NoNneighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('In-degree non-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




ratio_indeg = indeg_neighpacemaker_mean_network./indeg_NoNneighpacemaker_mean_network;

ax1 = figure(6)
b = imagesc(ratio_indeg);
set(b,'AlphaData',~isnan(ratio_indeg)) %~isnan(idx)
% caxis manual
% caxis([0 30]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(ratio_indeg,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(ratio_indeg,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Raio: in-deg neigh / in-deg non-neigh')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square






outdeg_neighpacemaker_NoOutput_network = flipud(outdeg_neighpacemaker_NoOutput_network).*Zero_diag;

ax1 = figure(7)
b = imagesc(outdeg_neighpacemaker_NoOutput_network);
set(b,'AlphaData',~isnan(outdeg_neighpacemaker_NoOutput_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_neighpacemaker_NoOutput_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_neighpacemaker_NoOutput_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion zero Out-degree among neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




outdeg_NoNneighpacemaker_NoOutput_network = flipud(outdeg_NoNneighpacemaker_NoOutput_network).*Zero_diag;

ax1 = figure(8)
b = imagesc(outdeg_NoNneighpacemaker_NoOutput_network);
set(b,'AlphaData',~isnan(outdeg_NoNneighpacemaker_NoOutput_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_NoNneighpacemaker_NoOutput_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_NoNneighpacemaker_NoOutput_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion zero Out-degree among non-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square







indeg_neighpacemaker_NoInput_mean_network = flipud(indeg_neighpacemaker_NoInput_mean_network).*Zero_diag;

ax1 = figure(9)
b = imagesc(indeg_neighpacemaker_NoInput_mean_network);
set(b,'AlphaData',~isnan(indeg_neighpacemaker_NoInput_mean_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_neighpacemaker_NoOutput_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_neighpacemaker_NoOutput_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion zero In-degree among neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




indeg_NoNneighpacemaker_NoInput_mean_network = flipud(indeg_NoNneighpacemaker_NoInput_mean_network).*Zero_diag;

ax1 = figure(10)
b = imagesc(indeg_NoNneighpacemaker_NoInput_mean_network);
set(b,'AlphaData',~isnan(indeg_NoNneighpacemaker_NoInput_mean_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(outdeg_NoNneighpacemaker_NoOutput_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(outdeg_NoNneighpacemaker_NoOutput_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion zero In-degree among NoN-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square








inHubs_neighpacemaker_mean_network = flipud(inHubs_neighpacemaker_mean_network).*Zero_diag;

ax1 = figure(11)
b = imagesc(inHubs_neighpacemaker_mean_network);
set(b,'AlphaData',~isnan(inHubs_neighpacemaker_mean_network)) %~isnan(idx)
caxis manual
caxis([0 0.5]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(inHubs_neighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(inHubs_neighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('proportion in-hubs among neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



inHubs_NoNneighpacemaker_mean_network = flipud(inHubs_NoNneighpacemaker_mean_network).*Zero_diag;

ax1 = figure(12)
b = imagesc(inHubs_NoNneighpacemaker_mean_network);
set(b,'AlphaData',~isnan(inHubs_NoNneighpacemaker_mean_network)) %~isnan(idx)
caxis manual
caxis([0 0.5]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(inHubs_NoNneighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(inHubs_NoNneighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('proportion in-hubs among non-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square






inHubsCD_neighpacemaker_mean_network = flipud(inHubsCD_neighpacemaker_mean_network).*Zero_diag;

ax1 = figure(13)
b = imagesc(inHubsCD_neighpacemaker_mean_network);
set(b,'AlphaData',~isnan(inHubsCD_neighpacemaker_mean_network)) %~isnan(idx)
caxis manual
caxis([0 0.5]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(inHubsCD_neighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(inHubsCD_neighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('proportion c-d Unit in-hubs among neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



inHubsCD_NoNneighpacemaker_mean_network = flipud(inHubsCD_NoNneighpacemaker_mean_network).*Zero_diag;

ax1 = figure(14)
b = imagesc(inHubsCD_NoNneighpacemaker_mean_network);
set(b,'AlphaData',~isnan(inHubsCD_NoNneighpacemaker_mean_network)) %~isnan(idx)
caxis manual
caxis([0 0.5]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(inHubsCD_NoNneighpacemaker_mean_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(inHubsCD_NoNneighpacemaker_mean_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('proportion c-d Unit in-hubs among non-neighbors initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





