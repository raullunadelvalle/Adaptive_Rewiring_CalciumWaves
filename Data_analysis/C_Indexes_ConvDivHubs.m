clear all
close all
clc


% Use this script to generate figures similar to those in Figure 2 and Figure 4 of:
%Luna, R., Li, J., Bauer, R., & van Leeuwen, C. (2023). Retinal waves in adaptive 
%rewiring networks orchestrate convergence and divergence in the visual system.
%Network Neuroscience




cd '..\scripts\Output'

fid = py.open('boost_nor_MatPlusCoords.pckl','rb'); 
coords = py.pickle.load(fid);
num_nodes = size(double(coords{1}{1}),1);


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

fid = py.open('mean_modularity.pckl','rb');  
data_mod = py.pickle.load(fid);

fid = py.open('all_cp.pckl','rb');
conNpairs = py.pickle.load(fid);

fid = py.open('all_hubs_15.pckl','rb');
Hubs = py.pickle.load(fid);

fid = py.open('all_hubs_nodeNo_15.pckl','rb');
Hubs_nodeNo = py.pickle.load(fid);

fid = py.open('all_pl.pckl','rb');
pathLength = py.pickle.load(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


mod = zeros(length(p_dist),length(p_boost));

prop_conNpairs_network = zeros(length(p_dist),length(p_boost));
prop_convHubs_network = zeros(length(p_dist),length(p_boost));
prop_divHubs_network = zeros(length(p_dist),length(p_boost));
avEff_network = zeros(length(p_dist),length(p_boost));

initiator_isDivHub_network = zeros(length(p_dist),length(p_boost));
randNode_isDivHub_network = zeros(length(p_dist),length(p_boost));


ii = 0;
for i = p_dist
    ii = ii+1;
    jj = 0;
    for j = p_boost
        jj = jj+1;

        if (i+j) <= 1

            mod(ii,jj) = data_mod{tau}{0.5}{i}{j};

            prop_conNpairs_network_k = [];
            prop_convHubs_network_k = [];
            prop_divHubs_network_k = [];
            avEff_network_k = [];

            initiator_isDivHub_network_k = [];
            randNode_isDivHub_network_k = [];

            kk = 0;
            for k = 0:K-1
                kk = kk+1;

                prop_conNpairs_network_r = [];
                prop_convHubs_network_r = [];
                prop_divHubs_network_r = [];
                avEff_network_r = [];

                initiator_isDivHub_network_r = [];
                randNode_isDivHub_network_r = [];

                rr = 0;
                for r = 0:R-1
                    rr = rr+1;
                    prop_conNpairs = conNpairs{tau}{0.5}{i}{j}{k}{r};
                    prop_conNpairs = double(prop_conNpairs);
                    prop_conNpairs_network_r(rr) = sum(sum(prop_conNpairs))/(size(prop_conNpairs,1)*size(prop_conNpairs,2));

                    prop_convHubs = Hubs{'in'}{tau}{0.5}{i}{j}{k}{r};
                    prop_convHubs = double(prop_convHubs);
                    prop_convHubs_network_r(rr) = prop_convHubs/num_nodes;   

                    prop_divHubs = Hubs{'out'}{tau}{0.5}{i}{j}{k}{r};

                    prop_divHubs = double(prop_divHubs);
                    prop_divHubs_network_r(rr) = prop_divHubs/num_nodes;   

                    avEff = pathLength{tau}{0.5}{i}{j}{k}{r};
                    avEff = double(avEff);
                    avEff_network_r(rr) = (1/(num_nodes*(num_nodes-1))) * sum(sum(1./avEff));

                    % Does the boosted node belong in the list of divergent hubs?
                    Hubs_nodeNo_out = int64(Hubs_nodeNo{'out'}{tau}{0.5}{i}{j}{k}{r});
                    node = boosted_nodes{k};
                    node = double(node{1});


                    if isempty(Hubs_nodeNo_out) == 1
                    initiator_isDivHub_network_r(rr) = NaN;
                    randNode_isDivHub_network_r(rr) = NaN;

                    else

                    if ismember(node,Hubs_nodeNo_out) == 1
                        initiator_isDivHub = 1;
                    else
                        initiator_isDivHub = 0;
                    end

                    initiator_isDivHub_network_r(rr) = initiator_isDivHub;
                    %%%%%%%%


                    %random nodes. How often do they become divergent hubs?
                    randNode_isDivHub_g = [];
                    gg = 0;
                    for g = 1:1000
                        gg = gg+1;
                        rand_node = node;
                        while rand_node == node
                            rand_node = randi([0 99]);
                        end

                        if ismember(rand_node,Hubs_nodeNo_out) == 1
                            randNode_isDivHub = 1;
                        else
                            randNode_isDivHub = 0;
                        end
                        randNode_isDivHub_g(gg) = randNode_isDivHub;
                    end

                    randNode_isDivHub = nansum(randNode_isDivHub_g)/g;
                    randNode_isDivHub_network_r(rr) = randNode_isDivHub;

                    end

                end

                prop_conNpairs_network_k(kk) = nanmean(prop_conNpairs_network_r);
                prop_convHubs_network_k(kk) = nanmean(prop_convHubs_network_r);
                prop_divHubs_network_k(kk) = nanmean(prop_divHubs_network_r);
                avEff_network_k(kk) = nanmean(avEff_network_r);

                initiator_isDivHub_network_k(kk) = nanmean(initiator_isDivHub_network_r);
                randNode_isDivHub_network_k(kk) = nanmean(randNode_isDivHub_network_r);
            end

            prop_conNpairs_network(ii,jj) = nanmean(prop_conNpairs_network_k);
            prop_convHubs_network(ii,jj) = nanmean(prop_convHubs_network_k);
            prop_divHubs_network(ii,jj) = nanmean(prop_divHubs_network_k);
            avEff_network(ii,jj) = nanmean(avEff_network_k);

            initiator_isDivHub_network(ii,jj) = nansum(initiator_isDivHub_network_k) / (kk-sum(isnan(initiator_isDivHub_network_k))); %nanmean(initiator_isDivHub_network_k);
            randNode_isDivHub_network(ii,jj) = nansum(randNode_isDivHub_network_k) / (kk-sum(isnan(randNode_isDivHub_network_k))); %nanmean(randNode_isDivHub_network_k);
        end

    end
end


Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

mod = flipud(mod).*Zero_diag;

ax2 = figure(1)
b = imagesc(mod);
set(b,'AlphaData',~isnan(mod)) %~isnan(idx)
%caxis manual
%caxis([0 0.65]); 

colormap(ax2,cool)
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(mod,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(mod,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Modularity')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

avEff_network = flipud(avEff_network).*Zero_diag;

ax1 = figure(2)
b = imagesc(avEff_network);
set(b,'AlphaData',~isnan(avEff_network)) %~isnan(idx)
%caxis manual
%caxis([-10 45]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(avEff_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(avEff_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Average efficiency')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

prop_conNpairs_network = flipud(prop_conNpairs_network).*Zero_diag;

ax1 = figure(3)
b = imagesc(prop_conNpairs_network);
set(b,'AlphaData',~isnan(prop_conNpairs_network)) %~isnan(idx)
%caxis manual
%caxis([-10 45]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(prop_conNpairs_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(prop_conNpairs_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion of connected pairs')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

prop_convHubs_network = flipud(prop_convHubs_network).*Zero_diag;

ax1 = figure(4)
b = imagesc(prop_convHubs_network);
set(b,'AlphaData',~isnan(prop_convHubs_network)) %~isnan(idx)
caxis manual
caxis([0.04 0.1224]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(prop_convHubs_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(prop_convHubs_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion of convergent hubs')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square






Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

prop_divHubs_network = flipud(prop_divHubs_network).*Zero_diag;

ax1 = figure(5)
b = imagesc(prop_divHubs_network);
set(b,'AlphaData',~isnan(prop_divHubs_network)) %~isnan(idx)
caxis manual
caxis([0.04 0.1224]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(prop_divHubs_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(prop_divHubs_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Proportion of divergent hubs')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;


initiator_isDivHub_network = flipud(initiator_isDivHub_network).*Zero_diag;

ax1 = figure(6)
b = imagesc(initiator_isDivHub_network);
set(b,'AlphaData',~isnan(initiator_isDivHub_network)) %~isnan(idx)
caxis manual
caxis([0 0.7]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(initiator_isDivHub_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(initiator_isDivHub_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Succes rate initiator becomes div Hub')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;


randNode_isDivHub_network = flipud(randNode_isDivHub_network).*Zero_diag;

ax1 = figure(7)
b = imagesc(randNode_isDivHub_network);
set(b,'AlphaData',~isnan(randNode_isDivHub_network)) %~isnan(idx)
caxis manual
caxis([0 0.7]); %caxis([0.04 0.1224]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(randNode_isDivHub_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(randNode_isDivHub_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Succes rate non-initiator becomes div Hub')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square






