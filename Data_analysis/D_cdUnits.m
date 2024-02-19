clear all
close all
clc


% Use this script to generate figures similar to those in Figure 5 and Figure 6 of:
%Luna, R., Li, J., Bauer, R., & van Leeuwen, C. (2023). Retinal waves in adaptive 
%rewiring networks orchestrate convergence and divergence in the visual system.
%Network Neuroscience



cd '..\scripts\Output'

fid=py.open('simulation_params.pckl','rb'); 
params = py.pickle.load(fid);

p_dist = params(3);
p_dist = cell2mat(cell(p_dist{1}));

p_boost = params(4);
p_boost = cell2mat(cell(p_boost{1}));

K = params(6);
K = (double(K{1})); % It is a problem that I have to subtract 2
R = params(7);
R = double(R{1});

tau = params(1);
a = cell(tau{1});
tau = double(a{1});

fid = py.open('boosted_nodes.pckl','rb');
boosted_nodes = py.pickle.load(fid);





cd '..\..\Data_analysis\Metrics'

fid = py.open('all_cd_pairs_15.pckl','rb');  
cd_pairs = py.pickle.load(fid);

fid = py.open('all_interm_size_15.pckl','rb'); 
interm_size = py.pickle.load(fid);

fid = py.open('all_interm_density_15.pckl','rb'); 
interm_density = py.pickle.load(fid);

fid = py.open('all_periph_density_15.pckl','rb'); 
periph_density = py.pickle.load(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




num = zeros(length(p_dist),length(p_boost));
num_network = zeros(length(p_dist),length(p_boost));

num_cd_units = zeros(length(p_dist),length(p_boost));



it_size = zeros(length(p_dist),length(p_boost));
it_density = zeros(length(p_dist),length(p_boost));
pr_density = zeros(length(p_dist),length(p_boost));

it_size_network = zeros(length(p_dist),length(p_boost));
it_density_network = zeros(length(p_dist),length(p_boost));
pr_density_network = zeros(length(p_dist),length(p_boost));

ii = 0;
for i = p_dist
    ii = ii+1;
    jj = 0;
    for j = p_boost
        jj = jj+1;

        if (i+j) <= 1


            num_conv_div_network_k = [];
            num_conv_div_k = [];

            num_cd_units_k = [];

            interm_size_k = [];
            interm_density_k = [];
            periph_density_k = [];

            interm_size_network_k = [];
            interm_density_network_k = [];
            periph_density_network_k = [];
            

            for k = 0:K-1

                conv_div_network_r = [];
                num_conv_div_r = [];

                num_cd_units_r = [];

                interm_size_r = [];
                interm_density_r = [];
                periph_density_r = [];

                interm_size_network_r = [];
                interm_density_network_r = [];
                periph_density_network_r = [];

                for r = 0:R-1
                    node = boosted_nodes{k};
                    node = double(node{1});
        
                    conv_div = int64(cd_pairs{tau}{0.5}{i}{j}{k}{r});

                    number_cd_units = size(conv_div,2);

                    
                    % Network has at least 1 convergent-divergent unit,
                    % regardless of it including the pacemaker
                    if isempty(conv_div) == 1
                        conv_div_network = 0;

                        int_size_network = NaN;
                        int_density_network = NaN;
                        per_density_network = NaN;
                    else
                        conv_div_network = 1;

                        int_size = interm_size{tau}{0.5}{i}{j}{k}{r};
                        int_size = cell(int_size);
                        int_size = cellfun(@double,int_size);
                        int_size_network = int_size(find(conv_div(2,:) ~= node));

                        int_density = interm_density{tau}{0.5}{i}{j}{k}{r};
                        int_density = cell(int_density);
                        int_density = cellfun(@double,int_density);
                        int_density_network = int_density(find(conv_div(2,:) ~= node));

                        per_density = periph_density{tau}{0.5}{i}{j}{k}{r};
                        per_density = cell(per_density);
                        per_density = cellfun(@double,per_density);
                        per_density_network = per_density(find(conv_div(2,:) ~= node));
                    end

                    conv_div_network_r = [conv_div_network_r conv_div_network];

                    num_cd_units_r = [num_cd_units_r number_cd_units];


                    interm_size_network_r = [interm_size_network_r int_size_network];
                    interm_density_network_r = [interm_density_network_r int_density_network];
                    periph_density_network_r = [periph_density_network_r per_density_network];

                    % Convergent divergent unit that includes pacemaker,
                    % there is at least 1
                    if ismember(node,conv_div(2,:)) == 1
                        num_conv_div = 1;

                        int_size = interm_size{tau}{0.5}{i}{j}{k}{r};
                        int_size = cell(int_size);
                        int_size = cellfun(@double,int_size);
                        int_size_pacemaker = int_size(find(conv_div(2,:) == node));

                        int_density = interm_density{tau}{0.5}{i}{j}{k}{r};
                        int_density = cell(int_density);
                        int_density = cellfun(@double,int_density);
                        int_density_pacemaker = int_density(find(conv_div(2,:) == node));

                        per_density = periph_density{tau}{0.5}{i}{j}{k}{r};
                        per_density = cell(per_density);
                        per_density = cellfun(@double,per_density);
                        per_density_pacemaker = per_density(find(conv_div(2,:) == node));

                    else
                        num_conv_div = 0;

                        int_size_pacemaker = NaN;
                        int_density_pacemaker = NaN;
                        per_density_pacemaker = NaN;
                    end
                    num_conv_div_r = [num_conv_div_r num_conv_div];
                    interm_size_r = [interm_size_r int_size_pacemaker];
                    interm_density_r = [interm_density_r int_density_pacemaker];
                    periph_density_r = [periph_density_r per_density_pacemaker];
                    
                end

                % Network has at least 1 convergent-divergent unit,
                % regardless of it including the pacemaker
                num_conv_div_network = nansum(conv_div_network_r);
                num_conv_div_network_k = [num_conv_div_network_k num_conv_div_network];

                num_cd_units_network = nanmean(num_cd_units_r);
                num_cd_units_k = [num_cd_units_k num_cd_units_network];


                int_size_network = nanmean(interm_size_network_r);
                interm_size_network_k = [interm_size_network_k int_size_network];

                int_density_network = nanmean(interm_density_network_r);
                interm_density_network_k = [interm_density_network_k int_density_network];

                per_density_network = nanmean(periph_density_network_r);
                periph_density_network_k = [periph_density_network_k per_density_network];




                % Convergent divergent unit that includes pacemaker
                num_conv_div = nansum(num_conv_div_r);
                num_conv_div_k = [num_conv_div_k num_conv_div];

                int_size = nanmean(interm_size_r);
                interm_size_k = [interm_size_k int_size];

                int_density = nanmean(interm_density_r);
                interm_density_k = [interm_density_k int_density];

                per_density = nanmean(periph_density_r);
                periph_density_k = [periph_density_k per_density];
            end

            % Proportion of times networks develop cd units     Wrong:Number Convergent divergent units in the network
            num_network(ii,jj) = nanmean(num_conv_div_network_k);

            % Average Number of cd units in networks
            num_cd_units(ii,jj) = nanmean(num_cd_units_k);


            % Intermediate size Convergent divergent unit in the network
            it_size_network(ii,jj) = nanmean(interm_size_network_k);

            % Intermediate density Convergent divergent unit in the network
            it_density_network(ii,jj) = nanmean(interm_density_network_k);

            % Peripheral density Convergent divergent unit in the network
            pr_density_network(ii,jj) = nanmean(periph_density_network_k);




            % Number Convergent divergent unit that includes pacemaker
            num(ii,jj) = nanmean(num_conv_div_k);

            % Intermediate size Convergent divergent unit that includes pacemaker
            it_size(ii,jj) = nanmean(interm_size_k);

            % Intermediate density Convergent divergent unit that includes pacemaker
            it_density(ii,jj) = nanmean(interm_density_k);

            % Peripheral density Convergent divergent unit that includes pacemaker
            pr_density(ii,jj) = nanmean(periph_density_k);

        end

    end
end



Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;



num_cd_units = flipud(num_cd_units).*Zero_diag;



ax1 = figure(1)
b = imagesc(num_cd_units);
set(b,'AlphaData',~isnan(num_cd_units)) %~isnan(targeted)
% caxis manual
% caxis([0 50]); 

colormap(ax1,cool)
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(num_cd_units,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(num_cd_units,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Total number cd units')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



num_network = flipud(num_network/R).*Zero_diag;

figure(2)
b = imagesc(num_network);
set(b,'AlphaData',~isnan(num_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(num_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(num_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('c-d unit network success')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



num = flipud(num/R).*Zero_diag;

figure(3)
b = imagesc(num);
set(b,'AlphaData',~isnan(num)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap cool
%cob = colorbar;
c.LineWidth = 3;

xticks([1:size(num,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(num,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('c-d unit pacemaker success')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




figure(4)

it_size = flipud(it_size).*Zero_diag;

b = imagesc(it_size);
set(b,'AlphaData',~isnan(it_size)) %~isnan(idx)
caxis manual
caxis([0 100]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(it_size,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(it_size,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('pacemaker interm sub size')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




it_density = flipud(it_density).*Zero_diag;

figure(5)
b = imagesc(it_density);
set(b,'AlphaData',~isnan(it_density)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(it_density,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(it_size,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('pacemaker interm sub density')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





figure(6)

it_size_network = flipud(it_size_network).*Zero_diag;

b = imagesc(it_size_network);
set(b,'AlphaData',~isnan(it_size_network)) %~isnan(idx)
caxis manual
caxis([0 100]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(it_size_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(it_size_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('network interm sub size')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square


it_density_network = flipud(it_density_network).*Zero_diag;

figure(40)
b = imagesc(it_density_network);
set(b,'AlphaData',~isnan(it_density_network)) %~isnan(idx)
caxis manual
caxis([0 1]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(it_density_network,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(it_size_network,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('network interm sub density')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square



