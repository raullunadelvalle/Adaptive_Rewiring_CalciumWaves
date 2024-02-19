clear all
close all
clc

% Use this script to generate figures similar to those in Figure 3A of:
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



% Average in-degree targeted
In_degree_Targeted_runs = {};
% Average in-degree non-targeted
In_degree_nonTargeted_runs = {};
% Average out-degree targeted
Out_degree_Targeted_runs = {};
% Average out-degree non-targeted
Out_degree_nonTargeted_runs = {};



fid = py.open('boosted_nodes.pckl','rb');
boosted_nodes = py.pickle.load(fid);


cd '..\..\Data_analysis\Metrics'

fid = py.open('InOut_degree.pckl','rb');
Degree = py.pickle.load(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


In_degree_Targeted = zeros(length(p_dist),length(p_boost));
In_degree_nonTargeted = zeros(length(p_dist),length(p_boost));
Out_degree_Targeted = zeros(length(p_dist),length(p_boost));
Out_degree_nonTargeted = zeros(length(p_dist),length(p_boost));

ii = 0;
for i = p_dist
    ii = ii+1;
    jj = 0;
    for j = p_boost
        jj = jj+1;

        if (i+j) <= 1

            In_degree_Targeted_k = [];
            In_degree_nonTargeted_k = [];
            Out_degree_Targeted_k = [];
            Out_degree_nonTargeted_k = [];

            kk = 0;
            for k = 0:K-1
                kk = kk+1;

                In_degree_Targeted_r = [];
                In_degree_nonTargeted_r = [];
                Out_degree_Targeted_r = [];
                Out_degree_nonTargeted_r = [];

                rr = 0;
                for r = 0:R-1
                    rr = rr+1;

                    node = boosted_nodes{k};
                    node = double(node{1});
                    
                    Degree_in = double(Degree{'in'}{tau}{0.5}{i}{j}{k}{r});
                    deg_Targeted_in = Degree_in(2,find(Degree_in(1,:)==node));
                    In_degree_Targeted_r(rr) = deg_Targeted_in;
                    deg_nonTargeted_in = Degree_in(2,find(Degree_in(1,:)~=node));
                    In_degree_nonTargeted_r(rr) = mean(deg_nonTargeted_in);

                    Degree_out = double(Degree{'out'}{tau}{0.5}{i}{j}{k}{r});
                    deg_Targeted_out = Degree_out(2,find(Degree_out(1,:)==node));
                    Out_degree_Targeted_r(rr) = deg_Targeted_out;
                    deg_nonTargeted_out = Degree_out(2,find(Degree_out(1,:)~=node));
                    Out_degree_nonTargeted_r(rr) = mean(deg_nonTargeted_out);
                end

                In_degree_Targeted_k(kk) = nanmean(In_degree_Targeted_r);
                In_degree_nonTargeted_k(kk) = nanmean(In_degree_nonTargeted_r);
                Out_degree_Targeted_k(kk) = nanmean(Out_degree_Targeted_r);
                Out_degree_nonTargeted_k(kk) = nanmean(Out_degree_nonTargeted_r);
            end

            In_degree_Targeted(ii,jj) = nanmean(In_degree_Targeted_k);
            In_degree_nonTargeted(ii,jj) = nanmean(In_degree_nonTargeted_k);
            Out_degree_Targeted(ii,jj) = nanmean(Out_degree_Targeted_k);
            Out_degree_nonTargeted(ii,jj) = nanmean(Out_degree_nonTargeted_k);
        end

    end
end




%% In-degree

Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

In_degree_Targeted = flipud(In_degree_Targeted).*Zero_diag;

ax1 = figure(1)
b = imagesc(In_degree_Targeted);
set(b,'AlphaData',~isnan(In_degree_Targeted)) %~isnan(idx)
% caxis manual
% caxis([1 10]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(In_degree_Targeted,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(In_degree_Targeted,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('In-degree Initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

In_degree_nonTargeted = flipud(In_degree_nonTargeted).*Zero_diag;

ax1 = figure(2)
b = imagesc(In_degree_nonTargeted);
set(b,'AlphaData',~isnan(In_degree_nonTargeted)) %~isnan(idx)
% caxis manual
% caxis([1 10]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(In_degree_nonTargeted,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(In_degree_nonTargeted,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('In-degree nonInitiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




In_t_nt = In_degree_Targeted./In_degree_nonTargeted;

ax = figure(20)

b = imagesc(In_t_nt);
set(b,'AlphaData',~isnan(In_t_nt)) %~isnan(targeted)
% caxis manual
% caxis([0 10]); 

colormap(ax,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(In_t_nt,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(In_t_nt,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Ratio Fig1/Fig2') %Targeted / non-targeted

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





%% Out-degree


Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

Out_degree_Targeted = flipud(Out_degree_Targeted).*Zero_diag;

ax1 = figure(3)
b = imagesc(Out_degree_Targeted);
set(b,'AlphaData',~isnan(Out_degree_Targeted)) %~isnan(idx)
% caxis manual
% caxis([10 50]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(Out_degree_Targeted,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(Out_degree_Targeted,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Out-degree Initiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square




Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

Out_degree_nonTargeted = flipud(Out_degree_nonTargeted).*Zero_diag;

ax1 = figure(4)
b = imagesc(Out_degree_nonTargeted);
set(b,'AlphaData',~isnan(Out_degree_nonTargeted)) %~isnan(idx)
%caxis manual
%caxis([-10 45]); 

colormap(ax1,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(Out_degree_nonTargeted,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(Out_degree_nonTargeted,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Out-degree nonInitiator')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square





Out_t_nt = Out_degree_Targeted./Out_degree_nonTargeted;

ax = figure(40)

b = imagesc(Out_t_nt);
set(b,'AlphaData',~isnan(Out_t_nt)) %~isnan(targeted)
% caxis manual
% caxis([0 10]); 

colormap(ax,cool)%
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(Out_t_nt,2)])
xticklabels(p_boost)
xlabel('p_{wave}')
    
yticks([1:size(Out_t_nt,1)])
yticklabels(fliplr(p_dist))
ylabel('p_{proximity}')

title('Ratio Fig3/Fig4') %Targeted / non-targeted

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square

