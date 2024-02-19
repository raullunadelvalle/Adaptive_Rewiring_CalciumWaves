clear all
close all
clc


% Use this script to generate a figure similar to Figure 5C (lower panel) of:
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

fid = py.open('all_cd_pairs_15.pckl','rb');  
cd_pairs = py.pickle.load(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 0;
for i = p_dist
    ii = ii+1;
    jj = 0;
    for j = p_boost
        jj = jj+1;
        if (i+j) <= 1 
            rand_num_conv_div_k = [];
            kk = 0;
            for k = 0:K-1
                kk = kk+1;
                rand_num_conv_div_r = [];
                rr = 0;
                for r = 0:R-1
                    rr = rr+1;
                    node = boosted_nodes{k};
                    node = double(node{1});
        
                    conv_div = int64(cd_pairs{tau}{0.5}{i}{j}{k}{r});

                    
                    if isempty(conv_div(2,:)) == 1
                        rand_num_conv_div_r(rr) = NaN;
                    else

                    
                    %random nodes. How often do they become divergent nodes in c-ds?
                    rand_num_conv_div_g = [];
                    gg = 0;
                    for g = 1:1000
                        gg = gg+1;
                        rand_node = node;
                        while rand_node == node
                            rand_node = randi([0 99]);
                        end


                        if ismember(rand_node,conv_div(2,:)) == 1
                            rand_num_conv_div = 1;
                        else
                            rand_num_conv_div = 0;
                        end

                        rand_num_conv_div_g(gg) = rand_num_conv_div;

                    end

                    rand_num_conv_div = nansum(rand_num_conv_div_g)/gg;
                    rand_num_conv_div_r(rr) = rand_num_conv_div;
                    end 
                    
                end


                rand_num_conv_div = nanmean(rand_num_conv_div_r);  %nansum RETOMAR AQUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                rand_num_conv_div_k(kk) = rand_num_conv_div;

            end

            % Number Convergent divergent unit that includes random node
            rand_num(ii,jj) = nansum(rand_num_conv_div_k) / (kk-sum(isnan(rand_num_conv_div_k))); %nanmean(rand_num_conv_div_k);
        end
    end
end




Zero_diag = tril( ones(length(p_dist),length(p_boost)) ); 
Zero_diag(Zero_diag==0) = NaN;

rand_num = flipud(rand_num/1).*Zero_diag;



figure(1)
b = imagesc(rand_num);
set(b,'AlphaData',~isnan(rand_num)) 
caxis manual
caxis([0 0.7]); 

colormap cool
cob = colorbar;
c.LineWidth = 3;

xticks([1:size(rand_num,2)])
xticklabels(p_boost)
xlabel('p Wave')
    
yticks([1:size(rand_num,1)])
yticklabels(fliplr(p_dist))
ylabel('p Distance')

title('c-d unit rand node success')

set(gca,'linewidth',2)
set(gca,'FontSize',20)
set(gcf,'color','w');
box off

axis square
