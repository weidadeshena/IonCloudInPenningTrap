% a random script to plot some data for my project report
% ignore

% alpha = [0.01,0.05,0.1,0.15,0.2,0.3,0.45,0.6,0.85,0.99,1.05,1.2,1.6,2.2,3,5,10,50,200];
% N = size(alpha,2);
% 
% max_r = zeros(N,1);
% max_z = zeros(N,1);
% for i = 1:N
%     [omega_z_squared,omega_1_squared] = Greens(alpha(i));
%     final_positions = crystal_graphs_energy(50,omega_z_squared,omega_1_squared,0,1);
%     max_r(i) = max(final_positions(:,end));
%     max_z(i) = max(final_positions(:,3));
% end
    

% alpha = [0.2,1/3,0.5,0.75,0.999,4/3,2,3,5];

% alpha = 3;
% 
alpha = 0.1
[omega_z_squared,omega_1_squared] = Greens(alpha);
[final_positions] = crystal_graphs_energy(300,omega_z_squared,omega_1_squared,1,0);


% N_list = [10,20,30,40,50,60,80,100,150,200];
% s = size(N_list,2);
% coulomb_energy = zeros(s,1);
% trapping_energy = zeros(s,1);
% total_energy = zeros(s,1);
% for i = 1:s
%     [final_positions,energy_c,energy_t,total_energy] = crystal_graphs_energy(N_list(i),1/3,1/3,0,1);
%     coulomb_energy(i) = energy_c;
%     trapping_energy(i) = energy_t;
%     total_energy(i) = total_energy;
% end



