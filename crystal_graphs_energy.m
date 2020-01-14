function [final_positions] = crystal_graphs_energy(ions,om_z_squared,om_1_squared)

ion_positions = crystals(ions,om_z_squared,om_1_squared);

final_positions = reshape(ion_positions(end,:),3,size(ion_positions,2)/3)'; 
final_positions(:,4) = sqrt(final_positions(:,1).^2 + final_positions(:,2).^2 );

min_axis = min(min(final_positions(:,:)));
max_axis = max(max(final_positions(:,:)));

figure(1)
hold off;
for ion_number = 1:ions
    
    plot3(ion_positions(end,ion_number*3-2),ion_positions(end,ion_number*3-1),ion_positions(end,3*ion_number),'b*');
    
    hold on
end



axis([min_axis max_axis min_axis max_axis min_axis max_axis]);
axis square
grid on
view(0,0)

% plot3(final_positions(:,4),zeros(size(final_positions,1)),final_positions(:,3),'ro')

hold off

final_positions(:,[4 3]);

total_coulomb_energy=0;
total_trap_energy=0;


for ion_number_1 = 1:ions-1
for ion_number_2 = ion_number_1+1:ions
    
    xxx = ( final_positions(ion_number_1,1) - final_positions(ion_number_2,1))^2;
        yyy = ( final_positions(ion_number_1,2) - final_positions(ion_number_2,2))^2;
        zzz = ( final_positions(ion_number_1,3) - final_positions(ion_number_2,3))^2;
        ion_ion_distance = (xxx + yyy + zzz)^(0.5);
    
        
    total_coulomb_energy = total_coulomb_energy +  1/ion_ion_distance;   
    
end
end

for ion_number_1 = 1:ions
    total_trap_energy = total_trap_energy + 0.5 * (om_1_squared * (final_positions(ion_number_1,1)^2 + final_positions(ion_number_1,2)^2) + om_z_squared * final_positions(ion_number_1,3)^2);

    
end
hold on
total_energy = total_coulomb_energy+ total_trap_energy;
total_energy


%figure(2)
% [gaussian_filtered_image,simulated_image] = simulate_image(final_positions(:,4),final_positions(:,3));
% %axis square
% 
% 
% function [filtered_image,sim_image] = simulate_image(amplitude,z_pos);
% %Make a 2d image where the photo lives
% pixel_unit_scaling = 3.8;
% image_size = 101;%round(2* max(max(pixel_unit_scaling*amplitude),max(pixel_unit_scaling*abs(z_pos)))) + 10;
% 
% sim_image = zeros(round(image_size));
% image_parts = sim_image;
% extra_ion = sim_image;
% distance_from_centre = (1:image_size) - image_size/2 - 0.5;
% 
% %amplitude = 20;
% pixel_size = 1;
% %start with one particle oscillating with amplitude = 20
% %ion_positions = zeros(size(amplitude,1),5);
% 
% % ion_positions(:,1)=51 - pixel_unit_scaling .* z_pos;
% % ion_positions(:,2)=floor(51 - pixel_unit_scaling .* z_pos);
% % ion_positions(:,3)=1 + floor(51-pixel_unit_scaling.*z_pos) - (51-pixel_unit_scaling.*z_pos);
% % ion_positions(:,4)=ceil(51-pixel_unit_scaling.*z_pos);
% % ion_positions(:,5)=( - floor(51-pixel_unit_scaling.*z_pos) + (51-pixel_unit_scaling.*z_pos));
% 
% for k = 1:size(amplitude,1)
% 
%     
%     if amplitude(k) < 0.001 
%         
% %         image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),50) = ...
% %         (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) * 0.5;
% %         
% %         sim_image = sim_image + image_parts;
% % 
% %         image_parts(:,:) = 0;
%         
%         
%         image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),51) = ...
%         1- ((51-pixel_unit_scaling*z_pos(k)) - floor(51-pixel_unit_scaling*z_pos(k)) );   
%     
% 
%     
%         sim_image = sim_image + image_parts;
% 
%         image_parts(:) = 0;
%         
%         
% %         image_parts(floor(51-pixel_unit_scaling*z_pos(k)),50) = ...
% %         (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) * 0.5;
% %         
% %         sim_image = sim_image + image_parts;
% % 
% %         image_parts(:,:) = 0;
%         
%         
%         image_parts(ceil(51-pixel_unit_scaling*z_pos(k)),51) = ...
%         ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k)));     
%     
%         
%         
%         sim_image = sim_image + image_parts;
% 
%         image_parts(:) = 0;
%         
%     else
%         
%         
%         image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),:) = ...
%         (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) *...
%         real( asin((distance_from_centre + pixel_size/2)/(pixel_unit_scaling*amplitude(k)))...
%         - asin((distance_from_centre - pixel_size/2)/(pixel_unit_scaling*amplitude(k))) ) / pi;    
% 
% %   (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k)))% *...
% 
%         extra_ion = extra_ion + image_parts;
%         
%         image_parts(:) = 0;
% 
% 
% 
%         image_parts(ceil(51-pixel_unit_scaling*z_pos(k)),:) = ...
%         ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k))) *...
%         real( asin((distance_from_centre + pixel_size/2)/(pixel_unit_scaling*amplitude(k)))...
%         - asin((distance_from_centre - pixel_size/2)/(pixel_unit_scaling*amplitude(k))) ) / pi;   
% 
%     
%  %  ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k))) %* 
% 
% 
%         extra_ion = extra_ion + image_parts;
% 
%         image_parts(:) = 0;
%         sum(extra_ion(:));
%         %extra_ion = 1/sum(extra_ion(:));
%         
%         sim_image = sim_image + extra_ion;
%         extra_ion(:) = 0;
%     end
%     
% end
% 
% %normalising image matrix
% 
% %sum(sim_image(:));
% 
% %tic
% gaussian_filter = fspecial('gaussian',5,1.1);
% filtered_image = imfilter(sim_image,gaussian_filter);
%toc

%colorbar
%figure(1)
%imagesc(sim_image);
%colorbar
%figure(2)
%imagesc(filtered_image);
%colorbar

% amplitude
% ion_positions


function [positions] = crystals(num_of_ions,omegaz_squared,omega1_squared)

%Write this initial script for n ions
%Put constants here.
%num_of_ions = 3;
num_of_time_steps = 10000;
%omega1_squared = 1-0.5*omegaz_squared;
time_step_size = 0.1;

ion_forces = zeros(num_of_time_steps,num_of_ions*3,num_of_ions);
positions = zeros(num_of_time_steps,3*num_of_ions);
ion_ion_distance_cubed = zeros(num_of_ions, num_of_ions);
trap_forces = zeros(num_of_time_steps,3*num_of_ions);
total_force = zeros(num_of_time_steps,3*num_of_ions);

positions(1,:) = rand(1,3*num_of_ions);

tic

for time_step = 1:num_of_time_steps;


for ion_number = 1:num_of_ions
    
    for k = 1:num_of_ions
        xxx = ( positions(time_step,3*ion_number-2) - positions(time_step,3*k -2))^2;
        yyy = ( positions(time_step,3*ion_number-1) - positions(time_step,3*k -1))^2;
        zzz = ( positions(time_step,3*ion_number) - positions(time_step,3*k))^2;
        ion_ion_distance_cubed(ion_number,k) = (xxx + yyy + zzz)^1.5;
    end
    
end



for ion_number = 1:num_of_ions
    
    for l = 1:num_of_ions
    
        if l ~= ion_number
        
            ion_forces(time_step,l*3-2,ion_number) = (-positions(time_step,3*l -2) + positions(time_step,3*ion_number-2))...
            /ion_ion_distance_cubed(ion_number,l);
            %force from the x-direction
    
            ion_forces(time_step,l*3-1,ion_number) = (-positions(time_step,3*l -1) + positions(time_step,3*ion_number-1))...
            /ion_ion_distance_cubed(ion_number,l);
            %force from the y-direction
    
            ion_forces(time_step,l*3,ion_number) = (-positions(time_step,3*l) + positions(time_step,3*ion_number))...
            /ion_ion_distance_cubed(ion_number,l);
            %force from the z-direction
            
        end
    
    end
        
end

for ion_number = 1:num_of_ions
    
    trap_forces(time_step,ion_number*3-2) = -omega1_squared * positions(time_step,3*ion_number-2);
    trap_forces(time_step,ion_number*3-1) = -omega1_squared * positions(time_step,3*ion_number-1);
    trap_forces(time_step,ion_number*3)   = -omegaz_squared * positions(time_step,3*ion_number);
    
    
end



for ion_number = 1:num_of_ions
   
        positions(time_step+1,3*ion_number-2) = positions (time_step,3*ion_number-2) + time_step_size*(sum(ion_forces(time_step,1:3:end,ion_number)) + trap_forces(time_step,3*ion_number-2));
        positions(time_step+1,3*ion_number-1) = positions (time_step,3*ion_number-1) + time_step_size*(sum(ion_forces(time_step,2:3:end,ion_number)) + trap_forces(time_step,3*ion_number-1));
        positions(time_step+1,3*ion_number) = positions (time_step,3*ion_number) + time_step_size*(sum(ion_forces(time_step,3:3:end,ion_number)) + trap_forces(time_step,3*ion_number));

        total_force(time_step,ion_number*3-2) = sum(ion_forces(time_step,1:3:end,ion_number)) + trap_forces(time_step,3*ion_number-2);
        total_force(time_step,ion_number*3-1) = sum(ion_forces(time_step,2:3:end,ion_number)) + trap_forces(time_step,3*ion_number-1);
        total_force(time_step,ion_number*3) = sum(ion_forces(time_step,3:3:end,ion_number)) + trap_forces(time_step,3*ion_number);
end

end


toc