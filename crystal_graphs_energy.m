function [final_positions] = crystal_graphs_energy(ions,om_z_squared,om_1_squared,graph,energy)

% ion_positions = crystals(ions,om_z_squared,om_1_squared);
final_positions = crystals(ions,om_z_squared,om_1_squared);

% final_positions = reshape(ion_positions(end,:),3,size(ion_positions,2)/3)'; 
r = sqrt(final_positions(:,1).^2 + final_positions(:,2).^2 );

final_positions = horzcat(final_positions,r);
if graph
    disp('plotting graph')
    min_axis = min(min(final_positions(:,:)));
    max_axis = max(max(final_positions(:,:)));

    figure(1)
%     hold off;

    plot3(final_positions(:,1),final_positions(:,2),final_positions(:,3),'b*')


    axis([min_axis max_axis min_axis max_axis min_axis max_axis]);
    axis square
    grid on
    view(0,0)


    hold off

    figure(2)
    [gaussian_filtered_image,simulated_image] = simulate_image(final_positions(:,4),final_positions(:,3));
    
end
if energy
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
        total_trap_energy = total_trap_energy + 0.5 * (om_1_squared * (final_positions(ion_number_1,1)^2 + ...
        final_positions(ion_number_1,2)^2) + om_z_squared * final_positions(ion_number_1,3)^2);

    end
%     hold on
    total_energy = total_coulomb_energy+ total_trap_energy;
    total_energy
end


%axis square


function [filtered_image,sim_image] = simulate_image(amplitude,z_pos);
%Make a 2d image where the photo lives
pixel_unit_scaling = 3.8;
image_size = 101;%round(2* max(max(pixel_unit_scaling*amplitude),max(pixel_unit_scaling*abs(z_pos)))) + 10;

sim_image = zeros(round(image_size));
image_parts = sim_image;
extra_ion = sim_image;
distance_from_centre = (1:image_size) - image_size/2 - 0.5;

%amplitude = 20;
pixel_size = 1;
%start with one particle oscillating with amplitude = 20
%ion_positions = zeros(size(amplitude,1),5);

% ion_positions(:,1)=51 - pixel_unit_scaling .* z_pos;
% ion_positions(:,2)=floor(51 - pixel_unit_scaling .* z_pos);
% ion_positions(:,3)=1 + floor(51-pixel_unit_scaling.*z_pos) - (51-pixel_unit_scaling.*z_pos);
% ion_positions(:,4)=ceil(51-pixel_unit_scaling.*z_pos);
% ion_positions(:,5)=( - floor(51-pixel_unit_scaling.*z_pos) + (51-pixel_unit_scaling.*z_pos));

for k = 1:size(amplitude,1)

    
    if amplitude(k) < 0.001 
        
%         image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),50) = ...
%         (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) * 0.5;
%         
%         sim_image = sim_image + image_parts;
% 
%         image_parts(:,:) = 0;
        
        
        image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),51) = ...
        1- ((51-pixel_unit_scaling*z_pos(k)) - floor(51-pixel_unit_scaling*z_pos(k)) );   
    

    
        sim_image = sim_image + image_parts;

        image_parts(:) = 0;
        
        
%         image_parts(floor(51-pixel_unit_scaling*z_pos(k)),50) = ...
%         (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) * 0.5;
%         
%         sim_image = sim_image + image_parts;
% 
%         image_parts(:,:) = 0;
        
        
        image_parts(ceil(51-pixel_unit_scaling*z_pos(k)),51) = ...
        ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k)));     
    
        
        
        sim_image = sim_image + image_parts;

        image_parts(:) = 0;
        
    else
        
        
        image_parts(floor(51 - pixel_unit_scaling * z_pos(k)),:) = ...
        (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k))) *...
        real( asin((distance_from_centre + pixel_size/2)/(pixel_unit_scaling*amplitude(k)))...
        - asin((distance_from_centre - pixel_size/2)/(pixel_unit_scaling*amplitude(k))) ) / pi;    

%   (1 + floor(51-pixel_unit_scaling*z_pos(k)) - (51-pixel_unit_scaling*z_pos(k)))% *...

        extra_ion = extra_ion + image_parts;
        
        image_parts(:) = 0;



        image_parts(ceil(51-pixel_unit_scaling*z_pos(k)),:) = ...
        ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k))) *...
        real( asin((distance_from_centre + pixel_size/2)/(pixel_unit_scaling*amplitude(k)))...
        - asin((distance_from_centre - pixel_size/2)/(pixel_unit_scaling*amplitude(k))) ) / pi;   

    
 %  ( - floor(51-pixel_unit_scaling*z_pos(k)) + (51-pixel_unit_scaling*z_pos(k))) %* 


        extra_ion = extra_ion + image_parts;

        image_parts(:) = 0;
        sum(extra_ion(:));
        %extra_ion = 1/sum(extra_ion(:));
        
        sim_image = sim_image + extra_ion;
        extra_ion(:) = 0;
    end
    
end

%normalising image matrix

%sum(sim_image(:));

%tic
gaussian_filter = fspecial('gaussian',5,1.1);
filtered_image = imfilter(sim_image,gaussian_filter);
toc

% colorbar
% figure(1)
% imagesc(sim_image);
% colorbar
figure(2)
imagesc(filtered_image);
% colorbar
% 
% amplitude
% ion_positions


function [positions] = crystals(num_of_ions,omegaz_squared,omega1_squared)

%Write this initial script for n ions
%Put constants here.
%num_of_ions = 3;
num_of_time_steps = 10000;
%omega1_squared = 1-0.5*omegaz_squared;
time_step_size = 0.1;

% ion_forces = zeros(num_of_ions,num_of_ions*3);
r_cubed = zeros(num_of_ions, num_of_ions);
trap_forces = zeros(num_of_ions,3);
% total_force = zeros(num_of_time_steps,3*num_of_ions);


positions = rand(num_of_ions,3);


tic
for time_step = 1:num_of_time_steps

coulomb_forces = zeros(num_of_ions,3);
for ion_number = 1:num_of_ions
    for k = 1:num_of_ions
        r_cubed(ion_number,k) = sum((positions(ion_number,:)-positions(k,:)).^2)^1.5;
    end 
end



for ion_number = 1:num_of_ions
    
    for l = 1:num_of_ions
    
        if l ~= ion_number
  
            coulomb_forces(ion_number,1) = coulomb_forces(ion_number,1) + (-positions(l,1) +...
            positions(ion_number,1))/r_cubed(ion_number,l); 
            %force from the x-direction
            coulomb_forces(ion_number,2) = coulomb_forces(ion_number,2) + (-positions(l,2) +...
            positions(ion_number,2))/r_cubed(ion_number,l); 
            %force from the y-direction
            coulomb_forces(ion_number,3) = coulomb_forces(ion_number,3) + (-positions(l,3) +...
            positions(ion_number,3))/r_cubed(ion_number,l); 
            %force from the z-direction
            
        end
    
    end
        
end


% trap_forces = -omega1_squared*positions;
% trap_forces(:,3) = -omegaz_squared*positions(:,3);

for ion_number = 1:num_of_ions
    trap_forces(ion_number,1) = -omega1_squared * positions(ion_number,1);
    trap_forces(ion_number,2) = -omega1_squared * positions(ion_number,2);
    trap_forces(ion_number,3)   = -omegaz_squared * positions(ion_number,3);
end


positions = positions + time_step_size*(coulomb_forces+trap_forces);

% for ion_number = 1:num_of_ions
%    
%         positions(ion_number,1) = positions (ion_number,1) + time_step_size*(coulomb_forces(ion_number,1) + trap_forces(ion_number,1));
%         positions(ion_number,2) = positions (ion_number,2) + time_step_size*(coulomb_forces(ion_number,2) + trap_forces(ion_number,2));
%         positions(ion_number,3) = positions (ion_number,3) + time_step_size*(coulomb_forces(ion_number,3) + trap_forces(ion_number,3));
%         
%         
% 
% %         total_force(time_step,ion_number*3-2) = sum(ion_forces(1:3:end,ion_number)) + trap_forces(3*ion_number-2);
% %         total_force(time_step,ion_number*3-1) = sum(ion_forces(2:3:end,ion_number)) + trap_forces(3*ion_number-1);
% %         total_force(time_step,ion_number*3) = sum(ion_forces(3:3:end,ion_number)) + trap_forces(3*ion_number);
% end

end
toc

