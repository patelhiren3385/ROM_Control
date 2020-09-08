clear all
close all
clc
addpath('lib');
%% Normal distribution
nd.variance = 25 ; 
nd.mu = 0;
nd.sigma = sqrt(nd.variance);
nd.n = 5 ; 
distribution.delta = 1:1:nd.n ; 
%% Monte Carlo Simulation
for i = 1:length(distribution.delta)
    MC_s = 'Monte Carlo iteration no. %f \n' ; 
    fprintf(MC_s,i)
    %% Dimensions
    dim.len = 0.5 ; 
    dim.length = dim.len ; % in m
    dim.width = 0.03; % in m
    dim.depth = 0.002; % in m
    dim.support_condition = 'c'; % 'c' for cantilever
    %% Material
    % steel
    E = 210; % in Pa
    percentage_of_variation = 0.30 ; 
    Mean_value = E ; 
    no_of_points = nd.n ;
    mc.E = get_log_normal_distribution(percentage_of_variation,Mean_value,no_of_points) ; 
    nu = 0.3;
    rho = 7800; % in kg/m^3
    cdr = 0.0002; % critical damping ratio
    material = get_mechanical_properties(1e9.*mc.E.Emc(i), nu, rho, dim.width, -dim.depth/2, dim.depth/2);

    %% Actuator
    % piezo
    E_piezo = 139e9; % in Pa
    nu_piezo = 0.3;
    rho_piezo = 7500; % in kg/m^3

    pzt_depth = 0.2e-3; % in m

    piezo_act = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, dim.depth/2, dim.depth/2 + pzt_depth);
    piezo_sen = get_mechanical_properties(E_piezo, nu_piezo, rho_piezo, dim.width, -dim.depth/2-pzt_depth, -dim.depth/2);
    piezo_act.piezoelectric_constant = 15.29; % Cm^-2
    piezo_sen.piezoelectric_constant = 15.29; % Cm^-2
    piezo_act.dielectric_constant = 11e-9; % Fm^-2
    piezo_sen.dielectric_constant = 11e-9; % Fm^-2
    %% Element
    dof_per_node = 2;
  
    elements.sec1 = 6 ; % 4 cm - 0.08 element size
    elements.sec2 = 19 ; % 8 cm - 0.08 element size
    elements.sec3 = 1 ; % 2 cm - 0.4 element size
    elements.sec4 = 19 ; % 8 cm - 0.16 element size
    elements.sec5 = 27 ; % 16 cm - 0.8 element size
    elements.sec6 = 19 ; % 8 cm - 0.8 element size
    elements.sec7 = 9 ; % 4 cm - 0.8 element size
    
    
    elements.beam = elements.sec7 + elements.sec6 + elements.sec5 + elements.sec4 + elements.sec3 + elements.sec2 + elements.sec1 ;
    %% Geometry
    % 100 - Elements
    dim.length_sec1 = 0.03 ;
    dim.length_sec2 = 0.095 ;
    dim.length_sec3 = 0.005 ;
    dim.length_sec4 = 0.095 ;
    dim.length_sec5 = 0.1350 ;
    dim.length_sec6 = 0.095 ;
    dim.length_sec7 = 0.0450 ;
    
    
    beam = get_nodes_coords_connectivity(dim, elements, dof_per_node);
       
    force = [0;0];
    point_load = 0;

    wdof = beam.dofs(1,:); %Transverse deflection
    txdof = beam.dofs(2,:); %Theta_x

   
    element_act1 = elements.sec1 + 1 : elements.sec1 + elements.sec2 ; % 11 - 210
    element_act2 = element_act1(end) + elements.sec3 + 1 :elements.sec1 + elements.sec2 + elements.sec3 + elements.sec4 ; % 221 - 320
    element_act3 = element_act2(end) + elements.sec5 + 1 : elements.sec1 + elements.sec2 + elements.sec3 + elements.sec4 + elements.sec5 + elements.sec6  ; % 341 - 360
    
    % seperate actuators while assembly
    element_sen1 = element_act1;
    element_sen2 = element_act2;
    element_sen3 = element_act3;
    n_act = 3;
    n_sen = 3;

    %% Matrix construction
    K_global = zeros(beam.total_dofs, beam.total_dofs);
    M_global = zeros(beam.total_dofs, beam.total_dofs);
    F_global = zeros(beam.total_dofs,3); %External Load
    P_global = zeros(beam.total_dofs,n_act); %Actuator load
    C_global = zeros(n_sen,beam.total_dofs); %y = c*x
    K_beam = zeros(beam.total_dofs, beam.total_dofs);
    M_beam = zeros(beam.total_dofs, beam.total_dofs);

    for element = 1:beam.total_elements
        dof_address = node2dof(beam.connectivity(element,:),dof_per_node);
        el_connect = dof_address(:) ; 
        x1 = beam.element_coordinates(element,1) ; 
        x2 = beam.element_coordinates(element,2) ; 

        element_matrices = get_element_matrices(material.D, material.rho, x1,x2);

        K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + element_matrices.stiffness;
        M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + element_matrices.mass;
        F_global(el_connect,:) = F_global(el_connect,:) + element_matrices.force;

        K_beam(el_connect, el_connect) = K_beam(el_connect, el_connect) + element_matrices.stiffness;
        M_beam(el_connect, el_connect) = M_beam(el_connect, el_connect) + element_matrices.mass;

        
        if any(element_act1 == element)
            pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.piezoelectric_constant,beam.width, x1, x2, piezo_act.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            P_global(el_connect,1) = P_global(el_connect,1) + pzt_matrices.force;
        end
        
        if any(element_act2 == element)
            pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.piezoelectric_constant,beam.width, x1, x2, piezo_act.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            P_global(el_connect,2) = P_global(el_connect,2) + pzt_matrices.force;
        end
        
        if any(element_act3 == element)
            pzt_matrices = get_actuator_matrices(piezo_act.D, piezo_act.rho, piezo_act.piezoelectric_constant,beam.width, x1, x2, piezo_act.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            P_global(el_connect,3) = P_global(el_connect,3) + pzt_matrices.force;
        end
        
        
        if any(element_sen1 == element)
            pzt_matrices = get_sensor_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.piezoelectric_constant,piezo_sen.dielectric_constant,pzt_depth, x1, x2, piezo_sen.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            C_global(1,el_connect) = C_global(1,el_connect) + pzt_matrices.force;
        end
       
        if any(element_sen2 == element)
            pzt_matrices = get_sensor_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.piezoelectric_constant,piezo_sen.dielectric_constant,pzt_depth, x1, x2, piezo_sen.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            C_global(2,el_connect) = C_global(2,el_connect) + pzt_matrices.force;
        end
        
        if any(element_sen3 == element)
            pzt_matrices = get_sensor_matrices(piezo_sen.D, piezo_sen.rho, piezo_sen.piezoelectric_constant,piezo_sen.dielectric_constant,pzt_depth, x1, x2, piezo_sen.lever_arm);
            K_global(el_connect, el_connect) = K_global(el_connect, el_connect) + pzt_matrices.stiffness;
            M_global(el_connect, el_connect) = M_global(el_connect, el_connect) + pzt_matrices.mass;
            C_global(3,el_connect) = C_global(3,el_connect) + pzt_matrices.force;
        end
        
        
    end
    K = K_global(beam.free_dofs, beam.free_dofs);
    M = M_global(beam.free_dofs, beam.free_dofs);
    F = F_global(beam.free_dofs,:);
    P = P_global(beam.free_dofs,:);
    C = C_global(:,beam.free_dofs);

    Kb = K_beam(beam.free_dofs, beam.free_dofs);
    Mb = M_beam(beam.free_dofs, beam.free_dofs);
    
    F = zeros(length(F),1) ; 
    F(end-1) = 1 ; 
    %% Modal analysis
    [eigenvectors, eigenvalues] = eig(K,M);
    [eigenvalues, indices] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:,indices);
    frequencies = sqrt(eigenvalues);
    w1 = frequencies(1)/(2*pi); % in HZ
    w2 = frequencies(2)/(2*pi);
    bet = 2*cdr/(w1+w2);
    alp = w1*w2*bet;
    
    D = alp*M + bet*K; %Dmaping matrix
    freq = frequencies./(2*pi) ; 
    %% SEREP
    %%

%     rom.config = 'Best selection for 50 elements' ;
%     rom.act_1 = [7 ,8 ,9 ,10,11,13,14,15,17,18,19,21,22,23,24] ; %7:24
%     rom.act_2 = [27,28,29,30,31,33,34,35,36,37,38,39,40,41,42,43,44] ; %27:44
%     rom.act_3 = [73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90] ; %73:90
%     rom.sec_1 = [1,2,3,4,5,6] ; %1:6
%     rom.sec_2 = [25,26] ; %25,26
%     rom.sec_3 = [45,46,47,49,50,51,53,54,55,57,58,59,61,62,63,65,66,67,68,69,71,72] ; %45:72
%     rom.sec_4 = [91,92,93,95,96,97,99,100] ; %91:100
%     rom.dof = [rom.sec_1,rom.act_1,rom.sec_2,rom.act_2,rom.sec_3,rom.act_3,rom.sec_4] ; 
% 
%     rom.config = 'Best selection for 100 elements' ;
%     rom.act_1 = [ 13, 14, 15, 16, 19, 20, 23, 24, 27, 28, 29, 30, 33, 34, 36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50] ; %13:50
%     rom.act_2 = [ 53, 54, 55, 56, 59, 60, 63, 64, 67, 68, 69, 70, 73, 74, 76, 77, 78, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90] ; %53:90
%     rom.act_3 = [145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182] ; %145:182
%     rom.sec_1 = [1,2,3,4,5,7,8,9,10,11,12] ; %1:12
%     rom.sec_2 = [51,52] ; %51,52
%     rom.sec_3 = [ 91, 92, 93, 94, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,...
%                  123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144] ; %91:144 
%     rom.sec_4 = [183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200] ; %183:200 
    
%     rom.dof = [rom.sec_1,rom.act_1,rom.sec_2,rom.act_2,rom.sec_3,rom.act_3,rom.sec_4] ; 
    
%     rom.dof = [1,40,80,120,160,200] ; % 5 Element
%     rom.dof = [1,8,16,24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,160,168,176,184,192,200] ; %25 element
%     rom.dof = [1,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,...,
%                 116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200] ; %50 element
    rom.dof = [1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,...,
             64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,...,
             120,122,124,126,128,130,132,134,136,138,140,142,144,146,148,150,152,154,156,158,160,162,164,166,168,...,
             170,172,174,176,178,180,182,184,186,188,190,192,194,196,198,200] ; % 100 element
    
    rom.config = 'Hundred_dofs_selected' ;
    rom.sec_1 = [1:(element_act1(1)-1)*2 ] ; % 
    rom.sec_2 = [ element_act1(1)*2 - 1 : element_act1(end)*2 ] ; %Actuator 1 
    rom.sec_3 = [rom.sec_2(end)+1 : (element_act2(1)-1)*2] ; 
    rom.sec_4 = [ rom.sec_3(end)+1 : element_act2(end)*2] ; %Actuator 2 
    rom.sec_5 = [rom.sec_4(end)+1:(element_act3(1)-1)*2] ; 
    rom.sec_6 = [ rom.sec_5(end) + 1 : element_act3(end)*2] ; %Actuator 3 
    rom.sec_7 = [ rom.sec_6(end) + 1 : elements.beam*2 ] ;
    rom.sec_detail = [rom.sec_1(1),rom.sec_1(end) ; rom.sec_2(1),rom.sec_2(end) ; rom.sec_3(1),rom.sec_3(end) ; rom.sec_4(1),rom.sec_4(end) ; rom.sec_5(1),rom.sec_5(end) ; rom.sec_6(1),rom.sec_6(end) ; rom.sec_7(1),rom.sec_7(end) ] ; 
    
    %%
    rom.eig_vec = eigenvectors(rom.dof,rom.dof) ; 
    rom.gen_inv = (rom.eig_vec'*rom.eig_vec)\rom.eig_vec' ;
    rom.T_serep = eigenvectors(:,1:length(rom.dof))*rom.gen_inv ; 
    
    rom.Mb = rom.T_serep'*Mb*rom.T_serep ;
    rom.Kb = rom.T_serep'*Kb*rom.T_serep ;
    rom.M = rom.T_serep'*M*rom.T_serep ;
    rom.K = rom.T_serep'*K*rom.T_serep ;
    rom.D = rom.T_serep'*D*rom.T_serep ;
    rom.C = C*rom.T_serep ; 
    rom.F = rom.T_serep'*F ; 
    rom.P = rom.T_serep'*P ;
    
    [rom.eigenvectors, rom.eigenvalues] = eig(rom.K,rom.M);
    [rom.eigenvalues, rom.indices] = sort(diag(rom.eigenvalues));
    rom.eigenvectors = rom.eigenvectors(:,rom.indices);
    rom.frequencies = sqrt(rom.eigenvalues);
    rom.freq = rom.frequencies./(2*pi) ; 
    Mode_ROMt = rom.T_serep*rom.eigenvectors ;
     
    %% Dynamics - Full Order
    n = length(beam.free_dofs);
    
    Asys = [zeros(n,n), eye(n);-M\K, -M\D];
    Bext = [zeros(n,1); M\F];
    Bcont = [zeros(n,n_act); M\P];
    Csys = [C, zeros(n_sen,n)];
    sys_full = ss(Asys,Bcont,Csys,zeros(size(Csys,1),size(Bcont,2))) ; 
    x0 = zeros(n,1);
    xd0 = zeros(n,1);
    y0 = [x0;xd0];
    %% Dynamics - ROM using SEREP
    rom.n = length(rom.dof) ; 
    
    rom.Asys = [zeros(rom.n,rom.n) , eye(rom.n) ; -rom.M\rom.K , -rom.M\rom.D ] ; 
    rom.Bext = [zeros(rom.n,1) ; rom.M\rom.F] ;
    rom.Bcont = [zeros(rom.n,n_act) ; rom.M\rom.P] ;
    rom.Csys = [rom.C , zeros(n_sen,rom.n)] ; 
    
    rom.x0 = zeros(rom.n,1) ; 
    rom.xd0 = zeros(rom.n,1) ; 
    rom.y0 = [ rom.x0 ; rom.xd0 ] ; 
    
    sys_rom = ss(rom.Asys,rom.Bcont,rom.Csys,zeros(size(rom.Csys,1),size(rom.Bcont,2))) ;
    %% Simulation time
    tn = 100;
    tf = 10;
    %tspan = linspace(0,tf,tn);
    %% LQR - Controller
    %Full Order
    lqr_c.Q = Csys'*Csys ;
    lqr_c.R = (0.01.*(Bcont'*Bcont))\eye(n_act) ;
    
    f_lqr_st = cputime ;
    [lqr_c.fo_gain,~,~] = lqr(sys_full,lqr_c.Q,lqr_c.R);
    f_lqr_et = cputime - f_lqr_st ;
    
    %Reduced Order
    lqr_c.rom_Q = rom.Csys'*rom.Csys;
    lqr_c.rom_R = (0.01.*(rom.Bcont'*rom.Bcont))\eye(n_act) ;
    
    r_lqr_st = cputime ;
    [lqr_c.rom_gain,~,~] = lqr(sys_rom,lqr_c.rom_Q,lqr_c.rom_R);
    r_lqr_et = cputime - r_lqr_st ; 
    
    time_saved = 1 - (r_lqr_et/f_lqr_et) ; 
    per_time_saved = zeros(length(distribution.delta),1) ;
    per_time_saved(i) = time_saved ; 
    
    %Mapping ROM gain to Full order gain
    lqr_c.mfull_gain = [lqr_c.rom_gain(:,1:length(rom.dof))*rom.T_serep',lqr_c.rom_gain(:,length(rom.dof)+1:end)*rom.T_serep'] ; 
    %% ODE solution - LQR regulator problem
%--------------------------------------------------------------------------
%Dynamic solution using - lsim
%--------------------------------------------------------------------------
    statespace.ydot_full_uc = ss(Asys,Bext,Csys,zeros(size(Csys,1),size(Bext,2))) ;
    statespace.ydot_full_c =  ss(Asys - Bcont*(lqr_c.fo_gain),Bext,Csys,zeros(size(Csys,1),size(Bext,2))) ;
    statespace.ydot_rom_uc = ss(rom.Asys,rom.Bext,rom.Csys,zeros(size(rom.Csys,1),size(rom.Bext,2))) ;
    statespace.ydot_rom_c = ss(rom.Asys - rom.Bcont*(lqr_c.rom_gain),rom.Bext,rom.Csys,zeros(size(rom.Csys,1),size(rom.Bext,2))) ; 
    
    [u_ext,tspan] = gensig('square',5,tf,1/tn) ; 
    u_ext(252:end) = 0; 
     
    [~,~,y_full_uc] = lsim(statespace.ydot_full_uc, u_ext, tspan, y0) ; %Full system - Uncontrolled
    [~,~,y_full_c] = lsim(statespace.ydot_full_c , u_ext , tspan , y0) ; %Full system - Controlled
    [~,~,y_rom_uc] = lsim(statespace.ydot_rom_uc, u_ext , tspan , rom.y0) ; %Reduced system - Uncontrolled  
    [~,~,y_rom_c] = lsim(statespace.ydot_rom_c, u_ext , tspan , rom.y0) ; %Mapped control - Controlled
    
    y_rom_uc = [y_rom_uc(:,1:length(rom.dof))*rom.T_serep',y_rom_uc(:,length(rom.dof)+1:end)*rom.T_serep'] ;
    y_rom_c_mapped = [y_rom_c(:,1:length(rom.dof))*rom.T_serep',y_rom_c(:,length(rom.dof)+1:end)*rom.T_serep'] ;
    rom_Csys = [rom.C , zeros(n_sen,rom.n)] ; 
    
    post_processing.Uncontrolled.y_full_uc(:,:,i) = y_full_uc ; 
    post_processing.Controlled.y_full_c(:,:,i) = y_full_c ;
    post_processing.Uncontrolled.y_rom_uc(:,:,i) = y_rom_uc;
    
    post_processing.Controlled.y_rom_c(:,:,i) = y_rom_c_mapped ; 
    
    post_processing.Controlled.Csys(:,:,i) = Csys ;
    post_processing.Controlled.rom_Csys(:,:,i) = rom_Csys ; 
    
    post_processing.Controlled.fo_gain(:,:,i) = lqr_c.fo_gain ;
    post_processing.Controlled.rom_gain(:,:,i) = lqr_c.rom_gain ; 
    post_processing.Controlled.mfull_gain(:,:,i) = lqr_c.mfull_gain ;
    
    post_processing.Controlled.C(:,:,i) = C ;  
    post_processing.Controlled.rom_C(:,:,i) = rom.C ;
    post_processing.Controlled.y_rom_c_red(:,:,i) = y_rom_c ;
    
end


                                % Post processing
                                 
tn = length(tspan) ;
U_con = zeros(tn, beam.total_dofs,length(distribution.delta));
V_con = zeros(tn, beam.total_dofs,length(distribution.delta));
U_rcon = zeros(tn, beam.total_dofs,length(distribution.delta));
V_rcon = zeros(tn, beam.total_dofs,length(distribution.delta));
U_uncon = zeros(tn, beam.total_dofs,length(distribution.delta));
V_uncon = zeros(tn, beam.total_dofs,length(distribution.delta));
U_runcon = zeros(tn, beam.total_dofs,length(distribution.delta));
V_runcon = zeros(tn, beam.total_dofs,length(distribution.delta));

U_con(:,beam.free_dofs,:) = post_processing.Controlled.y_full_c(:,1:n,:);
V_con(:,beam.free_dofs,:) = post_processing.Controlled.y_full_c(:,(n+1):2*n,:);
U_rcon(:,beam.free_dofs,:) = post_processing.Controlled.y_rom_c(:,1:n,:);
V_rcon(:,beam.free_dofs,:) = post_processing.Controlled.y_rom_c(:,(n+1):2*n,:);
U_uncon(:,beam.free_dofs,:) = post_processing.Uncontrolled.y_full_uc(:,1:n,:);
V_uncon(:,beam.free_dofs,:) = post_processing.Uncontrolled.y_full_uc(:,(n+1):2*n,:);
U_runcon(:,beam.free_dofs,:) = post_processing.Uncontrolled.y_rom_uc(:,1:n,:);
V_runcon(:,beam.free_dofs,:) = post_processing.Uncontrolled.y_rom_uc(:,(n+1):2*n,:);

%% Input Output
Output = zeros(tn,n_act,length(distribution.delta)) ; 
Input = zeros(tn,n_act,length(distribution.delta)) ; 
for j = 1:length(distribution.delta)
    Output(:,:,j) = post_processing.Controlled.y_full_c(:,:,j)*post_processing.Controlled.Csys(:,:,j)' ;
    Input(:,:,j) = post_processing.Controlled.y_full_c(:,:,j)*(-post_processing.Controlled.fo_gain(:,:,j)'); 
end

%% ROM Input Output - Calculate reduced order output then taking it back to Full order

rOutput = zeros(tn,n_act,length(distribution.delta)) ; 
rInput = zeros(tn,n_act,length(distribution.delta)) ; 
for j = 1:length(distribution.delta)
    rOutput(:,:,j) = post_processing.Controlled.y_rom_c_red(:,:,j)*post_processing.Controlled.rom_Csys(:,:,j)' ;
    rInput(:,:,j) = post_processing.Controlled.y_rom_c_red(:,:,j)*(-post_processing.Controlled.rom_gain(:,:,j)') ; 
end



%--------------------------------------------------------------------------
%Saving Solution
%--------------------------------------------------------------------------
solution.Uncon_full = U_uncon;
solution.Uncon_rom = U_runcon;
solution.Vncon_full = V_uncon;
solution.Vncon_rom = V_runcon;

solution.Ucon_full = U_con;
solution.Ucon_rom = U_rcon;
solution.Vcon_full = V_con;
solution.Vcon_rom = V_rcon;

solution.Output_full = Output;
solution.Input_full = Input;
solution.Output_rom = rOutput;
solution.Input_rom = rInput;
solution.rom.dof = rom.dof ;

%--------------------------------------------------------------------------
%Generate the Path for saving RESULTS
%--------------------------------------------------------------------------
folder = sprintf('results/%s/Regulator_problem/%s/forced_vib_control_with_delta_%d_t_%0.1f_t_dof_%d_rom_dof_%d_time_saved_%0.1f',date,rom.config,length(distribution.delta),tf,beam.total_dofs,length(rom.dof),max(per_time_saved));  
mkdir(folder);
mkdir(folder,'PDFs') ;
Pdf_plots = fullfile(folder,'PDFs') ; 

fname = sprintf('%s/dynamics',folder);
filename = 'variables';
filename = sprintf('%s_%s.mat',fname,filename);
save(filename,'tspan','solution','beam');

fname_pdf = sprintf('%s/PDF',Pdf_plots) ; 

lga = [];
lgs = [];
for i=1:n_act
    lga = [lga;sprintf('Actuator_%i',i)];
end

for i=1:n_sen
    lgs = [lgs;sprintf('Sensor_%i',i)];
end

%--------------------------------------------------------------------------
%Calculation for Obtaining PDF 
%--------------------------------------------------------------------------
post_processing.pdf.Ucon_full = zeros(length(distribution.delta),beam.total_dofs) ; 
post_processing.pdf.Uncon_full = zeros(length(distribution.delta),beam.total_dofs) ; 
post_processing.pdf.Input_full = zeros(length(distribution.delta),n_act) ; 
post_processing.pdf.Output_full = zeros(length(distribution.delta),n_sen) ;

post_processing.pdf.Ucon_rom = zeros(length(distribution.delta),beam.total_dofs) ; 
post_processing.pdf.Uncon_rom = zeros(length(distribution.delta),beam.total_dofs) ; 
post_processing.pdf.Input_rom = zeros(length(distribution.delta),n_act) ; 
post_processing.pdf.Output_rom = zeros(length(distribution.delta),n_sen) ;

%--------------------------------------------------------------------------

for i = 1:length(distribution.delta)
    post_processing.pdf.Ucon_full(i,:) = max(abs(U_con(:,:,i))) ; 
    post_processing.pdf.Uncon_full(i,:) = max(abs(U_uncon(:,:,i))) ; 
    post_processing.pdf.Input_full(i,:) = max(abs(Input(:,:,i))) ; 
    post_processing.pdf.Output_full(i,:) = max(abs(Output(:,:,i))) ; 

    post_processing.pdf.Ucon_rom(i,:) = max(abs(U_rcon(:,:,i))) ; 
    post_processing.pdf.Uncon_rom(i,:) = max(abs(U_runcon(:,:,i))) ;
    post_processing.pdf.Input_rom(i,:) = max(abs(rInput(:,:,i))) ; 
    post_processing.pdf.Output_rom(i,:) = max(abs(rOutput(:,:,i))) ;
end

%--------------------------------------------------------------------------

post_processing.pdf.plt_pdf.Ucon_full = zeros(2,100,beam.total_nodes) ;
post_processing.pdf.plt_pdf.Uncon_full = zeros(2,100,beam.total_nodes) ;
post_processing.pdf.plt_pdf.Input_full = zeros(2,100,n_act) ;
post_processing.pdf.plt_pdf.Output_full = zeros(2,100,n_sen) ;

post_processing.pdf.plt_pdf.Ucon_rom = zeros(2,100,beam.total_nodes) ;
post_processing.pdf.plt_pdf.Uncon_rom = zeros(2,100,beam.total_nodes) ;
post_processing.pdf.plt_pdf.Input_rom = zeros(2,100,n_act) ;
post_processing.pdf.plt_pdf.Output_rom = zeros(2,100,n_sen) ;

%--------------------------------------------------------------------------

for i = 1:beam.total_nodes
    [pdf_Uncon_full,x_Uncon_full] = ksdensity(post_processing.pdf.Uncon_full(:,2*i-1));
    [pdf_Uncon_rom,x_Uncon_rom] = ksdensity(post_processing.pdf.Uncon_rom(:,2*i-1));
    [pdf_Ucon_full,x_Ucon_full] = ksdensity(post_processing.pdf.Ucon_full(:,2*i-1));
    [pdf_Ucon_rom,x_Ucon_rom] = ksdensity(post_processing.pdf.Ucon_rom(:,2*i-1));
    post_processing.pdf.plt_pdf.Ucon_full(:,:,i) = [x_Ucon_full;pdf_Ucon_full] ; 
    post_processing.pdf.plt_pdf.Ucon_rom(:,:,i) = [x_Ucon_rom;pdf_Ucon_rom] ;
    post_processing.pdf.plt_pdf.Uncon_full(:,:,i) = [x_Uncon_full;pdf_Uncon_full] ; 
    post_processing.pdf.plt_pdf.Uncon_rom(:,:,i) = [x_Uncon_rom;pdf_Uncon_rom] ;
end

for i = 1:n_act
    [pdf_Input_full,x_Input_full] = ksdensity(post_processing.pdf.Input_full(:,i));
    [pdf_Input_rom,x_Input_rom] = ksdensity(post_processing.pdf.Input_rom(:,i));
    post_processing.pdf.plt_pdf.Input_full(:,:,i) = [x_Input_full;pdf_Input_full] ; 
    post_processing.pdf.plt_pdf.Input_rom(:,:,i) = [x_Input_rom;pdf_Input_rom] ; 
end

for i = 1:n_sen
    [pdf_Output_full,x_Output_full] = ksdensity(post_processing.pdf.Output_full(:,i));
    [pdf_Output_rom,x_Output_rom] = ksdensity(post_processing.pdf.Output_rom(:,i));
    post_processing.pdf.plt_pdf.Output_full(:,:,i) = [x_Output_full;pdf_Output_full] ; 
    post_processing.pdf.plt_pdf.Output_rom(:,:,i) = [x_Output_rom;pdf_Output_rom] ;
end
%--------------------------------------------------------------------------
% Displacement vs PDF of Last node
%--------------------------------------------------------------------------

figure;
plot(post_processing.pdf.plt_pdf.Uncon_full(1,:,end-1),post_processing.pdf.plt_pdf.Uncon_full(2,:,end-1),'r','LineWidth',2) ;
hold on
plot(post_processing.pdf.plt_pdf.Uncon_rom(1,:,end-1),post_processing.pdf.plt_pdf.Uncon_rom(2,:,end-1),'k--','LineWidth',2) ;
grid on
legend('FEA','SEREP')
xlabel(sprintf('Displacement of last node - (Uncontrolled)'))
ylabel(sprintf('PDF of last node'))
filename = sprintf('of last node - Uncontrolled') ; 
filename = sprintf('%s_%s',fname_pdf,filename) ; 
saveas(gcf,filename,'fig') ;
saveas(gcf,filename,'png') ;

figure;
plot(post_processing.pdf.plt_pdf.Ucon_full(1,:,end-1),post_processing.pdf.plt_pdf.Ucon_full(2,:,end-1),'r','LineWidth',2) ;
hold on
plot(post_processing.pdf.plt_pdf.Ucon_rom(1,:,end-1),post_processing.pdf.plt_pdf.Ucon_rom(2,:,end-1),'k--','LineWidth',2) ; 
grid on
legend('FEA','SEREP')
xlabel(sprintf('Displacement of last node - (Controlled)'))
ylabel(sprintf('PDF of last node'))
filename = sprintf('of last node - Controlled') ; 
filename = sprintf('%s_%s',fname_pdf,filename) ; 
saveas(gcf,filename,'fig') ;
saveas(gcf,filename,'png') ;

% --------------------------------------------------------------------------
% Comparision of FEA and SEREP - Uncontrolled
% --------------------------------------------------------------------------
figure;
plot(tspan, U_uncon(:,end-1,1),'r','LineWidth',2);
hold on
plot(tspan, U_runcon(:,end-1,1),'k--','LineWidth',2);

grid on
legend('FEA','SEREP');
xlabel('Time (s)');
ylabel('Deformation (m) - Uncontrolled');
filename = 'End_node_def_Comp_uncont';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

% --------------------------------------------------------------------------
% Comparision of FEA and SEREP - Controlled
% --------------------------------------------------------------------------

figure;
plot(tspan, U_con(:,end-1,1),'r','LineWidth',2);
hold on
plot(tspan, U_rcon(:,end-1,1),'k--','LineWidth',2);

grid on
legend('FEA','SEREP');
xlabel('Time (s)');
ylabel('Deformation (m) - Controlled');
filename = 'End_node_def_Comp_cont';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

err_uncon = relError(post_processing.pdf.plt_pdf.Uncon_full(1,:,end-1),post_processing.pdf.plt_pdf.Uncon_rom(1,:,end-1))  ; 
err_con = relError(post_processing.pdf.plt_pdf.Ucon_full(1,:,end-1),post_processing.pdf.plt_pdf.Ucon_rom(1,:,end-1))  ; 

FullDofs = beam.total_dofs ; 
Romdofs = length(rom.dof) ; 
ErrorPdfUncont = err_uncon(1,2) ; 
ErrorPdfCont = err_con(1,2) ; 
TimeSaved = max(per_time_saved) ;
Text_table = table(FullDofs,Romdofs,ErrorPdfUncont,ErrorPdfCont,TimeSaved) ; 
save_table= fname ;
table_path_format = [save_table 'OutputData.txt'];
writetable(Text_table,table_path_format,'Delimiter','\t');
%--------------------------------------------------------------------------