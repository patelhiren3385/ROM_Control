                                
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
folder = sprintf('results/%s/Regulator_problem/%s/forced_vib_control_with_delta_%d_t_%0.1f_t_dof_%d_rom_dof_%d',date,rom.config,length(distribution.delta),tf,beam.total_dofs,length(rom.dof));  
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
% 
% for i = 1:length(distribution.delta)
%     post_processing.pdf.Ucon_full(i,:) = (-min(U_con(:,:,i)) < max(U_con(:,:,i))).*(max(U_con(:,:,i)) - min(U_con(:,:,i))) + min(U_con(:,:,i)) ; 
%     post_processing.pdf.Uncon_full(i,:) = (-min(U_uncon(:,:,i)) < max(U_uncon(:,:,i))).*(max(U_uncon(:,:,i)) - min(U_uncon(:,:,i))) + min(U_uncon(:,:,i)) ; 
%     post_processing.pdf.Input_full(i,:) = (-min(Input(:,:,i)) < max(Input(:,:,i))).*(max(Input(:,:,i)) - min(Input(:,:,i))) + min(Input(:,:,i)) ; 
%     post_processing.pdf.Output_full(i,:) = (-min(Output(:,:,i)) < max(Output(:,:,i))).*(max(Output(:,:,i)) - min(Output(:,:,i))) + min(Output(:,:,i)) ; 
%     
%     post_processing.pdf.Ucon_rom(i,:) = (-min(U_rcon(:,:,i)) < max(U_rcon(:,:,i))).*(max(U_rcon(:,:,i)) - min(U_rcon(:,:,i))) + min(U_rcon(:,:,i)) ; 
%     post_processing.pdf.Uncon_rom(i,:) = (-min(U_runcon(:,:,i)) < max(U_runcon(:,:,i))).*(max(U_runcon(:,:,i)) - min(U_runcon(:,:,i))) + min(U_runcon(:,:,i)) ; 
%     post_processing.pdf.Input_rom(i,:) = (-min(rInput(:,:,i)) < max(rInput(:,:,i))).*(max(rInput(:,:,i)) - min(rInput(:,:,i))) + min(rInput(:,:,i)) ; 
%     post_processing.pdf.Output_rom(i,:) = (-min(rOutput(:,:,i)) < max(rOutput(:,:,i))).*(max(rOutput(:,:,i)) - min(rOutput(:,:,i))) + min(rOutput(:,:,i)) ; 
%     
% end


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
%Visulization of system
%--------------------------------------------------------------------------
figure ; 
visulization.beam_total_dofs = 1:elements.beam ;
visulization.beam_active_dofs = ceil(rom.dof/2) ;
visulization.ele_length = (beam.length/elements.beam)*1000 ; 
visulization.length_of_sens_act = visulization.ele_length*size(element_act,2) ; 


plot(0:1,zeros(2,1),'k','linewidth',80) ;
hold on
plot(visulization.beam_total_dofs,zeros(length(visulization.beam_total_dofs),1),'y','linewidth',15) ;
hold on 
for i = 1:n_act
 plot(element_act(i,:),zeros(size(element_act,2),1),'r','linewidth',5)
 hold on
end
hold on
plot(visulization.beam_active_dofs,zeros(length(visulization.beam_active_dofs),1),'k*','MarkerSize',5) ; 
hold on
plot(visulization.beam_active_dofs,zeros(length(visulization.beam_active_dofs),1),'kd','MarkerSize',10) ; 
hold on
dim_ = [.2 .5 .3 .3];
str = sprintf(' Beam - %d mm \n No of elements - %d \n Length of element - %d mm \n Length of sensors & actuators - %d mm \n Marker indicates Active DOFs',1000.*beam.length,elements.beam,visulization.ele_length,visulization.length_of_sens_act);
annotation('textbox',dim_,'String',str,'FitBoxToText','on')
xlabel('Number of Elements')
filename = 'Visulization of System';
filename = sprintf('%s_%s',fname,filename);
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

% %--------------------------------------------------------------------------
% %PDF of all Nodes 
% %--------------------------------------------------------------------------
% % for i = 1:beam.total_nodes
% %     figure;
% %     plot(post_processing.pdf.plt_pdf.Uncon_full(1,:,i),post_processing.pdf.plt_pdf.Uncon_full(2,:,i),'LineWidth',2) ;
% %     xlabel(sprintf('Displacement of %ith node - FEA (Uncontrolled)',i))
% %     ylabel(sprintf('PDF of %ith node',i))
% %     filename = sprintf('of Node number %i - FEA (Uncontrolled)',i) ; 
% %     filename = sprintf('%s_%s',fname_pdf,filename) ; 
% %     saveas(gcf,filename,'fig') ;
% %     saveas(gcf,filename,'png') ;   
% % end
% % 
% % for i = 1:beam.total_nodes
% %     figure;
% %     plot(post_processing.pdf.plt_pdf.Uncon_rom(1,:,i),post_processing.pdf.plt_pdf.Uncon_rom(2,:,i),'LineWidth',2) ; 
% %     xlabel(sprintf('Displacement of %ith node - SEREP (Uncontrolled)',i))
% %     ylabel(sprintf('PDF of %ith node',i))
% %     filename = sprintf('of Node number %i - SEREP (Uncontrolled)',i) ; 
% %     filename = sprintf('%s_%s',fname_pdf,filename) ; 
% %     saveas(gcf,filename,'fig') ;
% %     saveas(gcf,filename,'png') ;   
% % end
% % 
% % 
% % for i = 1:beam.total_nodes
% %     figure;
% %     plot(post_processing.pdf.plt_pdf.Ucon_full(1,:,i),post_processing.pdf.plt_pdf.Ucon_full(2,:,i),'LineWidth',2) ;
% %     xlabel(sprintf('Displacement of %ith node - FEA (Controlled)',i))
% %     ylabel(sprintf('PDF of %ith node',i))
% %     filename = sprintf('of Node number %i - FEA (Controlled)',i) ; 
% %     filename = sprintf('%s_%s',fname_pdf,filename) ; 
% %     saveas(gcf,filename,'fig') ;
% %     saveas(gcf,filename,'png') ;   
% % end
% % 
% % for i = 1:beam.total_nodes
% %     figure;
% %     plot(post_processing.pdf.plt_pdf.Ucon_rom(1,:,i),post_processing.pdf.plt_pdf.Ucon_rom(2,:,i),'LineWidth',2) ; 
% %     xlabel(sprintf('Displacement of %ith node - SEREP (Controlled)',i))
% %     ylabel(sprintf('PDF of %ith node',i))
% %     filename = sprintf('of Node number %i - SEREP (Controlled)',i) ; 
% %     filename = sprintf('%s_%s',fname_pdf,filename) ; 
% %     saveas(gcf,filename,'fig') ;
% %     saveas(gcf,filename,'png') ;   
% % end
% % 

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

%--------------------------------------------------------------------------
% Frequencies comparison of Full and SEREP
%--------------------------------------------------------------------------
fr = 1:1:200;
omega = fr*2*pi;

KK = rom.T_serep*rom.K*rom.T_serep' ; 
MM = rom.T_serep*rom.M*rom.T_serep' ;

F_full = zeros(size(K,1),size(K,1),length(fr)) ;
F_se = zeros(size(rom.K,1),size(rom.K,1),length(fr)) ;
for i = 1:length(omega)
    F_full(:,:,i) = (K - omega(i).^2.*M)\eye(size(K,1)) ;
    F_se(:,:,i) = (rom.K - omega(i).^2.*rom.M)\eye(size(rom.K,1)) ;  
end
F_1(1,:) = abs(F_full(1,1,:));
F_1_se(1,:) = abs(F_se(1,1,:)); 

figure;
semilogy(fr,F_1,'r','LineWidth',1);
hold on
semilogy(fr,F_1_se,'--k','LineWidth',1);
xlabel('Frequency, rad/s'); 
ylabel('FRF');
legend('FEA','SEREP')
filename = sprintf('Frequency comparison') ; 
filename = sprintf('%s_%s',fname,filename) ; 
saveas(gcf,filename,'fig') ;
saveas(gcf,filename,'png') ;  
saveas(gcf,filename,'eps') ;  

%--------------------------------------------------------------------------
% MAC
%--------------------------------------------------------------------------

for I=1:24
    for J=1:24
        mac(I,J)=Mac(Mode_ROMt(:,I),eigenvectors(:,J));
    end
end
figure;
bar3(mac);
xlabel(sprintf('Modes - FEA'))
ylabel(sprintf('Modes - SEREP'))
zlabel(sprintf('MAC'))
filename = sprintf('MAC') ; 
filename = sprintf('%s_%s',fname,filename) ; 
saveas(gcf,filename,'fig') ;
saveas(gcf,filename,'png') ; 
saveas(gcf,filename,'eps') ;  

%--------------------------------------------------------------------------
% PDF of Sensor Output and Actuator Input - FEA & SEREP                   
%--------------------------------------------------------------------------

% figure;
% for i = 1:n_act
%     plot(post_processing.pdf.plt_pdf.Input_full(1,:,i),post_processing.pdf.plt_pdf.Input_full(2,:,i),'LineWidth',2) ;
%     hold on
%     grid on
%     legend(lga)
%     xlabel(sprintf('Voltage input for actuators - FEA '))
%     ylabel(sprintf('PDF of actuators input'))
%     filename = sprintf('of actuator - FEA') ; 
%     filename = sprintf('%s_%s',fname_pdf,filename) ; 
%     saveas(gcf,filename,'fig') ;
%     saveas(gcf,filename,'png') ;   
% end
% 
% figure;
% for i = 1:n_act
%     plot(post_processing.pdf.plt_pdf.Input_rom(1,:,i),post_processing.pdf.plt_pdf.Input_rom(2,:,i),'LineWidth',2) ;
%     hold on
%     grid on
%     legend(lga)
%     xlabel(sprintf('Voltage input for actuators - SEREP '))
%     ylabel(sprintf('PDF of actuators input')) 
%     filename = sprintf('of actuator - SEREP') ; 
%     filename = sprintf('%s_%s',fname_pdf,filename) ; 
%     saveas(gcf,filename,'fig') ;
%     saveas(gcf,filename,'png') ;   
% end
% 
% figure;
% for i = 1:n_sen
%     plot(post_processing.pdf.plt_pdf.Output_full(1,:,i),post_processing.pdf.plt_pdf.Output_full(2,:,i),'LineWidth',2) ; 
%     hold on
%     grid on
%     legend(lgs)
%     xlabel(sprintf('Voltage output of sensors - FEA'))
%     ylabel(sprintf('PDF of sensors output'))
%     filename = sprintf('of sensor - FEA') ; 
%     filename = sprintf('%s_%s',fname_pdf,filename) ; 
%     saveas(gcf,filename,'fig') ;
%     saveas(gcf,filename,'png') ;   
% end
% 
% figure;
% for i = 1:n_sen
%     plot(post_processing.pdf.plt_pdf.Output_rom(1,:,i),post_processing.pdf.plt_pdf.Output_rom(2,:,i),'LineWidth',2) ;
%     hold on
%     grid on
%     legend(lgs)
%     xlabel(sprintf('Voltage output of sensors - SEREP '))
%     ylabel(sprintf('PDF of sensors output'))
%     filename = sprintf('of sensor - SEREP') ; 
%     filename = sprintf('%s_%s',fname_pdf,filename) ; 
%     saveas(gcf,filename,'fig') ;
%     saveas(gcf,filename,'png') ;
% end
% 
% % %--------------------------------------------------------------------------
% % %Uncontrolled & Controlled Surface Plot for FEA and SEREP
% % %--------------------------------------------------------------------------
% % % figure;
% % % surf(real(U_uncon(:,wdof,1)));
% % % xlabel('Length (m) - FEA)');
% % % ylabel('Time (s)');
% % % zlabel('Deformation (m) - Uncontrolled');
% % % filename = 'space_time_UC_FEA';
% % % filename = sprintf('%s_%s',fname,filename);
% % % saveas(gcf,filename,'fig');
% % % saveas(gcf,filename,'png');
% % % 
% % % figure;
% % % surf(real(U_runcon(:,wdof,1)));
% % % xlabel('Length (m) - SEREP');
% % % ylabel('Time (s)');
% % % zlabel('Deformation (m) - Uncontrolled');
% % % filename = 'space_time_UC_SEREP';
% % % filename = sprintf('%s_%s',fname,filename);
% % % saveas(gcf,filename,'fig');
% % % saveas(gcf,filename,'png');
% % % 
% % % figure;
% % % surf(real(U_con(:,wdof,1)));
% % % xlabel('Length (m) - FEA');
% % % ylabel('Time (s)');
% % % zlabel('Deformation (m) - Controlled');
% % % filename = 'space_time_C_FEA';
% % % filename = sprintf('%s_%s',fname,filename);
% % % saveas(gcf,filename,'fig');
% % % saveas(gcf,filename,'png');
% % % 
% % % figure;
% % % surf(real(U_rcon(:,wdof,1)));
% % % xlabel('Length (m) - SEREP');
% % % ylabel('Time (s)');
% % % zlabel('Deformation (m) - Controlled');
% % % filename = 'space_time_C_SEREP';
% % % filename = sprintf('%s_%s',fname,filename);
% % % saveas(gcf,filename,'fig');
% % % saveas(gcf,filename,'png');

%--------------------------------------------------------------------------
%Comparision of FEA and SEREP - Uncontrolled
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
%Comparision of FEA and SEREP - Controlled
%--------------------------------------------------------------------------

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

%--------------------------------------------------------------------------
%Time series plot of Sensor Output and Actuator Input
%--------------------------------------------------------------------------

% figure;
% plot(tspan, Output(:,:,1),'LineWidth',2);
% grid on
% legend(lgs);
% xlabel('Time (s) - FEA');
% ylabel('Voltage (V)');
% filename = 'sensor_FEA';
% filename = sprintf('%s_%s',fname,filename);
% saveas(gcf,filename,'fig');
% saveas(gcf,filename,'png');
% 
% figure;
% plot(tspan, Input(:,:,1),'LineWidth',2);
% grid on
% legend(lga);
% xlabel('Time (s) - FEA');
% ylabel('Voltage (V)');
% filename = 'actuator_FEA';
% filename = sprintf('%s_%s',fname,filename);
% saveas(gcf,filename,'fig');
% saveas(gcf,filename,'png');
% 
% 
% figure;
% plot(tspan, rOutput(:,:,1),'LineWidth',2);
% legend(lgs);
% grid on
% xlabel('Time (s) - SEREP');
% ylabel('Voltage (V)');
% filename = 'sensor - SEREP';
% filename = sprintf('%s_%s',fname,filename);
% saveas(gcf,filename,'fig');
% saveas(gcf,filename,'png');
% 
% figure;
% plot(tspan, rInput(:,:,1),'LineWidth',2);
% grid on
% legend(lga);
% xlabel('Time (s) - SEREP');
% ylabel('Voltage (V)');
% filename = 'actuator - SEREP';
% filename = sprintf('%s_%s',fname,filename);
% saveas(gcf,filename,'fig');
% saveas(gcf,filename,'png');

%--------------------------------------------------------------------------