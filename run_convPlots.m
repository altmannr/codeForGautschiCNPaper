function f1 = run_convPlots
%
% Creates convergence plots (in M and A norms).
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
format long;

T=1;

%%% kinetic bc
%export_txt='PDAE_Euler'; 
%export_txt='PDAE_CNIMEX'; 
export_txt='PDAE_Gautschi'; 

% spatial refinements -- (0:9)
n_vect = [4];
% temporal refinements -- (-0:-1:-8)
tau_vect=2.^(-1:-1:-9);

% legend location
location_text='SouthEast';
% axis limit and reference line height
k=1;
factor=10^(1-k);
ylimits=[10^-6 10^0];

%% reading errors
for i=1:length(n_vect)
    meshno=n_vect(i);

    data=load(['meshes/disk_meshdata',num2str(meshno),'.txt']);
    DOF{i}=['$h=',num2str(data(1,1)),'$'];
    
    for j=1:length(tau_vect)
        tau=tau_vect(j);
            if exist(['errors_',export_txt,'/error_',export_txt,'_n',num2str(meshno),'_tau',num2str(tau),'.txt'], 'file')==2
                err = load(['errors_',export_txt,'/error_',export_txt,'_n',num2str(meshno),'_tau',num2str(tau),'.txt']);
                
                % L^\infty errors on [0,T]
                I = ismember(err(1,:),T);
                Points = (1:length(I));
                T_n = Points(I);
                %%% kinetic bc without splitting 
                error_Om_L2(i,j) = err(1);  % L2-L2 error
                error_Om_H1(i,j) = err(2);  % energy error
            else  
                error_Om_L2(i,j)=NaN;
                error_Om_H1(i,j)=NaN; % semi-norm
            end
    end
end

%% plot
f1=figure('position',[15 55 1200 500]);

% x 0
subplot(1,2,1)
loglog_conv(tau_vect,error_Om_L2)
hold on;
plot(tau_vect,tau_vect.^2/factor,': red','LineWidth',1);
plot(tau_vect,tau_vect./factor,'-. blue','LineWidth',1);
%plot(tau_vect,tau_vect.^0.5/factor,'-- black','LineWidth',1);
hold off; 
h1=title('$\|u-(u^n)^\ell\|_{L^\infty(L^2)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
h1=ylabel(['\verb|',export_txt,'| errors'],'Interpreter','latex');
set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)

% 0 x
subplot(1,2,2)
loglog_conv(tau_vect,error_Om_H1)
loglog_conv(tau_vect,sqrt(error_Om_L2.^2 + error_Om_H1.^2))
hold on;
plot(tau_vect,tau_vect.^2/factor,': red','LineWidth',1);
plot(tau_vect,tau_vect./factor,'-. blue','LineWidth',1);
%plot(tau_vect,tau_vect.^0.5/factor,'-- black','LineWidth',1);
hold off; 
h1=title('$\|u-(u^n)^\ell\|_{L^\infty(H^1)}$');
set(h1,'Interpreter','latex','FontSize',18);
h1=xlabel('step size ($\tau$)');
set(h1,'Interpreter','latex','FontSize',16);
% h1=ylabel('errors','Interpreter','latex');
% set(h1,'Interpreter','latex','FontSize',16);
h2=legend(DOF,'Location',location_text,'FontSize',10);
set(h2,'Interpreter','latex');
xlim([(1/1.2)*min(tau_vect) 1.2*max(tau_vect)])
ylim(ylimits)
end
    
%% 
function loglog_conv(vect,err)
    [n1 n2]=size(err);
    
    symbols='sox.d+^*v><';
    ms=[6 6 6 8 6 6 6 6 6 6 6];
    gr=(linspace(.66,0,n1))';
    colors=[gr gr gr];
    
    for jj=1:n1
        loglog(vect,err(jj,:), ...
               'LineWidth',1,...
               'Marker',symbols(jj),...
               'MarkerSize',ms(jj),...
               'Color', colors(jj,:));%'Color', 'black'); %
        if jj==1
            hold on;
        end
    end
    hold off;
end
