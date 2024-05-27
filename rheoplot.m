function [] = rheoplot(timetype,rheodata,vemodel,flowtype,plottype)

    % extract stresses from rheodata for convenience
    sxx = rheodata.stress(1,:); sxy = rheodata.stress(2,:); sxz = rheodata.stress(3,:);
    syy = rheodata.stress(4,:); syz = rheodata.stress(5,:); szz = rheodata.stress(6,:);
    time = rheodata.time;

    if flowtype == 1 % simple shear

        if strcmp(timetype,'startup')

            if strcmp(plottype,'visc')
    
                figure; 
                plot(time,sxy/rheodata.rate_for_startup,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient shear viscosity $\eta(t)$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \eta $','interpreter', 'LaTeX','FontSize',28); % y-axis label
        
                figure; 
                plot(time,(sxx-syy)/rheodata.rate_for_startup^2,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient first normal stress coefficient $\Psi_{1}(t)$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \Psi_{1}$','interpreter', 'LaTeX','FontSize',28); % y-axis label  
        
                figure; 
                plot(time,(sxx-syy)./sxy,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient stress ratio $S(t)=N_1(t)/\tau_{xy}(t)$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$S$','interpreter', 'LaTeX','FontSize',28); % y-axis label
        
                figure; 
                plot(time,(syy-szz)/rheodata.rate_for_startup^2,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient second normal stress coefficient $\Psi_{2}(t)$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \Psi_{2}$','interpreter', 'LaTeX','FontSize',28); % y-axis label

            elseif strcmp(plottype,'stress')

                figure; 
                plot(time,sxy,'LineWidth',2); hold on;                
                plot(time,sxx-syy,'LineWidth',2);
                plot(time,syy-szz,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient stress','Interpreter','LaTeX','FontSize',24)                
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \tau $','interpreter', 'LaTeX','FontSize',28); % y-axis label
                legend('$\tau_{xy}$','$N_1 = \tau_{xx}-\tau_{yy}$','$N_2 = \tau_{yy}-\tau_{zz}$','Interpreter','LaTeX')

            end
    
        elseif strcmp(timetype,'steady')
    
            if strcmp(plottype,'visc')

                figure; 
                plot(rheodata.rates,sxy./rheodata.rates,'LineWidth',2)
                set(gca,'FontSize',16); set(gca,'xscale','log'); set(gca,'yscale','log')
                title('Steady shear viscosity $\eta(\dot{\gamma})$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \eta $','interpreter', 'LaTeX','FontSize',28); % y-axis label
        
                figure; 
                plot(rheodata.rates,(sxx-syy)./rheodata.rates.^2,'LineWidth',2)
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady first normal stress coefficient $\Psi_{1}(\dot{\gamma})$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \Psi_{1}$','interpreter', 'LaTeX','FontSize',28); % y-axis label  
        
                figure; 
                plot(rheodata.rates,(sxx-syy)./sxy,'LineWidth',2)
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady stress ratio $S(\dot{\gamma})=N_1(\dot{\gamma})/\tau_{xy}(\dot{\gamma})$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$S$','interpreter', 'LaTeX','FontSize',28); % y-axis label
        
                figure; 
                plot(rheodata.rates,(syy-szz)./rheodata.rates.^2,'LineWidth',2)
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady second normal stress coefficient $\Psi_{2}(\dot{\gamma})$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \Psi_{2}$','interpreter', 'LaTeX','FontSize',28); % y-axis label

            elseif strcmp(plottype,'stress')

                figure; 
                plot(rheodata.rates,sxy,'LineWidth',2); hold on;                
                plot(rheodata.rates,sxx-syy,'LineWidth',2);
                if any(syy-szz<0)
                    plot(rheodata.rates,-(syy-szz),'LineWidth',2);
                else
                    plot(rheodata.rates,syy-szz,'LineWidth',2);                    
                end
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady-state stress','Interpreter','LaTeX','FontSize',24)                
                xlabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$\tau$','interpreter', 'LaTeX','FontSize',28); % y-axis label
                if any(syy-szz<0)
                    legend('$\tau_{xy}$','$N_1 = \tau_{xx}-\tau_{yy}$', ...
                        '$-N_2 = \tau_{zz}-\tau_{yy}$','Interpreter','LaTeX')
                else
                    legend('$\tau_{xy}$','$N_1 = \tau_{xx}-\tau_{yy}$', ...
                        '$N_2 = \tau_{yy}-\tau_{zz}$','Interpreter','LaTeX')
                end

            end
    
        elseif strcmp(timetype,'startup_stress')
    
            if strcmp(plottype,'strain')
            
                figure; 
                plot(rheodata.time,rheodata.strain,'LineWidth',2)
                set(gca,'FontSize',16);
                set(gca,'xscale','linear')
                set(gca,'yscale','linear')
                title('Transient shear strain for imposed shear stress','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$\gamma$','interpreter', 'LaTeX','FontSize',28); % y-axis label

            elseif strcmp(plottype,'rate')

                figure; 
                plot(rheodata.time,rheodata.rates,'LineWidth',2)
                set(gca,'FontSize',16);
                set(gca,'xscale','linear')
                set(gca,'yscale','linear')
                title('Transient strain-rate for imposed shear stress','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$\dot{\gamma}$','interpreter', 'LaTeX','FontSize',28); % y-axis label                

            end

        end

    elseif flowtype == 2 || flowtype == 3  % uni/bi-extensional flow

        if strcmp(timetype,'startup')

            if strcmp(plottype,'visc')

                figure; 
                plot(time,(sxx-syy)/rheodata.rate_for_startup,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient extensional viscosity $\eta_{\rm E}(t)$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \eta_{\rm E} $','interpreter', 'LaTeX','FontSize',28); % y-axis label

            elseif strcmp(plottype,'stress')

                figure; 
                plot(time,sxx-syy,'LineWidth',2)
                set(gca,'xscale','linear'); set(gca,'yscale','linear'); set(gca,'FontSize',16);
                title('Transient normal stress $N_1(t)=\tau_{xx}-\tau_{yy}$','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ N_1 $','interpreter', 'LaTeX','FontSize',28); % y-axis label

            end

        elseif strcmp(timetype,'steady')

            if strcmp(plottype,'visc')

                figure; 
                plot(rheodata.rates,(sxx-syy)./rheodata.rates,'LineWidth',2)
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady extensional viscosity $\eta_{\rm E}(\dot{\epsilon})$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\epsilon}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$ \eta_{\rm E} $','interpreter', 'LaTeX','FontSize',28); % y-axis label

            elseif strcmp(plottype,'stress')

                figure; 
                plot(rheodata.rates,(sxx-syy)./rheodata.rates,'LineWidth',2)
                set(gca,'xscale','log'); set(gca,'yscale','log'); set(gca,'FontSize',16);
                title('Steady normal stress $N_1(\dot{\epsilon})=\tau_{xx}-\tau_{yy}$','Interpreter','LaTeX','FontSize',24)
                xlabel('$\dot{\epsilon}$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$N_1$','interpreter', 'LaTeX','FontSize',28); % y-axis label

            end

         elseif strcmp(timetype,'startup_stress')
    
            if strcmp(plottype,'strain')
             
                figure; 
                plot(rheodata.time,rheodata.strain,'LineWidth',2)
                set(gca,'FontSize',16);
                set(gca,'xscale','linear')
                set(gca,'yscale','linear')
                title('Transient elongational strain for imposed normal stress','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$\epsilon$','interpreter', 'LaTeX','FontSize',28); % y-axis label
            
            elseif strcmp(plottype,'rate')
             
                figure; 
                plot(rheodata.time,rheodata.rates,'LineWidth',2)
                set(gca,'FontSize',16);
                set(gca,'xscale','linear')
                set(gca,'yscale','linear')
                title('Transient elongational rate for imposed normal stress','Interpreter','LaTeX','FontSize',24)
                xlabel('$t$','interpreter', 'LaTeX','FontSize',28); % x-axis label
                ylabel('$\dot{\epsilon}$','interpreter', 'LaTeX','FontSize',28); % y-axis label

            end

         end
    end
