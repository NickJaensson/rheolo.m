function [] = rheoplot(type,rheodata,vemodel)

    if strcmp(type,'startup')

        figure; 
        plot(rheodata.time,rheodata.stress(2,:)/rheodata.rate_for_startup,'LineWidth',2)
        set(gca,'FontSize',16);
        set(gca,'xscale','linear')
        set(gca,'yscale','linear')
        title('Transient shear viscosity $\eta(t)$','Interpreter','LaTeX','FontSize',24)
        x = xlabel('$t$','FontSize',28); % x-axis label
        y = ylabel('$ \eta $','FontSize',28); % y-axis label
        set(x, 'interpreter', 'LaTeX')
        set(y, 'interpreter', 'LaTeX')

    elseif strcmp(type,'steady')

        figure; 
        plot(rheodata.rates,rheodata.stress(2,:)./rheodata.rates,'LineWidth',2)
        set(gca,'FontSize',16);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        title('Steady shear viscosity $\eta(\dot{\gamma})$','Interpreter','LaTeX','FontSize',24)
        x = xlabel('$\dot{\gamma}$','FontSize',28); % x-axis label
        y = ylabel('$ \eta $','FontSize',28); % y-axis label
        set(x, 'interpreter', 'LaTeX')
        set(y, 'interpreter', 'LaTeX')

    elseif strcmp(type,'startup_stress')

        figure; 
        plot(rheodata.time,rheodata.strain,'LineWidth',2)
        set(gca,'FontSize',16);
        set(gca,'xscale','linear')
        set(gca,'yscale','linear')
        title('Transient strain for imposed stress','Interpreter','LaTeX','FontSize',24)
        x = xlabel('$\dot{\gamma}$','FontSize',28); % x-axis label
        y = ylabel('$ \eta $','FontSize',28); % y-axis label
        set(x, 'interpreter', 'LaTeX')
        set(y, 'interpreter', 'LaTeX')
    

    end

end
