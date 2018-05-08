%% Model Overview
%%%%%% Assumptions/Inputs/Outputs %%%%%%
% 1) Each sensor has a different spring constant
% 2) The "other" FA linkage (for the sake of argument this might be talin)
%    has a single spring constant
% 3) Knowns (inputs)
%           a) Nv: number of vinculin molecules in a single FA
%           b) Nt: number of talin molecules in a single FA
%           c) kv: estimates of the entropic spring constants for our three
%                  sensors (~1.54, 1.09, 0.88 pN/nm for 9, 7, 5 repeats)
%           c) kt: entropic spring constant of talin
%                  choose to vary this up to 2^4 times greater and
%                  less than kv
% 4) Unknowns (outputs)
%           a) Fv: force per vinculin molecule
%           b) Ft: force per talin molecule
%           c) Dv: extension per vinculin molecule
%           d) Dt: extension per talin molecule
%
%%%%%% Brief Model Description %%%%%%
% Model performs a brute force search for model parameters that fulfill the
% strain-control observation

%% Clean Up and Set Up
clc
clear
close all

% Where to save data and filenames
type = 'Hybrid'; % options: 'Hybrid' or 'Parallel'
prot1 = 'Sensor';
prot2 = 'Linker';
if exist(fullfile(pwd,type),'dir')~=7
    mkdir(fullfile(pwd,type));
end
addpath(genpath(pwd));

% make custom red-green colormap
a = linspace(0,0.9,128)';
r = horzcat(a,zeros([128,1]),zeros([128,1]));
g = horzcat(zeros([128,1]),a,zeros([128,1]));
rg = vertcat(flipud(r),g);
rg(:,[1,2,3]) = rg(:,[1,3,2]); % actually switch to red-blue

clear a r g
jet1 = jet;
jet2(:,1) = interp1(1:64,jet1(:,1),0.25:0.25:64);
jet2(:,2) = interp1(1:64,jet1(:,2),0.25:0.25:64);
jet2(:,3) = interp1(1:64,jet1(:,3),0.25:0.25:64);
jet2(end,:) = 0.2;
jet2(1:28,:) = [];

%% Input Variables
% Bulk input (BI) to the entire structure (either a force or an extension
if strcmpi(type,'Hybrid')
    BI = 1.75; % nm
    BIunit = 'nm';
elseif strcmpi(type,'Parallel')
    BI = 300; % pN
    BIunit = 'pN';
end
% Stiffness of Sensor Element
% kv = [1.54, 1.09, 0.88]; % entropic spring constant of VinTS-GGS5,7,9 (pN/nm)
kv = 0.05:0.05:2.3; % same mean so we can compare relative stiffness still

% Constants (unknown)
Nv = 2.^(3.3219:0.5714:7.3219); % number of vinculin molecules in a single FA (sum = 160)
Nt = 2.^(3.3219:0.5714:7.3219); % number of talin molecules in a single FA (sum = 160)
kt = 2.^(-3.77:1.1425:4.25); % entropic spring constant for talin (pN/nm) -- somewhere between 16x less than and 16x greater than that of vinculin
% Nv = 2.^(3.3219:0.005005:7.3219); % number of vinculin molecules in a single FA (sum = 160)
% Nt = 2.^(3.3219:0.005005:7.3219); % number of talin molecules in a single FA (sum = 160)
% kt = 2.^(-3.77:0.01003:4.25); % entropic spring constant for talin (pN/nm) -- somewhere between 16x less than and 16x greater than that of vinculin

% Create all combinations of vectors
Nv_vec = repmat(Nv,[1,length(kt)]);
Nt_vec = repmat(fliplr(Nv),[1,length(kt)]);
kt_vec = repelem(kt,length(Nv));
rel_stiffnes_vec = log2(kt_vec./(mean(kv))); % kt_kv_ratio Metric (y axis)

n = length(Nv_vec);
plottype = 'simple';
mkrsz = 60;
for z = 1:length(BI)
    Dv = zeros(n,3);
    Dt = zeros(n,3);
    Dv_COV = zeros(n,1);
    F = zeros(n,3);
    Fv = zeros(n,3);
    Ft = zeros(n,3);
    Fv_COV = zeros(n,1);
    log2_Fv_Dv_ratio = zeros(n,1); % OK large = more likely extension controlled system (PLOT as dot color)
    log2_Nt_Nv_ratio = zeros(n,1); % OK Relative molecular ratio of talin to vinculin (PLOT as x axis)
    for i = 1:n
        for j = 1:length(kv)
            if strcmpi(type,'Hybrid')
                % Consitutive Equations (solved for Hybrid Model)
                F = (BI(z)/(((kv(j)*Nv_vec(i))^-1+(kt_vec(i)*Nt_vec(i))^-1)));
                Fv(i,j) = F/Nv_vec(i);
                Ft(i,j) = F/Nt_vec(i);
                Dv(i,j) = F/(kv(j)*Nv_vec(i));
                Dt(i,j) = F/(kt_vec(i)*Nt_vec(i));
            elseif strcmpi(type,'Parallel')
                % Consitutive Equations (solved for Parallel Model)
                Dv(i,j) = BI(z)/((kv(j)*Nv_vec(i))+(kt_vec(i)*Nt_vec(i)));
                Dt(i,j) = BI(z)/((kv(j)*Nv_vec(i))+(kt_vec(i)*Nt_vec(i)));
                Fv(i,j) = kv(j)*Dv(i,j);
                Ft(i,j) = kt_vec(i)*Dv(i,j);
            end
        end
        Fv_COV(i,1) = std(Fv(i,:))./mean(Fv(i,:)); % coefficient of variation for forces
        Dv_COV(i,1) = std(Dv(i,:))./mean(Dv(i,:)); % coefficient of variation for forces across vinculin
        log2_Fv_Dv_ratio(i,1) = log2(Fv_COV(i,1)./Dv_COV(i,1)); % Control Metric
        log2_Nt_Nv_ratio(i,1) = log2(Nt_vec(i)./Nv_vec(i)); % Nt_Nv_ratio Metric (x axis)
    end
    
    % Colors for plotting
    c = jet;
    cr = 1:floor(length(c)./length(kt))+1:length(c);
    
    % For each relative stiffness, plot 6 different Nv:Nt ratios as lines
    for k = 1:length(kt)
        figure('Position',[50 500 1000 412]);
        hold on
        r1 = kt_vec == kt(k);
        for m = 1:length(Nv)
            r2 = Nv_vec == Nv(m);
            row = r1 & r2;
            
            subplot(1,2,1)
            hold on
            p1=plot(kv,Fv(row,:),'Color',c(cr(m),:),'LineWidth',2);
            axis([0 2.5 0 5]);
            xlabel('Spring Constant,  {\bf{\itk_{Sj}}} (pN/nm)')
            ylabel('Force,  {\bf{\itF_{Sj}}} (pN)');
            PrettyPlot
        
            subplot(1,2,2)
            hold on
            p2=plot(kv,Dv(row,:),'Color',c(cr(m),:),'LineWidth',2);
            axis([0 2.5 0 5]);
            xlabel('Spring Constant,  {\bf{\itk_{Sj}}} (pN/nm)')
            ylabel('Extension,  {\bf{\it\delta_{Sj}}} (nm)');
            PrettyPlot
        
            a = round(log2(kt(k)/mean(kv)),1);
            b = round(log2(Nt_vec(m)/Nv_vec(m)),1);
            t1 = a==4.0 & b==-0.6; %%%%%%%%%%%%%%%%%%%%%%%%%%%
            t2 = a==2.9 & b==0.6; %%%%%%%%%%%%%%%%%%%%%%%%%%%
            t3 = a==1.7 & b==1.7; %%%%%%%%%%%%%%%%%%%%%%%%%%%
            t4 = a==0.6 & b==2.9; %%%%%%%%%%%%%%%%%%%%%%%%%%%
            t5 = a==-0.6 & b==4.0; %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum([t1,t2,t3,t4,t5])==1
                p1.Color(4) = 1.0;
                p2.Color(4) = 1.0;
            else
                p1.Color(4) = 0.25;
                p2.Color(4) = 0.25;
            end
            formatSpec1 = 'Relative Stiffness = %0.1f';
            title(sprintf(formatSpec1,log2(kt(k)/mean(kv))));
            formatSpec2 = '%0.1f';
            l{m} = sprintf(formatSpec2,log2(Nt_vec(m)/Nv_vec(m)));
        end
        hleg = legend(l,'Location','NE');
        htitle = get(hleg,'Title');
        set(htitle,'String','Relative Abundance');
        set(hleg,'FontSize',10);
        
        % Saving Plot 0a
        outname0a = strcat(plottype,'_', type, '_StiffnessVersusForce_RelativeStiffness', num2str(k),'_bulk_input', strrep(num2str(BI(z)),'.','p') ,BIunit);
        saveas(gcf,fullfile(pwd,type,outname0a),'png')
    end
    
    % For each relative abundance, plot 6 different kv:kt ratios as lines
    for k = 1:length(Nv)
        figure('Position',[650 500 1000 412]);
        hold on
        r1 = Nv_vec == Nv(k);
        for m = 1:length(kt)
            r2 = kt_vec == kt(m);
            row = r1 & r2;
            
            subplot(1,2,1)
            hold on
            p3 = plot(kv,Fv(row,:),'Color',c(cr(m),:),'LineWidth',2);
            axis([0 2.5 0 5]);
            xlabel('Spring Constant,  {\bf{\itk_{Sj}}} (pN/nm)')
            ylabel('Force,  {\bf{\itF_{Sj}}} (pN)');
            PrettyPlot
            
            subplot(1,2,2)
            hold on
            p4 = plot(kv,Dv(row,:),'Color',c(cr(m),:),'LineWidth',2);
            axis([0 2.5 0 5]);
            xlabel('Spring Constant,  {\bf{\itk_{Sj}}} (pN/nm)')
            ylabel('Extension,  {\bf{\it\delta_{Sj}}} (nm)');
            PrettyPlot
            
            a = round(log2(kt(m)/mean(kv)),1);
            b = round(log2(Nt_vec(k)/Nv_vec(k)),1);
            t1 = a==4.0 & b==-0.6;
            t2 = a==2.9 & b==0.6;
            t3 = a==1.7 & b==1.7;
            t4 = a==0.6 & b==2.9;
            t5 = a==-0.6 & b==4.0;
            if sum([t1,t2,t3,t4,t5])==1
                p3.Color(4) = 1.0;
                p4.Color(4) = 1.0;
            else
                p3.Color(4) = 0.25;
                p4.Color(4) = 0.25;
            end
            formatSpec1 = 'Relative Abundance = %0.1f';
            title(sprintf(formatSpec1,log2(Nt_vec(k)/Nv_vec(k))));
            formatSpec2 = '%0.1f';
            l{m} = sprintf(formatSpec2,log2(kt(m)/mean(kv)));
        end
        hleg = legend(l,'Location','NW');
        htitle = get(hleg,'Title');
        set(htitle,'String','Relative Stiffness')
        
        % Saving Plot 0b
        outname0b = strcat(plottype,'_', type, '_StiffnessVersusForce_RelativeAbundance', num2str(k),'_delta', strrep(num2str(BI(z)),'.','p'),'nm');
        saveas(gcf,fullfile(pwd,type,outname0b),'png')
    end
    
    
    %% Plot log2Fv_delta_ratio
    figure('Position',[500 500 490 412])
    scatter(log2_Nt_Nv_ratio,rel_stiffnes_vec,mkrsz,log2_Fv_Dv_ratio,'Filled');
    PrettyPlot
    xlabel({['N_{' prot1(1) '}> N_{' prot2(1) '}                        N_{' prot1(1) '}= N_{' prot2(1) '}                        N_{' prot1(1) '}< N_{' prot2(1) '}'];['log_{2}(N_{' prot2 '}/N_{' prot1 '})']}); % vinculin:talin molecular ratio
    ylabel({['log_{2}(k_{' prot2 '}/k_{' prot1 '})'];['k_{' prot1(1) '}> k_{' prot2(1) '}                    k_{' prot1(1) '}= k_{' prot2(1) '}                    k_{' prot1(1) '}< k_{' prot2(1) '}']}); % vinculin:talin stiffness ratio
    axis([-4.5 4.5 -4.6 4.4]);
    colormap(rg);
    c = colorbar;
    c.Label.String = 'Force-Controlled                    Extension-Controlled';
    set(gca,'FontSize',11);
    c.Label.FontSize = 11;
    c.FontSize = 11;
    caxis([-5,5]);
    title(['Extension = ' num2str(BI(z)) 'nm']);
    
    % Saving Plot 1
    outname1 = strcat(plottype,'_', type, '_ControlMetric_', strrep(num2str(BI(z)),'.','p') ,'nm');
    saveas(gcf,fullfile(pwd,type,outname1),'png')
    
    %% Clear a bunch of variables for the next force
    close all
    clear Fv Ft Dv Dt Dv_COV Fv_COV i j log2_Fv_Dv_ratio log2_Nt_Nv_ratio outname0 outname1 p1 p2 p3 p4
end
rmpath(genpath(pwd))