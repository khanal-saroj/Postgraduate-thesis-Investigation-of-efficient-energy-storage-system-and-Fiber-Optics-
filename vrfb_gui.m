function vrfb_gui()

app = struct();
app.params = init_params();
app.stopRequested = false;
app.isRunning = false;
app.sim = [];


app.fig = uifigure('Name','VRFB Simulator','Position',[100 100 1250 750]);

% panel controls
ctrlPanel = uipanel(app.fig,'Title','Controls','Position',[10 10 300 730]); 

uilabel(ctrlPanel,'Position',[10 670 120 22],'Text','Current (A):');
app.sldI = uislider(ctrlPanel,'Position',[10 660 280 3],'Limits',[-1000 1000],'Value',app.params.I_cmd);
app.edtI = uieditfield(ctrlPanel,'numeric','Position',[10 620 120 24],'Value',app.params.I_cmd);

uilabel(ctrlPanel,'Position',[10 580 120 22],'Text','Flow (L/min):');
app.sldQ = uislider(ctrlPanel,'Position',[10 570 280 3],'Limits',[10 500],'Value',app.params.Q_flow_Lpm);
app.edtQ = uieditfield(ctrlPanel,'numeric','Position',[10 530 120 24],'Value',app.params.Q_flow_Lpm);

uilabel(ctrlPanel,'Position',[10 490 120 22],'Text','Sim time (min):');
app.edtT = uieditfield(ctrlPanel,'numeric','Position',[10 460 120 24],'Value',app.params.t_end_min);
uilabel(ctrlPanel,'Position',[150 490 120 22],'Text','dt (s):');
app.edtDt = uieditfield(ctrlPanel,'numeric','Position',[150 460 120 24],'Value',app.params.dt);

app.btnStart = uibutton(ctrlPanel,'push','Text','Start','Position',[10 410 80 30],'ButtonPushedFcn',@(btn,event) onStart());
app.btnStop  = uibutton(ctrlPanel,'push','Text','Stop','Position',[110 410 80 30],'ButtonPushedFcn',@(btn,event) onStop());
app.btnReset = uibutton(ctrlPanel,'push','Text','Reset','Position',[210 410 80 30],'ButtonPushedFcn',@(btn,event) onReset());

app.lblStatus = uilabel(ctrlPanel,'Position',[10 370 280 30],'Text','Status: Idle');

uilabel(ctrlPanel,'Position',[10 330 200 22],'Text','Number of cells:');
app.lblNcell = uilabel(ctrlPanel,'Position',[10 310 200 22],'Text',num2str(app.params.Ncell));
uilabel(ctrlPanel,'Position',[10 280 200 22],'Text','Cell area (cm^2):');
app.lblAcell = uilabel(ctrlPanel,'Position',[10 260 200 22],'Text',num2str(app.params.A_cell*1e4));


axW = 880; axH = 710;
app.ax1 = uiaxes(app.fig,'Position',[320 530 axW/2 axH/3]); title(app.ax1,'Stack Voltage (V)'); xlabel(app.ax1,'Time (min)');
app.ax2 = uiaxes(app.fig,'Position',[320 270 axW/2 axH/3]); title(app.ax2,'Species (Tank)');   xlabel(app.ax2,'Time (min)');
app.ax3 = uiaxes(app.fig,'Position',[320+axW/2 530 axW/2 axH/3]); title(app.ax3,'Stack Current (A)'); xlabel(app.ax3,'Time (min)');
app.ax4 = uiaxes(app.fig,'Position',[320+axW/2 270 axW/2 axH/3]); title(app.ax4,'SOC (neg tank)');    xlabel(app.ax4,'Time (min)');
app.ax5 = uiaxes(app.fig,'Position',[320 10 axW/2 axH/3]);        title(app.ax5,'Stack Temperature (FBG)'); xlabel(app.ax5,'Time (min)');
app.ax6 = uiaxes(app.fig,'Position',[320+axW/2 10 axW/2 axH/3]);  title(app.ax6,'FBG Wavelength (nm)');      xlabel(app.ax6,'Time (min)');

app.sldI.ValueChangedFcn = @(s,e) syncIfromSlider();
app.edtI.ValueChangedFcn = @(s,e) syncIfromEdit();
app.sldQ.ValueChangedFcn = @(s,e) syncQfromSlider();
app.edtQ.ValueChangedFcn = @(s,e) syncQfromEdit();

drawnow;


    function syncIfromSlider(), app.edtI.Value = app.sldI.Value; end
    function syncIfromEdit(),   app.sldI.Value = app.edtI.Value; end
    function syncQfromSlider(), app.edtQ.Value = app.sldQ.Value; end
    function syncQfromEdit(),   app.sldQ.Value = app.edtQ.Value; end

    function onStart()
        if app.isRunning
            uialert(app.fig,'Simulation already running.','Info'); return;
        end
        app.params.I_cmd       = double(app.edtI.Value);
        app.params.Q_flow_Lpm  = double(app.edtQ.Value);
        app.params.Q_flow_Ls   = app.params.Q_flow_Lpm / 60;
        app.params.t_end       = double(app.edtT.Value) * 60;
        app.params.dt          = double(app.edtDt.Value);
        app.stopRequested = false; app.isRunning = true;
        app.lblStatus.Text = 'Status: Running...';
        run_sim();   % nested function below
    end

    function onStop()
        if ~app.isRunning, app.lblStatus.Text = 'Status: Idle'; return; end
        app.stopRequested = true;
        app.lblStatus.Text = 'Status: Stop requested...';
    end

    function onReset()
        if app.isRunning
            uialert(app.fig,'Stop simulation first to reset.','Warning'); return;
        end
        app.params = init_params();
        app.sldI.Value = app.params.I_cmd;    app.edtI.Value = app.params.I_cmd;
        app.sldQ.Value = app.params.Q_flow_Lpm; app.edtQ.Value = app.params.Q_flow_Lpm;
        app.edtT.Value = app.params.t_end/60; app.edtDt.Value = app.params.dt;
        app.lblNcell.Text = num2str(app.params.Ncell);
        app.lblAcell.Text = num2str(app.params.A_cell*1e4);
        app.lblStatus.Text = 'Status: Reset to defaults'; drawnow;
    end

   

    function run_sim()
        p = app.params; dt = p.dt; Nt = floor(p.t_end/dt)+1; tvec = (0:Nt-1)'*dt;

        % Initial concentrations (mol/L)
        C_V2_t = p.SOC0*p.C_tot; C_V3_t = (1-p.SOC0)*p.C_tot;
        C_V5_t = p.SOC0*p.C_tot; C_V4_t = (1-p.SOC0)*p.C_tot;
        C_V2_e = C_V2_t; C_V3_e = C_V3_t; C_V5_e = C_V5_t; C_V4_e = C_V4_t;

      
        T_stack = p.T0;

       
        V_hist=nan(Nt,1); I_hist=nan(Nt,1); SOC_hist=nan(Nt,1);
        C_t_hist=nan(Nt,4); T_hist=nan(Nt,1); Lambda_hist=nan(Nt,1);

        isCharging=(p.I_cmd<0); mode='CC';
        cla(app.ax1); cla(app.ax2); cla(app.ax3); cla(app.ax4); cla(app.ax5); cla(app.ax6);
        hold(app.ax1,'on'); hold(app.ax2,'on'); hold(app.ax3,'on'); hold(app.ax4,'on'); hold(app.ax5,'on'); hold(app.ax6,'on');

        for k=1:Nt
            if app.stopRequested, app.lblStatus.Text='Status: Stopped by user.'; break; end

          
            Ecell = compute_Ecell(C_V2_e,C_V3_e,C_V4_e,C_V5_e,p);

           
            if strcmp(mode,'CC')
                I_target = p.I_cmd;
                eta = find_eta_for_I(I_target,C_V2_e,C_V3_e,p);
                Vstack = p.Ncell*(Ecell + eta) - I_target*p.R_ohmic; % R_ohmic is STACK-level
                if isCharging && Vstack>=p.V_cut_charge,      mode='CV'; end
                if ~isCharging && Vstack<=p.V_cut_discharge,  mode='CV'; end
            else
                if isCharging, Vtarget=p.V_cut_charge; else, Vtarget=p.V_cut_discharge; end
                I_target = find_I_for_Vtarget(Vtarget,Ecell,C_V2_e,C_V3_e,p);
                Vstack = Vtarget; eta=0;
            end
            if ~isfinite(I_target), I_target=0; eta=0; end

           
            mol_rate = I_target/p.F;  % mol/s (n_e=1)
            dC = -mol_rate/p.V_elec_L; % mol/L/s into electrolyte

            dCv2_e = p.k_mass*(C_V2_t-C_V2_e) + dC;
            dCv3_e = p.k_mass*(C_V3_t-C_V3_e) - dC;
            dCv5_e = p.k_mass*(C_V5_t-C_V5_e) - dC;
            dCv4_e = p.k_mass*(C_V4_t-C_V4_e) + dC;

            dCv2_t = (p.Q_flow_Ls/p.V_tank_L)*(C_V2_e-C_V2_t);
            dCv3_t = (p.Q_flow_Ls/p.V_tank_L)*(C_V3_e-C_V3_t);
            dCv5_t = (p.Q_flow_Ls/p.V_tank_L)*(C_V5_e-C_V5_t);
            dCv4_t = (p.Q_flow_Ls/p.V_tank_L)*(C_V4_e-C_V4_t);

            
            C_V2_e=max(C_V2_e+dCv2_e*dt,0); C_V3_e=max(C_V3_e+dCv3_e*dt,0);
            C_V5_e=max(C_V5_e+dCv5_e*dt,0); C_V4_e=max(C_V4_e+dCv4_e*dt,0);
            C_V2_t=max(C_V2_t+dCv2_t*dt,0); C_V3_t=max(C_V3_t+dCv3_t*dt,0);
            C_V5_t=max(C_V5_t+dCv5_t*dt,0); C_V4_t=max(C_V4_t+dCv4_t*dt,0);

            
            Q_gen = abs(I_target*Vstack); 
            dT = (Q_gen - p.hA*(T_stack-p.Tamb))/(p.C_th_stack);
            T_stack = T_stack + dT*dt;

         
            lambda_nm = p.lambda0 + p.alpha_FBG*(T_stack - p.T0);

            
            V_hist(k)=Vstack; I_hist(k)=I_target;
            SOC_hist(k)=C_V2_t/(C_V2_t+C_V3_t);
            C_t_hist(k,:)=[C_V2_t,C_V3_t,C_V4_t,C_V5_t];
            T_hist(k)=T_stack;
            Lambda_hist(k)=lambda_nm;

            
            if mod(k,round(1/dt))==0 || k==1
                plot(app.ax1,tvec(1:k)/60,V_hist(1:k),'b-'); ylabel(app.ax1,'V_{stack}');
                plot(app.ax2,tvec(1:k)/60,C_t_hist(1:k,:)); ylabel(app.ax2,'Conc (mol/L)');
                legend(app.ax2,{'V2','V3','V4','V5'},'Location','best');
                plot(app.ax3,tvec(1:k)/60,I_hist(1:k),'r-'); ylabel(app.ax3,'I (A)');
                plot(app.ax4,tvec(1:k)/60,SOC_hist(1:k),'k-'); ylabel(app.ax4,'SOC');
                plot(app.ax5,tvec(1:k)/60,T_hist(1:k),'m-'); ylabel(app.ax5,'T (°C)');
                plot(app.ax6,tvec(1:k)/60,Lambda_hist(1:k),'g-'); ylabel(app.ax6,'\lambda (nm)');
                drawnow limitrate;
            end
        end

        % Save to base workspace
        sim.t=tvec; sim.V=V_hist; sim.I=I_hist;
        sim.SOC=SOC_hist; sim.C_t=C_t_hist; sim.T=T_hist; sim.lambda=Lambda_hist;
        sim.params=p;
        assignin('base','vrfb_gui_sim',sim);
        app.sim=sim; app.isRunning=false;
        if ~app.stopRequested
            app.lblStatus.Text='Status: Finished normally. Results saved.';
        end
    end
end


function params=init_params()
    

    params.F=96485; params.R=8.314; params.T=298.15; params.n_e=1;

    
    params.Ncell=200; 
    params.A_cell=0.40; 
    params.A_stack=params.A_cell*params.Ncell;

    
    params.V_tank_L=16700;          % per tank (~16.7 m^3)
    params.V_elec_L=50;             % in channels
    params.Q_flow_Lpm=250; 
    params.Q_flow_Ls=params.Q_flow_Lpm/60;

    
    params.I_cmd=890; 
    params.dt=1; 
    params.t_end=120*60; 
    params.t_end_min=120;

   
    params.C_tot=1.6;               % mol/L
    params.SOC0=0.5; 
    params.E0_cell=1.40;            % OCV at 50% SOC, 25°C
    params.R_ohmic=0.06;            % Ω, STACK-level (~200 * 0.3 mΩ)
    params.k0=3e-3; 
    params.alpha=0.5; 
    params.k_pol=0;

    
    params.V_cut_charge    = params.Ncell*1.65;
    params.V_cut_discharge = params.Ncell*0.90;

    
    params.k_mass = params.Q_flow_Ls/params.V_elec_L;  % s^-1

    
    params.T0=25; 
    params.C_th_stack=1.5e6;        % J/K
    params.hA=640;                  % W/K (h~80 W/m^2K, A~8 m^2)
    params.Tamb=25;

    
    params.lambda0 = 1530;          % nm
    params.alpha_FBG = 0.010;       % nm/°C
end

function E=compute_Ecell(Cv2,Cv3,Cv4,Cv5,p)
    eps=1e-12;
    E = p.E0_cell + (p.R*p.T/(p.n_e*p.F)) * ...
        log((max(Cv5,eps).*max(Cv2,eps))./(max(Cv4,eps).*max(Cv3,eps)));
end

function eta=find_eta_for_I(I_target,C_red,C_ox,p)
    j_target = I_target / p.A_stack;
    f = @(eta) BV_current_density(eta,C_red,C_ox,p) - j_target;
    try
        eta = fzero(f,0);
    catch
        eta = sign(I_target)*0.01;
    end
end

function I_sol=find_I_for_Vtarget(Vtarget,Ecell,C_red,C_ox,p)
    fI = @(I) I - (BV_current_density(Vtarget/p.Ncell - Ecell + (I*p.R_ohmic)/p.Ncell, C_red, C_ox, p) .* p.A_stack);
    try
        I_sol = fzero(fI,0);
    catch
        I_sol = 0;
    end
end

function j=BV_current_density(eta,C_red,C_ox,p)
    Cr=max(C_red,1e-12); Co=max(C_ox,1e-12);
    j0=p.k0*p.F*(Cr.^p.alpha).*(Co.^(1-p.alpha));
    j=j0.*(exp((p.alpha*p.F*eta)/(p.R*p.T)) - exp((-(1-p.alpha)*p.F*eta)/(p.R*p.T)));
end