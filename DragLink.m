classdef DragLink < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        LengthwithcouplermmEditField    matlab.ui.control.NumericEditField
        LengthwithcouplermmEditFieldLabel  matlab.ui.control.Label
        LengthwithcrankmmEditField      matlab.ui.control.NumericEditField
        LengthwithcrankmmEditFieldLabel  matlab.ui.control.Label
        ThetawithcrankdegEditField      matlab.ui.control.NumericEditField
        ThetawithcrankdegEditFieldLabel  matlab.ui.control.Label
        ThetawithcouplerdegEditField    matlab.ui.control.NumericEditField
        ThetawithcouplerdegEditFieldLabel  matlab.ui.control.Label
        ForceoncrankNEditField          matlab.ui.control.NumericEditField
        ForceoncrankNEditFieldLabel     matlab.ui.control.Label
        ForceoncouplerNEditField        matlab.ui.control.NumericEditField
        ForceoncouplerNEditFieldLabel   matlab.ui.control.Label
        AngvelocityradsecEditField      matlab.ui.control.NumericEditField
        AngvelocityradsecEditFieldLabel  matlab.ui.control.Label
        ANIMATIONOFDRAGLINKMECHANISMLabel  matlab.ui.control.Label
        TorqueNmEditField               matlab.ui.control.NumericEditField
        TorqueNmEditFieldLabel          matlab.ui.control.Label
        ResetButton                     matlab.ui.control.Button
        exitButton                      matlab.ui.control.Button
        LengthwithRockermmEditField     matlab.ui.control.NumericEditField
        LengthwithRockermmEditFieldLabel  matlab.ui.control.Label
        StaticButton                    matlab.ui.control.Button
        ThetawithRockerdegEditField     matlab.ui.control.NumericEditField
        ThetawithRockerdegEditFieldLabel  matlab.ui.control.Label
        ForceonRockerNEditField         matlab.ui.control.NumericEditField
        ForceonRockerNEditFieldLabel    matlab.ui.control.Label
        Theta2degEditField              matlab.ui.control.NumericEditField
        Theta2degEditFieldLabel         matlab.ui.control.Label
        TimevariationsecEditField       matlab.ui.control.NumericEditField
        TimevariationsecEditFieldLabel  matlab.ui.control.Label
        StaticForceAnalysisLabel        matlab.ui.control.Label
        saveButton                      matlab.ui.control.Button
        MechanismEditField              matlab.ui.control.EditField
        MechanismEditFieldLabel         matlab.ui.control.Label
        MinAngleEditField               matlab.ui.control.NumericEditField
        MinAngleEditFieldLabel          matlab.ui.control.Label
        MaxAngleEditField               matlab.ui.control.NumericEditField
        MaxAngleEditFieldLabel          matlab.ui.control.Label
        TransmissionButton              matlab.ui.control.Button
        TransmissionAnglesLabel         matlab.ui.control.Label
        AccelerationButton              matlab.ui.control.Button
        velocityButton                  matlab.ui.control.Button
        AnimationButton                 matlab.ui.control.Button
        CalculateButton                 matlab.ui.control.Button
        EnterLinkLengthsLabel           matlab.ui.control.Label
        L3mmEditField                   matlab.ui.control.NumericEditField
        L3mmEditFieldLabel              matlab.ui.control.Label
        L2mmEditField                   matlab.ui.control.NumericEditField
        L2mmEditFieldLabel              matlab.ui.control.Label
        L1mmEditField                   matlab.ui.control.NumericEditField
        L1mmEditFieldLabel              matlab.ui.control.Label
        L0mmEditField                   matlab.ui.control.NumericEditField
        L0mmEditFieldLabel              matlab.ui.control.Label
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
      theta_max
      theta_min
      theta2
      theta3
      theta4
      theta1
      theta
      r1
      r2
      r3
      r4
      Lengths
      sorted_lengths
      L_max
      L_min
      l1
      l2
      time
      t
      Torque
      omega3
      omega4
      alpha3
      alpha4
      delt
      timer
      ang_speed
      
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.r1=app.L0mmEditField.Value;
            app.r2=app.L1mmEditField.Value;
            app.r3=app.L2mmEditField.Value;
            app.r4=app.L3mmEditField.Value;
        
        end

        % Button pushed function: TransmissionButton
        function TransmissionButtonPushed(app, event)
            app.r1=app.L0mmEditField.Value;
            app.r2=app.L1mmEditField.Value;
            app.r3=app.L2mmEditField.Value;
            app.r4=app.L3mmEditField.Value;
            max_d=app.r1+app.r2
            min_d=app.r1-app.r2  
            app.theta_max=acosd((app.r3^2+app.r4^2-max_d^2)/(2*app.r3*app.r4))
            app.theta_min=acosd((app.r3^2+app.r4^2-min_d^2)/(2*app.r3*app.r4))
            app.MinAngleEditField.Value=app.theta_min;
            app.MaxAngleEditField.Value=app.theta_max;
            
        end

        % Button pushed function: AnimationButton
        function AnimationButtonPushed(app, event)

             app.delt=app.TimevariationsecEditField.Value;
             app.r1=app.L0mmEditField.Value;
             app.r2=app.L1mmEditField.Value;
             app.r3=app.L2mmEditField.Value;
             app.r4=app.L3mmEditField.Value;
             app.Lengths = [app.r1 app.r2 app.r3 app.r4];
             app.sorted_lengths = sort(app.Lengths);
             app.L_max = app.sorted_lengths(end);
             app.L_min = app.sorted_lengths(1);
             app.l1 = app.sorted_lengths(2);
             app.l2 = app.sorted_lengths(3);
          
            if (app.L_max + app.L_min <app.l1+app.l2)&&(app.L_min==app.Lengths(1))
              app.ang_speed =app.AngvelocityradsecEditField.Value; 
              app.timer=(2*3.14)/(app.ang_speed);
           
              app.t=0:app.delt:app.timer;
              app.theta2 = app.ang_speed*app.t;
              app.theta1 = deg2rad(0);
             for i=1:length(app.theta2)

              theta2_prime(i) =app.theta2(i) - app.theta1;
              delta(i) = sqrt(app.r1^2 + app.r2^2 -2*app.r1*app.r2*cos(theta2_prime(i)));
              beta(i) = acos( (app.r1^2 + delta(i)^2 - app.r2^2) / (2*app.r1*delta(i)));
              psi(i) = acos( (app.r3^2 + delta(i)^2 - app.r4^2) / (2*app.r3*delta(i)));
              lambda(i) = acos( (app.r4^2 + delta(i)^2 - app.r3^2) / (2*app.r4*delta(i)));

             if(theta2_prime<=pi)
              app.theta3(i) = psi(i)-(beta(i)-app.theta1);
              app.theta4(i) = pi-lambda(i)-(beta(i)-app.theta1);
              gamma(i)= acos( (app.r3^2+app.r4^2-delta(i)^2) / (2*app.r3*app.r4)) - pi/2;
             else
              app.theta3(i) = psi(i)+(beta(i)+app.theta1);
              app.theta4(i) = pi-lambda(i)+(beta(i)+app.theta1);
              gamma(i)= acos( (app.r3^2+app.r4^2-delta(i)^2) / (2*app.r3*app.r4)) - pi/2;
             end
             Ax(i) = app.r2*cos(app.theta2(i));
             Ay(i) = app.r2*sin(app.theta2(i));
             Bx(i) = app.r2*cos(app.theta2(i))+app.r3*cos(app.theta3(i));
             By(i) = app.r2*sin(app.theta2(i))+app.r3*sin(app.theta3(i));
             Box(i) = app.r1*cos(app.theta1);
             Boy(i) = app.r1*sin(app.theta1);
             plot(app.UIAxes,[0 Ax(i)], [0 Ay(i)],'ro-','LineWidth',5);
             hold(app.UIAxes,"on"); %r2
             plot(app.UIAxes,[Ax(i) Bx(i)], [Ay(i) By(i)], 'go-','LineWidth',5);
             hold(app.UIAxes,"on"); %r3
             plot(app.UIAxes,[Bx(i) Box(i)], [By(i) Boy(i)], 'bo-','LineWidth',5); 
             hold(app.UIAxes,"on"); %r4
             plot(app.UIAxes,[Box(i) 0], [Boy(i) 0], 'co-','LineWidth',5);
             hold(app.UIAxes,"off");%r1
             %grid(app.UIAxes,"on")
             legend(app.UIAxes,'crank[L2]','Coupler[L3]','Rocker[L4]','Ground[L0]');
             xlim(app.UIAxes,[-800,800]);
             ylim(app.UIAxes,[-800,800]);
             pause(0.001);
             end 
             app.theta3
             app.theta4
             app.omega3 = (app.theta3(2:end) -app.theta3(1:end-1))./0.5;
             app.alpha3 = (app.omega3(2:end) -app.omega3(1:end-1))./0.5;
             app.omega4 = (app.theta4(2:end) -app.theta4(1:end-1))./0.5;
             app.alpha4 = (app.omega4(2:end) -app.omega4(1:end-1))./0.5;
            end
        end

        % Button pushed function: velocityButton
        function velocityButtonPushed(app, event)
            comet(app.t(2:end),app.omega3);
            xlabel('Time(seconds)');
            ylabel('Angular velocity of Link 3 (rad/sec)');
            
        end

        % Button pushed function: CalculateButton
        function CalculateButtonPushed(app, event)
            app.r1=app.L0mmEditField.Value;
            app.r2=app.L1mmEditField.Value;
            app.r3=app.L2mmEditField.Value;
            app.r4=app.L3mmEditField.Value;
            Lengths = [app.r1 app.r2 app.r3 app.r4];
            sorted_lengths = sort(Lengths);
            L_max = sorted_lengths(end);
            L_min = sorted_lengths(1);
            l1 = sorted_lengths(2);
            l2 = sorted_lengths(3);
            if L_max >= L_min+l1+l2
                   msgbox("Linkage is not possible","error");
            return
            else
               app.MechanismEditField.Value="Linkage is possibe"
    
            end
            if L_max + L_min > l1+l2
              msgbox("TRIPLE ROCKER (NON-GRASHOF)","error")
            else
              app.MechanismEditField.Value="GRASHOF MECHANISM"
            if L_max + L_min == l1+l2
               msgbox("CHANGE POINT",'error');
            elseif L_min == Lengths(2)
               msgbox("CRANK ROCKER",'error');
            elseif L_min == Lengths(3)
              msgbox("DOUBLE ROCKER",'error');
            else 
              app.MechanismEditField.Value="DRAG LINK";
                
            end
         end
        end

        % Button pushed function: AccelerationButton
        function AccelerationButtonPushed(app, event)
         comet(app.t(3:end),app.alpha3);
          xlabel('Time(seconds)');
          ylabel('Angular Acceleration of Link3(rad/sec^2)');
          %grid(app.UIAxes3,"on")
        end

        % Button pushed function: saveButton
        function saveButtonPushed(app, event)
            app.r1=app.L0mmEditField.Value;
            app.r2=app.L1mmEditField.Value;
            app.r3=app.L2mmEditField.Value;
            app.r4=app.L3mmEditField.Value;
            [file,path] = uiputfile('*.txt');
            filename = fullfile(path,file);
            fileID = fopen(file,'w+');
            s='Input Link Lengths';
            su='Static Force Analysis';
            Tr='Transmission angles are';
            fprintf(fileID,'\n %s\n\n',s);
            fprintf(fileID,'Length of ground = %d mm\n',app.r1);
            fprintf(fileID,'Length of crank  = %d mm\n',app.r2);  
            fprintf(fileID,'Length of coupler= %d mm\n',app.r3);
            fprintf(fileID,'Length of rocker = %d mm\n\n\n\n',app.r4);
            fprintf(fileID,'%s\n',Tr);
            fprintf(fileID,'Maximum Angle is %12.6f\n',app.theta_max);
            fprintf(fileID,'Manimum Angle is %12.6f\n\n',app.theta_min);
            fprintf(fileID,'%s\n',su);
            fprintf(fileID,'Torque acting the Link2 is :%12.6f N/m\n\n',app.Torque/1000);
          
            t1=deg2rad(0);
            app.ang_speed =app.AngvelocityradsecEditField.Value; 
            app.timer=(2*3.14)/(app.ang_speed);
            app.t=0:app.delt:app.timer;
            app.theta2=app.t*app.ang_speed;
            alpha2 = 0;
            for i = 1:length(app.theta2)
            syms T3 T4 o3 o4 a3 a4
            T2 = app.theta2(i);
            eq1 = app.r1 + app.r4*cos(T4) == app.r2*cos(T2) + app.r3*cos(T3);
            eq2 = app.r4*sin(T4) == app.r2*sin(T2) + app.r3*sin(T3);
            eq3 = app.r4*o4*sin(T4) == app.r2*app.ang_speed*sin(T2) + app.r3*o3*sin(T3);
            eq4 = app.r4*o4*cos(T4) == app.r2*app.ang_speed*cos(T2) + app.r3*o3*cos(T3);
            eq5 = -app.r4*a4*sin(T4)-app.r4*power(o4,2)*cos(T4) == -app.r2*alpha2*sin(T2) - app.r2*power(app.ang_speed,2)*cos(T2) - app.r3*a3*sin(T3) - app.r3*power(o3,2)*cos(T3);
            eq6 = -app.r4*a4*cos(T4)-app.r4*power(o4,2)*sin(T4) == app.r2*alpha2*cos(T2) - app.r2*power(app.ang_speed,2)*sin(T2) + app.r3*a3*cos(T3) - app.r3*power(o3,2)*sin(T3);
            [T3,T4,o3,o4,a3,a4] = vpasolve(eq1,eq2,eq3,eq4,eq5,eq6,[T3,T4,o3,o4,a3,a4]);
            t3(i)=T3(1);
            t4(i) = T4(1);
            t3dot(i) = o3(1);
            t4dot(i) = o4(1);
            t3ddot(i) = a3(1);
            t4ddot(i) = a4(1);
            end
            omg=[app.t; app.theta2; t3dot; t3ddot; t4dot; t4ddot];
            fprintf(fileID,'%6s %12s %12s %12s %12s %12s\n','Time(seconds)','theta2','omega3(rad/sec)','omega4(rad/sec)','alpha3(rad/sec^2)','alpha4(rad/sec^2)');
            fprintf(fileID,'%12.8f    %6.2f        %12.8f    %12.8f      %12.8f     %12.8f\n',omg);
            fclose(fileID);
            
            
        end

        % Button pushed function: StaticButton
        function StaticButtonPushed(app, event)
            syms FR Fc Fcr TR Tc Tcr LR Lc Lcr fbc fdc fab theta_3 theta_4 t3 t4
            app.r1=app.L0mmEditField.Value;
            app.r2=app.L1mmEditField.Value;
            app.r3=app.L2mmEditField.Value;
            app.r4=app.L3mmEditField.Value;
            theta_2=app.Theta2degEditField.Value;
            FR=app.ForceonRockerNEditField.Value;
            LR=app.LengthwithRockermmEditField.Value;
            TR=app.ThetawithRockerdegEditField.Value;
             Lc=app.LengthwithcouplermmEditField.Value;
             Fc=app.ForceoncouplerNEditField.Value;
            Tc=app.ThetawithcouplerdegEditField.Value;
            Fcr=app.ForceoncrankNEditField.Value;
            Lcr=app.LengthwithcrankmmEditField.Value;
            Tcr=app.ThetawithcrankdegEditField.Value;
            t1=0;
            eqn1=app.r2*cosd(theta_2)+app.r3*cosd(theta_3)-app.r1*cosd(t1)-app.r4*cosd(theta_4)
            eqn2=app.r2*sind(theta_2)+app.r3*sind(theta_3)-app.r1*sind(t1)-app.r4*sind(theta_4)
            [t4,t3]=vpasolve([eqn1,eqn2],[theta_4,theta_3]);
            t4
            t3
            %Analysis on link4
            syms fbc
            FBC=fbc*[cosd(t3) sind(t3) 0];
            rcd=[app.r4*cosd(t4) app.r4*sind(t4) 0 ];
            rde=LR*[cosd(t4) sind(t4) 0 ];
            F1=FR*[cosd(TR) sind(TR)  0];
            eqn1=cross(rcd,FBC);
            eqn2=cross(rde,F1);
            eqn=eqn1+eqn2==0;
            [Fbc]=vpasolve(eqn,fbc);
            rab=[app.r2*cosd(theta_2) app.r2*sind(theta_2) 0 ];
            FBC=Fbc*[cosd(t3) sind(t3) 0];
            Torque1=double(cross(rab,FBC));
            
            %Analysis on Link3
            FDC=fdc*[cosd(t4) sind(t4) 0];
            rbc=[app.r2*cosd(t3) app.r3*sind(t3) 0 ];
            rbe=[Lc*cosd(t3) Lc*sind(t3) 0 ];
            F2=Fc*[cosd(Tc) sind(Tc)  0];
            eqn3=cross(rbc,FDC);
            eqn4=cross(rbe,F2);
            eqn=eqn3+eqn4==0;
            [Fdc]=vpasolve(eqn,fdc);
            rab=[app.r2*cosd(theta_2) app.r2*sind(theta_2) 0 ];
            FDC=Fdc*[cosd(t4) sind(t4) 0];
            Torque2=double(cross(rab,FDC));
            
            %Analysis on Link2
            syms fbc
            FBC=fbc*[cos(t3) sin(t3) 0];
            F3=Fcr*[cosd(Tcr) sind(Tcr) 0];
            rab=[app.r2*cosd(theta_2) app.r2*sind(theta_2) 0 ];
            rae=[Lcr*cosd(theta_2) Lcr*sind(theta_2) 0 ];
            eqn5=cross(rae,F3);
            eqn6=cross(rab,FBC);
            eqn=eqn5+eqn6==0;
            [Fbc]=vpasolve(eqn,fbc);
            FBC=Fbc*[cos(t3) sin(t3) 0];
            rab=[app.r2*cosd(theta_2) app.r2*sind(theta_2) 0 ];
            Torque3=double(cross(rab,FBC));
            
            app.Torque=(Torque1(3)+Torque2(3)+Torque3(3));
            app.TorqueNmEditField.Value=app.Torque/1000;
            message = sprintf('Torque value of the link2 =%d N/m',app.Torque);
            uiwait(msgbox(message, 'Torque')); 
        end

        % Button pushed function: exitButton
        function exitButtonPushed(app, event)
            closereq();
        end

        % Button pushed function: ResetButton
        function ResetButtonPushed(app, event)
            app.L0mmEditField.Value=0;
            app.L1mmEditField.Value=0;
            app.L2mmEditField.Value=0;
            app.L3mmEditField.Value=0;
            app.TimevariationsecEditField.Value=0;
            app.ThetawithRockerdegEditField.Value=0;
            app.Theta2degEditField.Value=0;
            app.ForceonRockerNEditField.Value=0;
            app.LengthwithRockermmEditField.Value=0;
            app.MechanismEditField.Value='';
            app.MaxAngleEditField.Value=0;
            app.MinAngleEditField.Value=0;
            app.AngvelocityradsecEditField.Value=0;
            app.ForceonRockerNEditField.Value=0;
            app.ForceoncouplerNEditField.Value=0;
            app.ForceoncrankNEditField.Value=0;
            app.ThetawithRockerdegEditField.Value=0;
            app.ThetawithcouplerdegEditField.Value=0;
            app.ThetawithcrankdegEditField.Value=0;
            app.LengthwithRockermmEditField.Value=0;
            app.LengthwithcouplermmEditField.Value=0;
            app.LengthwithcrankmmEditField.Value=0;
            app.TorqueNmEditField.Value=0;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1380 801];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [296 117 747 617];

            % Create L0mmEditFieldLabel
            app.L0mmEditFieldLabel = uilabel(app.UIFigure);
            app.L0mmEditFieldLabel.HorizontalAlignment = 'right';
            app.L0mmEditFieldLabel.FontSize = 17;
            app.L0mmEditFieldLabel.Position = [82 584 64 22];
            app.L0mmEditFieldLabel.Text = 'L0(mm)';

            % Create L0mmEditField
            app.L0mmEditField = uieditfield(app.UIFigure, 'numeric');
            app.L0mmEditField.FontSize = 17;
            app.L0mmEditField.Position = [191 584 72 22];

            % Create L1mmEditFieldLabel
            app.L1mmEditFieldLabel = uilabel(app.UIFigure);
            app.L1mmEditFieldLabel.HorizontalAlignment = 'right';
            app.L1mmEditFieldLabel.FontSize = 17;
            app.L1mmEditFieldLabel.Position = [81 552 64 22];
            app.L1mmEditFieldLabel.Text = 'L1(mm)';

            % Create L1mmEditField
            app.L1mmEditField = uieditfield(app.UIFigure, 'numeric');
            app.L1mmEditField.FontSize = 17;
            app.L1mmEditField.Position = [191 552 72 22];

            % Create L2mmEditFieldLabel
            app.L2mmEditFieldLabel = uilabel(app.UIFigure);
            app.L2mmEditFieldLabel.HorizontalAlignment = 'right';
            app.L2mmEditFieldLabel.FontSize = 17;
            app.L2mmEditFieldLabel.Position = [78 520 64 22];
            app.L2mmEditFieldLabel.Text = 'L2(mm)';

            % Create L2mmEditField
            app.L2mmEditField = uieditfield(app.UIFigure, 'numeric');
            app.L2mmEditField.FontSize = 17;
            app.L2mmEditField.Position = [191 522 72 22];

            % Create L3mmEditFieldLabel
            app.L3mmEditFieldLabel = uilabel(app.UIFigure);
            app.L3mmEditFieldLabel.HorizontalAlignment = 'right';
            app.L3mmEditFieldLabel.FontSize = 17;
            app.L3mmEditFieldLabel.Position = [81 490 64 22];
            app.L3mmEditFieldLabel.Text = 'L3(mm)';

            % Create L3mmEditField
            app.L3mmEditField = uieditfield(app.UIFigure, 'numeric');
            app.L3mmEditField.FontSize = 17;
            app.L3mmEditField.Position = [191 490 72 22];

            % Create EnterLinkLengthsLabel
            app.EnterLinkLengthsLabel = uilabel(app.UIFigure);
            app.EnterLinkLengthsLabel.FontSize = 18;
            app.EnterLinkLengthsLabel.FontWeight = 'bold';
            app.EnterLinkLengthsLabel.FontColor = [1 0.0745 0.651];
            app.EnterLinkLengthsLabel.Position = [63 645 168 40];
            app.EnterLinkLengthsLabel.Text = 'Enter Link Lengths';

            % Create CalculateButton
            app.CalculateButton = uibutton(app.UIFigure, 'push');
            app.CalculateButton.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.CalculateButton.FontSize = 17;
            app.CalculateButton.Position = [105 363 99 29];
            app.CalculateButton.Text = 'Calculate';

            % Create AnimationButton
            app.AnimationButton = uibutton(app.UIFigure, 'push');
            app.AnimationButton.ButtonPushedFcn = createCallbackFcn(app, @AnimationButtonPushed, true);
            app.AnimationButton.Position = [854 48 125 31];
            app.AnimationButton.Text = 'Animation';

            % Create velocityButton
            app.velocityButton = uibutton(app.UIFigure, 'push');
            app.velocityButton.ButtonPushedFcn = createCallbackFcn(app, @velocityButtonPushed, true);
            app.velocityButton.Position = [1059 50 107 32];
            app.velocityButton.Text = 'velocity';

            % Create AccelerationButton
            app.AccelerationButton = uibutton(app.UIFigure, 'push');
            app.AccelerationButton.ButtonPushedFcn = createCallbackFcn(app, @AccelerationButtonPushed, true);
            app.AccelerationButton.Position = [1213 50 108 35];
            app.AccelerationButton.Text = 'Acceleration';

            % Create TransmissionAnglesLabel
            app.TransmissionAnglesLabel = uilabel(app.UIFigure);
            app.TransmissionAnglesLabel.FontSize = 18;
            app.TransmissionAnglesLabel.FontWeight = 'bold';
            app.TransmissionAnglesLabel.FontColor = [1 0.0745 0.651];
            app.TransmissionAnglesLabel.Position = [63 291 186 22];
            app.TransmissionAnglesLabel.Text = 'Transmission Angles';

            % Create TransmissionButton
            app.TransmissionButton = uibutton(app.UIFigure, 'push');
            app.TransmissionButton.ButtonPushedFcn = createCallbackFcn(app, @TransmissionButtonPushed, true);
            app.TransmissionButton.FontSize = 17;
            app.TransmissionButton.Position = [93 159 116 29];
            app.TransmissionButton.Text = 'Transmission';

            % Create MaxAngleEditFieldLabel
            app.MaxAngleEditFieldLabel = uilabel(app.UIFigure);
            app.MaxAngleEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxAngleEditFieldLabel.FontSize = 17;
            app.MaxAngleEditFieldLabel.Position = [50 240 91 22];
            app.MaxAngleEditFieldLabel.Text = ' Max-Angle';

            % Create MaxAngleEditField
            app.MaxAngleEditField = uieditfield(app.UIFigure, 'numeric');
            app.MaxAngleEditField.FontSize = 17;
            app.MaxAngleEditField.Position = [155 240 70 22];

            % Create MinAngleEditFieldLabel
            app.MinAngleEditFieldLabel = uilabel(app.UIFigure);
            app.MinAngleEditFieldLabel.HorizontalAlignment = 'right';
            app.MinAngleEditFieldLabel.FontSize = 17;
            app.MinAngleEditFieldLabel.Position = [59 208 82 22];
            app.MinAngleEditFieldLabel.Text = 'Min-Angle';

            % Create MinAngleEditField
            app.MinAngleEditField = uieditfield(app.UIFigure, 'numeric');
            app.MinAngleEditField.FontSize = 17;
            app.MinAngleEditField.Position = [157 208 69 22];

            % Create MechanismEditFieldLabel
            app.MechanismEditFieldLabel = uilabel(app.UIFigure);
            app.MechanismEditFieldLabel.HorizontalAlignment = 'right';
            app.MechanismEditFieldLabel.FontSize = 18;
            app.MechanismEditFieldLabel.Position = [335 50 114 33];
            app.MechanismEditFieldLabel.Text = 'Mechanism';

            % Create MechanismEditField
            app.MechanismEditField = uieditfield(app.UIFigure, 'text');
            app.MechanismEditField.FontSize = 16;
            app.MechanismEditField.Position = [461 40 286 46];

            % Create saveButton
            app.saveButton = uibutton(app.UIFigure, 'push');
            app.saveButton.ButtonPushedFcn = createCallbackFcn(app, @saveButtonPushed, true);
            app.saveButton.Position = [23 70 97 30];
            app.saveButton.Text = 'save';

            % Create StaticForceAnalysisLabel
            app.StaticForceAnalysisLabel = uilabel(app.UIFigure);
            app.StaticForceAnalysisLabel.FontSize = 19;
            app.StaticForceAnalysisLabel.FontWeight = 'bold';
            app.StaticForceAnalysisLabel.FontColor = [1 0.0745 0.651];
            app.StaticForceAnalysisLabel.Position = [1159 705 198 29];
            app.StaticForceAnalysisLabel.Text = 'Static Force Analysis';

            % Create TimevariationsecEditFieldLabel
            app.TimevariationsecEditFieldLabel = uilabel(app.UIFigure);
            app.TimevariationsecEditFieldLabel.HorizontalAlignment = 'right';
            app.TimevariationsecEditFieldLabel.FontSize = 17;
            app.TimevariationsecEditFieldLabel.Position = [13 454 151 22];
            app.TimevariationsecEditFieldLabel.Text = 'Time-variation(sec)';

            % Create TimevariationsecEditField
            app.TimevariationsecEditField = uieditfield(app.UIFigure, 'numeric');
            app.TimevariationsecEditField.FontSize = 17;
            app.TimevariationsecEditField.Position = [191 454 72 22];

            % Create Theta2degEditFieldLabel
            app.Theta2degEditFieldLabel = uilabel(app.UIFigure);
            app.Theta2degEditFieldLabel.HorizontalAlignment = 'right';
            app.Theta2degEditFieldLabel.FontSize = 19;
            app.Theta2degEditFieldLabel.Position = [1156 653 114 24];
            app.Theta2degEditFieldLabel.Text = 'Theta 2(deg)';

            % Create Theta2degEditField
            app.Theta2degEditField = uieditfield(app.UIFigure, 'numeric');
            app.Theta2degEditField.FontSize = 19;
            app.Theta2degEditField.Position = [1285 653 72 24];

            % Create ForceonRockerNEditFieldLabel
            app.ForceonRockerNEditFieldLabel = uilabel(app.UIFigure);
            app.ForceonRockerNEditFieldLabel.HorizontalAlignment = 'right';
            app.ForceonRockerNEditFieldLabel.FontSize = 19;
            app.ForceonRockerNEditFieldLabel.Position = [1098 622 172 24];
            app.ForceonRockerNEditFieldLabel.Text = 'Force on Rocker(N)';

            % Create ForceonRockerNEditField
            app.ForceonRockerNEditField = uieditfield(app.UIFigure, 'numeric');
            app.ForceonRockerNEditField.FontSize = 19;
            app.ForceonRockerNEditField.Position = [1285 622 72 24];

            % Create ThetawithRockerdegEditFieldLabel
            app.ThetawithRockerdegEditFieldLabel = uilabel(app.UIFigure);
            app.ThetawithRockerdegEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetawithRockerdegEditFieldLabel.FontSize = 19;
            app.ThetawithRockerdegEditFieldLabel.Position = [1069 497 204 24];
            app.ThetawithRockerdegEditFieldLabel.Text = 'Theta with Rocker(deg)';

            % Create ThetawithRockerdegEditField
            app.ThetawithRockerdegEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThetawithRockerdegEditField.FontSize = 19;
            app.ThetawithRockerdegEditField.Position = [1285 497 72 24];

            % Create StaticButton
            app.StaticButton = uibutton(app.UIFigure, 'push');
            app.StaticButton.ButtonPushedFcn = createCallbackFcn(app, @StaticButtonPushed, true);
            app.StaticButton.FontSize = 19;
            app.StaticButton.Position = [1190 228 133 31];
            app.StaticButton.Text = 'Static';

            % Create LengthwithRockermmEditFieldLabel
            app.LengthwithRockermmEditFieldLabel = uilabel(app.UIFigure);
            app.LengthwithRockermmEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthwithRockermmEditFieldLabel.FontSize = 19;
            app.LengthwithRockermmEditFieldLabel.Position = [1057 282 213 24];
            app.LengthwithRockermmEditFieldLabel.Text = 'Length with Rocker(mm)';

            % Create LengthwithRockermmEditField
            app.LengthwithRockermmEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengthwithRockermmEditField.FontSize = 19;
            app.LengthwithRockermmEditField.Position = [1285 282 72 24];

            % Create exitButton
            app.exitButton = uibutton(app.UIFigure, 'push');
            app.exitButton.ButtonPushedFcn = createCallbackFcn(app, @exitButtonPushed, true);
            app.exitButton.Position = [82 20 94 22];
            app.exitButton.Text = 'exit';

            % Create ResetButton
            app.ResetButton = uibutton(app.UIFigure, 'push');
            app.ResetButton.ButtonPushedFcn = createCallbackFcn(app, @ResetButtonPushed, true);
            app.ResetButton.Position = [145 70 103 31];
            app.ResetButton.Text = 'Reset';

            % Create TorqueNmEditFieldLabel
            app.TorqueNmEditFieldLabel = uilabel(app.UIFigure);
            app.TorqueNmEditFieldLabel.HorizontalAlignment = 'right';
            app.TorqueNmEditFieldLabel.FontSize = 19;
            app.TorqueNmEditFieldLabel.Position = [1087 162 111 24];
            app.TorqueNmEditFieldLabel.Text = 'Torque(N/m)';

            % Create TorqueNmEditField
            app.TorqueNmEditField = uieditfield(app.UIFigure, 'numeric');
            app.TorqueNmEditField.FontSize = 19;
            app.TorqueNmEditField.Position = [1212 161 108 24];

            % Create ANIMATIONOFDRAGLINKMECHANISMLabel
            app.ANIMATIONOFDRAGLINKMECHANISMLabel = uilabel(app.UIFigure);
            app.ANIMATIONOFDRAGLINKMECHANISMLabel.FontSize = 24;
            app.ANIMATIONOFDRAGLINKMECHANISMLabel.FontWeight = 'bold';
            app.ANIMATIONOFDRAGLINKMECHANISMLabel.FontColor = [0 0 1];
            app.ANIMATIONOFDRAGLINKMECHANISMLabel.Position = [447 732 476 45];
            app.ANIMATIONOFDRAGLINKMECHANISMLabel.Text = 'ANIMATION OF DRAGLINK MECHANISM';

            % Create AngvelocityradsecEditFieldLabel
            app.AngvelocityradsecEditFieldLabel = uilabel(app.UIFigure);
            app.AngvelocityradsecEditFieldLabel.HorizontalAlignment = 'right';
            app.AngvelocityradsecEditFieldLabel.FontSize = 17;
            app.AngvelocityradsecEditFieldLabel.Position = [11 422 165 22];
            app.AngvelocityradsecEditFieldLabel.Text = 'Ang-velocity(rad/sec)';

            % Create AngvelocityradsecEditField
            app.AngvelocityradsecEditField = uieditfield(app.UIFigure, 'numeric');
            app.AngvelocityradsecEditField.FontSize = 17;
            app.AngvelocityradsecEditField.Position = [193 422 70 22];

            % Create ForceoncouplerNEditFieldLabel
            app.ForceoncouplerNEditFieldLabel = uilabel(app.UIFigure);
            app.ForceoncouplerNEditFieldLabel.HorizontalAlignment = 'right';
            app.ForceoncouplerNEditFieldLabel.FontSize = 19;
            app.ForceoncouplerNEditFieldLabel.Position = [1096 584 175 24];
            app.ForceoncouplerNEditFieldLabel.Text = 'Force on coupler(N)';

            % Create ForceoncouplerNEditField
            app.ForceoncouplerNEditField = uieditfield(app.UIFigure, 'numeric');
            app.ForceoncouplerNEditField.FontSize = 19;
            app.ForceoncouplerNEditField.Position = [1285 584 72 24];

            % Create ForceoncrankNEditFieldLabel
            app.ForceoncrankNEditFieldLabel = uilabel(app.UIFigure);
            app.ForceoncrankNEditFieldLabel.HorizontalAlignment = 'right';
            app.ForceoncrankNEditFieldLabel.FontSize = 19;
            app.ForceoncrankNEditFieldLabel.Position = [1111 540 159 24];
            app.ForceoncrankNEditFieldLabel.Text = 'Force on crank(N)';

            % Create ForceoncrankNEditField
            app.ForceoncrankNEditField = uieditfield(app.UIFigure, 'numeric');
            app.ForceoncrankNEditField.FontSize = 19;
            app.ForceoncrankNEditField.Position = [1285 540 72 24];

            % Create ThetawithcouplerdegEditFieldLabel
            app.ThetawithcouplerdegEditFieldLabel = uilabel(app.UIFigure);
            app.ThetawithcouplerdegEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetawithcouplerdegEditFieldLabel.FontSize = 19;
            app.ThetawithcouplerdegEditFieldLabel.Position = [1064 457 206 24];
            app.ThetawithcouplerdegEditFieldLabel.Text = 'Theta with coupler(deg)';

            % Create ThetawithcouplerdegEditField
            app.ThetawithcouplerdegEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThetawithcouplerdegEditField.FontSize = 19;
            app.ThetawithcouplerdegEditField.Position = [1285 457 72 24];

            % Create ThetawithcrankdegEditFieldLabel
            app.ThetawithcrankdegEditFieldLabel = uilabel(app.UIFigure);
            app.ThetawithcrankdegEditFieldLabel.HorizontalAlignment = 'right';
            app.ThetawithcrankdegEditFieldLabel.FontSize = 19;
            app.ThetawithcrankdegEditFieldLabel.Position = [1076 414 190 24];
            app.ThetawithcrankdegEditFieldLabel.Text = 'Theta with crank(deg)';

            % Create ThetawithcrankdegEditField
            app.ThetawithcrankdegEditField = uieditfield(app.UIFigure, 'numeric');
            app.ThetawithcrankdegEditField.FontSize = 19;
            app.ThetawithcrankdegEditField.Position = [1285 414 72 24];

            % Create LengthwithcrankmmEditFieldLabel
            app.LengthwithcrankmmEditFieldLabel = uilabel(app.UIFigure);
            app.LengthwithcrankmmEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthwithcrankmmEditFieldLabel.FontSize = 19;
            app.LengthwithcrankmmEditFieldLabel.Position = [1066 366 200 24];
            app.LengthwithcrankmmEditFieldLabel.Text = 'Length with crank(mm)';

            % Create LengthwithcrankmmEditField
            app.LengthwithcrankmmEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengthwithcrankmmEditField.FontSize = 19;
            app.LengthwithcrankmmEditField.Position = [1285 366 72 24];

            % Create LengthwithcouplermmEditFieldLabel
            app.LengthwithcouplermmEditFieldLabel = uilabel(app.UIFigure);
            app.LengthwithcouplermmEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthwithcouplermmEditFieldLabel.FontSize = 19;
            app.LengthwithcouplermmEditFieldLabel.Position = [1054 324 216 24];
            app.LengthwithcouplermmEditFieldLabel.Text = 'Length with coupler(mm)';

            % Create LengthwithcouplermmEditField
            app.LengthwithcouplermmEditField = uieditfield(app.UIFigure, 'numeric');
            app.LengthwithcouplermmEditField.FontSize = 19;
            app.LengthwithcouplermmEditField.Position = [1285 324 72 24];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DragLink

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end