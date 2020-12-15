classdef ght_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        imageButton        matlab.ui.control.Button
        shapeButton        matlab.ui.control.Button
        findButton         matlab.ui.control.Button
        Label              matlab.ui.control.Label
        Label2             matlab.ui.control.Label
        Label3             matlab.ui.control.Label
        scoreLabel         matlab.ui.control.Label
        yLabel             matlab.ui.control.Label
        xLabel             matlab.ui.control.Label
        UIAxes             matlab.ui.control.UIAxes
        EdgeDropDownLabel  matlab.ui.control.Label
        EdgeDropDown       matlab.ui.control.DropDown
        sigmaSpinnerLabel  matlab.ui.control.Label
        sigmaSpinner       matlab.ui.control.Spinner
    end

    
    properties (Access = private)
        %Property  Description
%         image=null;
%         shape=null;
%         x=null;
%         y=null;
%         score=null;
        image=[];
        imageGray=[];
        imageEdg=[];
        shape=[];
        x;
        y;
        score;
    end
    
    methods (Access = private)
        
        function [score,  y, x ] = Generalized_hough_transform(app,Is,Iedg,Itm) 
            %Find template/shape Itm in greyscale image Is using generalize hough trasform
            %Use generalized hough transform to find Template/shape binary image given in binary image Itm in grayscale image Is (greyscale image)
            %Is is greyscale  picture were the template Itm should be found 
            %Itm is bool edge image of the template with edges markedd ones
            %Return the x,y location  cordniates  which gave the best match 
            %Also return the score of each this point (number of point matching)
            %The x,y are the cordinates in image Is in which the  the top left edge of image Itm (1,1) should be positioned in order to give the best match
            %Is=imread('');
            %Itm=imread('');
            %if nargin<3 thresh=1;end;
            %--------------------------create edge and system edge images------------------------------------------------------------------------------------------------------------------------
            %Is=rgb2gray(Is);
            %Iedg=edge(Is,'canny'); % Take canny edge images of Is with automatic threshold
            %}
            %--------------------------------------------------------------------------------------------------------------------------------------
            [y, x]=find(Itm>0); % find all y,x cordinates of all points equal 1 inbinary template image Itm
            nvs=size(x);% number of points in the  template image
            %-------------------Define Yc and Xc ----------------------------------------------
            Cy=1;%round(mean(y));% find object y center, note that any reference point will do so the origin of axis hence 1 could be used just as well
            Cx=1;%round(mean(x));% find object z center, note that any reference point will do so the origin of axis hence 1 could be used just as well
            %------------------------------create gradient map of Itm, distrobotion between zero to pi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            GradientMap = gradient_direction( app,Itm );
            %%%%%%%%%%%%%%%%%%%%%%%Create an R-Table of Itm gradients to  parameter space in parameter space.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %---------------------------create template descriptor array------------------------------------
            MaxAngelsBins=30;% devide the angel space to MaxAngelsBins uniformed space bins
            MaxPointsPerangel=nvs(1);% maximal amount of points corresponding to specific angel
            PointCounter=zeros(MaxAngelsBins);% counter for the amount of edge points associate with each angel gradient
            Rtable=zeros(MaxAngelsBins,MaxPointsPerangel,2); % assume maximum of 100 points per angle with MaxAngelsBins angles bins between zero and pi and x,y for the vector to the center of each point
            % the third adimension are vector from the point to the center of the vessel
            %------------------fill the  angel bins with points in the Rtable---------------------------------------------------------
            for f=1:1:nvs(1)
                bin=round((GradientMap(y(f), x(f))/pi)*(MaxAngelsBins-1))+1; % transform from continues gradient angles to discrete angle bins and 
                PointCounter(bin)=PointCounter(bin)+1;% add one to the number of points in the bin
                if (PointCounter(bin)>MaxPointsPerangel)
                    disp('exceed max bin in hugh transform');
                end
                Rtable(bin, PointCounter(bin),1)= Cy-y(f);% add the vector from the point to the object center to the bin
                Rtable(bin, PointCounter(bin),2)= Cx-x(f);% add the vector from the point to the object center to the bin
            end
            %plot(pc);
            %pause;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%create and populate hough space%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %-----------------------------use the array in previous image to identify the template in the main image Is----------------------------------------
            [y, x]=find(Iedg>0); % find all edg point in the in edge image Iedg of the main image Is
            np=size(x);% find number of edge points Is edge image
            GradientMap=gradient_direction(app,Is); % create gradient direction  map of the Is
            Ss=size(Is); % Size of the main image Is
            houghspace=zeros(size(Is));% the hough space assume to be in size of the image but it should probably be smaller
                for f=1:1:np(1)
                      bin=round((GradientMap(y(f), x(f))/pi)*(MaxAngelsBins-1))+1; % transform from continues gradient angles to discrete angle bins and 
                      for fb=1:1:PointCounter(bin)
                          ty=Rtable(bin, fb,1)+ y(f);
                          tx=Rtable(bin, fb,2)+ x(f);
                           if (ty>0) && (ty<Ss(1)) && (tx>0) && (tx<Ss(2))  
                               houghspace(Rtable(bin, fb,1)+ y(f), Rtable(bin, fb,2)+ x(f))=  houghspace(Rtable(bin, fb,1)+ y(f), Rtable(bin, fb,2)+ x(f))+1; % add point in were the center of the image should be according to the pixel gradient
                           end
                      end
                end
            %{
            %====================================show The Hough Space in color==================================================================================================
            imtool(houghspace);
            imshow(houghspace,[]);
            colormap jet
            colorbar
            pause
            %}
            %============================================Find best match in hough space=========================================================================================
            %---------------------------------------------------------------------------normalized according to template size (fraction of the template points that was found)------------------------------------------------------------------------------------------------
            Itr=houghspace;%./(sum(sum(Itm))); % Itr become the new score matrix
            %---------------------------------------------------------------------------find  the location best score all scores which are close enough to the best score
            %imtool(Itr,[]);
            mx=max(max(Itr));% find the max score location
            [y,x]=find(Itr==mx,1,'first');
            %[y,x]=find(Itr>=thresh*mx,  10, 'first'); % find the location first 10 best matches which their score is at least thresh percents of the maximal score and pot them in the x,y array
            score=Itr(y,x); % find max score in the huogh space 
        end
        
        function [ Is ] = gradient_direction(~, i3 )
            %return map of the absolute direction from -pi/2 to pi/2  of gradient in every point of the gradient  i3(only half circle does not have negative directionss
            %-------------------------------------------------------------------
            Dy=imfilter(double(i3),[1; -1],'same');%x first derivative  sobel mask
            Dx=imfilter(double(i3),[1  -1],'same');% y sobel first derivative
            %Is=atan2(Dy,Dx)+pi();
            Is=mod(atan2(Dy,Dx)+pi(), pi());%atan(Dy/Dx);%note that this expression can reach infinity if dx is zero mathlab aparently get over it but you can use the folowing expression instead slower but safer: 
            %mod(atan2(Dy,Dx)+pi(), pi());%gradient direction map going from 0-180
            %--------------------show the image-----------------------------------------------
            %{
            imshow(Is,[]);% the ,[]  make sure the display will be feeted to doube image
            colormap jet
            colorbar
            pause;
            %}
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.findButton.Enable='off';
            app.UIAxes.Visible='off';
            app.sigmaSpinner.Enable='off';
        end

        % Button pushed function: imageButton
        function imageButtonPushed(app, event)
            [file,path]=uigetfile('*.jpg;*.bmp;*.gif;*.tiff;*.png');
            if isequal(file,0)
                return
            else
                a=imread(fullfile(path,file));
                %app.Image.ImageSource=a;
                app.image=a;
                imtool(a);
                a=rgb2gray(a);
                %app.image=double(a)/255;
                app.imageGray=a;
                value = app.EdgeDropDown.Value;
                app.imageEdg=edge(app.imageGray,value);
                imshow(app.imageEdg,'parent',app.UIAxes);
                %app.Image.ImageSource=app.image;
                if ~isempty(app.shape)
                    app.findButton.Enable='on';
                end
            end
        end

        % Button pushed function: shapeButton
        function shapeButtonPushed(app, event)
            [file,path]=uigetfile('*.jpg;*.bmp;*.gif;*.tiff;*.png');
             if isequal(file,0)
                return
            else
                a=imread(fullfile(path,file));
                app.shape=(a(:,:,1)==0);
                imtool(app.shape);
                if ~isempty(app.imageGray)
                    app.findButton.Enable='on';
                end
             end
        end

        % Button pushed function: findButton
        function findButtonPushed(app, event)
            [app.score,app.y,app.x]=Generalized_hough_transform(app,app.imageGray,app.imageEdg,app.shape);
            app.scoreLabel.Text=string(app.score);
            app.xLabel.Text=string(app.x);
            app.yLabel.Text=string(app.y);
            [a,b,~]=size(app.image);
            temp=zeros([a,b,3],'uint8');
            for i=1:a
                for j=1:b
                    if (i-app.x+1>0 && j-app.y+1>0 && i-app.x+1<=size(app.shape,1) && j-app.y+1<=size(app.shape,2) && app.shape(i-app.x+1,j-app.y+1))
                        temp(i,j,1)=255;
                    else
                        temp(i,j,:)=app.image(i,j,:);
                    end
                end
            end
            
            imshow(temp,'parent',app.UIAxes);
        end

        % Value changed function: EdgeDropDown
        function EdgeDropDownValueChanged(app, event)
            value = app.EdgeDropDown.Value;
            if(value=="log")
                app.sigmaSpinner.Value=2;
                app.sigmaSpinner.Enable='on';
                sigma=app.sigmaSpinner.Value;
                app.imageEdg=edge(app.imageGray,value,0,sigma);
            elseif(value=="Canny")
                app.sigmaSpinner.Value=1.4;
                app.sigmaSpinner.Enable='on';
                sigma=app.sigmaSpinner.Value;
                app.imageEdg=edge(app.imageGray,value,[],sigma);
            else
                app.sigmaSpinner.Enable='off';
                app.imageEdg=edge(app.imageGray,value);
            end
            imshow(app.imageEdg,'parent',app.UIAxes);
        end

        % Value changed function: sigmaSpinner
        function sigmaSpinnerValueChanged(app, event)
            value = app.EdgeDropDown.Value;
            if(value=="log")
                sigma=app.sigmaSpinner.Value;
                app.imageEdg=edge(app.imageGray,value,0,sigma);
            elseif(value=="Canny")
                sigma=app.sigmaSpinner.Value;
                app.imageEdg=edge(app.imageGray,value,[],sigma);
            end
            imshow(app.imageEdg,'parent',app.UIAxes);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create imageButton
            app.imageButton = uibutton(app.UIFigure, 'push');
            app.imageButton.ButtonPushedFcn = createCallbackFcn(app, @imageButtonPushed, true);
            app.imageButton.Position = [18 388 100 22];
            app.imageButton.Text = 'image';

            % Create shapeButton
            app.shapeButton = uibutton(app.UIFigure, 'push');
            app.shapeButton.ButtonPushedFcn = createCallbackFcn(app, @shapeButtonPushed, true);
            app.shapeButton.Position = [407 388 100 22];
            app.shapeButton.Text = 'shape';

            % Create findButton
            app.findButton = uibutton(app.UIFigure, 'push');
            app.findButton.ButtonPushedFcn = createCallbackFcn(app, @findButtonPushed, true);
            app.findButton.Position = [517 388 100 22];
            app.findButton.Text = 'find';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [70 64 35 22];
            app.Label.Text = 'score';

            % Create Label2
            app.Label2 = uilabel(app.UIFigure);
            app.Label2.Position = [161 64 25 22];
            app.Label2.Text = 'y';

            % Create Label3
            app.Label3 = uilabel(app.UIFigure);
            app.Label3.Position = [241 64 25 22];
            app.Label3.Text = 'x';

            % Create scoreLabel
            app.scoreLabel = uilabel(app.UIFigure);
            app.scoreLabel.Position = [104 64 25 22];
            app.scoreLabel.Text = '';

            % Create yLabel
            app.yLabel = uilabel(app.UIFigure);
            app.yLabel.Position = [185 64 44 22];
            app.yLabel.Text = '';

            % Create xLabel
            app.xLabel = uilabel(app.UIFigure);
            app.xLabel.Position = [265 64 44 22];
            app.xLabel.Text = '';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
            app.UIAxes.Position = [1 85 640 291];

            % Create EdgeDropDownLabel
            app.EdgeDropDownLabel = uilabel(app.UIFigure);
            app.EdgeDropDownLabel.HorizontalAlignment = 'right';
            app.EdgeDropDownLabel.Position = [128 388 34 22];
            app.EdgeDropDownLabel.Text = 'Edge';

            % Create EdgeDropDown
            app.EdgeDropDown = uidropdown(app.UIFigure);
            app.EdgeDropDown.Items = {'Sobel', 'Prewitt', 'Roberts', 'log', 'Canny'};
            app.EdgeDropDown.ValueChangedFcn = createCallbackFcn(app, @EdgeDropDownValueChanged, true);
            app.EdgeDropDown.Position = [177 388 100 22];
            app.EdgeDropDown.Value = 'Sobel';

            % Create sigmaSpinnerLabel
            app.sigmaSpinnerLabel = uilabel(app.UIFigure);
            app.sigmaSpinnerLabel.HorizontalAlignment = 'right';
            app.sigmaSpinnerLabel.Position = [289 388 33 22];
            app.sigmaSpinnerLabel.Text = 'sigma';

            % Create sigmaSpinner
            app.sigmaSpinner = uispinner(app.UIFigure);
            app.sigmaSpinner.Step = 0.1;
            app.sigmaSpinner.Limits = [0 Inf];
            app.sigmaSpinner.ValueChangedFcn = createCallbackFcn(app, @sigmaSpinnerValueChanged, true);
            app.sigmaSpinner.Position = [334 388 63 22];
            app.sigmaSpinner.Value = 2;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ght_exported

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