classdef ght_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure           matlab.ui.Figure
        imageButton        matlab.ui.control.Button
        shapeButton        matlab.ui.control.Button
        findButton         matlab.ui.control.Button
        Label              matlab.ui.control.Label
        UIAxes             matlab.ui.control.UIAxes
        EdgeDropDownLabel  matlab.ui.control.Label
        EdgeDropDown       matlab.ui.control.DropDown
        sigmaSpinnerLabel  matlab.ui.control.Label
        sigmaSpinner       matlab.ui.control.Spinner
        stepSpinnerLabel   matlab.ui.control.Label
        stepSpinner        matlab.ui.control.Spinner
        findrotateButton   matlab.ui.control.Button
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
            Itr=houghspace/(sum(sum(Itm))); % Itr become the new score matrix
            %---------------------------------------------------------------------------find  the location best score all scores which are close enough to the best score
            %imtool(Itr,[]);
            mx=max(max(Itr));% find the max score location
            [y,x]=find(Itr==mx);
            %[y,x]=find(Itr==mx,1,'first');
            %[y,x]=find(Itr>=thresh*mx,  10, 'first'); % find the location first 10 best matches which their score is at least thresh percents of the maximal score and pot them in the x,y array
            %score=Itr(y,x); % find max score in the huogh space 
            score=mx;
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
        
        function [ItmAng, nAng, BestScore]= ght_with_rotate(app,Is,Iedg,Itm)
            %{
            Find an object that fit Template Itm in image Is.
            The orientation of the template and the object in the image does not have to be the same as that as the template. 
            The template Itm is matched to the image Is in various of rotations and the best match is chosen. 
            The function use  Generalized Hough transforms to match the template to the
            image.
            Input (Essential):
            Is: Color image with the object to be found.
            Itm: Template of the object to be found. The template is written as binary image with the boundary of the template marked 1(white) and all the rest of the pixels marked 0. 
            Output
            Ismarked: The image (Isresize) with the template marked upon it in the location of and size of the best match.
            Iborders: Binary image of the borders of the template/object in its scale and located in the place of the best match on the image. 
            Ybest Xbest: location on the image (in pixels) were the template were found to give the best score (location the top left pixel of the template Itm in the image).
            ItmAng: The rotation angle of  the template in degrees that give the best match
            BestScore: Score of the best match found in the scan (the score of the output).
            Algorithm:
            The function rotate the template Itm in various of angles and for each rotation search for the template in the image. 
            The angle and location in the image that gave the best match for the template are chosen.
            %}
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initialize optiona paramters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if (nargin<1)  
%                 Is=imread('Is.jpg');  
%             end %Read image
%             if (nargin<2)  
%                 Itm=imread('Itm.tif');
%             end %Read template image
            %Is=rgb2gray(Is);
            %Itm=logical(Itm);% make sure Itm is boolean image
            BestScore=-100000;
            nAng=1;
            ItmAng=zeros(1,100);
%             close all;
%             imtool close all;
            %%%%%%%%%%%%%%%%Some parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555555555
    %         St=size(Itm);
    %         Ss=size(Is);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Main Scan  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for Ang=0:app.stepSpinner.Value:359 % rotate the template Itm 1 degree at the time and look for it in the image Is
                
                %disp([num2str((Ang)/3.6) '% Scanned']);
              Itr=Rotate_binary_edge_image(app,Itm,Ang);
            %----------------------------------------------------------------------------------------------------------------------------------------- 
             % the actuall recogniton step of the resize template Itm in the orginal image Is and return location of best match and its score can occur in one of three modes given in search_mode
                         [score,  ~,~ ]=Generalized_hough_transform(app,Is,Iedg,Itr);% use generalized hough transform to find the template in the image
                 %--------------------------if the correct match score is better then previous best match write the paramter of the match as the new best match------------------------------------------------------
                 if (score>BestScore) % if item  result scored higher then the previous result
                    BestScore=score;% remember best score
                    nAng=1;
                    ItmAng(:)=0;
                    ItmAng(1)=Ang;
                 elseif(score==BestScore)
                    nAng=nAng+1;
                    ItmAng(nAng)=Ang;
                 end
            %-------------------------------mark best found location on image----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%show   best match optional part can be removed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if BestScore>-100000% Display best match
%                  Itr=Rotate_binary_edge_image(app,Itm,ItmAng);
%                         [yy,xx] =find(Itr);
%                          Ismarked=set2(app,Is,[yy,xx],255,Ybest,Xbest);%Mark best match on image
%                         imshow(Ismarked,'Parent',app.UIAxes);
%                         Iborders=false(size(Is));
%                    Iborders=set2(app,Iborders,[yy,xx],1,Ybest,Xbest);
%                
%             else % if no match 
%                disp('Error no match founded');
%                 Ismarked=0;% assign arbitary value to avoid 
%                    %Iborders=0;
%                    Iborders=0;
%                 
%             end
        end
        
        function [mat]=Rotate_binary_edge_image(app,I,Ang)
        %{
        Rotate binary (logical) edge image (I) in (Ang) degrees 
        The rotated output image will also also be a binary edge image.
        The connectivity/topology of all edges/curves in the input image (I) will be maintained and the line thickness of the curves
         in the output image (mat) will remain 1 pixel.
        The center of rotation is the center of the image
        The dimensions of the output image (mat) will be different from the 
        input image (I) and will be set such that the rotated image is fully within the image frame.
        Input
        I: Binary  edge image (logical type) consist of lines and curves with 
        thickness of one pixel (such as curves, contour line, template, or edge images)
        Ang: Rotation angle of the image in Degrees
         
        Output
        mat: Rotated version of the input image (I), also binary edge image, the connectivity/topology of the edges/curves in input image  (I)
        is maintained and also the line thickness remain one pixel.
        x' = x*cos(theta) - y*sin (theta); rotated coordinates
        y' = x*sin(theta) + y*cos (theta);
        %}
        [Hight,Width]=size(I);
        CntX=(Width+1)/2;
        CntY=(Hight+1)/2; %find center of orignal image
        theta=Ang/180*pi; % convert angle from degree  to radians
        %============================find size of new image===========================================================================================
        Rot=[   cos(theta) -sin(theta)
                sin(theta)  cos(theta)];
        Corners=[ 1-CntY     1-CntX
                  1-CntY     Width-CntX
                  Hight-CntY 1-CntX
                  Hight-CntY Width-CntX];  %Corners of the old image  center at zero
        NewCorn=Corners*Rot;   
        NMaxX=max(NewCorn(:,2));% corners of new image
        NMinX=min(NewCorn(:,2));
        NMaxY=max(NewCorn(:,1));
        NMinY=min(NewCorn(:,1));
        NSizeY=NMaxY-NMinY; %Size of new image Y
        NSizeX=NMaxX-NMinX; %Size of new image Y
        Dx=-NMinX;          % Tanslation of new image x
        Dy=-NMinY;          % Tanslation of new image y
        mat = false(round(NSizeY+1),round(NSizeX+1));% create canvas for rotated image
        %===============================================rotate image==================================================================================================
        [y,x]=find(I);
        [n,~]=size(x);
        for f=1:n %scan all points on image
             p1=[y(f)-CntY x(f)-CntX]*Rot+[Dy+1 Dx+1];
             mat(round(p1(1)),round(p1(2)))=1;% mark rotated point
            
            
            if x(f)>1 && I(y(f),x(f)-1)==1 % Connect horizontal neighbour 
                   Nx=x(f)-1; Ny=y(f);
                   p1=[y(f)-CntY, x(f)-CntX]*Rot+[Dy+1 Dx+1];
                   p2=[Ny-CntY, Nx-CntX]*Rot+[Dy+1 Dx+1];
                   mat=ConnectPoints(app,mat,round(p1(2)),round(p1(1)),round(p2(2)),round(p2(1)));
            end
            if y(f)>1 && I(y(f)-1,x(f))==1 %connect vertical neighbour
                    Nx=x(f); Ny=y(f)-1;  
                    p1=[y(f)-CntY, x(f)-CntX]*Rot+[Dy+1 Dx+1];
                    p2=[Ny-CntY, Nx-CntX]*Rot+[Dy+1 Dx+1];
                    mat=ConnectPoints(app,mat,round(p1(2)),round(p1(1)),round(p2(2)),round(p2(1)));
            end
            if y(f)>1  && x(f)>1 &&  I(y(f)-1,x(f))==0 && I(y(f),x(f)-1)==0 && I(y(f)-1,x(f)-1)==1 % connect diagonal neighbor
                    Nx=x(f)-1; Ny=y(f)-1;  
                    p1=[y(f)-CntY, x(f)-CntX]*Rot+[Dy+1 Dx+1];
                    p2=[Ny-CntY, Nx-CntX]*Rot+[Dy+1 Dx+1];
                    mat=ConnectPoints(app,mat,round(p1(2)),round(p1(1)),round(p2(2)),round(p2(1)));
            end
            if y(f)>1  && x(f)<Width &&  I(y(f)-1,x(f))==0 && I(y(f),x(f)+1)==0 && I(y(f)-1,x(f)+1)==1 % connect second diagonal neighbor
                    Nx=x(f)+1; Ny=y(f)-1;  
                    p1=[y(f)-CntY, x(f)-CntX]*Rot+[Dy+1 Dx+1];
                    p2=[Ny-CntY, Nx-CntX]*Rot+[Dy+1 Dx+1];
                    mat=ConnectPoints(app,mat,round(p1(2)),round(p1(1)),round(p2(2)),round(p2(1)));
            end 
            %}
          
        end
        %imshow(mat);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [mat]=ConnectPoints(~,mat,x1,y1,x2,y2)
            % draw line between x1,y1 to x2,y2 on binary image mat
            %make sure cordinates are inside image mat
            SizeMat=size(mat);
            x1=max(x1,1);x2=max(x2,1);y1=max(y1,1);y2=max(y2,1);
            y1=min(y1,SizeMat(1));y2=min(y2,SizeMat(1));x1=min(x1,SizeMat(2));x2=min(x2,SizeMat(2));
            if x1==x2 % for horizontal line
               yy=linspace(y1,y2,abs(y1-y2)*5);
               xx = ones(size(yy))*x1;
            else % for none horizontal line
               xx = linspace(x1,x2,abs(x1-x2)*5+abs(y1-y2)*5);   
               yy = y1+(xx-x1)*(y2-y1)/(x2-x1);                     
            end
            round(yy);
           index = sub2ind(size(mat),round(yy),round(xx));  
            mat(index) = 1;   
           %imshow(mat)
           %pause(0.1);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            app.findButton.Enable='off';
            app.findrotateButton.Enable='off';
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
                %imtool(a);
                a=rgb2gray(a);
                %app.image=double(a)/255;
                app.imageGray=a;
                value = app.EdgeDropDown.Value;
                app.imageEdg=edge(app.imageGray,value);
                imshow(app.imageEdg,'parent',app.UIAxes);
                %app.Image.ImageSource=app.image;
                if ~isempty(app.shape)
                    app.findButton.Enable='on';
                    app.findrotateButton.Enable='on';
                end
                app.Label.Text="image added";
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
                %imtool(app.shape);
                if ~isempty(app.imageGray)
                    app.findButton.Enable='on';
                    app.findrotateButton.Enable='on';
                end
                app.Label.Text="shape added";
             end
        end

        % Button pushed function: findButton
        function findButtonPushed(app, event)
            app.Label.Text="shape search in progress";
            tic
            [score,y,x]=Generalized_hough_transform(app,app.imageGray,app.imageEdg,app.shape);
            time=toc;
            s=size(x,1);
            [a,b,~]=size(app.image);
            temp=zeros([a,b,3],'uint8');
            temp(:,:,:)=app.image(:,:,:);
            for k=1:s
                for i=1:size(app.shape,1)
                    for j=1:size(app.shape,2)
                        if (app.shape(i,j))
                            temp(i+y(k),j+x(k),1)=255;
                            temp(i+y(k),j+x(k),2)=0;
                            temp(i+y(k),j+x(k),3)=0;
                        end
                    end
                end
            end
            app.Label.Text="found "+string(s)+" shapes, score "+string(score)+', time: '+string(time);
            imshow(temp,'Parent',app.UIAxes);
        end

        % Value changed function: EdgeDropDown
        function EdgeDropDownValueChanged(app, event)
            app.Label.Text="edge search in progress";
            value = app.EdgeDropDown.Value;
            if(value=="log")
                app.sigmaSpinner.Value=2;
                app.sigmaSpinner.Enable='on';
                sigma=app.sigmaSpinner.Value;
                app.imageEdg=edge(app.imageGray,value,[],sigma);
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
            app.Label.Text="edges found";
        end

        % Value changed function: sigmaSpinner
        function sigmaSpinnerValueChanged(app, event)
            app.Label.Text="edge search in progress";
            value = app.EdgeDropDown.Value;
            sigma=app.sigmaSpinner.Value;
            app.imageEdg=edge(app.imageGray,value,[],sigma);
            imshow(app.imageEdg,'parent',app.UIAxes);
            app.Label.Text="edges found";
        end

        % Button pushed function: findrotateButton
        function findrotateButtonPushed(app, event)
            app.Label.Text="shape search in progress";
            tic
            [Ang,nAng,score]=ght_with_rotate(app,app.imageGray,app.imageEdg,app.shape);
            time=toc;
            s=0;
            [a,b,~]=size(app.image);
            temp=zeros([a,b,3],'uint8');
            temp(:,:,:)=app.image(:,:,:);
            for l=1:nAng
                Itr=Rotate_binary_edge_image(app,app.shape,Ang(l));
                [~,y,x]=Generalized_hough_transform(app,app.imageGray,app.imageEdg,Itr);
                s1=size(x,1);
                s=s+s1;
%                 zmienna="iteracja"
%                 Ang
%                 size(Ang)
%                 score
%                 y
%                 x
%                 s
%                 s1
                for k=1:s1
                    for i=1:size(Itr,1)
                        for j=1:size(Itr,2)
                            if (Itr(i,j))
                                temp(i+y(k),j+x(k),1)=255;
                                temp(i+y(k),j+x(k),2)=0;
                                temp(i+y(k),j+x(k),3)=0;
                            end
                        end
                    end
                end
            end
            app.Label.Text="found "+string(s)+" shapes, score: "+string(score)+', nAng: '+string(nAng)+', time: '+string(time);
            imshow(temp,'Parent',app.UIAxes);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 826 530];
            app.UIFigure.Name = 'MATLAB App';

            % Create imageButton
            app.imageButton = uibutton(app.UIFigure, 'push');
            app.imageButton.ButtonPushedFcn = createCallbackFcn(app, @imageButtonPushed, true);
            app.imageButton.Position = [18 480 100 22];
            app.imageButton.Text = 'image';

            % Create shapeButton
            app.shapeButton = uibutton(app.UIFigure, 'push');
            app.shapeButton.ButtonPushedFcn = createCallbackFcn(app, @shapeButtonPushed, true);
            app.shapeButton.Position = [391 480 100 22];
            app.shapeButton.Text = 'shape';

            % Create findButton
            app.findButton = uibutton(app.UIFigure, 'push');
            app.findButton.ButtonPushedFcn = createCallbackFcn(app, @findButtonPushed, true);
            app.findButton.Position = [498 480 100 22];
            app.findButton.Text = 'find';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [42 37 388 22];
            app.Label.Text = '';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, '')
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
            app.UIAxes.Position = [1 82 826 373];

            % Create EdgeDropDownLabel
            app.EdgeDropDownLabel = uilabel(app.UIFigure);
            app.EdgeDropDownLabel.HorizontalAlignment = 'right';
            app.EdgeDropDownLabel.Position = [117 480 34 22];
            app.EdgeDropDownLabel.Text = 'Edge';

            % Create EdgeDropDown
            app.EdgeDropDown = uidropdown(app.UIFigure);
            app.EdgeDropDown.Items = {'Sobel', 'Prewitt', 'Roberts', 'log', 'Canny'};
            app.EdgeDropDown.ValueChangedFcn = createCallbackFcn(app, @EdgeDropDownValueChanged, true);
            app.EdgeDropDown.Position = [166 480 100 22];
            app.EdgeDropDown.Value = 'Sobel';

            % Create sigmaSpinnerLabel
            app.sigmaSpinnerLabel = uilabel(app.UIFigure);
            app.sigmaSpinnerLabel.HorizontalAlignment = 'right';
            app.sigmaSpinnerLabel.Position = [276 480 33 22];
            app.sigmaSpinnerLabel.Text = 'sigma';

            % Create sigmaSpinner
            app.sigmaSpinner = uispinner(app.UIFigure);
            app.sigmaSpinner.Step = 0.1;
            app.sigmaSpinner.Limits = [0.1 Inf];
            app.sigmaSpinner.ValueChangedFcn = createCallbackFcn(app, @sigmaSpinnerValueChanged, true);
            app.sigmaSpinner.Position = [321 480 63 22];
            app.sigmaSpinner.Value = 2;

            % Create stepSpinnerLabel
            app.stepSpinnerLabel = uilabel(app.UIFigure);
            app.stepSpinnerLabel.HorizontalAlignment = 'right';
            app.stepSpinnerLabel.Position = [605 480 28 22];
            app.stepSpinnerLabel.Text = 'step';

            % Create stepSpinner
            app.stepSpinner = uispinner(app.UIFigure);
            app.stepSpinner.Limits = [1 360];
            app.stepSpinner.Position = [648 480 63 22];
            app.stepSpinner.Value = 2;

            % Create findrotateButton
            app.findrotateButton = uibutton(app.UIFigure, 'push');
            app.findrotateButton.ButtonPushedFcn = createCallbackFcn(app, @findrotateButtonPushed, true);
            app.findrotateButton.Position = [716 480 100 22];
            app.findrotateButton.Text = 'find rotate';

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