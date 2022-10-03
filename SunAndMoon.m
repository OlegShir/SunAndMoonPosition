function Sun_and_Moon
SCR=get(0,'ScreenSize');
Lat_start_E=0; Lon_start_E=0;SunAz_start=90; SunEl_start=0;
%% Forms of source infomation
ab=figure('Position',[(SCR(3)-830)/2 (SCR(4)-480)/2 360 480],'Resize',...
    'off','NumberTitle','off','Name', ...
    'Position of the sun and moon','toolbar','none',...
    'MenuBar','none');

panel_left=uipanel(ab,'Title','Initial parameters','TitlePosition',...
    'Centertop','Units','pixels','Position',[10 240 340 240]);
%% Date
uicontrol(panel_left,'Style','text','Position',[10 197 200 20],...
    'String','Date of observation:','HorizontalAlignment','left');

data_string=uicontrol(panel_left,'Style','Edit','Position',...
    [120 202 110 20],'HorizontalAlignment','left',...
    'Enable','off','Tag','dateEdit',...
    'HorizontalAlignment','center');

uicontrol(panel_left,'Style','PushButton','Position',...
    [240 202 90 20],'String','Choose',...
    'HorizontalAlignment','center','Callback',@Calendar);
%% Time
uicontrol(panel_left,'Style','text','Position',[10 170 110 20],'String',...
    'Observation time:','HorizontalAlignment', 'left');

time_hour=uicontrol(panel_left,'Style','Edit','Position',[120 175 30 20],...
    'String',hour(now),'HorizontalAlignment','left','Enable', ...
    'on','HorizontalAlignment','center');

uicontrol(panel_left,'Style','text','Position',[160 170 20 20],...
    'String','h.','HorizontalAlignment', 'left');

time_minutes=uicontrol(panel_left,'Style','Edit','Position',...
    [180 175 30 20],'String',minute(now),'HorizontalAlignment',...
    'left','Enable','on','HorizontalAlignment','center');

uicontrol(panel_left,'Style','text','Position',[220 170 30 20],...
    'String','min.','HorizontalAlignment', 'left');

uicontrol(panel_left,'Style','text','Position',[260 170 30 20],...
    'String','UTC:','HorizontalAlignment', 'left');

UTC=uicontrol(panel_left,'Style','Edit','Position',[300 175 30 20],...
    'String','3','HorizontalAlignment','left','Enable', ...
    'on','HorizontalAlignment','center');

%% Coordinate object
uicontrol(panel_left,'Style','text','Position',[10 143 200 20],...
    'String','Observation point location coordinates:',...
    'HorizontalAlignment', 'left');
% Laties
panel_lat= uipanel(panel_left,'Title','Latitude','TitlePosition',...
    'Centertop','Units','pixels','Position',[10 95 320 45]);

lat_hour=uicontrol(panel_lat,'Style','Edit','Position',[5 5 30 20],...
    'HorizontalAlignment','left','Enable', ...
    'on','String','60','HorizontalAlignment','center');

uicontrol(panel_lat,'Style','text','Position',[40 3 30 20],...
    'String','deg.','HorizontalAlignment', 'left');

lat_min=uicontrol(panel_lat,'Style','Edit','Position',[75 5 30 20],...
    'HorizontalAlignment','left','Enable', ...
    'on','String','0','HorizontalAlignment','center');

uicontrol(panel_lat,'Style','text','Position',[110 3 30 20],...
    'String','min.','HorizontalAlignment', 'left');

lat_sec=uicontrol(panel_lat,'Style','Edit','Position',[145 5 30 20],...
    'HorizontalAlignment','left','Enable', ...
    'on','String','0','HorizontalAlignment','center');

uicontrol(panel_lat,'Style','text','Position',[180 3 30 20],...
    'String','sec.','HorizontalAlignment', 'left');

north_lat=uicontrol(panel_lat,'Style','radiobutton','Enable','on',...
    'Position',[215 7 20 20],'HorizontalAlignment','center',...
    'Value',1, 'Callback',@N_ratio);

uicontrol(panel_lat,'Style','text','Position',[235 3 30 20],...
    'String','north','HorizontalAlignment', 'left')

south_lat=uicontrol(panel_lat,'Style','radiobutton','Enable','on',...
    'Position',[270 7 20 20],'HorizontalAlignment','center',...
    'Callback',@S_ratio);

uicontrol(panel_lat,'Style','text','Position',[290 3 30 20],...
    'String','south','HorizontalAlignment','left');

% Longitude
panel_lon=uipanel(panel_left,'Title','Долгота','TitlePosition','Centertop',...
    'Units','pixels','Position',[10 40 320 45]);

lon_hour=uicontrol(panel_lon,'Style','Edit','Position',[5 5 30 20],...
    'String','30','HorizontalAlignment','left','Enable',...
    'on','HorizontalAlignment','center');

uicontrol(panel_lon,'Style','text','Position',[40 3 30 20],...
    'String','deg.','HorizontalAlignment', 'left');

lon_min=uicontrol(panel_lon,'Style','Edit','Position',[75 5 30 20],...
    'String','0','HorizontalAlignment','left','Enable', ...
    'on','HorizontalAlignment','center');

uicontrol(panel_lon,'Style','text','Position',[110 3 30 20],...
    'String','min.','HorizontalAlignment', 'left');

lon_sec=uicontrol(panel_lon,'Style','Edit','Position',[145 5 30 20],...
    'HorizontalAlignment','left','Enable',...
    'on','String','0','HorizontalAlignment','center');

uicontrol(panel_lon,'Style','text','Position',[180 3 30 20],...
    'String','sec.','HorizontalAlignment', 'left');

east_lon=uicontrol(panel_lon,'Style','radiobutton','Enable','on',...
    'Position',[215 7 20 20],'HorizontalAlignment','center',...
    'Value',1, 'Callback',@E_ratio);

uicontrol(panel_lon,'Style','text','Position',[235 3 30 20],...
    'String','east','HorizontalAlignment', 'left')

west_lon=uicontrol(panel_lon,'Style','radiobutton','Enable','on',...
    'Position',[270 7 20 20],'HorizontalAlignment','center',...
    'Callback',@W_ratio);

uicontrol(panel_lon,'Style','text','Position',[290 3 30 20],...
    'String','west','HorizontalAlignment','left');
%% Buttom's solve
uicontrol(panel_left,'Style','PushButton','Position',[200 10 130 25],...
    'String','Calculation','HorizontalAlignment','center',...
    'Callback',@Solve)
%% Solve
panel_down=uipanel(ab,'Title','Results','TitlePosition',...
    'Centertop','Units','pixels','Position',[10 10 340 230]);


%% Answer Sun

uicontrol(panel_down,'Style','text','Position',[10 194 110 20],...
    'String','Azimuth of the Sun:','HorizontalAlignment','left');
SunAzAns=uicontrol(panel_down,'Style','text','Position',[10 180 70 20],...
    'String','','HorizontalAlignment','left');

uicontrol(panel_down,'Style','text','Position',[10 166 110 20],...
    'String','Sun elevation angle:','HorizontalAlignment','left');
SunElAns=uicontrol(panel_down,'Style','text','Position',[10 152 70 20],...
    'String',"",'HorizontalAlignment','left');

uicontrol(panel_down,'Style','text','Position',[10 138 140 20],...
    'String','Distance to the Sun:','HorizontalAlignment','left');
SunDistAns=uicontrol(panel_down,'Style','text','Position',[10 124 120 20],...
    'String','','HorizontalAlignment','left');

%% Answer panel divide
uipanel(panel_down,'Units','pixels','Position',[0 124 140 2]);

%% Answer Moon
uicontrol(panel_down,'Style','text','Position',[10 98 120 20],...
    'String','Azimuth of the Moon:','HorizontalAlignment','left');
MoonAzAns=uicontrol(panel_down,'Style','text','Position',[10 84 70 20],...
    'String','','HorizontalAlignment','left');

uicontrol(panel_down,'Style','text','Position',[10 70 120 20],...
    'String','Moon elevation angle:','HorizontalAlignment','left');
MoonElAns=uicontrol(panel_down,'Style','text','Position',[10 56 100 20],...
    'String','','HorizontalAlignment','left');

uicontrol(panel_down,'Style','text','Position',[10 42 120 20],...
    'String','Distance to the Moon:','HorizontalAlignment','left');
MoonDistAns=uicontrol(panel_down,'Style','text','Position',[10 28 120 20],...
    'String','','HorizontalAlignment','left');
%% Visualization Earth
axesEarth= axes(panel_down,'Units','pixel','CameraViewAngleMode', 'manual','Position',[125 -120 230 460]);
axis equal
axesm('globe','galt',0)
gridm('glinestyle','-')
m = matfile('topo.mat');
topo = m.topo;
topolegend=m.topolegend;
geo=geoshow(topo,topolegend,'DisplayType','texturemap');
demcmap(topo);
land = shaperead('landareas.shp','UseGeoCoords',true);
linem([land.Lat],[land.Lon])

geo.AmbientStrength = 0.2;
geo.DiffuseStrength = 0.8;
geo.SpecularStrength = 0.5;
geo.SpecularExponent = 1;
geo.BackFaceLighting = 'unlit';
axis off
%camlight(SunAz_start,SunEl_start)
Sun_light=light;
lightangle(Sun_light,90,0);
camposm(Lat_start_E,Lon_start_E,22);
point=geoshow(0,0, 'Visible', 'off');
%% Divide panel_down
uipanel(panel_down,'Units','pixels','Position',[140 -5 2 250]);

%% Visualization
%axes_3D = axes(panel_right,'Units','pixel','Position',[10 10 440 450]);

%%
    function Solve (~,~)
        
        D=day(data_string.String);
        M=month(data_string.String);
        Y=year(data_string.String);
        for k = length(M):-1:1
            if ( M(k) <= 2 ) % january & february
                Y(k)  = Y(k) - 1.0;
                M(k) = M(k) + 12.0;
            end
        end
        UT=str2double(time_hour.String)+str2double(time_minutes.String)/60;
        %number of days in the Gregorian calendar
        jd = floor( 365.25*(Y + 4716.0)) + floor( 30.6001*( M + 1.0)) + 2.0 - ...
            floor( Y/100.0 ) + floor( floor( Y/100.0 )/4.0 ) + D - 1524.5 +(UT)/24 ;
        d = jd - 2451543.5;
        % longitude
        Lon=str2double(lon_hour.String)+str2double(lon_min.String)/60+str2double(lon_sec.String)/3600;
        if west_lon.Value == 1
            Lon=-Lon;
        end
        % latitude
        Lat=str2double(lat_hour.String)+str2double(lat_min.String)/60+str2double(lat_sec.String)/3600;
        if south_lat.Value == 1
            Lat=-Lat;
        end
        oblect = 23.4393-3.563e-7*d;%
        
        %% Find Sun Azimuth and Elevation angles
        %%Keplerian Elements of the Sun
        deg=pi/180;
        w=282.9404+4.70935e-5*d;
        e=0.016713+1.151e-9;
        M=356.0470+0.9856002585*d;
        rv=rev(M);
        L=rev(w+M);
        E=M+(180/pi)*e*sind(M)*(1+e*cosd(M));
        
        %
        x=cosd(E)-e;  y=sind(E)*sqrt(1-e*e); r=sqrt(x*x+y*y);
        v=(arctan2(y,x))*180/pi; lon=rev(v+w);%Longitude of the Sun
        x1=r*cosd(lon); y1=r*sind(lon); z1=0;
        xequart=x1; yequart=y1*cosd(oblect)-z1*sind(oblect);
        zequart=y1*sind(oblect)+z1*cosd(oblect);
        r1=sqrt(xequart^2+yequart^2+zequart^2);
        RA=(arctan2(yequart,xequart))*180/pi; RA=RA/15;
        Decl=(arctan2(zequart,sqrt(xequart^2+yequart^2)))*180/pi;
        
        %Calculate siderial time and hour angle
        G = str2double(UTC.String);% calculation of the difference with greenwich
        GSM0 = L/15 + 12 - G;
        SIDTIME=GSM0+UT+Lon/15;% sidereal time of revolution
        HA=(SIDTIME-RA)*15; %Hour Angle Calculation
        
        %Find the El and Az at the current SIDTIME
        SunEl_end=asind(sind(Decl).*sind(Lat) + cosd(Decl).*cosd(Lat).*cosd(HA));
        SunAz_end=acosd((sind(Decl) - sind(SunEl_end).*sind(Lat))./(cosd(SunEl_end).*cosd(Lat)));
        
        %Answer
        SunAzAns.String=unit_of_value(SunAz_end,1);
        SunElAns.String=unit_of_value(SunEl_end,1);
        format long
        SunDistAns.String=unit_of_value(fix(r1*149597870.700),2);
        
        %Find Moon Azimuth and Elevation angles%
        %Keplerian Elements of the Moon
        Nm=rev(125.1228-0.0529538083*d);%dlinnii uzel
        im=5.1454;
        wm=rev(318.0634+0.1643573223*d);
        am=60.2666;
        em=0.054900;
        Mm=rev(115.3654+13.0649929509*d);
        EarthRadEq = 6378.1370;
        
        %Compute E, the eccentric anomaly
        E0=Mm+(180/pi)*em*sind(Mm)*(1 + em * cosd(Mm));
        E1=E0-(E0-(180/pi)*em*sind(E0)- Mm)/(1-em*cosd(E0));
        Ls=L;
        Lm=rev(Nm + wm + Mm);
        D=rev(Lm-Ls);
        F=Lm-Nm;
        
        %Calculate Lunar perturbations
        longMA=-1.274*sind(Mm-2*D)+0.658*sind(2*D)-0.186*sind(rv)...
            -0.059*sind(2*Mm-2*D)-0.057*sind(Mm-2*D+rv)+0.053*sind(Mm+...
            2*D)+0.046*sind(2*D-rv)+0.041*sind(Mm-rv)...
            -0.035*sind(D)-0.031*sind(Mm+rv)-0.015*sind(2*F-2*D)...
            +0.011*sind(Mm-4*D);
        latMA=-0.173*sind(F-2*D)-0.055*sind(Mm-F-2*D)-0.046*sind(Mm+F-2*D)...
            +0.033*sind(F+2*D)+0.017*sind(2*Mm+F);
        distMA=-0.58*cosd(Mm - 2*D)-0.46*cosd(2*D);
        
        %Compute rectangular coordinates (x,y) in the plane of the lunar orbit
        xm=am*(cosd(E1)-em);
        ym=am*sqrt(1 - em^2)*sin(E1*pi/180);
        
        %convert this to distance and true anomaly
        rm=sqrt(xm^2+ym^2);
        vm=(arctan2(ym,xm))*180/pi;
        
        %Compute moon's position in ecliptic coordinates
        xeclip=rm*(cosd(Nm)*cosd(vm+wm)-sind(Nm)*sind(vm+wm)*cosd(im));
        yeclip=rm*(sind(Nm)*cosd(vm+wm)+cosd(Nm)*sind(vm+wm)*cosd(im));
        zeclip=rm*sind(vm+wm)*sind(im);
        
        %Add the calculated lunar perturbation terms to increase model fidelity
        [longM, latM, distM] = cart2sph(xeclip,yeclip,zeclip);
        [xeclipM, yeclipM, zeclipM] = sph2cart(longM + longMA*deg, ...
            latM + latMA*deg, ...
            distM + distMA);
        clear longM latM;
        
        %Compute ecliptic to equitorial
        xequatM=xeclipM ;
        yequatM=yeclipM*cosd(oblect)-zeclipM*sind(oblect);
        zequatM=yeclipM*sind(oblect)+zeclipM*cosd(oblect);
        RAM=rev(atan2(yequatM,xequatM)/deg);
        DeclM=atan2(zequatM,sqrt(xequatM^2+yequatM^2))/deg;
        HAM=rev((SIDTIME*15-RAM));
        
        %Find the El and Az at the current SIDTIME
        MoonEl=asind(sind(DeclM).*sind(Lat) + cosd(DeclM).*cosd(Lat).*cosd(HAM));
        MoonAz=acosd((sind(DeclM) - sind(MoonEl).*sind(Lat))./(cosd(MoonEl).*cosd(Lat)));
        horParal = 8.794/(rm*6379.14/149.59787e6);
        p = asin(cos(MoonEl.*(pi/180))*sin((horParal/3600).*(pi/180))).*(180/pi);
        MoonEl = MoonEl-p;
        if sind(HAM) >= 0
            MoonAz = 360-MoonAz;
        end
        
        %Answer
        MoonAzAns.String=unit_of_value(MoonAz,1);
        MoonElAns.String=unit_of_value(MoonEl,1);
        MoonDist=fix(EarthRadEq*(distM + distMA));
        MoonDistAns.String=unit_of_value(MoonDist,2);
        %% Geoshow
        
        axes(axesEarth);
        Lat_end=Lat;
        Lon_end=Lon;
        %if Lat_end~=Lat_start&&Lon_end~=Lon_start
        delete(point)
        point=geoshow(Lat_end,Lon_end, 'DisplayType','point','Marker','o','MarkerEdgeColor','r',...
            'MarkerFaceColor', 'r','MarkerSize', 6);
        
        %%
        i_Lat=1;i_Lon=1;i_SanAz=1;i_SanEl=1;
        if SunEl_end<SunEl_start
            i_SanEl=-1;
        end
        if SunAz_end<SunAz_start
            i_SanAz=-1;
        end
        if Lat_end<Lat_start_E
            i_Lat=-1;
        end
        if Lon_end<Lon_start_E
            i_Lon=-1;
        end
        K_Lat=abs(Lat_end-Lat_start_E)/360;
        K_Lon=abs(Lon_end-Lon_start_E)/360;
        K_SunAz=abs(SunAz_end-SunAz_start)/360;
        K_SunEl=abs(SunEl_end-SunEl_start)/360;
        for drive=1:2:360
            Lat_drive=Lat_start_E+i_Lat*K_Lat*drive;
            Lon_drive=Lon_start_E+drive*i_Lon*K_Lon;
            SunAz_drive=SunAz_start+drive*i_SanAz*K_SunAz;
            SunEl_drive=SunEl_start+drive*i_SanEl*K_SunEl;
            camposm(Lat_drive,Lon_drive,22);
            lightangle(Sun_light,SunAz_drive,SunEl_drive);
            drawnow;
        end
        SunAz_start=SunAz_end;
        SunEl_start=SunEl_end;
        Lat_start_E=Lat_end;
        Lon_start_E=Lon_end;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Переключение выбора между Север-Юг и Запад-Восток
    function N_ratio(~,~)
        if north_lat.Value==1
            south_lat.Value=0;
        end
        if north_lat.Value==0
            north_lat.Value=1;
        end
    end
    function S_ratio(~,~)
        if south_lat.Value==1
            north_lat.Value=0;
        end
        if south_lat.Value==0
            south_lat.Value=1;
        end
    end
    function E_ratio(~,~)
        if east_lon.Value==1
            west_lon.Value=0;
        end
        if east_lon.Value==0
            east_lon.Value=1;
        end
    end
    function W_ratio(~,~)
        if west_lon.Value==1
            east_lon.Value=0;
        end
        if west_lon.Value==0
            west_lon.Value=1;
        end
    end
%%
    function Calendar (~,~)
        geo=findobj('Tag','dateEdit');
        uicalendar('Weekend', [0 0 0 0 0 1 1], ...
            'SelectionType', 1, 'DestinationUI', ...
            geo, 'OutputDateFormat', 'dd mmmm yyyy');
    end
%% lawing unit of value
    function [value]=unit_of_value(value, unit)
        if unit==1
            degree_value=degrees2dms(value);
            value=string(degree_value(1)+string(char(176))+' '+degree_value(2)+"' "+fix(degree_value(3))+'"');
        elseif unit==2
            value=string(value)+" км.";
        else
        end
    end
%%



end