%Code written by Sophie Rothman
%Supplied as part of the supplement to "Waterfalls alter reach-scale
%fluvial erosion rates: Evidence from field data and process modeling" by
%Sophie D. Rothman, Joel S. Scheingross and Scott W. McCoy, in review at
%JGR:Earth Surface.
function[allwidth, equations, allsymnrm, P, cuts_index, WF_concatenated, wf_freq, P_savg, width_avg, wf_top, wf_bot, dist_abv_in_wf, dist_blw_wf, dist_in_wf]=width_calc_tribs_interp_nochute(cDEM, FD, S, elev_measure, interval, windowsize, k, avgby, thresh_slope, thresh_dh, long_top, lat_top, long_bot, lat_bot )%S is a stream object, must be trunk only right now
%INPUTS

%cDEM=DEM (1 m lidar) of the area of interest
%FD=TopoToolbox-derived FLOWobj of the area of interest
%S=STREAMobj of the stream of itnerest (may have tributaries)
%elev_measure=1D array of the elevations above the thalweg that the user
%   wants to measure
%interval=density of measurements, minimum is 1 (where 1 creates a 
%   measurement every node), whole integers only
%windowsize=integer, controlls the window created for analyzing each cross
%   section - 50 is often good, and will create a square of 100x100. If the
%   stream is small, a smaller integer can be used that will increase the
%   speed of the code, but may result in more NaN values if the width of the
%   channel exceeds the window
%long_top=longitude in UTM of the top of the reach in question. A null
%   value will trigger measurement of the full river
%lat_top=lattitude in UTM of the top of the reach in question. A null
%   value will trigger measurement of the full river
%long_bot=longitude in UTM of the bottom of the reach in question. A null
%   value will trigger measurement of the full river
%lat_bot=lattitude in UTM of the bottom of the reach in question. A null
%   value will trigger measurement of the full river

%OUTPUTS
y = smooth(S,S.y,'k',50,'nstribs',true); %smooth the planform by a 50 m window
x = smooth(S,S.x,'k',50,'nstribs',true);
n=5; %.5 interval across smoothing

ixc = getnal(S);
ixc(S.ix) = S.ixc;
A = flowacc(FD).*(FD.cellsize^2); %make A give drainage area in m2 instead of # of pixels
%putting data in a struct format
[lat,long,ixc2,elev,dist, area, smoothx, smoothy] = STREAMobj2XY(S,ixc,cDEM,S.distance, A, x, y); %can make a P out of this so we don't do it twice I did this one to incorporate curvatre
slope = zeros(size(dist));
slope(:) = NaN;

for i = 2:(length(dist)-1)
    slope(i) = (elev(i-1)-elev(i+1))./(dist(i-1)-dist(i+1));
end
%%
% Separate the different tributaries
NaN_locations = find(isnan(elev)); % these are the boundaries between the tributaries

%make sure to set the slope at the tributary and one cell above and below 
%to NaN, otherwise, you might accidentally associate these with waterfalls.
slope(NaN_locations) = NaN;  
slope(NaN_locations(1:end-1)+1) = NaN;
slope(NaN_locations-1) = NaN;

number_tribs = length(NaN_locations); 
NaN_locations = [1; NaN_locations];

trib_id = zeros(length(dist),1);
for i = 1:(number_tribs)
    trib_id(NaN_locations(i):NaN_locations(i+1))=i;
end

trib_id(NaN_locations) = NaN;
for jjjj=1:5
    trib_id(NaN_locations(1:end-1)+jjjj) = NaN;
    trib_id(NaN_locations(2:end)-jjjj) = NaN;
end
allnan=find(isnan(trib_id));
%% making a struct with channel data
P.x =lat; %x coordinate
P.smoothx=smoothx;
P.smoothy=smoothy;
P.y = long;
%P.z=nan(size(lat));
P.z= elev; %elevation
P.trib_id=trib_id;
P.ntribs=number_tribs;
P.d = dist; %distance along channel
P.s = nan(size(P.z));
P.s_deg = nan(size(P.z));
P.a_km = area./1e6; %drainage area along profile in km2
P.a_m = area; %drainage area along profile in m2

P.s(1:end) = slope; %slope (m/m)
P.s_deg(1:end) = atand(slope); %convert slope to degrees
avgby=250;

P_savg=zeros(size(P.d));  %constructing avg slope
for i =avgby+1:(size(P.d)-avgby)
    P_savg(i)=mean(P.s((i-avgby):(i+avgby)), 1, "omitnan");
end
%%
%making a grid - to build future fake DEMs off of
%note the if statements are because sometimes the spatial reference is
%documented slightly differently and it makes the columns and rows 1 short
tuo_grid.grid=cDEM.Z;
tuo_grid.ncols=length(cDEM.Z(1,:));
tuo_grid.nrows=length(cDEM.Z(:,1));
tuo_grid.x=zeros(1, tuo_grid.ncols);
if length(cDEM.georef.SpatialRef.XWorldLimits(1):cDEM.georef.SpatialRef.XWorldLimits(2))==tuo_grid.ncols

    tuo_grid.x(1,:)=cDEM.georef.SpatialRef.XWorldLimits(1):cDEM.georef.SpatialRef.XWorldLimits(2);
end
if length(cDEM.georef.SpatialRef.XWorldLimits(1):cDEM.georef.SpatialRef.XWorldLimits(2))~=tuo_grid.ncols

    tuo_grid.x(1,:)=cDEM.georef.SpatialRef.XWorldLimits(1):(cDEM.georef.SpatialRef.XWorldLimits(2)-1);
end
tuo_grid.y=zeros(tuo_grid.nrows, 1);
if length(cDEM.georef.SpatialRef.YWorldLimits(1):cDEM.georef.SpatialRef.YWorldLimits(2))==tuo_grid.nrows
    tuo_grid.y(:,1)=cDEM.georef.SpatialRef.YWorldLimits(2):-1:cDEM.georef.SpatialRef.YWorldLimits(1);
end
if length(cDEM.georef.SpatialRef.YWorldLimits(1):cDEM.georef.SpatialRef.YWorldLimits(2))~=tuo_grid.nrows
    tuo_grid.y(:,1)=(cDEM.georef.SpatialRef.YWorldLimits(2)-1):-1:(cDEM.georef.SpatialRef.YWorldLimits(1));
end
%saveas(gcf, strcat('E:/san_gab_tifs/resultsfigs/wf_', fname ,'.png'))

%%
%S=trunk(S);
if ~isnan(long_bot)
    [xds,yds] = snap2stream(S,long_bot,lat_bot);
    [xus,yus] = snap2stream(S,long_top,lat_top);
    mds=find(P.x==xds);
    nds=find(P.y==yds);
    [ibot,pos1]=intersect(mds,nds);
    mus=find(P.x==xus);
    nus=find(P.y==yus);
    [itop,pos2]=intersect(mus,nus);
end

if isnan(long_bot)
    itop=n;
    ibot=length(P.z)-n;
end

%NEED TO FIDDLE WITH ANGLES - currently just using the linear relationship,
%so could take out all the polynomial stuff



upperbound=itop+n; %% can lower down to 0 if you want - sets how far from the top you stop calculating
lowerbound=ibot-n;%Same but for the bottom - note the distance parameter is zero at the headwaters

riversmooth=zeros(length((P.x)), 3);
riversmooth_slope=zeros(length((P.x)), 1);
riverperp_poly=zeros(length((P.x)), 1);

riversmo_lin=zeros(length((P.x)), 2);
for i=upperbound:(lowerbound) %this started as i=n+1:(length(P.x)-n), but to look at a smaller section I'll do
    riversmooth(i,:)=polyfit(P.smoothx(i-n:i+n), P.smoothy(i-n:i+n), 2);
    riversmo_lin(i,:)=polyfit(P.smoothx(i-n:i+n), P.smoothy(i-n:i+n), 1);
    riversmooth_slope(i,:)=(2*riversmooth(i,1)*P.smoothx(i))+riversmooth(i,2);
    riverperp_poly(i)=-1/riversmooth_slope(i);
    
end
%riverperp=zeros((length(P.x)), 2);

%riverperp contains the information needed to draw a line in the DEM the
%first column has the slope- m, while the second is the intercept b= ( (-m*
%x1) +y1)
%riverperp=-1./(2*riversmooth(:,1)+riversmooth(:,2)); %this is the perpendicular line to the river
riverperp=-1./(riversmo_lin(:,1));
%angle=atan(1./riverperp_poly);
angle=atan(1./riverperp);
angle=rad2deg(angle);
stream_curve=2*riversmooth(:,1);

badangle=atan(1/riverperp);
%% making a bunch of arrays to use later - most of these I want to change into multidimensional arrays
notmessedup=0;  % this counts how many spots were skipped because the data didn't work
%now in inputs   % this controls the distance between cuts - i.e. if you want to calculate width every 1,2,5,10 m (given a 1m dataset) 
winsize=windowsize;    %how big the window is that you use to look for widths - i.e. if the width is more than winsize it will skip it
messedup=[];
num_messed=0;
fixedup=[];
num_fixed=0;
num_right=0;
num_left=0;

cuts_index=[upperbound:interval:(lowerbound)];
cuts_index=cuts_index(~ismember(cuts_index, allnan));
avgprof=nan(length(P.x),1 );
cuts_count=length(cuts_index);


cross_profs=NaN(length((P.x)),winsize*2 ); %this will record the swath profiles, with each index corresponding to a point on the elevation profile (only one every 50 is filled)
cross_profs_x=NaN(length((P.x)),winsize*2); %this will record the x coord  of swath profiles, with each row corresponding to a point on the river elevation profile
cross_profs_y=NaN(length((P.x)),winsize*2); %this will record the y coord of swath profiles, with each row corresponding to a point on the elevation profile (only one every 50 is filled)
%cross_profs_smooth=NaN(length((P.x)),winsize*2); 
min_swath_elevation=NaN(length(P.z), 1); %totest whether pulling a profile like this gets a much diff profile than the river profile
channel_elev=NaN(length(P.z),1);
channel_mid_col=NaN(length(P.z),1);
channel_mid_row=NaN(length(P.z),1);
lengthP=length(P.z);
lengthmeas=length(elev_measure);
%ALSO WANT TO BUILD WIDTHS WHILE WE"RE ZOOMED IN


allwidth=NaN(lengthP, lengthmeas);
prof_clipped_allxtot=NaN(lengthP, winsize*2, lengthmeas);
prof_clipped_allytot=NaN(lengthP, winsize*2, lengthmeas);
allsymnrm=NaN(lengthP, lengthmeas);

for i=cuts_index%cuts_index%cuts_index               %

    %i=16000; %JUST FOR TESTING
    %FIRST figure out where your profile point is in the DEM

    x1=find(int32(tuo_grid.x)==int32(P.x(i)));
    y1=find(int32(tuo_grid.y)==int32(P.y(i)));
    z1=tuo_grid.grid(y1, x1);
    %we check that we have the data to do so, and then make a smaller DEM
    %centered around that profile point
    %dist_array=GRIDobj(tuo_grid.x(x1-(winsize-1):x1+winsize), tuo_grid.y(y1-(winsize-1):y1+winsize), zsmooth_tuogrid(y1-(winsize-1):y1+winsize,x1-(winsize-1):x1+winsize));
    if (x1-(winsize-1)>0) && (x1+winsize < length(tuo_grid.x)) && (y1-(winsize-1)>0) && (y1+winsize< length(tuo_grid.y))

        %pgon = polyshape([0 0 (x1-(winsize-1)) (x1+winsize)],[(y1-(winsize-1)) 0 0 (y1+winsize]);
        dist_array=GRIDobj(tuo_grid.x(x1-(winsize-1):x1+winsize), tuo_grid.y(y1-(winsize-1):y1+winsize), tuo_grid.grid(y1-(winsize-1):y1+winsize,x1-(winsize-1):x1+winsize));

        [dist_array_coord_x, dist_array_coord_y]=getcoordinates(dist_array);
        D=dist2line(dist_array, dist_array_coord_x(winsize), dist_array_coord_y(winsize),angle(i));  %shows the distance of each point in the window from the perpendicular line to the river
        mindist_coord_x=NaN(1, 2*winsize); %array to put in the coordinates of the point in each row that is closest to that line
        mindist_coord_y=NaN(1, 2*winsize);
        mindist_index_y=NaN(1, 2*winsize);
        mindist_index_x=NaN(1, 2*winsize);
        mindist_elev=NaN(2*winsize,1); %elevations along that swath
        mindist_elev_smooth=NaN(2*winsize,1); %elevations along that swath
        array_onprof=zeros(2*winsize,2*winsize);
        channel_elev_loc=NaN(2,1);
    

        
        %GO THROUGH THE WINDOW to get an elevation profile
        for j=1:(2*winsize) %this loop goes through every row to find the ppoint on the row with the least
            % distance from the and puts that point on the cross profile
            %the problem is that if the angle is between 45 and 135 degrees, it
            %will strike horizontally, and if its between 0-45 or 135-180 it
            %will strike vertically so need to go through either rows or
            %columns depending on what this angle is. only for 45 it doesn't
            %matter
            if abs(angle(i))>45
                mindist_coord_y(j)=min(dist_array_coord_y(D.Z(:,j)==min(D.Z(:,j))));
        
                mindist_coord_x(j)=dist_array_coord_x(j);
                mindist_elev(j)=min(dist_array.Z( find(D.Z(:,j)==min(D.Z(:,j))),j ));
                %mindist_elev_smooth(j)=zsmooth_tuogrid(dist_array_coord_y==mindist_coord_y(j), dist_array_coord_x==mindist_coord_x(j) );
                mindist_index_x(j)=find(dist_array_coord_x==mindist_coord_x(j));
                mindist_index_y(j)=find(dist_array_coord_y==mindist_coord_y(j));
                array_onprof(mindist_index_y(j), mindist_index_x(j))=1;
            end
            if abs(angle(i))<=45
                mindist_coord_x(j)=min(dist_array_coord_x(D.Z(j,:)==min(D.Z(j,:))));
        
                mindist_coord_y(j)=dist_array_coord_y(j);
                mindist_elev(j)=min(dist_array.Z(j, find(D.Z(j,:)==min(D.Z(j,:)) )));
                mindist_index_x(j)=find(dist_array_coord_x==mindist_coord_x(j));
                mindist_index_y(j)=find(dist_array_coord_y==mindist_coord_y(j));
                %mindist_elev_smooth(j)=zsmooth_tuogrid(dist_array_coord_y==mindist_coord_y(j), dist_array_coord_x==mindist_coord_x(j) );
                array_onprof(mindist_index_y(j), mindist_index_x(j))=1;
            end
        end
        channel_elev(i)=min(mindist_elev((winsize-2):(winsize+2)));

        % these are the coordinates of the actual lowest part of the channel for if I want to center the window for channel width around it
        channel_elev_locx=find(dist_array_coord_x==min(mindist_index_x(mindist_elev==channel_elev(i))));%This is problematic because often it captures more than one river cross section
        channel_elev_locy=find(dist_array_coord_y==min(mindist_index_y(mindist_elev==channel_elev(i))));
       % channel_mid_col(i)=dist_array_coord_x(dist_array.Z==min(mindist_elev((winsize-10):(winsize+10))))
        %FROM HERE WE CAN JUST GO THROUGH OUR ARRAY OF ELEVS TO MEASURE AT
            %NEW STUFF
        for elevn=1:length(elev_measure)
                undern=(channel_elev(i)+elev_measure(elevn)>=mindist_elev);%channel_elev(i) is the min channel elev, 
                % elev_measure(elevn) is the distance above the channel that we a re measureing and
                % mindist_elev is an array of the channel cross cut
    
                underndif=undern(1:end-1)-undern(2:end);
                undernneg=find(underndif==-1); %pulls the index within the cross section of the left edge of under x m
                undernpos=find(underndif==1);
                locnl=max(undernneg((winsize-undernneg)>0)); %makes sure that we are grabbing the closest one to the left side its on the correct side of the profile (not grabbing a different cur farther upstream or downstream 
                locnr=min(undernpos((undernpos-winsize)>0));% makes sure we are grabbing the closest to the right side
    
    
    %            %sanity check figure
    %             figure
    %             plot(mindist_elev)
    %             hold on
    %             yline(channel_elev(i)+elev_measure(elevn))
    %             scatter([locnl, locnr],[mindist_elev(locnl), mindist_elev(locnr)] )
    %             scatter([locnl-1, locnr+1],[mindist_elev(locnl-1), mindist_elev(locnr+1)] )
    %             plot(undern*(channel_elev(i)+elev_measure(elevn)))
    %             scatter([locnl-extral, locnr+extrar], [mindist_elev(locnl)+locl_elev_below, mindist_elev(locnr)+locr_elev_below])
    %             ylim([min(mindist_elev)-3, max(mindist_elev)+1])
    
    
                if (~isempty(locnl) && ~isempty(locnr))
    
                    %INTERPOLATION
                    locnl=locnl+1;
                    locl_elev_below=channel_elev(i)+elev_measure(elevn)-mindist_elev(locnl);
                    locr_elev_below=channel_elev(i)+elev_measure(elevn)-mindist_elev(locnr);
                    slopel=(mindist_elev(locnl-1)-mindist_elev(locnl))/sqrt((mindist_coord_x(locnl-1)-mindist_coord_x(locnl))^2+(mindist_coord_y(locnl-1)-mindist_coord_y(locnl))^2);
        
                    sloper=(mindist_elev(locnr+1)-mindist_elev(locnr))/sqrt((mindist_coord_x(locnr+1)-mindist_coord_x(locnr))^2+(mindist_coord_y(locnr+1)-mindist_coord_y(locnr))^2);
                    extral=locl_elev_below/slopel;
                    extrar=locr_elev_below/sloper;
    
    
                    prof_clipped=mindist_elev(locnl:locnr);
                    prof_clipped_allxtot(i, 1:length(locnl:locnr), elevn)=mindist_coord_x(locnl:locnr);
                    prof_clipped_allytot(i, 1:length(locnl:locnr), elevn)=mindist_coord_y(locnl:locnr);
                    allwidth(i, elevn)=sqrt((mindist_coord_x(locnl)-mindist_coord_x(locnr))^2+(mindist_coord_y(locnl)-mindist_coord_y(locnr))^2)+extral+extrar;
                    centernl=sqrt((mindist_coord_x(locnl)-x1)^2 + (mindist_coord_y(locnl)-y1)^2)+extral;
                    centernr=sqrt((mindist_coord_x(locnr)-x1)^2 + (mindist_coord_y(locnr)-y1)^2)+extrar;
                    allsymnrm(i, elevn)=(centernl-centernr)./allwidth(i);
                end
    
    
           
    
                if (isempty(locnl) || isempty(locnr))
                    stopit=0;
                    num_messed=num_messed+1;
                    messedup(num_messed)=(i);
                    for error=1:4
                        if (isempty(locnl) & ~isempty(locnr) ) 
                            %maximum index that is a left boundary, within
                            %error of where the thalweg should be = 
                            locnl=max(undernneg(((winsize+error)-undernneg)>=0)); %the main error is that the channel is all to the left or right of point winsize
                        end 
                        if (~isempty(locnl) & isempty(locnr))
                            locnr=min(undernpos((undernpos-(winsize-error))>=0)); %the main error is that the channel is all to the left or right of point winsize
                        end
                       
    
                        if (~isempty(locnl) && ~isempty(locnr) & stopit==0)
                            %+1 to
                            %correct for marking error
                            locnl=locnl+1;
                            %INTERPOLATION
                            locl_elev_below=channel_elev(i)+elev_measure(elevn)-mindist_elev(locnl);
                            locr_elev_below=channel_elev(i)+elev_measure(elevn)-mindist_elev(locnr);
                            slopel=(mindist_elev(locnl-1)-mindist_elev(locnl))/sqrt((mindist_coord_x(locnl-1)-mindist_coord_x(locnl))^2+(mindist_coord_y(locnl-1)-mindist_coord_y(locnl))^2);
                
                            sloper=(mindist_elev(locnr+1)-mindist_elev(locnr))/sqrt((mindist_coord_x(locnr+1)-mindist_coord_x(locnr))^2+(mindist_coord_y(locnr+1)-mindist_coord_y(locnr))^2);
                            extral=locl_elev_below/slopel;
                            extrar=locr_elev_below/sloper;
    
                            allwidth(i, elevn)=sqrt((mindist_coord_x(locnl)-mindist_coord_x(locnr))^2+(mindist_coord_y(locnl)-mindist_coord_y(locnr))^2)+extral+extrar;
                            prof_clipped=mindist_elev(locnl:locnr);
        
                            centernl=sqrt((mindist_coord_x(locnl)-x1)^2 + (mindist_coord_y(locnl)-y1)^2)+extral;
                            centernr=sqrt((mindist_coord_x(locnr)-x1)^2 + (mindist_coord_y(locnr)-y1)^2)+extrar;
                            allsymnrm(i, elevn)=(centernl-centernr)./allwidth(i);
                            prof_clipped_allxtot(i, 1:length(locnl:locnr), elevn)=mindist_coord_x(locnl:locnr);
                            prof_clipped_allxtot(i, 1:length(locnl:locnr), elevn)=mindist_coord_y(locnl:locnr);
    %                                   %sanity check figure
    %                         figure
    %                         plot(mindist_elev)
    %                         hold on
    %                         yline(channel_elev(i)+elev_measure(elevn))
    %                         scatter([locnl, locnr],[mindist_elev(locnl), mindist_elev(locnr)] )
    %                         scatter([locnl-1, locnr+1],[mindist_elev(locnl-1), mindist_elev(locnr+1)] )
    %                         plot(undern*(channel_elev(i)+elev_measure(elevn)))
    %                         scatter([locnl-extral, locnr+extrar], [mindist_elev(locnl)+locl_elev_below, mindist_elev(locnr)+locr_elev_below])
    %                         ylim([min(mindist_elev)-3, max(mindist_elev)+1])
                        end
                        if (~isempty(locnl) && ~isempty(locnr))
                            stopit=1;
                        end
                    end
                end
        end%end if the data is in range loop
    end%end cuts_index loop
end

%% adding waterfalls
        % going trib by trib for assessing wf freq and average slope/ksn
profile_plots=0;
for i = 1:P.ntribs

% Create a temporary structure that's just for the tributary
    i1 = P.trib_id == i;
    P_temp = P;
    P_temp.x = P.x(i1);
    P_temp.y = P.y(i1);
    P_temp.z = P.z(i1);
    P_temp.a_km = P.a_km(i1);
    P_temp.a_m = P.a_m(i1);
    P_temp.d = P.d(i1);
    P_temp.s = P.s(i1);
    P_temp.s_deg = P.s_deg(i1);
    
    % make sure there are no NaN values in the temporary structure
    i2 = isnan(P_temp.d);
    P_temp.x = P_temp.x(~i2);
    P_temp.y = P_temp.y(~i2);
    P_temp.z = P_temp.z(~i2);
    P_temp.a_km = P_temp.a_km(~i2);
    P_temp.a_m = P_temp.a_m(~i2);
    P_temp.d = P_temp.d(~i2);
    P_temp.s = P_temp.s(~i2);
    P_temp.s_deg = P_temp.s_deg(~i2);
    fname=string(k);
    % Run the WF code
    eval(['[WF' num2str(i) '] = wf_finder(P_temp,thresh_slope,thresh_dh,fname,profile_plots);'])

    % Concatenate all the waterfall metrics from each tributary into a
    % single vector.
    if i>1
        eval(['WF_concatenated.count_wfs = [WF_concatenated.count_wfs; WF' num2str(i) '.count_wfs];'])
        eval(['WF_concatenated.UTME_wf_top = [WF_concatenated.UTME_wf_top; WF' num2str(i) '.UTME_wf_top];'])
        eval(['WF_concatenated.UTMN_wf_top = [WF_concatenated.UTMN_wf_top; WF' num2str(i) '.UTMN_wf_top];'])
        eval(['WF_concatenated.dist_bot = [WF_concatenated.dist_bot; WF' num2str(i) '.dist_bot];'])
        eval(['WF_concatenated.dist_top = [WF_concatenated.dist_top; WF' num2str(i) '.dist_top];'])
        eval(['WF_concatenated.z_bot = [WF_concatenated.z_bot; WF' num2str(i) '.z_bot];'])
        eval(['WF_concatenated.z_top = [WF_concatenated.z_top; WF' num2str(i) '.z_top];'])
        eval(['WF_concatenated.A = [WF_concatenated.A; WF' num2str(i) '.A];'])
        eval(['WF_concatenated.Lu = [WF_concatenated.Lu; WF' num2str(i) '.Lu];'])
        eval(['WF_concatenated.Hwf = [WF_concatenated.Hwf; WF' num2str(i) '.Hwf];'])
        eval(['WF_concatenated.Hwf2Lu = [WF_concatenated.Hwf2Lu; WF' num2str(i) '.Hwf2Lu];'])
        eval(['WF_concatenated.Su = [WF_concatenated.Su; WF' num2str(i) '.Su];'])
        eval(['WF_concatenated.Sr = [WF_concatenated.Sr; WF' num2str(i) '.Sr];'])
    elseif i==1
        WF_concatenated = WF1;
    end
end

% 
% figure
% plot(P.d,P.z,'k-','LineWidth',1,'DisplayName','profile')
% hold on
% plot(WF_concatenated.dist_top,WF_concatenated.z_top,'ko','markerfacecolor','w','DisplayName','Waterfall Lip','markersize',6)
% plot(WF_concatenated.dist_bot,WF_concatenated.z_bot, 'ko','markerfacecolor','r','DisplayName','Waterfall Base','markersize',6)
% title(fname)
% xlabel('Distance from Outlet (m)')
% ylabel('Elevation (m)')
% legend
% 
% saveas(gcf, strcat('F:/sangabs_lidar/v6_nochutes_ksn/' ,fname ,'_wfprof.png'))


%% making all the averages
P_savg=zeros(size(P.d));  %constructing avg slope
width_avg=zeros(size(P.d));
wf_freq=zeros(size(P.d));
wf_top=zeros(size(P.d));
wf_bot=zeros(size(P.d));
dist_abv_in_wf=zeros(size(P.d));
dist_blw_wf=zeros(size(P.d));
dist_in_wf=zeros(size(P.d));
updist=400;
downdist=400;

for jj = 1:P.ntribs
    WF=eval(['WF' num2str(jj)]);
    %the k= are for removing waterfalls in certain basins
    if k==185
                if length(WF.dist_bot(WF.dist_bot<2200 & WF.dist_bot>1800))>0
                    indexbad=find(WF.dist_bot<2200 & WF.dist_bot>1800)
                        WF.dist_bot(indexbad)=[];
                        k
                        jj
                end
                if length(WF.dist_top(WF.dist_top<2200 & WF.dist_top>1800))>0
                    indexbad=find(WF.dist_top<2200 & WF.dist_top>1800)
                        WF.dist_top(indexbad)=[];
                        k
                end 
      
    end

    if k==212
                if length(WF.dist_bot(WF.dist_bot<975))>0
                    indexbad=find(WF.dist_bot<975)
                        WF.dist_bot(indexbad)=[];
                        k
                        jj
                end
                if length(WF.dist_top(WF.dist_top<955))>0
                    indexbad=find(WF.dist_top<955)
                        WF.dist_top(indexbad)=[];
                        k
                end 
      
    end
     if k==159
                if length(WF.dist_bot(WF.dist_bot<584))>0
                    indexbad=find(WF.dist_bot<584)
                        WF.dist_bot(indexbad)=[];
                        k
                        jj
                end
                if length(WF.dist_top(WF.dist_top<584))>0
                    indexbad=find(WF.dist_top<584)
                        WF.dist_top(indexbad)=[];
                        k
                end 
      
     end


%       if k==210
%                 if length(WF.dist_bot(WF.dist_bot<584))>0
%                     indexbad=find(WF.dist_bot<584)
%                         WF.dist_bot(indexbad)=[];
%                         k
%                         jj
%                 end
%                 if length(WF.dist_top(WF.dist_top<584))>0
%                     indexbad=find(WF.dist_top<584)
%                         WF.dist_top(indexbad)=[];
%                         k
%                 end 
%       
%      end

      if k==220
                if length(WF.dist_bot(WF.dist_bot<1500))>0
                    indexbad=find(WF.dist_bot<1500)
                        WF.dist_bot(indexbad)=[];
                        k
                        jj
                end
                if length(WF.dist_top(WF.dist_top<1500))>0
                    indexbad=find(WF.dist_top<1500)
                        WF.dist_top(indexbad)=[];
                        k
                end 
      
    end
    if k==25
                if length(WF.dist_bot(WF.dist_bot<1150))>0
                    indexbad=find(WF.dist_bot<1505)
                        WF.dist_bot(indexbad)=[];
                        k
                        jj
                end
                if length(WF.dist_top(WF.dist_top<1150))>0
                    indexbad=find(WF.dist_top<1505)
                        WF.dist_top(indexbad)=[];
                        k
                end 
      
    end






    if~isempty(WF.dist_bot)
        if WF.dist_bot(1)>WF.dist_top(1) %some of them seem to have an erroneous wf bot on top
            WF.dist_bot(1:end-1)=WF.dist_bot(2:end);
            WF.dist_bot(end)=NaN;
            WF.dist_top(end)=NaN;
        end
    end
    for i =min(find(P.trib_id==jj))+avgby:((max(find(P.trib_id==jj)))-avgby)
        wf_freq(i)=length(WF.dist_bot(((WF.dist_bot<P.d((i-avgby)) & WF.dist_bot>=P.d(i+avgby)))));%.*400/(P.d((i-avgby))-P.d(i+avgby)); %length(WF.dist_bot((WF.dist_bot>(i-n) & WF.dist_bot<=(i+n))));
        width_avg(i)=mean(allwidth((i-avgby):(i+avgby), 1), 1, 'omitnan');
        P_savg(i)=mean(P.s((i-avgby):(i+avgby)), 1, "omitnan");

        if ismember(P.d(i), WF.dist_top)
            wf_top(i)=1;%length(WF.dist_top(((WF.dist_top<P.d((i-1)) & WF.dist_top>=P.d(i+1)))));%.*400/(P.d((i-avgby))-P.d(i+avgby)); %length(WF.dist_bot((WF.dist_bot>(i-n) & WF.dist_bot<=(i+n))));
    
        end

        if ismember(P.d(i), WF.dist_bot)
            wf_bot(i)=1;%length(WF.dist_bot(((WF.dist_bot<P.d((i-1)) & WF.dist_bot>=P.d(i+1)))));%.*400/(P.d((i-avgby))-P.d(i+avgby)); %length(WF.dist_bot((WF.dist_bot>(i-n) & WF.dist_bot<=(i+n))));
        end
    end
    dist_abv_in_wf_trib=zeros(size(P.d(P.trib_id==jj)));
    dist_in_wf_trib=zeros(size(P.d(P.trib_id==jj)));

    dist_blw_wf_trib=zeros(size(P.d(P.trib_id==jj)));
    for i=1:length(WF.dist_bot(~isnan(WF.dist_bot)))
        % we set all the points within the waterfall to 1

        inwf=find(P.d(P.trib_id==jj)<=WF.dist_top(i) & P.d(P.trib_id==jj)>=WF.dist_bot(i));%indices of waterfall sites
        dist_in_wf_trib(inwf)=ones(length(inwf), 1);

        if i>1 %if there isn't a waterfall within updist m upstream, label those points
            if (find(P.d(P.trib_id==jj)==WF.dist_bot(i))-find(P.d(P.trib_id==jj)==WF.dist_bot(i-1)))>=updist
                dist_abv_in_wf_trib( (max(inwf)-updist+1): (max(inwf)))=flip(linspace(1, updist, updist)*-1);
            end
            %if there is, only label that many
            if (find(P.d(P.trib_id==jj)==WF.dist_bot(i))-find(P.d(P.trib_id==jj)==WF.dist_bot(i-1)))<updist
                spacing=(find(P.d(P.trib_id==jj)==WF.dist_bot(i))-find(P.d(P.trib_id==jj)==WF.dist_bot(i-1)));
                dist_abv_in_wf_trib( (max(inwf)-spacing+1):(max(inwf)))=flip(linspace(1, spacing, spacing)*-1);

            end
        end
        if i==1 %special case for if its the top waterfall
            if (find(P.d(P.trib_id==jj)==WF.dist_top(i))-find(P.d(P.trib_id==jj)==max(P.d(P.trib_id==jj))))>=updist
                dist_abv_in_wf_trib(max(inwf)-updist+1: max(inwf))=flip(linspace(1, updist, updist)*-1);
            end
            if (find(P.d(P.trib_id==jj)==WF.dist_top(i))-find(P.d(P.trib_id==jj)==max(P.d(P.trib_id==jj))))<updist
                spacing=(find(P.d(P.trib_id==jj)==WF.dist_top(i))-find(P.d(P.trib_id==jj)==max(P.d(P.trib_id==jj))));

                dist_abv_in_wf_trib((max(inwf)-spacing+1): max(inwf))=flip(linspace(1, spacing, spacing)*-1);
            end
        end
        dist_abv_in_wf(P.trib_id==jj)=dist_abv_in_wf_trib;
        dist_in_wf(P.trib_id==jj)=dist_in_wf_trib;

        %now moving on to making this for below waterfalls
        if i==length(WF.dist_bot)%the last waterfall will need to check the bottom of the array
            if (find(P.d(P.trib_id==jj)==min(P.d(P.trib_id==jj)))-find(P.d(P.trib_id==jj)==WF.dist_bot(i)))>=downdist
                dist_blw_wf_trib( (max(inwf)+1): (max(inwf)+downdist))=linspace(1, downdist, downdist);
                
            end        
            if (find(P.d(P.trib_id==jj)==min(P.d(P.trib_id==jj)))-find(P.d(P.trib_id==jj)==WF.dist_bot(i)))<downdist
                spacing=(find(P.d(P.trib_id==jj)==min(P.d(P.trib_id==jj)))-find(P.d(P.trib_id==jj)==WF.dist_bot(i)));
                dist_blw_wf_trib( (max(inwf)+1): (max(inwf)+spacing))=linspace(1, spacing, spacing);
                
            end  
    

        end
        if i<length(WF.dist_top)
            if (find(P.d(P.trib_id==jj)==WF.dist_top(i+1))-find(P.d(P.trib_id==jj)==WF.dist_top(i)))>=downdist
                dist_blw_wf_trib( (max(inwf)+1): (max(inwf)+downdist))=linspace(1, downdist, downdist);
                
            end

            if (find(P.d(P.trib_id==jj)==WF.dist_top(i+1))-find(P.d(P.trib_id==jj)==WF.dist_top(i)))<downdist
                    spacing=(find(P.d(P.trib_id==jj)==WF.dist_top(i+1))-find(P.d(P.trib_id==jj)==WF.dist_top(i)))+4;
                if length(dist_blw_wf_trib)>=(max(inwf)+spacing)
                    dist_blw_wf_trib( (max(inwf)+1): (max(inwf)+spacing))=linspace(1, spacing, spacing);
                end
                if length(dist_blw_wf_trib)<(max(inwf)+spacing)
                    spacing=(find(P.d(P.trib_id==jj)==WF.dist_top(i+1))-find(P.d(P.trib_id==jj)==WF.dist_top(i)));
                    dist_blw_wf_trib( (max(inwf)+1): (max(inwf)+spacing))=linspace(1, spacing, spacing);
                end
            end 

        end
        dist_blw_wf(P.trib_id==jj)=dist_blw_wf_trib;
    end

end

% figure
% plot(P.d, wf_top)
% hold on
% plot(P.d, wf_bot)
% plot(P.d, wf_freq)
% 
% plot(P.d, dist_abv_in_wf)
% plot(P.d,  dist_blw_wf)
% xlabel('distance (m)')
% legend('top node', 'bottom node','30m window', 'above and in wf', 'below wf')




%%
equations=NaN(length(P.z), 4);
% for j =cuts_index
%     widthsforcut=allwidth(j, :);
%     isnotnan=~isnan(widthsforcut);
%         fitobj=fit(transpose(elev_measure(isnotnan)), transpose(widthsforcut(isnotnan)), 'power1');
%         fitwidths=coeffvalues(fitobj);
%         equations(j, 1)=fitwidths(1);
%         equations(j,2)=fitwidths(2);
% end
picplot=0;
count=1;
for i=length(elev_measure):-1:1
    
    elev_reord(count)=elev_measure(i);
    count=count+1;
end
%makes a plot of the width calc on the river- computationally very
%expensive 
if picplot==1
    figure
    imageschs(cDEM, cDEM, 'colorbarylabel', 'Elevation (m)')
    hold on
    cmap=parula(16);
    smallguy=0;
    lengthofcut=nan(length(P.d), length(elev_measure));
    for l=cuts_index
        for iiii=length(elev_measure):-1:1
            lengthofcut(l, iiii)=length(find(~isnan(prof_clipped_allxtot(l,:, iiii))));
            plot(prof_clipped_allxtot(l,1:lengthofcut(l, iiii), iiii), prof_clipped_allytot(l,1:lengthofcut(l, iiii), iiii),color=cmap(iiii*2, :) , LineWidth=3)
            hold on
        end
    end 
    scatter(WF_concatenated.UTME_wf_top, WF_concatenated.UTMN_wf_top, 'o', 'red')
    
    
    str_elev_reord=string(elev_reord);
    legend(str_elev_reord, 'Location', 'eastoutside')
    %scatter(WF.UTME_wf_top, WF.UTMN_wf_top, 'r')
    hold on
    xlabel('UTM Easting (m)')
    %title(fname)
    ylabel('UTM Northing (m)')
    %saveas(gcf, strcat('E:/san_gab_tifs/resultsfigs/hs_', fname ,'.png'))
    hold off
end
%% results analysis
%sum(isnan(allwidth(:, 2)))
%thirdnan=find(isnan(allwidth(:, 1)));
%halfnan=find(isnan(allwidth(:, 2)));
%third_messup=thirdnan(find(~ismember(thirdnan, halfnan)));
% for i=cuts_index(1:10)
%     fitobj=fit(transpose(elev_measure(~isnan(allwidth(i, 1:4)))), transpose(allwidth(i, 1:4)), 'power1')
%     fitobj2=fit(transpose(elev_measure( 5:8)), transpose(allwidth(i, 5:8)), 'poly3')
%     figure
%     scatter(elev_measure, allwidth(i, :))
%     hold on
%     %plot(elev_measure, equations(i, 1)*elev_measure.^3+ equations(i, 2).*elev_measure.^2 +equations(i, 3).*elev_measure+equations(i, 4));
%     plot(fitobj)
%     hold on
%     plot(fitobj2)
% end
% elevcat=categorical(elev_measure);
% legend(elevcat)
% xlabel('elev above bottom (m)')
% ylabel('width (m)')
end%end function