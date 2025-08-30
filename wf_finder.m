function [WF] = wf_finder(P,thresh_slope,thresh_dh,name,profile_plots)
%function that finds waterfalls and their metrics
%Written by Erika Groh, if used, site Groh and Scheingross (2022)
%% Meet the first condition (slope > threshold value):
thresh_slope_i = P.s_deg > thresh_slope; %this creates an index with values 0 = false and 1 = true if slope is over threshold for each location along the profile

diff_slope_i = diff(thresh_slope_i); %this finds the difference of each value in the index minus the value before it which means there will be a 1 where the index 
%goes from 0 to 1 and a -1 where the index goes from 1 to 0 so
%groups of 1 values in the slope index will be bounded by a 1 and -1 (1 at
%the top of the wf and -1 at the bottom of the wf) 

diff_slope_i = [NaN; diff_slope_i];

%% define the top and bottom of waterfalls as wherever there are multiple points together over the threshold slope
i_bottom = diff_slope_i == -1;
i_top = diff_slope_i == 1; %this selects the values in the index (of diffence between slopes over threshold)that are either 1 or -1

z_bottom = P.z(i_bottom);
z_top = P.z(i_top); %this assigns the locations of the 1 and -1 values to the actual elevation on the profile

dist_bottom = P.d(i_bottom);
dist_top = P.d(i_top); %this assigns the locations of the 1 and -1 values to the distance from outlet on the profile

%% define waterfall height
if length(z_bottom)==length(z_top)
    dh = abs((z_bottom) - (z_top)); %this takes the elevation point at the top of the waterfall minus the bottom to get the change in height between them (i.e. height of waterfall)
end
if length(z_bottom)>length(z_top)
    dh=abs((z_bottom(1:end-1)) - (z_top));
    dist_bottom=dist_bottom(1:end-1);
    z_bottom=z_bottom(1:end-1);
    i_bottom=i_bottom(1:end-1);
end
if length(z_bottom)<length(z_top)
    dh=abs((z_bottom) - (z_top(2:end)));
    dist_top=dist_top(2:end);
    z_top=z_top(2:end);
    i_top=i_top(2:end);
end
%% Meet the second condition (drop height > threshold value):
i_dhover = dh > thresh_dh;
wf_dh = dh(i_dhover); %this is saying there is a waterfall by given definition where the drop height is greater than the threshold value input

%% Find the top (lip) and bottom (base) of the waterfall

wf_dist_bot = dist_bottom(i_dhover);
wf_dist_top = dist_top(i_dhover); 

wf_z_bot = z_bottom(i_dhover);
wf_z_top = z_top(i_dhover);

%% Get the drainage area at each waterfall
DA_top = P.a_m(i_top);
Awf = DA_top(i_dhover); %drainage area at each wf top in m^2

%% Calculate number of waterfalls on each channel
count_wfs = length(wf_z_top);

%% find the UTM locations of the waterfall tops
x_top = P.x(i_top);
y_top = P.y(i_top);

wf_x_top = x_top(dh>thresh_dh); %x is UTM East
wf_y_top = y_top(dh>thresh_dh); %y is UTM North


%% Calculate the metrics of waterfall morphology

% set up some empty vectors to fill
Su = zeros(size(wf_dist_bot));
Su(:) = NaN;
Lu = Su;
Hwf = Su;
Sr = Su;

for i=2:length(wf_dist_bot)
    Su(i) = (wf_z_bot(i-1) - wf_z_bot(i)) ./ (wf_dist_bot(i-1) - wf_dist_bot(i)); %Slope of the waterfall unit: calculated from the base of a waterfall to the base of the next downstream waterfall
    Sr(i) = (wf_z_bot(i-1) - wf_z_top(i)) ./ (wf_dist_bot(i-1) - wf_dist_top(i)); %River slope (excluding waterfalls): calculated from the base of a waterfall to the lip of the next downstream waterfall
    Lu(i) = wf_dist_bot(i-1) - wf_dist_bot(i); %Length of the waterfall unit: Includes the waterfall and section above it to the base of the upstream waterfall
    Hwf(i) = wf_dh(i); %waterfall height (excludes the upstream-most waterfall)
end

Hwf2Lu = Hwf ./ Lu;

%%
%Combine the variables into one output
WF.count_wfs = count_wfs;
WF.UTME_wf_top = wf_x_top;
WF.UTMN_wf_top = wf_y_top;
WF.dist_bot = wf_dist_bot;
WF.dist_top = wf_dist_top;
WF.z_bot = wf_z_bot;
WF.z_top = wf_z_top;
WF.A = Awf; %m2
WF.Lu = Lu;
WF.Hwf = Hwf;
WF.Hwf2Lu = Hwf2Lu;
WF.Su = Su;
WF.Sr = Sr;
WF.dhover=i_dhover;

%% Plot the profiles with the waterfall lip and base
if profile_plots == 1
    figure
    plot(P.d,P.z,'k-','LineWidth',1,'DisplayName','profile')
    hold on
    plot(wf_dist_top,wf_z_top,'ko','markerfacecolor','w','DisplayName','Waterfall Lip','markersize',6)
    plot(wf_dist_bot,wf_z_bot, 'ko','markerfacecolor','r','DisplayName','Waterfall Base','markersize',6)
    title(name)
    xlabel('Distance from Outlet (m)')
    ylabel('Elevation (m)')
    legend
end