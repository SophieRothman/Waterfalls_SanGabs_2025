%This script calcualtes flow width. depth, and velocity given a cross
%sectional form estimated by 6 measurements at different elevations above
%the thalweg
%We use Equation 20 from Ferguson, (2007), where we set a1=6.5, a2=2.5
%Code written by Sophie Rothman
%Part of the supplement for "Waterfalls perturb channel and hillslope 
% dynamics in the San Gabriel Mountains, California" by Rothman,
% Scheingross and McCoy

%%   Required inputs
%elev_measure=a 1xm array of the m elevations above the thalweg that were
            %measured
%allwidths_bf=a Lxm array with the widths calculated at each elevation in
%           elev_array, across every cross section
%alldas=Lx1 array of drainage area at each cross section
%allslope_avg=Lx1 array of slope at each cross section, calculated over a 15 m
%       rolling average


R = 1.65; %normalized submerged sediment density (non-dimensional)
g = 9.81; %gravitational acceleration
a1=6.5;   %%empirical constants
a2=2.5;    %%empirical constants
detail=500;
Hrange=linspace(.01, 2.5, detail);  %Array of all the flow heights we're inputting to see to what Q it corresponds
Vrange=NaN(1, detail);  %array of the resultant velocities for the flow heights
Qrange=NaN(1, detail);  %Array of the resultant discharges for the flow hegihts
Wrange=NaN(length(dist), detail);



discharge=[33.4]; %SET DISCHARGE TO WHATEVER
elevmeasure=[.1, .25, .5, 1, 1.5, 2.5];

%Create arrays for the results
h_result=NaN(length(alldas), length(discharge)); 
v_result=NaN(length(alldas), length(discharge)); 
w_result=NaN(length(alldas), length(discharge)); 
q_result=NaN(length(alldas), length(discharge)); 

h_result=zeros((length(allwidths_bf(:,1))), length(discharge));
w_result=zeros((length(allwidths_bf(:,1))), length(discharge));
q_result=zeros((length(allwidths_bf(:,1))), length(discharge));
v_result=zeros((length(allwidths_bf(:,1))), length(discharge));

for i=1: length(discharge) %currently discharge is length 1 so loop isn't really needed

    for j=1:(length(allwidths_bf(:,1)))
        discharge=33.4*(alldas(j)/42); %adjust discharge to location discharge=33.4*(median_da./42); %adjust discharge to location

        D85=(-0.0076755*alldas(j)) + 0.23381;%In this case I used D50 data but the variable is called D85
        isnotnan=~isnan(allwidths_bf(j, :));
        

        if length(allwidths_bf(j, isnotnan))>3
            %We fit a linear interpolation scheme to these measurements to
            %asses widths between measurement points
            [fitobj, goodfit, o]=fit(transpose(elevmeasure( isnotnan)), transpose(allwidths_bf(j, isnotnan)), 'linearinterp');
            Wrange(j, :)=transpose(fitobj(Hrange));
    
            %This equation tells us what velocity/discharge a given flow height and
            %roughness has, so we just go through a range of flow depths to see what
            %produces the closest match to the discharge we want
            %while Q<q_baseline
                for m=1:detail
                     ushear=sqrt(g*Hrange(m)*allslope_avg(j));
                     Vrange(m)=ushear*(a1*a2*(Hrange(m)/D85))/(sqrt((a1^2)+((a2^2)*((Hrange(m)/D85)^(5/3)))));
                     Qrange(m)=Vrange(m)*Hrange(m)*Wrange(j, m);
                     if Qrange(m)<discharge(i)   %we want to know the H when Q=Qgage
                        H=Hrange(m); %Best estimate of flow depth
                        V=Vrange(m); %Best estimate of velocity
                        Q=Qrange(m); % Q 
                        W=Wrange(j, m);
                  
                     end
                end
            %end
            h_result(j, i)=H;
            v_result(j, i)=V;
            w_result(j, i)=W;
            q_result(j, i)=Q;
            
        end
        if length(allwidths_bf(j, isnotnan))<=3 %if not enough measurements were succesfull result is NaN
            h_result(j, i)=NaN;
            v_result(j, i)=NaN;
            w_result(j, i)=NaN;
            q_result(j, i)=NaN;

        end
    end
end