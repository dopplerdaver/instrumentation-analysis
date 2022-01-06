% Filename: SKEWT_F_DFDT_LWC_SNOWIE.m
%
% Purpose: plots T (temp_cor_c) and Td (temp_frost_c) 
%          from SLWsonde in skew-t format
%
% Usage: run 'SKEWT_F_DFDT_LWC_SNOWIE' after running 'loadplot_slwsonde'
%        uses pres, temp_cor_c and humidity from SLWsonde input from loadplot_slwsonde
%
% original code from Mike King's SKEWT.m
%


%clear; clc

%% Just change this to the name of the CSV file you want to import.
%% Make sure that Matlab is running from the same folder that the
%% CSV file is in, or it will not find it.
%FILENAME = 'ma001_20140301.csv' ;
%%FILENAME = 'ma001_20140308.csv' ;
%
%FID = fopen(FILENAME) ;
%
%STRING = [ '%s %s %f %f %f %f %f %f %f %f %f %f %f %f ' ...
%           '%f %f %f %f %f %f %f %f %f %f %f %f %f %f ' ...
%           '%f %f %f %s %f %f %f %f %f %f %f %f %*[^\n]' ] ;
%
%SDATA = textscan(FID,STRING,10000,'headerlines',17,'delimiter',',') ;
%
%fclose(FID) ;
%
%NDATA = dlmread('ma001_NIRSS_LWC.txt') ;
%

% Figure Parameters
TL  = [-20 20] ;                % Skew-T Temperature Limits 
PL  = [400 1000] ;              % Skew-T Pressure Limits
AL  = [0 10] ;                  % Ascent Rate Plot Limits            
FL  = [35 45] ;                 % Frequency Plot Limits
DFDTL  = [-0.02 0.02] ;               % DFrequencyDT Plot Limits
SL  = [0.0 1.5] ;               % SLWC Plot Limits
TT  = -20:10:20 ;               % Skew-T Temperature Axis Tick Values
PT  = [100:100:1000,1050] ;     % Skew-T Pressure Axis Tick Values
AT  = 0:1:10 ;                  % Ascent Rate Plot Ascent Rate Tick Values    
FT  = 40:1:50 ;                 % Frequency Plot Frequency Axis Tick Values 
ST  = 0.0:0.2:1.4 ;             % SLWC Plot SLWC Axis Tick Values

% Skew-T Isotherm Line Parameters
AR  = [1 2 1] ;                 % Skew-T Aspect Ratio [x y z]
Z   = 2.5 ;                     % Skew-T Aspect Ratio Quick Scaling Factor
DT  = TL(2)-TL(1) ;             % Skew-T Distance Term 

% Guesses? WAGS? Maybe SWAGS? Meh...?
MVD = 20 ;

% 3005 Range Gate Interval
RG = 98.4252 ; % [ft]

% Colors!
br  = [0.8 0.5 0.0] ;           % Brown-ish
gr  = [0.0 0.5 0.0] ;           % Green-ish

% Constants
G     = 1.4 ;
MB2P  = 1.E2 ;
M2FT  = 3.28 ;
A1    = 17.625 ;
B1    = 243.04 ;

% Figure Subplot Titles
X1T = {'Temperature' ; sprintf('[ %cC ]',char(176))} ;
X2T = {'Supercooled Liquid' ; 'Water Content' ; '[ g/m3 ]'} ;
X3T = {'Frequency' ; '[ hz ]'} ;
Y1T = {'Pressure' ; '[ mbar ]'} ;
Y2T = {'Altitude' ; '[ 1000 ft MSL ]'} ;
Y3T = {'Time After Release' ; '[ seconds ]'} ;


% Skew-T Tick Labels
Y1L = {'100','200','300','400','500','600',...
       '700','800','900','1000','1050'} ;
Y2L = {'53','39','30','24','18','14','10','6','3','0.4',' '} ;

%% Release Determination
%TM_  = SDATA{5} ; I2 = length(TM_) ; i = 2 ; 
%while TM_(i) > TM_(i-1) ; i = i + 1 ; end

%% Bracket Sounding Data
%I1 = i ; N = I2-I1+1 ;

%% Extracting Sounding Data 
RH = humidity;
T  = temp_cor_c;
P  = pres;
F  = freq;
TS = secs_pmid;
A  = alt;
%TS_ = SDATA{ 4}(I1:I2) ; TS(1) = TS_(1) ;
%TM_ = SDATA{ 5}(I1:I2) ; TM(1) = TM_(1) ; 
%P_  = SDATA{ 7}(I1:I2) ; P(1)  = P_(1)  ; 
%T_  = SDATA{ 8}(I1:I2) ; T(1)  = T_(1)  ;
%RH_ = SDATA{10}(I1:I2) ; RH(1) = RH_(1) ;
%A_  = SDATA{17}(I1:I2) ; A(1)  = A_(1)  ;
%F_  = SDATA{39}(I1:I2) ; F(1)  = F_(1)  ;

%% Remove Repeated Frequency Telemetry and Corresponding Sounding Telemetry

%j = 0 ;
%k = 1 ;
%for i = 2:N  
%    j = j+1 ;   
%    if F_(i) ~= F_(i-1)      
%       k     = k+1    ;  
%       j     = 0      ;     
%       TS(k) = TS_(i) ; 
%       TM(k) = TM_(i) ; 
%       P(k)  = P_(i)  ; 
%       T(k)  = T_(i)  ;
%       RH(k) = RH_(i) ; 
%       A(k)  = A_(i)  ; 
%       F(k)  = F_(i)  ;   
%    elseif j == 3      
%       k     = k+1    ; 
%       j     = 0      ;
%       TS(k) = TS_(i) ; 
%       TM(k) = TM_(i) ; 
%       P(k)  = P_(i)  ; 
%       T(k)  = T_(i)  ;
%       RH(k) = RH_(i) ; 
%       A(k)  = A_(i)  ; 
%       F(k)  = F_(i)  ;     
%    end          
%end

% Total Records
%N = k ;
N = size(temp_cor_c);

% Set Seconds Time for Release
%TS = TS(:) - TS(1) ;

% Isotherms 
TI = -400:10:100 ;
PI = isotpres(TI,PL,DT,Z,AR) ;

% Skew-T,Log-P Subplot and SLWC Subplot
F1 = figure(1) ;

hold on

% Skew-T,Log-P Subplot
subplot_tight(1,4,4) ;

AX1 = gca ;
set(AX1,   'YScale','log',...                    % P-Axis Logarithmic Scale
           'YDir','reverse',...                  % Reversed P-Axis
           'XLim',TL,...                         % T-Axis Limits
           'YLim',PL,...                         % P-Axis Limits
           'XColor','k',...                      % T-Axis Color
           'YColor','k',...                      % P-Axis Color
           'box','on',...                        % Plot Box, set, ON! :P
           'PlotBoxAspectRatioMode','manual',... % Manual Box Aspect Ratio
           'PlotBoxAspectRatio',AR,...           % Box Aspect Ratio Value
           'XTick',TT,...                        % T-Axis Tick Mark Values
           'YTick',PT,...                        % P-Axis Tick Mark Values
           'TickDir','out',...                   % Inny or Outy? Right? :P
           'YTickLabel',Y1L) ;                   % P-Axis Tick Labels

AX2 = axes('Position',get(AX1,'Position'),...    % Overlay 2nd Axis on 1st
           'XAxisLocation','top',...             % 2nd Horz.-Axis Location
           'YAxisLocation','right',...           % Altitude Axis Location
           'Color','none',...                    % Axes Color: INVISIBLE!
           'XColor','k',...                      % T-Axis Color
           'YColor','k',...                      % A-Axis Color
           'YScale','log',...                    % A-Axis Logarithmic Scale
           'YDir','reverse',...                  % Reversed A-Axis
           'XLim',TL,...                         % T-Axis Limits
           'YLim',PL,...                         % A-Axis Limits
           'PlotBoxAspectRatioMode','manual',... % Manual Box Aspect Ratio
           'PlotBoxAspectRatio',AR,...           % Box Aspect Ratio Value
           'XTick',TT,...                        % T-Axis Tick Mark Values
           'YTick',PT,...                        % A-Axis Tick Mark Values
           'XTickLabel',{},...                   % T-Axis Tick Labels
           'TickDir','out',...                   % Inny or Outy? Dumb Joke.
           'YTickLabel',{}) ;                    % A-Axis Tick Labels

% Skew-T Axis Label Declaration       
xlabel(AX1,X1T) ; ylabel(AX1,Y1T) ; 

% Isobars (modifying GridLineStyle had no effect to achieve solid lines) 
for i = 2:10 ; line(TL,[PT(i) PT(i)],'Parent',AX1,'Color',br) ; end

% Isotherms
for i = 1:51 ; line([TI(i) 100],[1050 PI(i)],'Parent',AX1,'Color',br) ; end

% Dry Adiabats
for j = 1:15 ;   
    hold on
    DALR  = 9.8 ;
    PD    = flip(100:10:1050) ;
    TD(6) = TL(2)-10*j+50 ;
    ZD(6) = p2a(PD(6)*MB2P) ;   
    for i = 1:length(PD)        
        if PD(i) <= 1000.       
           ZD(i) = p2a(PD(i)*MB2P) ;
           TD(i) = TD(i-1)-DALR*(ZD(i)-ZD(i-1))/1E3 ;        
        else       
           k     = 6-i ;
           ZD(k) = p2a(PD(k)*MB2P) ;
           TD(k) = TD(k+1)+DALR*(ZD(k+1)-ZD(k))/1E3 ;                      
        end          
    end    
    for i = 1:length(PD) ; TDI(i) = isottemp(PD(i),TD(i),PL,DT,Z,AR) ; end    
    line(TDI,PD,'Parent',AX1,'Color',br,'LineStyle','--') ;        
    clear TDI TD PD ZD   
end

% Wet Adiabats 
for j = 1:15 ;    
    hold on
    PW    = flip(200:10:1050) ;
    TW(6) = 44-4*j ;
    PS(6) = vappres(TW(6)) ;
    WS(6) = mixrat(PW(6),PS(6)) ;
    ZW(6) = p2a(PW(6)*MB2P) ;
    for i = 1:length(PW)        
        if PW(i) <= 1000.            
           SALR  = salr(c2k(TW(i-1)),WS(i-1)) ;
           ZW(i) = p2a(PW(i)*MB2P) ;
           TW(i) = k2c(c2k(TW(i-1))-SALR*(ZW(i)-ZW(i-1))) ;
           PS(i) = vappres(TW(i)) ;
           WS(i) = mixrat(PW(i),PS(i)) ;           
        else            
           k     = 6-i ;       
           SALR  = salr(c2k(TW(k+1)),WS(k+1)) ;       
           ZW(k) = p2a(PW(k)*MB2P) ;
           TW(k) = k2c(c2k(TW(k+1))+SALR*(ZW(k+1)-ZW(k))) ;
           PS(k) = vappres(TW(k)) ;       
           WS(k) = mixrat(PW(k),PS(k)) ;           
        end
    end       
    for i = 1:length(PW) ; TWI(i) = isottemp(PW(i),TW(i),PL,DT,Z,AR) ; end     
    line(TWI,PW,'Parent',AX1,'Color',gr,'LineStyle',':') ;    
    clear TWI TW PW  ZW PS WS    
end

% Saturation Mixing Ratio Lines
WS = [0.1 0.4 1 2 5 10 20 40 60] ;
for j = 1:length(WS)    
    PW = flip(200:10:1050) ;    
    for i = 1:length(PW)        
        PS    = PW(i)*(WS(j)/(621.97+WS(j))) ;
        TW(i) = ((273.15*log10(PS/6.11))/(7.5-log10(PS/6.11))) ;                
    end    
    for i = 1:length(PW) ; TWI(i) = isottemp(PW(i),TW(i),PL,DT,Z,AR) ; end     
    line(TWI,PW,'Parent',AX1,'Color',gr,'LineStyle',':') ;   
    if TWI(3)>TL(1) && TWI(3)<TL(2)     
       text(TWI(3),1025,num2str(WS(j)),'Color',gr,'BackgroundColor','w',...
            'FontSize',5)     
    end     
    clear TWI TW PW       
end

% Sonde Skew-T,Log-P Data Transformation
for i = 1:N;
    DP(i)  = dewpnt(T(i),RH(i)) ;
    TI(i)  = isottemp(P(i),T(i),PL,DT,Z,AR) ;
    DPI(i) = isottemp(P(i),DP(i),PL,DT,Z,AR) ;
end

% Plot Skew-T,Log-P Data
line(TI,P,'Parent',AX1,'Color','b','LineWidth',1.5) ;
line(DPI,P,'Parent',AX1,'Color','r','LineWidth',1.5) ;
   
hold on


% SLWC Subplot
subplot_tight(1,4,3)

AX3 = gca ;
set(AX3,'YScale','log',...                       % P-Axis Logarithmic Scale
        'YDir','reverse',...                     % Reversed P-Axis
        'PlotBoxAspectRatioMode','manual',...    % Manual Box Aspect Ratio
        'PlotBoxAspectRatio',[1 2 1],...         % Box Aspect Ratio Value        
        'XLim',SL,...                            % SLWC-Axis Limits
        'YLim',PL,...                            % P-Axis Limits
        'XColor','k',...                         % SLWC-Axis Color
        'YColor','k',...                         % P-Axis Color
        'box','on',...                           % Plot Box, set, ON! :P
        'XTick',ST,...                           % SLWC-Axis Tick Values
        'YTick',PT,...                           % P-Axis Tick Mark Values
        'YTickLabel',{},...                      % P-Axis Tick Labels: none
        'TickDir','out');                        % Inny or Outy?
 
AX4 = axes('Position',get(AX3,'Position'),...    % Overlay 2nd Axis on 1st
           'XAxisLocation','top',...             % 2nd Horz.-Axis Location
           'YAxisLocation','right',...           % Altitude Axis Location
           'Color','none',...                    % Axes Color: INVISIBLE!
           'XColor','k',...                      % SLWC-Axis Color
           'YColor','k',...                      % A-Axis Color
           'YScale','log',...                    % A-Axis Logarithmic Scale
           'YDir','reverse',...                  % Reversed A-Axis
           'PlotBoxAspectRatioMode','manual',... % Manual Box Aspect Ratio
           'PlotBoxAspectRatio',[1 2 1],...      % Box Aspect Ratio Value           
           'XLim',SL,...                         % SLWC-Axis Limits
           'YLim',PL,...                         % A-Axis Limits
           'XTick',ST,...                        % SLWC-Axis Tick Values
           'YTick',PT,...                        % A-Axis Tick Mark Values
           'XTickLabel',{},...                   % SLWC-Axis Tick Labels
           'TickDir','out',...                   % Outy? Still Dumb Joke :P
           'YTickLabel',Y2L) ;                   % A-Axis Tick Labels

% Skew-T Axis Label Declaration       
xlabel(AX3,X2T) ; ylabel(AX4,Y2T) ;

% Isobars (modifying GridLineStyle had no effect to achieve solid lines) 
for i = 2:length(PT)-1 ; line([0 1],[PT(i) PT(i)],'Parent',AX3,'Color',br) ; end

% Supercooled Liquid Water Content Lines
for i = 2:length(ST)-1 ; line([ST(i) ST(i)],PL,'Parent',AX3,'Color',br) ; end

% Rejection Array 
B(1:N) = 0 ; B = B' ;

% Rejection Cutoff and Rejection Array 
NF = 44. ; TH = 1.1 ; CO = TH*NF ; B(F > CO) = 2 ;

% Averaging Interval
AI  = 30. ;  

% Initialize Average Value Data Array Sets
FA = F(1) ; TSA = TS(1) ; PA  = P(1) ; AA  = A(1) ;

% Determine Minimum Frequency Marking Maximum SLWC
[FML,FI] = min(F) ;

% Defining Values at Maximum SLWC
TSML = TS(FI) ; PML  = P(FI) ; AML  = A(FI) ;

% Defining Values at Maximum Telemetry
TSMT = max(TS(B == 0)) ; FMT = max(F(B == 0)) ; 
PMT  = max(P(B == 0))  ; AMT = max(A(B == 0)) ;

% Define Number of Intervals Prior to Maximum SLWC
IML = floor((TSML-TSA(1))/AI) ;

% Develop Time Interval Matrix up to Time at Maximum SLWC
TIM = [ TS(1) ; TS(1)+AI ] ; for i = 2:IML ; TIM(:,i) = TIM(:,i-1)+AI ; end

% Time Interval Matrix at Maximum SLWC 
TIM(1:2,end+1) = [ TIM(2,end) ; TSML ] ; 
TIM(1:2,end+1) = TSML ;
TIM(1:2,end+1) = [ TSML ; TSML+AI ] ;

% Define Number of Intervals Between Maximum SLWC and Maximum Telemetry
IMT = floor((TSMT-TSML)/AI) ;

% Develop Time Interval Matrix Between Maximum SLWC and Maximum Telemetry
for i = length(TIM)+1:length(TIM)+IMT ; TIM(:,i) = TIM(:,i-1)+AI ; end

% Time Interval Matrix at Maximum Telemetry
TIM(2,end) = TSMT ; TIM(1:2,end+1) = TSMT ;

% Averaging Accepted Telemetry
for i = 2:length(TIM)-1    
    if TIM(1:2,i) == TSML       
       TSA(i) = TSML ;
       FA(i)  = FML ;
       PA(i)  = PML ;
       AA(i)  = AML ;
    elseif TIM(1:2,i) == TSMT       
       TSA(i) = TSMT ;
       FA(i)  = FMT ; 
       PA(i)  = PMT ;
       AA(i)  = AMT ;
    else    
       %disp(i)
       TSA(i) = mean(TS(TS > TIM(1,i) & TS <= TIM(2,i) & B == 0)) ;
       FA(i)  = mean(F (TS > TIM(1,i) & TS <= TIM(2,i) & B == 0)) ;
       PA(i)  = mean(P (TS > TIM(1,i) & TS <= TIM(2,i) & B == 0)) ; 
       AA(i)  = mean(A (TS > TIM(1,i) & TS <= TIM(2,i) & B == 0)) ;
    end   
end

% Remove Any NaN Values in Average Arrays Resulting from mean(Null)
FA(isnan(TSA)) = [] ; PA(isnan(TSA))  = [] ;
AA(isnan(TSA)) = [] ; TSA(isnan(TSA)) = [] ;

% Frequency Time Derivative
for i = 1:length(FA)-1 ; DFDT(i) = (FA(i+1)-FA(i))/(TSA(i+1)-TSA(i)) ; end

% Setting the Frequency Time Derivative for Maximum Telemetry
DFDT(end+1) = DFDT(end) ;

% Collection Efficiency Look-up Table (Anasphere Final Report NNX11CD74P)  
CE = [ 0    0.03 ; 3    0.16 ; 4    0.32 ; 5    0.43 ; 6    0.52 ;
       7    0.59 ; 8    0.65 ; 9    0.70 ; 10   0.74 ; 11   0.77 ;
       12   0.80 ; 13   0.82 ; 13.5 0.85 ; 14   0.86 ; 15   0.87 ;
       16   0.88 ; 17   0.89 ; 18   0.90 ; 19   0.91 ; 20   0.92 ;
       21   0.93 ; 22   0.94 ; 25   0.95 ; 30   0.96 ; 35   0.97 ;
       40   0.98 ; 45   0.99 ; 55   0.99 ; 100  1.00 ] ;

% Interpolate Collection Efficiency from MVD Guess
E = interp1(CE(:,1),CE(:,2),MVD) ;

% SLWC Calculation Parameters
F0 = 44.0    ; % Nominally
B0 = 6.8528  ; % Hill and Wolfinden (1980)
B1 = 44.366*100  ;
L  = 0.070  ; % [m] Anasphere Final Report NNX11CD74P p.20  
D  = 0.00063 ; % [m] Anasphere Final Report NNX11CD74P p.20

% SLWC Calculations
for i = 1:length(DFDT)-1 ; 
    %disp(i)
    SLWC(i)  = -(2*F0^2)/(E*B0*L*D*AA(i)*FA(i)^3)*DFDT(i) ;
end

% Negative SLWC Values are NOT physical! Begone! :P
SLWC(SLWC < 0) = 0. ;

% Noise Threshold (Need to Justify this but it makes things pretty for now)
SLWC(PA > 900) = 0.    ;
SLWC(PA < 500) = 0.    ;
SLWC(SLWC < 0.02) = 0. ;

% NIRSS LWC Profile
% see how NIRSS data files are loaded above
NDATA(2,1) = 225.*3.28084 ; 
for i = 2:485 ; NDATA(2,i) = NDATA(2,i-1)+RG ; end
for i = 1:485 ; NDATA(3,i) = a2p(NDATA(2,i)/3.28084)/100 ; end

% Plot SLWC Data
%patch([SLWC(12:15) 0],[PA(12:15) PA(15)],'b','Parent',AX3) ;
line(SLWC,PA,'Parent',AX3,'Color','b','LineWidth',1.5) ;
%line(SLWC2,PA,'Parent',AX3,'Color','k','LineWidth',1.5) ;
line(NDATA(1,:),NDATA(3,:),'Parent',AX3,'Color','r','LineWidth',1.5) ;



%  freq Subplot
subplot_tight(1,4,1)

AX3 = gca ;
set(AX3,'YScale','log',...                       % P-Axis Logarithmic Scale
        'YDir','reverse',...                     % Reversed P-Axis
        'PlotBoxAspectRatioMode','manual',...    % Manual Box Aspect Ratio
        'PlotBoxAspectRatio',[1 2 1],...         % Box Aspect Ratio Value        
        'XLim',FL,...                            % Freq-Axis Limits
        'YLim',PL,...                            % P-Axis Limits
        'XColor','k',...                         % SLWC-Axis Color
        'YColor','k',...                         % P-Axis Color
        'box','on',...                           % Plot Box, set, ON! :P
        'XTick',ST,...                           % SLWC-Axis Tick Values
        'YTick',PT,...                           % P-Axis Tick Mark Values
        'YTickLabel',{},...                      % P-Axis Tick Labels: none
        'TickDir','out');                        % Inny or Outy?
        
line(F(1:5:end-1),PA,'Parent',AX3,'Color','b','LineWidth',1.5) ;


%  DFDT Subplot
subplot_tight(1,4,2)
AX3 = gca ;
set(AX3,'YScale','log',...                       % P-Axis Logarithmic Scale
        'YDir','reverse',...                     % Reversed P-Axis
        'PlotBoxAspectRatioMode','manual',...    % Manual Box Aspect Ratio
        'PlotBoxAspectRatio',[1 2 1],...         % Box Aspect Ratio Value        
        'XLim',DFDTL,...                         % DFDT-Axis Limits
        'YLim',PL,...                            % P-Axis Limits
        'XColor','k',...                         % SLWC-Axis Color
        'YColor','k',...                         % P-Axis Color
        'box','on',...                           % Plot Box, set, ON! :P
        'XTick',ST,...                           % SLWC-Axis Tick Values
        'YTick',PT,...                           % P-Axis Tick Mark Values
        'YTickLabel',{},...                      % P-Axis Tick Labels: none
        'TickDir','out');                        % Inny or Outy?

line(DFDT,PA,'Parent',AX3,'Color','b','LineWidth',1.5) ;

subplot_tight(1,4,4) ;
% Skew-T Axis Label Declaration       
xlabel(AX1,X1T) ; ylabel(AX1,Y1T) ; 

tightfig ;

%C = floor(FML)-1-exp(1) ;
 
%TTICK = 0:100:1200 ;

%FTICK = log([floor(FML)-1 44 50:10:100]-C) ;

%FTLABEL = {floor(FML)-1 44 50:10:100} ;

%figure(2)

%hold on

%subplot(1,2,1)

%AX5 = gca ;

%set(AX5,   'XLim',[log(exp(1)) log(100-C)],...         % T-Axis Limits
%           'YLim',[0 1200],...                         % P-Axis Limits
%           'XColor','k',...                      % T-Axis Color
%           'YColor','k',...                      % P-Axis Color
%           'box','on',...                        % Plot Box set, ON! :P
%           'XTick',FTICK,...                     % T-Axis Tick Mark Values
%           'YTick',TTICK,...                     % P-Axis Tick Mark Values
%           'TickDir','out',...                   % Inny or Outy? Right? :P
%           'XTickLabel',FTLABEL) ;                   % P-Axis Tick Labels

%% Skew-T Axis Label Declaration       
%xlabel(AX5,X3T) ; ylabel(AX5,Y3T) ;

%for i = 3:length(FTICK)-1 ; 
%    line([FTICK(i) FTICK(i)],[0 1200],'Parent',AX5,'Color',br) ; 
%end
%for i = 2:length(TTICK)-1 ; 
% line([FTICK(1) FTICK(end)],[TTICK(i) TTICK(i)],'Parent',AX5,'Color',br) ; 
%end
       
%line(log(F(B==0)-C),TS(B == 0),'Color','b','MarkerSize',4,'LineStyle','-.')
%line(log(F(B~=0)-C),TS(B ~= 0),'Color','b','MarkerSize',4,'LineStyle','o')
%line(log(FA-C),TSA,'Color','r','LineWidth',1.5)
 
