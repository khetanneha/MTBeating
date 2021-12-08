% ----------------------------------------------------------------
% Neha Khetan, July 2021
% Athale lab, for Filament beating analysis from experiments of SY
% 1. Comparison b/w different methods - interpolation and smoothing
% 2. Input File types: 
%     a. Excel file - SAY: NeuronJ
%     b. Mat file -   DK:  tracking
% ----------------------------------------------------------------
clear all, clc, clf, close all;

ipath       = './' %'./../exptdata/input_files/';


%opath       = './exptplots/exp10/';
mkdir('./exptplots/');
opath = './exptplots/'

%% WHether data is from SY or DK
% Choose 0 = fro data from SY: excel - sheet (.xlsx) format from SY
%        1 = for data from DK: Filament tracking - .mat format
selectType = 0;

% To Interpolate points: for improving data quality ?
%   Methods: Linear, polynomial, spline
InterPolatePolyfit  = 0  ; % for polynomial based
InterPolateSpline   = 0  ;
InterPolateLinear   = 0 ;

% INPUT FILES from SY: Excel sheets: NeuronJ
% filename      = 'Coordinates.xlsx'; % col 5 and col 6
% filename      = 'antibody.xlsx';
% filename      = 'BiotinStrept.xlsx'
filename        = 'FinalBS.xlsx'

% if Dhruv data
DK_data         = load( './coordinates_v2.mat' );
 

% ----------------------------------------------------------
DoIKnowLenExpt = 0; % No: use au. vals - modify later
inputLen      =  4.469 %7.5; % um;
PivLength     =  0 ;


tanval ={}; curv ={}; BeatWidth = []; CounterFil = 0;
maxN   =[]; sumCurv = []; FreeTipAngle =[]; MTip =[]; End2EndDist = []; ds = [];

nFrame      = 31;
px2mu       = 1; %106/1000; % 1 px is 106 nm
dT          = 10; % 10 s
Time        = nFrame * dT;
Timepts     = [1:nFrame].*dT;




if (selectType == 0)
    
    yColNo = 5;
    xColNo = 6;
    % ---------------------
    FN                  = strcat( ipath , filename  );
    [ status ,sheets]   = xlsfinfo( FN );
    nFrame              = length( sheets );
    px2mu               = 1;
    PivLength           = 0 ;
    Timepts             = [1:nFrame].*dT;
    Time                = nFrame * dT;
    % ---------------------
    colc = colormap( turbo(nFrame + 1)); cc = 1;
    for k = 1:nFrame
        
        fn          = ['Sheet', sprintf('%i', k )];
        data        = xlsread( [ ipath , filename] , fn  );
        
        yy   = data( : , yColNo ).*px2mu ;
        xx   = data( : , xColNo ).*px2mu ;
        
        if k==1
            x0  = xx; y0 = yy;
            [ xx , yy ] = dataTransform( xx , yy   );
            %[ xx , yy ]  = dataTransformXAxis( xx , yy  , x0 , y0 );
            V0          = [ ( xx(end) - xx(1)) , ( yy(end) -yy(1))  ];
            
        else
            [ xx , yy ] = dataTransform( xx , yy   );
            %[ xx , yy ]  = dataTransformXAxis( xx , yy  , x0 , y0 );
        end
        
        %% Increase data points along contour -
        % work around the coarsed data-set in experiments
        if (InterPolateSpline)
            %Spline interpolation
            xUnique= handleXDuplicate( xx );
            xNew =[ 0:0.5:max(xUnique) ];
            xNew  = linspace(min(xUnique) , max(xUnique) , length(xx)) ;
            yNew = interp1( xUnique, yy , xNew , 'spline');
        end
        
        if  (InterPolateLinear)
            [ xx , yy , flag ] = dataInterpolate( xx , yy  );
            if flag
                [ xx , yy , flag ] = dataInterpolate( xx , yy  );
            end
        end
        
        if InterPolatePolyfit
            % polynomial fitting
            xNew = linspace( min(xx) , max(xx) , 2*length(xx) );
            %xNew = [min(xx):0.2: max(xx)]; % , length(xx) );
            P    = polyfit( xx , yy , 3 );
            yNew = polyval( P, xNew );
            yy = yNew';
            yy = yy(:)-yy(1);
            xx = xNew(:)-xNew(1);
        end
        
        %% ---
        BeatWidth  =[ BeatWidth ; ( abs(min(xx )) + abs(max(xx)) ) , ( abs(min(yy)) + abs(max(yy)) ) ];
        
        totalXX = xx + PivLength;
        
        
        if rem( CounterFil , 3 ) == 0
            figure(1),...
                plot( totalXX , yy , 'color' , colc(CounterFil + 1,:),  'linewidth' , 1.6   ), hold on,...
                
        end
    
              
            
    
    CounterFil = CounterFil + 1;

    
    %% ---- compute tanAngles
    % Vect: Fil. seg from pivot to free tip
    V1    = [ xx , yy ];
    tmp   = tanAngle( V1  );
    tanval{k} = tmp';
    maxN      = [ maxN ; length( tmp ) ];
    % curv = 1/R ~ dtheta/dl
    tmp2    = diff( tmp );
    dl      = sqrt(    (V1(2:end,1) - V1(1:end-1,1) ).^2  + (V1(2:end,2) - V1(1:end-1,2) ).^2 );
    curv{k} = (tmp2./dl(1:end-1))';
    sumCurv = [ sumCurv ; k , sum( (tmp2./dl(1:end-1)) ) ];
    ds = [ ds; dl];
    
    
    % for tip angles  % compute angles free end wrt pivot
    VFree        = [  ( V1(end, 1) - V1(1, 1) ) , ( V1(end,2) - V1(1,2))  ];
    FreeTipAngle = [ FreeTipAngle ; k , tangentAngle( V0 , VFree )*180/pi ];
    %FreeTipAngle = [ FreeTipAngle ; k , AnglePosX( VFree )*180/pi ];
    
    
    % -----
    MTip        = [ MTip ; k , V1(end,1) , V1( end , 2 ) ];
    End2EndDist = [ End2EndDist; sqrt( V1(end,1)^2 + V1( end , 2 )^2 ) ];
    end
    
    
    
    
else
    %% Dhruv Data Format
    inputLen = inputLen*px2mu;
    d  = DK_data;   
       fn = d.coordinates(:,3);
    xx = d.coordinates(:,1);
    yy = d.coordinates(:,2);
    % -------------------------------------------------------
    data = [ fn , xx , yy ];
    ti = min( fn );
    tf = max( fn );
    % --------------------------------------------------------
    dT          = 30; % 30 s
    Frames      = unique( fn );
    nFrame      = length( Frames )
    Timepts     = [1:nFrame].*dT;
    px2mu       = 106/1000; % 1 px is 106 nm
    PivLength   = 0;  % Since, in DK method of tracking - pivot is also included
    CounterFil  = 0;
    
    colc = colormap( turbo(nFrame + 1)); cc = 1;
    for k = 1:nFrame
        
        indx = data(:,1) == Frames(k);
        
        yy   = data( indx , 3 ).*px2mu ;
        xx   = data( indx , 2 ).*px2mu ;
        
        if k==1
            x0 = xx; y0 = yy;
            [ xx , yy ] = dataTransform( xx , yy  );
            V0          = [ ( xx(end) - xx(1)) , ( yy(end) -yy(1))  ];
        else
            [ xx , yy ] = dataTransform( xx , yy   );
        end
        %% Increase data points along contour -
        % work around the coarsed data-set in experiments
        if (InterPolateSpline)
            %Spline interpolation
            xUnique= handleXDuplicate( xx );
            xNew =[ 0:0.5:max(xUnique) ];
            xNew  = linspace(min(xUnique) , max(xUnique) , length(xx)) ;
            yNew = interp1( xUnique, yy , xNew , 'spline');
        end
        
        if  (InterPolateLinear)
            [ xx , yy , flag ] = dataInterpolate( xx , yy  );
            if flag
                [ xx , yy , flag ] = dataInterpolate( xx , yy  );
            end
        end
        
        if InterPolatePolyfit
            % polynomial fitting
            xNew = linspace( min(xx) , max(xx) , 2*length(xx) );
            %xNew = [min(xx):0.2: max(xx)]; % , length(xx) );
            P    = polyfit( xx , yy , 3 );
            yNew = polyval( P, xNew );
            yy = yNew';
            yy = yy(:)-yy(1);
            xx = xNew(:)-xNew(1);
        end
        
        BeatWidth  =[ BeatWidth ; ( abs(min(xx )) + abs(max(xx)) ) , ( abs(min(yy)) + abs(max(yy)) ) ];
        
        totalXX = xx + PivLength;
        % for comparison at the interval of 30 s
        %         if rem( CounterFil , 3 ) == 0
        %             figure(1),...
        %                 plot( totalXX , yy ,'linewidth' , 1.4  ), hold on,...
        %         end
        
        figure(1),...
           plot( totalXX , yy , 'color' , colc(CounterFil + 1,:),  'linewidth' , 1.6   ), hold on,...
           
        
        CounterFil = CounterFil + 1;
        
        
        
        %% ---- compute tanAngles
        % Vect: Fil. seg from pivot to free tip
        V1    = [ xx , yy ];
        tmp   = tanAngle( V1  );
        
        tanval{k} = tmp';
        maxN      = [ maxN ; length( tmp ) ];
        
        % curv = 1/R ~ dtheta/dl
        tmp2    = diff( tmp );
        dl      = sqrt(    (V1(2:end,1) - V1(1:end-1,1) ).^2  + (V1(2:end,2) - V1(1:end-1,2) ).^2 );
        curv{k} = (tmp2./dl(1:end-1))';
        sumCurv = [ sumCurv ; k , sum( (tmp2./dl(1:end-1)) ) ];
        
        ds = [ ds; dl];
        
        
        % for tip angles  % compute angles free end wrt pivot
        VFree        = [  ( V1(end, 1) - x0(1) ) , ( V1(end,2) - y0(1))  ];
        FreeTipAngle = [ FreeTipAngle ; k , tangentAngle( V0 , VFree )*180/pi ];
        %FreeTipAngle  = [ FreeTipAngle ; k , AnglePosX( VFree )*180/pi ];
        % -----
        MTip        = [ MTip ; k , V1(end,1) , V1( end , 2 ) ];
        End2EndDist =  [ End2EndDist; sqrt( V1(end,1)^2 + V1( end , 2 )^2 ) ];
    end
end

%% -----------------------------------------------------------------------------------
%
%---------------------------------------------------------------------------------------


filLength = ( round(max(End2EndDist)) + PivLength);

figure(1),...
    %xlim([ 0 filLength ]),...
    c = colorbar;
    c.Label.String = 'Time (s)',...
    caxis( [ 0 Time ]),...
    xlabel('X (\mum)'),...
    ylabel('Y (\mum)'),...
    set( gca , 'fontsize' , 24 ),...
    axis equal,...
    print( gcf, '-dpdf', [ opath , 'exptXY_', sprintf('%s', filename) , '.pdf'], '-r300' )
%%  ------------------------------------------------------------------------------------------------------------------------------------



Fh        = 2; % reshaping
maxN      = max( maxN );
fL        = maxN;
Mulfactor = round(max(End2EndDist))/fL;
xval      = [0:fL]%.*Mulfactor;

MatTanAng   = zeros( nFrame , maxN );
MatTanAng(:,:) = NaN;
[ m1 , n1 ] = size( tanval );
for kk=1:n1
    tmp = tanval{1,kk };
    MatTanAng( kk , 1:length( tmp ) ) = tmp;
end



XVAL = (xval + PivLength );
xMaxLabel = ( xval(end) + PivLength);

if (DoIKnowLenExpt)
    mulfac = inputLen/xMaxLabel;
else
    mulfac = 1;
end


figure( Fh ),...
    imagesc( XVAL.*mulfac , Timepts , MatTanAng , 'AlphaData',~isnan(MatTanAng)),...
    colorbar,...
    xlabel('L_{MT} (a.u., in terms data points)'),...% (\mum)'),...
    ylabel('Time (s)'),...
    set( gca , 'fontsize' , 24 ),...
    %title('Tangent angle \psi'),...
xlim([ 0 xMaxLabel ] ),...
    print( gcf, '-dpdf', [ opath , 'PSI_', sprintf('%s', filename) , '.pdf'], '-r300' )

% ------------------------------------------------------------------------
MatCurv     = zeros( nFrame , maxN );
[ m2 , n2 ] = size( curv );
for k2=1:n2
    tmp = curv{1,k2 };
    %ind = ((tmp(:) <=-10 )  );
    %tmp(ind) = 0;
    MatCurv( k2 , 1:length( tmp ) ) = tmp;
end
% [ mean(MatCurv(:)) , std(MatCurv(:)) ,  median(MatCurv(:)) ]
figure(Fh+1),...
    imagesc( XVAL.*mulfac ,Timepts , (MatCurv) ),...
    colorbar,...
    xlabel('L_{MT} (a.u., in terms data points)'),...% (\mum)'),...
    ylabel('Time (s)'),...
    set( gca , 'fontsize' , 24 ),...
    title('Curvature'),...
    % (rad/\mum),...
xlim([ 0 xMaxLabel ] ),...
 print( gcf,  '-dpdf', [ opath , 'Curv_', sprintf('%s', filename) , '.pdf'],  '-r300' )


figure( Fh+2 ),...
    plot( FreeTipAngle(:,1).*dT , FreeTipAngle(:,2) , 'k-', 'linewidth' , 2 ),...
    xlabel('Time (s)'),...
    set( gca , 'fontsize' , 24 ),...
    xlim([ 0 Time  ]),...
    ylabel('\phi '),...
    print( gcf,  '-dpdf', [ opath , 'FreeTipAng_', sprintf('%s', filename) , '.pdf'], '-r300' )


%% ========================================================================
MTip(:,2) = ( MTip(:,2) -MTip(1,2) );
MTip(:,3) = ( MTip(:,3) -MTip(1,3) );
MTip(:,1) = [ MTip(:,1).*dT   ];
tmpX = max(   abs( min( MTip(:,2) )) , max( MTip(:,2) ) );
tmpY = max(   abs( min( MTip(:,3) )) , max( MTip(:,3) ) );
tmpMax = round(max( tmpX , tmpY ));
figure( Fh + 3 ),...
    h1 = plot( MTip(:,1) , MTip(:,2) , 'linewidth' , 2 ),hold on,...
    h2 = plot( MTip(:,1) , MTip(:,3)  , 'linewidth' , 2 ),...
    ylabel( 'Position(\mum)'),...
    xlabel( 'Time (s)'),...
    xlim( [ 0 MTip(end,1)]),...
    ylim( [ -tmpMax tmpMax ]),...
    set( gca , 'fontsize' , 24),...
    legend( [ h1, h2  ] , 'X', 'Y'  , 'location' , 'best' , 'fontsize' , 10);
    print( gcf,  '-dpdf', [ opath , 'FreeTipXY_', sprintf('%s', filename) , '.pdf'], '-r300' )

figure( Fh + 4 ),...
    [ yF , xB ] = hist( ds ),...
    bar( xB , yF./sum( yF) , 1 ),...
    xlabel('\delta s (\mum)'),...
    ylabel('Frequency'),...
    set( gca, 'fontsize', 18 ),...
    print( gcf,  '-dpdf', [ opath , 'dsExpt_', sprintf('%s', filename) , '.pdf'], '-r300' )



QuantVals = [ max(  BeatWidth(:,1) ) ,  max(  BeatWidth(:,2) ) , max( abs(diff(FreeTipAngle(:,2)))) ];
dlmwrite( [ opath , 'OutStats_', sprintf('%s', filename) , '.out' ] ,QuantVals ,  'delimiter', '\t' );



%% ================================================================================================================================
function [ xx , yy ] = dataTransform( x , y , x0 , y0)
xx = x(:) -x(1);
yy = y(:) -y(1);
end

function [ xx , yy ] = dataTransformXAxis( x , y , x0 , y0)
xx = x(:) -x0(:)';
yy = y(:) -y0(:)';
end


function   [ X , Y , DOi ] = dataInterpolate( fX , fY   )
DOi = 0;


dx = fX(2:end) -fX(1:end-1);
dy = fY(2:end) -fY(1:end-1);

dr = sqrt( dx.^2 + dy.^2 );

DistCutoff =  0.1; %200;
X = fX(1); Y =fY(1);
for rk= 1:length( dr )
    kk = rk + 1;
    
    if dr(rk) > DistCutoff
        xnew  = 0.5*( fX( kk ) + fX(kk-1) );
        ynew  = 0.5*( fY( kk ) + fY(kk-1) );
        
        X = [ X ;  xnew  ; fX( kk ) ];
        Y = [ Y ;  ynew  ; fY( kk ) ];
    else
        X = [ X ; fX( kk )];
        Y = [ Y ; fY( kk )];
    end
end
dx = X(2:end) -X(1:end-1);
dy = Y(2:end) -Y(1:end-1);
dr2 = sqrt( dx.^2 + dy.^2 );

if( dr2(dr2 > DistCutoff ) )
    DOi = 1;
end

end

function theta = tanAngle( vv  )

dy = [ vv(2:end , 2 ) - vv(1:end-1, 2)];
dx = [ vv(2:end , 1 ) - vv(1:end-1, 1)];
theta =[];
for k = 1:length( dy )
    theta = [ theta ; atan2(dy(k), dx(k) ) ];
end
end




function ang = tangentAngle( V1 , V2 )

ang = [];

num = ( V1(:,1).*V2(:,2) ) - ( V2(:,1).*V1(:,2) );
den = ( V1(:,1).*V2(:,1) ) + ( V1(:,2).*V2(:,2) );
for i = 1:length( num )
    ang = [ ang ; atan2(  num(i) , den(i) ) ];
    
end
end

function ang = AnglePosX( V1  )

ang = atan2( V1(2), V1(1) );
end

function V= handleXDuplicate( V )


n = length( V );

tt = [];

for k = 1:n
    v1 = V(k)
    if( ismember( v1 , V ))
        indx = find( V == v1 );
        
        for i = 2:length( indx )
            V( indx(i) ) = V( indx(i) ) + 0.01;
        end
    end
end

end
