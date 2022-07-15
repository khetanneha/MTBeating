% Neha Khetan, 2021 -2022
% Computes the following:
% 1] Tangent angle: computes, plots kymograph, writs into file
% 2] Curvature : computes, plots kymograph and writes into file
% 3] Tip-angle : Angle at time 't' wrt to t=0 of the free-tip with the
% fixed end
function [ VectPiv , PivCoor , TipAngVel , BWInfo ,  TipAng] = tang_curv2(  Fh ,  f1 , filSeg , dT , fL ,  ST , ET , PivLen , filename , opath )


% Fiber coors: % Time identity      posX      posY    forceX    forceY   tension
%                   1                3          4       5         6 %      7
tottime = ET - ST;
Timepts = [ ST:dT:ET ];

nSeg = fL/filSeg;
nSeg1um = 1/filSeg;
iPiv = (PivLen/filSeg)+1;


tanval = {};  curv ={};
maxN = [];
sumCurv = [];
FreeTipAngle =[];


xx      = [  PivLen:0.1:(fL)   ];

ForceKymX = []; ForceKymY = []; NetForce = [];

colormap turbo
ds = [];
% -------------------------------------------------------------------------
plotTime =[ Timepts(:) - Timepts(1)]; % Time starts from '0'


for k=1:length( Timepts )
    
    Key = Timepts( k );
    data = f1(f1( :, 1) == Key , 3:7);
    
    
    fXt  = data( : , 1 );
    fYt  = data( : , 2 );
    
    if (k == 1)
        minX = abs( min( fXt) );
    end
    
    
    ForceX = data(:, 3); ForceY = data(:, 4);
    
    % Swapping rep. i.e. ( 0, 0 ) is plus end ;
    % direction is plus -minus end
    fXt    = -fXt; fXt = fXt + minX;
    fXt    = fXt( end:-1:1); fYt  = fYt( end:-1:1);
    ForceX = ForceX(end:-1:1); ForceY = ForceY(end:-1:1)
    
    
    fXt = fXt( iPiv:end) - fXt( iPiv);
    fYt = fYt( iPiv:end) - fYt( iPiv);
    
    % Vect: Fil. seg from pivot to free tip
    V1     = [ fXt , fYt ];
    tmp    = tanAngle( V1  );
    
    tanval{k} = tmp';
    maxN      = [ maxN ; length( tmp ) ];
    
    % curv = 1/R ~ dtheta/dl
    tmp2    = diff( tmp );
    dl      = sqrt(    (V1(2:end,1) - V1(1:end-1,1) ).^2  + (V1(2:end,2) - V1(1:end-1,2) ).^2 );
    curv{k} = (tmp2./dl(1:end-1))';
    sumCurv = [ sumCurv ; plotTime(k) , sum( (tmp2./dl(1:end-1)) ) ];
    ds      = [ ds ; dl' ];
    
    % for tip angles
    %-------
    if ( k == 1 )
        VectPiv = [     ( fXt( end ) - fXt( iPiv ) ) , ( fYt( end ) - fYt( iPiv ) ) ];
        PivCoor = [     fXt( iPiv )  , fYt( iPiv ) ];
        %         FreeTipAngle = [ FreeTipAngle ; plotTime(k) , 0 , fXt(end) , fYt(end) ];
        %     end
    end
    
    
    % compute angles free end wrt pivot
    VFree        = [  fXt(end) , fYt(end)  ];
    FreeTipAngle = [ FreeTipAngle ; plotTime(k) , AnglePosX(VFree)*180/pi , fXt(end) , fYt(end) ];
    % AngleV1V2( VectPiv , VFree )*180/pi];
    % Since, in sim. at t =0; filament is +ve X-axis
    % Angle , at the tip of immobile filament as in EXPT: we do not have
    % info on pivot end;
    % If: in expt we have info - the angle subtended at the pivot site can
    % be reverted to
    
    
    % Matrix for Motor forces on MT
    tmp       = sqrt( ForceX.^2  +  ForceY.^2 );
    ForceKymX = [ ForceKymX ; (ForceX)']; %[ ForceKymX ; ForceX'];
    ForceKymY = [ ForceKymY ; (ForceY)']; %[ ForceKymY ; ForceY'];
    NetForce  = [ NetForce ; plotTime( k ) , sqrt( sum(ForceX).^2 + sum(ForceY).^2 ) , length( tmp(tmp~= 0) ) ];
    
    
    
end

% reshaping
maxN = max( maxN );
MatTanAng   = zeros( length( Timepts ) , maxN );
[ m1 , n1 ] = size( tanval );
for kk=1:n1
    tmp = tanval{1,kk };
    MatTanAng( kk , 1:length( tmp ) ) = tmp;
end

figure( Fh ),...
    imagesc( xx , plotTime , MatTanAng ),...
    colorbar,...
    xlabel('L_{MT}'),...
    ylabel('Time (s)'),...
    set( gca , 'fontsize' , 18 ),...
    title('\psi'),...
    xlim([ 0 fL ]),...
    export_fig( gcf, [ opath , 'tangentAngle_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');

% ------------------------------------------------------------------------
MatCurv     = zeros( length( Timepts ) , maxN );
[ m2 , n2 ] = size( curv );
for k2=1:n2
    tmp = curv{1,k2 };
    %ind = ((tmp(:) <=-10 )  );
    %tmp(ind) = 0;
    MatCurv( k2 , 1:length( tmp ) ) = tmp;
end

figure(Fh+1),...
    imagesc( xx , plotTime , (MatCurv) ),...
    colorbar,...
    xlabel('L_{MT} (\mum)'),...
    ylabel('Time (s)'),...
    set( gca , 'fontsize' , 18 ),...
    title('Curvature'),... % (rad/\mum),...
    xlim([ 0 fL ]),...
    export_fig( gcf, [ opath , 'Curvature_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');


figure(Fh+1000),...
    imagesc( xx./fL , plotTime , (MatCurv) ),...
    colorbar,...
    xlabel('L_{MT} (\mum)'),...
    ylabel('Time (s)'),...
    set( gca , 'fontsize' , 18 ),...
    title('Curvature'),... % (rad/\mum),...
    xlim([ 0 1 ]),...
    export_fig( gcf, [ opath , 'NormCurvature_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');


%Kymographs
figure(Fh+2),...
    imagesc( xx , plotTime , (ForceKymX) ),...
    colorbar,...
    xlabel('L_{MT} (\mum)'),...
    ylabel('F_{x} (pN)'),...
    set( gca , 'fontsize' , 18 ),...
    xlim([ 0 fL ]),...
    %export_fig( gcf, [ opath , 'ForceX_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');

figure(Fh+3),...
    imagesc( xx , plotTime , (ForceKymY) ),...
    colorbar,...
    xlabel('L_{MT} (\mum)'),...
    ylabel('F_{y} (pN)'),...
    set( gca , 'fontsize' , 18 ),...
    xlim([ 0 fL ]),...
    %export_fig( gcf, [ opath , 'ForceY_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');

% figure(Fh+4),...
%     yyaxis left,...
%     plot( NetForce(:,1) , sum(ForceKymX,2), 'linewidth' , 2 ), hold on,...
%     xlabel('Time (s)'),...
%     ylabel('F_{x} (pN)'),...
%     set( gca , 'fontsize' , 18 ),...
%     yyaxis right,...
%     plot( NetForce(:,1) , sum(ForceKymY,2), 'linewidth' , 2 ), hold on,...
%     xlabel('Time (s)'),...
%     ylabel('F_{y} (pN)'),...
%     set( gca , 'fontsize' , 18 ),...
%     xlim([ 0 NetForce(end,1)  ]),...
%     %set(gca,'ColorScale','log')
% export_fig( gcf, [ opath , 'ForceXY_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');

figure( Fh + 5 ),...
    yyaxis left,...
    plot( NetForce(:,1) , NetForce(:,2) , 'linewidth' , 2 ), hold on,...
    xlabel('Time (s)'),...
    ylabel('Force (pN)'),...
    set( gca , 'fontsize' , 18 ),...
    xlim( [ 0 tottime ] ),...
    yyaxis right,...
    plot( NetForce(:,1) , NetForce(:,3) , 'linewidth' , 2 ), hold on,...
    xlabel('Time (s)'),...
    ylabel('N'),...
    set( gca , 'fontsize' , 18 ),...
    xlim( [ 0 tottime ] ),...
    export_fig( gcf, [ opath , 'NetForce_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');

% these files can be huge with density; save if necessary else SKIP
dlmwrite(  [ opath , 'MotForceNetNum', sprintf('%s', filename ) ,'.out'] , NetForce ,  'delimiter' , '\t', 'precision','%.6f' );
% dlmwrite(  [ opath , 'DatKymForceX', sprintf('%s', filename ) ,'.out'] , ForceKymX ,  'delimiter' , '\t', 'precision','%.6f' );
%dlmwrite(  [ opath , 'DatKymForceY', sprintf('%s', filename ) ,'.out'] , ForceKymY ,  'delimiter' , '\t', 'precision','%.6f' );


%- -------------------------------------------------------------------
Fh = Fh + 7;
%% Plot Free tip of the filament
[ Fh , TipAngVel , BWInfo  ] = FreeTipFil( FreeTipAngle ,  (Fh) , filename , opath );
%% FFT of tangent angle
Fh = computeFFT( xx , plotTime , MatTanAng , nSeg1um , (Fh + 1 ) , filename , opath );
%% Compute energy on filament
% Fh = computeEnergy( xx , plotTime , MatCurv , ds ,   nSeg1um , (Fh + 1 ) , filename , opath );

TipAng  = FreeTipAngle( : , 2 );
dlmwrite([opath, 'data_psi_', sprintf('%s', filename), '.out'] ,  MatTanAng , 'delimiter' , '\t');
dlmwrite([opath,'data_kappa_', sprintf('%s', filename), '.out'] ,  MatCurv , 'delimiter' , '\t');


% figure( Fh*40 ),...
%     [ yF , xB ] = hist( ds ),...
%     bar( xB , yF./sum( yF) , 1 ),...
%     xlabel('\delta s (\mum)'),...
%     ylabel('Frequency'),...
%     set( gca, 'fontsize', 18 ),...
% export_fig( gcf, [ opath , 'dsSim_', sprintf('%s', filename) , '.png'], '-r300' , '-transparent');


end

function ang = AngleV1V2( V1 , V2 )

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

function theta = tanAngle( vv  )

dy = [ vv(2:end , 2 ) - vv(1:end-1, 2)];
dx = [ vv(2:end , 1 ) - vv(1:end-1, 1)];
theta =[];
for k = 1:length( dy )
    theta = [ theta ; atan2(dy(k), dx(k) ) ];
end
end