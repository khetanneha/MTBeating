% Neha Khetan
% Sep 2021; Extending it for multiple runs
% July 2021: MT beating and filament dynamics


clear; clc; close all;
par1      = 6  % Var
par2      = 10  % No Runs
PivLen    = 2; % 1 um
% -----------------------
filLen        = [ 10 , 10 , 10 , 10 , 10 ] ; % Filament length
% filLen      = [ 15 , 15 , 15 , 15 , 15 ] ; % Filament length
% filLen      = [  5 , 10 , 15 , 20 , 25 ] ; % Filament length
% ---------------------

StartTime = 1;
EndTime   = 1200;
TotTime   = (EndTime - StartTime) ;%1200; % in s
dTQuant   = 1;  % tim res
dTViz     = 30; % for visual
dTQuant2  = 10; % expt comparison
filSeg    = 0.1; % 0.1 um
tim       = [ StartTime:dTQuant:EndTime]';
TR        = par1*par2;

Dyn       = 1;

for ol = 1:par1
    Fname    = ['./analyzed/p', sprintf('%d' , ol) ];
    mkdir( Fname ) ;
end


for ol = 1:par1
    CC = 1;
    opath = ['./analyzed/p', sprintf('%d' , ol ),'/' ]
    EnergY    = []; tipAngVel = []; BW     = [] ; Phi    =[];
    
    for il = ol:par1:TR        
        clf;close all;
        Fname    = [ 'p', sprintf('%d' , ol ) , '_nR_', sprintf('%d', CC ) ];
        fL       = filLen( ol );
        %ipath    = [ './input/test/run', sprintf('%04d', ( il - 1 ) ),'/' ]
        %ipath     = ['/media/nehakhetan/My Passport/Backup_IdeapadS145/Dynein_param/cytosimOutput_DynVel/run', sprintf('%04d', ( il - 1 ) ),'/' ];
         ipath     = ['/media/nehakhetan/My Passport/Backup_IdeapadS145/cytosimOutput_VISCOSITY/run', sprintf('%04d', ( il - 1 ) ),'/' ];
       
        
        
        % Fiber coors: % Time identity      posX      posY    forceX    forceY   tension
        f1          = load( [ ipath , 'fiber_posforce.out' ] );
        plot_mtXY(  1 ,  f1 , filSeg , dTViz , fL , StartTime, EndTime , PivLen , Fname , opath , Dyn );
        
        [ PivVect , PivXY , tipOmega  , tmp , tmp2 ]  = tang_curv2(  2 ,  f1 , filSeg , dTQuant , fL , StartTime, EndTime  , PivLen , Fname , opath );
        tipAngVel  = [ tipAngVel  , tipOmega ];
        BW         = [ BW ; tmp ];
        Phi         = [ Phi , tmp2 ];
        % ====  Fiber energy
        %         f3          = load( [ ipath , 'fiber_energy.out' ] );
        %         EnergY      = [ EnergY , f3(StartTime:EndTime ,3).*(10^-3) ]; % pN    
        CC = CC +1;
    end
    
    % dlmwrite( [ opath , 'CytosimEnergy_', sprintf('%s', Fname) , '.out' ] ,  EnergY , 'delimiter', '\t' );
    dlmwrite( [ opath , 'TipVel_', sprintf('%s', Fname) , '.out'] ,  tipAngVel , 'delimiter', '\t' );
    dlmwrite( [ opath , 'BW_', sprintf('%s', Fname) , '.out'] ,  BW , 'delimiter', '\t' );
    dlmwrite( [ opath , 'PHI_', sprintf('%s', Fname) , '.out'] ,  Phi , 'delimiter', '\t' );
end