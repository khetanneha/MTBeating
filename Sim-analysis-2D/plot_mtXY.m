% Neha Khetan, 2021
function plot_mtXY(  Fh ,  f1 , filSeg , dT , fL , ST , ET , PivLen , filename , opath , Dyn )


% Fiber coors: % Time identity      posX      posY    forceX    forceY   tension
%                   1                3          4       5         6 %      7
tottime = ET - ST;
Timepts = [ ST:dT:ET ];

nSeg = fL/filSeg;
iS  = 1; iE  = nSeg;


V1 = VideoWriter([ opath, 'xyMov', sprintf('%s', filename) , '.avi']); %,  'Uncompressed AVI');
V1.FrameRate = 10;
V1.Quality = 100;
open(V1);


for k=1:length( Timepts )
    
    Key = Timepts( k );
    data = f1(f1( :, 1) == Key , 3:7);
    
    
    fXt  = data( : , 1 );
    fYt  = data( : , 2 );
    
    if k == 1
        minX = abs( min( fXt) );
    end
    
    % Swapping rep. i.e. ( 0, 0 ) is plus end and the coordinate at the end
    % of the pivot
    
    if Dyn
        fXt  = -fXt; fXt = fXt + minX;
        fXt  = fXt( end:-1:1); fYt  = fYt( end:-1:1);
    else
      fXt = fXt + minX;
    end
    
    
    figure(1),...
        plot( fXt , fYt , 'linewidth' ,1.5 ), hold on;
        colormap( jet(length(Timepts) )); 
        c = colorbar,...
        c.Label.String = 'Time (s)',...
        caxis( [ 0 tottime]),...  
        xlabel(' X (\mum)'),...
        ylabel('Y (\mum)'),...
        set( gca , 'fontsize' , 18 ),...
      
    
   
    writeVideo(V1 , getframe( gcf ) );
end
close(V1);


export_fig( gcf, [ opath , 'xy_', sprintf('%s', filename) , '.pdf'], '-r300' , '-transparent');
end




