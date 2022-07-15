function nF = computeFFT( fx  , TIME , DATA  , unitLen , nF , Fname , OP )

% DX = unitLen; %*(5-1);
% pt = [ fx(1),  fx(end) ]
% fx
% unitLen
[ r , pt ] = size( DATA );
%pt = round(linspace( fx(1) , fx(end) , 10 ));

out = [];
[ m , n  ] = size( DATA );

% -----------------------------
NP  = 1;
tmp = zeros(NP,1);
% -----------------------------

t = [0:1:length(TIME)-1]';

% pt = [ 11 , 51 , 101 , 140 ];
tmpLeg ={}


for k=1:pt
    
    y         = smoothdata( DATA( :, (k) ));
    [ a , f ] = compute_fft( t , y , NP , nF , k )
    tmp(:)    = k;
    out       = [ out; tmp , f , a ];
    
end
dlmwrite(  [ OP , 'FFT_AmpFreq', sprintf('%s', Fname ) ,'.out'] , out ,  'delimiter' , '\t', 'precision','%.6f' );


h1 = figure(nF),...
    tmpLeg ={};
    for jj = 1:length(pt)
        tmpLeg{jj} = num2str( pt(jj) );
    end
    figure(nF),legend(tmpLeg , 'location','best' ,'fontsize', 10);
    export_fig( h1, [ OP , 'rawfft_', sprintf('%s', Fname) , '.png'], '-r300' , '-transparent');
h2 = figure(nF+1),...
    export_fig( h2, [ OP , 'fft_Mag_', sprintf('%s', Fname) , '.png'], '-r300' , '-transparent');
h3 = figure(nF+2),...
    export_fig( h3, [ OP , 'fft_Power_', sprintf('%s', Fname) , '.png'], '-r300' , '-transparent');
nF = nF + 2;

end

function [ Amp , Freq ] = compute_fft( t, y  , NP , nF , FI  )



% t =[0:0.01:1];
%w = 2*pi*5;
%y = 100.*sin( w*t);
%plot( t , y ),hold on;


figure(nF),...
    plot( t , y , 'linewidth' , 1.5 ),hold on,...
    set( gca , 'fontsize' , 18),hold on,...
 

%% Normalizing the values
dT = ( t(end) -t(1) )/length(t);
Fs = 1/dT;

y1   = fft( y )/length(y);
n    = length(y1);
df   = Fs/n;

ff   = [0:1:(n-1)].*(df);
yval = abs( y1 );

figure(nF + 1 ),...
    subplot(1,2,1),plot( ff(1:n/2) , (yval(1:n/2)) ), hold on,...
    xlabel('Frequency (Hz)'),ylabel('Magnitude'),set( gca , 'fontsize' , 18),...
    subplot(1,2,2),plot( ff(1:n/100) , (yval(1:n/100)) ), hold on,...
    xlabel('Frequency (Hz)'),ylabel('Magnitude'),set( gca , 'fontsize' , 18);

figure(nF + 2 ),...
    subplot(1,2,1),plot( ff(1:n/2) , (yval(1:n/2).^2/df)), hold on,...
    xlabel('Frequency (Hz)'),ylabel('Power'),set( gca , 'fontsize' , 18),...
    subplot(1,2,2),plot( ff(1:n/100) , (yval(1:n/100).^2/df)), hold on,...
    xlabel('Frequency (Hz)'),ylabel('Power'),set( gca , 'fontsize' , 18); % https://in.mathworks.com/matlabcentral/answers/706048-normalization-of-power-spectral-density


[ peak1, loc1 ]=findpeaks(  (yval(1:n/2)), ff(1:n/2) , 'Npeaks' , NP , 'SortStr','descend' );
Amp  = 2.*peak1;
Freq = loc1';
[ Amp , Freq ];

end