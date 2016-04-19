function h = dispGroupOverlay(LCMS, RT, RT_Adjustment, MZ, numRow, numCol,  GCF, LOG)

numSamples = size(LCMS,3);

if ((numRow*numCol) ~= numSamples)
    
    error('Number of samples not compatible with the display setup!')
    
end
    
if LOG == 0

    disp('The two input images are intensity maps. ')

    N = max(log(LCMS(:) + 1));

    for i = 1 : numSamples

        LCMS(:,:,i) = log(LCMS(:,:,i)+1) ./ N;

    end
    
else

    disp('The two input images are log-transformed.')
   
    N = max(LCMS(:));

    for i = 1 : numSamples

        LCMS(:,:,i) = LCMS(:,:,i) ./ N;

    end
 
end

MeanLogLCMS = mean(LCMS, 3);

Fig2Disp = cat(3, MeanLogLCMS, MeanLogLCMS, MeanLogLCMS);

marginWidth = .02;

marginHeight = .033;

imgWidth = (.76 - (numCol-1)*marginWidth)/numCol; % from .05 to .81

imgHeight = (.93 - (numRow-1)*marginHeight)/numRow; % from .05 to .95

figure(GCF)

clf

set(gcf, 'Position', [100, 100, 1280, 720])

h(1) = axes('position', [.83 marginHeight .12 .92]);
plot(RT_Adjustment'*60, RT*60)
%axis([-3.5*60 3.5*60 min(RT)*60 max(RT)*60])

axis([-3.5*60 3.5*60 670 4400])


for i = 1 : numSamples

    iRow = (floor((i-1) / numCol) + 1);
    
    iCol = i - numCol * (iRow - 1);

    Fig2Disp(:,:,2) = LCMS(:,:,i);
        
    h(i+1) = axes('position', [(marginWidth+imgWidth)*(iCol-1)+marginWidth (marginHeight+imgHeight)*(numRow - iRow)+marginHeight imgWidth imgHeight]); 

    imagesc(MZ, RT*60, 1 - Fig2Disp)
       
    grid on
    
    set(gca, 'YLim', [670 4400])
    
    axis off
        
    
end



%xlabel('m/z, Normalised', 'fontsize', 18)
%ylabel('RT, sec', 'fontsize', 18)
%title('LC-MS', 'fontsize', 20)