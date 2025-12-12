
%path to folder with Ca2+ imaging 

CNMFE_path = '' ; 
cd(CNMFE_path)

%filenames to include in analysis
sessions =  {} ; 

%paramters for Jaccard correlation analysis
inputPeakThreshold = 2.5 ; 
inputMicronsPerPixel = 1.85 ; % micronsPerPixel 2.5 = microns (inscopix), 1 (v3), 1.85 (v4)
inputMaxDist = 500 ; 
inputBinSize = 225 ; 
inputBStart = 50 ;
inputNumBins = 9; %9 for 50um Size
%%

% if analyzing a subset of frames input the paths to the csv files with
% those frame indicies here (csv files are output of the
% "alignEzTrackToMiniscope" jupyter notebooks)
% framesDir= 'F:\\JJM\\miniscope_analysis\\dSPNs\\clustering_analysis\\dSPNs_framesSubsetAnalysis\\' ;
framesDir = 'all_frames';

sizeSessions = size(sessions);

    if  strcmp(framesDir, 'all_frames')
        dirName = strcat('all_frames_', string(datetime('now', 'format', 'y_M_d_HH_mm-ss'), "yyyy-MM-dd-HH-mm-ss"), '_analysisOutput');
        mkdir(dirName);
        dirInput = dirName ;
    else
        dirInput = framesDir ;
    end
%%
for i=1:sizeSessions(1,2)
    session=sessions{1,i} ;
    analyzeJaccardsForSessionFn(session, dirInput, regExp, CNMFE_path, inputPeakThreshold, inputMicronsPerPixel, ...
        inputMaxDist, inputBinSize, inputBStart, inputNumBins) ; 

end


