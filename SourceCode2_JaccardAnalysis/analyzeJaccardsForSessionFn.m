function [outputFileArray, frameFileNames] = analyzeJaccardsForSessionFn(session, framesDir, regExp, CNMFE_path, ...
    inputPeakThreshold, inputMicronsPerPixel, inputMaxDist, inputBinSize, inputBStart, inputNumBins)

% input are analysis parameters, list of CNMFE files and parameters of # of
% velocity bins to analyze if performing a analysis on bin subsets

    if ~contains(framesDir, 'all_frames') 
        cd (framesDir) ; 
        mkdir('analysisOutput'); 
        frameFiles = dir(strcat(session, regExp)); 
        frameFileNames = vertcat(frameFiles.name);
        sizeFrameFileNames = size(frameFileNames);

        outputSessionArrayJaccards = {} ;
        outputSessionArraysignalPeaks = {} ;
        outputFilesPerVBin = cell(sizeFrameFileNames(1,1), 1);

        for i=1:sizeFrameFileNames(1,1)

            disp('loading session');
            disp(frameFileNames(i,:)); 
            outputFileArray_vbin = jaccardComputeJjmBinInputFrameSubsetFn(session, ...
                frameFileNames(i,:), framesDir, CNMFE_path, inputPeakThreshold, ...
                inputMicronsPerPixel, inputMaxDist, inputBinSize, inputBStart, inputNumBins) ; 

            
            outputFilesPerVBin{i, 1}=outputFileArray_vbin;
        end
        outputFileArray=outputFilesPerVBin;


    else 
        frameFileNames = {};
        %create save path here
        
        outputFileArray = jaccardComputeJjmBinInputFrameSubsetFn(session, ...
                'all_frames', framesDir, CNMFE_path, inputPeakThreshold, ...
                inputMicronsPerPixel, inputMaxDist, inputBinSize, inputBStart, inputNumBins); 

    end

end