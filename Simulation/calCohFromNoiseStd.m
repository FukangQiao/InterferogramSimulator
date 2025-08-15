% ===========================================================
% Filename:     calCohFromNoiseStd.m
% Date:   	 	2025-08-14
% Author:    	Fukang Qiao
% Description:  Calculate coherence noiseStd is the standard deviation of noise after multi-look
% The standard deviation after multi-look is equal to that before multi-look divided by sqrt(NL)
% ===========================================================
function coh = calCohFromNoiseStd(noiseStd,NLs,stdIsBeforeMultilook)
% note: noiseStd is the standard deviation of noise after multi-look
% The standard deviation after multi-look is equal to that before multi-look divided by sqrt(NL)

if nargin<1
    help calCohFromNoiseStd;
end
if nargin<2
    NLs = [1 1];
end
if nargin<3
    stdIsBeforeMultilook = 1;
end

NL = NLs(1)*NLs(2);
if stdIsBeforeMultilook
    noiseStd = multilook(noiseStd,NLs(1),NLs(2));
    noiseStd = noiseStd./sqrt(NL);
end

% Referred to the calculation formula in the PPT on the website:
% https://sarimggeodesy.github.io/teaching/SARDataAnalysis_Syllabus
% gamma = 1/Sqrt[2*NL*sigma^2 + 1];

coh = 1./sqrt(1+2*NL*noiseStd.^2);

