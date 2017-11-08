function [im1, im2, tu, tv] = read_flow_file(seqName, iSeq)

% read_image_flow reads in two pairs of frames and possibly the ground truth optical flow
% specific to the file location in Deqing Sun's computer
% 
% isColor~=0 : readin color image if possible; isColor==0: turn color into
% gray if necessary
% 
% called by: read_image_flow_tune_para, compute_flow_test,
%            generate_middelbury, compute_flow_tune_para
% 

% Authors: Deqing Sun, Department of Computer Science, Brown University 
% Contact: dqsun@cs.brown.edu
% $Date: $
% $Revision: $
%
% Copyright 2007-2010, Brown University, Providence, RI. USA
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the Technology Venture Office of Brown University
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.        

filePath = ['data' filesep]; % or the folder that you save the data

switch seqName
        
    case 'middle-eval' 
        
        filePrefex = [filePath 'eval-data/'];
        
        flowFolder = {'Army',  'Mequon', 'Schefflera', 'Wooden',  'Grove', 'Urban', ...
                      'Yosemite',  'Teddy', 'Basketball',  'Evergreen',  'Backyard',  'Dumptruck'};          
        
        imfilename1 = 'frame10.png';
        imfilename2 = 'frame11.png';
        
        im1=double(imread([filePrefex flowFolder{iSeq} filesep imfilename1]));
        im2=double(imread([filePrefex flowFolder{iSeq} filesep imfilename2]));

        % GT not available
        tu = nan(size(im1));     
        tv = tu; 
        
    case 'middle-other'
        % other sequence with GT provided
        
        imgFilePath     = [filePath 'other-data/'];        
        flowFilePath    = [filePath 'other-gt-flow/'];        

        subPath = {'Venus', 'Dimetrodon',   'Hydrangea',    'RubberWhale',...
                    'Grove2', 'Grove3', 'Urban2', 'Urban3', ...
                    'Walking', 'Beanbags',     'DogDance',     'MiniCooper'};
                
        
        img1Filename = [imgFilePath subPath{iSeq} filesep 'frame10.png'];
        img2Filename = [imgFilePath subPath{iSeq} filesep 'frame11.png'];        

        im1=double(imread(img1Filename));
        im2=double(imread(img2Filename));
        
        
        if iSeq <=8
            
            % Read in ground truth flow fields
            flowFilename = [flowFilePath subPath{iSeq} filesep 'flow10.flo'];
            flow = readFlowFile(flowFilename);
            tu = flow(:,:,1);
            tv = flow(:,:,2);
            
            % Set unknown values to nan
            UNKNOWN_FLOW_THRESH = 1e9; 
            tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
            tv (tv>UNKNOWN_FLOW_THRESH) = NaN;

        else
            
            % GT not available
            tu = nan(size(im1));
            tv = tu; 
            
        end;
        
    otherwise
        error('unknown image sequence!');
end;