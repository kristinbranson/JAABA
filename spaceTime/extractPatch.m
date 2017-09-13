function patch = extractPatch(im,loc1,loc2,theta,psz)

ni = size(im,3);
patch = padgrab(im,0,loc1-psz+1,loc1+psz,loc2-psz+1,loc2+psz,1,ni);
patch = imrotate(patch,theta*180/pi,'bicubic','crop');
patch = padgrab(patch,0,psz/2+1,psz+psz/2,psz/2+1,psz+psz/2,1,ni);