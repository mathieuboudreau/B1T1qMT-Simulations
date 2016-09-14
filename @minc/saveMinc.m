function [] = saveMinc(obj)
%SAVEMINC Saves object as a Minc file.
%   Uses NIAK to save into minc file. 
%   Calls custom niak file "niak_write_minc_ss" (modified by MB) so that it
%   handles single slice data correctly. If only interested in 3D/4D
%   volumes, change this function to "niak_write_minc_ss".
%
%   Date: June 7th 2014
%   Author: Mathieu Boudreau
%   Contact: mathieu.boudreau2@mail.mcgill.ca
%
%   Date last modified: June 7th 2014

    %Set filename
    obj.niak_hdr.file_name=obj.fileName;
    
    %Write out minc file
    niak_write_minc_ss(obj.niak_hdr, obj.niak_volume);
end

