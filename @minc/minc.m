classdef minc
    %MINC Class to store volume and HDR information for data to be saved as
    % minc files. Uses NIAK for I/O handling.
    %   
    % Date: June 7th 2014
    % Author: Mathieu Boudreau
    % Contact: mathieu.boudreau2@mail.mcgill.ca
    %
    % Date last modified: June 7th 2014
    %
    %
    %
    % ***
    
    properties
        fileName
    end
    
    properties (Access=private)
        niak_hdr
        niak_volume
        template='@minc/2Dminctemplate.mat';
    end

    
    methods        
        %Constructor
        function obj=minc()
            niak_dat = load('@minc/2Dminctemplate.mat');
            obj.niak_hdr = niak_dat.hdr;
            obj.niak_volume = niak_dat.vol.*0; %Set current niak volume to 0 
            obj.fileName = ['empty_' num2str(round(rand(1)*10000)) '.mnc'];

        end

        %Set Volume
        function obj=setVolume(obj,vol)
            obj.niak_volume=vol;
            
            if length(size(vol)) == 2
                obj.niak_hdr.dimensions=[size(vol) 1];
                obj.niak_hdr.info.dimensions=[size(vol) 1];
            else
                obj.niak_hdr.dimensions=size(vol);
                obj.niak_hdr.info.dimensions=size(vol);                
            end
        end
        %Get Volume
        function rvol=getVolume(obj)
            rvol=obj.niak_volume;
        end
        
        %Set Resolution
        function obj=setResolution(obj,voxel_size)
            obj.niak_hdr.info.voxel_size=voxel_size;             
        end

        %Get Resolution
        function voxRes=getResolution(obj)
            voxRes=obj.niak_hdr.info.voxel_size;             
        end

        %Get Dimensions
        function voxRes=getDimensions(obj)
            voxRes=obj.niak_hdr.info.dimensions;             
        end        
        
        %Set Repitition time
        function obj=setTR(obj,tr)
            display('TR must be set in seconds')
            obj.niak_hdr.details.acquisition.attvalue{getAttributeIndexNiak(obj.niak_hdr,'repetition_time')}=tr;
        end
        
        %Set Flip angle
        function obj=setFA(obj,fa)
            obj.niak_hdr.details.acquisition.attvalue{getAttributeIndexNiak(obj.niak_hdr,'flip_angle')}=fa;
        end
    end
    
end

