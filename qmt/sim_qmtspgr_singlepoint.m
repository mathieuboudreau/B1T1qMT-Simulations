function mt_measurement = sim_qmtspgr_singlepoint(protocol,protocol_index,TissueParameters, protocolFlag)

%% Simulate single point MTSPGR measurement
%
% PROJECT MANAGER: Mathieu Boudreau (_MJB)
% DESCRIPTION: 
%
%
% LAST EDITED BY: Mathieu Boudreau (_MJB)
% LAST MODIFIED ON: June 3rd 2014
%

    sled_flips{1} = protocol(protocol_index,1);
    sled_offsets = protocol(protocol_index,2);
    sled_pulse_durations = protocol(protocol_index,3);
    sled_trs = protocol(protocol_index,4);
    sled_excite_flips = protocol(protocol_index,5);
    %[trueParams(:,forIndex) simuParams(:,forIndex) simuError(:,forIndex)] = qmtspgr_simu_fit_error(MeasurementParameters, sled_excite_flips, TissueParameters, T1obs, protocolFlag);
    %forIndex = forIndex+1;
        %% Set Measurement Parameters
    %



    % Pre-allocate memory.
    my_pulse = cell(length(sled_flips),length(sled_offsets),length(sled_trs));

    % Pulse setup
    switch protocolFlag
        case 'custom'
            for forIndex = 1:length(sled_flips)

                temp_sled_flips = sled_flips{forIndex};
                %temp_sled_flips = cell2mat(temp_sled_flips);
                each_dim_sled_flips(forIndex) = length(temp_sled_flips);
                clear temp_sled_flips
            end
        dim_sled_flips = sum(each_dim_sled_flips);

    end

    for expt = 1:length(sled_trs)
      for del = 1:length(sled_offsets)

        switch protocolFlag
            case 'sled'
                fl_max = length(sled_flips);
            case 'custom'
                fl_max = each_dim_sled_flips(expt);
        end

        for fl = 1:fl_max
           switch protocolFlag
               case 'sled' 
                    my_pulse{fl,del,expt} = gaussian_hann(sled_flips(expt,fl),sled_pulse_durations(expt), sled_offsets(del),sled_trs(expt),200); 
               case 'custom' 
                    temp_sled_flips = sled_flips{expt};
                    my_pulse{fl,del,expt} = gaussian_hann(temp_sled_flips(fl),sled_pulse_durations(expt), sled_offsets(del),sled_trs(expt),200); 
                    clear temp_sled_flips;
           end
        end
      end
    end

    %% Set Tissue Parameters
    %
    %	p: 	    input parameters (1 = liquid/lorentzian, 2 = solid/specify)
    %		    [f r t2a t2b r1a r1b]

    load(TissueParameters)

    params_true(1) = m0(2)/m0(1);       % f
    params_true(2) = r;                 % r
    params_true(3) = t2(1);             % t2a
    params_true(4) = t2(2);             % t2b
    params_true(5) = r1(1);             % r1a
    params_true(6) = r1(2);             %r1b
    lineshape_true = 'superlrtz_line';


    %% Create Study object -> Cache Lineshape
    %

    study = struct('template','qmtspgr');
    switch protocolFlag
        case 'sled' 
            study.nominal_angles = reshape(sled_flips',1,4)
            study.nominal_offsets = {sled_offsets sled_offsets sled_offsets sled_offsets}
            study.nominal_offsets(1)
            study.pulse_type = {'gaussian_hann' 'gaussian_hann' 'gaussian_hann' 'gaussian_hann'}
            %study.pulse_type = {'gaussian' 'gaussian' 'gaussian' 'gaussian'};
            %study.pulse_type = {'cwoffres' 'cwoffres' 'cwoffres' 'cwoffres'};
            study.pulse_type

            study.pulse_duration = reshape([sled_pulse_durations' sled_pulse_durations']',1,4)
            study.TR = reshape([sled_trs' sled_trs']',1,4)
            study.flip = [sled_excite_flips(1) sled_excite_flips(1) sled_excite_flips(2) sled_excite_flips(2)]
            study.normalization_scale = [1 1 1 1];
            study.b0_shift = [0 0 0 0];
            study.b1_scale =[1 1 1 1];
        case 'custom'
            study = struct('template','qmtspgr');
            for forIndex = 1:length(sled_flips)

                temp_sled_flips = sled_flips{forIndex};
                %temp_sled_flips = cell2mat(temp_sled_flips);

                if forIndex == 1
                   study.nominal_angles = temp_sled_flips';
                else
                   study.nominal_angles = cat(1,study.nominal_angles,temp_sled_flips');
                end

                clear temp_sled_flips
            end
            dim_sled_flips = sum(dim_sled_flips);
            for dimIndex = 1:dim_sled_flips

                study.nominal_offsets{dimIndex} = sled_offsets
                study.nominal_offsets
                study.pulse_type{dimIndex} = 'gaussian_hann'
                study.pulse_type
                % study.pulse_duration = reshape([sled_pulse_durations' sled_pulse_durations']',1,4)
                % study.TR = reshape([sled_trs' sled_trs']',1,4)
                % study.flip = [sled_excite_flips(1) sled_excite_flips(1)
                % sled_excite_flips(2) sled_excite_flips(2)]
            end
                study.normalization_scale = ones(1,dim_sled_flips);
                study.b0_shift = zeros(1,dim_sled_flips);
                study.b1_scale = ones(1,dim_sled_flips);
    end

    lineshape_prep = cache_lineshape(study, lineshape_true);

    %% Run Simulation
    %

    % Pre-allocate memory.
    switch protocolFlag
        case 'sled' 
            mtspgrMeas = zeros(length(sled_flips), length(sled_offsets), length(sled_trs));
            mtspgrMeas_norm = zeros(length(sled_flips), length(sled_offsets), length(sled_trs));
        case 'custom'
            mtspgrMeas_norm = cell(length(sled_trs),1);
    end

    for expt = 1:length(sled_trs)
        % MT Measruement
        switch protocolFlag
            case 'sled'
                mtspgrMeas(:,:,expt) = mtspgr_rp2_model(params_true, sled_flips(expt,:), sled_offsets, my_pulse(:,:,expt), lineshape_prep, sled_excite_flips(expt), sled_trs(expt)); 
                mtspgrMeas_norm(:,:,expt) = mtspgrMeas(:,:,expt)/max(max(mtspgrMeas(:,:,expt)));
            case 'custom' 
                custom_sled_flips = sled_flips{expt}
                temp_mtspgrMeas_baseline = mtspgr_rp2_model_baseline(params_true, custom_sled_flips, sled_offsets, my_pulse(:,:,expt), lineshape_prep, sled_excite_flips(expt), sled_trs(expt));
                temp_mtspgrMeas = mtspgr_rp2_model(params_true, custom_sled_flips, sled_offsets, my_pulse(:,:,expt), lineshape_prep, sled_excite_flips(expt), sled_trs(expt));

                mtspgrMeas_norm{expt} = {temp_mtspgrMeas/temp_mtspgrMeas_baseline};
        end         


    end

   mt_measurement=mtspgrMeas_norm{1};

