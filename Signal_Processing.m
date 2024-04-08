classdef Signal_Processing

    properties
        Value {mustBeNumeric}
    end

    % List of methods for SignalProcessing Class
    methods
        function [LFP_PSD, auxF] = calc_PSD(~,LFP, Fs)
            %Fs = Frequency sample
            %TT = epoch size in seconds

            W = Fs; % Window size for pwelch (this can be adjusted empirically untill the desired frequency resolution is reached)
            NFFT=2*max(size(LFP(1,:))); % specifies the number of discrete Fourier transform (DFT)

            % Test the lenght of the resultant matrix
            lfp=LFP(1,:); % TT s epoch
            [Pxx,auxF] = pwelch(lfp,W,0,NFFT,Fs);

            % Pre allocate the PSD result
            LFP_PSD = nan(size(LFP,1),length(Pxx));

            for i=1:size(LFP,1)

                lfp=LFP(i,:); % TT s epoch
                W = Fs; % Window size for pwelch (this can be adjusted empirically untill the desired frequency resolution is reached)
                NFFT=2*max(size(lfp)); % specifies the number of discrete Fourier transform (DFT)

                clear A Pxx F A_norm

                [Pxx,auxF] = pwelch(lfp,W,0,NFFT,Fs);
                F=auxF;

                notch = F>58 & F<62;

                Pxx(notch)=NaN;

                A_norm=Pxx./nansum(Pxx);

                Fidx=find(F<200); %restricts the ploting up to 200 Hz

                LFP_PSD(i,:)=A_norm;
            end
        end

        % Basic pre-processing with blocked data
        function LFP = basic_preprocessing(~,LFP,Fs,band,filter)
            LFP = detrend(LFP);
            if filter
                LFP = eegfilt2(LFP,Fs,band(1),[]);
                LFP = eegfilt2(LFP,Fs,[],band(2));
            end
        end

        function [Call, phiall, confCall, f] = phase_coherence(~,LFP1,LFP2,Fs)

            params.Fs=Fs; % ATENÇÃO!!! sampling frequency
            params.fpass=[0 50]; % band of frequencies to be kept
            params.tapers=[3 5]; % taper parameters
            params.pad=0; % pad factor for fft
            params.err=[2 0.05]; % 1=theoretical; 2=Jacknife
            params.trialave=0;

            win=2;
            nblocks=size(LFP1,1);

            for j=1:nblocks

                data1=LFP1(j,:)';
                data2=LFP2(j,:)';

                [C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencysegc(data1,data2,win,params);

                Call(:,j)=C;
                phiall(:,j)=phi;
                confCall(:,j)=confC;

            end
        end

        % Spectral coherence
        function [Cxy_all, C_band, F2] = spectral_coherence(~,LFP1,LFP2,Fs,band_list,band_names)
            % Parameters
            NFFT=2^12;
            WINDOW= [];
            NOVERLAP=[];

            % Pre-allocate variables
            cxy_teste = mscohere(LFP1(1,:),LFP2(1,:),WINDOW,NOVERLAP,NFFT,Fs);
            Cxy_all = nan(length(cxy_teste),size(LFP1,1));
            % Create a C_band field for every band
            for ii = 1:length(band_names)
                C_band.(band_names{ii}) = nan(1,size(LFP1,1));
            end

            % Main loop for each block
            for i=1:size(LFP1,1)
                % Calculate spectral coherence
                [Cxy,F2]=mscohere(LFP1(i,:),LFP2(i,:),WINDOW,NOVERLAP,NFFT,Fs);

                notch= F2>58 & F2<62;
                Cxy_all(:,i)=Cxy;
                aux_Cxy=Cxy;
                aux_Cxy(notch)=NaN;
                Cxy_all(:,i)=aux_Cxy;
                
                % Get the band value for each epoch
                for ii = 1:length(band_names)
                    f = F2>=band_list(ii,1) & F2<=band_list(ii,2);
                    C_band.(band_names{ii})(i) = nanmean(aux_Cxy(f),1);
                end               

            end

        end

    end
end