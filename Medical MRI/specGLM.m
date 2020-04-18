function w = specGLM(spec,Freqshift,T2,Relax_time,DwellTime,model)
%T2(1) and Freqshift(1) should be WATER


%%% Chemical shift converted to frequency
shift = (4.7 - Freqshift)*100;

Nsamples = length(spec);
Nmetabs = length(Freqshift);
A = zeros(Nsamples,Nmetabs);

switch model
    case 1
        %%%% Model 1: Kronecker Delta at peaks. Shouldnt be used
        
        for I = 2:Nmetabs
            A(round(shift(I)),I) = 1;
        end
        A(1,1) = 1; %Work-around
        
        % Now solve A*weight=spec. These weights would be RELATIVE
        w = (A'*A)^(-1)*A'*spec; %Got very close, but only before using T2
        % After using T2 also got close
        
    case 2
        % Model2: columns of design matrix should be constructed by frequency
        % shifting the unsuppressed water signal
        unsup_water_sig = FID(4.7,Relax_time,Nsamples,DwellTime,1);
        t = linspace(0,Nsamples*DwellTime,Nsamples)';
        FIDs = unsup_water_sig.*exp(1i*2*pi*shift.*t);
        
        % Peaks should be width-scaled, as below. T2(1) should be WATER
        FIDs = FIDs./exp(-t/T2(1)).*exp(-t./T2');
        
        % Then the fIDs should be converted to frequencies
        A = fft(FIDs);
        A = A/Nsamples;
        A = real(A);
        A = [A(end,:);A(1:end-1,:)];
        
        w = (A.'*A)^(-1)*A.'*spec;
        
        
        
    case 3 %Including lipids, baseline, non-suppressed water
        unsup_water_sig = FID(4.7,Relax_time,Nsamples,DwellTime,1);
        t = linspace(0,Nsamples*DwellTime,Nsamples)';
        
        % Include lipids:
        T2(length(T2)+1) = 0.085;
        shift(length(shift)+1) = (4.7-1.3)*100;
        
        FIDs = unsup_water_sig.*exp(1i*2*pi*shift.*t);
        FIDs = FIDs./exp(-t/T2(1)).*exp(-t./T2');
        
        A = fft(FIDs);
        A = A/Nsamples;
        A = real(A);
        A = [A(end,:);A(1:end-1,:)];
        
        %Baseline: constant, linear, quadratic
        A(:,size(A,2)+1) = ones(1,Nsamples);
        A(:,size(A,2)+1) = linspace(0,1,Nsamples);
        A(:,size(A,2)+1) = A(:,size(A,2)).^2;
        
        % Derivatives
        
        %Still need non-suppressed water???
        w = (A.'*A)^(-1)*A.'*spec;
        
end



