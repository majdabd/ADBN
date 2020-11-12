
   function [w_corr,variance_matrix,std_matrix] = sliding_windows(time_courses,w,w_shape,ss,leonardi_filter)

% This function outputs windowed correlation time series of fMRI time
% courses measured in different Regions of Interest (ROIs).

% Input variables
%%%%%%%%%%%%%%%%%

% time_courses          2-D matrix of size (N,T) where N is the number of ROIs
%                       and T is the number of observations.
% w                     Width of the window used to get windowed correlations, 
%                       default value is 30 (times TR).
% leonardi_filter       Boolean taking 'on' or 'off' (default) value and indicating whether 
%                       the highpass filtering suggested in [1] should be applied or not. 
% w_shape               Indicates the shape of the window used to compute
%                       the sliding correlation. Default values is
%                       'squared' but 'tapered', as suggested in [2] to
%                       reduce spurious correlations is also available.

% Output variable
%%%%%%%%%%%%%%%%% 

% w_corr                3-D matrix of size (N,N,T-w+1) containing T-w+1
%                       square correlation matrices of size (N,N).
% variance_matrix       2-D matrix of size (N,N) encoding the variance of
%                       the fluctuations of the correlation between each 
%                       pair of variables.
% variance_vector       Vector of size (1,N). Element i (0<i<N+1) encodes the 
%                       average level of variance of fluctuations in the 
%                       correlations between ROI i and all other ROIs.
%                       The advantage of this way of encoding variance is
%                       that its elements can be plotted on a brain surface. 

% References
%%%%%%%%%%%%

% [1] N. Leonardi and D. Van de Ville, 2015. On spurious and real fluctuations
%     of dynamic functional connectivity during rest. Neuroimage, Vol 104,
%     pp. 430-436.
% [2] F. Pozzi, T. Di Matteo and T. Aste, 2012. Exponential smoothing
%     weighted correlations, Eur Phys J B, Vol. 85, pp. 1-21.
% [3] A. Zalesky et al., 2014. Time-resolved resting-state brain networks.
%     PNAS, Vol. 111, pp. 10341-10346. 


%% Setting default values, if needed

switch nargin
    case 4
        ss=1;
    case 3
        w_shape         = 'squared';
        ss=1;
    case 2
        w_shape         = 'squared';
        leonardi_filter = 'off';
        ss=1;
    case 1
        w_shape         = 'squared';
        leonardi_filter = 'off';
        w               = 30;
        ss=1;
end

%% Get dimension of data

N = size(time_courses,1); % Number of ROIs
T = size(time_courses,2); % Number of observations


%% Display verification message

disp(['Computing sliding correlation from data containing ' num2str(T) ' observations in ' num2str(N) ' ROIs, using a ' w_shape ' window of width ' num2str(w) ' TRs and filtering set to ''' leonardi_filter '''.'])

%% Filtering data, if required


switch leonardi_filter
    
    case 'on' % Filtering at 1/w as suggested in [1]. Its implementation might be handwavy, it should be cross-checked.
        
        fft_tc      = (fft(time_courses'))';
        k_cut       = round(T/w); % cutoff frequency
        filtered_TC = zeros(N,T);
        
        for k=k_cut:round(T/2) % keep only high frequency components 
            add         = (2*(abs(fft_tc(:,k))*ones(1,T)).*cos(ones(N,1)*(2*pi*(k-1)*[0:(T-1)]/T)+angle(fft_tc(:,k))*ones(1,T)))/T;
            filtered_TC = filtered_TC + add;          
        end
        
        TC=filtered_TC;
        
    case 'off' % No filtering
        
        TC = time_courses;
        
    otherwise
        
        disp('Error: value of third input should be ''on'' or ''off'', verify value and/or order of inputs.')
        return
        
end
%% Build the window shape



switch w_shape
    case 'squared'
        shape = ones(w,1)./w;
        tspan=fix((T-w+1)/ss)-1 ;
        l= length(shape);
        w_corr = zeros(N,N,tspan+1);
    case 'tukey'
        shape = tukeywin(w,1/2);
        shape=shape./sum(shape);
        l= length(shape);
        tspan  =   fix((T-l+1)/ss)-1;
        w_corr =   zeros(N,N,tspan+1);
    case 'exponential' % Implementation as in [3]
        theta = floor(w/3);
        w0    = (1 - exp(-1/theta))/(1-exp(-w/theta));
        shape = w0*exp(((1:w)-w)/theta);
        l=length(shape);
        tspan =  fix((T-w+1)/ss)-1 ;
        w_corr = zeros(N,N,tspan+1);
    case 'gaussian'
        sigma=2;
        f = @(m,s) exp(-0.5*m.^2/s)/sqrt(2*pi*s);
        window_rect = [zeros(1,sigma-1) ones(1,w) zeros(1,sigma-1)];
        arb = 3e3; % An arbitrary large number
        window_temp = conv(window_rect,f(-arb:arb,sigma^2));
        x = 1:(w+2*sigma-2); y = window_temp(x+arb);
        % Resample weights of window to make window(central_sample)=max(window); may be unnecessary
        xx = (1+0.5):((w+2*sigma-2)-0.5); yy = spline(x,y,xx);
        shape = yy/sum(yy);
        l=length(shape);
        tspan =  fix((T-l+1)/ss)-1 ;
        w_corr = zeros(N,N,tspan+1); 
    otherwise
        disp('Error: value of fourth input should be ''squared'' or ''tapered'', verify value and/or order of inputs.')
end
%% Computation of the sliding windows correlation
for t=0:tspan     
            start=t*ss+1;
            stop=start+l-1;
            span=start:stop;
            windowed_series = TC(:,span)';
            windowed_series = bsxfun(@minus, windowed_series, mean(windowed_series, 1));
            windowed_series = bsxfun(@times, windowed_series, 1./sqrt(sum(windowed_series.^2, 1)));
            w_corr(:, :, t+1) = weightedcorrs(windowed_series,shape); 
            w_corr(:,:,t+1)= w_corr(:,:,t+1)-diag(diag(w_corr(:,:,t+1)));
            w_corr(:,:,t+1)= atanh(w_corr(:,:,t+1));
end

%% Computation of variance of fluctuations in correlation
variance_matrix = var(w_corr,[],3);    
std_matrix=     std(w_corr,[],3);   
end       
