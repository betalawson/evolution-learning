function ELoptions = addDefaultOptions(options, N_feat)
%
%     ELoptions = addDefaultOptions(ELoptions)
%     ELoptions = addDefaultOptions(ELoptions, N_feat)
%
% This function fills out an options struct with default values, wherever
% the values for options are not already provided. If the number of
% features, N_feat, is not provided, it is set to 0 (meaningless, but
% useful for specifying options that do not depend on feature interactions
% or names)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    N_feat = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SET DEFAULT OPTION VALUES BELOW %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ELoptions.library_orders = 0:2;             % Maximum order for polynomial library
ELoptions.interactions = true(N_feat);      % All species interact by default
ELoptions.symmetric_payoff = false;         % Specifies whether to target learning a symmetric payoff matrix

ELoptions.regression_type = 'grad';         % Use 'grad' for gradient matching regression, or 'standard' for traditional least squares regression

ELoptions.Xdash_data = [];                  % By default, learn derivatives from data (option exists so user can provide derivative data)             

% Parameters used for function fitting
ELoptions.deriv_method = 'gp';              % Method for derivative estimation (smoothing)
ELoptions.visualise_fit = false;            % Visualise the function fits used for derivative estimation
ELoptions.kernel_bandwidth = -1;            % Local polynomial regression: bandwidth. Set to -1 for "auto" (2.5 times average point separation)
ELoptions.polynomial_order = 4;             % Local polynomial regression: polynomial order
ELoptions.polyfit_regularisation = 1e-7;    % Local polynomial regression: ridge parameter (used to regularise ill-defined regression problems)

% Parameters used for determination of library coefficients
ELoptions.force_positive = true;            % Force coefficients to be positive
ELoptions.shrinkage = 1e-6;                 % Extent of shrinkage to apply (L2 penalisation on non-zero library term coefficients). Set to 0 to apply no shrinkage. 
ELoptions.nlin_silent = false;              % If 'true', calls to nlinfit inside coefficient determination routines will supress common warnings

% Default names
for k = 1:N_feat
    ELoptions.names{k} = ['x',num2str(k)];
end

% Loop over all provided options and overwrite the defaults
if isstruct(options)
    
    % Get list of provided field names which are possible option overwrites
    provided_names = fieldnames(options);
    
    % Apply each overwrite (or just write useless fields if invalid)
    for k = 1:length(provided_names)
        ELoptions.(provided_names{k}) = options.(provided_names{k});
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ensure shrinkage parameter is valid
if ELoptions.shrinkage < 0
    warning('A negative value for shrinkage/ridge regression loss was provided. Using no shrinkage.');
    ELoptions.shrinkage = 0;
end