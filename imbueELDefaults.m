function ELoptions = imbueELDefaults(options, N_feat)
%
%     ELoptions = imbueELDefaults(ELoptions)
%     ELoptions = imbueELDefaults(ELoptions, N_spec)
%
% This function applies default options for evolution learning or for the
% subfunctions used within the evolution learning methodology. Where the
% user provides input options to this function, they are retained, so that
% only unspecified options receive their default options. If the number of
% features is not provided, it is set to 0 (meaningless, but useful for
% specifying options that do not depend on feature interactions or names)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    N_feat = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% SET DEFAULT OPTION VALUES BELOW %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ELoptions.library_orders = 0:2;             % Maximum order for polynomial library
ELoptions.interactions = true(N_feat);      % All species interact by default
ELoptions.symmetric_payoff = true;          % Specifies whether to target learning a symmetric payoff matrix

ELoptions.Xdash_data = [];                  % By default, learn derivatives from data (option exists so user can provide derivative data)             

% Parameters used for function fitting
ELoptions.evolutionary_derivative = true;   % Enforces derivatives summing to zero across components
ELoptions.deriv_method = 'gp';              % Method for derivative estimation (smoothing)
ELoptions.visualise_fit = false;            % Visualise the function fits used for derivative estimation
ELoptions.kernel_bandwidth = -1;            % Local polynomial regression: bandwidth. Set to -1 for "auto" (2.5 times average point separation)
ELoptions.polynomial_order = 4;             % Local polynomial regression: polynomial order
ELoptions.polyfit_regularisation = 1e-7;    % Local polynomial regression: ridge parameter (used to regularise ill-defined regression problems)

% Parameters used for determination of library coefficients
ELoptions.force_positive = true;            % Force coefficients to be positive
ELoptions.shrinkage = 1e-3;                 % Extent of shrinkage to apply (L2 penalisation on non-zero library term coefficients). Set to 0 to apply no shrinkage. 

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