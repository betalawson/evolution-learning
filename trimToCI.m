function x = trimToCI(x,p)
%
%    x = trimToCI(x,p)
%
% This function takes the input data 'x', which is a three-dimensional
% array intended for use in a boxplot, and looks over the "data replicates"
% in the third dimension, setting all extreme values (those falling outside
% of the central proportion 'p' of the data), to a value of NaN.

switch ndims(x)

    case 2
        
        % Extract the sizes of the data
        [Nx,Ny] = size(x);
        
        % Loop over each column
        for j = 1:Ny
            
            % Extract the data here
            dat = x(:,j);
            
            % Sort the data
            dat = sort(dat);
            
            % Trim out the most extreme values
            start_val = floor((0.5-p/2)*Nx + 1);
            end_val = ceil((0.5+p/2)*Nx);
            dat(1:start_val-1) = NaN;
            dat(end_val+1:end) = NaN;
            
            % Re-store the data in the vector
            x(:,j) = dat;
                
        end
            
    
    case 3

        % Extract the sizes of the data
        [Nx,Ny,Nz] = size(x);
        
        % Loop over each separate box object to be drawn
        for i = 1:Nx
            for j = 1:Ny
                
                % Extract the data here
                dat = squeeze(x(i,j,:));
                
                % Sort the data
                dat = sort(dat);
                
                % Trim out the most extreme values
                start_val = floor((0.5-p/2)*Nz + 1);
                end_val = ceil((0.5+p/2)*Nz);
                dat(1:start_val-1) = NaN;
                dat(end_val+1:end) = NaN;
                
                % Re-store the data in the vector
                x(i,j,:) = dat;
                
            end
        end
    
end