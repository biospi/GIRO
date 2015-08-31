function [IBspl, DIBspl] = BsplIntegral(t, x, d, NFactor )

% This is to evaluate the degree d B-spline functions and create a matrix. 
% of basis functions. 
% Input: t is the knots, x is the sampling positions, and d is the degree 
% of the B-spline. 
% Output: B is a numSampleSpl (length of x) -by- numCP (length of t - d - 1)
% matrix with rows being the dimension of the function to represent and the
% column being the dimension of the number of control points.

if size(t,1) < size(t,2)
    
    t = [t t(end)+1:t(end)+d+1]';
    
end

if size(x,1) < size(x,2)
    
    x = x';
    
end

numKnot = length(t);

numSampleSpl = length(x);

Bspl= zeros(numSampleSpl, numKnot-1);  % B splines

for j = 0 : d+1

numCP = numKnot - j - 1; 

for i = 1 : numCP
    
    if j == 0
    
        Bspl( ( x>=t(i) ) & ( x<t(i+j+1) ), i) = 1;
        
    else          
        
        Bspl( (x>=t(i))   & (x<t(i+j)) , i) = Bspl( (x>=t(i))   & (x<t(i+j)), i)        .* (x( x>=t(i) & (x<t(i+j))) - t(i)) / (t(i+j) - t(i));
        
        Bspl( (x>=t(i+1)) & (x<t(i+j+1)), i) = Bspl( (x>=t(i+1)) & (x<t(i+j+1)), i) + Bspl( (x>=t(i+1)) & (x<t(i+j+1)), i+1) .*( t(i+j+1) - x((x>=t(i+1)) & (x<t(i+j+1)))) / (t(i+j+1) - t(i+1));
        
    end
    
end

if j == d

    DBspl= zeros(size(Bspl)); % Differentiation of the B splines

    for i = 1 : numCP-1

        DBspl( (x>=t(i))  & (x<t(i+j+1)) , i) = Bspl( (x>=t(i))   & (x<t(i+j+1)), i)     * (j+1) / (t(i+j+1) - t(i));
        
        DBspl( (x>=t(i+1)) & (x<t(i+j+2)), i) = Bspl( (x>=t(i+1)) & (x<t(i+j+2)), i) + Bspl( (x>=t(i+1)) & (x<t(i+j+2)), i+1) * (-j-1) / (t(i+j+2) - t(i+1));
        
    end
    
end

% if j == d-1
% 
%     HBspl= zeros(size(B)); % Differentiation of the B splines
% 
%     for i = 1 : numCP-1
% 
%         HBspl( (x>=t(i))  & (x<t(i+j+1)) , i) = Bspl( (x>=t(i))   & (x<t(i+j+1)), i)        * (j+1) / (t(i+j+1) - t(i));
%         
%         HBspl( (x>=t(i+1)) & (x<t(i+j+2)), i) = Bspl( (x>=t(i+1)) & (x<t(i+j+2)), i) + Bspl( (x>=t(i+1)) & (x<t(i+j+2)), i+1) * (-j-1) / (t(i+j+2) - t(i+1));
%         
%     end
%     
% end

% Compute the integration of the degree d B-spline.
if (j == d+1)
    
    IBspl= zeros(numSampleSpl, numCP+4); % Integration of the B splines
            
    DIBspl= zeros(numSampleSpl, numCP+4); % Derivative of integration of the B splines
    
    % Construct the coef vector:
    CoefIB = (t((1 : numCP)+j)-t(1 : numCP))/j;
    
    for i = 1 : numCP
    
        IBspl((x>=t(i)) & (x<=t(i+j)),i) = Bspl((x>=t(i)) & (x<=t(i+j)) , i); % This is the look-up-table for the dictionary using potentially non-uniform knots.
        
        DIBspl((x>=t(i)) & (x<=t(i+j)),i) = DBspl((x>=t(i)) & (x<=t(i+j)) , i); 
        
    end
    
    for i = 1 : numCP
        
        IBspl((x>=t(i)) & (x<=t(i+j)),i) = sum(IBspl((x>=t(i)) & (x<=t(i+j)),i:end),2) * CoefIB(i);
        
        IBspl((x>=t(i)) & (x<=t(i+j)),i) = [diff(IBspl((x>=t(i)) & (x<=t(i+j)),i)); 0];
       
        DIBspl((x>=t(i)) & (x<=t(i+j)),i) = sum(DIBspl((x>=t(i)) & (x<=t(i+j)),i:end),2) * CoefIB(i);
        
        DIBspl((x>=t(i)) & (x<=t(i+j)),i) = [diff(DIBspl((x>=t(i)) & (x<=t(i+j)),i)); 0];
       
        
    end
    
    IBspl= IBspl(:,1:numCP+1);
    
    DIBspl= DIBspl(:,1:numCP+1);

end

end 

% NFactor is for normalising the Basis functions for the final warp:
if NFactor == 1 % L1 normalising
    
    for i = 1
        
        N = sum(Bspl(:,i));
        
        Bspl(:,i) = Bspl(:,i)/N;

        DIBspl(:,i) = DIBspl(:,i)/N;
        
%        HBspl(:,i) = HBspl(:,i)/N;
        
    end
        
end

IBspl = sparse(IBspl);
% 
DIBspl = sparse(DIBspl);

% HBspl = sparse(HBspl);

% varargout(3) = HBspl;



