clear all; close all; clc;

%%%%  Project: Blob Statistics

%% === STEP 1 - IMPORTING THE IMAGE ===
% I = imread('IMG_20200317_151602.jpg');  % Fish
% I = imread('IMG_20200317_151618.jpg');  % Lightning + circle
% I = imread('IMG_20200317_151636.jpg');  % screws
I = imread('IMG_20200317_151708.jpg');  % star
J = imresize(I, 0.1);

%% === STEP 2 - RGB>GRAYSCALE CONVERSION ===
R = J(:,:,1);
G = J(:,:,2);
B = J(:,:,3);

grayscale = (0.3*R) + (0.59*G) + (0.1140*B);    % Weighted Method

%% === STEP 3 - IMAGE BINARIZATION ===
T = 155;        % Threshold 

black = (grayscale < T);        % Splitting black & white bits
white = (grayscale >= T);

blobImg = zeros( size(grayscale) );     % Pre-allocation
blobImg(black) = 0;
blobImg(white) = 1;

%% === STEP 4 - CONNECTIVITY ANALYSIS ===
nRows = length(blobImg(:,1));   % Finding total number of rows in image
nCols = length(blobImg(1,:));   % Finding total number of rows in image

% Pre-allocation
r = 1;          
c = 1;
label = -2;     % label to be assigned to each different component

for i = 1:(nRows - 1)
    for j = 1:(nCols - 1)
        window = blobImg(r:r+1,c:c+1);  % 2x2 sqaure matrix for scanning bits
        
        %% CONNECTIVITY RULES
        if ~all( window(1) == window(:) )   % Rule 1
            % Variable declaration
            A = window(2,2);
            B = window(1,2);
            C = window(1,1);
            D = window(2,1);
            
            % Rule 2
            if (C<1) && (B==D) && (A<1)
                A = min(C, A);
                C = A;
                blobImg(r:r+1,c:c+1) = [C B; D A];  % updating values back
            end
            % Rule 3
            if (C==A) && (B<1) && (D<1)
                D = min(B, D);
                blobImg(blobImg==D) = B;    % all m's are equal to n
                
                B = D;
                blobImg(r:r+1,c:c+1) = [C B; D A];
            end
            % Rule 4
            if (C==1) && (B==1) && (D==1) && (A<1) 
                A = label;          % assigining a label
                blobImg(r:r+1,c:c+1) = [C B; D A];
                
                label = label - 1;  % creating a new label
            end
            % Rule 5
            if (B<1) && (D==1) && (A<1)
                A = min(B, A);
                B = A;
                blobImg(r:r+1,c:c+1) = [C B; D A];
            end
            % Rule 6
            if (B==1) && (D<1) && (A<1)
                A = min(D, A);
                blobImg(r:r+1,c:c+1) = [C B; D A];
            end
            % Rule 7
            if (B<1) && (D<1) && (A<1) && (B==D)
                A = min(D, A);
                B = A;
                blobImg(r:r+1,c:c+1) = [C B; D A];
            end
            % Rule 8
            if (B<1) && (D<1) && (A<1) && ~(B==D)
                A = max(B, D);
                blobImg(blobImg==D) = A;        % all m's are equal to n
                blobImg(blobImg==B) = A;        % all m's are equal to n
                
                D = A;
                B = A;
                blobImg(r:r+1,c:c+1) = [C B; D A];
            end
        end
        c = c+1;    % moving onto scanning next column
    end
    c = 1;          % resetting column value for new row
    r = r+1;        % moving onto scanning next row
end


%% === STEP 5 - CENTER OF AREA ===
blobLabels = unique(blobImg(blobImg<1));    % finding no. of blobs
m00 = zeros(size(blobLabels));      % Blob Area pre-allocation

for n = 1: length(blobLabels)
    m00(n) = sum(blobImg(:) == blobLabels(n) ); % finding area of each blob
end

% Pre-allocation
m01 = zeros(size(m00));     
m10 = m01;

m11 = m01;
m20 = m01;
m02 = m01;


% Center of Area
for n = 1:length(blobLabels)       % Sum of Columns
    [rows,cols] = find( blobImg == blobLabels(n) ); % blob coordinates
    
    % Finding variables for area & orientation formulas
    m01(n) = sum( cols );
    m10(n) = sum( rows );
    
    m11(n) = sum( rows.*cols );
    m20(n) = sum( rows.^2 );
    m02(n) = sum( cols.^2 );
end

% Center-point coordinates
i0 = m10./m00;
j0 = m01./m00;

%% === STEP 6 - ORIENTATION ===
% Pre-allocation
theta = zeros(size(m00));   % angle of orientation in degrees
j1 = theta;     % x point on axis of rotation
i1 = theta;     % y point on axis of rotation

for n = 1:length(blobLabels)    % total no. of blobs
    y = ( 2*( m00(n)*m11(n) - m10(n)*m01(n) ) );
    x = ( ( m00(n)*m20(n) - m10(n)^2 )-( m00(n)*m02(n) - m01(n)^2 ) );
    theta(n) = 0.5*atand(y/x);  % orientation mathematical formula
    
    % Conditions for correcting ambiguity in arc-tangent
    if ( y < 0 ) && ( x < 0 )
        theta(n) = theta(n) - 180;
    end
    if ( y > 0 ) && ( x < 0 )
        theta(n) = theta(n) - 180;
    end
    if ( y < 0 ) && ( x > 0 )
        theta(n) = theta(n) - 90;
    end
    if ( y > 0 ) && ( x > 0 )
        theta(n) = theta(n) - 90;
    end
    
    % New points on axis of rotation
    j1(n) = 45*cosd(theta(n)) + j0(n);
    i1(n) = 45*sind(-theta(n)) + i0(n);
end

%% ====== DISPLAYING RESULTS ======
figure(1)
imshow(blobImg);
axis on
hold on;
for n = 1:length(m00)
    plot( j0(n),i0(n) , 'r+','MarkerSize', 10,'LineWidth',1 );
    line( [j0(n) j1(n)], [i0(n) i1(n)], 'Color', 'r', 'LineWidth',0.5 );
end

print('-f1','PlotStar','-dpdf')

% ---------------------------
% END OF CODE
