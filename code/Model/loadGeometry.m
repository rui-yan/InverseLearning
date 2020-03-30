function tr = loadGeometry(filename)
    switch filename
        case 'LFP_P3.mat'
            epsilon = 30e-2;
        case 'LFP_P1.mat'
            epsilon = 25e-2;
        case 'LFP_P2.mat'
            epsilon = 10e-2;
        otherwise 
            epsilon = 9e-2;
    end
    data = load(filename);
    [row,col] = find( (data.FP + data.LFP) > epsilon);
    boundaries = boundary(row, col);
    X = col(boundaries);
    Y = row(boundaries);
    pgon = polyshape(X,Y);
    tr = triangulation(pgon);
end