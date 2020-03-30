function tr = loadGeometry(filename)
epsilon = 9e-2;
data = load(filename);
[row,col] = find( data.FP + data.LFP > epsilon);
boundaries = boundary(row, col);
X = col(boundaries);
Y = row(boundaries);
pgon = polyshape(X,Y);
tr = triangulation(pgon);
end