format long
fileID = "Mandelbrot_C.txt";
data = load(fileID, '-ascii');
x = data(:,1);
y = data(:,2);
z = data(:,3);
%contourf(data(:,1),data(:,2),data(:,3))
scatter(x,y,30,z,'filled','MarkerFaceAlpha',.75)