files=dir('*.dat');
nfiles=length(files);
for nfile=1:nfiles
    fprintf('%3i %s \n',nfile,files(nfile).name);
end
nfile=input('Input file number ');
dd=importdata(files(nfile).name);
fprintf('%s \n',...
' YR  MO DY UTSEC   X(KM)       Y(KM)       Z(KM)     VX(KM/S)   VY(KM/S)   VZ(KM/S)  DLAT(DEG) DLON(DEG) DALT(KM)')
fprintf('%4i %2i %2i %2i   %10.4f  %10.4f  %10.4f %10.4f %10.4f %10.4f  %8.4f  %8.4f  %8.4f\n',...
         dd.data(1,1:13))
