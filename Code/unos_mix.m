function A=unos_mix(name, ncol, nhedln, list_alfa)

% A=unos_mix(name, ncol, nhedln, list_alfa)
%
% Input of data from textual file of any extension.
% The data are of mixed (numeric and alphanumeric type). The file has only
% header lines, not column labels.
% Variables: 'name' = full name to file including extension and path;
%            'ncol' = number of columns in the file;
%            'nhedln' = number of header lines;
%            'list_alfa' = an array with indexes of alphanumeric columns.

a=fopen(name,'r');
format=repmat('f',ncol,1);
for i=1:length(list_alfa)
    format(list_alfa(i))='s';
end
format=[repmat('%',ncol,1) format repmat(' ',ncol,1)]';
format=reshape(format, 1,ncol*3);
A = textscan(a, format, 'headerlines', nhedln);
fclose(a);