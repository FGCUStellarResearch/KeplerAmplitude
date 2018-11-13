%Read in data from McQuillan first
%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Volumes/My Passport for Mac/Stellar Research /kic_mcquillan_xref2.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.


%% Initialize variables.
dirname = "/Users/Tesla/Documents/Derek Research";
filename = "/Users/Tesla/Documents/Derek Research/kic_mcquillan_xref.xlsx";
delimiter = ',';
startRow = 2;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
RAJ2000 = dataArray{:, 1};
DEJ2000 = dataArray{:, 2};
KIC = dataArray{:, 3};
Teff = dataArray{:, 4};
logg = dataArray{:, 5};
Mass = dataArray{:, 6};
Prot = dataArray{:, 7};
e_Prot = dataArray{:, 8};
Rper = dataArray{:, 9};
LPH = dataArray{:, 10};
w = dataArray{:, 11};
dPP = dataArray{:, 12};


%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

sample_length = 27.0;

%figure out start names (EPIC) for each file
ts1_files = dir("/Users/Tesla/Documents/Derek Research/kic_mcquillan_xref.xlsx");
num_files = numel({ts1_files.name});

fname = string;
for ii = 1:num_files
    temp = ts1_files(ii).name;
    temp = strsplit(temp,"-");
    fname(ii) = temp{1};
end
fname = unique(fname);
num_stars = numel(fname);

star_t = zeros(num_stars,2);
star_t = struct('time',[],'flux',[]);
star = zeros(num_stars,2);
star = struct('time',[],'flux',[],'range',[],'teff',[],'prot',[],'mcq_range',[]);
fstar = struct('range',[],'teff',[],'prot',[],'mcq_range',[],'tot_range',[]);

for ii = 1:num_stars
    temp = char(strcat(dirname,fname(ii),'*'));
    starfile = dir(temp);
    for jj = 1:numel(starfile)
        readme = strcat(dirname,starfile(jj).name);
        data = fitsread(readme,'binarytable');
        time = data{1};
        flux = data{8};
        star_t(jj).time = time(~isnan(flux));
        star_t(jj).flux = flux(~isnan(flux));
        star_t(jj).flux = star_t(jj).flux/median(star_t(jj).flux);
    end
    star(ii).time = vertcat(star_t(:).time);
    star(ii).flux = vertcat(star_t(:).flux);
    tot_range(ii) = prctile(star(ii).flux,95)-prctile(star(ii).flux,5);
    tot_range(ii) = 1e3*tot_range(ii);
    %next up is to sample in 27-day samples...start with 10 samples
    for kk = 1:100
        time0 = star(ii).time(1);
        timeE = star(ii).time(end)-sample_length;
        %randomly select start time
        start_time = rand(1)*(timeE-time0)+time0;
        %find index closest to this
        [val, start_index] = min(abs(star(ii).time-start_time));
        [val, end_index] = min(abs(star(ii).time-(star(ii).time(start_index)+sample_length)));
        %now calculate range within those values
        star(ii).range(kk) = prctile(star(ii).flux(start_index:end_index),95)-prctile(star(ii).flux(start_index:end_index),5);
        star(ii).range(kk) = 1e3*star(ii).range(kk);
    end
    S = fname(ii);
    S = char(S);
    S = str2num(S(5:end)); %now have KIC number for star ii
    mcq_index = find(KIC == S);
    if isempty(mcq_index)
        star(ii).teff = -99;
        star(ii).prot = -99;
        star(ii).mcq_range = -99;
    else
        star(ii).teff = Teff(mcq_index);
        star(ii).prot = Prot(mcq_index);
        star(ii).mcq_range = Rper(mcq_index);
    end
end
temp = struct2cell(star);
nums = numel(temp(1,:));
%clean up if stars don't match
jj=1;
for ii=1:nums
    if star(ii).teff > 0
        fstar(jj).range = star(ii).range;
        fstar(jj).teff = star(ii).teff;
        fstar(jj).prot = star(ii).prot;
        fstar(jj).mcq_range = star(ii).mcq_range;
        fstar(jj).tot_range = tot_range(ii);
        jj = jj+1;
    end
end

%Try making sample plot showing overall range and 100 samples
num_files = jj-1;

for ii=1:num_files
    full_range(ii) = fstar(ii).tot_range;
    period(ii) = fstar(ii).prot(1);
    rrange(ii,1:100) = fstar(ii).range;
    tstar(ii) = fstar(ii).teff(1);
    rel_range(ii,1:100) = (fstar(ii).range-fstar(ii).tot_range)/fstar(ii).tot_range;
end

figure(1)
%hold off
semilogy(period,full_range,'or','markersize',3,'markerfacecolor','r')
hold on

for ii=1:num_files
    semilogy(period(ii),rrange(ii,:),'.k','markersize',1)
    drange(ii) = std(rrange(ii,:));
end
xlabel('Period (d)')
ylabel('Photometric Range (ppt)')
hold on
semilogy(period,full_range,'or','markersize',3,'markerfacecolor','r')
axis([0 45 0.1 1e3])
saveas(gcf,"full90.png")
%saveas(gcf,"full90.fig")
hold off


figure(2)
clear drange
idx = find(tstar<4500);
snumfiles = numel(idx);
speriod = period(tstar<4500);
sfull_range = full_range(tstar<4500);

for ii=1:snumfiles
    srrange(ii,1:100) = rrange(ii,:);
end

semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
hold on

for ii=1:snumfiles
    semilogy(speriod(ii),srrange(ii,:),'.k','markersize',1)
    drange(ii) = std(srrange(ii,:));
end
xlabel('Period (d)')
ylabel('Photometric Range (ppt)')
hold on
semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
axis([0 45 0.1 1e3])
saveas(gcf,"t1_90.png")
%saveas(gcf,"t1_90.fig")
hold off

figure(3)
clear drange
idx = find(tstar>4500 & tstar<5500);
snumfiles = numel(idx);
speriod = period(tstar>4500 & tstar<5500);
sfull_range = full_range(tstar>4500 & tstar<5500);

for ii=1:snumfiles
    srrange(ii,1:100) = rrange(ii,:);
end

semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
hold on

for ii=1:snumfiles
    semilogy(speriod(ii),srrange(ii,:),'.k','markersize',1)
    drange(ii) = median(srrange(ii,:));
end
xlabel('Period (d)')
ylabel('Photometric Range (ppt)')
hold on
semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
%semilogy(speriod,drange,'ob','markersize',4,'markerfacecolor','b')
axis([0 45 0.1 1e3])
saveas(gcf,"t2_90.png")
%saveas(gcf,"t2_90.fig")
hold off


figure(4)
clear drange
idx = find(tstar>5500);
snumfiles = numel(idx);
speriod = period(tstar>5500);
sfull_range = full_range(tstar>5500);

for ii=1:snumfiles
    srrange(ii,1:100) = rrange(ii,:);
end

semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
hold on

for ii=1:snumfiles
    semilogy(speriod(ii),srrange(ii,:),'.k','markersize',1)
    drange(ii) = median(srrange(ii,:));
end
xlabel('Period (d)')
ylabel('Photometric Range (ppt)')
hold on
semilogy(speriod,sfull_range,'or','markersize',3,'markerfacecolor','r')
%semilogy(speriod,drange,'xb','markersize',4,'markerfacecolor','b')
axis([0 45 0.1 1e3])
saveas(gcf,"t3_90.png")
%saveas(gcf,"t3_90.fig")
hold off

