clear


%[file_0, pathname,filterindex] = uigetfile({'*.asc'},'Choisir le repertoire à traiter')
selpath = uigetdir
cd(selpath)

Liste=dir("*.asc");
Repertoire_arriver="../Data_corriger_appareil/";
mkdir(Repertoire_arriver)

%figure
%hold on
label=[]

for ifile=1:length(Liste)
    file_0=Liste(ifile).name
    
    Data = fileread(file_0); % remplacement des virgule par des points
    Data = strrep(Data, ',', '.');
    FID = fopen(file_0, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);
    
    sansech = importdata(file_0);

    corr=sansech(:,3)./sansech(:,2);    % fonction de correction
    %ab=-log10(abs(tr))  % absorbance

    nm=sansech(:,1);

    k=10^7./sansech(:,1);               % energie en cm^-1
% 
%     figure
%     title('en nm')
%     plot(nm, corr)
%     grid on

%     figure
%     plot(k,corr)
%     legend(file_0(1:end-4))
%     xlim([12000, 24000])
%     grid on

    BK2_HT=[nm, corr];

    writematrix(BK2_HT,(Repertoire_arriver+file_0(1:end-4)+".csv"), 'Delimiter', ';')

end