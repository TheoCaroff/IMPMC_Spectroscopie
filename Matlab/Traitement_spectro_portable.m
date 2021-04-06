clear

%[file_0, pathname,filterindex] = uigetfile({'*.asc'},'Choisir le repertoire à traiter')
selpath = uigetdir
cd(selpath)
    
Liste=dir("*.txt");
Repertoire_arriver="../Data_trait/";
mkdir(Repertoire_arriver)


for ifile=1:length(Liste)
    file_0=Liste(ifile).name
    
    NEWFILE=Repertoire_arriver+file_0(1:end-4)+".csv"
    
    ENTETE= ['Transmission spectro portable;', newline, 'nm;%T', newline]
    FID = fopen(NEWFILE, 'w');
    fwrite(FID, ENTETE, 'char');
    fclose(FID);
    
    Data = fileread(file_0); % remplacement des virgule par des points
    Data = strrep(Data, ',', '.');
    FID = fopen(file_0, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);
    
    DATA = importdata(file_0);
    nm=DATA.data(:,1);
    Tr=DATA.data(:,2)./100;
    
    DATAOK = [nm, Tr];

    writematrix(DATAOK, NEWFILE, 'Delimiter', ';', 'WriteMode','append')
end

Fichiertemp = tempdir+"repspectro";
FID = fopen(Fichiertemp, 'w')
fwrite(FID, selpath+"\..", 'char')
fclose(FID)
%https://www.developpez.net/forums/d1074908/autres-langages/python/general-python/fonction-renvoie-chemin-repertoire-temporaire/
