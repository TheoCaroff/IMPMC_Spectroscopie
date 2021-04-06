clear

selpath = uigetdir
cd(selpath)

[file_0, pathname,filterindex] = uigetfile({'*.asc'},'Choisir référence (ou blanc)')
[file_100, pathname,filterindex] = uigetfile({'*.asc'},'Choisir le fichier avec echantillon')

Repertoire_arriver="../Data_trait/";
mkdir(Repertoire_arriver)

VIRGULE = 1; % mettre 1 si besoin de remplacer les virgules par des points.
BACKGROUND = 1;

if BACKGROUND
    [file_noir, pathname,filterindex] = uigetfile({'*.asc'},'Choisir le fichier noir (ou background)')
end

if VIRGULE
    Data = fileread(file_0); % remplacement des virgule par des points
    Data = strrep(Data, ',', '.');
    FID = fopen(file_0, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);

    Data = fileread(file_100); % remplacement des virgule par des points
    Data = strrep(Data, ',', '.');
    FID = fopen(file_100, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);
    
    if BACKGROUND
        Data = fileread(file_noir); % remplacement des virgule par des points
        Data = strrep(Data, ',', '.');
        FID = fopen(file_noir, 'w');
        fwrite(FID, Data, 'char');
        fclose(FID);
    end
end

sansech = importdata(file_0);
avecech = importdata(file_100);

if BACKGROUND
    noir = importdata(file_noir)
    sansech(:, 2:3)= sansech(:, 2:3) - noir(:, 2:3);
    avecech(:, 2:3) = avecech(:, 2:3) - noir(:, 2:3);
end
    ref=sansech(:,2);


tr=avecech(:,2)./ref;               % transmission
ab=-log10(abs(tr));                 % absorbance

k=10^7./sansech(:,1);               % energie en cm^-1

figure
plot(k, tr)
title('Tr')

figure
plot(k,ab)
title('Absorbance')
%xlim([12000, 24000])
grid on

tr=10.^(-(ab)) %pour comparer avec Perkin, vu que double épaisseur on divse l'abs par 2.

BK2_HT=[sansech(:,1), tr];
if BACKGROUND
        writematrix(BK2_HT,Repertoire_arriver+file_100(1:end-4)+"_SANSREF.csv", 'Delimiter', ';')
        disp('ATTENTION TRAITEMENT SANS PRENDRE EN COMPTE LA LIGNE REFERENCE (cas spectro mono faisceau)')
    else
        writematrix(BK2_HT,Repertoire_arriver+file_100(1:end-4)+"_SANSNOIR_SANSREF.csv", 'Delimiter', ';')
        disp('ATTENTION TRAITEMENT SANS PRENDRE EN COMPTE LA LIGNE REFERENCE (cas spectro mono faisceau)')
end

Fichiertemp = tempdir+"repspectro";
FID = fopen(Fichiertemp, 'w')
fwrite(FID, [pathname filesep '..'], 'char')
fclose(FID)
