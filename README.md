# Muon-Decay
Le code le plus à jour est muon_decay_FFA_argparse_merged.py. Voici quelques instructions sur sont utilisation:

Ce code permet l'analyse de données ainsi que la fusion de fichiers «t_decay» déjà analysés.

## Analyse:

L'analyse nécessite plusieurs arguments essentiels: 
- La date de l'analyse doit être écrite entre guillements après l'argument **--date** ou **-d**. 
- Le chemin du dossier où se trouve les données à analyser sur votre ordinateur doit également être précisé entre guillemets dans l'argument **--folder_analyse** ou **-fa**. 
- Le numéro du scintilateur d'où provient les données doit être précisé grâce à l'argument **--scint**. Pour cet argument, vous avez le choix entre écrire 1 ou 2. 
- Le nom que vous voulez donner au fichier contenant les différences de temps doit être écrit entre guillemets dans l'argument **--tdID_analyse** ou **-tdID_a**. 
    - Si vous ne voulez pas enregistrer les temps de désintégration, écrivez rien pour l'argument **-tdID_analyse** ou **-tdID_a** et écrivez 0 pour l'argument **--save_times** ou **-tdsave**.

Voici un exemple de ce à quoi votre ligne dans votre ligne de commande devrait ressembler (l'ordre des arguments n'a pas d'importance):

**>python muon_decay_FFA_argparse_merged.py -d "07-07" -fa "C:\Users\p123456\Desktop\Donées muons\muon_decay_07-07" --scint 1 -tdID_a "t_decay_07-07"**

Il y a aussi quelques arguments optionnels qui peuvent êtres pertinents: 
- Si vous voulez analyser des données préselectionnées, il faut écrire 1 pour l'argument **--selected_files** ou **-sel_f**. 
- Les premières lettres devant le numéro de l'acquisition dans le nom des données est "acq" automatiquement, mais vous pouvez le changer à l'aide de l'argument **--file_prefixe** ou **-f_p**. 
- Si le nom des données fini autrement que part ".txt", précisez le dans l'argument **--file_ext** ou **-f_e**. 
- Le programme décide, selon le scintillateur utilisé, le seuil vertical au-dessus duquel il cherche des pics, appelé seuil, ainsi que le seuil horizontal par-delà il enregistre des pics, appelé dp_min.  Si vous voulez choisir vous-mêmes ces valeurs, précisez les dans l'argument **--seuil** ou **-s** pour le seuil vertical et **--dp_min** pour le seuil vertical.
   - Les valeurs associées au scintillateur 1 sont:
       - seuil = 0.03;
       - dp_min = 500.
   - Les valeurs associées au scintillateur 2 sont:
       - seuil = 0.045;
       - dp_min = 300.
   - Il est important de remarquer que l'unité de mesure du seuil est en mV, tandis que le dp_min est exprimée en nombre de canaux, c'est-à-dire den nombre de points individuels de l'acquisition.
- Le code enregistre et ferme automatiquement les figures des désintégrations. Si vous voulez voir les figures, l'argument -**-fshow** devrait être suivi de 1 et si vous ne voulez pas enregistrer les figures, l'argument **--fsave** devrait être suivi d'un 0.

## Fusion:

La fusion de fichiers de temps de désintégration peut être faite entre plusieurs fichiers déjà analysés, mais aussi entres des fichiers analysés et un fichier à analyser. Si vous voulez effectuer l'analyse d'une prise de donnée, puis la fusionner avec d'autres fichiers, lisez les paragraphes ci-dessus, puis ajoutez les arguments appropriés pour la section fusion.

La fusion de fichiers requiert certains arguments:
- La date de la fusion doit être écrite entre guillements après l'argument **--date** ou **-d**. 
- Le chemin du dossier où se trouve les fichers des temps de désintégration à fusionner sur votre ordinateur doit également être précisé entre guillemets dans l'argument **--folder_merge** ou **-fm**. 
- Le nom du fichier des temps de désintégration créé dans la fusion doit être précisé dans l'argument **--tdID_merge** ou **-tdID_m**. Il faut l'écrire entre guillemets
    - Si vous ne voulez pas enregistrer les temps de désintégration, écrivez rien pour l'argument **-tdID_merge** ou **-tdID_m** et écrivez 0 pour l'argument **--save_times** ou **-tdsave**.
- Les noms des fichiers à fusionner doivent être écrit sans guillemets et avec un espace les séparant dans l'argument **--t_decays** ou **-tds**. Vous pouvez ajouter autant de fichiers que vous voulez. 
   - Si vous voulez fusionner un fichier qui sera analysé par la section analyse avec d'autres fichiers de temps de désintégration, n'ajoutez *pas le nom du fichier des temps de désintégration qui sera créé dans cet argument*. Écrivez le simplement dans l'argument **--tdID_analyse** ou **-tdID_a**.
   - Il est important d'écrire le nom complet du fichier, avec la terminaison ".txt" ou toute autre terminaison.

Voici un exemple de ce à quoi votre ligne dans votre ligne de commande devrait ressembler (l'ordre des arguments n'a pas d'importance):

**>python muon_decay_FFA_argparse_merged.py -d "07-07" -fm "C:\Users\p123456\Desktop\Donées muons\t_decays" -tdID_m "t_decay_07-07" -tds t_decay_05-07-a.txt t_decay_06-07.txt**
