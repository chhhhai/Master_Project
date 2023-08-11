
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///E:/WinArl35/WinArl35/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (3.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 16/02/23 at 15:01:46", "3.xml#16_02_23at15_01_46"))
	insDoc(aux1, gLnk("R", "Settings", "3.xml#16_02_23at15_01_46_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "3.xml#16_02_23at15_01_46_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "3.xml#16_02_23at15_01_46_pop_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "3.xml#16_02_23at15_01_46_pop_exct_tests"))
