
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:///E:/WinArl35/WinArl35/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (2.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 16/02/23 at 14:32:58", "2.xml#16_02_23at14_32_58"))
	insDoc(aux1, gLnk("R", "Settings", "2.xml#16_02_23at14_32_58_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "2.xml#16_02_23at14_32_58_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "2.xml#16_02_23at14_32_58_pop_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "2.xml#16_02_23at14_32_58_pop_exct_tests"))
	aux1 = insFld(foldersTree, gFld("Run of 16/02/23 at 14:36:22", "2.xml#16_02_23at14_36_22"))
	insDoc(aux1, gLnk("R", "Settings", "2.xml#16_02_23at14_36_22_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "2.xml#16_02_23at14_36_22_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "2.xml#16_02_23at14_36_22_pop_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "2.xml#16_02_23at14_36_22_pop_exct_tests"))
	aux1 = insFld(foldersTree, gFld("Run of 16/02/23 at 14:44:01", "2.xml#16_02_23at14_44_01"))
	insDoc(aux1, gLnk("R", "Settings", "2.xml#16_02_23at14_44_01_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "2.xml#16_02_23at14_44_01_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "2.xml#16_02_23at14_44_01_pop_pairw_diff"))
		insDoc(aux2, gLnk("R", "Exact tests", "2.xml#16_02_23at14_44_01_pop_exct_tests"))
