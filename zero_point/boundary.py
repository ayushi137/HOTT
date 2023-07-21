##############################################################################
'''
This script is used to store the values needed to create the mask using
'creating_csr.py' and boundary values to make the plots look presentable in 
'zero_point_correction.py'  

This script doesn't require any module. However, it is used in other scipts:
creating_csr.py
zero_point_correction.py

To run: 
should be called appropriatly by other scripts

By: Ayushi Singh
Last Modified: 16 March 2022
'''
##############################################################################


# purpose = grid, edge

def boundaries(name, purpose ):
	if name == 'OrionA_C':
		if purpose == 'edge':
			top =  		[1341, 3026]
			left =		[88, 2175]
			bottom = 	[1592, 402]
			right = 	[2810, 1280]

		if purpose == 'grid':
			bottom = 190
			top = 290
			left = 170
			right = 290		

	elif name == 'OrionA_N':
		if purpose == 'edge':
			top =  		[2395, 4474]
			left =		[561, 1910]
			bottom = 	[2844, 448]
			right = 	[4677, 3051]

		if purpose == 'grid':	
			bottom = 230
			top = 390
			left = 190
			right = 350		

	elif name == 'OrionA_S':
		if purpose == 'edge':
			top =  		[1699, 3930]
			left =		[127, 2392]
			bottom = 	[2439, 415]
			right = 	[3942, 2034]

		if purpose == 'grid':	
			bottom = 100
			top = 230
			left = 110
			right = 270			

	elif name == 'OrionB_N':
		if purpose == 'edge':
			top    =  	[1760, 3040]
			left   =    [560, 1930]
			bottom = 	[2030, 400]
			right  = 	[3330, 1470]

		if purpose == 'grid':
			bottom = 470
			top = 580
			left = 70
			right = 200
			
	elif name == 'OrionB_NN':
		if purpose == 'edge':
			top    =  	[1770, 3580]
			left   =    [100, 2000]
			bottom = 	[1900, 430]
			right  = 	[3490, 2000]

		if purpose == 'grid':
			bottom = 520
			top = 670
			left = 10
			right = 160
	
	elif name == 'OrionB_S':
		if purpose == 'edge':
			top    =  	[2600, 4880]
			left   =    [470, 3060]
			bottom = 	[3100, 440]
			right  = 	[5190, 2190]

		if purpose == 'grid':
			bottom = 330
			top = 510
			left = 80
			right = 270

	elif name == 'CepL1157':
		if purpose == 'edge':
			top =  		[1056, 1736]
			left =		[378, 1040]
			bottom =	[1229, 82]
			right = 	[1913, 802]

		if purpose == 'grid':
			bottom = 85
			top = 160
			left = 378
			right = 453

	elif name == 'CepL1172':
		if purpose == 'edge':
			top =  		[776, 1794]
			left =		[527, 429]
			bottom = 	[1758, 187]
			right = 	[1970, 1581]

		if purpose == 'grid':
			bottom = 100
			top = 170
			left = 310
			right = 375

	elif name == 'CepL1228':
		if purpose == 'edge':
			top =   	[315, 1928]
			left = 		[293, 560]
			bottom =  	[1454, 462]
			right =   	[1492, 1809]

		if purpose == 'grid':
			bottom = 450
			top = 510
			left = 330
			right = 380	

	elif name == 'CepL1241':
		if purpose == 'edge':
			top =  		[1977, 2149]
			left =		[515, 2050]
			bottom = 	[685, 283]
			right = 	[2183, 333]

		if purpose == 'grid':
			bottom = 420
			top = 510
			left = 190
			right = 280		
		
	elif name == 'CepL1251':
		if purpose == 'edge':
			top =  		[2054, 1347]
			left =		[523, 1097]
			bottom = 	[700, 317]
			right = 	[2256, 555]

		if purpose == 'grid':
			bottom = 400
			top = 470
			left = 125
			right = 180		

	elif name == 'IC5146':
		if purpose == 'grid':
			bottom = 30
			top = 130
			left = None
			right = None	

		if purpose == 'edge':
			top = 		[2338, 2311] 	
			left =		[186, 1446]
			bottom = 	[656, 418]
			right = 	[2645, 1242]

	elif name == 'Perseus_W':
		if purpose == 'grid':
			top = 140
			bottom = 12
			left = 270
			right = 415
			
		if purpose == 'edge':
			top =  		[2700, 4040]
			left =		[600, 2750]
			bottom =	[2070, 400]
			right = 	[4280, 1670]	


	elif name == 'Perseus_E':
		if purpose == 'grid':
			bottom = 40
			top = 180
			left = 180
			right = 310

		if purpose == 'edge':		
			top =  		[1940, 4140]
			left =		[510, 3320]
			bottom = 	[2420, 400]
			right = 	[3940, 1200]

	
	elif name == 'Tau_B18':
		if purpose == 'grid':
			bottom = 80
			top = 170
			left = 280
			right = 380

		if purpose == 'edge':
			top =  		[2710, 3110]
			left =		[600, 1570]
			bottom = 	[1400, 460]
			right = 	[3610, 2020]

	elif name == 'OphL1688':
		if purpose == 'grid':
			bottom = 20
			top = 190
			left = 120
			right = 300	

		if purpose == 'edge':		
			top    =  	[2700, 5090]
			left   =    [540, 2070]
			bottom = 	[2970, 500]
			right  = 	[5180, 3450]

	elif name == 'OphL1712':
		if purpose == 'grid':
			bottom = 70
			top = 140
			left = 60
			right = 140	

		if purpose == 'edge':		
			top    = [1000, 1880]
			left   = [520, 1280]
			bottom = [1775, 397]
			right  = [2275, 1004]

	elif name == 'TauFill':
		if purpose == 'grid':
			bottom = 250
			top = 340
			left = 320
			right = 430	

		if purpose == 'edge':		
			top    =  	[2690, 2430]
			left   =    [540, 850]
			bottom = 	[900, 400]
			right  = 	[3100, 1950]
			
	elif name == 'TauL1489':
		if purpose == 'grid':
			bottom = 310
			top = 370
			left = 580
			right = 630	

		if purpose == 'edge':		
			top    =  	[600, 1230]
			left   =    [100, 900]
			bottom = 	[516, 400]
			right  = 	[1000, 745]

	elif name == 'TauL1517':
		if purpose == 'grid':
			bottom = 450
			top = 520
			left = 180
			right = 270	

		if purpose == 'edge':		
			top    =  	[1373, 1976]
			left   =    [113, 1105]
			bottom = 	[863, 409]
			right  = 	[1954, 1247]
			

	elif name == 'TauL1521':
		if purpose == 'grid':
			bottom = 300
			top = 460
			left = 380
			right = 590	

		if purpose == 'edge':		
			top    = [4728, 4957]
			left   = [725, 2117]
			bottom = [1861, 570]
			right  = [6049, 3438]

	elif name == 'TauL1539':
		if purpose == 'grid':
			bottom = 300
			top = 380
			left = 180
			right = 270	

		if purpose == 'edge':		
			top    =  	[1167, 2205]
			left   =    [87, 1425]
			bottom = 	[993, 376]
			right  = 	[2090, 1155]
	
	elif name == 'TauL1544':
		if purpose == 'grid':
			bottom = 280
			top = 380
			left = 100
			right = 210

		if purpose == 'edge':		
			top    =  	[2319, 2791]
			left   =    [178, 1191]
			bottom = 	[795, 467]
			right  = 	[2959, 2113]

	elif name == 'TauL1551':
		if purpose == 'grid':
			bottom = 20
			top = 100
			left = 370
			right = 450	

		if purpose == 'edge':		
			top    =  	[1761, 2055]
			left   =    [66, 997]
			bottom = 	[504, 436]
			right  = 	[2188, 1488]			

	elif name == 'TauS1':
		if purpose == 'grid':
			bottom = 160
			top = 250
			left = 320
			right = 410	

		if purpose == 'edge':		
			top    =  	[1758, 2405]
			left   =    [716, 1575]
			bottom = 	[1614, 402]
			right  = 	[2664, 1287]

	elif name == 'TauS2':
		if purpose == 'grid':
			bottom = 200
			top = 310
			left = 330
			right = 460	

		if purpose == 'edge':		
			top    =  	[2726, 3116]
			left   =    [564, 1561]
			bottom = 	[1402, 467]
			right  = 	[3597, 2014]
	
	elif name == 'TauS3':
		if purpose == 'grid':
			bottom = 250
			top = 330
			left = 430
			right = 520	

		if purpose == 'edge':		
			top    =  	[1645, 2053]
			left   =    [98, 1113]
			bottom = 	[583, 459]
			right  = 	[2118, 1434]

	elif name == 'TauT3':
		if purpose == 'grid':
			bottom = 240
			top = 340
			left = 500
			right = 600	

		if purpose == 'edge':		
			top    =  	[1968, 2476]
			left   =    [668, 1424]
			bottom = 	[1519, 292]
			right  = 	[2858, 1324]
			
	elif name == 'TauTMC':
		if purpose == 'grid':
			bottom = 250
			top = 390
			left = 280
			right = 420	

		if purpose == 'edge':		
			top    =  	[3023, 4075]
			left   =    [265, 2058]
			bottom = 	[1388, 553]
			right  = 	[4146, 2646]
			
	elif name == 'TauTMC_E':
		if purpose == 'grid':
			bottom = 250
			top = 330
			left = 220
			right = 310

		if purpose == 'edge':		
			top    =  	[1051, 2284]
			left   =    [73, 1596]
			bottom = 	[1070, 407]
			right =    [2024,1088] 																											
	elif name == 'CygnusX':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1238, 1435]
			left   = [437, 991]
			bottom = [952, 253]
			right  = [1707, 701]

	elif name == 'M16':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1316, 1537]
			left   = [472, 726]
			bottom = [915, 308]
			right  = [1755, 1119]

	elif name == 'M16_2':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':	
			top    = [1900, 3507]
			left   = [50, 2100]
			bottom = [1665, 384]
			right  = [3455, 1830]
	
	elif name == 'M17':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1222, 2400]
			left   = [467, 1680]
			bottom = [1855, 328]
			right  = [2600, 1030]

	elif name == 'M17_2':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1882, 3528]
			left   = [44, 2064]
			bottom = [1695, 375]
			right  = [3430, 1821]

	elif name == 'NGC2264':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [917, 1780]
			left   = [44, 887]
			bottom = [704, 294]
			right  = [1573, 1206]

	elif name == 'NCG7538':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1238, 1435]
			left   = [450, 1000]
			bottom = [939, 244]
			right  = [1703, 697]

	elif name == 'W3':
		if purpose == 'grid':
			bottom = 50
			top = 120
			left = 50
			right = 130	

		if purpose == 'edge':		
			top    = [2119, 1829]
			left   = [438, 1882]
			bottom = [506, 174]
			right  = [2172, 153]

	elif name == 'W3_2':
		if purpose == 'grid':
			bottom = 40
			top = 150
			left = None
			right = 130	

		if purpose == 'edge':		
			top    = [2320, 3322]
			left   = [70, 2698]
			bottom = [888, 528]
			right  = [3022, 1210]

	elif name == 'Spider':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [435, 4825]
			left   = [70, 517]
			bottom = [4403, 442]
			right  = [4806, 4699]

	elif name == 'W48':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1886, 3152]
			left   = [470, 1712]
			bottom = [2038, 368]
			right  = [3390, 1784]

	elif name == 'CygX-N':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [641, 2993]
			left   = [509, 357]
			bottom = [3114, 110]
			right  = [3292, 2713]
			
	elif name == 'CygX-S':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1926, 3365]
			left   = [91, 2893]
			bottom = [794, 447]
			right  = [2653, 910]

	elif name == 'MonOB1':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [699, 1462]
			left   = [24, 868]
			bottom = [590, 300]
			right  = [1253, 904]
	
	elif name == 'MonR2':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [750, 1481]
			left   = [40, 1000]
			bottom = [560, 333]
			right  = [1240, 814]

	elif name == 'Rosette':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [900, 2010]
			left   = [46, 866]
			bottom = [945, 294]
			right  = [1800, 1416]

	elif name == 'Vela':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1810, 3610]
			left   = [500, 830]
			bottom = [1500, 420]
			right  = [2828, 3200]

	elif name == 'Aquila':
		if purpose == 'grid':
			bottom = 80
			top = 260
			left = 90
			right = 270	

		if purpose == 'edge':
			top    = [3370, 5122]
			left   = [467, 3322]
			bottom = [2563, 406]
			right  = [5438, 2166]

	elif name == 'Ser_Main':
		if purpose == 'grid':
			bottom = 190
			top = 350
			left = 20
			right = 180	

		if purpose == 'edge':		
			top    = [1866, 4612]
			left   = [486, 3736]
			bottom = [2983, 362]
			right  = [4399, 1263]

	elif name == 'Serpens':
		if purpose == 'grid':
			bottom = 230
			top = 340
			left = 130
			right = 240	

		if purpose == 'edge':		
			top    = [2234, 3190]
			left   = [470, 2043]
			bottom = [1764, 392]
			right  = [3546, 1522]

	elif name == 'Aquila_W':
		if purpose == 'grid':
			bottom = 30
			top = 180
			left = 320
			right = 500	

		if purpose == 'edge':		
			top    = [1607, 4013]
			left   = [516, 3018]
			bottom = [3212, 366]
			right  = [4335, 1382]

	elif name == 'Cha_I':
		if purpose == 'grid':
			bottom = 110
			top = 240
			left = 240
			right = None	

		if purpose == 'edge':		
			top    = [2755, 3417]
			left   = [387, 3013]
			bottom = [589, 139]
			right  = [3012, 497]

	elif name == 'Cha_II':
		if purpose == 'grid':
			bottom = 120
			top = 230
			left = 20
			right = 120	

		if purpose == 'edge':		
			top    = [1040, 3160]
			left   = [95, 946]
			bottom = [1967, 343]
			right  = [2859, 2530]

	elif name == 'Cha_III':
		if purpose == 'grid':
			bottom = 50
			top = 150
			left = 50
			right = 180	

		if purpose == 'edge':		
			top    = [875, 3320]
			left   = [87, 1040]
			bottom = [2831, 372]
			right  = [3582, 2615]

	elif name == 'Csack_Glob1':
		if purpose == 'grid':
			bottom = 25
			top = 70
			left = 140
			right = 190	

		if purpose == 'edge':		
			top    = [488, 1146]
			left   = [111, 634]
			bottom = [669, 308]
			right  = [1013, 823]

	elif name == 'Csack_Glob2':
		if purpose == 'grid':
			bottom = 25
			top = 70
			left = 150
			right = 180	

		if purpose == 'edge':		
			top    = [284, 938]
			left   = [160, 338]
			bottom = [798, 262]
			right  = [908, 866]

	elif name == 'Coalsack':
		if purpose == 'grid':
			bottom = 100
			top = 160
			left = None
			right = 100	

		if purpose == 'edge':		
			top    = [576, 1542]
			left   = [550, 463]
			bottom = [2748, 575]
			right  = [2788, 1688]

	elif name == 'CrA_N':
		if purpose == 'grid':
			bottom = 30
			top = 160
			left = 90
			right = 230	

		if purpose == 'edge':		
			top    = [2672, 3864]
			left   = [526, 2025]
			bottom = [2077, 392]
			right  = [4213, 2190]

	elif name == 'CrA_S':
		if purpose == 'grid':
			bottom = 10
			top = 150
			left = None
			right = 140	

		if purpose == 'edge':		
			top    = [2108, 4071]
			left   = [474, 2832]
			bottom = [2620, 376]
			right  = [4212, 1572]

	elif name == 'Lupus_I':
		if purpose == 'grid':
			bottom = 300
			top = 430
			left = 200
			right = 320	

		if purpose == 'edge':		
			top    = [1747, 3517]
			left   = [79, 1642]
			bottom = [1800, 358]
			right  = [3434, 2270]

	elif name == 'Lupus_III':
		if purpose == 'grid':
			bottom = 140
			top = 220
			left = 20
			right = 100	

		if purpose == 'edge':		
			top    = [941, 2209]
			left   = [84, 1377]
			bottom = [1262, 359]
			right  = [2064, 1204]

	elif name == 'Lupus_IV-SP1':
		if purpose == 'grid':
			bottom = 30
			top = 110
			left = 60
			right = 180	

		if purpose == 'edge':		
			top    = [1233, 2522]
			left   = [451, 1628]
			bottom = [2211, 385]
			right  = [2986, 1300]

	elif name == 'Lupus_IV-SP2':
		if purpose == 'grid':
			bottom = 50
			top = 120
			left = 30
			right = 100	

		if purpose == 'edge':		
			top    = [1150, 1779]
			left   = [504, 1052]
			bottom = [1390, 341]
			right  = [2067, 1088]

	elif name == 'Musca':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1277, 3486]
			left   = [435, 2950]
			bottom = [2005, 107]
			right  = [2866, 691]
		
	elif name == 'North_Streamer':
		if purpose == 'grid':
			bottom = 110
			top = 255
			left = 10
			right = 200	

		if purpose == 'edge':		
			top    = [1461, 4558]
			left   = [135, 2927]
			bottom = [3724, 433]
			right  = [5050, 2027]

	elif name == 'Pipe_C':
		if purpose == 'grid':
			bottom = 90	
			top = 150
			left = 130
			right = 200	

		if purpose == 'edge':		
			top    = [1157, 1826]
			left   = [453, 1101]
			bottom = [1397, 350]
			right  = [2091, 1085]

	elif name == 'Pipe_E':
		if purpose == 'grid':
			bottom = 100
			top = 220
			left = 20
			right = 140	

		if purpose == 'edge':		
			top    = [1064, 3551]
			left   = [125, 2602]
			bottom = [2792, 381]
			right  = [3674, 1349]

	elif name == 'B59':
		if purpose == 'grid':
			bottom = 70
			top = 150
			left = 180
			right = 270	

		if purpose == 'edge':		
			top    = [1365, 2201]
			left   = [490, 1215]
			bottom = [1544, 340]
			right  = [2475, 1338]

	elif name == 'B68':
		if purpose == 'grid':
			bottom = 200
			top = 270
			left = 110
			right = 190	

		if purpose == 'edge':		
			top    = [1389, 2202]
			left   = [501, 1258]
			bottom = [1555, 346]
			right  = [2492, 1326]

	elif name == 'Pipe_fill_1':
		if purpose == 'grid':
			bottom = 30
			top = 140
			left = 60
			right = 190	

		if purpose == 'edge':		
			top    = [1754, 3510]
			left   = [60, 1920]
			bottom = [1838, 385]
			right  = [3541, 2004]
	
	elif name == 'Pipe_fill_2':
		if purpose == 'grid':
			bottom = 100
			top = 160
			left = 110
			right = 160	

		if purpose == 'edge':		
			top    = [1286, 1633]
			left   = [455, 779]
			bottom = [974, 345]
			right  = [1795, 1232]

	elif name == 'Pipe_fill_3':
		if purpose == 'grid':
			bottom = 90
			top = 130
			left = 160
			right = 210	

		if purpose == 'edge':		
			top    = [869, 1204]
			left   = [457, 744]
			bottom = [987, 317]
			right  = [1418, 799]

	elif name == 'Polaris':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [0,0]
			left   = [0,0]
			bottom = [0,0]
			right  = [0,0]

	elif name == 'Polaris_2':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [0,0]
			left   = [0,0]
			bottom = [0,0]
			right  = [0,0]

	elif name == 'Draco':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1131,5404]
			left   = [462,1158]
			bottom = [4722,685]
			right  = [5446,4930]
	elif name == 'California_E':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1680, 3273]
			left   = [559,2468]
			bottom = [2134,360]
			right  = [3341,1121]

	elif name == 'California_W':
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [1130,2000]
			left   = [541,1560]
			bottom = [1494,349]
			right  = [2162,817]

	else:
		if purpose == 'grid':
			bottom = None
			top = None
			left = None
			right = None	

		if purpose == 'edge':		
			top    = [0,0]
			left   = [0,0]
			bottom = [0,0]
			right  = [0,0]
	

	return top, left, bottom, right			
