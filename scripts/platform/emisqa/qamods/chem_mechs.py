molecDct = dict()

molecDct['cmaq_cb05'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
	  "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, \
	  "PSO4": 1, "PMC": 1,"ALD2": 46.8, "ALDX": 38.6, "CH4": 16.0,\
	  "ETH": 33.1, "ETHA": 30.1,"ETOH": 45.6, "FORM": 29.9, \
	  "IOLE": 55.3, "ISOP": 68.1, "MEOH": 32.0,"NVOL": 15.9, \
	  "OLE": 32.3, "PAR": 16.9, "TERP": 134.4, "TOL": 97.7, \
	  "UNK": 352.4, "UNR": 25.9, "XYL": 109.0, "HGIIGAS": 200.59,\
	  "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
	  "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
	  "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,\
	  "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
	  "BUTADIENE13":54.0904, "BENZENE": 78.1118, "VOC_INV": 1.0, \
	  "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, \
	  "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
	  "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
	  "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
	  "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
	  "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
	   "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,\
           "NAPHTH_72": 128.1705, "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
	   "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
	   "PNA": 1, "PNCOM": 1, "PNH4": 1, "PSI": 1, "PTI": 1, 
	   "ARSENIC_C": 1, "ARSENIC_F": 1, \
	   "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00,
	   "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00,
	   "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00}

molecDct['cmaq_cb05_soa'] = molecDct['cmaq_cb05']

molecDct['cmaq_cb05v2_soa'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
	  "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, \
	  "PSO4": 1, "PMC": 1,"ALD2": 43.8, "ALDX": 37.58, "CH4": 16.04,\
	  "ETH": 28.05, "ETHA": 30.07,"ETOH": 44.73, "FORM": 30.026, \
	  "IOLE": 54.55, "ISOP": 68.12, "MEOH": 32.04,"NVOL": 1.00, \
	  "OLE": 33.09, "PAR": 17.15, "TERP": 136.24, "TOL": 97.98, \
	  "UNK": 352.4, "UNR": 26.29, "XYL": 109.71, "HGIIGAS": 200.59,\
	  "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
	  "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
	  "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,\
	  "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
	  "BUTADIENE13":54.0904, "BENZENE": 78.1118, "VOC_INV": 1.0, \
	  "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, \
	  "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
	  "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
	  "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
	  "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
	  "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
	   "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,\
           "NAPHTH_72": 128.1705, "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
	   "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
	   "PNA": 1, "PNCOM": 1, "PNH4": 1, "PSI": 1, "PTI": 1, 
	   "ARSENIC_C": 1, "ARSENIC_F": 1, \
	   "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00,
	   "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00,
	   "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00}

molecDct['cmaq_cb05v2_mplite'] = molecDct['cmaq_cb05v2_soa']

molecDct['cmaq_cb05_tx'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, 
	  "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, 
	  "PSO4": 1, "PMC": 1,"ALD2": 46.8, "ALDX": 38.6, "CH4": 16.0,
	  "ETH": 33.1, "ETHA": 30.1,"ETOH": 45.6, "FORM": 29.9, 
	  "IOLE": 55.3, "ISOP": 68.1, "MEOH": 32.0,"NVOL": 15.9, 
	  "OLE": 32.3, "PAR": 16.9, "TERP": 134.4, "TOL": 97.7, 
	  "UNK": 352.4, "UNR": 25.9, "XYL": 109.0, "HGIIGAS": 200.59,
	  "HGNRVA":200.59, "PHGI": 1, "DIESEL_PEC": 1, 
	  "DIESEL_PMFINE": 1, "DIESEL_PNO3":1, "DIESEL_POC": 1, 
	  "DIESEL_PSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,
	  "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, 
	  "BUTADIENE13":54.0904, "BENZENE": 78.1118, "VOC_INV": 46.0, 
	  "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, 
	  "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,
	  "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, 
	  "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,
	  "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, 
	  "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, 
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, 
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, 
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,
	   "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,
           "NAPHTH_72": 128.1705, "NH3_FERT": 17}

molecDct['cmaq_saprc07T'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, 
	  "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, 
	  "PSO4": 1, "PMC": 1,"ALD2": 46.8, "ALDX": 38.6, "CH4": 16.0,
	  "ETH": 33.1, "ETHA": 30.1,"ETOH": 45.6, "FORM": 29.9, 
	  "IOLE": 55.3, "ISOP": 68.1, "MEOH": 32.0,"NVOL": 15.9, 
	  "OLE": 32.3, "PAR": 16.9, "TERP": 134.4, "TOL": 97.7, 
	  "UNK": 352.4, "UNR": 25.9, "XYL": 109.0, "HGIIGAS": 200.59,
	  "HGNRVA":200.59, "PHGI": 1, "DIESEL_PEC": 1, 
	  "DIESEL_PMFINE": 1, "DIESEL_PNO3":1, "DIESEL_POC": 1, 
	  "DIESEL_PSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,
	  "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, 
	  "BUTADIENE13":54.0904, "BENZENE": 78.1118, 
	  "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, 
	  "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,
	  "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, 
	  "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,
	  "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, 
	  "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, 
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, 
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, 
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1, "HFC": 66.0999, "HCHO": 1, "CCHO": 1, "BENZ": 78.1118, 
	   "13BDE": 54.0904, "ACRO": 56.0633, "BACL": 86.0892, "ACET": 36.0, "IPRD":70.0898, "PRPE": 42.0797, "PACD": 74.0785,  
	   "CRES": 104.0141, "ETHE": 24.0, "ARO1": 159.5964, "ARO2": 147.3216, "OLE2": 91.6257, "OLE1": 74.8872, "RCHO": 66.7692, 
	    "PRD2": 106.6958, "VOC": 1, "GLY": 58.0361, "APIN": 136.234, "BPIN": 120.0, "B124": 127.5588, "NROG": 102.4983, "BALD": 115.0695, 
	    "MVK": 140.2227, "SESQ": 180.0, "MEK": 72.1046, "ALK4": 82.814, "ALK5": 99.9047, "ALK2": 54.3473, "ALK3": 64.5106, "ALK1": 59.8675, 
	    "AACD": 60.052, "MGLY": 72.0627, "MACR": 70.0898, "ACYE": 26.0373, "FACD": 46.0254, "HCHO_PRIMARY":30.026 }

molecDct['cmaq_cb6'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
          "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, \
          "PSO4": 1, "PMC": 1,"ALD2": 44.0526, "ALDX": 43.65, "CH4": 16.042,\
          "ETH": 28.053, "ETHA": 30.069,"ETOH": 46.0684, "FORM": 30.026, \
          "IOLE": 56.11, "ISOP": 68.117, "MEOH": 32.042,"NVOL": 1.0001, \
          "OLE": 27.65, "PAR": 14.43, "TERP": 136.234, "TOL": 92.138, \
          "UNK": 137.19, "UNR": 28.86, "XYL": 106.165, "HGIIGAS": 200.59,\
          "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
          "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
          "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,\
          "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
          "BUTADIENE13":54.0904, "BENZ": 78.1118, "BENZENE": 78.1118, "VOC_INV": 1.0, \
          "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, \
          "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
          "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
          "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
          "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
          "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
           "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,\
           "NAPHTH_72": 128.1705, "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
           "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
           "PNA": 1, "PNCOM": 1, "PNH4": 1, "PSI": 1, "PTI": 1,
           "ARSENIC_C": 1, "ARSENIC_F": 1, \
           "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00,
           "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00,
           "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00,
           "ACET": 58.079, "ETHY": 26.037, "KET": 28.82, "PRPA": 44.096,
           "SOAALK": 92.1006, "XYLMN": 106.165, "NAPH": 128.1705, "NMOG": 1,
           "ACETONITRILE": 41.05, "ACRYLICACID": 72.06, "ACRYLONITRILE": 53.06,
           "CARBSULFIDE": 60.07, "CHLOROPRENE": 88.54, "ETHYLBENZ": 106.165,
           "HEXANE": 86.175, "METHCHLORIDE": 50.49, "STYRENE": 104.15,
           "XYLENES": 106.165, "VOC_BEIS": 46, "APIN": 120, "BPIN": 120, 
	   "SESQ": 180, "NR": 24, "CO2": 44.01, "CO2_ORG": 44.01}

molecDct['cmaq_cb6ae7'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
          "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, \
          "PSO4": 1, "PMC": 1,"ALD2": 44.0526, "ALDX": 43.65, "CH4": 16.042,\
          "ETH": 28.053, "ETHA": 30.069,"ETOH": 46.0684, "FORM": 30.026, \
          "IOLE": 56.11, "ISOP": 68.117, "MEOH": 32.042,"NVOL": 1.0001, \
          "OLE": 27.65, "PAR": 14.43, "TERP": 136.234, "TOL": 92.138, \
          "UNK": 137.19, "UNR": 28.86, "XYL": 106.165, "HGIIGAS": 200.59,\
          "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
          "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
          "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,\
          "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
          "BUTADIENE13":54.0904, "BENZ": 78.1118, "BENZENE": 78.1118, "VOC_INV": 1.0, \
          "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, \
          "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
          "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
          "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
          "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
          "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
           "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,\
           "NAPHTH_72": 128.1705, "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
           "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
           "PNA": 1, "PNCOM": 1, "PNH4": 1, "PSI": 1, "PTI": 1,
           "ARSENIC_C": 1, "ARSENIC_F": 1, \
           "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00,
           "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00,
           "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00,
           "ACET": 58.079, "ETHY": 26.037, "KET": 28.82, "PRPA": 44.096,
           "SOAALK": 92.1006, "XYLMN": 106.165, "NAPH": 128.1705,
           "ACETONITRILE": 41.05, "ACRYLICACID": 72.06, "ACRYLONITRILE": 53.06,
           "CARBSULFIDE": 60.07, "CHLOROPRENE": 88.54, "ETHYLBENZ": 106.165,
           "HEXANE": 86.175, "METHCHLORIDE": 50.49, "STYRENE": 104.15,
           "XYLENES": 106.165, "VOC_BEIS": 46, "APIN": 136, "BPIN": 120, 
	   "SESQ": 180, "NR": 24, "AACD": 60.052, "FACD": 46.025, 
	   "IVOC": 125.9429, "BENZOAPYRNE": 252.316, "GLY": 58.036, 
	   "ISPD": 70.091, "PACD": 42.22425, "GLYD": 60.052, "MGLY": 72.063,
	   "NMOG": 1.0}

molecDct['cmaq_cb6ae8'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
          "SULF": 98, "PEC": 1, "POC": 1, "PMFINE": 1, "PNO3": 1, \
          "PSO4": 1, "PMC": 1,"ALD2": 44.0526, "ALDX": 43.65, "CH4": 16.042,\
          "ETH": 28.053, "ETHA": 30.069,"ETOH": 46.0684, "FORM": 30.026, \
          "IOLE": 56.11, "ISOP": 68.117, "MEOH": 32.042,"NVOL": 1.0001, \
          "OLE": 27.65, "PAR": 14.43, "TERP": 136.234, "TOL": 92.138, \
          "UNK": 137.19, "UNR": 28.86, "XYL": 106.165, "HGIIGAS": 200.59,\
          "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
          "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
          "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "FORM_PRIMARY":30.026,\
          "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
          "BUTADIENE13":54.0904, "BENZ": 78.1118, "BENZENE": 78.1118, "VOC_INV": 1.0, "VOC": 1.0, \
          "NAPHTHALENE": 128.1705, "CL2": 70.91,"CHROMHEX_C": 1, \
          "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
          "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
          "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
          "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
          "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
           "PEC_72": 1, "PMFINE_72": 1, "POC_72": 1, "PMC_72": 1, "OTHER": 1,\
           "NAPHTH_72": 128.1705, "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
           "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
           "PNA": 1, "PNCOM": 1, "PNH4": 1, "PSI": 1, "PTI": 1, 
           "ARSENIC_C": 1, "ARSENIC_F": 1, \
           "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00, \
           "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00, \
           "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00, \
           "ACET": 58.079, "ETHY": 26.037, "KET": 28.82, "PRPA": 44.096, \
           "SOAALK": 92.1006, "XYLMN": 106.165, "NAPH": 128.1705, \
           "ACETONITRILE": 41.05, "ACRYLICACID": 72.06, "ACRYLONITRILE": 53.06, \
           "CARBSULFIDE": 60.07, "CHLOROPRENE": 88.54, "ETHYLBENZ": 106.165, \
           "HEXANE": 86.175, "METHCHLORIDE": 50.49, "STYRENE": 104.15, \
           "XYLENES": 106.165, "NOX_INV": 46, "HF": 20.01, "NMOG": 1.0, \
           "VOC_BEIS": 46, "APIN": 136.234, "BPIN": 120, "SESQ": 180, "NR": 24, \
           "AACD": 60.052, "FACD": 46.025, "IVOC": 125.9429, "BENZOAPYRNE": 252.316, \
           "TOG_INV": 1, "GLY": 58.036, "GLYD": 60.052, "ISPD": 70.091, "IVOCP3": 148.7341, \
           "IVOCP4": 163.4046, "IVOCP5": 119.0065, "IVOCP5ARO": 175.4026, "IVOCP6": 112.4847, "IVOCP6ARO": 102.8674, \
           "MGLY": 72.063, "PNCOMN1": 1, "PNCOMN2": 1, "PNCOMP0": 1, "PNCOMP1": 1, "PNCOMP2": 1, 
           "POCN1": 1, "POCN2": 1, "POCP0": 1, "POCP1": 1, "POCP2": 1, "SVOCN1": 228.294,
           "SVOCP0": 147.2124, "SVOCP1": 251.4408, "SVOCP2": 278.1936, "PACD": 42.22425}
					 
molecDct['cmaq_cracmmv1'] = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
            "SULF": 98, "PEC": 1, "PMOCN2": 1, "PMFINE": 1, "PNO3": 1, \
            "PSO4": 1, "PMC": 1, "FORM_PRIMARY":30.026,\
            "ALD2_PRIMARY": 44.0526, "ACROLEIN": 56.0633, \
            "BUTADIENE13":54.0904, "VOC_INV": 1.0, "VOC": 1, \
            "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
            "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
            "PFE": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
            "PNA": 1, "PMNCOMN2": 1, "PNH4": 1, "PSI": 1, "PTI": 1, 
            "SOAALK": 92.1006, "NMOG": 1.0, "VOC_BEIS": 46, "TOG_INV": 1, \
            "ACD": 44.053, "ACE": 26.037, "ACRO": 56.063, "ACT": 58.079, \
            "ALD": 65.3665, "API": 136.234, "BALD": 115.0873, "BDE13": 54.09, \
            "BEN": 78.112, "CSL": 122.0033, "ECH4": 16.042, "EOH": 46.068, \
            "ETE": 28.053, "ETEG": 62.068, "ETH": 30.069, "FURAN": 72.2989, \
            "GLY": 58.036, "HC3": 58.4107, "HC5": 82.5977, "HC10": 103.2933, \
            "HCHO": 30.026, "ISO": 68.117, "KET": 106.8355, "LIM": 136.234, \
            "MEK": 72.106, "MOH": 32.0419, "MVK": 70.09, "NAPH": 145.4619, \
            "OLI": 93.33, "OLT": 65.5403, "ORA1": 46.025, "ORA2": 64.3665, \
            "PHEN": 94.111, "PROG": 76.094, "VROCIOXY": 111.4328, "ROCP0": 147.211, \
            "ROCP1": 251.467, "ROCP2": 292.8446, "ROCP3": 240.468, "ROCP4": 239.6231, \
            "ROCP5": 210.6373, "ROCP5ARO": 171.1934, "ROCP6": 140.5942, "ROCP6ARO": 113.636, \
            "ROH": 75.5254, "SLOWROC": 107.0321, "TOL": 92.138, "XYM": 111.6467, \
            "XYE": 115.5535, "DCB1": 100.116, "HKET": 116.158, "MACR": 70.09, \
            "MCT": 124.575, "MGLY": 72.063, "ONIT": 105.09, "ROCN1": 212.246, \
            "ROCN2": 292.243, "SESQ": 204.357, "UALD": 84.116, "UNKKOH": 70.91, \
            "UNKSMILES": 353.37, "PNCOMN1": 1, "PNCOMN2": 1, "PNCOMP0": 1, \
	    "PNCOMP1": 1, "PNCOMP2": 1, "POCN1": 1, "POCN2": 1, \
	    "POCP0": 1, "POCP1": 1, "POCP2": 1, \
	    "CL2": 70.91}
