#!/usr/bin/env python

import sys
import os.path
import pandas as pd
import xlsxwriter as xl

class QA(object):

    def __init__(self, case, run_group, aermod_path):
        self.case = case
        self.run_group = run_group
        self.path = aermod_path
        self.wkb_name = os.path.join(self.path, 'qa', 'spreadsheets', '%s_%s_qa.xlsx' %(case, run_group))
        print('Writing to %s' %self.wkb_name)
        self.tabs = {'counts': 'helper file counts', 'checks': 'helper file check', 
          'emis': 'source county emissions', 'poll': 'pollutant emissions', 
          'scc': 'source scc emissions', 'temporal': 'temporal check', 'auto': 'automated checks',
          'grid': 'grid cell check', 'locs': 'point location check',
          'pt_params': 'point source parameters','fug_params': 'fugitive parameters',
          'cmv_params': 'CMV source parameters',
          'temporal_hourly': 'hourly temporal check','onroad_grid': 'onroad grid cell check'}
        self.tabdesc = {'counts': 'Number of unique cells or facilities and unique source combinations per AERMOD helper file for this run group.',
          'checks': 'Check of whether a source found in any helper file is also found in all other helper files.',
          'emis': 'County and/or AERMOD source group level emissions comparison between FF10 and emissions helper file.',
          'poll': 'Pollutant level emissions comparison between FF10 and emissions helper file.',
          'scc': 'Helper file emissions by SCC.',
          'temporal': 'Sum of scalars in non-hourly temporal helper file by AERMOD source. Sum formula varies by qflag as noted in Table 4 of NONPOINT and ONROAD design document.',
          'auto': 'Results of miscellaneous automated checks that do not output to a spreadsheet. Generally the results should be "none".',
          'grid': 'Grid cell alignment checks',
          'locs': 'Point source locations',
          'pt_params': 'Stack source release parameters',
          'fug_params': 'Fugitive source release parameters',
          'cmv_params': 'CMV source release parameters',
          'temporal_hourly': 'Sum of scalars hourly temporal helper files by county.',
          'onroad_grid': 'Onroad grid cell alignment checks'}
        self.tabproc = {'counts': self.write_counts, 'checks': self.write_checks,
          'emis': self.write_emis, 'poll': self.write_poll, 'scc': self.write_scc,
          'temporal_hourly': self.write_hourly, 'temporal': self.write_temporal,
          'onroad_grid': self.write_onroad_grid, 'grid': self.write_grid,
          'auto': self.write_auto, 'locs': self.write_locs, 'pt_params': self.write_pt_params,
          'fug_params': self.write_fug_params, 'cmv_params': self.write_cmv_params}
        self.hourly = False
        self.runtype = 'area'
        self.tabs_order = ['counts','checks','emis','poll','scc','temporal','grid','auto']
        if self.run_group.startswith('point'):
            self.tabs_order = ['counts','checks','locs','pt_params','fug_params','emis','poll',
              'temporal','temporal_hourly','auto']
            self.runtype = 'point'
            self.hourly = True
        elif self.run_group.startswith('CMV'):
            self.tabs_order = ['counts','checks','locs','cmv_params','emis','poll','temporal_hourly','auto']
            self.runtype = 'cmv'
            self.hourly = True
        elif self.run_group.startswith('ptair') or self.run_group.startswith('air'):
            #self.tabs_order = ['counts','checks','line_locs','nr_locs','line_params','nr_params',
            #  'emis','poll','temporal']
            self.tabs_order = ['counts','checks','locs','emis','poll','temporal','auto']
            self.runtype = 'airports'
        elif self.run_group[:3] in ('HOT','HDO','LDO'):
            self.hourly = True
            self.tabs_order = ['counts','checks','emis','poll','scc','temporal_hourly','onroad_grid','auto']
        elif 'RWC' in self.run_group or 'AG' in self.run_group:
            self.hourly = True
            self.tabs_order = ['counts','checks','emis','poll','scc','temporal_hourly','grid','auto']

    def __enter__(self):
        self.writer = pd.ExcelWriter(self.wkb_name, engine='xlsxwriter')
        self.wkb = self.writer.book
        return self

    def __exit__(self, type, value, traceback):
        self.writer.save()

    def init_formats(self):
        self.bold = self.wkb.add_format({'bold': 1})
        self.emis = self.wkb.add_format({'num_format': '#,##0.00'}) # eg 1,412.23
        self.pct = self.wkb.add_format({'num_format': '0.00%'})  # eg 5.33%
        self.whole = self.wkb.add_format({'num_format': '#,##0'})  # eg 1,412
        self.small_emis = self.wkb.add_format({'num_format': '0.000'}) # eg 2.235

    def write_readme(self):
        wks = self.wkb.add_worksheet('README')
        wks.write(0, 0, '%s AERMOD Helper QA' %self.case, self.bold)
        wks.write(2, 0, '%s inventory' %self.run_group, self.bold)
        wks.write(5, 0, '%s helper files' %self.run_group, self.bold) 
        wks.write(5, 1, self.path)
        for row, helper in enumerate(['Locations','Parameters','Temporal','Emissions']):
            wks.write(6+row, 0, helper)
            if helper == 'Parameters':
                fn = [self.run_group, 'area_params']
            elif helper == 'Emissions':
                if self.run_group.endswith('12'):
                    res = '12'
                elif self.run_group.endswith('4'):
                    res = '4'
                else:
                    res = self.run_group[-3]
                fn = [res, self.run_group, 'emis']
            elif helper == 'Temporal' and self.hourly:
                fn = [self.run_group, '[STFIPS]', 'hourly']
            else:
                fn = [self.run_group, helper.lower()]
            wks.write(6+row, 1, '%s.csv' %'_'.join(fn))
        if self.hourly and self.runtype != 'point':
            wks.write(10, 0, 'County X-Walk')
            wks.write(10, 1, '%s_county-to-gridcell.csv' %self.run_group)
        note_a = 'Note - inventory sources are only selected from the FF10 for those sources that have NATA emissions > 0, that are within the domain, and meet all run group criteria.'
        wks.write(11, 1, note_a)
        if self.runtype not in ['point','airports']:
            note_b = 'For area sources the source_groups file is used to cross-reference the inventory SCC to the run group.'
            wks.write(12, 1, note_b)
        wks.write(13, 0, 'Tab', self.bold)
        wks.write(13, 1, 'Description', self.bold)
        for row, tab in enumerate(self.tabs_order):
            wks.write(14+row, 0, self.tabs[tab])
            wks.write(14+row, 1, self.tabdesc[tab])
        wks.set_column(0, 0, 25)

    def write_header(self, row, header, wks):
        for n, col in enumerate(header):
            wks.write(row, n, col, self.bold)

    def write_counts(self):
        '''
        Write the counts of facility and source
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_counts.csv' %self.run_group)
        df = pd.read_csv(fn)
        sheet_name = self.tabs['counts']
        print(df.head())
        if self.runtype == 'point':
            df.columns = ['Location facilities','Location sources','Fugitive parameter facilities',
              'Fugitive parameter sources','Point parameter facilities','Point parameter sources',
              'Total parameter facilities','Total parameter sources','Temporal facilities',
              'Temporal sources','Emissions facilities','Emissions sources','Inventory facilities'] 
        elif self.runtype == 'cmv':
            df.columns = ['Location facilities','Location sources','Parameter facilities',
              'Parameter sources','Temporal facilities','Temporal sources',
              'Emissions facilities','Emissions sources']
        elif self.runtype == 'airports':
            df.columns = ['Line location facs','Line location sources','Nonrunway location facs',
              'Nonrunway location sources','Total location facs','Total location sources',
              'Nonrunway param facs','Nonrunway param sources','Line param facs','Line param sources',
              'Total param facs','Total param sources','Temporal facs','Temporal sources','Emis facs','Inv facs'] 
        else:
            df.columns = ['Location cells', 'Location sources', 'Params cells', 'Parms sources', 
              'Temporal cells', 'Temporal sources', 'Emissions cells', 'Emissions sources']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(0, df.columns, wks)
        wks.set_column('A:%s' %chr(64+len(df.columns)), 16)

    def write_checks(self):
        '''
        Write the source file check
        ''' 
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_srcid_qa.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'facility_id': str})
        sheet_name = self.tabs['checks']
        if self.runtype == 'point':
            df.columns = ['facility_id','src_id','locations','parameters (point or fugitive)','temporal','emissions','xwalk']
        elif self.runtype == 'airports':
            df.columns = ['facility_id','src_id','locations','parameters','temporal']
        else:
            df.columns=['met_cell','src_id','locations','parameters','temporal','emissions']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=2, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, df.columns, wks)
        wks.set_column('A:%s' %chr(64+len(list(df.columns))), 12)
        wks.merge_range('C1:E1', 'Source in helper file (Y/N)')
        wks.autofilter('A2:%s2' %chr(64+len(list(df.columns))))
        wks.freeze_panes(2,2)
 
    def write_emis(self):
        '''
        Write the source level emissions qa check
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_emis_qa.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'region_cd': str, 'facility_id': str})
        if self.runtype == 'point':
            df.columns = ['Facility ID','Source ID','FF10','Helper']
        elif self.runtype == 'airports':
            df.columns = ['Facility ID','FF10','Helper']
            df.insert(1, 'Runway', 'All')
        else:
            df.columns = ['FIPS','Source Group','FF10','Helper']
        sheet_name = self.tabs['emis']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=2, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, df.columns, wks)
        rows = len(df)
        wks.write('E2', 'Diff', self.bold)
        wks.write_array_formula('E3:E%s' %(2+rows), '=D3:D%s-C3:C%s' %(2+rows,2+rows), self.emis)
        wks.write('F2', '% Diff', self.bold)
        wks.write_array_formula('F3:F%s' %(2+rows), '=E3:E%s/C3:C%s' %(2+rows,2+rows), self.pct)
        wks.set_column('A:A', 10)
        wks.set_column('B:B', 16)
        wks.set_column('C:E', 16, self.emis)
        wks.set_column('F:F', 16, self.pct)
        wks.merge_range('C1:E1', 'Annual emissions (tons/yr)')
        wks.autofilter('A2:F2')
        wks.freeze_panes(2,2)
 
    def write_poll(self):
        '''
        Write the pollutant emissions comparison
        '''
        usecols = ['smoke_name','ann_value_aer','ann_value_inv']
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_poll_qa.csv' %self.run_group)
        df = pd.read_csv(fn, usecols=usecols)
        df.columns = ['SMOKE Name','FF10','Helper']
        sheet_name = self.tabs['poll']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=2, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, df.columns, wks)
        rows = len(df)
        wks.write('D2', 'Diff', self.bold)
        wks.write_array_formula('D3:D%s' %(2+rows), '=C3:C%s-B3:B%s' %(2+rows,2+rows), self.emis)
        wks.write('E2', '% Diff', self.bold)
        wks.write_array_formula('E3:E%s' %(2+rows), '=D3:D%s/B3:B%s' %(2+rows,2+rows), self.pct)
        wks.set_column('A:A', 16)
        wks.set_column('B:D', 12, self.emis)
        wks.set_column('E:E', 12, self.pct)
        wks.merge_range('B1:D1', 'Annual emissions (tons/yr)')

    def write_scc(self, scc_desc='/work/EMIS/em_v8/ge_dat/smkreport/sccdesc_2017platform_22nov2019_v0.txt'):
        '''
        Write the SCC emissions check
        '''
        scc_desc = pd.read_csv(scc_desc, names=['scc','scc_desc'], dtype={'scc': str}, comment='#')
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_scc_qa.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'scc': str})
        df = pd.merge(df, scc_desc, on='scc', how='left')
        df = df[['source_group','scc','scc_desc','ann_value']].copy()
        df.columns = ['Source Group','SCC','SCC Description','FF10 Emissions selected for AERMOD (tons/yr)']
        sheet_name = self.tabs['scc']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(0, df.columns, wks)
        wks.set_column('A:A', 32)
        wks.set_column('B:B', 12)
        wks.set_column('C:C', 50)
        wks.set_column('D:D', 32, self.small_emis)
        wks.autofilter('A1:D1')
        wks.freeze_panes(1,3)

    def write_hourly(self):
        '''
        Write the hourly temporal check
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', 'temporal_hour_%s_check.csv' %self.run_group)
        if os.path.exists(fn):
            df = pd.read_csv(fn, dtype={'fips': str, 'facility_id': str})
            if self.runtype in ('point','cmv'):
                df.columns = ['Facility ID','Source ID','Sum of scalars','Count of Scalars']
            else:
                df.columns = ['FIPS','Sum of scalars','Count of Scalars']
        else:
            df = pd.DataFrame([['No hourly temporal files found',self.run_group,self.runtype],])
        sheet_name = self.tabs['temporal_hourly']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(0, df.columns, wks)
        wks.set_column('A:%s' %chr(64+len(df.columns)), 16)
        wks.autofilter('A1:%s1' %chr(64+len(df.columns)))
        wks.freeze_panes(1,1)

    def write_temporal(self):
        '''
        Write the regular temporal check
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_temporal_check.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'facility_id': str})
        print(df.head())
        if self.runtype in ('point','cmv','airports'):
            df.columns = ['Facility ID','Source ID','Qflag','Sum of scalars']
        else:
            df.columns = ['met_cell','src_id','Qflag','Sum of scalars']
        sheet_name = self.tabs['temporal']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(0, df.columns, wks)
        wks.set_column('A:D', 16)
        wks.autofilter('A1:D1')
        wks.freeze_panes(1,1)

    def write_grid(self):
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_grid_comparison.csv' %self.run_group)
        usecols = ['met_cell','region_cd','rep_value','ann_value','diff','pdiff']
        df = pd.read_csv(fn, usecols=usecols, dtype={'region_cd': str})
        sheet_name = self.tabs['grid']
        df = df[usecols].copy() # Make sure the order is set correctly
        df.insert(0, 'run_group', self.run_group)
        df.columns = ['Run Group','Met Cell','FIPS','SMOKE Reports','AERMOD Helpers','Diff','% Diff']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=2, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, df.columns, wks)
        wks.set_column('A:C', 12)
        wks.set_column('D:F', 16, self.small_emis)
        wks.set_column('G:G', 16, self.pct)
        wks.merge_range('D1:F1', 'Annual Emissions (tons/yr)')
        wks.autofilter('A2:G2')
        wks.freeze_panes(2,3)

    def write_onroad_grid(self):
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_grid_comparison.csv' %self.run_group)
        usecols = ['met_cell','region_cd','in_rep','in_helper']
        df = pd.read_csv(fn, usecols=usecols, dtype={'region_cd': str})
        df.insert(0, 'run_group', self.run_group)
        sheet_name = self.tabs['onroad_grid']
        df.columns = ['Run Group','Met Cell','FIPS','SMOKE Reports','AERMOD Helpers']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=2, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, df.columns, wks)
        wks.set_column('A:C', 12)
        wks.set_column('D:E', 18)
        wks.merge_range('D1:E1', 'Cell Associated with FIPS and Source Group (Y/N)')
        wks.autofilter('A2:E2')
        wks.freeze_panes(2,4)

    def write_auto(self):
        wks = self.wkb.add_worksheet(self.tabs['auto'])
        wks.set_column('A:A', 40)
        wks.set_column('B:B', 20)
        wks.write(0, 0, 'Type', self.bold) 
        wks.write(0, 1, 'Results', self.bold) 
        types = ['Met cells with multiple UTM zones', 'Wrong run group in helpers',
          'Wrong number of verts (!=4)', 'Non-matching rh', 'Non-matching sz',
          'Unaligned met cells from subset (utmx!=x4 from neighboring cell)',
          'Met cell range for cells with emissions', 'Cells not in both locations and parameters',
          'Sum of monthly values does not equal annual value']
        for row, cell in enumerate(types):
            wks.write(1+row, 0, cell)
        wks.write(11, 0, 'Visual grid cell alignment and distance check. UTM, Lambert, and lat/lon coordinates are projected on the Lambert CONUS grid.')
        wks.write(12, 0, 'Distances measured with GIS ruler tool. Each point type is checked to be at the same Lambert location.')
        wks.write(13, 0, 'Order of vertices is also verified within each projection type.')

    def write_locs(self):
        '''
        Write the locations file
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_locs_check.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'facility_id': str, 'fac_source_type': str})
        if self.runtype == 'airports':
            usecols = ['facility_id','src_id','longitude_aer','latitude_aer','utm_x','utm_y','utm_zone',
              'longitude_inv','latitude_inv','utm_zone_inv','utm_x_inv','utm_y_inv']
            df = df[usecols].copy()
            df.columns = ['Facility ID','Source ID','Helper Longitude','Helper Latitude','Helper UTM x',
              'Helper UTM y','Helper UTM Zone','Inventory Longitude',
              'Inventory Latitude','Inventory UTM Zone','Inventory UTM x','Inventory UTM y']
        else:
            usecols = ['facility_id','src_id','longitude_aer','latitude_aer','utm_x','utm_y','utm_zone',
              'fac_source_type','longitude_inv','latitude_inv','utm_zone_inv','x_inv','y_inv']
            df = df[usecols].copy()
            df.columns = ['Facility ID','Source ID','Helper Longitude','Helper Latitude','Helper UTM x',
              'Helper UTM y','Helper UTM Zone','Facility Source Type','Inventory Longitude',
              'Inventory Latitude','Inventory UTM Zone','Inventory UTM x','Inventory UTM y']
        sheet_name = self.tabs['locs']
        cols = list(df.columns) + ['UTM x Diff','UTM y Diff','Longitude Diff','Latitude Diff']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, cols, wks)
        # Difference
        for row in range(len(df)):
            wks.write_formula('N%s' %(row+3), '=E%s-L%s' %(row+3, row+3))
            wks.write_formula('O%s' %(row+3), '=F%s-M%s' %(row+3, row+3))
            wks.write_formula('P%s' %(row+3), '=C%s-I%s' %(row+3, row+3))
            wks.write_formula('Q%s' %(row+3), '=D%s-J%s' %(row+3, row+3))
        wks.set_column('A:B', 12)
        wks.set_column('C:G', 16, self.small_emis)
        wks.set_column('I:Q', 16, self.small_emis)
        wks.autofilter('A2:Q2')
        wks.freeze_panes(2,2)

    def write_cmv_params(self):
        '''
        Write the CMV parameters file
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_point_param_qa.csv' %self.run_group)
        print(fn)
        df = pd.read_csv(fn, dtype={'facid': str})
        sheet_name = self.tabs['cmv_params']
        df.columns = ['Facility ID','Source ID','Release Height Helper (ft)',
          'Sigma Z','Release Height Group (ft)','Sigma Z Group']
        cols = list(df.columns) + ['Release Heigh Diff','Sigma Z Diff']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, cols, wks)
        # Difference
        for row in range(len(df)):
            wks.write_formula('G%s' %(row+3), '=C%s-E%s' %(row+3, row+3))
            wks.write_formula('H%s' %(row+3), '=D%s-F%s' %(row+3, row+3))
        wks.set_column('A:B', 12)
        wks.set_column('C:H', 16, self.small_emis)
        wks.autofilter('A2:H2')
        wks.freeze_panes(2,2)

    def write_pt_params(self):
        '''
        Write the point parameters file
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_point_param_qa.csv' %self.run_group)
        usecols = ['facility_id','src_id','stkhgt_aer','stktemp_aer','stkvel_aer','stkdiam_aer',
          'stkhgt_inv','stktemp_inv','stkvel_inv','stkdiam_inv']
        print(fn)
        df = pd.read_csv(fn, usecols=usecols, dtype={'facility_id': str})
        sheet_name = self.tabs['pt_params']
        df.columns=['Facility ID','Source ID','Helper Stack Height (ft)','Helper Stack Temp (F)',
          'Helper Stack Vel (ft/s)','Helper Stack Diam (ft)','Inv Stack Height (ft)','Inv Stack Temp (F)',
          'Inv Stack Vel (ft/s)','Inv Stack Diam (ft)']
        cols = list(df.columns) + ['Height Diff','Temp Diff','Vel Diff','Diam Diff']
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, cols, wks)
        # Difference
        for row in range(len(df)):
            wks.write_formula('K%s' %(row+3), '=C%s-G%s' %(row+3, row+3))
            wks.write_formula('L%s' %(row+3), '=D%s-H%s' %(row+3, row+3))
            wks.write_formula('M%s' %(row+3), '=E%s-I%s' %(row+3, row+3))
            wks.write_formula('N%s' %(row+3), '=F%s-J%s' %(row+3, row+3))
        wks.set_column('A:B', 12)
        wks.set_column('C:N', 16, self.small_emis)
        wks.autofilter('A2:N2')
        wks.freeze_panes(2,2)

    def write_fug_params(self):
        '''
        Write the point parameters file
        '''
        fn = os.path.join(self.path, 'qa', 'csvs', '%s_fug_param_qa.csv' %self.run_group)
        df = pd.read_csv(fn, dtype={'facility_id': str})
        sheet_name = self.tabs['fug_params']
        df.columns=['Facility ID','Source ID','Helper Height (ft)','Helper Width (ft)','Helper Length (ft)',
          'Helper Angle','Inv Height (ft)','Inv Width (ft)','Inv Length (ft)','Inv Angle']
        cols = list(df.columns) + ['Height Diff','Width Diff','Length Diff'] 
        df.to_excel(self.writer, sheet_name=sheet_name, index=False, startrow=1, header=False)
        wks = self.writer.sheets[sheet_name]
        self.write_header(1, cols, wks)
        # Difference
        for row in range(len(df)):
            wks.write_formula('K%s' %(row+3), '=C%s-G%s' %(row+3, row+3))
            wks.write_formula('L%s' %(row+3), '=D%s-H%s' %(row+3, row+3))
            wks.write_formula('M%s' %(row+3), '=E%s-I%s' %(row+3, row+3))
        wks.set_column('A:B', 12)
        wks.set_column('C:M', 16, self.small_emis)
        wks.autofilter('A2:M2')
        wks.freeze_panes(2,2)

def main():
    '''
    aermod_path = '/work/EMIS/em_v7.2/2014platform/2014fd_nata_aermod/aermod'
    case = '2014fd_aermod'
    run_group = 'HDOFF3HI'
    '''
    if len(sys.argv) != 4:
        raise ValueError('./prog rungroup case aermod_path')
    run_group, case, aermod_path = sys.argv[1:]
    with QA(case, run_group, aermod_path) as qa:
        qa.init_formats()
        qa.write_readme()
        for tab in qa.tabs_order:
            qa.tabproc[tab]()

if __name__ == '__main__':
    main()
