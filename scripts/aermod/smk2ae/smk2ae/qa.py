from builtins import str
from builtins import object
import os.path

class QA(object):
    def __init__(self, work_path):
        self.work_path = work_path
        self.moved_facs = []

    def write_moved_facs(self, state):
        '''
        Write the list of facilities where some sources are moved from one met grid to another
        '''
        qa_name = os.path.join(self.work_path, 'qa', 'moved_facilities_%s.csv' %state)
        with open(qa_name,'w') as qa_file:
            qa_file.write('facility_id\n')
            for fac in self.moved_facs:
                qa_file.write('%s\n' %fac)

    def write_totals(self, state):
        with open(os.path.join(self.work_path, 'qa', 'general_qa_%s.csv' %state), 'w') as qa_file:
            qa_file.write('input_fac-rel-proc,uniq_fac_src,uniq_points,uniq_fugs,total_benzene,total_non_hourly,total_hourly\n')
            try:
                hrly = self.hourly
            except AttributeError:
                hrly = '0'
            line = [str(x) for x in [self.sources,self.uniq_srcs,self.points,self.fugs,self.total_benzene,
              self.n_temp,hrly]]
            qa_file.write('%s\n' %','.join(line))

