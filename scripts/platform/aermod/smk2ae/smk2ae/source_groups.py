#!/usr/bin/env python

from builtins import object
import pandas as pd

class SourceGroups(object):
    '''
    Read the source groups file
    This file is used to x-ref the SCCs into source groups for area sources
    '''
    def __init__(self, group_file):
        self.xref = self.read_groups(group_file)
        self.scc_list = list(self.xref['scc'].drop_duplicates())

    def read_groups(self, group_file):
        df = pd.read_csv(group_file, dtype={'scc': str, 'run_group': str}, 
          usecols=['scc','run_group','source_group'])
        return df

