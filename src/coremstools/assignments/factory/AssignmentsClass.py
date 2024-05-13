from pandas import DataFrame, read_csv, unique
from numpy import mean, std, array, zeros, shape, where, log10, average
from tqdm import tqdm

from coremstools.assignments.calc.AssignmentError import AssignmentError
from coremstools.assignments.calc.QualityControl import QualityControl
from coremstools.assignments.calc.Dispersity import Dispersity 
from coremstools.Parameters import Settings
import coremstools.assignments.calc.Helpers as lcmsfns

class Assignments(QualityControl, Dispersity, AssignmentError):


    def __init__(self, sample_df, t_interval = 2):
        
        self.t_int = t_interval
        self.sample_df = sample_df     # dataframe

    def add_mol_class(self):

        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            df = read_csv(fpath)
            heter = lcmsfns.get_heteroatoms(df)
            molclasses = lcmsfns.get_mol_class(heter)
            df2 = lcmsfns.assign_mol_class(df,molclasses)
            df2.to_csv(fpath, index = False)

    def run_internal_std_qc(self,timerange):

        temp = QualityControl.StandardQC(self.sample_df, timerange)

    def run_assignment_error_plot(self):

        for f in self.sample_df['File']:

            fpath = Settings.assignments_directory + f.split('.')[0] + Settings.csvfile_addend + '.csv'
            save_file = fpath.split('.')[0] + '_mz-error.jpg'
            AssignmentError.ErrorPlot(read_csv(fpath), save_file)
            save_file = fpath.split('.')[0] + '_rt-error.jpg'
            AssignmentError.RTAssignPlot(read_csv(fpath), save_file)

    def run_dispersity_calculation(self):

        print('running dispersity calculation on ...')

        for f in self.sample_df['File']:
            print('  ' + f)
            Dispersity.CalculateDispersity(f, self.t_int)






