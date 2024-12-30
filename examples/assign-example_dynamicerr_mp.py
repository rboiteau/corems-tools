import os
import warnings
import pandas as pd
import multiprocessing as mp
import seaborn as sns
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
import sys
sys.path.append('./')

module_dir = 'c:/Users/15087/Desktop/CoreMSdev/CoreMS'
os.environ['COREMS_DATABASE_URL'] = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
sys.path.append(module_dir)

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration


def assign_formula(file, times, interval, error, ref): 


	# peak picking and assignment
	MSParameters.mass_spectrum.min_picking_mz=50	# integer, minimum m/z that will be assigned
	MSParameters.mass_spectrum.max_picking_mz=800 # integer, maximum m/z that will be assigned
	MSParameters.mass_spectrum.noise_threshold_method = 'log'
	MSParameters.mass_spectrum.noise_threshold_log_nsigma = 10 # integer, should be ~7-10 for most cases (lower will keep more noise)
	MSParameters.mass_spectrum.noise_threshold_log_nsigma_bins = 500 # integer, should be kept around 500
	MSParameters.ms_peak.peak_min_prominence_percent = 0.02 # relative intensity of lowest vs highest peak in spectrum. Keep around 0.01 to 0.1 . 

	# assigment & database
	MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@corems-molformdb-1:5432/coremsapp'
	MSParameters.molecular_search.db_chunk_size = 500
	MSParameters.molecular_search.error_method = None
	MSParameters.mass_spectrum.noise_threshold_method = 'log'
	MSParameters.molecular_search.min_ppm_error = -0.25 # acceptable post-calibration minimum mass error for assigments (orbitraps, window should be 1-2ppm)
	MSParameters.molecular_search.max_ppm_error = 0.25 # acceptable post-calibration maximum mass error for assigments (orbitraps, window should be 1-2ppm)
	#MSParameters.molecular_search.mz_error_score_weight: float = 0.6
	#MSParameters.molecular_search.isotopologue_score_weight: float = 0.4
    #MSParameters.molecular_search.score_method = "prob_score"
    #MSParameters.molecular_search.output_score_method = "prob_score"
	MSParameters.ms_peak.legacy_resolving_power = False

	# calibration
	MSParameters.mass_spectrum.min_calib_ppm_error = error-2 # minimum ppm error for detecting m/z calibrant peaks
	MSParameters.mass_spectrum.max_calib_ppm_error = error+2 # maximum ppm error for detecting m/z calibrant peaks
	MSParameters.mass_spectrum.calib_pol_order = 2
	#MSParameters.mass_spectrum.calib_sn_threshold = 50

	MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@corems-molformdb-1:5432/coremsapp'
	MSParameters.molecular_search.score_method = "prob_score"
	MSParameters.molecular_search.output_score_method = "prob_score"
	print("Loading file: "+ file)
	
	parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
	parser.chromatogram_settings.scans = (-1, -1)
	tic=parser.get_tic(ms_type='MS',smooth=False,peak_detection=False)[0]
	tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
	
	ref_file = file.split('.')[0] +'_'+ ref
	refmasslist = ref_file

	results = []
	for timestart in times:
		print(timestart)
		scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
		mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)
		
		MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()
		mass_spectrum.molecular_search_settings.min_dbe = 0
		mass_spectrum.molecular_search_settings.max_dbe = 20
		mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 40)
		mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 80)
		mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
		mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 10)
		mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 2)
		mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)

		mass_spectrum.molecular_search_settings.isProtonated = True
		mass_spectrum.molecular_search_settings.isRadical = False
		mass_spectrum.molecular_search_settings.isAdduct = False
		mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'D': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'S': 2
                                                                        }
		
		SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
		mass_spectrum.percentile_assigned(report_error=True)
		
		assignments=mass_spectrum.to_dataframe()
		assignments['Time']=timestart
		results.append(assignments)
		
	results=pd.concat(results,ignore_index=True)
	results['file'] = file
	results['Molecular Class']=results['Molecular Formula'].str.replace('\d+', '').str.replace(' ', '')
	results['Molecular Class'][results['Heteroatom Class']=='unassigned']='unassigned'
	results['Molecular Class'][results['Is Isotopologue']==1]='Isotope'
	return(results)

def errorplot(LCMS_annotations,filename):
 
	#### Plot and save error distribution figure
	fig, ((ax1, ax2)) = plt.subplots(1,2)
	fig.set_size_inches(12, 6)
	sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular Class',data=LCMS_annotations,ax=ax1, edgecolor='none')
	ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
	ax1.set_title('a', fontweight='bold', loc='left')
	sns.kdeplot(x='m/z Error (ppm)',data=LCMS_annotations,hue='Time',ax=ax2,legend=False)
	ax2.set_title('b', fontweight='bold', loc='left')
	fig.tight_layout()
	fig.savefig(filename, dpi=200,format='jpg')   

def rt_assign_plot(LCMS_annotations,filename):

	#### Plot library assignments over time
	assign_summary=[]
	for time in LCMS_annotations['Time'].unique():
		current={}
		current['Time']=time
		for mol_class in LCMS_annotations['Molecular Class'].unique():
			current[mol_class]=len(LCMS_annotations[(LCMS_annotations['Molecular Class']==mol_class) & (LCMS_annotations['Time']==time)])
		assign_summary.append(current)

	df=pd.DataFrame(assign_summary)
	df=df.sort_values(by='Time')

	df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
	plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
	plt.savefig(filename, bbox_inches='tight',format='jpg')


# Define CoreMS LCMS functions
def run_assignment(time_min,time_max,interval,data_dir,file,refmasslist,error):
	times = list(range(time_min,time_max,interval))

	output=assign_formula(file, times, interval, error, refmasslist)
	fname = file.replace('.raw','')
	output.to_csv(data_dir+fname+'.csv')

	#errorplot(output,data_dir+fname+'_errorplot.jpg')
	#rt_assign_plot(output,data_dir+fname+'_rt_assign_plot.jpg')


if __name__ == '__main__':

    data_dir = 'C:/Users/15087/Desktop/CoreMS/BATS Data/'
    refmasslist = "calibrants.ref"
    samplelist_file = "BATS_samplelist.csv"
    
    samplelist=pd.read_csv(data_dir+samplelist_file)

    results = []

    interval = 2
    time_min = 10
    time_max = 14

    os.chdir(data_dir)
    samplelist=pd.read_csv(data_dir+samplelist_file)
	
    #Run first file alone to ensure that the docker database is generated.
    run_assignment(time_min,time_max,interval,data_dir,samplelist.loc[0,'File'],refmasslist,samplelist.loc[0,'m/z error (ppm)'])

    args=[]
    for index, row in samplelist.iloc[1:].iterrows():
        file=row.loc['File'] 
        error=row.loc['m/z error (ppm)']
        args.append((time_min,time_max,interval,data_dir,file,refmasslist,error))
    with mp.Pool(processes=4) as pool:
        results = pool.starmap(run_assignment, args)
