#!/usr/bin/env conda run -n ct_extract_env python

import pathlib

def get_release(self):
	data_dir = "{}/data".format(pathlib.Path.cwd().parents[0])

	# below commented block works for dir (result of extracting .zip)
	# dir_list = [folds for folds in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir,folds))]
	# dir_list = [i.split("_")[0] for i in dir_list if "_extracted" in i]
	# file_dates = [dt.datetime.strptime(date, '%Y%m%d').date() for date in test_dir_list] # convert all strings in list into datetime objects
	# latest_file_date = max(file_dates)
	# latest_file_date = latest_file_date.strftime('%m/%d/%Y')
	# return latest_file_date

	

	# below works for zip files (that are unzipped still)
	def get_dates(date):
	    try:
	        return dt.datetime.strptime(date, '%Y%m%d').date()
	    except ValueError:
	        pass
  

	file_list = [f for f in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir, f)) and f.endswith(".zip")]
	file_list = [i.split("_")[0] for i in dir_list]
	file_list = [get_dates(x) for x in file_list]
	file_list = list(filter(None, file_list))

	latest_file_date = max(file_list)
	latest_file_date = latest_file_date.strftime('%m/%d/%Y')
	return latest_file_date


 