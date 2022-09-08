if [ -e bs_results/spectrum_bs.h5 ]; then
	python ~/lalibe/tests/regression/py_scripts/compare_h5_files.py bs_results/spectrum_bs.h5 ~/andre_master/tests/bs_results/spectrum_bs.h5
	rm bs_results/spectrum_bs.h5
else echo "bs_results/spectrum_bs.h5 does not exist"
fi 
