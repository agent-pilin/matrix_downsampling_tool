# Count table downsampling tool

This tool is destinated to downsample count tables in order to obtain desired **RPIB** that can be **RRPIB** (mean **raw** read count per barcode) or **MRPIB** (mean **aligned** read count per barcode). Downsampled count tables can be used for correct sample comparison.

### Input

Input file is the yaml file with input parameters:

* **count_table_folders:** folder(s) containing count tables to downsample.

* **downsampling_folders:** folder(s) containing downsampling statistics per barcode;

* **whitelist_folders:** folder(s) containing whitelist(s) of desired barcodes;

* **count_prefix:** prefix of count table names;

* **downsampling_prefix:** prefix of downsampling statistics files;

* **whitelist_prefix:** prefix of whitelist files;

* **samples:** sample list of count tables to downsample;

* **output_folder:** folder where output count tables will be written. If the folder does not exist, it will be created;

* **output_prefix:** prefix for output filenames;

* **add_RPIB_prefix:** **1** if you want to add prefix *'rpib' + [RPIB value]*, **0** otherwise;

* **desired_RPIB:** desired value of RRPIB (or MRPIB) for downsample count tables. When putting the value **0**, the tool will use minimal RRPIB (or MRPIB) among selected samples as desired value. Desired value of RRPIB (or MRPIB) must be **not greater** than minimal RRPIB (or MRPIB among samples).

### Output

Output of the tool is four text files:

1) **downsampled count tables** -- count tables downsampled up to desired RRPIB or MRPIB;

2) **downsampling plot** -- plot with downsampling curves fitted by NLS model. Point represent values downsampling statistics files used to fit models, solid black vertical line represents selected value of mean reads per barcode, dotted horizontal lines represent values of median transcript per barcode used for downsample count tables;

<img src="https://github.com/scipio-bioscience/eval_tools/blob/master/matrix_downsampling_tool/test_result/output/ds_plot_ds_1000_rpib1000.png" width="70%">

3) **summary check** -- tsv file with estimated by NLS median transcript count per barcode and real transcript count per barcode computed after downsampling. This file is destinated to check whether the median transcript count per barcode of downsampled count tables are similar to desired and not biased by barcode heterogeneity nether by randomizer bias;

*Example:*

| sample | estimated_transcript_median | real_transcript_median |
| --- | --- | --- |
| X081460_08_1_0.05 | 1944.1 |	1946.5 |
| X081460_08_2_0.05 |	2474.3 |	2476 |

4) **summary stats** -- tsv file with summary statistics of transcript count per barcode distribution and RPIB value for each sample before and after downsampling.

*Example:*

| sample | transcript_per_barcode_min |	transcript_per_barcode_Q1 |	transcript_per_barcode_median |	transcript_per_barcode_mean |	transcript_per_barcode_Q3 |	transcript_per_barcode_max |	transcript_per_barcode_IQR\* |	transcript_per_barcode_stdev |	transcript_per_barcode_CV\*\* | 	transcript_per_barcode_rIQR\*\*\* |	RPIB | transcript_per_barcode_state |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| X081460_08_1_0.05	| 952 |	1360.25 |	1946.5 | 2846.884 |	3095.5 | 56799 | 1735.25 | 4039.873 |	1.419 | 0.891 | 6767.276 | before downsampling |
| X081460_08_2_0.05 |	1057 | 1750.5 |	2491.5 | 3192.838 |	4108 | 13431 | 2357.5 |	2060.772 | 0.645 | 0.946 | 6818.228 | before downsampling |
| X081460_08_1_0.05 |	163 |	244.25 | 345 | 503.896 | 555.75 |	10085	| 311.5 |	718.74 | 1.408 | 0.89 | 1000 | after downsampling |
| X081460_08_2_0.05 |	163 |	293 |	410 |	533.203 |	701.25 | 2211 |	408.25 | 345.059 | 0.644 | 0.911 | 1000 | after downsampling |

<img src="https://render.githubusercontent.com/render/math?math=^*transcript\_per\_barcode\_IQR\ (interquartile\ range) = transcript\_per\_barcode\_Q3 - transcript\_per\_barcode\_Q1">

<img src="https://render.githubusercontent.com/render/math?math=^{**}transcript\_per\_barcode\_CV = \frac{transcript\_per\_barcode\_stdev}{transcript\_per\_barcode\_mean}">

<img src="https://render.githubusercontent.com/render/math?math=^{***}transcript\_per\_barcode\_rIQR = \frac{transcript\_per\_barcode\_IQR}{transcript\_per\_barcode\_median}">


### Launch tool

The tool is R script **matrix_downsampling_tool.R.** You can modify and launch shell file **launch_script.sh** or write command directly in the terminal. Before launching shell file, make sure that the file is allowed to execution.

The command looks like:

*Rscript [path_to_matrix_downsampling_tool.R] [config_file]*

When the tool is correctly run, you will see the message **"Mission completed."** in the terminal. If the parameter **desired_mean_reads_per_barcode** exceeds minimal RRPIB (or MRPIB) of selected samples, the execution will be stopped with the following message:

**"Parameter desired_RPIB (XXX) is higher than minimal RPIB among samples (YYY)."**, where **XXX** -- selected RPIB in config file, **YYY** -- minimal RPIB among selected samples. Also RPIB per sample for initial count tables will be shown. You can use this value or lower when rerunning **matrix downsampling tool.**

Another errors of execution may indicate wrong format of input files or non-existing pathways to files.

Libraries used in the tool: **dplyr**, **data.table**, **tibble**, **yaml**, **ggplot2**, **RColorBrewer**.


### General information

**Matrix downsampling tool** uses information from downsampling statistics files for fitting NLS model of downsampling curve (separately for each sample). Then it estimates the value of median transcript count correponding to desired RRPIB (or MRPIB) for each sample respectively. After that the ratio between estimated and initial median transcript count per barcode is computed for each sample and is used as downsampling coefficient for count table. So, computed proportion of transcripts is randomly selected from each count table for building downsampled count table.

You can use **matrix downsampling tool** for downsampling count tables either from one folder (one project) and from different folders (different projects). If you need to downsample count tables from different folders, you should specify list of folders containing input data in the parameters **count_table_folders**, **downsampling_folders**, **whitelist_folders**. The samples of interest are specified in the parameter **samples**. Notice that you can specify just some samples from each folder. The parameters **count_table_folders**, **downsampling_folders**, **whitelist_folders** have not to be necessarily ordered by sample order in the parameter **samples**. The sample in **downsampling plot** and **summary table** is specified by the parameter **samples**.


Don't forget to put **"/"** at the end of all folder names in the config file. 

* The pathway to input count tables is concatenated by config parameters:

*[count_table_folders] + [samples] + [count_prefix]*

*Example:*

The following parameters:

***count_table_folders:** /network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/*

***samples:** X081460_08_1_0.05*

***count_table_prefix:** \_merged_counts.tsv*

give full pathway to the input file:

*/network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/X081460_08_1_0.05_X081460_08_1_0.05_merged_counts.tsv*

* The pathway to input downsampling statistics files is concatenated by config parameters:

*[downsampling_folders] + [samples] + [downsampling_prefix]*

*Example:*

The following parameters:

***downsampling_folders:** /network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/*

***samples:** X081460_08_1_0.05*

***downsampling_prefix:** \_ds_merged_raw.txt*

give full pathway to the input file:

*/network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/X081460_08_1_0.05_X081460_08_1_0.05_ds_merged_raw.txt*

* The pathway to input whitelist files is concatenated by config parameters:

*[whitelist_folders] + [samples] + [whitelist_prefix]*

*Example:*

The following parameters:

***whitelist_folders:** /network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/*

***samples:** X081460_08_1_0.05*

***whitelist_prefix:** \_ds_merged_raw.txt*

give full pathway to the input file:

*/network/lustre/dtlake01/scipio/scripts/github_clone_repo/eval_tools/matrix_downsampling_tool/test_result/input/X081460_08_1_0.05_X081460_08_1_0.05_auto_whitelist_distance.txt*

* The pathway to output files is concatenated by config parameters:

**Count tables:** *[output_folder] + [samples] + [output_prefix] + [RPIB_prefix] + '\_counts.tsv'*

**Downsamping plot:** *[output_folder] + 'ds_plot' + [output_prefix] + [RPIB_prefix] + '.png'*

**Summary check:** *[output_folder] + 'summary_check' + [output_prefix] + [RPIB_prefix] + '.tsv'*

**Summary stats:** *[output_folder] + 'summary_stats' + [output_prefix] + [RPIB_prefix] + '.tsv'*

The folder **test_result** contains test input and output files of **matrix downsampling tool.** The corresponding config file is **config_matrix_downsampling_test.yml**. You can use this config file as template for customer using. The folder also contains another config file **config_matrix_downsampling_complexe_example.yml** with example of more complexe use case of downsampling tool when count tables from different folders (different projects) are needed to downsample.
