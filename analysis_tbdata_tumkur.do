/*****************************************************************
* Author: S. Mohan
* Affiliation: Sciences Po, Paris
* Last update: 21 Feb 2019
* Description: This dofile helps run necessary diagnostics to narrow down on an appropriate regression model to 
	predict treatment outcomes among TB cases in Tumkur city, Karnataka (India). This analysis was carried out for
	the master's thesis project - "Non-Biological Determinants of Treatment Success among Pulmonary TB cases in Tumkur city"
*********************************************************************/


*********************************************************************************
* 0. Preparation *
*********************************************************************************
* Set File paths
	global source "/Users/sakshimohan/Dropbox (Personal)/Personal/Thesis/Tumkur Data/TB Data/"
	global raw "$source/Final/"
	global output "$source/Final/output"
	global pc "/Users/sakshimohan/Dropbox (Personal)/Personal/Thesis/thesis/Spatial Analysis"
	global iph "$source/Final/raw/iph"
	global cleaning "$source/Final/cleaning"
	global tex "$source/Final/tex/"
	global temp "/Users/sakshimohan/Documents/spatial_exp/temp/"
		
* Housekeeping
	set matsize 1000
	set more off
	
* Load merged dataset
	use "$output/clean_data.dta", clear

* Drop children aged 5 or under
	drop if age <= 5
	
* Sort by time	
	sort elapsed_time tb_no
	
* Global control varlist
	global controlvars female lnage sputum_pos  previous_treat hiv 
	
* Look at the distribution of all dependent and control variables
	cd "$temp/"
	foreach var of varlist treat_success pp_ref intervention_dumreg $controlvars age lnage lndist_dotop lndist_ref lndist_dtc clusterb_750_66{
		hist `var', saving("`var'", replace) percent xtitle(,size(small)) 
		local list `list' `var'.gph
		}
		gr combine `list'
		graph export "$output/histograms_allvars.png", replace

* Generate Log file
		capture log close
		log using "$output/spatial_autocorrelation_with weights correction.log", replace
	
log on
/*********************************************************************************
* I. Run Spatial Diagnosis *
**********************************************************************************

Spatial diagnostics are run using three types of characterisations of spatial autocorrelation:
1. Based on Standardized Binary Matrix (Only Spatial)
2. Based on Standardized Distance Matrix (Only Spatial)
3. Based on Mahalanobis Distance Matrix (Both Space and Time)

The spatial weight matrices generated are later employed in the regression analysis
**********************************************************************************/
	
	* 1) Based on Standardized Binary Matrix (Only Spatial) - Band = [0,1250 meters]
	*----------------------------------------------------------------------------------
		
		// I try to generate the spatial weights matrics for the different variable combinations that 
		// can be used in the analysis. Since the number of observations missing for each set of variables
		// varies, a different weights matrix has to be defined for each variable combination
		
		// Size of different matrices based on missing observations for different variable combinations
			
		* Generate globals containing 4 sets of variables
			global miss_fullreg treat_success ref_source intervention_dumreg $controlvars dot_op_type_clean lndist_dotop lndist_ref lndist_dtc  x y pulmonary 
			global miss_norefsource treat_success $controlvars dot_op_type_clean x y
			global miss_nodist treat_success ref_source intervention_dumreg $controlvars dot_op_type_clean pulmonary
			global miss_outcome treat_success x y 
		
		* Generate a variable for each variable set denoting whether any variable in the set is missing
			foreach missingtype in miss_fullreg  miss_outcome miss_norefsource miss_nodist{
				cap drop `missingtype'
				egen `missingtype' = rowmiss($`missingtype')
			}
				
		* Generate spatial weights matrix and global spatial autocorrelation statistics for the 4 variable sets
			foreach missingtype in miss_fullreg miss_norefsource miss_nodist miss_outcome{
				preserve
				keep if `missingtype' == 0  
			
			// Weights matrix
				spatwmat, name(Wstandbin_`missingtype') xcoord(xproj) ycoord(yproj) band(0 1250) binary
			
			// Use nearest neighbour criteria for features with no neighbours
				* Generating and eigen vector requires every feature to have a neighbour. 
				* Eigen vector is required for spatial regression
						count
						local levels = r(N)
						
						* Identify elements with zero neighbours
						mata: st_matrix("test", rowsum(st_matrix("Wstandbin_`missingtype'")))
						local j 1
						
						forval i = 1/`levels'{
							if test[`i',1] == 0{
								local zeroelement`j' = `i'
								local j = `j'+1
							}
						}
						
						qui{
							noisily di "`j' patients have 0 neighbours"
							local k = `j' - 1
							forval i = 1/`k'{
								noisily di "Processing `zeroelement`i'' - `i'"
								local outlier = `zeroelement`i''
								di "`outlier'"
								mat distmat`outlier' = J(1,`levels',.)
								
							* Generate a matrix of distances from all other patients for each patient with zero neighbours
								forval l = 1/`levels' {
									mat distmat`outlier'[1,`l'] = sqrt((xproj[`i'] - xproj[`l'])^2 + (yproj[`i'] - yproj[`l'])^2)
									}
								
								* Replace distance from self as missing
								mat distmat`outlier'[1,`i'] = .
								
								* Find out the patient at the minimum distance
								mata: st_matrix("A`outlier'", rowmin(st_matrix("distmat`outlier'")))
								
								* Replace xij = 1 for the nearest neighbour for the no neighbour patients
								forval l = 1/`levels' {
									if distmat`outlier'[1,`l'] == A`outlier'[1,1] {
										di "Nearest neighbour - `l'"
										mat Wstandbin_`missingtype'[`outlier',`l'] = 1
										mat Wstandbin_`missingtype'[`l',`outlier'] = 1
										mat A = Wstandbin_`missingtype'[`outlier',`l'...]
									}
								}
							}
						}
				
			// Global Spatial Autocorrelation measures (Moran's I and Geary's c)
				spatgsa $`missingtype', weights(Wstandbin_`missingtype') moran geary
				
			// Spatial Correllelogram
				spatcorr treat_success, bands(250(250)2000) xcoord(xproj) ycoord(yproj) graph
				spatcorr treat_success, bands(250(250)2000) xcoord(xproj) ycoord(yproj) cumulative 
				
				restore
			
				matrix symeigen X_`missingtype' lambda_`missingtype' = Wstandbin_`missingtype'	
		} 
			
	
	log off	

		
	log on
	* 2) Based on Standardized Distance Matrix (Only Spatial) - Band = [0,1250 meters]
	*----------------------------------------------------------------------------------

		// Size of different matrices based on missing observations for different variable combinations	
		// Note: the same globals as under 1) are used 
			
		* Generate spatial weights matrix and global spatial autocorrelation statistics for the 4 variable sets
			foreach missingtype in miss_fullreg miss_norefsource miss_nodist miss_outcome{
				preserve
				keep if `missingtype' == 0 
			
			// Weights matrix
				spatwmat, name(Wstanddis_`missingtype') xcoord(xproj) ycoord(yproj) band(0 1250) eigenval(eigen_`missingtype')
			
			// Global Spatial Autocorrelation measures (Moran's I and Geary's c)
				spatgsa $`missingtype', weights(Wstanddis_`missingtype') moran geary
				
			// Spatial Correllelogram
				spatcorr treat_success, bands(250(250)2000) xcoord(xproj) ycoord(yproj)
				
				restore
			} 
			
	log off	
	
	log on
	
	* 3) Based on Mahalanobis Distance Matrix (Both Space and Time)
	*----------------------------------------------------------------------------------
		* Generate spatial weights matrix and global spatial autocorrelation statistics for the 4 variable sets
			foreach missingtype in miss_fullreg miss_norefsource miss_nodist miss_outcome{
				preserve
				keep if `missingtype' == 0

				* Use projected coordinates and number of days since registration for treatment 
					local maha_vars xproj yproj elapsed_time
				
				* Find the central value of each variable
					sum `maha_vars'
					foreach i of local `maha_vars'{
						qui sum `i'
						cap drop `i'cen
						gen `i'cen = `i' - r(mean)
					}
				
				* Generate Mahalanobis Inverse Distance Matrix
					matrix accum SSCP = xproj yproj elapsed_time, nocons dev
					matrix S = SSCP/(r(N) - 1)
					matrix Sinv = invsym(S)
					mkmat xprojcen yprojcen elapsed_timecen, matrix(Xcen)
					matrix Xcentran = Xcen'
					matrix mahaldist_`missingtype' = Xcen * Sinv * Xcentran
					matrix distsvec = vecdiag(mahaldist_`missingtype')
					mat mahalanobis = distsvec'
					svmat mahalanobis, name(distances_`missingtype')
					
				* Export matrix in a csv file
					putexcel A1 = matrix(mahaldist_`missingtype') using "$output/mahadist_`missingtype'.csv", replace
				
				restore
			}
			
			* Save the above matrix as a datafile to be import by the spatwmat command
				foreach missingtype in miss_fullreg miss_norefsource miss_nodist miss_outcome{
					preserve
						import excel using "$output/mahadist_`missingtype'.csv", clear
						save "$output/mahadist_`missingtype'.dta", replace
					restore
				}
				
			* Run spatial diagnostics
				foreach missingtype in miss_fullreg miss_norefsource miss_nodist miss_outcome{
					preserve
					keep if `missingtype' == 0 
						// Import mahalanobis distance weights matrix
						spatwmat using "$output/mahadist_`missingtype'.dta", name(mahadist_`missingtype')  // xcoord(xproj) ycoord(yproj) band(0 2000) standardize 
						matrix symeigen eigen_`missingtype' lambdamaha_`missingtype' = mahadist_`missingtype'
						
						// Global Index of Spatial Autocorrelation
						spatgsa $`missingtype', weights(mahadist_`missingtype') moran geary
						
						// Spatial Correllelogram 
						*spatcorr treat_success, bands(250(250)2000) xcoord(xproj) ycoord(yproj) cumulative 
						*spatcorr treat_success, bands(250(250)2000) xcoord(xproj) ycoord(yproj)
			
					restore
				}
			
		log off
		
		* Edit mahalanobis distance matrix values to 0 if the spatial distance between two cases is greater than 750 meters
			
			mat mahadist_miss_fullreg_adj = mahadist_miss_fullreg
			
			local changes = 0
			forval i = 1/402{
				forval j = 1/402{
					local a = Wstandbin_miss_fullreg[`i',`j']
					if `a'== 0{
					di "mat[`i',`j'] changed"
					mat mahadist_miss_fullreg_adj[`i',`j'] = 0
					mat mahadist_miss_fullreg_adj[`j',`i'] = 0
					local changes = `changes'+1
					}
					else{
					 di "mat[`i',`j'] not changed"
					}
				}
			}
			di "`changes' out of 403*403 elements changed to 0"
			
			count
			local levels = 403
						
			* Identify elements with zero neighbours
				mata: st_matrix("test", rowsum(st_matrix("mahadist_miss_fullreg_adj")))
				local j 1
				
				forval i = 1/`levels'{
					if test[`i',1] == 0{
						local zeroelement`j' = `i'
						local j = `j'+1
					}
				}
				
				qui{
					noisily di "`j' patients have 0 neighbours"
					local k = `j' - 1
					forval i = 1/`k'{
						noisily di "Processing `zeroelement`i'' - `i'"
						local outlier = `zeroelement`i''
						di "`outlier'"
						mat distmat`outlier' = mahadist_miss_fullreg[`outlier',1..403]
						
						* Find out the patient at the minimum distance
						mata: st_matrix("A`outlier'", rowmin(st_matrix("distmat`outlier'")))
						
						* Replace xij = 1 for the nearest neighbour for the no neighbour patients
						forval l = 1/`levels' {
							if distmat`outlier'[1,`l'] == A`outlier'[1,1] {
								di "Nearest neighbour - `l'"
								mat mahadist_miss_fullreg_adj[`outlier',`l'] = mahadist_miss_fullreg[`outlier',`l']
								mat mahadist_miss_fullreg_adj[`l',`outlier'] = mahadist_miss_fullreg[`l',`outlier']
								mat A = mahadist_miss_fullreg_adj[`outlier',`l'...]
								noisily mat list A
							}
						}
					}
					}
		
	
/*********************************************************************************
* I. Test Regression Models  *
**********************************************************************************

Four types of regression models are analysed:
A. Non-spatial regression
B. Spatial Regression With Spatial Weights Matrix
C. Spatial Regression With Mahalanobis Weights Matrix
D. Space-time Scan Clusters

Wald tests are carried out for 5 models under each of the above (consisting of an incremental set of regressors). 
Results of this analysis will be used to determine the most appropriate model to be used for the final analysis. 
**********************************************************************************/	

* Set locals to indicate if the winsorized and logged versions of veriable should be used
* (To use winsorized variables - "w", to used logged variables - "ln")
* Note: These are applied to all 4 types of models 
local winsor = ""
local ln = "ln"
local clustervar clusterb_750_66

	* A. Non-spatial regression
	**********************************************************************************	
	* i) Linear regression (Preliminary)
		*-------------------------------
						
			global prelim i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean
			global lndist_nodtc i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean lndist_dotop lndist_ref			
			global lndist i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean lndist_dotop lndist_ref lndist_dtc
			global intlndist_nodtc i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean female##(c.lndist_dotop c.lndist_ref) c.lnage##(c.lndist_dotop c.lndist_ref) i.ref_source##(c.lndist_dotop c.lndist_ref)			
			global intlndist i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean female##(c.lndist_dotop c.lndist_ref c.lndist_dtc) c.lnage##(c.lndist_dotop c.lndist_ref c.lndist_dtc) i.ref_source##(c.lndist_dotop c.lndist_ref c.lndist_dtc)
			
				local test_count 14
				local varset_count 5
				local row_count `test_count' 
				local j 1
				
				mat fullmat = J(`row_count',`varset_count'*2,.)
				
			foreach varset in prelim lndist intlndist{

				
				reg treat_success $`varset' if class == "P", robust
				est sto `varset'
				
				cap testparm lndist*
				estadd scalar distance_testparm_p = r(p) 
				cap testparm i.ref_source intervention_dum
				estadd scalar refsource_testparm_p = r(p)
				cap testparm female##(c.lndist_ref c.lndist_dotop c.lndist_dtc) c.lnage##(c.lndist_ref c.lndist_dotop c.lndist_dtc) i.ref_source##(c.lndist_ref c.lndist_dotop c.lndist_dtc)
				estadd scalar distint_testparm_p = r(p)
				
				* List of Tests
				*----------------
					* General variable combinations 
						local test1 i.ref_source
						local test2 i.ref_source intervention_dumreg
						local test3 $controlvars
						local test4 i.dot_op_type_clean
						local test5 i.ref_source i.dot_op_type_clean
						local test6 i.ref_source intervention_dumreg i.dot_op_type_clean
						local test7 i.ref_source $controlvars
					* Distance variable combinations
						local test8 `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test9 $controlvars `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test10 i.ref_source intervention_dumreg `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
					* Distance variable combinations with interactions
						local test11 female female##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test12 c.lnage##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test13 i.ref_source##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')						
						local test14 female female##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor') c.lnage##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor') i.ref_source##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')						
						
						* Fill table of results
					forval i = 1/`test_count'{
						cap testparm `test`i''
						mat fullmat[`i',`j'] = r(p)
						mat fullmat[`i',`j'+1] = r(df)
						}
						local j `j'+2
					
				}

					mat rownames fullmat = "1 Ref Source" "2 Ref Source + Intervention" "3 Demo \& Bio Controls" "4 DOTS" "5 Ref + DOTS" "6 Ref + Intervention + DOTS" "7 Ref + Demo \& Bio Controls" "8 Distances" "9 Dist + Bio \& Dem Controls" "10 Dist + Ref + Intervention" ///
						"11 Female ** Distances" "12 Age ** Distances" "13 Ref Source ** Distances" "14 All Interactions with Dist"
	
					frmttable using "$tex/tables/testparm_results_full.tex", statmat(fullmat) varlabels replace ///
						ctitle("", p-value(i), df(i), p-value(ii), df(ii), p-value(iii), df(iii), p-value(iv), df(iv), p-value(v), df(v)) tex ///
						hlines(1) vlines(110101010101) 
					
					esttab prelim lndist_nodtc lndist intlndist_nodtc intlndist  using "$tex/tables/linear_nonspatial.tex", replace title(Linear Regressions - Without Spatial Correction) p scalars(N) r2 ///
						mtitle("Preliminary(i)" "Distances I(ii)" "Distances II(iii)" "Interactions I(iv)" "Interactions II(v)") star(* .10 ** .05 *** .01) nobaselevels  ///
						label noconstant addnote("Note:distint_testparm_p tests the joint significance of all the interactions with distance") stats(r2 N distance_testparm_p refsource_testparm_p distint_testparm_p) ///
						coeflabels(1.dot_op_type_clean "DOT Operator: AWW" 2.dot_op_type_clean  "DOT Operator: Asha" 3.dot_op_type_clean "DOT Operator: Community DOT Provider" 4.dot_op_type_clean "DOT Operator: PHC" 5.dot_op_type_clean "DOT Operator: PP" 6.dot_op_type_clean "DOT Operator: TBHV" 2.ref_source "Referred by PP(Base:PHC)" 3.ref_source "Self-referred(Base:PHC)")  ///
						order(2.ref_source 3.ref_source intervention_dumreg female lnage sputum_pos previous_treatment hiv 1.dot_op_type_clean 2.dot_op_type_clean 3.dot_op_type_clean 4.dot_op_type_clean 5.dot_op_type_clean lndist_ref lndist_dotop lndist_dtc 1.female#lndist_ref 1.female#lndist_dotop 1.female#lndist_dtc lnage#lndist_ref lnage#lndist_dotop lnage#lndist_dtc)
		 
		 
	* B. Spatial Regression With Spatial Weights Matrix
	**********************************************************************************	
			global prelim i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean
			global lndist_nodtc i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean lndist_dotop lndist_ref			
			global lndist i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean lndist_dotop lndist_ref lndist_dtc
			global intlndist_nodtc i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean intlndist_dotop_female intlndist_ref_female intlndist_dotop_lnage intlndist_ref_lnage // intlndist_dotop_refpp intlndist_ref_refpp intlndist_dotop_refself intlndist_ref_refself  			
			//global intlndist i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean intlndist_dotop_female intlndist_ref_female intlndist_dtc_female intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage // intlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself 
			
				local test_count 14
				local varset_count 4
				local row_count `test_count'
				local j 1
				
				mat fullmat_spat = J(`row_count',`varset_count'*2,.)
			
			foreach varset in prelim lndist_nodtc lndist intlndist_nodtc{ // intlndist

				preserve
				
				if "`varset'" == "prelim"{
					keep if miss_nodist == 0
					xi: spatreg treat_success $`varset', weights(Wstandbin_miss_nodist) eigenval(X_miss_nodist) model(error) robust
					est sto sp`varset'
				}
				else{
					keep if miss_fullreg == 0 & class == "P"
					xi: spatreg treat_success $`varset', weights(Wstandbin_miss_fullreg) eigenval(X_miss_fullreg) model(error) robust
					est sto sp`varset'
				}
				

				
				* List of Tests
				*---------------
					
					* General variable combinations 
						local test1 i.ref_source
						local test2 i.ref_source intervention_dumreg
						local test3 $controlvars
						local test4 i.dot_op_type_clean
						local test5 i.ref_source i.dot_op_type_clean
						local test6 i.ref_source intervention_dumreg i.dot_op_type_clean
						local test7 i.ref_source $controlvars
					* Distance variable combinations
						local test8 `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test9 $controlvars `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test10 i.ref_source intervention_dumreg `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
					* Distance variable combinations with interactions
						local test11 female `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_female intlndist_ref_female intlndist_dtc_female
						local test12 lnage `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage 
						local test13 i.ref_source i`ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' ntlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself
						local test14 female lnage i.ref_source `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_female intlndist_ref_female intlndist_dtc_female intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage intlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself						
						
						* Fill table of results
					forval i = 1/`test_count'{
						cap testparm `test`i''
						mat fullmat_spat[`i',`j'] = r(p)
						mat fullmat_spat[`i',`j'+1] = r(df)
						}
						local j `j'+2
					
					restore
				}

					mat rownames fullmat_spat = "1 Ref Source" "2 Ref Source + Intervention" "3 Demo \& Bio Controls" "4 DOTS" "5 Ref + DOTS" "6 Ref + Intervention + DOTS" "7 Ref + Demo \& Bio Controls" "8 Distances" "9 Dist + Bio \& Dem Controls" "10 Dist + Ref + Intervention" ///
						"11 Female ** Distances" "12 Age ** Distances" "13 Ref Source ** Distances" "14 All Interactions with Dist"
	
					frmttable using "$tex/tables/testparm_results_full_spatial.tex", statmat(fullmat_spat) varlabels replace ///
						ctitle("", p-value(i), df(i), p-value(ii), df(ii), p-value(iii), df(iii), p-value(iv), df(iv)) tex ///
						hlines(1) vlines(1101010101) 
					
					esttab spprelim splndist_nodtc splndist spintlndist_nodtc  using "$tex/tables/linear_spatial.tex", replace title(Linear Regressions - With Spatial Correction) p scalars(N) r2 ///
						mtitle("Preliminary(i)" "Distances I(ii)"  "Distances II(iii)" "Interactions I(iv)") star(* .10 ** .05 *** .01) drop(_cons) ///
						label noconstant addnote() ///
						coeflabels(1.dot_op_type_clean "DOT Operator: AWW" _Idot_op_ty_2  "DOT Operator: Asha" _Idot_op_ty_3 "DOT Operator: Community DOT Provider" _Idot_op_ty_4 "DOT Operator: PHC" _Idot_op_ty_5 "DOT Operator: PP" _Idot_op_ty_6 "DOT Operator: TBHV" _Iref_sourc_2 "Referred by PP(Base:PHC)" _Iref_sourc_3 "Self-referred(Base:PHC)")  ///
						order(_Iref_sourc_2 _Iref_sourc_3 intervention_dumreg female lnage sputum_pos previous_treatment hiv _Idot_op_ty_2 _Idot_op_ty_3 _Idot_op_ty_4 _Idot_op_ty_5 _Idot_op_ty_6 lndist_ref lndist_dotop lndist_dtc int*female int*lnage)

		
		
	* C. Spatial Regression With Mahalanobis Weights Matrix
	**********************************************************************************	
	
			global prelim i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean
			global lndist_nodtc i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean lndist_dotop lndist_ref			
			global lndist i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean lndist_dotop lndist_ref lndist_dtc
			global intlndist_nodtc i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean intlndist_dotop_female intlndist_ref_female intlndist_dotop_lnage intlndist_ref_lnage // intlndist_dotop_refpp intlndist_ref_refpp intlndist_dotop_refself intlndist_ref_refself  			
			//global intlndist i.ref_source intervention_dumreg $controlvars i.dot_op_type_clean intlndist_dotop_female intlndist_ref_female intlndist_dtc_female intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage // intlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself 
			
				local test_count 14
				local varset_count 4
				local row_count `test_count' 
				local j 1
				
				mat fullmat_maha = J(`row_count',`varset_count'*2,.)
				
		foreach varset in prelim lndist_nodtc lndist intlndist_nodtc{ // intlndist

				preserve
				
				if "`varset'" == "prelim"{
					keep if miss_nodist == 0
					xi: spatreg treat_success $`varset', weights(mahadist_miss_nodist) eigenval(eigen_miss_nodist) model(error) robust
					est sto sp`varset'_maha
				}
				else{
					keep if miss_fullreg == 0 & class == "P"
					xi: spatreg treat_success $`varset', weights(mahadist_miss_fullreg) eigenval(eigen_miss_fullreg) model(error) robust
					est sto sp`varset'_maha
				}
				

				
				* List of Tests
				*----------------
					
					* General variable combinations 
						local test1 i.ref_source
						local test2 i.ref_source intervention_dumreg
						local test3 $controlvars
						local test4 i.dot_op_type_clean
						local test5 i.ref_source i.dot_op_type_clean
						local test6 i.ref_source intervention_dumreg i.dot_op_type_clean
						local test7 i.ref_source $controlvars
					* Distance variable combinations
						local test8 `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test9 $controlvars `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test10 i.ref_source intervention_dumreg `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
					* Distance variable combinations with interactions
						local test11 female `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_female intlndist_ref_female intlndist_dtc_female
						local test12 lnage `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage 
						local test13 i.ref_source i`ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' ntlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself
						local test14 female lnage i.ref_source `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' intlndist_dotop_female intlndist_ref_female intlndist_dtc_female intlndist_dotop_lnage intlndist_ref_lnage intlndist_dtc_lnage intlndist_dotop_refpp intlndist_ref_refpp intlndist_dtc_refpp intlndist_dotop_refself intlndist_ref_refself intlndist_dtc_refself						
						
						* Fill table of results
					forval i = 1/`test_count'{
						cap testparm `test`i''
						mat fullmat_maha[`i',`j'] = r(p)
						mat fullmat_maha[`i',`j'+1] = r(df)
						}
						local j `j'+2
					
					restore
				}

					mat rownames fullmat_maha = "1 Ref Source" "2 Ref Source + Intervention" "3 Demo \& Bio Controls" "4 DOTS" "5 Ref + DOTS" "6 Ref + Intervention + DOTS" "7 Ref + Demo \& Bio Controls" "8 Distances" "9 Dist + Bio \& Dem Controls" "10 Dist + Ref + Intervention" ///
						"11 Female ** Distances" "12 Age ** Distances" "13 Ref Source ** Distances" "14 All Interactions with Dist"
	
					frmttable using "$tex/tables/testparm_results_full_spatial_maha.tex", statmat(fullmat_maha) varlabels replace ///
						ctitle("", "", p-value(i), df(i), p-value(ii), df(ii), p-value(iii), df(iii), p-value(iv), df(iv)) tex ///
						hlines(1) vlines(110101010101) 
						
					esttab spprelim_maha splndist_nodtc_maha splndist_maha spintlndist_nodtc_maha  using "$tex/tables/linear_spatial_maha.tex", replace title(Linear Regressions - With Spatial Correction) p scalars(N) r2 ///
						mtitle("Preliminary(i)" "Distances I(ii)"  "Distances II(iii)" "Interactions I(iv)") star(* .10 ** .05 *** .01) drop(_cons) ///
						label noconstant addnote() ///
						coeflabels(1.dot_op_type_clean "DOT Operator: AWW" _Idot_op_ty_2  "DOT Operator: Asha" _Idot_op_ty_3 "DOT Operator: Community DOT Provider" _Idot_op_ty_4 "DOT Operator: PHC" _Idot_op_ty_5 "DOT Operator: PP" _Idot_op_ty_6 "DOT Operator: TBHV" _Iref_sourc_2 "Referred by PP(Base:PHC)" _Iref_sourc_3 "Self-referred(Base:PHC)")  ///
						order(_Iref_sourc_2 _Iref_sourc_3 intervention_dumreg female lnage sputum_pos previous_treatment hiv _Idot_op_ty_2 _Idot_op_ty_3 _Idot_op_ty_4 _Idot_op_ty_5 _Idot_op_ty_6 lndist_ref lndist_dotop lndist_dtc int*female int*lnage)
			
	
	* D) Space-time Scan Clusters
	*-----------------------------------------------------

			global prelim i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean `clustervar'##(female c.lnage previous_treatment hiv sputum_pos intervention_dum i.ref_source)
			global lndist_nodtc i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean lndist_dotop lndist_ref	`clustervar'		
			global lndist i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean lndist_dotop lndist_ref lndist_dtc `clustervar'
			global intlndist_nodtc i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean female##(c.lndist_dotop c.lndist_ref) c.lnage##(c.lndist_dotop c.lndist_ref) i.ref_source##(c.lndist_dotop c.lndist_ref)	clusterb_750_66 `clustervar'##(c.lndist_ref c.lndist_dotop)		
			global intlndist i.ref_source intervention_dumreg $controlvars ib6.dot_op_type_clean female##(c.lndist_dotop c.lndist_ref c.lndist_dtc) c.lnage##(c.lndist_dotop c.lndist_ref c.lndist_dtc) i.ref_source##(c.lndist_dotop c.lndist_ref c.lndist_dtc) `clustervar'##(c.lndist_ref c.lndist_dotop c.lndist_dtc) `clustervar'##(female c.lnage previous_treatment hiv sputum_pos intervention_dum i.ref_source)
			
				local test_count 16
				local varset_count 5
				local row_count `test_count'  // *`varset_count' + `varset_count'
				local j 1
				
				mat fullmat = J(`row_count',`varset_count'*2,.)
				
			foreach varset in prelim lndist_nodtc lndist intlndist_nodtc intlndist{

				
				reg treat_success $`varset' if class == "P", robust
				est sto `varset'
				
				cap testparm lndist*
				estadd scalar distance_testparm_p = r(p) 
				cap testparm i.ref_source intervention_dum
				estadd scalar refsource_testparm_p = r(p)
				cap testparm female##(c.lndist_ref c.lndist_dotop c.lndist_dtc) c.lnage##(c.lndist_ref c.lndist_dotop c.lndist_dtc) i.ref_source##(c.lndist_ref c.lndist_dotop c.lndist_dtc) clusterb_750_66##(c.lndist_ref c.lndist_dotop c.lndist_dtc)
				estadd scalar distint_testparm_p = r(p)
				cap testparm clusterb_750_66##(c.lndist_ref c.lndist_dotop c.lndist_dtc)
				estadd scalar cluster_testparm_p = r(p)
				
				* List of Tests
				*----------------
					
					* General variable combinations 
						local test1 i.ref_source
						local test2 i.ref_source intervention_dumreg
						local test3 $controlvars
						local test4 i.dot_op_type_clean
						local test5 i.ref_source i.dot_op_type_clean
						local test6 i.ref_source intervention_dumreg i.dot_op_type_clean
						local test7 i.ref_source $controlvars
					* Distance variable combinations
						local test8 `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test9 $controlvars `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
						local test10 i.ref_source intervention_dumreg `ln'dist_dtc`winsor' `ln'dist_ref`winsor'  `ln'dist_dotop`winsor' 
					* Distance variable combinations with interactions
						local test11 female female##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test12 c.lnage##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test13 i.ref_source##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test14 female female##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor') c.lnage##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor') i.ref_source##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')	
						local test15 `clustervar'##(c.`ln'dist_dtc`winsor' c.`ln'dist_ref`winsor'  c.`ln'dist_dotop`winsor')
						local test16 `clustervar'##($controlvars)
						
						* Fill table of results
					forval i = 1/`test_count'{
						cap testparm `test`i''
						mat fullmat[`i',`j'] = r(p)
						mat fullmat[`i',`j'+1] = r(df)
						}
						local j `j'+2
					
				}

					mat rownames fullmat = "1 Ref Source" "2 Ref Source + Intervention" "3 Demo \& Bio Controls" "4 DOTS" "5 Ref + DOTS" "6 Ref + Intervention + DOTS" "7 Ref + Demo \& Bio Controls" "8 Distances" "9 Dist + Bio \& Dem Controls" "10 Dist + Ref + Intervention" ///
						"11 Female ** Distances" "12 Age ** Distances" "13 Ref Source ** Distances" "14 All Interactions with Dist" "15 Cluster X Distances" "16 Cluster X Controls"
	
					 frmttable using "$tex/tables/testparm_results_full.tex", statmat(fullmat) varlabels replace ///
						ctitle("", p-value(i), df(i), p-value(ii), df(ii), p-value(iii), df(iii), p-value(iv), df(iv), p-value(v), df(v)) tex ///
						hlines(1) vlines(110101010101) 
					
					esttab prelim lndist_nodtc lndist intlndist_nodtc intlndist  using "$tex/tables/linear_scancluster.tex", replace title(Linear Regressions - Without Spatial Correction) p scalars(N) r2 ///
						mtitle("Preliminary(i)" "Distances I(ii)" "Distances II(iii)" "Interactions I(iv)" "Interactions II(v)") star(* .10 ** .05 *** .01) nobaselevels  ///
						label noconstant addnote("Note:distint_testparm_p tests the joint significance of all the interactions with distance") stats(r2 N distance_testparm_p refsource_testparm_p distint_testparm_p cluster_testparm_p) ///
						coeflabels(1.dot_op_type_clean "DOT Operator: AWW" 2.dot_op_type_clean  "DOT Operator: Asha" 3.dot_op_type_clean "DOT Operator: Community DOT Provider" 4.dot_op_type_clean "DOT Operator: PHC" 5.dot_op_type_clean "DOT Operator: PP" 6.dot_op_type_clean "DOT Operator: TBHV" 2.ref_source "Referred by PP(Base:PHC)" 3.ref_source "Self-referred(Base:PHC)")  ///
						order(2.ref_source 3.ref_source intervention_dumreg female lnage sputum_pos previous_treatment hiv 1.dot_op_type_clean 2.dot_op_type_clean 3.dot_op_type_clean 4.dot_op_type_clean 5.dot_op_type_clean lndist_ref lndist_dotop lndist_dtc 1.female#lndist_ref 1.female#lndist_dotop 1.female#lndist_dtc lnage#lndist_ref lnage#lndist_dotop lnage#lndist_dtc)
		 

					esttab prelim lndist_nodtc lndist intlndist_nodtc intlndist  ,  p scalars(N) r2 ///
						mtitle("Preliminary(i)" "Distances I(ii)" "Distances II(iii)" "Interactions I(iv)" "Interactions II(v)") star(* .10 ** .05 *** .01) nobaselevels  ///
						label noconstant addnote("Note:distint_testparm_p tests the joint significance of all the interactions with distance") stats(r2 N distance_testparm_p refsource_testparm_p distint_testparm_p cluster_testparm_p) 
						
**************************************************
* Edit all summary files for LaTeX Compatibility *
**************************************************

cd "$tex/tables/"
 
local files : dir . files "testparm_*.tex"

foreach texfile of local files{
	filefilter "`texfile'" "new_`texfile'", from("\BSdocumentclass[]{article}\BSpagestyle{empty}\BSbegin{document}") to("") replace
	filefilter "new_`texfile'" "`texfile'", from("\BSend{document}") to("") replace
	rm "new_`texfile'"
	sleep 500
	}	

