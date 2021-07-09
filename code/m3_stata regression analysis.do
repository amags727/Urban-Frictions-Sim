clear 
set more off 
global WD "/Volumes/GoogleDrive/My Drive/Github/Urban-Frictions-Sim/output/"
cd "/Volumes/GoogleDrive/My Drive/Github/Urban-Frictions-Sim/output/"

//Set the parameters
local d_shock_min = 90
local d_shock_max = 99
local d_shock_length =10

local xi_min =0
local xi_max= 100
local xi_length= 5

//Import and clean the results from the simulations without frictions
import delimited "no_friction_results.csv", numericcols(1) clear 
local varlist `"d_shock_09 d_shock_091 d_shock_092 d_shock_093 d_shock_094 d_shock_095 d_shock_096 d_shock_097  d_shock_098 "'
local shock = `d_shock_min'
foreach var in `varlist' d_shock_099{
	preserve
	keep og_distance og_popularity `var'
	rename `var' log_welfare_improvement
	gen d_shock = `shock'
	
	tempfile file_`var'
	save `file_`var'', replace
	restore
	local shock = `shock' + (`d_shock_max' - `d_shock_min')/ (`d_shock_length'-1)
} 
use `file_d_shock_099', clear 


foreach var in `varlist'{
	append using `file_`var''
}
gen xi_level=`xi_max'
tempfile no_friction
save `no_friction'

//Import and clean the results from the simulations with frictions
import delimited "friction_results.csv", numericcols(1) clear 
local varlist 
local xi= `xi_min'
local varlist `"xi_level_0 xi_level_025 xi_level_05 xi_level_075 "'

foreach var in `varlist' xi_level_1{
	preserve
	keep og_distance og_popularity `var'
	rename `var' log_welfare_improvement
	gen xi_level= `xi'
	local xi= `xi'+ (`xi_max'- `xi_min')/ (`xi_length'-1)
	tempfile file_`var'
	save `file_`var'', replace
	restore
}

use `file_xi_level_1', clear
foreach var in `varlist'{
	append using `file_`var''
}
gen d_shock=`d_shock_min'
append using `no_friction'

capture mkdir "stata output"
save "stata output/cleaned results.dta",replace 
rename og_distance  log_distance
rename og_popularity  log_popularity


//Conduct the regressions 
capture mkdir "stata output/regression tables"
reg log_welfare_improvement i.xi_level##c.log_distance
outreg2 using "stata output/regression tables", replace ctitle(Log welfare improvement)


