version 17
clear all
set more off

tempfile results
tempname handle
postfile `handle' str12 estimator str12 statistic str32 case str100 term ///
    double estimate std_error using `results', replace

capture program drop PostCurrentEstimates
program define PostCurrentEstimates
    syntax, Handle(name) Estimator(string) Statistic(string) Case(string)
    tempname coefficients covariance
    matrix `coefficients' = e(b)
    matrix `covariance' = e(V)
    local coefficient_names : colfullnames `coefficients'
    forvalues column = 1/`=colsof(`coefficients')' {
        local term : word `column' of `coefficient_names'
        post `handle' (`"`estimator'"') (`"`statistic'"') ///
            (`"`case'"') (`"`term'"') ///
            (`coefficients'[1, `column']) ///
            (sqrt(`covariance'[`column', `column']))
    }
end

quietly use "data/cattaneo2.dta", clear

quietly teffects ra (bweight prenatal1 mmarried mage fbaby) (mbsmoke)
post `handle' ("ra") ("ate") ("default") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])
post `handle' ("ra") ("ate") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects ra (bweight prenatal1 mmarried mage fbaby) (mbsmoke), atet
post `handle' ("ra") ("att") ("default") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])
post `handle' ("ra") ("att") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects ra (bweight prenatal1 mmarried mage fbaby) (mbsmoke), pomeans
post `handle' ("ra") ("pomeans") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])
post `handle' ("ra") ("pomeans") ("default") ("treated_pom") ///
    (_b[POmean:1.mbsmoke]) (_se[POmean:1.mbsmoke])

quietly teffects ipw (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit)
post `handle' ("ipw") ("ate") ("default") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])
post `handle' ("ipw") ("ate") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects ipw (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), atet
post `handle' ("ipw") ("att") ("default") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])
post `handle' ("ipw") ("att") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects ipwra (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit)
post `handle' ("ipwra") ("ate") ("default") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])
post `handle' ("ipwra") ("ate") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects ipwra (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), pomeans
post `handle' ("ipwra") ("pomeans") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])
post `handle' ("ipwra") ("pomeans") ("default") ("treated_pom") ///
    (_b[POmean:1.mbsmoke]) (_se[POmean:1.mbsmoke])

quietly teffects aipw (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit)
post `handle' ("aipw") ("ate") ("default") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])
post `handle' ("aipw") ("ate") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

quietly teffects aipw (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), pomeans
post `handle' ("aipw") ("pomeans") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])
post `handle' ("aipw") ("pomeans") ("default") ("treated_pom") ///
    (_b[POmean:1.mbsmoke]) (_se[POmean:1.mbsmoke])

quietly teffects aipw (bweight fbaby mage mmarried prenatal1) ///
    (mbsmoke fbaby foreign medu mmarried, probit), atet
post `handle' ("aipw") ("att") ("default") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])
post `handle' ("aipw") ("att") ("default") ("control_pom") ///
    (_b[POmean:0.mbsmoke]) (_se[POmean:0.mbsmoke])

/* Nearest-neighbor matching. cattaneo3 adds magecat to cattaneo2. */
quietly use "data/cattaneo3.dta", clear

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke)
post `handle' ("nnmatch") ("ate") ("default_ate") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke), atet
post `handle' ("nnmatch") ("att") ("default_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke), vce(iid)
post `handle' ("nnmatch") ("ate") ("iid_ate") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke), atet vce(iid)
post `handle' ("nnmatch") ("att") ("iid_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke), nneighbor(3)
post `handle' ("nnmatch") ("ate") ("neighbors3") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch ///
    (bweight mage prenatal1 mmarried fbaby) (mbsmoke), ///
    vce(robust, nn(4))
post `handle' ("nnmatch") ("ate") ("vce_neighbors4") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch (bweight mage) (mbsmoke), ///
    ematch(prenatal1 mmarried fbaby) metric(euclidean)
post `handle' ("nnmatch") ("ate") ("exact_euclidean") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch (bweight mage fage) (mbsmoke), ///
    ematch(prenatal1 mmarried fbaby) biasadj(mage fage)
post `handle' ("nnmatch") ("ate") ("bias_adjusted") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects nnmatch (bweight mage) (mbsmoke), ///
    metric(euclidean) caliper(50)
post `handle' ("nnmatch") ("ate") ("caliper50") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

/* Propensity-score matching. */
quietly use "data/cattaneo2.dta", clear

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu)
post `handle' ("psmatch") ("ate") ("logit_ate") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), atet
post `handle' ("psmatch") ("att") ("logit_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), vce(iid)
post `handle' ("psmatch") ("ate") ("logit_iid_ate") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), atet vce(iid)
post `handle' ("psmatch") ("att") ("logit_iid_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit)
post `handle' ("psmatch") ("ate") ("probit_ate") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), atet
post `handle' ("psmatch") ("att") ("probit_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), nneighbor(4)
post `handle' ("psmatch") ("ate") ("neighbors4") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), ///
    vce(robust, nn(4))
post `handle' ("psmatch") ("ate") ("vce_neighbors4") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), caliper(0.1)
post `handle' ("psmatch") ("ate") ("caliper_0_1") ("effect") ///
    (_b[ATE:r1vs0.mbsmoke]) (_se[ATE:r1vs0.mbsmoke])

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), ///
    atet vce(iid) caliper(0.03)
post `handle' ("psmatch") ("att") ("caliper_0_03_att") ("effect") ///
    (_b[ATET:r1vs0.mbsmoke]) (_se[ATET:r1vs0.mbsmoke])

/* telasso */
quietly telasso ///
    (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried mage c.mage#c.mage fbaby medu), ///
    selection(plugin) nolog
PostCurrentEstimates, handle(`handle') estimator(telasso) statistic(ate) ///
    case(default)

/*
Manual example census. These commands reproduce every displayed teffects
command in the dedicated manual entries.
All coefficients and standard errors are posted under stable manual_* cases.
*/

/* teffects intro: six-estimator tour */
quietly use "data/cattaneo2.dta", clear

quietly teffects ra (bweight mmarried mage prenatal1 fbaby) (mbsmoke)
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_intro_ra)

quietly teffects ipw (bweight) ///
    (mbsmoke mmarried mage prenatal1 fbaby, probit)
PostCurrentEstimates, handle(`handle') estimator(ipw) statistic(ate) ///
    case(manual_intro_ipw)

quietly teffects ipwra (bweight mmarried mage prenatal1 fbaby) ///
    (mbsmoke mmarried mage prenatal1 fbaby)
PostCurrentEstimates, handle(`handle') estimator(ipwra) statistic(ate) ///
    case(manual_intro_ipwra)

quietly teffects aipw (bweight mmarried mage prenatal1 fbaby) ///
    (mbsmoke mmarried mage prenatal1 fbaby)
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_intro_aipw)

quietly teffects nnmatch (bweight mage fage) (mbsmoke), ///
    ematch(prenatal1 mmarried fbaby) biasadj(mage fage)
PostCurrentEstimates, handle(`handle') estimator(nnmatch) statistic(ate) ///
    case(manual_intro_nnmatch)

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried mage prenatal1 fbaby, probit)
PostCurrentEstimates, handle(`handle') estimator(psmatch) statistic(ate) ///
    case(manual_intro_psmatch)

/* teffects aipw */
quietly use "data/cattaneo2.dta", clear

quietly teffects aipw (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), ///
    pomeans aequations
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(pomeans) ///
    case(manual_aipw_2)

quietly teffects aipw (bweight prenatal1 mmarried fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, ///
    hetprobit(c.mage)), aequations
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_aipw_3)

quietly teffects aipw (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), wnls
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_aipw_4)

quietly teffects aipw (bweight fbaby mage mmarried prenatal1) ///
    (mbsmoke fbaby foreign medu mmarried, probit), atet
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(att) ///
    case(manual_aipw_5)

quietly teffects aipw (bweight fbaby mage mmarried prenatal1) ///
    (msmoke fbaby foreign medu mmarried), atet
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(att) ///
    case(manual_aipw_6)

/* teffects ipw */
quietly use "data/cattaneo2.dta", clear

quietly teffects ipw (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), coeflegend
PostCurrentEstimates, handle(`handle') estimator(ipw) statistic(ate) ///
    case(manual_ipw_3)

/* teffects ipwra */
quietly use "data/cattaneo2.dta", clear

quietly teffects ipwra (bweight prenatal1 mmarried mage fbaby) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, probit), ///
    pomeans aequations
PostCurrentEstimates, handle(`handle') estimator(ipwra) statistic(pomeans) ///
    case(manual_ipwra_2)

quietly teffects ipwra (bweight prenatal1 mmarried fbaby c.mage) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu, ///
    hetprobit(c.mage##c.mage)), aequations
PostCurrentEstimates, handle(`handle') estimator(ipwra) statistic(ate) ///
    case(manual_ipwra_3)

/* teffects multivalued */
quietly use "data/bdsianesi5.dta", clear

quietly teffects ra (wage london eastern paed math7, poisson) (ed)
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_multivalued_1_ra)

quietly teffects aipw (wage london eastern paed math7, poisson) ///
    (ed math7 read7 maed paed)
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_multivalued_1_aipw)

quietly teffects, coeflegend
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_multivalued_2_replay)

quietly teffects aipw (wage london eastern paed math7, poisson) ///
    (ed math7 read7 maed paed), control(A) coeflegend
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(ate) ///
    case(manual_multivalued_2_control_a)

quietly teffects ipwra (wage london eastern paed math7, poisson) ///
    (ed math7 read7 maed paed), atet control(A) tlevel(H)
PostCurrentEstimates, handle(`handle') estimator(ipwra) statistic(att) ///
    case(manual_multivalued_3_att)

quietly teffects aipw (wage london eastern paed math7, poisson) ///
    (ed math7 read7 maed paed), pom
PostCurrentEstimates, handle(`handle') estimator(aipw) statistic(pomeans) ///
    case(manual_multivalued_4_pom)

/* teffects nnmatch */
quietly use "data/cattaneo3.dta", clear

quietly teffects nnmatch (bweight mage) (mbsmoke), ///
    ematch(prenatal1 mmarried fbaby) metric(euclidean)
PostCurrentEstimates, handle(`handle') estimator(nnmatch) statistic(ate) ///
    case(manual_nnmatch_2)

quietly teffects nnmatch (bweight mage fage) (mbsmoke), ///
    ematch(prenatal1 mmarried fbaby) biasadj(mage fage)
PostCurrentEstimates, handle(`handle') estimator(nnmatch) statistic(ate) ///
    case(manual_nnmatch_3)

quietly teffects nnmatch (bweight) (mbsmoke), ///
    ematch(i.mmarried i.magecat)
PostCurrentEstimates, handle(`handle') estimator(nnmatch) statistic(ate) ///
    case(manual_nnmatch_4)

quietly teffects ra (bweight i.mmarried##i.magecat) (mbsmoke)
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_nnmatch_4_ra_compare)

/* teffects psmatch */
quietly use "data/cattaneo2.dta", clear

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), caliper(0.1)
PostCurrentEstimates, handle(`handle') estimator(psmatch) statistic(ate) ///
    case(manual_psmatch_2_caliper)

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), ///
    atet vce(iid) caliper(0.03)
PostCurrentEstimates, handle(`handle') estimator(psmatch) statistic(att) ///
    case(manual_psmatch_2_att)

quietly teffects psmatch (bweight) ///
    (mbsmoke mmarried c.mage##c.mage fbaby medu), nneighbor(4)
PostCurrentEstimates, handle(`handle') estimator(psmatch) statistic(ate) ///
    case(manual_psmatch_3)

/* teffects ra */
quietly use "data/cattaneo2.dta", clear

quietly teffects ra ///
    (bweight prenatal1 mmarried mage fbaby) (mbsmoke), ///
    pomeans aequations
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(pomeans) ///
    case(manual_ra_3)

quietly teffects ra ///
    (bweight prenatal1 mmarried mage fbaby) (mbsmoke), coeflegend
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_ra_4)

quietly teffects ra ///
    (bweight prenatal1 mmarried mage fbaby, poisson) (mbsmoke)
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_ra_5)

quietly use "data/pollution.dta", clear
quietly teffects ra ///
    (pollution rainfall i.traffic industrial i.train, fprobit) (guzzler)
PostCurrentEstimates, handle(`handle') estimator(ra) statistic(ate) ///
    case(manual_ra_6)

postclose `handle'
quietly use `results', clear
sort estimator statistic case term
export delimited using "teffects_reference.csv", replace
