#!/usr/bin/env Rscript

BAYES = FALSE
MAXIMAL = FALSE

library(lme4)
if (BAYES) {
    library(rstanarm)
    library(bridgesampling)
}

pr = function(x, buffer=NULL) {
    if (is.null(buffer)) {
        buffer = stderr()
    }
    cat(paste0(x, '\n'), file=buffer, append=TRUE)
}

fit_stan = function(form, df, name, mm=TRUE, iter=2000) {
    cat(paste0('Regressing model ', name, '...\n'))
    #m = lmer(f, REML=F, data=df_)
    diagnostic_file = paste0('models_001/', test, hemi, suffix, '.blme.diag.csv')
    if (mm) {
        m = stan_glmer(form, data=df, iter=iter, diagnostic_file=file.path(tempdir(), "df.csv"))
    } else {
        m = stan_glm(form, data=df, iter=iter, diagnostic_file=file.path(tempdir(), "df.csv"))
    }
    sink(paste0('models_001/', name, '.blme.summary.txt'))
    pr("Formula:", stdout())
    pr(form, stdout())
    pr("Log marginal likelihood:", stdout())
    pr(logml(m), stdout())
    print(summary(m))
    sink()

    readline("")

    return(m)
}

fit_lme4 = function(form, df, name, mm=TRUE) {
    cat(paste0('Regressing model ', name, '...\n'))
    if (mm) {
        m = lmer(f, REML=F, data=df)
    } else {
        m = lm(f, data=df)
    }
    sink(paste0('models_001/', name, '.lme.summary.txt'))
    pr("Formula:", stdout())
    pr(form, stdout())
    print(summary(m))
    sink()

    return(m)
}

fit = function(...) if (BAYES) fit_stan(...) else fit_lme4(...)

compare_stan = function(m1, m0, name) {
    s1 = bridge_sampler(m1)
    s0 = bridge_sampler(m0)

    b = bayes_factor(s1, s0)
    print(logml(s1))
    print(logml(s0))
    print(b)
    b = b[[1]]
    print(b)
    if (b < 0.3333) {
        result = 'Evidence favors the null hypothesis'
    } else if (b > 3) {
        result = 'Evidence favors the alternative hypothesis'
    } else {
        result = 'Evidence inconclusive'
    }

    sum_path = paste0('tests_001/', name, '.bf.summary.txt')

    sink(sum_path)
    pr('==================\nBayes factor analysis\n', stdout())
    pr(paste0('Test name:       \n', name), stdout())
    print(paste0('Bayes factor: ', b))
    print(paste0('Result: ', result))
    pr('------------------\nFull model\n', stdout())
    print(summary(m1))
    pr('\n\n', stdout())
    pr('------------------\nAblated model\n', stdout())
    print(summary(m0))
    sink()
}

compare_lme4 = function(m1, m0, name) {
    lrt = anova(m1, m0)

    sum_path = paste0('tests_001/', name, '.lrt.summary.txt')

    sink(sum_path)
    pr('==================\nLikelihood ratio test\n', stdout())
    pr(paste0('Test name:       \n', name), stdout())
    print(lrt)
    pr('------------------\nFull model\n', stdout())
    print(summary(m1))
    pr('\n\n', stdout())
    pr('------------------\nAblated model\n', stdout())
    print(summary(m0))
    sink()
}

compare = function(...) if (BAYES) compare_stan(...) else compare_lme4(...)

save_fn = function(m, name) {
    if (BAYES) {
        save(m, file=paste0('models_001/', name, '.blme.Rdata'))
    } else {
        save(m, file=paste0('models_001/', name, '.lme.Rdata'))
    }
}

paths = list(

    # Lang

    SENTvNONW_lang=c('../../data/fMRI_ToM/001/langfROIs_langlocSN_n149_001_data.csv'),
    BELIvPHOT_lang=c('../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv'),
    MENTvPHYS_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'),
    TOMVvTOMN_lang=c(
        '../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv',
        '../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'
    ),
    SENTvBELI_lang=c(
        '../../data/fMRI_ToM/001/langfROIs_langlocSN_n149_001_data.csv',
        '../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv'
    ),
    SENTvPHOT_lang=c(
        '../../data/fMRI_ToM/001/langfROIs_langlocSN_n149_001_data.csv',
        '../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv'
    ),
    BELIvNONW_lang=c(
        '../../data/fMRI_ToM/001/langfROIs_langlocSN_n149_001_data.csv',
        '../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv'
    ),
    PHOTvNONW_lang=c(
        '../../data/fMRI_ToM/001/langfROIs_langlocSN_n149_001_data.csv',
        '../../data/fMRI_ToM/001/langfROIs_ToMshort_n149_001_data.csv'
    ),
    MENTvSOCL_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'),
    MENTvPAIN_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'),
    PHYSvPAIN_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'),
    PHYSvSOCL_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv'),
    SOCLvPAIN_lang=c('../../data/fMRI_ToM/001/langfROIs_Cloudy_001_data.csv')

)

tests = names(paths)

options(mc.cores = parallel::detectCores())

if (!dir.exists('models_001')) {
    dir.create('models_001', recursive = TRUE)
}
if (!dir.exists('tests_001')) {
    dir.create('tests_001', recursive = TRUE)
}

for (test in tests) {
    df = data.frame()
    for (path in paths[[test]]) {
        df_ = read.table(path, header=TRUE, sep=',')
        df = rbind(df, df_)
    }
    df$Participant = substr(df$Subject, 1, 3)
    df$ismental = df$Effect %in% c('bel', 'ment')
    df$isverbal = df$Effect %in% c('bel', 'pho')
    conds = c()
    verbal_v_nonverbal = FALSE
    if (grepl('SENT', test, fixed=TRUE)) conds = append(conds, 'S')
    if (grepl('NONW', test, fixed=TRUE)) conds = append(conds, 'N')
    if (grepl('BELI', test, fixed=TRUE)) conds = append(conds, 'bel')
    if (grepl('PHOT', test, fixed=TRUE)) conds = append(conds, 'pho')
    if (grepl('MENT', test, fixed=TRUE)) conds = append(conds, 'ment')
    if (grepl('PHYS', test, fixed=TRUE)) conds = append(conds, 'phys')
    if (grepl('SOCL', test, fixed=TRUE)) conds = append(conds, 'reln')
    if (grepl('PAIN', test, fixed=TRUE)) conds = append(conds, 'pain')
    if (grepl('TOMV', test, fixed=TRUE)) {
        verbal_v_nonverbal = TRUE
        conds = append(conds, 'bel')
        conds = append(conds, 'pho')
    }
    if (grepl('TOMN', test, fixed=TRUE)) {
        verbal_v_nonverbal = TRUE
        conds = append(conds, 'ment')
        conds = append(conds, 'phys')
    }
    df = df[df$Effect %in% conds,]
    if (verbal_v_nonverbal) {
        subjs = unique(df[df$Effect == 'ment',]$Participant)
        df = df[df$Participant %in% subjs,]
    }
    if (substr(test, 11, nchar(test)) == 'lang') {
        hemis = c('RH', 'AngG')
    } else {
        hemis = c('')
    }
    for (hemi in hemis) {
        m1 = NULL
        m0 = NULL
        if (hemi == 'LH') {
            df_ = df[df$ROI < 6,]
        } else if (hemi == 'RH') {
            df_ = df[(df$ROI > 6) & (df$ROI < 12),]
        } else if (hemi == 'AngG') {
            df_ = df[(df$ROI == 6) | (df$ROI == 12),]
        } else {
            df_ = df
        }
        for (suffix in c('', '_h0')) {
            name = paste0(test, hemi, suffix)
            if (suffix == '_h0') {
                if (MAXIMAL) {
                    if (verbal_v_nonverbal) {
                        f = 'EffectSize ~ ismental + isverbal + (ismental * isverbal | Participant) + (ismental * isverbal | ROI)'
                    } else {
                        f = 'EffectSize ~ (Effect | Participant) + (Effect | ROI)'
                    }
                } else {
                    if (verbal_v_nonverbal) {
                        f = 'EffectSize ~ ismental + isverbal + (1 | Participant) + (1 | ROI)'
                    } else {
                        f = 'EffectSize ~ (1 | Participant) + (1 | ROI)'
                    }
                }
            } else {
                if (MAXIMAL) {
                    if (verbal_v_nonverbal) {
                        f = 'EffectSize ~ ismental * isverbal + (ismental * isverbal | Participant) + (ismental * isverbal | ROI)'
                    } else {
                        f = 'EffectSize ~ Effect + (Effect | Participant) + (Effect | ROI)'
                    }
                } else {
                    if (verbal_v_nonverbal) {
                        f = 'EffectSize ~ ismental * isverbal + (1 | Participant) + (1 | ROI)'
                    } else {
                        f = 'EffectSize ~ Effect + (1 | Participant) + (1 | ROI)'
                    }
                }
            }

            m = fit(f, df_, name, mm=TRUE)
            if (suffix == '_h0') {
                m0 = m
            } else {
                m1 = m
            }
            save_fn(m, name)
        }
        compare(m1, m0, paste0(test, hemi))

        for (roi in unique(df_$ROI)) {
            df__ = df_[df_$ROI == roi,]
            m1 = NULL
            m0 = NULL
            for (suffix in c('', '_h0')) {
                name = paste0(test, hemi, '_roi', roi, suffix)
                if (suffix == '_h0') {
                    if (MAXIMAL) {
                        if (verbal_v_nonverbal) {
                            f = 'EffectSize ~ ismental + isverbal + (ismental * isverbal | Participant)'
                        } else {
                            f = 'EffectSize ~ (1 | Participant)'
                        }
                    } else {
                        if (verbal_v_nonverbal) {
                            f = 'EffectSize ~ ismental + isverbal'
                        } else {
                            f = 'EffectSize ~ 1'
                        }
                    }
                } else {
                    if (MAXIMAL) {
                        if (verbal_v_nonverbal) {
                            f = 'EffectSize ~ ismental * isverbal + (ismental * isverbal | Participant)'
                        } else {
                            f = 'EffectSize ~ Effect + (1 | Participant)'
                        }
                    } else {
                        if (verbal_v_nonverbal) {
                            f = 'EffectSize ~ ismental * isverbal'
                        } else {
                            f = 'EffectSize ~ Effect'
                        }
                    }
                }

                m = fit(f, df__, name, mm=MAXIMAL)
                if (suffix == '_h0') {
                    m0 = m
                } else {
                    m1 = m
                }
                save_fn(m, name)
            }
            compare(m1, m0, paste0(test, hemi, '_roi', roi))
        }
    }
}
