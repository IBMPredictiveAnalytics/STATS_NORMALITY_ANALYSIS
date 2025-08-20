#/***********************************************************************
# *
# * (C) Copyright Jon K Peck, 2024
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# version 1.0.0

# history
# 7-aug-2024 original version
# 16-sep-2024 finishing touches including variable case and TO support
# 25-sep-2024 adjustments to Mardia output
# 14-jul-2025 major rewrite to deal with problems in mvn


# helpers
gtxt <- function(...) {
    return(gettext(...,domain="STATS_NORMALITY_ANALYSIS"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_NORMALITY_ANALYSIS"))
}

loadmsg = "The R %s package is required but could not be loaded."
tryCatch(suppressWarnings(library(MVN, warn.conflicts=FALSE)), error=function(e){
    stop(gtxtf(loadmsg,"MVN"), call.=FALSE)
}
)
tryCatch(suppressWarnings(library(energy, warn.conflicts=FALSE)), error=function(e){
    stop(gtxtf(loadmsg,"energy"), call.=FALSE)
}
)
tryCatch(suppressWarnings(library(MASS, warn.conflicts=FALSE)), error=function(e){
    stop(gtxtf(loadmsg,"MASS"), call.=FALSE)
}
)
tryCatch(suppressWarnings(library(nortest, warn.conflicts=FALSE)), error=function(e){
    stop(gtxtf(loadmsg,"nortest"), call.=FALSE)
}
)
tryCatch(suppressWarnings(library(moments)), error=function(e){
    stop(gtxtf(loadmsg,"moments"), call.=FALSE)
}
)
tryCatch(suppressWarnings(library(boot)), error=function(e){
    stop(gtxtf(loadmsg,"boot"), call.=FALSE)
}
)


mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment
    
    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.
        
        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 
        
        if (is.null(msg) || dostop) {
            spssdata.CloseDataConnection()
            lcl$display(inproc)  # display messages and end procedure state

            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any
        
        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
                procok = TRUE
            }
        } else {
            procok = inproc
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    procok = TRUE
                },
                error = function(e) {
                    prockok = FALSE
                }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings", isSplit=FALSE) # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row,
                                               gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory,
                                                spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

casecorrect = function(vlist, warns) {
    # correct the case of variable names
    # vlist is a list of names, possibly including TO and ALL
    # unrecognized names are returned as is as the GetDataFromSPSS api will handle them

    tryCatch(
        {
        dictnames = spssdictionary.GetDictionaryFromSPSS()["varName",]
        }, error = function(e) {warns$warn(gtxt("The active dataset has no variables"),
                 dostop=TRUE)}
    )
    names(dictnames) = tolower(dictnames)
    dictnames['all'] = "all"
    dictnames['to'] = "to"
    correctednames = list()
    for (item in vlist) {
        lcitem = tolower(item)
        itemc = dictnames[[lcitem]]
        if (is.null(itemc)) {
            warns$warn(gtxtf("Invalid variable name: %s", item), dostop=TRUE)
        }
        correctednames = append(correctednames, itemc)
    }

    return(correctednames)
}

reorderbyvar = function(ut, variables) {
    # return the data frame in order of the variables
    # ut is a data frame with the first column containing variable names
    # variables is a list of variable names in the desired order
    
    # mvn returns its output with a leading blank 
    # on the variable names column values, so remove them
    # make the first column a factor with levels in order of variables
    uu = ordered(trimws(ut[[1]]), levels=variables)
    udfbind = cbind(uu, ut[-1])
    # sort the data frame by varnames
    udfbind = udfbind[order(udfbind[[1]]),]
    return(udfbind)
}


procname=gtxt("Normality Analysis")
warningsprocname = gtxt("Normality Analysis")
omsid="STATSNORMALITY"
warns = Warn(procname=warningsprocname,omsid=omsid)

univartests = c(sw="SW", cvm="CVM", lillie="Lillie", sf="SF", ad="AD")

# main worker
domvn<-function(idvar=NULL, variables, mvntests=NULL, univariatetests=NULL, 
    bootstrapreps=1000, uniplots=TRUE, scatterplots=FALSE, nscatterplotvars=10, 
    scatterplotsize=100, noutliers=0, outlierdetection="quan", 
    scaledata=FALSE, desc=FALSE) {

    domain<-"STATS_NORMALITY_ANALYSIS"
    setuplocalization(domain)
    # DEBUG
    ###sink(file="c:/temp/normout.log", type="output")
    ###f = file("c:/temp/normsgs.log", open="w")
    ###sink(file=f, type="message")

    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(gtxt("Case weights are not supported by this procedure and will be ignored"), dostop=FALSE)
    }
    ut = list()
    for (item in univariatetests) {ut[[item]] = univartests[[item]]}   # case correct test keywords
    univariatetests = ut
    if (length(mvntests) > 0) {
        for (item in 1:length(mvntests)) {
            if (mvntests[[item]] == "dh") {
                mvntests[item] = "doornik_hansen"
            }
        }
    }

    # correct variable name case, including the id variable, if any
    variables = casecorrect(c(variables, idvar), warns)  # get data api requires case match
    if (!is.null(idvar)) {
        idvar = variables[[length(variables)]]
        variables = variables[-length(variables)]
    }
    
    splitvars = spssdata.GetSplitVariableNames()
    nsplitvars = length(splitvars)
    if (length(intersect(tolower(variables), tolower(splitvars))) > 0) {
        warns$warn(gtxt("Split variables cannot be included in the list of variables to analyze"), dostop=TRUE)
    }

    # nvars can be an underestimate if TO or ALL is used, but TO would have to involve
    # at least two variables, and ALL is very unlikely to be used.
    # TO or ALL in splitvars might mess things up but that would be very unlikely usage.
    
    nvars = length(variables)
    if (nvars < 2) {
        warns$warn(gtxt("At least two variables must be specified"), dostop=TRUE)
    }
    # splitlist will hold the split number for each plot file generated
    # in order to assign splits to plots later
    splitlist = list()
    if (length(splitvars) > 0) {
        needsplittbl = TRUE
        splittbl = data.frame(matrix(ncol = nsplitvars, nrow=0))
        colnames(splittbl) = list(splitvars)
    } else {
        needsplittbl = FALSE
    }
    splitnumber = 0
    caption = gtxtf("Tests computed by R MVN package, version %s", packageVersion("mvn"))  

    # get all the variables, including split vars but drop split vars from dta after extracting the split values.
    # factor mode is "labels" in order to pick up labelled split values, if there are splits.
    # Procedure will only use complete cases.
    # SPSS date variables are not converted via rDate=POSIXct, because mvn code does not handle
    # R date variables correctly
    
    #  plotlistfile will hold the file names for all the generated plots
    
    if (uniplots || scatterplots) {
        plotlistfile = tempfile("plotlist", tmpdir=tempdir(), fileext=".txt")
        pf = file(plotlistfile, open="w")   # accumulate list of generated plot files
        ###print(sprintf("plot list file: %s, %s", plotlistfile, pf))
    }
    spsspkg.StartProcedure(gtxt("Normality Analysis"),"STATS NORMALITY ANALYSIS")
    
    # main analysis loop
    
    while (!spssdata.IsLastSplit()) {
        if (is.null(idvar)) {
            dta = spssdata.GetSplitDataFromSPSS(paste(c(variables, splitvars), collapse=" "), 
                missingValueToNA=TRUE, factorMode="labels")
        } else {
            # this will fail if id values are not unique
            tryCatch(
            {
            dta = spssdata.GetSplitDataFromSPSS(paste(c(variables, splitvars), collapse=" "), 
                missingValueToNA=TRUE, factorMode="labels", row.label=idvar)
            },
            error = function(e) {
                warns$warn(e, dostop=TRUE)
            } 
            )
        }
        splitnumber = splitnumber + 1
        # save split values
        if (needsplittbl) {
            dd = data.frame(dta[1, (nvars + 1):ncol(dta)])
            names(dd) = names(splittbl)
            splittbl = rbind(splittbl, dd)
            dta = dta[1:nvars]   # remove split vars
            splitprefix = gtxtf("(Split %d)", splitnumber)   # for labelling charts
        } else {
            splitprefix = ""
        }
        if (any(sapply(dta, is.factor))) {
            warns$warn(gtxtf("Categorical variables cannot be used in this procedure"), dostop=TRUE)
        }
        if (scaledata) {
            dta = scale(dta)
            scaledata2 = TRUE
        } else {
            scaledata2 = FALSE
        }

        ncases = nrow(dta)
        if (("royston" %in% mvntests || "SW" %in% univariatetests) && (ncases > 5000 || ncases < 3)) {  # missing value cases will be discarded later, which could invalidate this test
            warns$warn(gtxt("The Royston and Shapiro-Wilk tests cannot be used with more than 5000 or fewer than 3 cases"), dostop=TRUE)
        }

        if (desc) {
            dodesc(dta, scaledata2)
        }
        ###save(dta, desc, nvars, variables, mvntests, univariatetests, uniplots, scatterplots, file="c:/temp/dodesc.rdata")
        # mvn insists on a univariate and a multivariate test or raises an error
        # Univariate normality tests
        if (length(univariatetests) > 0) {
            douniv(dta, univariatetests, scaledata2, caption, variables)
        }

        if (length(mvntests) > 0) {
            domv(dta, mvntests, bootstrapreps, scaledata2, caption)
        }
       newsplits = dographics(dta, needsplittbl, splittbl, splitvars, splitnumber, nvars,  
            scaledata, splitprefix, uniplots, 
            scatterplots, nscatterplotvars, scatterplotsize, noutliers, idvar, outlierdetection, pf)
       splitlist = c(splitlist, newsplits)
    }
    
    spssdata.CloseDataConnection()
    close(pf)
    drawgraphs(plotlistfile, needsplittbl, splittbl, splitvars, splitlist)
    
    # print messages and clean up
    warns$display(inproc=FALSE)
    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})
    # DEBUG
    ###sink(file=NULL, type="output")
    ###sink(file=NULL, type="message")
}


dodesc = function(dta, scaleddata) {
    # return descriptive statistics for dta
    
    # scaleddata indicates whether dta has  been standardized
    ###tryCatch({desctable = suppressWarnings(mvn(dta, scale=FALSE, desc=TRUE))$Descriptives},
    tryCatch({desctable = suppressWarnings(mvn(dta, scale=FALSE, descriptives=TRUE))$descriptives},
             error=function(e){warns$warn(
                 gtxt("Cannot compute descriptives due to data conditions.\n  Perhaps too few complete cases, too little variance, or too highly correlated variables"), dostop=TRUE)
             }
    )
             
    if (scaleddata) {
        caption = gtxt("Variables are standardized")
    }
    else {
        caption = gtxt("Variables are not standardized")
    }
    row.names(desctable) = colnames(dta)
    desctable = desctable[-1]
    colnames(desctable) = c(gtxt("n"), gtxt("Mean"), gtxt("Std.Dev."), gtxt("Median"), gtxt("Min"),
                            gtxt("Max"), gtxt("25th"), gtxt("75th"), gtxt("Skewness"), gtxt("Kurtosis"))

    spsspivottable.Display(desctable, 
                           title=gtxt("Descriptive Statistics"),
                           templateName="UNIVARIATESTATS",
                           isSplit=TRUE,
                           rowdim = gtxt("Variables"), 
                           hiderowdimlabel=FALSE, 
                           hidecoldimtitle=TRUE,
                           format=formatSpec.GeneralStat,
                           caption=gtxt(caption)
    )
}


douniv = function(dta, univariatetests, scaledata, caption, variables) {
    ###ut = data.frame('Variable'=NULL, 'Test'=NULL, 'Statistic'=NULL, 'p Value'=NULL)

    ut =  data.frame(Variable=character(), Test=character(), Statistic=double(), 'p value'=double())
    for (item in univariatetests) {
        tryCatch({
            ###res = suppressWarnings(mvn(dta, univariateTest=item, scale=FALSE, desc=FALSE, mvnTest="mardia"))   #
            res = suppressWarnings(mvn(dta, univariate_test=item, scale=FALSE, 
                descriptives=FALSE, mvn_test="hz"))   
        #ignore the mandatory mvnTest
            ut = rbind(ut, res$univariate_normality[c(2, 1, 3, 4)])
        }, 
        error = function(e) {warns$warn(gtxtf("Test %s cannot be calculated", item), dostop=FALSE)}
        )
    }
    # for translation...
    if (nrow(ut) > 0) {
        ut = reorderbyvar(ut, variables)
        colnames(ut) = c(gtxt("Variable"), gtxt("Test"), gtxt("Statistic"), gtxt("P Value"))

        spsspivottable.Display(ut,
                               title=gtxt("Univariate Tests"),
                               templateName="UNIVARIATENORMALITY",
                               hiderowdimtitle=TRUE, 
                               hidecoldimtitle=TRUE,
                               ###rowlabels=as.character(1:nrow(ut)),
                               format=formatSpec.GeneralStat,
                               caption=caption)
    }
}

domv = function(dta, mvntests, bootstrapreps, scaledata2, caption) {
    # multivariate tests
    tryCatch({
        res = mvtestresults(dta, mvntests, bootstrapreps, scaledata2)
    },
    error=function(e) {
        warns$warn(gtxt("Multivariate tests cannot be calculated.\n Perhaps there are too few complete cases, too little variance, or too highly correlated variables"), dostop=TRUE)}
    )
    mt = res[[1]]
    
    # for translation
    colnames(mt) = c(gtxt("Test", "Statistic", "P Value"))
    if (!is.null(res[[2]])) {
        # what to do for split files d.fs?
        caption = sprintf(gtxtf("Doornik-Hansen degrees of freedom: %s\n%s", res[2], caption))
    }

    spsspivottable.Display(mt,
       title = gtxt("Multivariate Normality Tests"),
       templateName="MULTIVARIATENORMALITY",
       rowdim = gtxt("Tests"), 
       hiderowdimlabel=TRUE, 
       hidecoldimtitle=TRUE,
       format=formatSpec.GeneralStat,
       caption = caption)
}

doboxplot = function(dta, nvars) {
    pfilebox = tempfile("box", tmpdir=tempdir(), fileext = ".png")
    if (nvars <= 10) {
        width = 600
    } else {
        width = 1000
    }
    ff = png(pfilebox, width=width, height=400, units="px")
    boxplot(dta)
    dev.off()
    return(pfilebox)
}

dographics = function(dta, needsplittbl, splittbl, splitvars, splitnumber, nvars, scaledata, splitprefix,  uniplots, scatters, nscatterplotvars, scatterplotsize, noutliers, idvar, outlierdetection, pf) {
    # produce univariate and multivariate png files as requested and return
    # file listing the files but do not draw graphs
    # Also produce the outlier table if requested
    # return the splitnumbers for any graphs generated
    
    newsplitlist = c()
    if (uniplots) {  # boxplot and histogram-qq pairs
        # boxplot
        ff = doboxplot(dta, nvars)
        newsplitlist = c(newsplitlist, splitnumber)
        writeLines(text=ff, con=pf)
        # histogram and qq plots
        pfile = comboplot(dta)
        newsplitlist = c(newsplitlist, splitnumber)
        writeLines(text=pfile, con=pf)
        ###print("done with uniplots")
    }

    # do scatterplots as SPLOM
    if (nvars >= 2 && scatters) {
        ###print("doing scatterplots")
        dtalimit = min(nscatterplotvars, ncol(dta))
        ###save(dta, dtalimit, file="c:/temp/dtaetc.rdata")
        pfile = scatplot(data.frame(dta[, 1:dtalimit]), scatterplotsize)
        if (!is.null(pfile)) {
            newsplitlist = c(newsplitlist, splitnumber)
        }
        writeLines(text=pfile, con=pf)
    }
        
    if (noutliers > 0) {
            tryCatch(
                {res = 
                    mvn(dta, subset=NULL, scale=FALSE, 
                    multivariate_outlier_method = outlierdetection)
                },
                error = function(e) {print(e)
                    warns$warn(gtxt("Outlier analysis cannot be completed.\n Perhaps there are too few complete cases, too little variance, or too highly correlated variables"), dostop=TRUE)
                }
            )
        ntoshow = min(noutliers, nrow(res$multivariate_outliers))
        if (ntoshow == 0) {
            warns$warn(gtxt("There are no outliers"))
        } else {
            oo = data.frame(res$multivariate_outliers)
            oo = head(oo, noutliers)
            colnames(oo) = c(idvar, gtxt("Mahalanobis Distance"))
            ###save(oo, res, file="c:/temp/oo.rdata")
            spsspivottable.Display(
                oo[2],
                title = gtxt("Top Outliers"),
                templateName = "NORMALITYOUTLIERS",
                rowdim = idvar,
                rowlabels = oo[1],
                hiderowdimtitle=FALSE,
                hidecoldimtitle=TRUE,
                caption = gtxtf("Detection Method: %s.  Outlier Limit: %s", outlierdetection, noutliers)
            )
        }
    }
    ###print("done dographics")
    return(newsplitlist)
}


drawgraphs = function(pf, needsplittbl, splittbl, splitvars, splitlist) {
    # display split table if needed and draw graphs
    if (needsplittbl && !is.null(pf)) {
        names(splittbl) = splitvars
        spsspivottable.Display(
            splittbl,
            title = gtxtf("Table of Splits"),
            templateName = "SPLITTABLE",
            isSplit = FALSE,
            rowdim = gtxt("Split"),
            hiderowdimtitle=FALSE,
            hidecoldimtitle=TRUE,
            rowlabels=as.character(seq(1:nrow(splittbl))),
            caption = gtxtf("Use this table to identify splits in plots")
        )
    }
    spsspkg.EndProcedure()
    if (is.null(pf)) {
        return()
    }
    # insert the graphs as listed in pf
    # roundabout through INSERT, because I18N characters are not handled properly
    # in labels
    
    outlinelabel = gtxt("Normality Analysis")
    labelparm = sprintf("%s", paste(splitlist, collapse=" "))
    if (needsplittbl) {
        lbparm = sprintf("LABELPARM= %s", labelparm)
        #lbparm = paste(sprintf('"(Split %s)"', splitlist, sep= " "))
    } else {
        lbparm = c()
    }

    pf2 = pf
    cmd = c("* Encoding: UTF-8.",
            sprintf("STATS INSERT CHART "),
            lbparm,
            sprintf("CHARTLIST='%s'", pf),
            "HEADER='Normality PlotS'", 
            sprintf("OUTLINELABEL='%s '", outlinelabel)
    )
    ###print(sprintf("in drawgraphs.  cmd: %s", cmd))
    syntemp = tempfile("synplt", tmpdir=tempdir(), fileext=".sps")
    writeLines(text=cmd, con=syntemp, useBytes=TRUE)
    
    # STATS INSERT CHART's xml file is not read on first invocation after installation

    tryCatch(
        {
            spsspkg.Submit(sprintf("INSERT FILE='%s' ENCODING='UTF8'", syntemp))
            unlink(syntemp)
        },
        error = function(e) {
            print(e)
            print("Please restart SPSS Statistics to complete installation of this command")
        }
    )
}

comboplot = function(dta, width=600, height=300) {
    # create png file of pairs of histograms and qqnormal plots

    pfile = tempfile("histqq", tmpdir=tempdir(), 
        fileext=".png")
    dta = data.frame(dta)
    rowsofplots = ncol(dta)
    height = max(height, 200 * rowsofplots)
    width = max(width, 400)
    cols = 2
    varnames = names(dta)

    tryCatch(
        {
            png(pfile, width=width, height=height, units="px")
            par(mfrow=c(rowsofplots, cols))
            for (v in 1:rowsofplots) {
                main = varnames[[v]]
                hist(dta[[v]], main=main, xlab=varnames[v])
                qqnorm(dta[[v]], main=main, xlab=varnames[v])
                qqline(dta[[v]])
            }
            dev.off()
        }, error=function(e) {print(e)}
    )
    
    return(pfile)
}

scatplot = function(dta, scatterplotsize=100, title=NULL) {
    # create png file of scatterplots for dta variables
    
    pfile = tempfile("splom", tmpdir=tempdir(), fileext=".png")
    nvars = ncol(dta)
    size = scatterplotsize
    tryCatch(
        {
        png(pfile, width = nvars * size, height = nvars * size)
        pairs(dta, main=title)
        dev.off()
        }, error = function(e) {print(e)
            ###save(e, size, dta, file="c:/temp/scatter.rdata")
            warns$warn("Could not draw scatterplot matrix.  Perhaps too many variables",
                dostop=TRUE)
        return(NULL)
        }
    )
    return(pfile)
}


mvtestresults <- function(dta, mvntests, bootstrapreps, scaledata) {
    # return a list with a data frame containing the multivariate test results and the d.f. for dh or NULL
    
    # dta is the data frame to analyze
    # mvn is a list of the mv test parameters
    # bootstrapreps is the number of replications for the energy e test
    
    # Each result is a list of test name, statistic, and p value
    # test result structure returned by mvn vary by test type :-(
    
    mt = data.frame(Test=numeric(), Statistic=numeric(), "p value"= numeric())
    row = 1
    dhdegf = NULL # for Doornik-Hansen
    
    for (item in mvntests) {
        mvt = suppressWarnings(mvn(dta, mvn_test=item, scale=FALSE,   
            descriptives=FALSE, B=bootstrapreps, univariate_test="AD", 
            subset=NULL)) # ignore univariate

        if (item == "doornik-hansen") {
            dhdegf = mvt$multivariate_normality[1, 3]
            mvt$multivariate_normality[1,3] = mvt$multivariate_normality[1,4] # squeeze out df value
        }
        if (item != "mardia") {
            mt[row,] = mvt$multivariate_normality[1, c(1,2,3)]
            mt[row, 2] = round4(mt[[row, 2]])
            mt[row, 3] = round4(mt[[row, 3]])
            if (mt[row, 3] == 0) {
                mt[row, 3] = "<.001"
            }
            row = row + 1
        } else {  # mardia has two values
            mtg = getmardia(mvt$multivariate_normality)
            mt = rbind(mt, mtg)
            row = row + 2
        }

    }
    return(list(mt, dhdegf))
}


round4 = function(x) {
    if (is.numeric(x)) {
        return (round(x, 4))
    }
    return(x)
}

getmardia = function(mvnormality) {  # mardia has test and p value as factors
    
    mt = data.frame(Test=numeric(), Statistic=numeric(), "p value"=numeric())
    mt[1, 1] = mvnormality$Test[1]
    mt[1, 2] = round4(as.numeric(mvnormality$Statistic)[[1]])
    respv = round4(mvnormality$"p.value")
    pvs = respv
    # mvn only returns one sig value if skewness and kurtosis sig values are the same
    # this many happen at the extremes where values are both 0.
    if (length(pvs) == 1) {
        pvs[2] = pvs[1]
    }
    mt[1, 3] = round4(pvs[1])  #skewness sig
    
    if (mt[1, 3] == 0) {
        mt[1, 3] = "<.001"
    }
    
    mt[2, 1] = mvnormality$Test[2]  # Kurtosis
    mt[2, 2] = round4(as.numeric(mvnormality$Statistic[[2]]))
    mt[2, 3] = round4(pvs[[2]])  # kurtosis sig
    if (mt[2, 3] == 0) {
        mt[2, 3] = "<.001"
    }
    return(mt)
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    if (!is.null(fpath)) {
        bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
    }
} 


Run<-function(args){

    cmdname = args[[1]]
    args <- args[[2]]
    ###options(error = function() traceback())

    # note SW, CVM Lillie, SF, AD
    
    # variable keywords are typed as varname instead of existingvarlist in
    # order to allow for case correction of names later, since the data fetching apis are
    # case sensitive
    
    oobj <- spsspkg.Syntax(templ=list(
        spsspkg.Template("VARIABLES", subc="", ktype="varname", var="variables", islist=TRUE),
        
        spsspkg.Template("MVNTESTS", subc="OUTPUT", ktype="str", var="mvntests", 
            vallist=list("mardia", "hz", "royston", "dh", "energy"), islist=TRUE),
        spsspkg.Template("UNIVARIATETESTS", subc="OUTPUT", ktype="str", var="univariatetests", 
            vallist=list("sw", "cvm", "lillie", "sf", "ad"), islist=TRUE),
        spsspkg.Template("BOOTSTRAPREPS", subc="OUTPUT", ktype="int", var="bootstrapreps", islist=FALSE),
        spsspkg.Template("SCALEDATA", subc="OUTPUT", ktype="bool", var="scaledata", islist=FALSE),
        spsspkg.Template("UNIPLOTS", subc="OUTPUT", ktype="bool", var="uniplots",
            islist=FALSE),
        spsspkg.Template("SCATTERPLOTS", subc="OUTPUT", ktype="bool", var="scatterplots",
             islist=FALSE),
        spsspkg.Template("NSCATTERPLOTVARS", subc="OUTPUT", ktype="int", 
            var="nscatterplotvars", islist=FALSE, vallist=list(1,1000)),
        spsspkg.Template("SCATTERPLOTSIZE", subc="OUTPUT", ktype="int", var="scatterplotsize",
            islist=FALSE, vallist=list(50,1000)),
        spsspkg.Template("DESCRIPTIVES", subc="OUTPUT", ktype="bool", var="desc", islist=FALSE),

        spsspkg.Template("IDVAR", subc="OUTLIERS", ktype="varname", var="idvar", islist=FALSE),
        spsspkg.Template("NOUTLIERS", subc="OUTLIERS", ktype="int", var="noutliers", islist=FALSE),
        spsspkg.Template("OUTLIERDETECTION", subc="OUTLIERS", ktype="str", var="outlierdetection",
            vallist=list("quan", "adj"),  islist=FALSE)
        ))

    if ("HELP" %in% attr(args,"names"))
        helper(cmdname)
    else {
        res <- spsspkg.processcmd(oobj, args, "domvn")
    }
}


helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
    if (exists("spsspkg.helper")) {
    assign("helper", spsspkg.helper)
}
