#/***********************************************************************
# *
# * (C) Copyright Jon K Peck, 2024
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# version 1.0.0

# history
# 31-jul-2024 original version

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

# helpers
gtxt <- function(...) {
    return(gettext(...,domain="STATS_NORMALITY_ANALYSIS"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_NORMALITY_ANALYSIS"))
}


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
    
    ### dictnames = spssdictionary.GetDictionaryFromSPSS()["varName", ]

    dictnames = spssdictionary.GetDictionaryFromSPSS()["varName",]
    names(dictnames) = tolower(dictnames)
    correctednames = list()
    for (item in vlist) {
        lcitem = tolower(item)
        if (lcitem %in% dictnames) {
            itemc = dictnames[[lcitem]]
            correctednames = append(correctednames, itemc)
        } else {
            if (!(lcitem %in% list("to", "all"))) {
                warns$warn(gtxtf("Invalid variable name: %s", item), dostop=TRUE)
            }
        }

    }
    return(correctednames)
}


procname=gtxt("Normality Analysis")
warningsprocname = gtxt("Normality Analysis")
omsid="STATSNORMALITY"
warns = Warn(procname=warningsprocname,omsid=omsid)

# main worker
univartests = c(sw="SW", cvm="CVM", lillie="Lillie", sf="SF", ad="AD")


domvn<-function(idvar=NULL, variables, mvntests=NULL, univariatetests=NULL, bootstrapreps=1000,
    uniplots=NULL, multivarplots=NULL, noutliers=0, outlierdetection="quan", 
    scaledata=FALSE, desc=FALSE) {

    domain<-"STATS_NORMALITY_ANALYSIS"
    setuplocalization(domain)

    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(gtxt("Case weights are not supported by this procedure and will be ignored"), dostop=FALSE)
    }
    ut = list()
    for (item in univariatetests) {ut[[item]] = univartests[[item]]}   # case correct test keywords
    univariatetests = ut
    spsspkg.StartProcedure(gtxt("Normality Analylsis"),"STATS NORMALITY ANALYSIS")
    variables = casecorrect(variables, warns)  # get data api requires case match
    splitvars = spssdata.GetSplitVariableNames()
    nsplitvars = length(splitvars)
    if (length(intersect(tolower(variables), tolower(splitvars))) > 0) {
        warns$warn(gtxt("Split variables cannot be included in the list of variables to analyze"), dostop=TRUE)
    }

    nvars = length(variables)
    nallvars = nvars
    if (length(splitvars) > 0) {
        needsplittbl = TRUE
        splittbl = data.frame(matrix(ncol = nsplitvars, nrow=0))
        colnames(splittbl) = list(splitvars)
    } else {
        needsplittbl = FALSE
    }
    splitnumber = 0
    if (noutliers > 0) {
        if (length(variables) < 2) {
            warns$warn(gtxt("At least two variables must be named for outlier analysis"), dostop=TRUE)
        }
    }
    if (desc) {
        if (nvars < 2) {
            warns$warn(gtxt("Descriptives requires at least two variables"), dostop=FALSE)
        }
    }
    if ("contour" %in% multivarplots || "persp" %in% multivarplots) {
        warns$warn(gtxt("Contour and Perspective plots are limited to two variables.  The first two will be plotted"),
            dostop=FALSE)
    }

    # get all the variables, including split vars but drop split vars from dta after extracting the split values.
    # factor mode is "labels" in order to pick up labelled split values, if there are splits.
    # procedure will only use complete cases.
    while (!spssdata.IsLastSplit()) {
        # if (is.null(idvar)) {
        #     dta = spssdata.GetSplitDataFromSPSS(c(idvar, variables, splitvars), missingValueToNA=TRUE, factorMode="labels")
        # } else {
        #     dta = spssdata.GetSplitDataFromSPSS(c(idvar, variables), missingValueToNA=TRUE, factorMode="labels",
        #     row.label=idvar)
        # }
        if (is.null(idvar)) {
            dta = spssdata.GetSplitDataFromSPSS(c(variables, splitvars), missingValueToNA=TRUE, factorMode="labels")
        } else {
            dta = spssdata.GetSplitDataFromSPSS(c(variables, splitvars), missingValueToNA=TRUE, factorMode="labels",
                                                row.label=idvar)
        }

        splitnumber = splitnumber + 1
        # save split values
        if (needsplittbl) {
            dd = data.frame(dta[1, (nallvars + 1):ncol(dta)])
            names(dd) = names(splittbl)
            splittbl = rbind(splittbl, dd)
            dta = dta[1:nallvars]   # remove split vars
            splitprefix = gtxtf("(Split %d)", splitnumber)   # for labelling charts
        } else {
            splitprefix = ""
        }
        
        if (any(sapply(dta, is.factor))) {
              warns$warn(gtxtf("Categorical variables cannot be used in this procedure"), dostop=TRUE)
        }
        ncases = nrow(dta)
        if (("royston" %in% mvntests || "SW" %in% univariatetests) && (ncases > 5000 || ncases < 3)) {  # missing value cases will be discarded later
            warns$warn(gtxt("The Royston and Shapiro-Wilk tests cannot be used with more than 5000 or fewer than 3 cases"), dostop=TRUE)
        }
    
        caption = gtxtf("Computed by R MVN package, version %s", packageVersion("mvn"))
        
        if (desc && nvars >= 2) {
                desctable = suppressWarnings(mvn(dta, scale=scaledata, desc=TRUE)$Descriptives)
                if (scaledata) {
                    caption = gtxt("Variables are standardized")
                }
                else {
                    caption = gtxt("Variables are not standardized")
                }
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
    
        # mvn insists on a univariate and a multivariate test or raises an error
        # Univariate normality tests
        if (length(univariatetests) > 0) {
            ut = data.frame('Variable'=NULL, 'Test'=NULL, 'Statistic'=NULL, 'p Value'=NULL)
            for (item in univariatetests) {
                tryCatch( {
                    res = suppressWarnings(mvn(dta, univariateTest=item, scale=scaledata, desc=FALSE, mvnTest="mardia"))   # ignore the mandatory mvnTest
                    ut = rbind(ut, res$univariateNormality[c(2, 1, 3, 4)])
                }, error = function(e) {warns$warn(gtxtf("Test %s cannot be calculated", item), dostop=FALSE)}
                )
            }
            # for translation...
            colnames(ut) = c(gtxt("Variable"), gtxt("Test"), gtxt("Statistic"), gtxt("P Value"))
            spsspivottable.Display(ut,
                  title=gtxt("Univariate Tests"),
                  templateName="UNIVARIATENORMALITY",
                  isSplit=TRUE,
                  rowdim = gtxt("Variables"), 
                  hiderowdimlabel=TRUE, 
                  hidecoldimtitle=TRUE,
                  format=formatSpec.GeneralStat,
                  caption=caption)
            }
            
            if (length(mvntests) > 0) {
                res = mvtestresults(dta, mvntests, bootstrapreps, scaledata)
                mt = res[[1]]
                # for translation
                colnames(mt) = c(gtxt("Test", "Statistic", "P Value"))
    
                if (!is.null(res[[2]])) {
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

        # outliers
        if (noutliers > 0) {
            spssRGraphics.SetGraphicsLabel(gtxtf("%s Normality Analysis - Outliers", splitprefix))
            ntoshow = min(noutliers, nrow(res$multivariateOutliers))
            res = suppressWarnings(mvn(dta, subset=NULL, scale=scaledata, 
                multivariateOutlierMethod = outlierdetection, showOutliers=TRUE))
            oo = data.frame(res$multivariateOutliers[2])
            oo = head(oo, noutliers)
            colnames(oo) = gtxt("Mahalanobis Distance")
            spsspivottable.Display(
                oo,
                title = gtxtf("Top %s Outliers Sorted by Mahalanobis Distance", noutliers),
                templateName = "NORMALITYOUTLIERS",
                rowdim = idvar,
                hiderowdimtitle=FALSE,
                hidecoldimtitle=TRUE,
                caption = gtxtf("Detection Method: %s", outlierdetection)
            )
        }
        
        # graphics
        # if (("contour" %in% multivarplots || "scatter" %in% multivarplots) && nvars > 2) {
        #         warns$warn(gtxt("Contour and scatter plots cannot be used with more than two variables"), dostop=TRUE)
        # }
        if (needsplittbl) {
            names(splittbl) = splitvars
            spsspivottable.Display(
                splittbl,
                title = gtxtf("Table of Splits"),
                templateName = "SPLITTABLE",
                isSplit = FALSE,
                hidecoldimtitle=TRUE,
                rowlabels=as.character(seq(1:nrow(splittbl))),
                caption = gtxtf("Use this table to identify splits in plots")
            )
        }

        basename = gtxtf("%s Normality Univariate Plot", splitprefix)
        spssRGraphics.SetGraphicsLabel(basename)
        for (pl in uniplots) {
            res = suppressWarnings(mvn(dta, univariatePlot = pl ))
        }
        
        # multivariate plots
        basename = gtxtf("%s Normality Multivariate Plot", splitprefix)
        spssRGraphics.SetGraphicsLabel(basename)
        if (nvars > 1) {
            for (pl in multivarplots) {
                if (pl == "qq") {
                    res = mvn(data=dta, scale=scaledata, multivariatePlot="qq")
                } else if (pl == "contour") {
                    # if (nvars > 2) {
                    #     warns$warn(gtxt("Contour plots are limited to two variables.  Plotting first two"), dostop = FALSE)
                    # }
                    res = mvn(dta[1:2], scale=scaledata, multivariatePlot="contour")
                } else if (pl == "persp") {
                    # if (nvars > 2) {
                    #     warns$warn(gtxt("Perspective plots are limited to two variables.  Plotting first two"), dostop=FALSE)
                    # }
                    res =mvn(dta[1:2], scale=scaledata, multivariatePlot = "persp")
                }
            }
        }
    }
    spssdata.CloseDataConnection()


    warns$display(inproc=TRUE)
    spsspkg.EndProcedure()
    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})

}

mvtestresults <- function(dta, mvntests, bootstrapreps, scaledata) {
    # return a list with a data frame containing the multivariate test results and the d.f. for dh or NULL
    
    # dta is the data frame to analyze
    # mvn is a list of the mv test parameters
    # bootstrapreps is the number of replications for the energy e test
    
    # Each result is a list of test name, statistic, and p value
    # test result structure returned by mvn vary by test type :-(
    
    mt = data.frame(Test=0, Statistic=0, "p value"= 0)
    row = 1
    dhdegf = NULL # for Doornik-Hansen
    for (item in mvntests) {
        mvt = suppressWarnings(mvn(dta, item, scale=scaledata, desc=FALSE, R=bootstrapreps, univariateTest="AD", subset=NULL)) # ignore univariate
        # can't use rbind because column names vary by test

        if (item == "dh") {
            dhdegf = mvt$multivariateNormality[1, 3]
            mvt$multivariateNormality[1,3] = mvt$multivariateNormality[1,4] # squeeze out df value
        }
        if (item != "mardia") {
            mt[row,] = mvt$multivariateNormality[1, c(1,2,3)]
            mt[row, 2] = round4(mt[[row, 2]])
            mt[row, 3] = round4(as.numeric(mt[[row, 3]]))
            if (mt[row, 3] == 0) {
                mt[row, 3] = "<.001"
            }
            row = row + 1
        } else {  # mardia has test and p value as factors
            mt[row, 1] = mvt$multivariateNormality$Test[1]
            mt[row, 2] = round4(as.numeric(levels(mvt$multivariateNormality$Statistic))[[1]])
            mt[row, 3] = round4(as.numeric(levels(mvt$multivariateNormality$"p value"))[[1]])
            if (mt[row, 3] == 0) {
                mt[row, 3] = "<.001"
            }
            row = row + 1
            mt[row, 1] = mvt$multivariateNormality$Test[2]
            mt[row, 2] =round4(as.numeric(levels(mvt$multivariateNormality$Statistic[[1]])[[2]]))
            mt[row, 3] = "--"
            row = row + 1
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
    oobj <- spsspkg.Syntax(templ=list(
                spsspkg.Template("VARIABLES", subc="", ktype="varname", var="variables", islist=TRUE),
                
                spsspkg.Template("MVNTESTS", subc="OUTPUT", ktype="str", var="mvntests", 
                    vallist=list("mardia", "hz", "royston", "dh", "energy"), islist=TRUE),
                spsspkg.Template("UNIVARIATETESTS", subc="OUTPUT", ktype="str", var="univariatetests", 
                    vallist=list("sw", "cvm", "lillie", "sf", "ad"), islist=TRUE),
                spsspkg.Template("BOOTSTRAPREPS", subc="OUTPUT", ktype="int", var="bootstrapreps", islist=FALSE),
                spsspkg.Template("SCALEDATA", subc="OUTPUT", ktype="bool", var="scaledata", islist=FALSE),
                spsspkg.Template("UNIPLOTS", subc="OUTPUT", ktype="str", var="uniplots",
                    vallist=list("qq", "histogram", "box", "scatter"),  islist=TRUE),
                spsspkg.Template("MULTIVARPLOTS", subc="OUTPUT", ktype="str", var="multivarplots",
                    vallist=list("qq", "persp", "contour"),  islist=TRUE),
                spsspkg.Template("DESCRIPTIVES", subc="OUTPUT", ktype="bool", var="desc", islist=FALSE),

                spsspkg.Template("IDVAR", subc="OUTLIERS", ktype="existingvarlist", var="idvar", islist=FALSE),
                spsspkg.Template("NOUTLIERS", subc="OUTLIERS", ktype="int", var="noutliers", islist=FALSE),
                spsspkg.Template("OUTLIERDETECTION", subc="OUTLIERS", ktype="str", var="outlierdetection",
                    vallist=list("quan", "adj"),  islist=FALSE)
                ))

    if ("HELP" %in% attr(args,"names"))
        #writeLines(helptext)
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

#}
