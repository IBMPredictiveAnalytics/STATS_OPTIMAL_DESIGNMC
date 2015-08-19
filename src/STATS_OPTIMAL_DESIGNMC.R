#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2015
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.0

# history
# 02-feb-2015 - original version



optdesmc <- function(varnames=NULL,  frml=NULL, factors="no",
    mixtures="no", lows=-Inf, highs=Inf, centers=NULL,
    nlevels=NULL, roundtos=NULL, ncand=NULL,
    constraintfunc=NULL, mixturesum=1,
	model="linear", constant=TRUE,
	ntrials=NULL, criterion="d", center=FALSE, 
	initial="random", repeats=5, designalg="exact",
	outputdataset, confounding=TRUE) {

    setuplocalization("STATS_OPTIMAL_DESIGNMC")	

    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Optimal Design - Monte Carlo")
    warningsprocname = gtxt("Optimal Design - Monte Carlo: Warnings")
    omsid="STATSOPTDESIGNMC"
    warns = Warn(procname=warningsprocname,omsid=omsid)
    
    # The bool type with lists is broken prior to Statistics 23, so patch this up with str workaround
    factors = fixtype(factors, warns)
    mixtures = fixtype(mixtures, warns)

    tryCatch(library(AlgDesign), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","AlgDesign"),
            dostop=TRUE)
    }
    )
    hasformula = !is.null(frml)
	# validation
    spec = validate(varnames, frml, factors, nlevels, lows, highs, centers, roundtos, mixtures, warns)
    variables = spec[,1]
    factors = spec[,7]
    vlevels = spec[,5]
    if (any(as.logical(mixtures))) {
        constant = FALSE
    }

	alldatasets = spssdata.GetDataSetList()
	if ("*" %in% alldatasets) {
		warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
			dostop=TRUE)
	}
	if (tolower(outputdataset) %in% tolower(alldatasets)) {
		warns$warn(gtxt("The output dataset name is already in use"),
			dostop=TRUE)
	}
    # If using a constraint function, load and check for required name "dfilter"
    # source into a separate environment to minimize conflict possibilities
    if (!is.null(constraintfunc)) {
        funcenv = new.env()
        tryCatch(sys.source(constraintfunc, envir=funcenv),
            error=function(e) {
                warns$warn(e$msg, dostop=TRUE)
            }
        )

        if (!exists("dfilter", envir=funcenv)) {
            warns$warn(gtxt("The constraint function source file does not define a function named dfilter"),
                dostop=TRUE)
        } else {
            dfilter = get("dfilter", envir=funcenv)
        }
    }
	criterion = toupper(criterion)
	

	factorlist = buildfactorlist(variables, factors, warns)
	frml = genfrml(variables, frml, model, constant, factorlist)
    arglist = list(
        frml=as.formula(frml),
        data = spec,
        nTrials = ntrials,
        approximate = designalg == "approx",
        evaluateI = FALSE,
        mixtureSum=mixturesum,
        criterion = criterion,
        RandomStart = ifelse(initial == "random", TRUE, FALSE),
        nRepeats = repeats,
        DFrac=1, 
        CFrac=1
    )
    if (is.null(ntrials)) {
        arglist["nTrials"] = NULL # optMorteCarlo won't tolerate null value
    }
    # constraint function for montecarlo
    if (!is.null(constraintfunc)) {
        arglist[["constraints"]] = dfilter
    }
    if (!is.null(ncand)) {
        arglist["nCand"] = ncand
    }

    res = tryCatch(do.call(optMonteCarlo, arglist),
           error = function(e) {
               warns$warn(e, dostop=TRUE)
           }
        )

    if (!any(sapply(res$design[, 2:ncol(res$design)], is.numeric))) {
        center = FALSE
    }

    # print results
    # 
    displayresults(res, spec, hasformula, variables, vlevels, model, constraintfunc,
		constant, criterion, center, initial, repeats, designalg, ncand,
		outputdataset, frml, confounding, warns)

	# create dataset
	spsspkg.EndProcedure()
	gendataset(res, outputdataset, dsdict, variables, factorlist, warns)
	warns$display(inproc=FALSE)

	
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}

validate = function(varnames, frml, factors, nlevels, lows, highs, centers, roundtos, mixtures, warns) {
    # return data frame of design specifications
    # varnames is a list of variable names
    # frml is a formula for the design
    # exactly one of varnames and frml must be supplied
    # types is an optional list of types (default=factor)
    # low, high, center, and roundto are variable parameters
    # mixture is TRUE/FALSE for each variable
    
    if (is.null(varnames) && is.null(frml)) {
        warns$warn(gtxt("Either VARNAMES or FORMULA must be specified."),
            dostop=TRUE)
    }
    if (!is.null(frml)) {
        if (substr(frml, 1, 1) != "~") {
            frml = paste("~", frml, sep="")
        }
        varnamesx = tryCatch(all.vars(formula(frml)),
            error = function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )
        if (!is.null(varnames)) {
            if (length(setdiff(varnames, varnamesx)) > 0) {
                warns$warn(gtxt("Formula variable names are inconsistent with VARNAMES specification"),
                    dostop=TRUE)
            }
        } else {
            varnames = varnamesx
        }
    }
    nvars = length(varnames)
    lcvarnames = tolower(varnames)
    if (length(union(lcvarnames, lcvarnames)) != nvars) {
        warns$warn(gtxt("Duplicate variable names were found"), dostop=TRUE)
    }
    

    factors = unlist(fixup(nvars, factors, gtxt("factors"), warns))
    lows = fixup(nvars, lows, gtxt("low values"), warns)
    highs = fixup(nvars, highs, gtxt("high values"), warns)
    centers = fixup(nvars, centers, gtxt("center values"), warns)
    nlevels = fixup(nvars, nlevels, gtxt("number of levels"), warns)
    roundtos = fixup(nvars, roundtos, gtxt("round values"), warns)
    mixtures = unlist(fixup(nvars, mixtures, gtxt("mixture specification"), warns))
    
    for (i in 1:nvars) {
        if (lows[[i]] >= highs[[i]] || (centers[[i]] < lows[[i]] || centers[[i]] > highs[[i]])) {
            warns$warn(gtxtf("Invalid low, high, or center specifications for variable %s", varnames[[i]]),
                       dostop=TRUE)
        }
    }

    return(data.frame(var=unlist(varnames), lows, highs, centers, 
        nlevels, roundtos, factors, mixtures))
}

fixup = function(nvars, item, label, warns) {
    # check and adjust item
    # nvars is the required length
    # item is the list to process
    # label is for use in an error message

    if (length(item) == 1) {
        item = rep(item, nvars)
    }
    if (length(item) != nvars) {
        warns$warn(gtxtf("An input item has the wrong length: %s", label), dostop=TRUE)
    }
    return(item)
}

genfrml = function(variables, frml, model, constant, factorlist) {
	# generate and return a formula for the design
	# variables are the variables in the model
	# model indicates polynomial order - applied to all
    # nonfactor variables
	# constant is TRUE if a constant term should be included
	# factorlist provides numbers of the terms that are factors

    if (!is.null(frml)) {
        if (substr(frml,1,1) != "~") {
            frml = paste("~", frml, sep="")
        }
        return(frml)
    }
	if (model == "linear") {
		part = paste(variables, collapse="+")
	} else {
		factors = variables[factorlist]
		nonfactors = setdiff(variables, factors)
		if (length(nonfactors) > 0) {
			part = paste(nonfactors, collapse=",")
            if (model == "cubics") {
                model="cubicS"
            }
			part = paste(model,"(", part, ")", sep="", collapse="")
		} else {
			part = NULL
		}

		if (length(factors) > 0) {
			f = paste(factors, collapse="+")
			if (is.null(part)) {
				part = f
			} else {
				part = paste(part, f, sep="+", collapse="")
			}
		}
	}
	if (!constant) {
		part = paste(part, "-1", collapse="")
	}

	return(paste("~", part, collapse=""))
}

buildfactorlist = function(variables, factors, warns) {
	# return vector of indexes of factors
	

	cvars = c()
    count = 1
    for (i in 1:length(variables)) {
        if (factors[[i]]) {
            cvars[count] = i
            count = count + 1
        }
    }
	return(cvars)
}
			

rounder <- function(s) {
    if (!is.null(s)) {
        s = round(s, 5)
    } else {
        s = "."
    }
    return(s)
}

gendataset = function(res, outputdataset, dsdict, variables, factorlist, warns) {
	# create output dataset holding the design
    # besides the design variable values, it has a repetitions variable

	varspec = list()
    varnames = names(res$design)
    
    # fix (within reason) any duplicate name with the one added by the calculation
    # duplicate names for design variables are caught earlier by
    # the function
    
	varnameslc = lapply(varnames, tolower)
    while (varnameslc[1] %in% varnameslc[2: length(varnameslc)]) {
        varnames[1] = paste(varnames[1],".", sep="")
        varnameslc[1] = tolower(varnames[1])
    }
    

	for (v in 1:ncol(res$design)) {
        if (is.numeric(res$design[[v]])) {
            format = "F8.2"
            mlevel = "scale"
            len = 0
        } else {
            ll = levels(res$design[[v]])
            len = max(sapply(ll, nchar),1) * 3  # worst case Unicode expansion
            format = paste("A", len, sep="")
            mlevel = "nominal"
        }

		vspec = c(varnames[v], "", len, format, mlevel)
		varspec[[v]] = vspec
	}

	dsdict = do.call(spssdictionary.CreateSPSSDictionary, varspec)

	spssdictionary.SetDictionaryToSPSS(outputdataset, dsdict)
	spssdata.SetDataToSPSS(outputdataset, res$design)
	spssdictionary.EndDataStep()

}


displayresults = function(res, data, hasformula, variables, vlevels, model, constraintfunc,
		constant, criterion, center, initial, repeats, designalg, ncand,
		outputdataset, frml, confounding, warns) {

    StartProcedure(gtxt("Optimal Design"), "STATSOPTDESIGN")
    # summary results
    scaption = gtxt("Computations done by R package AlgDesign")
    lbls = c(
        gtxt("Number of Candidate Points"),
		gtxt("Actual Number of Trials"),
        gtxt("Model"),
		gtxt("Expanded Model"),
		gtxt("Include Constant"),
        gtxt("Initial Design"), 
        gtxt("Criterion"), 
		gtxt("Center Data for Evaluation"),
        gtxt("Design Generation Repeats"), 
        gtxt("Algorithm"), 
        gtxt("Output Dataset")
    )
	    expandedmodel = frml
    if (hasformula) {
        model = gtxt("formula")
    }
    ntrials = nrow(res$design)

    vals = c(
        ifelse(is.null(ncand), gtxt("not specified"), ncand),
		ntrials,
        model,
		expandedmodel,
        ifelse(constant, gtxt("Yes"), gtxt("No")),
		initial,
		criterion,
        ifelse(center, gtxt("Yes"), gtxt("No")),
		repeats,
		designalg,
		outputdataset
    )

    # settings and result summary
    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Settings and Results Summary"),
        collabels=c(gtxt("Summary")), templateName="OPTDESSUMMARYMC", outline=gtxt("Summary"),
        caption = gtxt("Results computed by AlgDesign package"))
      
    names(data) = c(
        gtxt("Variable"), gtxt("Low"), gtxt("High"), gtxt("Center"),
        gtxt("Levels"), gtxt("Round to Decimals"), gtxt("Factor"), gtxt("Mixture Variable")
    )
    data[7] = sapply(data[7], function(x) {ifelse(x, gtxt("Yes"), gtxt("No"))})
    data[8] = sapply(data[8], function(x) {ifelse(x, gtxt("Yes"), gtxt("No"))})
    caption = ifelse(!is.null(constraintfunc),
        gtxtf("Constraint function from file: %s", constraintfunc),
        gtxt("Constraint function: None")
    )
    spsspivottable.Display(
        data, 
        title=gtxt("Variable Specifications"),
        caption=caption
        )
    # evaluation
	# eval.design does not work if all terms are factors and centering is specified

    emsg = ""
	designeval = 
        tryCatch(eval.design(frml=frml, design=res$design, confounding=TRUE, variances=TRUE,
		center = center),
        error= function(e) {
            emsg <<- e$message
            warns$warn(e$message, dostop=FALSE)
        }
        )

	lbls = c("D",
		"A (Average Coefficient Variance)",
		"Ge (Minimax Normalized Variance Efficiency)",
		"Dea (Lower Bound on D Efficiency)",
		gtxt("Determinant"),
		gtxt("Diagonality"),
		gtxt("Geometric Mean of Coefficient Variances")
	)
	vals = c(res$D,
		res$A,
		res$Ge,
		res$Dea,
		designeval$determinant,
		designeval$diagonality,
		designeval$gmean.variances
	)
    lbls = lbls[1:length(vals)]  # singular design loses some statistics
	spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), 
		title = gtxt("Design Evaluation"), caption=emsg,
	collabels=c(gtxt("Statistics")), templateName="OPTDESEVAL", outline=gtxt("Evaluation"))

	# confounding
	if (confounding && !is.null(designeval)) {
		condf = data.frame(designeval$confounding)
		names(condf) = row.names(condf)
		spsspivottable.Display(condf, title=gtxt("Confounding Matrix"),
			templateName="OPTDESCONF", outline=gtxt("Confounding"),
			caption= gtxt("Each column gives the coefficients of that variable regressed on the others")
		)
	}

}

fixtype = function(alist, warns) {
    # return a list of logical values
    # workaround function for ECM211559
    # alist is a list of str values in "yes", "no", "true", "false"
    
    isvalid = all(alist %in% list("yes", "true", "no", "false"))
    if (!isvalid) {
        warns$warn(gtxt("A yes/no keyword has an invalid value"), dostop=TRUE)
    }
    return(lapply(alist, function(x) {x == "yes"||x=="true"}))

}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_OPTIMAL_DESIGNMC"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_OPTIMAL_DESIGNMC"))
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
    lcl = list2env(lcl) # makes this list into an environment

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
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        } else {
            procok = TRUE
        }
        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
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

Run<-function(args){
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("VARNAMES", subc="", ktype="varname", var="varnames", islist=TRUE),
        spsspkg.Template("FORMULA", subc="", ktype="literal", var="frml"),
        spsspkg.Template("FACTORS", subc="", ktype="str", var="factors", islist=TRUE), # should be bool
        spsspkg.Template("LEVELS", subc="", ktype="int", var="nlevels", islist=TRUE,
            vallist=list(2)),
        spsspkg.Template("MIXTURES", subc="", ktype="str", var="mixtures", islist=TRUE), # should be bool
        spsspkg.Template("LOWS", subc="", ktype="float", var="lows", islist=TRUE),
        spsspkg.Template("HIGHS", subc="", ktype="float", var="highs", islist=TRUE),
        spsspkg.Template("CENTERS", subc="", ktype="float", var="centers", islist=TRUE),
        spsspkg.Template("ROUNDTOS", subc="", ktype="int", var="roundtos", islist=TRUE),
        spsspkg.Template("CONSTRAINTFUNC", subc="", ktype="literal", var="constraintfunc"),
        spsspkg.Template("MIXTURESUM", subc="", ktype="float", var="mixturesum"),
        
        spsspkg.Template("MODEL", subc="",  ktype="str", var="model", 
			vallist=list("linear", "quad", "cubic", "cubics")),
        spsspkg.Template("TRIALS", subc="", ktype="int", var="ntrials"),
		spsspkg.Template("CONSTANT", subc="", ktype="bool", var="constant"),
		
        spsspkg.Template("CRITERION", subc="OPTIONS", ktype="str", 
            var="criterion", vallist=list("d", "a")),   # "i" does not work - bug in package
        spsspkg.Template("NUMCAND", subc="OPTIONS", ktype="int",
            var="ncand", vallist=list(2)),
		spsspkg.Template("CENTER", subc="OPTIONS", ktype="bool", 
            var="center"),
		spsspkg.Template("INITIAL", subc="OPTIONS", ktype="str", var="initial",
			vallist=list("random","nullification")),
        spsspkg.Template("REPEATS", subc="OPTIONS", ktype="int", 
            var="repeats", vallist = list(1)), 
        spsspkg.Template("DESIGNALG", subc="OPTIONS", ktype="str", 
            var="designalg", vallist=list("approx", "exact")),
			
		spsspkg.Template("CONFOUNDING", subc="OUTPUT", ktype="bool", var="confounding"),
			
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", 
        var="outputdataset"),		
        
        #spsspkg.Template("CONSTRAINTS", subc="OUTPUT", ktype="bool", var="constraintsds"),
        #spsspkg.Template("ALLVARS", subc="OUTPUT", ktype="bool", var="allvars"),
		
		spsspkg.Template("", subc="OBJBOUNDS", ktype="literal", var="objbounds", islist=TRUE)
    ))        
if ("HELP" %in% attr(args,"names"))
    helper(cmdname)
else
    res <- spsspkg.processcmd(oobj,args,"optdesmc")
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
