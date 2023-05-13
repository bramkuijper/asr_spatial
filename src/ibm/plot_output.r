#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("tidyverse", warn.conflicts=F))
suppressPackageStartupMessages(library("jsonlite", warn.conflicts=F))
suppressPackageStartupMessages(library("patchwork", warn.conflicts=F))

# from a list of values like x1, x2, x3
# create a reasonable variable name, like x
make.var.name <- function(vars) {
    
    var1 <- vars[[1]]

    return(gsub(pattern="[_0-1]",replacement="",x=var1))
}

find.params <- function(filename) {

    f <- readLines(filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

xvar <- "time"

# structure of the graph in JSON
plots_json <- paste0('[
    {"xvar" : "',xvar,'",
    "yvar" : ["mean_T_f","mean_T_m"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["ss_T_f","ss_T_m"]
    },
    {"xvar" : "',xvar,'",
    "yvar" : ["mean_Tb_f","mean_Tb_m"]
    },
    {"xvar" : "',xvar,'",
    "yvar" : "mean_phen"
    },
    {"xvar" : "',xvar,'",
    "yvar" : "var_phen"
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["ss_Tb_f","ss_Tb_m"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["n_care_f","n_care_m"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["n_mate_f","n_mate_m"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["n_juv_f","n_juv_m"],
        "ylim" : [0,250]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : "juv_sr",
        "ylim" : [ -0.05, 1.05]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["n_ad_tot_f","n_ad_tot_m"]
    },
    {
        "xvar" : "',xvar,'",
        "yvar" : ["n_tot_f","n_tot_m"]
    }
]')

#file.name <- "sim_asr_20230509_101417900895_0"

# if the file name not provided raise error
if (!exists("file.name"))
{
    
    # get command line arguments
    args = commandArgs(trailingOnly=TRUE)
    
    # give an error message if you do not provide it with a simulation file name
    if (length(args) < 1)
    {
        print("provide a simulation file name")
        stop()
    }
    
    file.name <- args[[1]]
}

param.line <- find.params(file.name)

data.tibble <- read_delim(file=file.name
        ,delim=";"
        ,n_max=param.line-1
        ,col_names=T)

if (nrow(data.tibble) > 50000)
{
    data.tibble <- data.tibble[data.tibble$time %% 10 == 0,]
}

# get the parameters
data.tibble.params <- read_delim(file=file.name
        ,delim=";"
        ,skip=param.line
        ,col_names=c("name","value")
        )

# transpose the tibble with the parameters
params <- data.tibble.params %>% pivot_wider(
        names_from = name
        ,values_from = value)

data.tibble <- data.tibble %>% mutate(
        n_ad_tot_f= n_mate_f + n_care_f,
        n_ad_tot_m= n_mate_m + n_care_m,
        n_tot_f= n_ad_tot_f + n_juv_f,
        n_tot_m= n_ad_tot_m + n_juv_m,
        juv_sr = n_juv_m / (n_juv_f + n_juv_m)
        )


plot.structure <- fromJSON(plots_json, simplifyVector = F)

plot.structure.l <- length(plot.structure)

# list with all the plots
plot.list <- list(rep(NA,times=plot.structure.l))

plot.list.idx <- 1

# first plot stuff for the tibble as a whole
for (plot_struct_idx in 1:plot.structure.l)
{
    # get the (potential list of) y variable(s)
    # as this is a list and hence highly structured
    # hence, try to flatten it
    yvar <- unlist(plot.structure[[plot_struct_idx]]$yvar)

    if (length(yvar) > 1)
    {
        yvar_name <- make.var.name(yvar)
        yvar_values <- paste0(yvar_name,"_values")

        sub.data <- pivot_longer(data=data.tibble
                ,cols=yvar
                ,names_to=yvar_name
                ,values_to=yvar_values)

        plot.list[[plot.list.idx]] <- ggplot(data=sub.data
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=yvar_values)) + geom_line(mapping=aes_string(colour=yvar_name))
    } else {
        plot.list[[plot.list.idx]] <- ggplot(data=data.tibble
                ,mapping=aes_string(x=plot.structure[[plot_struct_idx]]$xvar
                        ,y=plot.structure[[plot_struct_idx]]$yvar)) + geom_line()
    }

    # add ylim
    if ("ylim" %in% names(plot.structure[[plot_struct_idx]]))
    {
        plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + ylim(
                unlist(
                        plot.structure[[plot.list.idx]]$ylim)
                )
    }
    
    plot.list[[plot.list.idx]] <- plot.list[[plot.list.idx]] + theme_classic()
    
    plot.list.idx <- plot.list.idx + 1
}

title <- ""

if (exists("params") && "mu_juv_0_f" %in% params)
{
    title <- paste0(
            "mu_juv_0_f: ",params["mu_juv_0_f"],
            ", mu_juv_1_f: ",params["mu_juv_1_f"],
            ", mu_juv_0_m: ",params["mu_juv_0_m"],
            ", mu_juv_1_m: ",params["mu_juv_1_m"],
            ", sigma: ",params["sigma"])
}


wrap_plots(plot.list,ncol=1) + plot_annotation(
        title=title)

file.name <- paste0("graph_",basename(file.name),".pdf")

ggsave(file.name,height= 3 * plot.structure.l)
