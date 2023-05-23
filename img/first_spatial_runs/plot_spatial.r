library("tidyverse")


the.data <- read_delim(file="summary.csv", delim=";")

find.params <- function(filename) {

    f <- readLines(con=filename)

    seq.rev <- rev(seq(1,length(f),1))

    for (line_i in seq.rev)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }
}

any.mu <- function(params)
{
    mus <- params[,c("mu_juv_0_f","mu_juv_0_m","mu_juv_1_f","mu_juv_1_m")]

    sum_mu_f <- mus$mu_juv_0_f + mus$mu_juv_1_f
    sum_mu_m <- mus$mu_juv_0_m + mus$mu_juv_1_m

    return(c(sum_mu_f,sum_mu_m))
}

get.all.data <- function(data.set)
{
    all.data <- NULL

    for (row in 1:nrow(data.set))
    {
        file.name <- data.set[row,] %>% pull("file")

        param.line <- find.params(filename=file.name)

        long.data <- read_delim(file=file.name,n_max=param.line-1,delim=";")

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

        # extract summed mortality values
        muvals <- any.mu(params=params)

        long.data <- long.data %>% mutate(
                mu_f = muvals[1]
                ,mu_m = muvals[2]
                ,replicate = row)

        all.data <- bind_rows(all.data,long.data)
    }

    return(all.data)
} # end function get.all.data()

all.the.bloody.data <- get.all.data(data.set=the.data)

all.the.bloody.data <- all.the.bloody.data %>% mutate(
        replicate_f = factor(replicate)
        ,time_is_of_the_essence=time
        )

all.the.bloody.data <- arrange(all.the.bloody.data, replicate_f, time)

select.dat <- filter(all.the.bloody.data,time==max(all.the.bloody.data$time))

ggplot(data=all.the.bloody.data
        ,mapping=aes(x=mean_T_f,y=mean_T_m)) +
    geom_path(mapping=aes(colour=replicate_f)) +
    geom_point(data=select.dat,
            mapping=aes(x=mean_T_f, y=mean_T_m)) +
    facet_grid(mu_m~mu_f, labeller=label_both) +
    xlab("Female care") +
    ylab("Male care") +
    theme(legend.position="none")



ggsave("plot_T.pdf")

ggplot(data=all.the.bloody.data
        ,mapping=aes(x=mean_phen_f,y=mean_phen_m)) +
    geom_path(mapping=aes(colour=replicate_f)) +
    facet_grid(mu_m~mu_f)


ggsave("plot_phen.pdf")
select.dat <- filter(all.the.bloody.data,time==max(all.the.bloody.data$time))

ggplot(data=all.the.bloody.data
        ,mapping=aes(x=mean_Tb_f,y=mean_Tb_m)) +
    geom_path(mapping=aes(colour=replicate_f)) +
    facet_grid(mu_m~mu_f) +
    geom_point(data=select.dat,
            mapping=aes(x=mean_Tb_f, y=mean_Tb_m)) 

ggsave("plot_Tb.pdf")




