#' Create a new ppfunnel plot
#'
#' `ppfunnel()` initializes a ppfunnel plot. It can be used to declare the input
#' data frame, outcome type, and indicator of interest.
#' @param data A data frame including summary information used for creating a
#' graphic. See `indiv.data` for details on required columns.
#' @param indiv.data logical. If `TRUE` (default) then each row of `data`
#' corresponds to an individual (e.g., patient or hospital discharge), and columns
#' `"O"` (observed outcome), `"E"` (expected outcome), `"ID"` (provider ID) are
#' required. If `FALSE` then each row of `data` corresponds to a provider (e.g.,
#' hospital or facility), and `"O"` (observed outcome sum) and `"E"` (expected
#' outcome sum) are required. When `outcome = "bernoulli"`, an additional column
#' `"probs"`, each element of which is a list of expected outcomes for a certain
#' provider, is required.
#' @param outcome A character string indicating the outcome type.
#' If `"poisson"` then the outcome is a count (i.e., 0, 1, 2, ...).
#' If `"bernoulli"` then the outcome is a binary indicator (i.e., 0 or 1).
#' Abbreviation supported.
#' @param indicator A character string indicating the indicator of interest.
#' `"ISR"` (indirectly standardized ratio, default) and `"Prop"` (proportion) are
#' supported.
#' @param target numeric. It specifies the desired expectation for providers
#' considered "in-control". If `indicator = "ISR"` then `target` defaults to 1.
#' @param test A character string indicating the type of statistical test used
#' for calculating the precision and control limits. `"score"` (score test,
#' default) and `"exact"` (exact test) are supported.
#' @param method A character string indicating the method used for overdispersion
#' due to incomplete risk adjustment. Options include `"FE"` (fixed effects
#' modeling omitting overdispersion), `"MO"` (multiplicative overdispersion without
#' Winsorization), `"MOW"` (multiplicative overdispersion with Winsorization),
#' `"MODW"` (multiplicative overdispersion with debiased Winsorization), `"AO"`
#' (additive overdispersion), `"indivEN.exact"` (individualized empirical null
#' without Taylor approximation), `"indivEN.approx"` (individualized empirical
#' null with Taylor approximation), `"indivEN.0meanapprox"` (individualized
#' empirical null with Taylor approximation and zero mean assumption).
#' @param q.winsor A numeric vector specifying the lower and upper thresholds for
#' Winsorization, with `c(0.05, 0.95)` as the default.
#' @param levels A numeric vector with 1 to 5 elements specifying the p-value(s)
#' for displaying control limits, possibly multiple pairs when the length of
#' `levels` is greater than 1.
#' @param flags A vector of three character strings for flagging providers, with
#' `"c("lower", "expected", "higher")"` as the default.
#' @param colors A vector of three character strings specifying the colors for
#' flagging.
#' @param xlab,ylab Titles for the `x` and `y` axes, respectively.
#' @param xlim,ylim Numeric vectors of length two providing ranges of the axes
#' passed on to `scale_x_continuous` and `scale_y_continuous` of `ggplot2`.
#' `NULL` uses the default scale range.
#' @param xexpand,yexpand Numeric vectors of length two for adding some padding
#' to data to ensure that they are placed some distance away from the axes.
#' `xexpand` and `yexpand` will be passed on to `scale_x_continuous` and
#' `scale_y_continuous` of `ggplot2`.
#' @param xbreaks,ybreaks Numeric vectors passed on to `scale_x_continuous` and
#' `scale_y_continuous` of `ggplot2`.
#' @param ptsize A positive number specifying the size of the points, passed on
#' to `geom_point` of `ggplot2`.
#' @param lwd A positive number specifying the line width, passed on to
#' `geom_line` of `ggplot2`.
#' @param title A character string specifying the text of the title, passed on
#' to `labs` of `ggplo2`.
#' @param ... A list of additional arguments passed on to `theme` of `ggplot2`.
#'
#' @details
#' `ppfunnel()` is used to create funnel plots for profiling health care
#' providers. The default `method = "FE"` for constructing control limits is
#' based on a model with fixed provider effects (see
#' [He et al., 2013](https://doi.org/10.1007/s10985-013-9264-6) for an example
#' of Bernoulli outcomes). Methods `"MO"`, `"MOW"`, `"MODW"`, `"AO"` are based on
#' [Spiegelhalter (2004)](https://doi.org/10.1002/sim.1970). Methods based on
#' the individualized empirical null approach were discussed in
#' (cite online preprint).
#'
#' @return A `ggplot` object of class `"gg"` and `"ggplot"`.
#'
#' @examples
#' library(ggplot2)
#' mytheme <- theme(legend.justification=c(1,1), legend.position=c(1,1),
#' legend.box="horizontal", legend.text=element_text(size=14),
#' axis.title=element_text(size=14),
#' axis.text=element_text(size=14),
#' plot.title=element_text(hjust=0), text=element_text(size=13))
#' # Poisson outcome with method = "FE"
#' SHR %>% ppfunnel(indiv.data=F, legend.justification=c(1,1),
#' legend.position=c(1,1), legend.box="horizontal",
#' legend.text=element_text(size=14), axis.title=element_text(size=14),
#' axis.text=element_text(size=14), plot.title=element_text(hjust=0),
#' text=element_text(size=13))
#' # alternative use
#' SHR %>% ppfunnel(indiv.data=F) + mytheme
#' # Poisson outcome with method = "indivEN.0meanapprox"
#' SHR %>% ppfunnel(indiv.data=F, method="indivEN.0meanapprox") + mytheme
#'
#' @author Wenbo Wu and Kevin He
#'
#' @references
#' Spiegelhalter DJ (2004) "Funnel plots for comparing institutional
#' performance", *Stat Med* 24(8): 1185--1202 <https://doi.org/10.1002/sim.1970>.
#'
#' He K, Kalbfleisch JD, Li Y, and Li Y (2013) "Evaluating hospital readmission
#' rates in dialysis facilities; adjusting for hospital effects",
#' *Lifetime Data Anal* 19(4): 490--512 <https://doi.org/10.1007/s10985-013-9264-6>.
ppfunnel <- function(data, indiv.data=TRUE, outcome="poisson", indicator="ISR",
                     target=1, test="score", method="FE", q.winsor=c(0.05, 0.95),
                     levels=c(0.05, 0.01, 0.1),
                     flags=c("lower", "expected", "higher"),
                     colors=c('blue', '#E69F00', 'red'), xlab=NULL, ylab=NULL,
                     xlim=NULL, ylim=NULL, xexpand=waiver(), yexpand=waiver(),
                     xbreaks=waiver(), ybreaks=waiver(), ptsize=0.8, lwd=0.6,
                     title=waiver(), ...) {
  if (!is.numeric(target)) stop("'target' is not numeric!")
  if (!is.character(flags)) stop("'flags' is not character-valued!")
  if (length(flags)!=3) stop("length of 'flags' unequal to 3!")
  if (method %in% c("MOW", "MODW") &&
      (!is.numeric(q.winsor) ||
       length(q.winsor)!=2 ||
       any(range(q.winsor)>1) ||
       any(range(q.winsor)<0))) stop("'q.winsor' is invalid!")
  if (indiv.data) {
    if (!"ID"%in%names(data)) stop("'ID' not in 'data'!")
    m <- length(unique(data %>% pull(ID)))
  } else m <- NROW(data) # provider count
  level <- levels[1L] # primary level for flagging
  levels <- sort(levels) # levels for funnels
  n.level <- length(levels) # level count
  if (n.level>5) stop("'levels' greater than 5!")
  # Poisson outcome
  if (grepl("^pois", tolower(outcome))) { # lowercased 'outcome' starts with "pois"
    if (!"data.frame"%in%class(data)) stop("Invalid class of 'data'!")
    if (indicator=="ISR") { # indirectly standardized ratio
      if (!all(c("O","E") %in% names(data))) stop("'O' or 'E' not in 'data'!")
      if (indiv.data)
        data <- data %>% group_by(ID) %>%
          summarise(O=sum(O), E=sum(E), .groups="drop")
      if (target<0) stop("'target' is negative!")
      if (test=="score") {
        if (method=="FE") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E),
                   prec=E, prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                               flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               dev=data %>% pull(E) %>% `^`(-1) %>% `*`(target) %>% sqrt %>%
                                 rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0), ulim=target+dev,
                   level=factor(level))
        } else if (method=="indivEN.exact") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          res.indivEN <-
            indivEN.pois.exact(tbl.flag$stat, tbl.flag$E, target=target)
          tbl.flag <- tbl.flag %>%
            mutate(mean=sqrt(target)*res.indivEN[2]*sqrt(E),
                   sd=sqrt(res.indivEN[3]+
                             target*(res.indivEN[3]^4-res.indivEN[3]^2)*E),
                   stat=(stat-mean)/sd, prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]), flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               const=data %>% pull(E) %>% `^`(-1) %>% `*`(target) %>% sqrt %>%
                                 rep(each=n.level),
                               mean=tbl.flag %>% pull(mean) %>% rep(each=n.level),
                               sd=tbl.flag %>% pull(sd) %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(llim=pmax(target+const*(mean-c.val*sd),0),
                   ulim=target+const*(mean+c.val*sd),
                   level=factor(level))
        } else if (method=="indivEN.approx") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          res.indivEN <-
            indivEN.pois.approx(tbl.flag$stat, tbl.flag$E, target=target)
          tbl.flag <- tbl.flag %>%
            mutate(mean=res.indivEN[2],
                   sd=sqrt(1+target*res.indivEN[3]*E),
                   stat=(stat-mean)/sd, prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]), flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               const=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target) %>% sqrt %>% rep(each=n.level),
                               mean=tbl.flag %>% pull(mean) %>% rep(each=n.level),
                               sd=tbl.flag %>% pull(sd) %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(llim=pmax(target+const*(mean-c.val*sd),0),
                   ulim=target+const*(mean+c.val*sd),
                   level=factor(level))
        } else if (method=="indivEN.0meanapprox") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          res.indivEN <-
            indivEN.pois.0meanapprox(tbl.flag$stat, tbl.flag$E, target=target)
          tbl.flag <- tbl.flag %>%
            mutate(sd=sqrt(1+target*res.indivEN[2]*E),
                   stat=stat/sd, prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               const=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target) %>% sqrt %>% rep(each=n.level),
                               sd=tbl.flag %>% pull(sd) %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(llim=pmax(target-const*c.val*sd,0),
                   ulim=target+const*c.val*sd,
                   level=factor(level))
        } else if (method=="MO") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          phi <- pmax(tbl.flag %>% pull(stat) %>% `^`(2) %>% mean, 1)
          tbl.flag <- tbl.flag %>%
            mutate(stat=stat/sqrt(phi), prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               dev=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target*phi) %>% sqrt %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0), ulim=target+dev,
                   level=factor(level))
        } else if (method=="MOW") {
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          q.stat <- tbl.flag %>% pull(stat) %>% quantile(sort(q.winsor), T)
          tbl.flag <- tbl.flag %>%
            mutate(stat.w=pmax(pmin(stat, q.stat[2L]), q.stat[1L]))
          phi <- pmax(tbl.flag %>% pull(stat.w) %>% `^`(2) %>% mean, 1)
          tbl.flag <- tbl.flag %>%
            mutate(stat=stat/sqrt(phi), prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                               flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               dev=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target*phi) %>% sqrt %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0), ulim=target+dev,
                   level=factor(level))
        } else if (method=="MODW") {
          q.winsor <- sort(q.winsor)
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          q.stat <- tbl.flag %>% pull(stat) %>% quantile(q.winsor, T)
          tbl.flag <- tbl.flag %>%
            mutate(stat.w=pmax(pmin(stat, q.stat[2L]), q.stat[1L]))
          phi <- tbl.flag %>% pull(stat.w) %>% `^`(2) %>% mean
          q.lower <- qnorm(q.winsor[1L])
          q.upper <- qnorm(q.winsor[2L])
          var.inctrl <-
            1 + (q.lower*dnorm(q.lower)-q.upper*dnorm(q.upper)) /
            (pnorm(q.upper)-pnorm(q.lower)) -
            (dnorm(q.lower)-dnorm(q.upper))^2/(pnorm(q.upper)-pnorm(q.lower))^2
          var.outofctrl <-
            (q.lower^2-q.lower)*(q.winsor[1L])/sum(q.winsor[1L]+1-q.winsor[2L]) +
            (q.upper^2-q.upper)*(1-q.winsor[2L])/sum(q.winsor[1L]+1-q.winsor[2L])
          var <- (q.winsor[2L]-q.winsor[1L])*var.inctrl +
            sum(q.winsor[1L]+1-q.winsor[2L])*var.outofctrl
          phi <- pmax(phi/var, 1)
          tbl.flag <- tbl.flag %>%
            mutate(stat=stat/sqrt(phi), prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                               flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               dev=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target*phi) %>% sqrt %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0), ulim=target+dev,
                   level=factor(level))
        } else if (method=="AO") {
          sum.weight <- data %>% pull(E) %>% sum %>% `/`(target)
          sum.weightsq <- data %>% pull(E) %>% `^`(2) %>% sum %>% `/`(target^2)
          tbl.flag <- data %>%
            mutate(ind=O/E, stat=(O-target*E)/sqrt(target*E))
          tausq <- tbl.flag %>% pull(stat) %>% `^`(2) %>% sum %>% `-`(m-1) %>%
            `/`(sum.weight-sum.weightsq/sum.weight)
          tbl.flag <- tbl.flag %>%
            mutate(stat=(ind-target)/sqrt(tausq+target/E), prec=E,
                   prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                               flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                               dev=data %>% pull(E) %>% `^`(-1) %>%
                                 `*`(target) %>% `+`(tausq) %>%
                                 sqrt %>% rep(each=n.level),
                               level=levels %>% rep(times=m),
                               c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0), ulim=target+dev,
                   level=factor(level))
        }
      } else if (test=="exact") {
        tbl.flag <- data %>%
          mutate(ind=O/E, prec=E,
                 prob.lower=ppois(O,target*E)-0.5*dpois(O,target*E),
                 p=2*pmin(prob.lower, 1-prob.lower),
                 flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                             flags[2L]),
                 flag=factor(flag, flags))
        tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                             E=prec, target.E=target*E,
                             level=levels %>% rep(times=m)) %>%
          mutate(o.lower=qpois(level/2,target.E),
                 o.lower=ifelse(ppois(o.lower-1,target.E)+
                                  0.5*dpois(o.lower,target.E)>=level/2,
                                o.lower,o.lower+1),
                 lambda.lower=(dpois(o.lower,target.E)+
                                 2*ppois(o.lower-1,target.E)-level)/
                   (dpois(o.lower-1,target.E)+dpois(o.lower,target.E)),
                 llim=pmax(o.lower-lambda.lower,0)/E,
                 o.upper=qpois(1-level/2,target.E),
                 o.upper=ifelse(ppois(o.upper-1,target.E)+
                                  0.5*dpois(o.upper,target.E)>=1-level/2,
                                o.upper,o.upper+1),
                 lambda.upper=(dpois(o.upper,target.E)+
                                 2*ppois(o.upper-1,target.E)-2+level)/
                   (dpois(o.upper-1,target.E)+dpois(o.upper,target.E)),
                 ulim=(o.upper-lambda.upper)/E,
                 level=factor(level))
      }
      if (is.null(ylim)) ylim <- c(0, max(tbl.flag %>% pull(ind)))
    } else if (indicator=="DSR") { # directly standardized ratio
      if (target<0) stop("'target' is negative!")
    } else if (indicator=="Prop") {

    }
  } else if (grepl("^bern", tolower(outcome))) { # Bernoulli outcome
    if (!"data.frame"%in%class(data)) stop("Invalid class of 'data'!")
    if (indicator=="ISR") { # indirectly standardized ratio
      if (!all(c("O","E") %in% names(data))) stop("'O' or 'E' not in 'data'!")
      if (indiv.data) {
        ls.probs <- split(data %>% pull(E), data %>% pull(ID)) %>% unname
        data <- data %>% group_by(ID) %>%
          summarise(O=sum(O), E=sum(E), .groups="drop") %>%
          mutate(probs=ls.probs)
      } else {
        if (!"probs" %in% names(data)) stop("'probs' not in 'data'!")
        if (!"ID" %in% names(data))
          data <- data %>% mutate(ID=as.character(1:n()))
      }
      if (test=="score") {
        if (method=="FE") {
          tbl.flag <- data %>% rowwise() %>%
            mutate(probs=list(target*unlist(probs)),
                   V=max(sum(unlist(probs)*(1-unlist(probs))), 0.001),
                   N=length(unlist(probs)),
                   ind=O/E, stat=(O-target*E)/sqrt(V),
                   prec=E^2/V, prob.lower=pnorm(stat),
                   p=2*pmin(prob.lower, 1-prob.lower),
                   flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                               flags[2L]),
                   flag=factor(flag, flags))
          tbl.ctrlim <-
            tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                   N=tbl.flag %>% pull(N) %>% rep(each=n.level),
                   dev=1/sqrt(prec),
                   level=levels %>% rep(times=m),
                   c.val=qnorm(1-level/2)) %>%
            mutate(dev=c.val*dev, llim=pmax(target-dev,0),
                   ulim=pmin(target+dev,N),
                   level=factor(level))
        } else if (method=="indivEN.exact") {

        } else if (method=="indivEN.approx") {

        } else if (method=="indivEN.0meanapprox") {

        } else if (method=="MO") {

        } else if (method=="MOW") {

        } else if (method=="MODW") {

        } else if (method=="AO") {

        }
      } else if (test=="exact") {
        tbl.flag <- data %>% rowwise() %>%
          mutate(probs=list(target*unlist(probs)),
                 V=max(sum(unlist(probs)*(1-unlist(probs))), 0.001),
                 N=length(unlist(probs)),
                 ind=O/E, prec=E^2/V,
                 prob.lower=ppoibin(O,unlist(probs))-0.5*dpoibin(O,unlist(probs)),
                 p=2*pmin(prob.lower,1-prob.lower),
                 flag=ifelse(p<level,ifelse(ind>target,flags[3L],flags[1L]),
                             flags[2L]),
                 flag=factor(flag, flags))
        tbl.ctrlim <- tibble(prec=tbl.flag %>% pull(prec) %>% rep(each=n.level),
                             E=tbl.flag %>% pull(E) %>% rep(each=n.level),
                             N=tbl.flag %>% pull(N) %>% rep(each=n.level),
                             probs=tbl.flag %>% pull(probs) %>% rep(each=n.level),
                             level=levels %>% rep(times=m)) %>% rowwise() %>%
          mutate(o.lower=qpoibin(level/2,unlist(probs)),
                 o.lower=ifelse(ppoibin(o.lower-1,unlist(probs))+
                                  0.5*dpoibin(o.lower,unlist(probs))>=level/2,
                                o.lower,o.lower+1),
                 lambda.lower=(dpoibin(o.lower,unlist(probs))+
                                 2*ppoibin(o.lower-1,unlist(probs))-level)/
                   (dpoibin(o.lower-1,unlist(probs))+dpoibin(o.lower,unlist(probs))),
                 llim=pmax(o.lower-lambda.lower,0)/E,
                 o.upper=qpoibin(1-level/2,unlist(probs)),
                 o.upper=ifelse(ppoibin(o.upper-1,unlist(probs))+
                                  0.5*dpoibin(o.upper,unlist(probs))>=1-level/2,
                                o.upper,o.upper+1),
                 lambda.upper=(dpoibin(o.upper,unlist(probs))+
                                 2*ppoibin(o.upper-1,unlist(probs))-2+level)/
                   (dpoibin(o.upper-1,unlist(probs))+dpoibin(o.upper,unlist(probs))),
                 ulim=pmin(o.upper-lambda.upper,N)/E,
                 level=factor(level))
      }
    } else if (indicator=="DSR") { # directly standardized ratio

    } else if (indicator=="Prop") {

    }
  }
  labs.colorshape <- tbl.flag %>% pull(flag) %>% table %>%
    paste0(flags, "|", ., "|", round(prop.table(.)*100), "%")
  labs.linetype <- paste0((1-levels)*100,"%")
  vals.linetype <- character(n.level)
  vals.linetype[levels==level] <- 'solid'
  vals.linetype[levels!=level] <-
    c('dashed','dotted','dotdash','longdash','twodash')[1:(n.level-1)]
  ggplot() +
    geom_point(data=tbl.flag, aes(x=prec, y=ind, color=flag, shape=flag),
               size=ptsize) +
    scale_shape_manual(name=bquote("flagging "~(alpha==.(level))),
                       labels=labs.colorshape, values=c(15,17,19)) +
    scale_color_manual(name=bquote("flagging "~(alpha==.(level))),
                       labels=labs.colorshape, values=colors) +
    geom_line(data=tbl.ctrlim, aes(x=prec, y=llim, group=level, linetype=level),
              size=lwd) +
    geom_line(data=tbl.ctrlim, aes(x=prec, y=ulim, group=level, linetype=level),
              size=lwd) +
    scale_linetype_manual(name="ctrl limits", values=vals.linetype,
                          labels=labs.linetype) +
    scale_x_continuous(name=xlab, limits=xlim, breaks=xbreaks, expand=xexpand) +
    scale_y_continuous(name=ylab, limits=ylim, breaks=ybreaks, expand=yexpand) +
    guides(shape=guide_legend(order=1), color=guide_legend(order=1),
           linetype=guide_legend(reverse=T, order=2)) +
    geom_hline(yintercept=target, size=lwd, linetype="dashed") +
    labs(title=title) + theme_classic() + theme(...)
}
