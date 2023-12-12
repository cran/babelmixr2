## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(babelmixr2)
one.compartment <- function() {
  ini({
    tka <- log(1.57); label("Ka (1/hr)")
    tcl <- log(2.72); label("Cl (L/hr)")
    tv <- log(31.5); label("V (L)")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7; label("additive residual error (mg/L)")
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    vc <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / vc * center
    cp <- center / vc
    cp ~ add(add.sd)
  })
}

## ----model-update-------------------------------------------------------------
prepared <-
  nlmixr2(
    one.compartment,
    data = theo_sd,
    est = "pknca",
    control = pkncaControl(concu = "ng/mL", doseu = "mg", timeu = "hr", volumeu = "L")
  )

## ----examine-update-----------------------------------------------------------
prepared$ui

knitr::knit_print(
  summary(prepared$nca)
)

## ----fit----------------------------------------------------------------------
fit <- nlmixr(prepared, data = theo_sd, est = "focei", control = list(print = 0))

fit

## ----model-ncaData------------------------------------------------------------
# Choose a subset of the full dataset for NCA
dNCA <- theo_sd[theo_sd$ID <= 6, ]

preparedNcaData <-
  nlmixr2(
    one.compartment,
    data = theo_sd,
    est = "pknca",
    control = pkncaControl(concu = "ng/mL", doseu = "mg", timeu = "hr", volumeu = "L", ncaData = dNCA)
  )

preparedNcaData$ui

