# lm_eqn is a copy from https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
              r2 = format(sqrt(summary(m)$r.squared), digits = 3)))
    as.character(as.expression(eq));
}

do_goset <- funciton(universe, selected) {}

