##' Plot MBNL1 dosing data and fit the Hill equation (PRISM - http://www.graphpad.com/guides/prism/6/curve-fitting/index.htm?reg_hill_slope.htm) Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
##'
##' Input dataframe(s) should have the following columns: Dox, PSI. Dox - doxycycline treatment in ng/ul; PSI - percent spliced in
##' @title FitDosingData
##' @param df1 data.frame object
##' @param df2 data.frame object
##' @param output_prefix prefix to use for plot output file
##' @return plotted data in eps format and the fit parameters: EC50, HillSlope
##' @author Adam Struck
FitDosingData <- function(df1, df2 =  NULL, output_prefix) {
    require(nls2)
    require(ggplot2)
    require(graphics)

    HillEq <- PSI ~ min(PSI) + (max(PSI) - min(PSI)) / (1 + 10^( (LogEC50 - logDox) * HillSlope ))

    name.df1 <- as.character(as.list(match.call())[[2]])
    names(df1) <- c("Dox", "PSI")
    df1$logDox <- log10(df1$Dox)
    df1 <- df1[is.finite(df1$logDox),]
    df1.plot <- data.frame(logDox = unique(df1$logDox),
                           PSI = sapply(unique(df1$logDox), function(x) { mean(df1[which(df1$logDox == x),"PSI"], na.rm = TRUE) }),
                           SD = sapply(unique(df1$logDox), function(x) { sd(df1[which(df1$logDox == x),"PSI"], na.rm = TRUE) }))
    
    st1 <- expand.grid(LogEC50 = seq(-0.5, 1.5, len = 10), HillSlope = seq(-10, 10, len = 10))
    mod1 <- nls2(HillEq, start = st1, algorithm = "brute-force", data = df1.plot) 
    fit <- nls2(HillEq, start = mod1, algorithm = "default", data = df1.plot)
    logdox.fit <- seq( min(df1$logDox), max(df1$logDox), length = length(df1$logDox) )
    PSI.fit <- predict( fit, data.frame(logDox = logdox.fit) )

    if( !is.null(df2) ) {
        name.df2 <- as.character(as.list(match.call())[[3]])
        names(df2) <- c("Dox", "PSI")
        df2$logDox <- log10(df2$Dox)
        df2 <- df2[is.finite(df2$logDox),]
        df2.plot <- data.frame(logDox = unique(df2$logDox),
                               PSI = sapply(unique(df2$logDox), function(x) { mean(df2[which(df2$logDox == x),"PSI"], na.rm = TRUE) }),
                               SD = sapply(unique(df2$logDox), function(x) { sd(df2[which(df2$logDox == x),"PSI"], na.rm = TRUE) }))
        
        mod2 <- nls2(HillEq, start = st1, algorithm = "brute-force", data = df2.plot)
        fit2 <- nls2(HillEq, start = mod2, algorithm = "default", data = df2.plot)
        logdox.fit2 <- seq( min(df2$logDox), max(df2$logDox), length = length(df2$logDox) )
        PSI.fit2 <- predict( fit2, data.frame(logDox = logdox.fit2) )
        df12 <- data.frame(df1.plot, df2.plot)
        
        plot.result <- ggplot(data = df12, environment=environment()) +
            geom_point( aes( x=logDox, y=PSI ), colour="blue", size=2.5 ) +
            geom_line( aes( x=logdox.fit, y=PSI.fit, colour="name.df1") ) +
            geom_errorbar( aes(x=logDox, ymax=PSI+SD, ymin=PSI-SD, width=0.025 ), colour="blue" ) +
            geom_point( aes( x=logDox.1, y=PSI.1 ), colour="red", size=2.5 ) +
            geom_line( aes( x=logdox.fit2, y=PSI.fit2, colour="name.df2") ) +
            geom_errorbar( aes(x=logDox.1, ymax=PSI.1+SD.1, ymin=PSI.1-SD.1, width=0.025 ), colour="red" ) +
            ylim(0,1) + xlab("log( Dox )") + ylab("PSI") + theme_gray(base_size=16) +
            scale_color_manual(labels = c(name.df1, name.df2), values = c("name.df1" = "blue", "name.df2" = "red"), name = "")

        result2 <- cbind(as.data.frame(rbind(coef(fit), coef(fit2))),
                         data.frame(DeltaPSI=c(df1.plot$PSI[length(df1.plot$PSI)] - df1.plot$PSI[1], df2.plot$PSI[length(df2.plot$PSI)]-df2.plot$PSI[1]),
                                    Bottom=c(df1.plot$PSI[1], df2.plot$PSI[1]),
                                    Top=c(df1.plot$PSI[length(df1.plot$PSI)], df2.plot$PSI[length(df2.plot$PSI)]) ))
        result2$EC50 <- 10^result2$LogEC50
        rownames(result2) <- c(name.df1, name.df2)
        result <- list(name.df1=coef(summary(fit)), name.df2=coef(summary(fit2)))
        names(result) <-  c(name.df1, name.df2)
        
    } else {
        plot.result <- ggplot(data=df1.plot, environment=environment()) +
            geom_point( aes( x=df1.plot$logDox, y=df1.plot$PSI ), colour="blue", size=2.5 ) +
            geom_line( aes( x=logdox.fit, y=PSI.fit, colour="name.df1") ) +
            geom_errorbar( aes(x=df1.plot$logDox, ymax=df1.plot$PSI+df1.plot$SD, ymin=df1.plot$PSI-df1.plot$SD, width=0.025 ), colour="blue" ) +
            ylim(0,1) + xlab("log( Dox )") + ylab("PSI") + theme_gray(base_size=16) +
            scale_color_manual(labels = c(name.df1), values = c("name.df1" = "blue"), name = "")

        result2 <- cbind(as.data.frame(t(coef(fit))),
                         data.frame(DeltaPSI=c(df1.plot$PSI[length(df1.plot$PSI)] - df1.plot$PSI[1]),
                                    Bottom=df1.plot$PSI[1],
                                    Top=df1.plot$PSI[length(df1.plot$PSI)]) )
        result2$EC50 <- 10^result2$LogEC50
        rownames(result2) <- c(name.df1)
        result <- list(name.df1=coef(summary(fit)))
        names(result) <-  c(name.df1)
    }
    postscript(sprintf("%s.ps", output_prefix), width=9, height=7)
    print( plot.result )
    graphics.off()
    return(result2)
}
