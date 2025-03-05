# animations for slides

library(EvolutionaryGames)
library(tidyverse)

get_payoff <- function(
        # Strat_r: Receiving strat, Strat_p: Providing strat
    Strat_r, Strat_p, # expects rows of data frame with columns Eff, AIP, AIS
    p = list(
        # Parameters
        c = 4, # cost
        f = 5, # receive benefit - c <> b_r might depend on knowledge gap??
        b = 2, # self benefit
        q = 0.8, # LLM is 0.8 of decent feedback
        k = 0.1, # reduces cost of giving fb by this
        a = 1,
        gamma = 0
    )) {
    # Relabelling to match model formula
    Theta_R = Strat_r$Eff
    Theta_P = Strat_p$Eff
    # Gamma_R = Strat_r$AIS
    # if (!is.null(global_gamma)) {
    #     Gamma_R = global_gamma
    # }
    # Pi_R = Strat_r$AIP
    # Pi_P = Strat_p$AIP
    f <- p$f
    b <- p$b
    c <- p$c
    a <- p$a
    q <- p$q
    k <- p$k
    # global AI use
    Pi_R = p$gamma
    Pi_P = p$gamma
    Gamma_R = p$gamma
    
    V.f <- (1-Pi_P + q * Pi_P ) * Theta_P  + Gamma_R * q
    V.bc <- ((1 - Pi_R) + Pi_R * k) * Theta_R
    V.a <- 1 - abs(Theta_R - Theta_P)
    
    V = V.f * f + 
        V.bc * (b - c) + 
        V.a * a + 
        Gamma_R * (q * f - k * c)
    return(V)
}


# custom version of phase plot - VERY minor adjustment on colour scale from EvolutionaryGames package
# It's a hack to try and get the colours the same across the various plots
phase3d <- function (A, dynamic, params = NULL, trajectories = NULL, contour = FALSE, 
          vectorField = FALSE, strategies = c("1", "2", "3"), max_lim = 0.8) 
{
    if (!is.matrix(A) || !is.numeric(A)) {
        stop("A must be a numeric matrix.")
    }
    else if (nrow(A) != 3 || ncol(A) != 3) {
        stop("A must be of size 3x3.")
    }
    else if (!is.null(params) && !is.numeric(params)) {
        stop("params must be numeric.")
    }
    else if (!is.null(trajectories) && !is.numeric(trajectories)) {
        stop("trajectories must be a numeric matrix.")
    }
    else if (!is.null(trajectories) && ncol(trajectories) != 
             3) {
        stop("trajectories must be of size mx3.")
    }
    else if (!is.null(strategies) && length(strategies) != 3) {
        stop("Number of strategies does not match the number of columns of A.")
    }
    parameters <- as.vector(A)
    x <- y <- NULL
    if (!is.null(params)) {
        p <- as.vector(params)
        parameters <- c(parameters, p)
    }
    times <- seq(0, 100, by = 0.01)
    refSimp <- triangle()$coords
    odeData <- arrowData <- c()
    if (!is.null(trajectories)) {
        for (i in 1:nrow(trajectories)) {
            out <- deSolve::ode(y = trajectories[i, ], times = times, 
                                func = dynamic, parms = parameters)
            odeData <- rbind(odeData, out[, -1])
        }
        odeData <- geometry::bary2cart(refSimp, odeData)
        odeData <- data.frame(x = odeData[, 1], y = odeData[, 
                                                            2])
        arrNum <- 20
        dist <- length(times)/(arrNum + 1)
        step <- c()
        for (i in seq(1, nrow(odeData), length(times))) {
            for (j in 1:arrNum) {
                step <- c(step, i + j * dist)
            }
        }
        step2 <- step + 1
        xend <- yend <- NULL
        arrowData <- data.frame(x = c(odeData[step, 1]), y = c(odeData[step, 
                                                                       2]), xend = c(odeData[step2, 1]), yend = c(odeData[step2, 
                                                                                                                          2]))
    }
    maxVelocity <- 0
    density <- c()
    x <- y <- z <- c()
    for (i in seq(0, 1, by = 0.1)) {
        for (j in seq(0, 1, by = 0.1)) {
            if (i + j > 1) {
                break
            }
            x <- c(x, i)
            y <- c(y, j)
            z <- c(z, 1 - i - j)
            dX <- dynamic(state = c(i, j, 1 - i - j), parameters = parameters)[[1]]
            dist <- sqrt(dX[1]^2 + dX[2]^2 + dX[3]^2)
            density <- c(density, dist)
            maxVelocity <- max(dist, maxVelocity)
        }
    }
    contourData <- cbind(x, y, z)
    x1 <- y1 <- z1 <- c()
    x2 <- y2 <- z2 <- c()
    for (i in seq(0, 1, by = 0.05)) {
        for (j in seq(0, 1, by = 0.05)) {
            if (i + j > 1) {
                break
            }
            x1 <- c(x1, i)
            y1 <- c(y1, j)
            z1 <- c(z1, 1 - i - j)
            dX <- dynamic(state = c(i, j, 1 - i - j), parameters = parameters)[[1]]
            dist <- sqrt(dX[1]^2 + dX[2]^2 + dX[3]^2)
            x2 <- c(x2, dX[1] * dist)
            y2 <- c(y2, dX[2] * dist)
            z2 <- c(z2, dX[3] * dist)
        }
    }
    vecData1 <- cbind(x1, y1, z1)
    vecData2 <- cbind(x2, y2, z2)
    density <- density/maxVelocity
    contourData <- geometry::bary2cart(refSimp, contourData)
    vecData1 <- geometry::bary2cart(refSimp, vecData1)
    vecData2 <- geometry::bary2cart(refSimp, vecData2)
    vecData <- cbind(vecData1, vecData2)
    vecData <- data.frame(x = vecData[, 1], y = vecData[, 2], 
                          xend = vecData[, 1] + vecData[, 3] * 0.5, yend = vecData[, 
                                                                                   2] + vecData[, 4] * 0.5)
    if (any(vecData < 0)) {
        vecData$xend = vecData$x
        vecData$yend = vecData$y
    }
    resol <- 300
    contourData <- interp::interp(contourData[, 1], contourData[, 
                                                                2], density, seq(0, 1, length = resol), seq(0, 1, length = resol))
    contourData <- reshape2::melt(contourData$z)
    contourData[, 1:2] <- contourData[, 1:2]/resol
    contourData <- data.frame(x = contourData[, 1], y = contourData[, 
                                                                    2], z = contourData[, 3])
    p <- triangle(strategies)$canvas
    pal <- (grDevices::colorRampPalette(c("blue", "cyan", "green", 
                                          "yellow", "red")))(5)
    if (contour) {
        p <- p + ggplot2::geom_tile(data = contourData, ggplot2::aes(x = x, 
                                                                     y = y, fill = z)) + ggplot2::scale_fill_gradientn(
                                                                         colours = pal, 
                                                                         na.value = NA,
                                                                         limits = c(0, max_lim - maxVelocity)) + ggplot2::labs(fill = "Velocity")
    }
    if (!is.null(trajectories)) {
        p <- p + ggplot2::geom_point(data = odeData, ggplot2::aes(x = x, 
                                                                  y = y), size = 0.1, shape = 16) + ggplot2::geom_segment(data = arrowData, 
                                                                                                                          ggplot2::aes(x = x, y = y, xend = xend, yend = yend), 
                                                                                                                          arrow = ggplot2::arrow(length = ggplot2::unit(0.25, 
                                                                                                                                                                        "cm"), type = "closed"), size = 0)
    }
    if (vectorField) {
        p <- p + ggplot2::geom_segment(data = vecData, ggplot2::aes(x = x, 
                                                                    y = y, xend = xend, yend = yend), arrow = ggplot2::arrow(length = ggplot2::unit(0.1, 
                                                                                                                                                    "cm"), type = "closed"), size = 0.4)
    }
    print(p)
}

p = list(
    c = 4, # cost
    f = 5, # receive benefit
    b = 2, # self benefit
    q = 0.8, # Value from AI
    k = 0.1, # Cost of using AI
    a = 1 # assortment
)


# game modelling
all_strategies <- tribble(
    ~Label, ~Eff, ~AIS, ~AIP,
    "CN", 1,   0, 0,
    "TN", 0.5, 0, 0,
    "FN", 0,   0, 0,
    "CS", 1,   1, 0,
    "TS", 0.5, 1, 0,
    "FS", 0,   1, 0,
    "CO", 1,   0, 1,
    "TO", 0.5, 0, 1,
    "FO", 0,   0, 1,
    "CB", 1,   1, 1,
    "TB", 0.5, 1, 1,
    "FB", 0,   1, 1,
    "CO.8", 1, 0, 0.8,
    "CO.6", 1, 0, 0.6,
    "CO.4", 1, 0, 0.4,
    "CO.2", 1, 0, 0.2
)
build_payoff_matrix <- function(
        S,
        # Parameters
        p = list(
            c = 2, # cost
            f = 2, # receive benefit
            b = 1, # self benefit
            a = 1,
            q = 0.8, # LLM is 0.8 of decent feedback
            k = 0.1 # reduces cost of giving fb by this, but also reduces self benefit
        )
) {
    N_strategies <- nrow(S)
    payoffs <- array(rep(NA_real_, N_strategies*N_strategies), 
                     dim = c(N_strategies, N_strategies)) 
    for (R in 1:N_strategies) {
        for (P in 1:N_strategies) {
            payoffs[R,P] <- get_payoff(slice(S, R), slice(S, P), p)
        }
    }
    return(payoffs)
}

print_payoff_matrix <- function(M, stgys) {
    rownames(M) <- stgys
    colnames(M) <- stgys
    print(M)
}

fetch_payoff_matrix <- function(strategies, p, as_df = FALSE,
                                add_exp_V = FALSE) {
    # strategies: character vector, p: list
    S.df <- tibble(Label = strategies) |> # preserves order of strategies
        inner_join(all_strategies, by = "Label")
    payoff_matrix <- build_payoff_matrix(S.df, p)
    if (as_df) {
        colnames(payoff_matrix) <- strategies
        payoff_df <- as_tibble(payoff_matrix) |> 
            mutate(`Payoff to strategy:` = strategies) |> 
            select(`Payoff to strategy:`, everything())
        if (add_exp_V) {
            payoff_df <- payoff_df |> 
                rowwise() |> 
                mutate(ExpV = mean(c_across(where(is.numeric))))
        }
        return(payoff_df)
    } else {
        return(payoff_matrix)
    }
}

state <- matrix(c(0.4, 0.3, 0.3), 1, 3, byrow=TRUE)

# Produces ggplot!!
plot_phases <- function(p.p, heat = TRUE,
                        genAI = TRUE, max_lim = 5) {
    if (genAI) {
        strats <- c("CB", "TB", "FB")
    } else {
        strats <- c("CN", "TN", "FN")
    }
    payoff <- fetch_payoff_matrix(
        strats, 
        p = p.p)
    phase3d(payoff, Replicator, NULL, NULL, heat, TRUE,
               strategies = c("Cooperator", "Token-effort", "Free-rider"),
            max_lim = max_lim) +
        guides(fill = "none") +
        geom_text(aes(x = 0, y = 0.9, label = "Environment\nparameters")) +
        geom_text(aes(x = -0.08, y = 0.8, label = "Cost:")) +
        geom_linerange(aes(y = 0.8, xmin = 0.02, xmax = 0.02+p.p$c / 30),
                       linewidth = 3, color = "darkred") +
        geom_text(aes(x = -0.12, y = 0.7, label = "Self benefit:")) +
        geom_linerange(aes(y = 0.7, xmin = 0.02, xmax = 0.02+p.p$b / 30),
                       linewidth = 3, color = "darkblue") +
        geom_text(aes(x = -0.14, y = 0.6, label = "GenAI use:")) +
        geom_linerange(aes(y = 0.6, xmin = 0.02, xmax = 0.02+p.p$gamma / 5),
                       linewidth = 3, color = "darkgreen")
        
}    


# Creating a list of parameter settings
N <- 31
Increment <- 4 / (N-1) 
pstart <- list(
    c = 4, # cost
    f = 5, # receive benefit
    b = 2, # self benefit
    q = 2, # Value from AI
    k = 0.1, # Cost of using AI
    a = 1.5, # assortment
    gamma = 0.0 # switching AI use to global param
)
plist <- list()
for (i in 1:N) { # up b
    plist[[i]] <- pstart
    pstart$b <- pstart$b + Increment
}
pstart <- plist[[N]]
for (i in 1:N) { # up gamma
    plist[[N+i]] <- pstart
    pstart$gamma <- pstart$gamma + 1 / (N+1)
}
pstart <- plist[[2*N]]
for (i in 1:N) { # down b
    plist[[i + 2*N]] <- pstart
    pstart$b <- pstart$b - Increment
}
pstart <- plist[[3*N]]
for (i in 1:N) { # down gamma
    plist[[i + 3*N]] <- pstart
    pstart$gamma <- pstart$gamma - 1 / (N+1)
}

plot_phases(plist[[1]], heat = T, genAI = T, max_lim = 1.9)

animate_parameter_list <- function(plist, filename, 
                                   heat = T, genAI = T,
                                   fps = 10,
                                   MAX_lim = 2.1) {
    N <- length(plist)
    ggplotlist <- list()
    plotlist <- list()
    for (i in 1:N) { # slow bit
        img_path <- paste0("temp_plot_", i, ".png")
        ggsave(img_path, 
               plot_phases(plist[[i]], heat = heat, genAI = genAI, max_lim = MAX_lim),
               width = 5, height = 4, dpi = 300)
        plotlist[[i]] <- img_path
    }
    
    img_list <- image_read(as.character(plotlist))
    animation <- image_animate(image_join(img_list), fps = fps,
                               optimize = T)
    image_write(animation, filename)
    # removing plots
    file.remove(as.character(plotlist))
}

# animate_parameter_list(plist, "annoai.gif", T, F)

animate_parameter_list(plist, "angenai.gif", T, T, fps = 25, 2)
animate_parameter_list(plist, "angenai.long.gif", T, T, fps = 20, 1.9)
