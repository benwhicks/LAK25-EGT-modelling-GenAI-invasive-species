# animations for slides

library(EvolutionaryGames)
library(tidyverse)

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
            q = 0.8, # LLM is 0.8 of decent feedback
            k = 0.1 # reduces cost of giving fb by this, but also reduces self benefit
        )
) {
    N_strategies <- nrow(S)
    payoffs <- array(rep(NA_real_, N_strategies*N_strategies), dim = c(N_strategies, N_strategies)) 
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
                        genAI = TRUE) {
    if (genAI) {
        strats <- c("CO", "TO", "FO")
    } else {
        strats <- c("CN", "TN", "FN")
    }
    payoff <- fetch_payoff_matrix(
        strats, 
        p = p.p)
    phaseDiagram3S(payoff, Replicator, NULL, NULL, heat, TRUE,
               strategies = c("Cooperator", "Token-effort", "Free-rider")) +
        guides(fill = "none") +
        geom_text(aes(x = 0, y = 0.9, label = "Environment\nparameters")) +
        geom_text(aes(x = -0.08, y = 0.8, label = "Cost:")) +
        geom_linerange(aes(y = 0.8, xmin = 0, xmax = p.p$c / 30),
                       linewidth = 3, color = "darkred") +
        geom_text(aes(x = -0.12, y = 0.7, label = "Self benefit:")) +
        geom_linerange(aes(y = 0.7, xmin = 0, xmax = p.p$b / 30),
                       linewidth = 3, color = "blue") +
        geom_text(aes(x = -0.14, y = 0.6, label = "Social benefit:")) +
        geom_linerange(aes(y = 0.6, xmin = 0, xmax = p.p$a / 30),
                       linewidth = 3, color = "darkgreen")
        
}    
plot_phases(plist[[20]], F, F)

# Creating a list of parameter settings
N <- 41
Increment <- 0.1
pstart <- p
plist <- list()
for (i in 1:N) { # up b
    plist[[i]] <- pstart
    pstart$b <- pstart$b + Increment
}
pstart <- plist[[N]]
for (i in 1:N) { # up a
    plist[[N+i]] <- pstart
    pstart$a <- pstart$a + Increment/2
}
pstart <- plist[[2*N]]
for (i in 1:N) { # down b
    plist[[i + 2*N]] <- pstart
    pstart$b <- pstart$b - Increment
}
pstart <- plist[[3*N]]
for (i in 1:N) { # down a
    plist[[i + 3*N]] <- pstart
    pstart$a <- pstart$a - Increment/2
}

animate_parameter_list <- function(plist, filename, 
                                   heat = T, genAI = T) {
    N <- length(plist)
    ggplotlist <- list()
    plotlist <- list()
    for (i in 1:N) { # slow bit
        img_path <- paste0("temp_plot_", i, ".png")
        ggsave(img_path, 
               plot_phases(plist[[i]], heat = heat, genAI = genAI),
               width = 5, height = 4, dpi = 300)
        plotlist[[i]] <- img_path
    }
    
    img_list <- image_read(as.character(plotlist))
    animation <- image_animate(image_join(img_list), fps = 5,
                               optimize = T)
    image_write(animation, filename)
    # removing plots
    file.remove(as.character(plotlist))
}

animate_parameter_list(plist, "annn.gif", T, F)
animate_parameter_list(plist, "anai.gif", T, T)

library(magick)



saveGIF(
    expr = {
        walk(
            plotlist,
            ~ plot(get(.))
        )
    },
    movie.name = "implicit_my3.gif"
)

saveGIF(
    expr = plotlist,
    movie.name = "basic.gif"
)

plot_phases(p, T, F)    

    geom_text(aes(label = "Text", x = 1, y = 0.5)) + 
    geom_poi

