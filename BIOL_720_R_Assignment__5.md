Question 2
==========

    diploid_sim <- function(p0 = 0.1, w1 = 1, w2 = 0.95, w3 = 0.9,  n = 500) {
      
        p <- rep(NA,n)  
        w_bar <- rep(NA, n)
         
        p[1] <- p0 

          w_bar[1] <- (p[1]^2*w1) + ((2*p[1]*(1-p[1]))*w2) + ((1-p[1])^2)*w3
        
        for (i in 2:n) {
            w_bar[i - 1] <- (p[i-1]^2*w1) + ((2*p[i-1]*(1-p[i-1]))*w2) + ((1-p[i-1])^2)*w3 
            p[i] <- ((p[i - 1]^2) * (w1/w_bar[i - 1])) + (p[i-1]*(1-p[i-1]))*(w2/w_bar[i-1])
        }
          
        if (p > 0.9999) {
        fixation <- min(which.max(p > 0.9999))
        print("fixation for A1 occurs approximately at generation", fixation )  
        } else {
            maxAlleleFreq <- max(p)
            print("fixation of A1 does not occur, max. allele frequency is")
            print(maxAlleleFreq, digits = 2) 
        }
          
    plot(x = 1:n, y = p, 
         xlab="Generations", 
         ylab="Allele frequency (p)", 
         pch = 20, col = "red", cex.lab = 1.5)
    }


    diploid_sim()

    ## Warning in if (p > 0.9999) {: the condition has length > 1 and only the
    ## first element will be used

    ## [1] "fixation of A1 does not occur, max. allele frequency is"
    ## [1] 1

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-1-1.png)

Question 3
==========

    gen_drift <- function(a, HD, ind, gen) {
      alleles <- c(sample(c("a", "A"), 
                          size = HD*(ind), 
                          replace = TRUE, 
                          prob = c(a, 1 - a)))

      allelefreq <- matrix(data = NA, nrow = gen, ncol = 2)

      for(i in 1:gen){
        alleles <- sample(alleles, HD*ind, replace=TRUE)
        allelefreq[i,2] <- (length(alleles[alleles =="A"]))/(HD*(ind))
        allelefreq[i,1] <- (length(alleles[alleles =="a"]))/(HD*(ind))
      }

     a <- allelefreq[,1]
     A <- allelefreq[,2]

     plot(1, type="n", xlim=c(1,gen), ylim=c(0,1), xlab="generation", ylab="Frequency Of Alleles")
     lines(x = A, type = "l", col = "red")
     lines(x = a, type = "l", col = "blue")
     legend("topleft", legend=c("A Frequency", "a Frequency"),
           col=c("red", "blue"), lty=1, cex=0.8)
    }

    gen_drift(0.5, 2, 20, 100)

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-2-1.png)

Question 4
==========

    allele_lost <- function(A, gen, runs) {
      
    Aprop1 <- matrix(data = NA, nrow = runs, ncol = 1)

    for(j in 1:runs){
    alleles <- c(sample(c("a", "A"), size = 400, replace = TRUE, prob = c(1 - A, A)))
    Aprop <- integer()


    for(i in 1:gen) {
    alleles <- sample(alleles, 400, replace=TRUE)
    Afreq <- (length(alleles[alleles == "A"]))/(400)
    Aprop[i] <- Afreq
    }

    Aprop1[j,1] <- tail(Aprop, n = 1)
    }

    zerocount <- (length(which(Aprop1[,1] == "0")))
    print("Proportion of simulation where pA = 0")
    print(zerocount/runs)
    }


    allele_lost(0.5, 100, 1000)

    ## [1] "Proportion of simulation where pA = 0"
    ## [1] 0.006

Question 5
==========

    allele_traject <- function(A, gen, runs) {

    plot(1, type = "n", xlim = c(1,gen), ylim = c(0,1), xlab = "Generation", ylab = "Allele Frequency", main = "Influence of Genetic Drift on Allele Frequencies")

    for(j in 1:runs){
    alleles <- c(sample(c("a", "A"), size = 400, replace = TRUE, prob = c(1 - A, A)))
    Aprop <- integer()


    for(i in 1:gen) {
    alleles <- sample(alleles, 400, replace=TRUE)
    Afreq <- (length(alleles[alleles =="A"]))/(400)
    Aprop[i] <- Afreq
    }
    lines(x = Aprop, type = "l", col = sample(rainbow(n= runs)), lwd=2)

    }

    }

    allele_traject(0.5, 100, 100)

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-4-1.png)

Question 6
==========

from question (stochastic):
===========================

x &lt;- seq(from =1, to = 10, length.out = 20) \# length.out is how many
observations we will have a &lt;- 0.5 \# intercept b &lt;- 0.1 \# slope
y\_deterministic &lt;- a + b\*x

y\_simulated &lt;- rnorm(length(x), mean = y\_deterministic, sd = 2)

mod\_sim &lt;- lm(y\_simulated ~ x) p\_val\_slope &lt;-
summary(mod\_sim)$coef\[2,4\] \# extracts the p-value p\_val\_slope

    p_val <- function(a, b, n, rse){
      
    set.seed(720)
    x <- seq(from =1, to = 10, length.out = n)
    y_deterministic <- a + b*x
    pvals <- integer(length = 1000)

    y_simulated <-rnorm(length(x), mean = y_deterministic, sd = rse)
    mod_sim <- lm(y_simulated ~ x)
    print(summary(mod_sim)$coef[2,4])

    }

    p_val(0.5, 0.1, 20, 3)

    ## [1] 0.4640427

    # to confirm it works, put set.seed() in the function

    # Run 1000 times

    hist_thousand <- function(a, b, n, rse){
      
    x <- seq(from =1, to = 10, length.out = n)
    y_deterministic <- a + b*x
    pvals <- integer(length = 1000)

    for(i in 1:1000) {
    y_simulated <-rnorm(length(x), mean = y_deterministic, sd = rse)
    mod_sim <- lm(y_simulated ~ x)
    p_val_slope <- summary(mod_sim)$coef[2,4]
    pvals[i] <- p_val_slope 
    }

    hist(pvals, ylim = c(0,200), main = "Frequency of P Values", xlab = "P Values")
    punder0.05 <- which(pvals < 0.05)
    print("Proportion of P-Values Less than 0.05")
    print(length(punder0.05)/1000)
    }

    hist_thousand(0.5, 0.1, 20, 1.5)

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    ## [1] "Proportion of P-Values Less than 0.05"
    ## [1] 0.131

    # Change to slope = 0

    hist_thousand(0.5, 0, 20, 1.5)

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    ## [1] "Proportion of P-Values Less than 0.05"
    ## [1] 0.048

    # The histogram is more evenly distributed. Around 50/50 for above or below P = 0.05
    # There was a much greater frequency of P values under 0.05 when slope = 0.1
    # This is because the slope & intercept make the points for the line fitting the curve. Setting slope = 0 fits the points to a straight line & no relationship between x & y is observed. The P values show nothing significant in this case.

    # slope = 0.1, a = 0.5, rse = 1.5
    # for loop

    pval_loop <- function(a, b, rse){
      
    pvals <- integer()
    forplot <- matrix(data = NA, nrow = 19, ncol = 2)
    n=10

    for(j in 1:10000){

    x <- seq(from =1, to = 10, length.out = n)
    y_deterministic <- a + b*x
      
    for(i in 1:100) {
    y_simulated <-rnorm(length(x), mean = y_deterministic, sd = rse)
    mod_sim <- lm(y_simulated ~ x)
    p_val_slope <- summary(mod_sim)$coef[2,4]
    pvals[i] <- p_val_slope
    }

    punder0.05 <- which(pvals < 0.05) # count ones under 0.05
    forplot[j,2] <- (length(punder0.05))/100
    forplot[j,1] <- n
    n =  n + 5
    if( n > 100){ 
    break
    }
    }
    print(forplot)
    plot( x = forplot[,1], y = forplot[,2], type = "p", xlab = "Sample Size", ylab = "Proportion < 0.05", main = "Effect of Sample size on P Value", col = "red", pch = 16)
    }

    pval_loop(0.5, 0.1, 1.5)

    ##       [,1] [,2]
    ##  [1,]   10 0.08
    ##  [2,]   15 0.09
    ##  [3,]   20 0.15
    ##  [4,]   25 0.10
    ##  [5,]   30 0.12
    ##  [6,]   35 0.22
    ##  [7,]   40 0.22
    ##  [8,]   45 0.22
    ##  [9,]   50 0.24
    ## [10,]   55 0.31
    ## [11,]   60 0.29
    ## [12,]   65 0.39
    ## [13,]   70 0.23
    ## [14,]   75 0.36
    ## [15,]   80 0.23
    ## [16,]   85 0.43
    ## [17,]   90 0.31
    ## [18,]   95 0.41
    ## [19,]  100 0.41

![](BIOL_720_R_Assignment__5_files/figure-markdown_strict/unnamed-chunk-8-1.png)
